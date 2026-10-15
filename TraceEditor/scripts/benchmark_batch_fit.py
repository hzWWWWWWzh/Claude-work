#!/usr/bin/env python
"""Benchmark batch fitting: sequential loop vs core.batch_fit.run_batch's
ProcessPoolExecutor, on synthetic ROIs sized like real data.

Gives a reproducible, re-runnable wall-clock comparison instead of a one-off
manual timing, and doubles as a regression guard: if the pooled path ever
accidentally serializes itself (e.g. a bug drops it back to `max_workers=1`,
or the O(n^2) `/Parameters` write bottleneck this was built to avoid comes
back), the reported speedup will drop toward — or below — 1x here, which
correctness tests alone wouldn't catch (they don't assert on wall-clock
time).

Two scenarios, both real models from force_spectroscopy.py:

- `ExponentialDrift`: cheap per-fit (~ms). At realistic batch sizes the pool
  LOSES to the sequential loop — each worker process pays a real, if now
  reduced, one-time import cost that a cheap fit can't amortize. Included so
  this stays honest: multiprocessing is not a universal win here.
- `zHMf2`: the heaviest model in force_spectroscopy.py (calls zHMf four
  times per evaluation for bilinear interpolation over non-integer
  n_dimer/n_gap). At real_data.h5's actual scale (2701 points/ROI, 51 ROIs,
  5 free params) a single fit takes several seconds, so the pool wins
  decisively — this is the scenario that justifies run_batch existing.

    uv run python scripts/benchmark_batch_fit.py [scenario] [n_rois] [max_workers]
    uv run python scripts/benchmark_batch_fit.py zhmf2 51
    uv run python scripts/benchmark_batch_fit.py drift 150
"""

import sys
import time
from pathlib import Path

import h5py
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "tests"))

from conftest import make_synthetic_workflow, write_workflow_h5  # noqa: E402

from trace_editor.core import force_spectroscopy as fs  # noqa: E402
from trace_editor.core.batch_fit import run_batch, _fit_one_roi  # noqa: E402
from trace_editor.core.h5_data import prepare_fit  # noqa: E402


def _drift_params():
    _, params, names = prepare_fit("ExponentialDrift")
    params["Z0"].set(value=0.0, vary=True)
    params["dZ_dt"].set(value=0.0, vary=True)
    params["tau"].set(value=1e9, vary=False)
    return params, names


def _build_drift_file(tmp_dir: Path, n_rois: int) -> Path:
    beads = [make_synthetic_workflow(noise_um=0.05, seed=i, Z0=10.0 + i, dZ_dt=-1.0 - 0.01 * i)
             for i in range(n_rois)]
    return write_workflow_h5(tmp_dir / "benchmark_drift.h5", beads)


def _zhmf2_params():
    _, params, names = prepare_fit("zHMf2")
    names["x"] = "F (pN)"
    params["L"].set(value=3646, vary=False)
    params["A"].set(value=50, vary=False)
    params["S"].set(value=1000, vary=False)
    params["k"].set(value=0.5, vary=True, min=0.01, max=5)
    params["n_dimer"].set(value=90.0, vary=True, min=0, max=300)
    params["n_gap"].set(value=1.0, vary=True, min=0, max=100)
    params["g1"].set(value=1.9, vary=True, min=0, max=10)
    params["g2"].set(value=4.11, vary=True, min=0, max=10)
    return params, names


def _build_zhmf2_file(tmp_dir: Path, n_rois: int, n_samples: int = 2701) -> Path:
    """Write a file sized like real_data.h5 (2701 pts/ROI), with Z traces
    generated from the real zHMf2 model so fits have something to converge to.
    """
    path = tmp_dir / "benchmark_zhmf2.h5"
    force = np.linspace(0.5, 60, n_samples)
    with h5py.File(path, "w") as f:
        f.attrs["Date"] = "2026-07-13"
        f.attrs["User"] = "benchmark"
        f.attrs["Notes"] = "zHMf2 benchmark (synthetic, sized like real_data.h5)"
        f.attrs["axis"] = "Shift_mm"
        f.attrs["num_rois"] = n_rois
        f.attrs["roi_size"] = 64
        f.attrs["bin_roi_size"] = 64

        param_dtype = np.dtype([("ROI", np.int32), ("X (pix)", np.float32),
                                 ("Y (pix)", np.float32), (" selected", np.float64)])
        params = np.zeros(n_rois, dtype=param_dtype)
        for i in range(n_rois):
            params[i]["ROI"] = i
            params[i][" selected"] = 1.0
        f.create_dataset("Parameters", data=params)

        traces = f.create_group("Traces")
        traces.create_dataset("Time (s)", data=np.linspace(0, 60, n_samples))
        traces.create_dataset("Shift (mm)", data=np.linspace(10.8, 1.0, n_samples))
        traces.create_dataset("F (pN)", data=force)

        z_grp = traces.create_group("Z (um)")
        for i in range(n_rois):
            rng = np.random.default_rng(i)
            z_true = np.asarray(fs.zHMf2(force, L=3646, n_dimer=95.0 + i * 0.3,
                                          n_gap=2.0, k=0.5))
            z = z_true + rng.normal(0, 0.01, size=n_samples)
            z_grp.create_dataset(f"{i:03d}", data=z.astype(np.float32))

    return path


SCENARIOS = {
    "drift": {
        "model": "ExponentialDrift",
        "build_file": _build_drift_file,
        "params": _drift_params,
        "default_n": 150,
        "note": "cheap model: expect the pool to LOSE (import overhead dominates)",
    },
    "zhmf2": {
        "model": "zHMf2",
        "build_file": _build_zhmf2_file,
        "params": _zhmf2_params,
        "default_n": 51,   # matches real_data.h5's actual ROI count
        "note": "heaviest model, real_data.h5 scale: expect the pool to WIN",
    },
}


def main():
    import tempfile

    scenario_name = sys.argv[1] if len(sys.argv) > 1 else "zhmf2"
    scenario = SCENARIOS[scenario_name]
    n_rois = int(sys.argv[2]) if len(sys.argv) > 2 else scenario["default_n"]
    max_workers = int(sys.argv[3]) if len(sys.argv) > 3 else None

    with tempfile.TemporaryDirectory() as tmp:
        path = scenario["build_file"](Path(tmp), n_rois)
        tasks = [(str(path), f"{i:03d}") for i in range(n_rois)]
        params, names = scenario["params"]()
        model_name = scenario["model"]

        print(f"Scenario: {scenario_name} ({scenario['note']})")
        print(f"  model={model_name}, n_rois={n_rois}, max_workers={max_workers or 'default'}")

        t0 = time.perf_counter()
        for filepath, roi in tasks:
            _fit_one_roi(filepath, roi, model_name, params, names)
        t_seq = time.perf_counter() - t0
        print(f"  sequential: {t_seq:7.2f} s  ({t_seq / n_rois * 1000:.1f} ms/roi)")

        t0 = time.perf_counter()
        run_batch(tasks, model_name, params, names, max_workers=max_workers)
        t_pool = time.perf_counter() - t0
        print(f"  pooled:     {t_pool:7.2f} s")

        print(f"  speedup:    {t_seq / t_pool:5.2f}x")


if __name__ == "__main__":
    main()
