"""Unit tests for core.batch_fit (Qt-free — no QApplication/qtbot needed).

Exercises the multiprocessing batch-fit path directly, mirroring the same
sequential fit the GUI's single-ROI `_run_fit` performs, so these tests can
run without ever importing PySide6.
"""

import multiprocessing

import numpy as np
import pytest
from lmfit import Parameters

from trace_editor.core.batch_fit import (
    run_batch, RoiFitResult, BatchFitError, _fit_one_roi, model_wants_multiprocessing,
)
from trace_editor.core.h5_data import H5Data, prepare_fit

from conftest import write_workflow_h5, make_synthetic_workflow, TRUE_L_BP


pytestmark = pytest.mark.fitting


def _drift_params():
    _, params, names = prepare_fit("ExponentialDrift")
    params["Z0"].set(value=0.0, vary=True)
    params["dZ_dt"].set(value=0.0, vary=True)
    params["tau"].set(value=1e9, vary=False)
    return params, names


class TestFitOneRoi:
    """The picklable per-ROI worker function, called directly (no pool)."""

    def test_recovers_drift_params(self, tmp_path):
        """dZ_dt recovers exactly; Z0 legitimately absorbs the mean WLC height
        at the selected low-force points (see tests/README.md's z_bead
        invariant) so it is not expected to equal the injected Z0 here.

        Restricted to the low-force plateau (Shift (mm) >= p90) via explicit
        criteria, matching how the GUI actually drives this fit --
        ExponentialDrift's model only recovers cleanly in that low-force
        region (see TestDefaultSelectionCriteria).
        """
        criteria = [{"channel": "Shift (mm)", "min": "p90"}]

        beads = [make_synthetic_workflow(noise_um=0.0, seed=0, Z0=12.0, dZ_dt=-1.5)]
        path = write_workflow_h5(tmp_path / "one.h5", beads)
        params, names = _drift_params()

        result = _fit_one_roi(str(path), "000", "ExponentialDrift", params, names,
                               criteria=criteria)

        assert isinstance(result, RoiFitResult)
        assert result.param_values["dZ_dt (nm_s)"][0] == pytest.approx(-1.5, abs=1e-1)
        assert np.isfinite(result.param_values["Z0 (um)"][0])
        assert result.r2 > 0.5
        assert np.isfinite(result.r2_full)

    def test_raises_on_bad_roi(self, tmp_path):
        beads = [make_synthetic_workflow(noise_um=0.0, seed=0)]
        path = write_workflow_h5(tmp_path / "one.h5", beads)
        params, names = _drift_params()

        with pytest.raises(Exception):
            _fit_one_roi(str(path), "999", "ExponentialDrift", params, names)


class TestRunBatch:
    """run_batch: the ProcessPoolExecutor-backed batch entry point."""

    def _make_file(self, tmp_path, n_beads=4, noise_um=0.0):
        beads = [make_synthetic_workflow(noise_um=noise_um, seed=i,
                                         Z0=10.0 + i, dZ_dt=-1.0 - 0.1 * i)
                 for i in range(n_beads)]
        return write_workflow_h5(tmp_path / "batch.h5", beads)

    def test_empty_tasks_returns_immediately(self):
        assert run_batch([], "ExponentialDrift", *_drift_params()) == []

    def test_fits_all_rois_and_writes_to_disk(self, tmp_path):
        """dZ_dt (unaffected by the Z0/WLC-height offset ambiguity, see
        tests/README.md's z_bead invariant) recovers per-ROI; Z0 is written
        and finite but not compared to the injected value for that reason.
        """
        path = self._make_file(tmp_path, n_beads=3)
        tasks = [(str(path), f"{i:03d}") for i in range(3)]
        params, names = _drift_params()

        results = run_batch(tasks, "ExponentialDrift", params, names, max_workers=2)

        assert len(results) == 3
        assert {r.roi for r in results} == {"000", "001", "002"}

        # Results are written back to the h5 file by run_batch itself.
        h5 = H5Data(path)
        for i in range(3):
            roi = f"{i:03d}"
            assert np.isfinite(h5.get_parameter("Z0 (um)", roi))
            dz_dt = h5.get_parameter("dZ_dt (nm_s)", roi)
            assert dz_dt == pytest.approx(-1.0 - 0.1 * i, abs=0.2)

    def test_progress_callback_reaches_total(self, tmp_path):
        path = self._make_file(tmp_path, n_beads=4)
        tasks = [(str(path), f"{i:03d}") for i in range(4)]
        params, names = _drift_params()

        calls = []
        run_batch(tasks, "ExponentialDrift", params, names, max_workers=2,
                   progress_cb=lambda done, total: calls.append((done, total)))

        assert len(calls) == 4
        assert calls[-1] == (4, 4)
        assert all(total == 4 for _, total in calls)

    def test_error_cb_receives_bad_roi_and_others_still_complete(self, tmp_path):
        path = self._make_file(tmp_path, n_beads=2)
        tasks = [(str(path), "000"), (str(path), "999"), (str(path), "001")]
        params, names = _drift_params()

        errors = []
        results = run_batch(tasks, "ExponentialDrift", params, names, max_workers=2,
                             error_cb=errors.append)

        assert len(results) == 2
        assert {r.roi for r in results} == {"000", "001"}
        assert len(errors) == 1
        assert isinstance(errors[0], BatchFitError)
        assert errors[0].roi == "999"

    def test_stop_cb_halts_remaining_work(self, tmp_path):
        path = self._make_file(tmp_path, n_beads=8)
        tasks = [(str(path), f"{i:03d}") for i in range(8)]
        params, names = _drift_params()

        seen = {"done": 0}

        def stop_after_first(done, total):
            seen["done"] = done

        results = run_batch(tasks, "ExponentialDrift", params, names, max_workers=2,
                             progress_cb=lambda d, t: stop_after_first(d, t),
                             stop_cb=lambda: seen["done"] >= 1)

        # Stops promptly: far fewer than all 8 ROIs get fit.
        assert len(results) < 8

    def test_no_leaked_worker_processes_after_completion(self, tmp_path):
        path = self._make_file(tmp_path, n_beads=3)
        tasks = [(str(path), f"{i:03d}") for i in range(3)]
        params, names = _drift_params()

        run_batch(tasks, "ExponentialDrift", params, names, max_workers=2)

        assert multiprocessing.active_children() == []

    def test_no_leaked_worker_processes_after_early_stop(self, tmp_path):
        path = self._make_file(tmp_path, n_beads=8)
        tasks = [(str(path), f"{i:03d}") for i in range(8)]
        params, names = _drift_params()

        run_batch(tasks, "ExponentialDrift", params, names, max_workers=2,
                   stop_cb=lambda: True)

        assert multiprocessing.active_children() == []

    def test_writes_are_batched_per_file_not_per_roi(self, tmp_path, monkeypatch):
        """Regression guard: run_batch must write results via one
        set_parameters_bulk call per file, not one set_parameter call per
        (roi, field). Per-call writes rewrite the entire /Parameters array
        every time, which is O(n_rois^2) overall and was an ~8s bottleneck
        at 300 ROIs in practice -- dwarfing both the fit compute and the
        process-pool overhead this module exists to amortize.
        """
        import trace_editor.core.batch_fit as batch_fit_mod

        path = self._make_file(tmp_path, n_beads=5)
        tasks = [(str(path), f"{i:03d}") for i in range(5)]
        params, names = _drift_params()

        bulk_calls = []
        real_bulk = H5Data.set_parameters_bulk

        def spying_bulk(self, values_by_roi):
            bulk_calls.append(len(values_by_roi))
            return real_bulk(self, values_by_roi)

        monkeypatch.setattr(H5Data, "set_parameters_bulk", spying_bulk)
        monkeypatch.setattr(
            H5Data, "set_parameter",
            lambda *a, **kw: (_ for _ in ()).throw(
                AssertionError("run_batch must not call set_parameter per-ROI")))

        run_batch(tasks, "ExponentialDrift", params, names, max_workers=2)

        # One bulk call for the one file involved, covering all 5 ROIs.
        assert bulk_calls == [5]

    def test_matches_sequential_apply_fit(self, tmp_path):
        """Pooled results match what apply_fit gives sequentially for the same
        data, when the same point selection (_get_fit_selection) is applied
        on both sides -- this is what _fit_one_roi does internally.
        """
        from trace_editor.core.h5_data import apply_fit
        from trace_editor.core.batch_fit import _get_fit_selection

        path = self._make_file(tmp_path, n_beads=1, noise_um=0.05)
        tasks = [(str(path), "000")]
        params, names = _drift_params()

        # Read from a file untouched by run_batch's writes (fit on file A,
        # compare against fitting the same untouched data read via file B).
        h5_before = H5Data(path)
        t = h5_before.get_channel("Time (s)", "000")
        z = h5_before.get_channel("Z (um)", "000")
        sel = _get_fit_selection(h5_before, "000")
        sequential, _ = apply_fit("ExponentialDrift", t, z, selection=sel,
                                   initial_parameters=params)

        [pooled] = run_batch(tasks, "ExponentialDrift", params, names, max_workers=1)

        assert pooled.param_values["Z0 (um)"][0] == pytest.approx(
            sequential.params["Z0"].value, abs=1e-6)
        assert pooled.param_values["dZ_dt (nm_s)"][0] == pytest.approx(
            sequential.params["dZ_dt"].value, abs=1e-6)


class TestDefaultSelectionCriteria:
    """[model_selection.<model>] config-driven default point selection.

    Regression coverage for a real bug: percentile bounds (e.g.
    ExponentialDrift's `min = "p90"` on Shift (mm)) must be inclusive at the
    boundary. A channel that plateaus at its max value (real magnetic-
    tweezer Shift traces often do -- the stepper holds position at the top
    of its range) can have many points exactly equal to its 90th percentile;
    a strict `>` comparison excludes all of them, silently selecting zero
    points instead of "the top 10%".
    """

    def test_percentile_min_is_inclusive_on_plateaued_channel(self, tmp_path):
        from trace_editor.core.batch_fit import _evaluate_criteria

        beads = [make_synthetic_workflow(noise_um=0.0, seed=0)]
        path = write_workflow_h5(tmp_path / "plateau.h5", beads)

        # Force a plateau: clamp the top 20% of Shift (mm) to its max value,
        # mimicking a stepper that holds position at full retraction.
        import h5py
        with h5py.File(path, "r+") as f:
            shift = f["Traces"]["Shift (mm)"][:]
            cutoff = np.percentile(shift, 80)
            shift[shift >= cutoff] = shift.max()
            f["Traces"]["Shift (mm)"][:] = shift

        h5 = H5Data(path)
        criteria = [{"channel": "Shift (mm)", "min": "p90"}]
        sel = _evaluate_criteria(criteria, h5, "000")

        assert sel.sum() > 0, "percentile min should include the plateau, not exclude it"

    def test_exponential_drift_default_selection_nonempty_on_plateaued_shift(self, tmp_path):
        """End-to-end: _get_fit_selection's Default mode for ExponentialDrift
        must not silently select zero points on real-shaped data."""
        from trace_editor.core.batch_fit import _get_fit_selection

        beads = [make_synthetic_workflow(noise_um=0.0, seed=0)]
        path = write_workflow_h5(tmp_path / "plateau2.h5", beads)

        import h5py
        with h5py.File(path, "r+") as f:
            shift = f["Traces"]["Shift (mm)"][:]
            cutoff = np.percentile(shift, 80)
            shift[shift >= cutoff] = shift.max()
            f["Traces"]["Shift (mm)"][:] = shift

        criteria = [{"channel": "Shift (mm)", "min": "p90"}]

        h5 = H5Data(path)
        sel = _get_fit_selection(h5, "000", mode="default", criteria=criteria)

        assert sel is not None
        assert sel.sum() > 0


class TestModelWantsMultiprocessing:
    """The :multiprocessing: docstring tag that opts a model into pooling.

    Benchmarking showed the pool is a net loss for cheap models (import +
    IPC overhead exceeds the fit cost) and a large win for expensive ones
    (zHMf/zHMf2, where a single fit is on the order of a second or more) --
    so pooling must be opt-in per model, read from the model's own
    docstring rather than hardcoded here.
    """

    def test_heavy_models_opt_in(self):
        assert model_wants_multiprocessing("zHMf") is True
        assert model_wants_multiprocessing("zHMf2") is True

    def test_cheap_models_do_not_opt_in(self):
        assert model_wants_multiprocessing("ExponentialDrift") is False
        assert model_wants_multiprocessing("zWLC") is False
        assert model_wants_multiprocessing("zFiber") is False

    def test_unknown_model_does_not_opt_in(self):
        assert model_wants_multiprocessing("NoSuchModel") is False

    def test_tag_does_not_corrupt_prepare_fit_param_ranges(self):
        """Regression guard: prepare_fit's own docstring parser finds
        parameter ranges via `doc.split(f"{key}:")[-1]` (last occurrence,
        no exact-line matching), so a naively-placed :multiprocessing: line
        could in principle be mistaken for part of a :param <name>: line.
        Confirm zHMf/zHMf2 -- the two tagged models -- still parse every
        parameter's [min, max] range correctly after tagging.
        """
        for model_name in ("zHMf", "zHMf2"):
            _, params, names = prepare_fit(model_name)
            assert names["x"] == "F (pN)"
            assert names["y"] == "Z (um)"
            assert params["L"].min == 1.0
            assert params["L"].max == float("inf")
            assert params["Z0"].min == -2.0
            assert params["Z0"].max == 0.0
            assert params["angle"].min == 0.0
            assert params["angle"].max == 90.0

    def test_run_batch_auto_detect_uses_pool_for_tagged_model(self, tmp_path, monkeypatch):
        """run_batch must actually spawn a ProcessPoolExecutor when the model
        opts in, and must NOT when it doesn't -- verified by spying on
        ProcessPoolExecutor.__init__ rather than trusting the return value,
        since both paths produce equivalent RoiFitResults.
        """
        from concurrent.futures import ProcessPoolExecutor

        pool_created = []
        real_init = ProcessPoolExecutor.__init__

        def spying_init(self, *a, **kw):
            pool_created.append(True)
            return real_init(self, *a, **kw)

        monkeypatch.setattr(ProcessPoolExecutor, "__init__", spying_init)

        beads = [make_synthetic_workflow(noise_um=0.0, seed=0)]
        path = write_workflow_h5(tmp_path / "one.h5", beads)
        tasks = [(str(path), "000")]
        params, names = _drift_params()

        run_batch(tasks, "ExponentialDrift", params, names)
        assert pool_created == [], "cheap model must not spawn a pool"

        run_batch(tasks, "ExponentialDrift", params, names, use_pool=True)
        assert pool_created == [True], "use_pool=True must force a pool regardless of the model"

    def test_run_batch_use_pool_false_forces_sequential(self, tmp_path, monkeypatch):
        import trace_editor.core.batch_fit as batch_fit_mod
        from concurrent.futures import ProcessPoolExecutor

        pool_created = []
        real_init = ProcessPoolExecutor.__init__
        monkeypatch.setattr(
            ProcessPoolExecutor, "__init__",
            lambda self, *a, **kw: pool_created.append(True) or real_init(self, *a, **kw))

        # Force "opts in" even for a cheap model, to prove use_pool=False wins.
        monkeypatch.setattr(batch_fit_mod, "model_wants_multiprocessing", lambda name: True)

        beads = [make_synthetic_workflow(noise_um=0.0, seed=0)]
        path = write_workflow_h5(tmp_path / "one.h5", beads)
        tasks = [(str(path), "000")]
        params, names = _drift_params()

        run_batch(tasks, "ExponentialDrift", params, names, use_pool=False)
        assert pool_created == [], "use_pool=False must force sequential regardless of the model"
