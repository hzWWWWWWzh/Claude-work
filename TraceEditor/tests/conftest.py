"""Pytest configuration and fixtures for trace_editor tests.

The synthetic data here is built by driving the *real* ForceSpectroscopy
models (`MagnetForce`, `zWLC`, `ExponentialDrift`) so that tests exercise the
same maths as the app. See `make_synthetic_workflow` for the ground-truth
construction of the standard analysis workflow.
"""

import os
import sys
import pytest
import numpy as np
import h5py
import tempfile
from pathlib import Path

# Render Qt offscreen so GUI tests never pop up a real window (CI-safe). Must be
# set before QApplication is created, i.e. before importing QtWidgets below.
os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

from PySide6.QtWidgets import QApplication

# Make sibling helper modules (e.g. workflow_helpers) importable from tests.
sys.path.insert(0, str(Path(__file__).parent))

from trace_editor.core import force_spectroscopy as fs


# Directory where inspectable artifacts are written (git-ignored except real_data.h5).
DATA_DIR = Path(__file__).parent / "data"

# Ground-truth parameters for the synthetic workflow dataset. These are the
# values the full open→force→drift→WLC pipeline must recover.
TRUE_L_BP = 10000.0       # contour length (bp)
TRUE_A_NM = 50.0          # persistence length (nm) — zWLC default
TRUE_S_PN = 1000.0        # stretch modulus (pN) — zWLC default
TRUE_ZBEAD_UM = 0.0       # bead height offset (um)
FMAX_PN = 100.0           # MagnetForce fmax used by the app's calculate_force


@pytest.fixture(scope="session")
def qapp():
    """Create a QApplication instance for the entire test session."""
    app = QApplication.instance()
    if app is None:
        app = QApplication([])
    yield app


@pytest.fixture(autouse=True)
def isolated_qsettings(tmp_path_factory, monkeypatch):
    """Redirect every `QSettings("TraceEditor", "TraceEditor")` construction
    to a throwaway .ini file for the duration of each test, so tests never
    read or write the developer's real app preferences (last_dir,
    watch_folder, and notably active_pipeline -- which, unlike the others,
    actually changes fit behavior, e.g. which selection criteria apply --
    see TestPointSelection/TestFitThroughGui, which broke the first time
    this leaked).

    `QSettings.setPath`/`setDefaultFormat` don't redirect the plain
    `QSettings(org, app)` two-argument constructor on macOS (it resolves to
    NativeFormat/the real plist regardless), so this patches the class
    itself to always use the explicit-format constructor instead.
    """
    from PySide6.QtCore import QSettings

    ini_dir = str(tmp_path_factory.mktemp("qsettings"))
    real_init = QSettings.__init__

    def patched_init(self, *args, **kwargs):
        if args == ("TraceEditor", "TraceEditor"):
            real_init(self, QSettings.Format.IniFormat, QSettings.Scope.UserScope,
                       "TraceEditor", "TraceEditor")
        else:
            real_init(self, *args, **kwargs)

    monkeypatch.setattr(QSettings, "__init__", patched_init)
    QSettings.setPath(QSettings.Format.IniFormat, QSettings.Scope.UserScope, ini_dir)
    yield


# ---------------------------------------------------------------------------
# Synthetic workflow data (the core of the fitting/workflow tests)
# ---------------------------------------------------------------------------

def make_shift_sweep(n_samples: int, n_cycles: int = 2,
                     shift_min: float = 1.0, shift_max: float = 10.8,
                     dwell_frac: float = 0.12) -> np.ndarray:
    """Shift (mm) sweep with dwell times at top and bottom, mirroring the
    real instrument trajectory.

    High Shift = magnet retracted = LOW force; low Shift = HIGH force.

    The trace starts at the top, then for each cycle: ramps down, dwells at the
    bottom, ramps up, dwells at the top — ending at shift_max. `dwell_frac` is
    the fraction of each ramp+dwell segment spent dwelling at the extreme.

    With n_cycles=2 the trace goes down-and-back twice, ending at the top,
    matching real_data.h5.
    """
    # Segments per cycle: [dwell top, ramp down, dwell bottom, ramp up].
    # We prepend an initial top dwell and let each cycle end back at the top.
    n_ramp = max(2, int(round(n_samples * (1.0 - dwell_frac) / (2 * n_cycles))))
    n_dwell = max(1, int(round(n_samples * dwell_frac / (2 * n_cycles))))

    down = np.linspace(shift_max, shift_min, n_ramp)
    up = np.linspace(shift_min, shift_max, n_ramp)
    dwell_top = np.full(n_dwell, shift_max)
    dwell_bot = np.full(n_dwell, shift_min)

    segments = [dwell_top]  # start with a top dwell
    for _ in range(n_cycles):
        segments += [down, dwell_bot, up, dwell_top]

    sweep = np.concatenate(segments)

    # Resample to exactly n_samples (preserves dwell/ramp shape).
    xp = np.linspace(0.0, 1.0, len(sweep))
    x = np.linspace(0.0, 1.0, n_samples)
    return np.interp(x, xp, sweep)


def make_synthetic_workflow(n_samples: int = 1500, n_cycles: int = 2,
                            noise_um: float = 0.0,
                            Z0: float = 12.0, dZ_dt: float = -1.5,
                            tau: float = np.inf, seed: int = 0) -> dict:
    """Build one bead's ground-truth workflow data using the real fs models.

    Pipeline (matches the app):
      1. Shift (mm) triangle sweep.
      2. F (pN)   = MagnetForce(shift, fmax=FMAX_PN).
      3. z_wlc    = zWLC(F, L=TRUE_L_BP, z_bead=TRUE_ZBEAD_UM).   ground truth
      4. drift    = ExponentialDrift(time, Z0, dZ_dt, tau).
      5. Z (um)   = z_wlc + drift + optional noise.               observed

    Returns a dict with time, shift, force, z_wlc (truth), drift, z (observed)
    and the true parameters.
    """
    rng = np.random.default_rng(seed)
    time = np.linspace(0.0, 36.0, n_samples)
    shift = make_shift_sweep(n_samples, n_cycles=n_cycles)

    force = np.asarray(fs.MagnetForce(shift, fmax=FMAX_PN))
    z_wlc = np.asarray(fs.zWLC(force, L=TRUE_L_BP, A=TRUE_A_NM,
                               S=TRUE_S_PN, z_bead=TRUE_ZBEAD_UM))
    drift = np.asarray(fs.ExponentialDrift(time, Z0=Z0, dZ_dt=dZ_dt, tau=tau))

    z = z_wlc + drift
    if noise_um > 0.0:
        z = z + rng.normal(0.0, noise_um, size=n_samples)

    # A high, roughly constant bead amplitude so the drift-fit amplitude filter
    # (A >= 95th-pct-of-first-tenth / 10) keeps the low-force points.
    amp = 100.0 + rng.normal(0.0, 1.0, size=n_samples)

    return {
        "time": time,
        "shift": shift,
        "force": force,
        "z_wlc": z_wlc.astype(np.float64),
        "drift": drift,
        "z": z.astype(np.float64),
        "amp": amp,
        "Z0": Z0, "dZ_dt": dZ_dt, "tau": tau,
        "L_bp": TRUE_L_BP, "A_nm": TRUE_A_NM, "S_pN": TRUE_S_PN,
        "z_bead": TRUE_ZBEAD_UM,
    }


def write_workflow_h5(filepath: Path, beads: list, notes: str = "synthetic") -> Path:
    """Write a list of per-bead workflow dicts to a new-format .h5 file.

    Shared Time/Shift come from the first bead (identical across beads).
    Format matches real_data.h5: root attrs, /Parameters compound dataset,
    /Traces with shared datasets and per-channel groups keyed by "NNN".
    """
    filepath = Path(filepath)
    filepath.parent.mkdir(parents=True, exist_ok=True)
    n_rois = len(beads)
    n_samples = len(beads[0]["time"])

    with h5py.File(filepath, "w") as f:
        f.attrs["Date"] = "2026-07-07 10:00:00"
        f.attrs["User"] = "test"
        f.attrs["Notes"] = notes
        f.attrs["axis"] = "Shift_mm"
        f.attrs["num_rois"] = n_rois
        f.attrs["roi_size"] = 64
        f.attrs["bin_roi_size"] = 64

        param_dtype = np.dtype([
            ("ROI", np.int32),
            ("X (pix)", np.float32),
            ("Y (pix)", np.float32),
            (" selected", np.float64),
        ])
        params = np.zeros(n_rois, dtype=param_dtype)
        for i in range(n_rois):
            params[i]["ROI"] = i
            params[i]["X (pix)"] = 640.0
            params[i]["Y (pix)"] = 640.0
            params[i][" selected"] = 1.0
        f.create_dataset("Parameters", data=params)

        traces = f.create_group("Traces")
        traces.create_dataset("Time (s)", data=beads[0]["time"], dtype=np.float64)
        traces.create_dataset("Shift (mm)", data=beads[0]["shift"], dtype=np.float64)
        traces.create_dataset("Intensity",
                              data=np.full(n_samples, 100.0), dtype=np.float64)

        z_grp = traces.create_group("Z (um)")
        x_grp = traces.create_group("X (um)")
        y_grp = traces.create_group("Y (um)")
        a_grp = traces.create_group("A (a.u.)")
        for i, bead in enumerate(beads):
            key = f"{i:03d}"
            z_grp.create_dataset(key, data=bead["z"], dtype=np.float32)
            x_grp.create_dataset(key, data=np.full(n_samples, 640.0), dtype=np.float32)
            y_grp.create_dataset(key, data=np.full(n_samples, 640.0), dtype=np.float32)
            a_grp.create_dataset(key, data=bead["amp"], dtype=np.float32)

    return filepath


@pytest.fixture
def synthetic_bead_noiseless():
    """One bead of ground-truth workflow data, no noise."""
    return make_synthetic_workflow(noise_um=0.0, seed=0)


@pytest.fixture
def synthetic_bead_noisy():
    """One bead of ground-truth workflow data with small measurement noise."""
    return make_synthetic_workflow(noise_um=0.05, seed=1)


# Standard 2-bead synthetic dataset (bead 0 noiseless, bead 1 small noise).
# Used to build all three on-disk workflow files below.
def _standard_beads():
    return [
        make_synthetic_workflow(noise_um=0.0, seed=0, Z0=12.0, dZ_dt=-1.5),
        make_synthetic_workflow(noise_um=0.05, seed=1, Z0=11.0, dZ_dt=-1.0),
    ]


# Three on-disk workflow files with distinct lifecycles:
#   1. committed, permanent, raw (no fits)      -> synthetic_workflow.h5
#   2. regenerated every run, raw (no fits)     -> synthetic_workflow_generated.h5
#   3. regenerated every run, raw + workflow fits -> synthetic_workflow_fitted.h5
# Only file 1 is committed (see tests/data/.gitignore).
PERMANENT_WORKFLOW_H5 = DATA_DIR / "synthetic_workflow.h5"
GENERATED_WORKFLOW_H5 = DATA_DIR / "synthetic_workflow_generated.h5"
FITTED_WORKFLOW_H5 = DATA_DIR / "synthetic_workflow_fitted.h5"

_WORKFLOW_NOTES = "synthetic workflow: bead 0 noiseless, bead 1 small noise; no fits"


@pytest.fixture(scope="session")
def permanent_workflow_h5():
    """The committed, permanent raw synthetic file (no fits).

    Created only if missing so it stays stable in git; delete it and re-run to
    regenerate. Read-only for tests — never write fits into this one.
    """
    if not PERMANENT_WORKFLOW_H5.exists():
        write_workflow_h5(PERMANENT_WORKFLOW_H5, _standard_beads(), notes=_WORKFLOW_NOTES)
    return PERMANENT_WORKFLOW_H5


@pytest.fixture
def generated_workflow_h5():
    """A freshly regenerated raw synthetic file (no fits, git-ignored).

    Rewritten on every request so it always reflects the current generator.
    """
    write_workflow_h5(GENERATED_WORKFLOW_H5, _standard_beads(), notes=_WORKFLOW_NOTES)
    return GENERATED_WORKFLOW_H5


@pytest.fixture
def fitted_workflow_h5(generated_workflow_h5):
    """Copy of the generated file for the workflow test to write fits into.

    Git-ignored; left on disk after the run so it can be opened in the app to
    visually inspect the drift and zWLC fits.
    """
    import shutil
    shutil.copyfile(generated_workflow_h5, FITTED_WORKFLOW_H5)
    return FITTED_WORKFLOW_H5


# ---------------------------------------------------------------------------
# Generic small file for file-I/O tests
# ---------------------------------------------------------------------------

@pytest.fixture
def temp_h5_file():
    """Create a temporary H5 file with realistic test data (5 ROIs)."""
    with tempfile.NamedTemporaryFile(suffix=".h5", delete=False) as f:
        filepath = Path(f.name)

    n_rois = 5
    n_samples = 1000
    rng = np.random.default_rng(42)

    with h5py.File(filepath, "w") as f:
        f.attrs["Date"] = "2026-07-07 10:00:00"
        f.attrs["User"] = "test_user"
        f.attrs["Notes"] = "Test experiment data"
        f.attrs["axis"] = "Shift_mm"
        f.attrs["num_rois"] = n_rois
        f.attrs["roi_size"] = 64
        f.attrs["bin_roi_size"] = 64

        param_dtype = np.dtype([
            ("ROI", np.int32),
            ("X (pix)", np.float32),
            ("Y (pix)", np.float32),
            (" selected", np.float64),
        ])
        params = np.zeros(n_rois, dtype=param_dtype)
        for i in range(n_rois):
            params[i]["ROI"] = i
            params[i]["X (pix)"] = 32 + rng.random() * 10
            params[i]["Y (pix)"] = 32 + rng.random() * 10
            params[i][" selected"] = 1.0
        f.create_dataset("Parameters", data=params)

        settings_dtype = np.dtype([
            ("Group", "S50"), ("Key", "S50"), ("Value", "S100"),
        ])
        settings = np.array(
            [(b"Camera", b"Gain", b"100"), (b"Laser", b"Power", b"50")],
            dtype=settings_dtype
        )
        f.create_dataset("Settings", data=settings)

        traces = f.create_group("Traces")
        time = np.linspace(0, 100, n_samples)
        shift = 5 + 0.01 * time + 0.1 * np.sin(2*np.pi*time/50)
        intensity = 100 + 10*rng.standard_normal(n_samples)
        traces.create_dataset("Time (s)", data=time, dtype=np.float64)
        traces.create_dataset("Shift (mm)", data=shift, dtype=np.float64)
        traces.create_dataset("Intensity", data=intensity, dtype=np.float64)

        channels_per_roi = ["Z (um)", "X (um)", "Y (um)", "A (a.u.)"]
        for channel in channels_per_roi:
            channel_group = traces.create_group(channel)
            for i in range(n_rois):
                z0 = 10 + i
                drift_rate = 0.001 * (i + 1)
                if channel == "Z (um)":
                    data = z0 + drift_rate * time + 0.5 * np.sin(2*np.pi*time/50) + 0.05*rng.standard_normal(n_samples)
                elif channel == "X (um)":
                    data = 32 + i + 0.1*rng.standard_normal(n_samples)
                elif channel == "Y (um)":
                    data = 32 + 0.1*rng.standard_normal(n_samples)
                else:  # A (a.u.)
                    data = 50 + 5*rng.standard_normal(n_samples)
                channel_group.create_dataset(f"{i:03d}", data=data, dtype=np.float32)

    yield filepath
    filepath.unlink()


@pytest.fixture
def real_h5_file():
    """Path to the committed real-data reference file, or None if absent.

    Note: real_data.h5 is low quality (first ok-ish bead is 008); use it for
    format/smoke checks, not for asserting fit accuracy.
    """
    test_file = DATA_DIR / "real_data.h5"
    return test_file if test_file.exists() else None


# ---------------------------------------------------------------------------
# GUI helpers
# ---------------------------------------------------------------------------

def load_file_headless(window, path):
    """Load an .h5 into a TraceEditor without going through the file dialog.

    Mirrors the tail of `open_files`: register the H5Data, add it to the file
    combo, and activate it (which populates channels, the ROI list, and draws).
    """
    from trace_editor.core.h5_data import H5Data
    from pathlib import Path as _Path

    path = str(path)
    data = H5Data(path)
    window._open_files[path] = data
    window.fileComboBox.blockSignals(True)
    window.fileComboBox.addItem(_Path(path).name, userData=path)
    window.fileComboBox.blockSignals(False)
    window._activate_file(path)
    return data


@pytest.fixture
def editor(qapp):
    """A bare TraceEditor window (offscreen), closed on teardown."""
    from trace_editor.gui.main_window import TraceEditor
    win = TraceEditor()
    yield win
    win.close()


@pytest.fixture
def editor_with_file(qapp, generated_workflow_h5):
    """A TraceEditor with the synthetic workflow file loaded headlessly.

    Yields (window, h5data). Uses the git-ignored generated file so tests may
    freely write fits into it without touching committed data.
    """
    from trace_editor.gui.main_window import TraceEditor
    win = TraceEditor()
    data = load_file_headless(win, generated_workflow_h5)
    yield win, data
    win.close()
