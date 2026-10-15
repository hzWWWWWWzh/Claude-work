# TraceEditor Test Suite

Tests for file I/O, the fitting models, the full analysis workflow, and the GUI.

## Install & run

```bash
uv sync --extra test          # or: pip install -e '.[test]'
uv run pytest                 # all tests
uv run pytest -m workflow     # just the workflow integration tests
uv run pytest --cov=src/trace_editor --cov-report=html
```

Config lives in `pyproject.toml` under `[tool.pytest.ini_options]` (there is no
`pytest.ini`). Markers: `io`, `fitting`, `workflow`, `gui`, `slow`.

## Layout

| File | What it covers |
|------|----------------|
| `conftest.py` | Fixtures + synthetic-data builders (drives the real `fs` models) |
| `workflow_helpers.py` | The two-stage drift→WLC fit, wrapping the app's `apply_fit` |
| `test_file_io.py` | `H5Data` read/write, channels, parameters, disk round-trip; on the temp, synthetic, and real files |
| `test_fitting.py` | Real `ExponentialDrift` / `MagnetForce` / `zWLC` models + `apply_fit`/`prepare_fit` |
| `test_workflow.py` | Full open→force→drift→WLC pipeline on synthetic data |
| `test_gui.py` | Real open→select ROI→plot→fit flow, driven headlessly (pytest-qt) |
| `test_batch_fit.py` | `core.batch_fit`: sequential/pooled correctness, stop/error handling, no leaked worker processes, and the `:multiprocessing:` docstring opt-in (see `docs/architecture.md`) — all Qt-free |

## The workflow test (`test_workflow.py`)

This reproduces the exact analysis done in the app and checks the logic end to
end. Ground truth is built with the **real** ForceSpectroscopy models:

```
Shift sweep (dwell at top → down → dwell at bottom → up, ×2, ending at top)
  → MagnetForce(fmax=100)         = F (pN)
  → zWLC(L=10000, z_bead=0)       = z_wlc          (ground truth extension)
  → + ExponentialDrift(Z0,dZ_dt)  = Z (um)         (observed trace)
```

The workflow under test then inverts it:

1. **Fit `ExponentialDrift`** over the low-force points only (high Shift ≥ 90th
   percentile — mirroring the GUI's `_get_fit_selection`) and subtract it. This
   also removes the non-zero WLC height at low force.
2. **Fit `zWLC`** with free `L` and `z_bead` (A, S fixed). This restores the WLC
   shape and height.

Assertions (noiseless: tight; small-noise: loose):

- `L` returns to 10000 bp.
- WLC **shape** restored — offset-removed RMS is small.
- WLC **height** (peak-to-peak, offset-independent) restored.
- **`z_bead` invariant**: `z_bead ≈ -mean(zWLC(F_selected))` — the drift fit's
  `Z0` absorbs the mean WLC height over the selected low-force points, and
  `z_bead` cancels exactly that. (A constant vertical offset is genuinely
  degenerate between drift `Z0` and `z_bead`, so we do *not* assert absolute
  `z_bead`.)

## Test data files (`tests/data/`)

Three synthetic files with distinct lifecycles, plus the real reference:

| File | Committed? | Contents |
|------|-----------|----------|
| `real_data.h5` | yes | Real instrument data (low quality; first ok-ish bead is 008). Format reference / smoke tests only — do **not** assert fit accuracy on it. |
| `synthetic_workflow.h5` | yes | Permanent raw synthetic input (2 beads, no fits). Stable for inspection; delete + re-run to regenerate. |
| `synthetic_workflow_generated.h5` | no (git-ignored) | Raw synthetic, regenerated every run. Guards that generation is reproducible. |
| `synthetic_workflow_fitted.h5` | no (git-ignored) | The generated file + workflow fits written into `/Parameters` and `F (pN)` into Traces. **Open this in the app to visually inspect the fits.** |

`bead 0` is noiseless, `bead 1` has small noise. All synthetic data uses fixed
seeds, so the committed and regenerated files match byte-for-byte in their
traces.

## Fixtures (`conftest.py`)

- `qapp` — session QApplication for GUI tests.
- `temp_h5_file` — generic 5-ROI temp file for file-I/O tests.
- `real_h5_file` — path to `real_data.h5`, or `None`.
- `synthetic_bead_noiseless` / `synthetic_bead_noisy` — one bead of ground-truth
  workflow data (dict) with the true parameters.
- `permanent_workflow_h5` / `generated_workflow_h5` / `fitted_workflow_h5` — the
  three on-disk files above.

`make_synthetic_workflow(...)` and `write_workflow_h5(...)` are the underlying
builders if you need custom parameters in a new test.

## GUI tests

Qt runs **offscreen** (`QT_QPA_PLATFORM=offscreen`, set in `conftest.py` before
the QApplication is created), so no window ever appears and the tests are
CI-safe — nothing pops up on your screen.

Files are loaded **headlessly** via `conftest.load_file_headless(window, path)`,
which bypasses the native file dialog and mirrors the tail of `open_files`. This
lets the tests drive the real flow and assert on outcomes rather than just
constructing the window:

- Loading a file populates the file / ROI / channel combos and picks a current ROI.
- `redraw()` scatters the selected trace onto the axes (points land in
  `ax.collections`; `ax.lines` only holds a model curve when a fit is shown).
- `calculate_force()` writes the `F (pN)` shared dataset and flags each ROI.
- `_run_fit(...)` runs the real drift→WLC fit and writes recovered params
  (L ≈ 10000 bp, R²) into the h5.

Fixtures: `editor` (bare window) and `editor_with_file` (window + synthetic file
loaded, yields `(window, h5data)`).

## Batch fit benchmark (`scripts/benchmark_batch_fit.py`)

Not part of the pytest suite (no assertions, just timing output) — a
reproducible way to check whether batch-fit multiprocessing is actually
paying off, and a regression guard for the sequential-vs-pooled decision
itself (see `docs/architecture.md`'s "Sequential vs. multiprocessing"
section for the full explanation and the rule-of-thumb threshold for
tagging a model `:multiprocessing: true`):

```bash
uv run python scripts/benchmark_batch_fit.py zhmf2 51   # heavy model: pool should win
uv run python scripts/benchmark_batch_fit.py drift 150  # cheap model: pool should lose
```
