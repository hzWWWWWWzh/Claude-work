# TraceEditor

A PySide6 desktop app for reviewing and fitting single-molecule force
spectroscopy traces from a magnetic tweezers instrument. It reads the new
instrument `.h5` format and fits exponential drift and worm-like-chain (zWLC)
models to per-bead extension data.

It replaces an older tkinter app (kept for reference in `oldhdf_versions/`).

## Repo layout

- `src/trace_editor/` — the current app (PySide6 GUI, `.h5` format)
- `oldhdf_versions/` — legacy scripts for the old PyTables/`.hdf` format (kept for reference, not used by new work)
- `notebooks/`, `nb_scripts/` — analysis notebooks (paired via jupytext)
- `docs/` — architecture, file format, and generated API reference
- `tests/` — test suite (see [`tests/README.md`](tests/README.md))

## Documentation

- [`docs/architecture.md`](docs/architecture.md) — internals: UI layout, Z-correction logic, fit/batch/selection flow
- [`docs/file-format.md`](docs/file-format.md) — the `.h5` (and legacy `.hdf`) file format spec
- [`docs/api/`](docs/api/) — API reference generated from docstrings
- [`CONTRIBUTING.md`](CONTRIBUTING.md) — design philosophy and contributor rules

## Installation

### Install with uv (fastest):

Best installed with [uv](https://docs.astral.sh/uv/), a 10–100× faster alternative to pip/conda.

**Install uv:**

```bash
# macOS / Linux
curl -LsSf https://astral.sh/uv/install.sh | sh
```
```powershell
# Windows
powershell -ExecutionPolicy ByPass -c "irm https://astral.sh/uv/install.ps1 | iex"
```
See the [uv installation guide](https://docs.astral.sh/uv/getting-started/installation/) for more details.

Download the repository then run:

```bash
uv sync
```

This builds a virtual environment at `.venv` in the project folder with all base dependencies pinned to `uv.lock` for reproducibility.

**Optional extras:**

| Extra | Who it's for | Command |
|---|---|---|
| `oldhdf` | Anyone running the legacy scripts in `oldhdf_versions/` (**File > Import Legacy .hdf…** works without this — PyTables is a base dependency) | `uv sync --extra oldhdf` |

**To open trace_editor**
Simplest is to use `uv run trace_editor` which automatically selects the virtual environment.
Alternatively you can activate the virtual environment (like conda) with:
| Platform | Command |
|----------|---------|
| macOS/Linux | `source .venv/bin/activate` |
| Windows (cmd) | `.venv\Scripts\activate.bat` |
| Windows (PowerShell) | `.venv\Scripts\Activate.ps1` |
and then run with `trace_editor` or `python -m trace_editor`.

---

### Install with pip:

```bash
# base only
pip install -e .
# with legacy oldhdf support
pip install -e '.[oldhdf]'
```

### Install in conda env:

```bash
conda create -n trace_editor python=3.13
conda activate trace_editor
pip install -e '.[oldhdf]' # or just '.'
```
and then run with `trace_editor` or `python -m trace_editor`.

---

## Usage

### Run the app

```bash
uv run trace_editor
# or
uv run python -m trace_editor
```

### Typical workflow

1. Open `.h5` files (LUT files are skipped automatically).
2. Browse ROIs with Shift+F11 / Shift+F12; toggle selection with Shift+F10.
3. **Calculate Force** — converts Shift (mm) → F (pN), stored as a shared channel.
4. Fit **ExponentialDrift** (batch) — sets Z0, dZ_dt, tau per ROI; also X0/Y0.
5. Fit **zWLC** (batch) — fits DNA extension; stores L (bp), Zbead (um), R2 zWLC.
6. Select → Traces → Criteria: keep traces with `R2 zWLC` > 0.9.
7. Browse, inspect, export.

See [`docs/architecture.md`](docs/architecture.md) for keyboard shortcuts and the
full internals.

### Auto-opening + auto-fitting new files

If your tweezer acquisition software writes `.h5` files to a folder on the
same computer, TraceEditor can watch that folder and automatically fit +
open each new file as it's written — no external script needed:

1. **File > Watch Folder…** — pick the folder your acquisition software
   writes to. **File > Watching** toggles it on/off without forgetting the
   folder.
2. Pick a **pipeline** from the toolbar dropdown (below the model selector) —
   a named, ordered sequence of fit steps defined in `app_config.toml` (see
   "Fit pipelines" below).
3. **Parameters > Auto-fit on File Open** — when checked, every newly opened
   file (watched, manually opened, or legacy-imported) is fit automatically
   with the selected pipeline.

The toolbar's **Run** button runs the selected pipeline against all open
files on demand, independent of auto-fit-on-open — handy after editing the
pipeline config, or for files opened before auto-fit was turned on.

### Fit pipelines

A pipeline is a named, ordered sequence of fit steps — e.g. fit
`ExponentialDrift`, calculate Force, then fit `zWLC` — declared entirely in
`app_config.toml`. This is the single source of truth for both the
interactive Fit dialog and automatic fitting; see
`app_config.default.toml` for the full schema and a worked example, and
[`docs/architecture.md`](docs/architecture.md) ("Fit pipelines") for
internals.

### Inspecting `.h5` files

**File > Inspect H5…** opens the current file's raw HDF5 structure (all
groups/datasets/attributes, not just what TraceEditor's UI shows) in the
[H5Web](https://github.com/silx-kit/h5web) viewer, via the
[H5Web VS Code extension](https://marketplace.visualstudio.com/items?itemName=h5web.vscode-h5web).

Setup (one-time):

1. Install [VS Code](https://code.visualstudio.com/) if you don't have it.
2. Install the **H5Web** extension from the VS Code Extensions panel
   (search "H5Web", publisher `h5web`), or from the
   [marketplace page](https://marketplace.visualstudio.com/items?itemName=h5web.vscode-h5web).
3. Make sure the `code` CLI is on your `PATH` — in VS Code, open the Command
   Palette (Cmd/Ctrl+Shift+P) and run **Shell Command: Install 'code' command
   in PATH**.

If VS Code isn't set up, TraceEditor falls back to opening
[myhdf5.hdfgroup.org](https://myhdf5.hdfgroup.org) in your browser — a
zero-install HDF5 viewer where you manually pick the file (files are read
entirely client-side; nothing is uploaded).

### Per-user defaults (force calibration + fit pipelines)

Force-calibration constants (`L1`, `L2`, `Fmin`, `Fmax`, `alpha`, `offset` —
used as the fallback for `calculate_force` pipeline steps that don't set a
given key themselves) and fit pipelines (named, ordered sequences of fit
steps — each pipeline fully explicit, with no fallback between pipelines)
live in
[`src/trace_editor/gui/app_config.default.toml`](src/trace_editor/gui/app_config.default.toml).
These defaults suit a typical tweezer setup, but once you've calibrated your
own instrument (and know your usual DNA construct length, etc.) you generally
want to keep using those values without retyping them.

To set your own defaults: copy `app_config.default.toml` to `app_config.toml`
in the same directory (`src/trace_editor/gui/`) and edit your copy.
`app_config.toml` is gitignored, so it stays local to you and takes priority
over the checked-in defaults whenever it exists.
```bash
cp src/trace_editor/gui/app_config.default.toml src/trace_editor/gui/app_config.toml
```

### Point selection (Default vs. Manual)

Each pipeline step can define a point-selection criterion (`criteria` — see
"Fit pipelines" above), evaluated automatically when fitting — e.g.
`ExponentialDrift` restricts to the lowest-force points (`Shift (mm)` ≥ 90th
percentile) in the checked-in `default_pipeline`. The Fit dialog has a
**Default / Manual** toggle: Manual mode uses whatever point mask is currently persisted for that
ROI (**Select > Points > By Criteria…** or **Select > Points > In view**,
which builds criteria from the plot's current axis limits). Switching modes
just changes which mask is highlighted/fit — nothing is overwritten until
you actually run a fit. See `docs/architecture.md` ("Selection modes") for
details.

**Known limitation:** `ExponentialDrift`'s default criterion currently only
restricts by `Shift (mm)` percentile. The legacy app also required
sufficient bead tracking amplitude in the same fit — points from the first
tenth of the trace set a reference amplitude, and later points below ~1/10th
of that reference were excluded, catching tracking dropout/bead loss that a
pure force-based filter wouldn't. That amplitude sub-filter was not ported
to the new `[model_selection]` config and should be revisited if drift fits
regress on real (noisy/dropout-prone) data — see
`oldhdf_versions/TraceIO.py` `fit()` (lines 1005-1014) for the original logic.
Note the same file also special-cases `zMagnetForce` similarly (lines
1017-1023, high-extension + amplitude filter), which was never ported to
either the old or new TraceEditor's `_get_fit_selection`.

