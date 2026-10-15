# Contributing to TraceEditor

## Design philosophy (read before making changes)

- **Keep it small and simple.** Push back on requests that add complexity,
  abstraction layers, or features not clearly needed. Three simple lines beat a
  clever abstraction.
- **Match the original.** UI layout, fitting workflow, Z correction logic,
  parameter naming — all should match the legacy `TraceEditor.py` and its
  siblings unless there is a concrete reason to differ. When in doubt, read the
  original first.
- **Comments explain _why_, not _what_** — and only when non-obvious. Don't add
  docstrings to obvious methods; do document the public API and the non-trivial
  GUI entry points (see "Docstrings" below).
- **UI changes go in the `.ui` file only.** Edit `main_window.ui`, then
  regenerate `ui_main_window.py` with `pyside6-uic`. Never edit
  `ui_main_window.py` by hand.
- **Keep documentation in sync.** Any change that affects behaviour, the file
  format, the workflow, or the UI must be reflected in the relevant doc
  ([`docs/architecture.md`](docs/architecture.md),
  [`docs/file-format.md`](docs/file-format.md), or a docstring) before the task
  is considered done.

## Legacy codebase (reference — do not modify)

The old tkinter app is the canonical reference for layout, workflow, and fitting
logic. It lives in `oldhdf_versions/` at the repo root and must **not** be
modified. When implementing any feature, **read the corresponding original code
first** before writing anything.

| File | Role |
|---|---|
| `TraceEditor.py` | Main tkinter app — layout, draw_figure, fit workflow, keyboard nav |
| `TraceEditorGUIs.py` | tkinter dialogs: SetFitParametersWindow, SetSelectPointsWindow, SetSelectTracesWindow, HistogramGUI |
| `TraceIO.py` | Old HDF backend: `hdf_data` class, `fit()`, `select_points()`, `select_trace()`, `calculate_drift()` |
| `ForceSpectroscopy.py` | Physics models — **also used by the new app** (copied to `core/force_spectroscopy.py`) |
| `MultipProcessJob.py` | Multiprocessing batch runner (replaced by QThread in the new app) |
| `FitSteps.py` | Step-fitting utilities |

## Development setup

```bash
uv sync --extra test          # install with test deps
uv run trace_editor           # run the app
uv run pytest                 # run the test suite
```

## UI changes — mandatory process

1. Edit `src/trace_editor/gui/main_window.ui` (Qt Designer XML).
2. Regenerate the Python UI:
   ```bash
   pyside6-uic src/trace_editor/gui/main_window.ui \
     -o src/trace_editor/gui/ui_main_window.py
   ```
3. Never touch `ui_main_window.py` directly — it will be overwritten.

## Docstrings

Use **reStructuredText `:param:` / `:return:` style** throughout, for two reasons:

1. The API docs under [`docs/api/`](docs/api/) are generated from docstrings.
2. `prepare_fit` in `core/h5_data.py` **parses** the docstrings of the
   ForceSpectroscopy models to extract parameter display names, units, and fit
   bounds. Their exact format matters:

   ```
   :param x: F (pN) # Force
   :param L: L (bp) [1000, 20000] # Contour length
   :return: Z (um) # extension
   ```

   The text after the field name up to `#` becomes the display label; a
   `[min, max]` sets the fit bounds. **Do not reformat the `fs` model
   docstrings** without checking `prepare_fit` still parses them.

### Regenerating the API docs

```bash
uv sync --extra docs                        # once, to install pydoc-markdown
uv run python scripts/gen_api_docs.py       # writes docs/api/*.md
```

Run this after changing public docstrings so the agent-readable API docs stay
current.

A pre-commit hook does this automatically. Enable it once per clone:

```bash
git config core.hooksPath .githooks
```

When you commit a change to any `src/**/*.py`, the hook
([`.githooks/pre-commit`](.githooks/pre-commit)) regenerates `docs/api/` and
stages it into the same commit. It requires the `docs` extra installed
(`uv sync --extra docs`). Commits that don't touch `src/` skip the hook.

## What is NOT yet implemented

- **Histogram / parameter distribution plot** — see
  `TraceEditorGUIs.HistogramGUI` in the original.
- **Export to CSV/Excel** — see `TraceIO.export()` in the original.
- **Force calibration dialog** — currently uses hardcoded defaults; see the
  `TraceEditor.py` calibration dialog.
