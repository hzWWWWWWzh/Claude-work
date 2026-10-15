# Table of Contents

* [trace\_editor.gui.main\_window](#trace_editor.gui.main_window)
  * [FileBrowserDialog](#trace_editor.gui.main_window.FileBrowserDialog)
  * [SelectPointsDialog](#trace_editor.gui.main_window.SelectPointsDialog)
  * [SelectTracesDialog](#trace_editor.gui.main_window.SelectTracesDialog)
  * [FitDialog](#trace_editor.gui.main_window.FitDialog)
    * [result\_token](#trace_editor.gui.main_window.FitDialog.result_token)
  * [FitWorker](#trace_editor.gui.main_window.FitWorker)
  * [PipelineWorker](#trace_editor.gui.main_window.PipelineWorker)
    * [status](#trace_editor.gui.main_window.PipelineWorker.status)
  * [TraceEditor](#trace_editor.gui.main_window.TraceEditor)
    * [open\_files](#trace_editor.gui.main_window.TraceEditor.open_files)
    * [open\_single\_file](#trace_editor.gui.main_window.TraceEditor.open_single_file)
    * [redraw](#trace_editor.gui.main_window.TraceEditor.redraw)
    * [remove\_model\_viewed\_trace](#trace_editor.gui.main_window.TraceEditor.remove_model_viewed_trace)
    * [remove\_model\_all\_traces](#trace_editor.gui.main_window.TraceEditor.remove_model_all_traces)
    * [remove\_all\_processing](#trace_editor.gui.main_window.TraceEditor.remove_all_processing)
    * [calculate\_force](#trace_editor.gui.main_window.TraceEditor.calculate_force)

<a id="trace_editor.gui.main_window"></a>

# trace\_editor.gui.main\_window

<a id="trace_editor.gui.main_window.FileBrowserDialog"></a>

## FileBrowserDialog Objects

```python
class FileBrowserDialog(QDialog)
```

Open File dialog: a sortable/queryable table of `.h5` metadata for one
folder, so files can be triaged (axis, tags, selected/total beads, LUT
file used) before opening them.

Populated via `core.file_browser.scan_folder`, which reads only root
attributes, `/Settings`, and the small `/Parameters` array — no trace data
is loaded.

<a id="trace_editor.gui.main_window.SelectPointsDialog"></a>

## SelectPointsDialog Objects

```python
class SelectPointsDialog(QDialog)
```

Channel min/max/derivative/invert table with live preview on current ROI.

<a id="trace_editor.gui.main_window.SelectTracesDialog"></a>

## SelectTracesDialog Objects

```python
class SelectTracesDialog(QDialog)
```

Parameter min/max/invert table; OK applies to all selected ROIs.

<a id="trace_editor.gui.main_window.FitDialog"></a>

## FitDialog Objects

```python
class FitDialog(QDialog)
```

Modal dialog for setting fit parameters.

Columns: Name | Value | Type (Global/Local/Fit) | Min | Max
Buttons: Cancel | Fit Current | Fit All

On Fit Current → caller fits current ROI with these params.
On Fit All → caller fits all selected ROIs with these params.
On Cancel (or closing via the window's X) → live-preview writes made
while the dialog was open are reverted, restoring the H5 parameters
for this ROI to their pre-dialog state.

`live_preview_fn` is called whenever any value/type changes so the
main plot updates in real time (matching original apply_model behaviour).

<a id="trace_editor.gui.main_window.FitDialog.result_token"></a>

#### result\_token

"OK" | "Batch" | "Cancel"

<a id="trace_editor.gui.main_window.FitWorker"></a>

## FitWorker Objects

```python
class FitWorker(QObject)
```

Runs a batch fit on a background QThread by delegating to
`core.batch_fit.run_batch`, which parallelizes the per-ROI fits across a
`ProcessPoolExecutor`.

This class is a thin Qt adapter: it owns no fitting logic itself, only
translates `run_batch`'s plain `progress_cb`/`stop_cb`/`error_cb`
callables into Qt signals and the `stop_requested` flag, so the QThread
plumbing in `TraceEditor._start_batch` (and the crash-safety fixes in
`_clear_batch_thread`/`closeEvent`) doesn't need to change.

<a id="trace_editor.gui.main_window.PipelineWorker"></a>

## PipelineWorker Objects

```python
class PipelineWorker(QObject)
```

Runs a `core.auto_fit_pipeline.run_pipeline` call on a background
`QThread` for one or more files, one after another.

Thin Qt adapter, same shape as `FitWorker` — `run_pipeline`'s
`progress_cb`/`error_cb`/`status_cb` are plain callables, translated here
into Qt signals so `TraceEditor._run_pipeline_on_files` can reuse the
same thread-lifecycle handling `_start_batch` already has.

<a id="trace_editor.gui.main_window.PipelineWorker.status"></a>

#### status

model_name, step, total_steps, file, total_files

<a id="trace_editor.gui.main_window.TraceEditor"></a>

## TraceEditor Objects

```python
class TraceEditor(QMainWindow, Ui_MainWindow)
```

Main application window.

A 3-column PySide6 window (model/params | plot | axis controls) over one or
more open `H5Data` files. Holds the open files, the current file/ROI
selection, and the matplotlib canvas. Single-ROI fitting delegates to
`core.h5_data.apply_fit` directly (`_run_fit`); batch fitting delegates to
`core.batch_fit.run_batch` via the background `FitWorker`/`QThread`
(`_start_batch`), which parallelizes per-ROI fits across a process pool.
See `docs/architecture.md` for the layout, the Z-correction logic, and the
batch/selection flow.

<a id="trace_editor.gui.main_window.TraceEditor.open_files"></a>

#### open\_files

```python
def open_files()
```

Prompt for `.h5` files and load them.

Opens the `FileBrowserDialog` (a table of file metadata — axis, Notes
tag, selected/total beads — for the chosen folder), pre-highlighting
rows already open elsewhere in the app. Registers each newly selected
file as an `H5Data`, closes any deselected already-open file, adds
newly opened files to the file combo, and activates the last newly
added file. LUT files are skipped. Errors per file are shown in the
status bar, not raised. If auto-fit-on-open is enabled, the active
pipeline is run once across every newly added file.

<a id="trace_editor.gui.main_window.TraceEditor.open_single_file"></a>

#### open\_single\_file

```python
def open_single_file(path: str) -> bool
```

Register `path` as an open `H5Data` file and add it to the file combo.

Shared by `open_files()` (manual, interactive) and the in-app folder
watcher (`_check_watched_folder`). Does not activate the file — call
`_activate_opened_file(path)` after, if `path` was newly added.
Returns whether `path` is now open (`True` if newly added or already
open; `False` on error or if it's a LUT file).

<a id="trace_editor.gui.main_window.TraceEditor.redraw"></a>

#### redraw

```python
def redraw(*_)
```

Redraw the plot for the current file/ROI from current control values.

Scatters the selected X vs Y channels (with drift/Zbead correction
applied to Z/X/Y), overlays the current model curve if one is chosen, and
refreshes the parameter list. A no-op with no file/ROI selected. Called
on every control change and after each fit.

<a id="trace_editor.gui.main_window.TraceEditor.remove_model_viewed_trace"></a>

#### remove\_model\_viewed\_trace

```python
def remove_model_viewed_trace()
```

Remove the selected model's fit parameters from the currently viewed trace only.

<a id="trace_editor.gui.main_window.TraceEditor.remove_model_all_traces"></a>

#### remove\_model\_all\_traces

```python
def remove_model_all_traces()
```

Remove the selected model's fit parameters from every ROI in every open file.

<a id="trace_editor.gui.main_window.TraceEditor.remove_all_processing"></a>

#### remove\_all\_processing

```python
def remove_all_processing()
```

Remove all calculated parameters, F (pN), and manual point
selection masks from all open files -- everything derived from the
raw instrument data, which is left untouched.

<a id="trace_editor.gui.main_window.TraceEditor.calculate_force"></a>

#### calculate\_force

```python
def calculate_force()
```

Convert Shift (mm) → F (pN) for every open file via `MagnetForce`.

Uses the active pipeline's `calculate_force` step params (falling back
to `[force_defaults]`, via `resolve_force_defaults` — the same
resolution `run_pipeline` uses for a `calculate_force` step), writes
the result as the shared `F (pN)` dataset, and marks every ROI with
`F (pN) calculated = 1.0`. Runs across all open files with a progress
bar.

