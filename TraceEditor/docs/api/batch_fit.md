# Table of Contents

* [trace\_editor.core.batch\_fit](#trace_editor.core.batch_fit)
  * [model\_wants\_multiprocessing](#trace_editor.core.batch_fit.model_wants_multiprocessing)
  * [RoiFitResult](#trace_editor.core.batch_fit.RoiFitResult)
    * [param\_values](#trace_editor.core.batch_fit.RoiFitResult.param_values)
  * [BatchFitError](#trace_editor.core.batch_fit.BatchFitError)
  * [run\_batch](#trace_editor.core.batch_fit.run_batch)

<a id="trace_editor.core.batch_fit"></a>

# trace\_editor.core.batch\_fit

Batch fitting across many ROIs, parallelized with a process pool.

Qt-free: this module is imported by worker processes spawned via
`concurrent.futures.ProcessPoolExecutor`, so it (and everything it imports)
must not import PySide6 — on macOS/Windows the default "spawn" start method
re-imports this module fresh in every worker process, and pulling in Qt there
would be wasted work at best and unsafe at worst (Qt objects are not
fork/spawn-safe). The GUI layer (`gui/main_window.py`) adapts this module's
plain-callable interface (`progress_cb`/`stop_cb`) onto Qt signals.

`_correct_z_for_fit`, `_get_fit_selection`, `_store_fit_result` and
`_exponential_drift` live here (not in `gui/main_window.py`) for the same
reason: they run inside worker processes as part of `_fit_one_roi`, so they
must be reachable without importing Qt. `gui/main_window.py` imports them
from here for its single-ROI (non-batch) fit path, so there is one
implementation shared by both.

<a id="trace_editor.core.batch_fit.model_wants_multiprocessing"></a>

#### model\_wants\_multiprocessing

```python
def model_wants_multiprocessing(model_name: str) -> bool
```

Whether `model_name` opts into batch-fit multiprocessing via its docstring.

Benchmarking (see `scripts/benchmark_batch_fit.py`) showed the process
pool is only worth its fixed per-worker startup cost (spawn + import,
roughly 1s x worker count) for models expensive enough that a single fit
takes on the order of a second or more -- e.g. `zHMf`/`zHMf2`, which loop
or call sub-models multiple times per evaluation. For cheap models (e.g.
`ExponentialDrift`, `zWLC`) the pool is slower than a plain sequential
loop, sometimes by 5-10x, so pooling must be opt-in per model rather than
always-on.

A model marks itself as worth pooling with a `:multiprocessing: true`
line in its docstring (see `zHMf`/`zHMf2` in `force_spectroscopy.py`),
following the same `:tag:` convention `prepare_fit` already uses for
`:param:`/`:return:`. Looked up by *exact line match* (not the
substring/`str.split` approach `prepare_fit` uses for parameter ranges)
so this tag can never be mistaken for part of a `:param <name>:` line or
vice versa, regardless of what any parameter is named.

**Arguments**:

- `model_name`: a ForceSpectroscopy model name, e.g. `"zHMf2"`.

**Returns**:

`True` if the docstring contains a `:multiprocessing: true` line.

<a id="trace_editor.core.batch_fit.RoiFitResult"></a>

## RoiFitResult Objects

```python
@dataclass
class RoiFitResult()
```

One ROI's completed fit, ready to be written to its h5 file.

Picklable: carries only plain floats and strings across the process
boundary, not the `lmfit.model.ModelResult`/`H5Data` objects used to
compute it (those are not guaranteed to pickle cleanly and are not
needed by the caller, which only writes `param_values`/`r2`/`r2_full`
back to disk).

<a id="trace_editor.core.batch_fit.RoiFitResult.param_values"></a>

#### param\_values

label -> (value, stderr)

<a id="trace_editor.core.batch_fit.BatchFitError"></a>

## BatchFitError Objects

```python
@dataclass
class BatchFitError()
```

One ROI's fit that raised an exception, carrying enough to report it.

<a id="trace_editor.core.batch_fit.run_batch"></a>

#### run\_batch

```python
def run_batch(tasks,
              model_name,
              lmfit_params,
              names,
              max_workers=None,
              progress_cb=None,
              stop_cb=None,
              error_cb=None,
              use_pool=None,
              selection_mode="default",
              criteria=None)
```

Fit `model_name` to every (filepath, roi) task, and write results.

Runs in a `ProcessPoolExecutor` only when `model_name` opts in (see
`model_wants_multiprocessing`) or `use_pool=True` is forced; otherwise
falls back to a plain sequential loop over `_fit_one_roi` in this
process. This matters: benchmarking showed pooling is a net loss (import
+ IPC overhead exceeds the fit cost) for cheap models like
`ExponentialDrift`/`zWLC`, and a large win (3-8x) for expensive models
like `zHMf`/`zHMf2` where a single fit takes on the order of a second or
more. Callers should not have to know this per model — the model's own
docstring says whether it's worth pooling.

When pooling, the pool is created, used, and fully shut down (via
`with`) before this function returns, so callers never need to manage
worker-process lifetime themselves. Each ROI's fit runs in a worker
process (`_fit_one_roi`); result writes back to the h5 file happen here,
in the caller's process, so concurrent workers never write the same file
at the same time. The sequential path writes the same way, once at the
end, so callers see identical behavior either way.

Qt-free and independently testable: `progress_cb`/`stop_cb` are plain
callables, not signals, so this can be driven from a script or test with
no `QApplication`. `gui/main_window.py`'s `FitWorker.run()` is a thin
adapter that supplies callables which emit Qt signals / check a
`stop_requested` flag.

**Arguments**:

- `tasks`: list of `(filepath, roi)` pairs to fit.
- `model_name`: a ForceSpectroscopy model name, e.g. `"zWLC"`.
- `lmfit_params`: `lmfit.Parameters` with vary flags set by the caller;
deep-copied per ROI inside each worker so ROIs fit independently.
- `names`: channel/label map from `prepare_fit`.
- `max_workers`: worker process count when pooling; defaults to one
less than the machine's CPU count (minimum 1), leaving a core free
for the caller's own thread/event loop. Ignored when not pooling.
- `progress_cb`: optional `(done: int, total: int) -> None`, called
after each ROI completes (success or failure).
- `stop_cb`: optional `() -> bool`; polled between completions, and if
it returns `True`, remaining tasks are abandoned (pooled: not-yet-
started futures are cancelled, already-running ones are still
awaited so the pool can shut down cleanly; sequential: the loop just
stops) and `run_batch` returns early.
- `error_cb`: optional `(error: BatchFitError) -> None`, called for
each ROI whose fit raises an exception (mirrors the previous
sequential `FitWorker`'s per-ROI try/except -> status message).
- `use_pool`: override the docstring-driven auto-detection —
`True`/`False` forces pooled/sequential regardless of `model_name`'s
`:multiprocessing:` tag; `None` (default) auto-detects.
- `selection_mode`: `"default"` or `"manual"` — forwarded to
`_get_fit_selection` for every ROI (see there). In `"manual"` mode
each ROI uses its own persisted `Traces/Selection/<roi>` mask.
- `criteria`: selection criteria to evaluate in `"default"` mode
(already resolved by the caller from the active pipeline, if any) —
forwarded to `_get_fit_selection` for every ROI.

**Returns**:

list of `RoiFitResult` for ROIs that fit successfully. ROIs that
raised are reported via `error_cb`, not included here or raised.

