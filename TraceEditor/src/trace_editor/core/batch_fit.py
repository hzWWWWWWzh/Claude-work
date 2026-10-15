"""Batch fitting across many ROIs, parallelized with a process pool.

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
"""

import math
import os
from concurrent.futures import ProcessPoolExecutor, as_completed
from dataclasses import dataclass, field
from os import process_cpu_count

import numpy as np

# H5Data never holds a file open across calls (every read/write is its own
# short-lived open/close), so there is no window where two processes could
# corrupt an in-progress write -- HDF5's file locking only guards against
# that, and otherwise just rejects the legitimate concurrent short opens this
# module relies on (a worker process reading an ROI while the parent process
# writes back a different, already-completed ROI's result). Set unless the
# caller has deliberately overridden it. Environment variables are inherited
# by spawned worker processes, so this must be set before the pool starts.
os.environ.setdefault("HDF5_USE_FILE_LOCKING", "FALSE")

from . import force_spectroscopy as fs
from .h5_data import H5Data, apply_fit, r2_label, r2_full_label


def model_wants_multiprocessing(model_name: str) -> bool:
    """Whether `model_name` opts into batch-fit multiprocessing via its docstring.

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

    :param model_name: a ForceSpectroscopy model name, e.g. `"zHMf2"`.
    :return: `True` if the docstring contains a `:multiprocessing: true` line.
    """
    func = getattr(fs, model_name, None)
    if func is None:
        return False
    doc = func.__doc__ or ""
    for line in doc.splitlines():
        if line.strip().startswith(":multiprocessing:"):
            value = line.strip()[len(":multiprocessing:"):].strip().split("#")[0].strip()
            return value.lower() == "true"
    return False


def _exponential_drift(time, Z0=0.0, dZ_dt=0.0, tau=np.inf):
    """Evaluate the exponential drift model matching ForceSpectroscopy.ExponentialDrift.

    :param time: time array (s).
    :param Z0: extension offset at t=0 (um).
    :param dZ_dt: drift rate (nm/s).
    :param tau: exponential decay time constant (s); `inf` gives linear drift.
    :return: drift array (um), same shape as `time`.
    """
    if dZ_dt == 0.0:
        return np.zeros_like(time, dtype=float) + Z0
    tau = min(float(tau), 1e9)
    a = -dZ_dt / 1000.0 * tau   # convert nm/s to um/s, scale by tau
    return a * np.exp(-time / tau) + (Z0 - a)


def _correct_z_for_fit(y_data: np.ndarray, model_name: str, names: dict,
                        h5data: H5Data, roi: str,
                        zbead_baseline: float | None = None) -> np.ndarray:
    """Apply drift + Zbead correction to y before fitting, matching TraceIO.fit() line 1030-1031.

    For all models except ExponentialDrift and zMagnetForce, when the model y-channel
    is Z (um), the original subtracts drift and Zbead so the fit sees corrected Z.
    This keeps z_bead in zWLC fitting to a small residual offset (~0) rather than
    absorbing the full bead height (~2-8 um).

    :param y_data: raw y-channel data for this ROI.
    :param model_name: a ForceSpectroscopy model name, e.g. `"zWLC"`.
    :param names: channel/label map from `prepare_fit`.
    :param h5data: the file containing `roi`.
    :param roi: zero-padded ROI id.
    :param zbead_baseline: Zbead value to subtract, overriding the value
        stored in `h5data`. Pass the dialog's pre-edit snapshot during live
        preview, where `h5data`'s stored Zbead may already reflect an
        in-progress (uncommitted) edit — see `FitDialog._on_param_changed`.
    :return: corrected copy of `y_data` (or `y_data` unchanged if correction
        does not apply to this model/channel).
    """
    if model_name in ("ExponentialDrift", "zMagnetForce"):
        return y_data
    if names.get("y") != "Z (um)":
        return y_data

    y = y_data.copy()
    # Subtract exponential drift
    try:
        time = h5data.get_channel("Time (s)", roi)
        Z0    = h5data.get_parameter("Z0 (um)",      roi, default=0.0)
        dZ_dt = h5data.get_parameter("dZ_dt (nm_s)", roi, default=0.0)
        tau   = h5data.get_parameter("tau (s)",       roi, default=np.inf)
        if math.isnan(Z0):    Z0 = 0.0
        if math.isnan(dZ_dt): dZ_dt = 0.0
        if math.isnan(tau):   tau = np.inf
        y -= _exponential_drift(time, Z0, dZ_dt, tau)
    except Exception:
        pass

    # Subtract Zbead
    if zbead_baseline is not None:
        zbead = zbead_baseline
    else:
        zbead = h5data.get_parameter("Zbead (um)", roi, default=0.0)
        if math.isnan(zbead):
            zbead = 0.0
    y -= zbead

    return y


def _resolve_bound(value, trace: np.ndarray):
    """Resolve one criterion bound: a plain number, or a percentile string like `"p90"`.

    `"p90"` means "the 90th percentile of `trace`'s values for this ROI" —
    used by `[model_selection.<model>]` defaults (e.g. ExponentialDrift's
    "only the lowest-force points", expressed as `min = "p90"` on the
    Shift (mm) channel) so the bound adapts per-trace instead of being a
    fixed number that would only work for one instrument's calibration.

    :param value: `""` (unbounded), a numeric string, or `"p<0-100>"`.
    :param trace: the channel's data for this ROI, used for percentile bounds.
    :return: the resolved float bound, or `None` for `""`.
    """
    if value == "":
        return None
    if isinstance(value, str) and value.startswith("p"):
        return float(np.percentile(trace, float(value[1:])))
    return float(value)


def _evaluate_criteria(criteria: list, h5data: H5Data, roi: str) -> np.ndarray:
    """Evaluate point-level selection criteria against one ROI's data.

    Pure/read-only — mirrors TraceIO.select_points()'s masking logic without
    writing anything to `h5data`. Each criterion: {channel, min, max,
    derivative, invert}; `min`/`max` may be a plain number or a percentile
    string (see `_resolve_bound`). Z (um) is drift+zbead corrected before
    thresholding (matching original line 965-968).

    Used both by `_apply_point_selection` (which persists the result — the
    explicit SelectPointsDialog/"in view" action) and by `_get_fit_selection`
    (which re-evaluates live for preview/fitting and must not have side
    effects, since it runs on every plot redraw).

    :return: float array (1.0/0.0), same length as the ROI's channels.
    """
    channels = h5data.all_channels()
    n = len(h5data.get_channel(channels[0], roi))
    selection = np.ones(n, dtype=float)

    for c in criteria:
        lo = c.get("min", "")
        hi = c.get("max", "")
        if lo == "" and hi == "":
            continue

        trace = h5data.get_channel(c["channel"], roi).copy()

        if c["channel"] == "Z (um)":
            time = h5data.get_channel("Time (s)", roi)
            Z0    = h5data.get_parameter("Z0 (um)",      roi, default=0.0)
            dZ_dt = h5data.get_parameter("dZ_dt (nm_s)", roi, default=0.0)
            tau   = h5data.get_parameter("tau (s)",       roi, default=np.inf)
            zbead = h5data.get_parameter("Zbead (um)",   roi, default=0.0)
            if math.isnan(Z0):    Z0 = 0.0
            if math.isnan(dZ_dt): dZ_dt = 0.0
            if math.isnan(tau):   tau = np.inf
            if math.isnan(zbead): zbead = 0.0
            trace -= _exponential_drift(time, Z0, dZ_dt, tau) + zbead

        if c.get("derivative", False):
            trace = np.diff(trace, prepend=np.nan)

        lo_is_percentile = isinstance(lo, str) and lo.startswith("p")
        hi_is_percentile = isinstance(hi, str) and hi.startswith("p")
        lo_val = _resolve_bound(lo, trace)
        hi_val = _resolve_bound(hi, trace)
        lo = -np.inf if lo_val is None else lo_val
        hi =  np.inf if hi_val is None else hi_val

        # Percentile bounds are inclusive at the boundary (matching the
        # original TraceIO.fit()'s ">=", e.g. ExponentialDrift's Shift >=
        # p90): a plateaued channel can have many points exactly at its
        # percentile value (e.g. a trace that saturates at its max Shift),
        # and a strict "<"/">" would silently select nothing. Plain numeric
        # bounds keep the existing exclusive "<"/">" semantics used by
        # SelectPointsDialog's manual range entry.
        lo_mask = (trace >= lo) if lo_is_percentile else (trace > lo)
        hi_mask = (trace <= hi) if hi_is_percentile else (trace < hi)
        mask = lo_mask & hi_mask
        if c.get("invert", False):
            mask = ~mask
        selection *= mask.astype(float)

    return selection


def _apply_point_selection(criteria: list, h5data: H5Data, roi: str) -> np.ndarray:
    """Evaluate point-level selection criteria and persist the mask to h5.

    Used by the explicit `SelectPointsDialog`/"Select points in view"
    actions, where the user is deliberately setting the ROI's stored
    selection. See `_evaluate_criteria` for the masking logic itself.
    """
    selection = _evaluate_criteria(criteria, h5data, roi)
    h5data.set_parameter(" selected_points", float(np.sum(selection)), roi=roi)
    h5data.set_selection_mask(roi, selection)
    return selection


def _get_fit_selection(h5data: H5Data, roi: str, mode: str = "default",
                        criteria: list | None = None) -> np.ndarray | None:
    """Return a boolean mask selecting which data points to use for fitting.

    - `mode="default"`: evaluates `criteria` fresh against this ROI's data
      (percentile bounds like `"p90"` are per-trace, so this can't be
      precomputed once and reused across ROIs). No `criteria` (or an empty
      list) uses every point. Resolving *which* criteria apply for a given
      model/pipeline is the caller's job (see
      `core.auto_fit_pipeline.resolve_criteria`) — this function only
      evaluates whatever list it's given, and has no notion of "pipeline".
    - `mode="manual"`: uses this ROI's own persisted point-selection mask at
      `Traces/Selection/<roi>` (set via `SelectPointsDialog`/"Select points
      in view" — see `_apply_point_selection`), if any; otherwise every
      point. Each ROI in a batch run reads its own persisted mask, so
      "Fit All" in manual mode applies whatever selection each ROI already
      has stored, not one shared criteria list.

    :param h5data: the file containing `roi`.
    :param roi: zero-padded ROI id.
    :param mode: `"default"` or `"manual"`.
    :param criteria: selection criteria to evaluate in `"default"` mode
        (same shape as `SelectPointsDialog`'s table); `None`/`[]` means
        every point.
    :return: boolean mask, or `None` to use every point.
    """
    if mode == "manual":
        mask = h5data.get_selection_mask(roi)
        return mask == 1.0 if mask is not None else None

    if not criteria:
        return None

    sel = _evaluate_criteria(criteria, h5data, roi)
    return sel == 1.0


def _resolve_fit_param(p: str, names: dict, result, model_name: str,
                        h5data: H5Data, roi: str,
                        zbead_baseline: float | None = None):
    """Compute the (label, value, stderr) to store for one fitted parameter.

    Because `_correct_z_for_fit` subtracts stored Zbead before fitting, the
    fitted z_bead is a residual. The true Zbead = stored_Zbead + fitted_z_bead,
    so z_bead is accumulated against the currently stored value rather than
    replacing it, keeping repeated fit iterations stable.

    :param p: lmfit parameter key.
    :param names: channel/label map from `prepare_fit`.
    :param result: an `lmfit.model.ModelResult` (or compatible) with `.params`.
    :param model_name: a ForceSpectroscopy model name, e.g. `"zWLC"`.
    :param h5data: the file containing `roi`.
    :param roi: zero-padded ROI id.
    :param zbead_baseline: Zbead value to accumulate onto, overriding the
        value stored in `h5data`. Must match whatever `zbead_baseline` was
        passed to `_correct_z_for_fit` for this same fit, so the residual is
        added back onto the baseline it was measured against.
    :return: `(label, value, stderr)`, or `None` if `p` is not a storable
        parameter (the model's x/y channel keys).
    """
    label = names.get(p, p)
    if label in ("x", "y"):
        return None
    val = float(result.params[p].value)
    stderr = result.params[p].stderr or 0.0
    if p == "z_bead" and model_name not in ("ExponentialDrift", "zMagnetForce"):
        if zbead_baseline is not None:
            stored = zbead_baseline
        else:
            stored = h5data.get_parameter("Zbead (um)", roi, default=0.0)
            if math.isnan(stored):
                stored = 0.0
        val = stored + val
    return label, val, stderr


def _store_fit_result(p: str, names: dict, result, model_name: str,
                       h5data: H5Data, roi: str,
                       zbead_baseline: float | None = None):
    """Write one fit parameter to h5, accumulating z_bead into stored Zbead.

    :param p: lmfit parameter key.
    :param names: channel/label map from `prepare_fit`.
    :param result: an `lmfit.model.ModelResult` (or compatible) with `.params`.
    :param model_name: a ForceSpectroscopy model name, e.g. `"zWLC"`.
    :param h5data: the file containing `roi`.
    :param roi: zero-padded ROI id.
    :param zbead_baseline: forwarded to `_resolve_fit_param` — see there.
    """
    resolved = _resolve_fit_param(p, names, result, model_name, h5data, roi,
                                   zbead_baseline=zbead_baseline)
    if resolved is None:
        return
    label, val, stderr = resolved
    h5data.set_parameter(label, val, roi=roi)
    h5data.set_parameter(f"{label} Error", stderr, roi=roi)


@dataclass
class RoiFitResult:
    """One ROI's completed fit, ready to be written to its h5 file.

    Picklable: carries only plain floats and strings across the process
    boundary, not the `lmfit.model.ModelResult`/`H5Data` objects used to
    compute it (those are not guaranteed to pickle cleanly and are not
    needed by the caller, which only writes `param_values`/`r2`/`r2_full`
    back to disk).
    """
    filepath: str
    roi: str
    model_name: str
    param_values: dict = field(default_factory=dict)   # label -> (value, stderr)
    r2: float = float("nan")
    r2_full: float = float("nan")


@dataclass
class BatchFitError:
    """One ROI's fit that raised an exception, carrying enough to report it."""
    filepath: str
    roi: str
    message: str

    def __str__(self):
        return f"ROI {self.roi}: {self.message}"


def _fit_one_roi(filepath: str, roi: str, model_name: str, lmfit_params,
                  names: dict, selection_mode: str = "default",
                  criteria: list | None = None) -> RoiFitResult:
    """Fit `model_name` to one ROI of one file. Runs inside a worker process.

    Opens its own `H5Data(filepath)` rather than receiving one from the
    caller — `H5Data` does no in-memory caching (every read/write reopens the
    file), so this is cheap, and it avoids pickling a live file handle across
    the process boundary. Only *reads* the file; parameter writes are done by
    the caller in the main process (see `run_batch`) so concurrent worker
    processes never write the same h5 file at the same time.

    :param filepath: path to the `.h5` file containing `roi`.
    :param roi: zero-padded ROI id.
    :param model_name: a ForceSpectroscopy model name, e.g. `"zWLC"`.
    :param lmfit_params: `lmfit.Parameters` with vary flags set by the caller
        (deep-copied per ROI internally so concurrent calls don't share state).
    :param names: channel/label map from `prepare_fit`.
    :param selection_mode: `"default"` or `"manual"` — see `_get_fit_selection`.
    :param criteria: selection criteria to evaluate in `"default"` mode
        (already resolved by the caller from the active pipeline, if any —
        see `core.auto_fit_pipeline.resolve_criteria`) — see `_get_fit_selection`.
    :return: `RoiFitResult` with the fitted parameter values, R2, and
        whole-trace R2.
    :raises Exception: any error raised while reading data or fitting is
        propagated to the caller via the `Future` (see `run_batch`).
    """
    import copy

    data = H5Data(filepath)
    params = copy.deepcopy(lmfit_params)

    # Seed per-ROI stored values for non-varying *Local* params only. Global
    # params are a fixed value shared across every trace (see FitDialog's
    # Global/Local/Fit docs) -- re-seeding them from this ROI's own stored
    # value would silently discard whatever the dialog was set to for every
    # ROI except the one the live preview already wrote it to. The Global/
    # Local/Fit type travels in lmfit's `user_data` slot (set by FitDialog),
    # not a bare custom attribute -- a plain attribute doesn't survive the
    # `copy.deepcopy` above or the pickling across the process-pool boundary.
    stored = data.list_parameters(roi)
    for p in params:
        label = names.get(p, p)
        if label in ("x", "y"):
            continue
        if not params[p].vary and (params[p].user_data or "Local") == "Local":
            stored_val = stored.get(label, float("nan"))
            if not math.isnan(stored_val):
                params[p].value = stored_val

    x = data.get_channel(names["x"], roi)
    y = data.get_channel(names["y"], roi)
    y = _correct_z_for_fit(y, model_name, names, data, roi)
    sel = _get_fit_selection(data, roi, mode=selection_mode, criteria=criteria)
    result, out_names = apply_fit(model_name, x, y, selection=sel, initial_parameters=params)

    param_values = {}
    for p in result.params:
        resolved = _resolve_fit_param(p, out_names, result, model_name, data, roi)
        if resolved is None:
            continue
        label, val, stderr = resolved
        param_values[label] = (val, stderr)

    r2 = float(getattr(result, "r2", np.nan))
    r2_full = float(getattr(result, "r2_full", np.nan))
    return RoiFitResult(filepath=filepath, roi=roi, model_name=model_name,
                         param_values=param_values, r2=r2, r2_full=r2_full)


def run_batch(tasks, model_name, lmfit_params, names,
              max_workers=None, progress_cb=None, stop_cb=None, error_cb=None,
              use_pool=None, selection_mode="default", criteria=None):
    """Fit `model_name` to every (filepath, roi) task, and write results.

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

    :param tasks: list of `(filepath, roi)` pairs to fit.
    :param model_name: a ForceSpectroscopy model name, e.g. `"zWLC"`.
    :param lmfit_params: `lmfit.Parameters` with vary flags set by the caller;
        deep-copied per ROI inside each worker so ROIs fit independently.
    :param names: channel/label map from `prepare_fit`.
    :param max_workers: worker process count when pooling; defaults to one
        less than the machine's CPU count (minimum 1), leaving a core free
        for the caller's own thread/event loop. Ignored when not pooling.
    :param progress_cb: optional `(done: int, total: int) -> None`, called
        after each ROI completes (success or failure).
    :param stop_cb: optional `() -> bool`; polled between completions, and if
        it returns `True`, remaining tasks are abandoned (pooled: not-yet-
        started futures are cancelled, already-running ones are still
        awaited so the pool can shut down cleanly; sequential: the loop just
        stops) and `run_batch` returns early.
    :param error_cb: optional `(error: BatchFitError) -> None`, called for
        each ROI whose fit raises an exception (mirrors the previous
        sequential `FitWorker`'s per-ROI try/except -> status message).
    :param use_pool: override the docstring-driven auto-detection —
        `True`/`False` forces pooled/sequential regardless of `model_name`'s
        `:multiprocessing:` tag; `None` (default) auto-detects.
    :param selection_mode: `"default"` or `"manual"` — forwarded to
        `_get_fit_selection` for every ROI (see there). In `"manual"` mode
        each ROI uses its own persisted `Traces/Selection/<roi>` mask.
    :param criteria: selection criteria to evaluate in `"default"` mode
        (already resolved by the caller from the active pipeline, if any) —
        forwarded to `_get_fit_selection` for every ROI.
    :return: list of `RoiFitResult` for ROIs that fit successfully. ROIs that
        raised are reported via `error_cb`, not included here or raised.
    """
    if not tasks:
        return []

    should_pool = use_pool if use_pool is not None else model_wants_multiprocessing(model_name)

    total = len(tasks)
    done = 0
    results = []

    def handle_outcome(outcome) -> bool:
        """Record one completed/failed fit and report progress.

        :return: whether to stop taking on further work (`stop_cb` fired).
        """
        nonlocal done
        if isinstance(outcome, RoiFitResult):
            results.append(outcome)
        elif error_cb is not None:
            error_cb(outcome)
        done += 1
        if progress_cb is not None:
            progress_cb(done, total)
        return stop_cb is not None and stop_cb()

    if should_pool:
        if max_workers is None:
            max_workers = max(1, (process_cpu_count() or 2) - 1)

        with ProcessPoolExecutor(max_workers=max_workers) as pool:
            futures = {
                pool.submit(_fit_one_roi, filepath, roi, model_name, lmfit_params, names,
                            selection_mode, criteria): (filepath, roi)
                for filepath, roi in tasks
            }
            for future in as_completed(futures):
                filepath, roi = futures[future]
                try:
                    outcome = future.result()
                except Exception as e:
                    outcome = BatchFitError(filepath=filepath, roi=roi, message=str(e))

                if handle_outcome(outcome):
                    for f in futures:
                        f.cancel()
                    break
    else:
        for filepath, roi in tasks:
            try:
                outcome = _fit_one_roi(filepath, roi, model_name, lmfit_params, names,
                                        selection_mode, criteria)
            except Exception as e:
                outcome = BatchFitError(filepath=filepath, roi=roi, message=str(e))

            if handle_outcome(outcome):
                break

    _write_results(results)
    return results


def _write_results(results: list):
    """Write completed fits' parameter values back to their h5 files.

    Runs in the caller's (main) process, never inside a worker, so
    concurrent ROI fits never write the same file at the same time.

    Grouped and written one file at a time (not one `H5Data.set_parameter`
    call per ROI per field) because `set_parameter` reads, mutates, and
    rewrites the *entire* `/Parameters` array on every call -- with N ROIs
    and several fields per ROI that is O(N^2) I/O overall, which dominated
    wall time in early benchmarking (an 8s bottleneck at 300 ROIs, dwarfing
    both the fitting itself and the process-pool overhead this module exists
    to amortize). Batching to one open/read/mutate-all/write/close cycle per
    file turns that into O(N).

    :param results: `RoiFitResult`s to write, in any order/mix of files.
    """
    by_file: dict[str, list[RoiFitResult]] = {}
    for r in results:
        by_file.setdefault(r.filepath, []).append(r)

    for filepath, file_results in by_file.items():
        data = H5Data(filepath)
        data.set_parameters_bulk({
            outcome.roi: _flatten_param_values(outcome)
            for outcome in file_results
        })


def _flatten_param_values(outcome: RoiFitResult) -> dict:
    """Expand one ROI's (value, stderr) pairs into the flat {field: value}
    layout `H5Data.set_parameters_bulk` expects, plus the R2 and whole-trace
    R2 fields.

    :param outcome: a successfully computed `RoiFitResult`.
    :return: `{field_name: value}`, e.g. `{"L (bp)": 8000.0, "L (bp) Error": 12.0, ...}`.
    """
    fields = {}
    for label, (val, stderr) in outcome.param_values.items():
        fields[label] = val
        fields[f"{label} Error"] = stderr
    fields[r2_label(outcome.model_name)] = outcome.r2
    fields[r2_full_label(outcome.model_name)] = outcome.r2_full
    return fields
