"""
Read/write the new h5 instrument format (h5py).

Parameters and selection are stored by extending the /Parameters compound
dataset with additional float64 fields, one row per ROI.
"""

import h5py
import numpy as np
import lmfit
from pathlib import Path

from . import force_spectroscopy as fs
from lmfit import Model


SELECTED_FIELD = " selected"
POINT_SELECTION_CHANNEL = "Selection"


class H5Data:
    """Read/write access to one new-format `.h5` instrument file.

    Wraps an on-disk h5py file. Trace data is read on demand; fit parameters and
    selection are stored by extending the `/Parameters` compound dataset with
    float64 fields, one row per ROI. There is no in-memory cache — every getter
    reads from disk and every setter writes to disk, so instances stay
    consistent across processes.

    See `docs/file-format.md` for the on-disk layout.
    """

    def __init__(self, filepath):
        """
        :param filepath: path to an existing `.h5` file.
        :raises FileNotFoundError: if the file does not exist.
        """
        self.filepath = Path(filepath)
        if not self.filepath.exists():
            raise FileNotFoundError(self.filepath)

    # ------------------------------------------------------------------
    # File identity
    # ------------------------------------------------------------------

    def is_lut_file(self) -> bool:
        """LUT files use a .LUT extension, never .h5. This always returns False."""
        return False

    def num_rois(self) -> int:
        """Return the number of tracked beads (ROIs) in the file."""
        with h5py.File(self.filepath, "r") as f:
            return int(f.attrs["num_rois"])

    def get_axis(self) -> str:
        """Return the `axis` root attribute (`"Shift_mm"` or `"Focus_mm"`), or empty string."""
        with h5py.File(self.filepath, "r") as f:
            raw = f.attrs.get("axis", "")
        if isinstance(raw, bytes):
            raw = raw.decode("utf-8", errors="replace")
        return str(raw)

    def get_user(self) -> str:
        """Return the `User` root attribute, or empty string."""
        with h5py.File(self.filepath, "r") as f:
            raw = f.attrs.get("User", "")
        if isinstance(raw, bytes):
            raw = raw.decode("utf-8", errors="replace")
        return str(raw)

    def set_notes(self, notes: str):
        """Overwrite the `Notes` root attribute."""
        with h5py.File(self.filepath, "r+") as f:
            f.attrs["Notes"] = notes

    def get_lut_filename(self) -> str:
        """Return the LUT file used to calculate bead height, as a bare
        filename with no path or extension (e.g. `"data_013"`).

        Read from `/Settings` (`Imaging`/`lut_filename`). Returns empty string
        if `/Settings` is absent (e.g. legacy-imported files) or the key isn't
        set.
        """
        with h5py.File(self.filepath, "r") as f:
            if "Settings" not in f:
                return ""
            settings = f["Settings"][:]
        for row in settings:
            group = row["Group"]
            key = row["Key"]
            if isinstance(group, bytes):
                group = group.decode("utf-8", errors="replace")
            if isinstance(key, bytes):
                key = key.decode("utf-8", errors="replace")
            if group == "Imaging" and key == "lut_filename":
                value = row["Value"]
                if isinstance(value, bytes):
                    value = value.decode("utf-8", errors="replace")
                return Path(str(value)).stem
        return ""

    def selected_count(self) -> int:
        """Return the number of ROIs with ` selected` = 1.0 (defaults to all if unset)."""
        with h5py.File(self.filepath, "r") as f:
            arr = self._read_params_array(f)
            if SELECTED_FIELD in arr.dtype.names:
                return int(np.count_nonzero(arr[SELECTED_FIELD] == 1.0))
        return int(arr.shape[0])

    def list_rois(self) -> list:
        """Return ROI ids as zero-padded 3-digit strings: `["000", "001", ...]`."""
        n = self.num_rois()
        return [f"{i:03d}" for i in range(n)]

    def get_notes(self) -> str:
        """Return Notes attribute from root group, or empty string."""
        with h5py.File(self.filepath, "r") as f:
            raw = f.attrs.get("Notes", "")
        if isinstance(raw, bytes):
            raw = raw.decode("utf-8", errors="replace")
        return str(raw)

    def get_date(self) -> str:
        """Return Date attribute from root group, or empty string."""
        with h5py.File(self.filepath, "r") as f:
            raw = f.attrs.get("Date", "")
        if isinstance(raw, bytes):
            raw = raw.decode("utf-8", errors="replace")
        return str(raw)

    # ------------------------------------------------------------------
    # Channel discovery
    # ------------------------------------------------------------------

    def shared_channels(self) -> list:
        """Channels shared across all ROIs (e.g. `Time (s)`, `Shift (mm)`, `F (pN)`).

        These are datasets directly under `/Traces`.
        """
        with h5py.File(self.filepath, "r") as f:
            traces = f["Traces"]
            return [k for k, v in traces.items() if isinstance(v, h5py.Dataset)]

    def roi_channels(self) -> list:
        """Per-ROI channels (e.g. `Z (um)`, `X (um)`, `A (a.u.)`).

        These are groups under `/Traces`, each holding one dataset per ROI.
        Excludes `Selection`, which is a point-selection mask, not a plottable
        or fittable physical channel.
        """
        with h5py.File(self.filepath, "r") as f:
            traces = f["Traces"]
            return [k for k, v in traces.items()
                    if isinstance(v, h5py.Group) and k != POINT_SELECTION_CHANNEL]

    def all_channels(self) -> list:
        """All channel names: shared channels followed by per-ROI channels."""
        return self.shared_channels() + self.roi_channels()

    # ------------------------------------------------------------------
    # Trace data
    # ------------------------------------------------------------------

    def get_shared(self, channel: str) -> np.ndarray:
        """Return a shared channel's data as float64.

        :param channel: a name from `shared_channels()`.
        :return: 1-D array, length = number of samples.
        """
        with h5py.File(self.filepath, "r") as f:
            return f["Traces"][channel][:].astype(np.float64)

    def get_trace(self, channel: str, roi: str) -> np.ndarray:
        """Return one ROI's trace for a per-ROI channel, as float64.

        :param channel: a name from `roi_channels()`.
        :param roi: zero-padded ROI id, e.g. `"003"`.
        :return: 1-D array, length = number of samples.
        """
        with h5py.File(self.filepath, "r") as f:
            return f["Traces"][channel][roi][:].astype(np.float64)

    def get_channel(self, channel: str, roi: str) -> np.ndarray:
        """Return data for any channel — shared or per-ROI — dispatching on type.

        :param channel: any name from `all_channels()`.
        :param roi: zero-padded ROI id; ignored for shared channels.
        :return: 1-D float64 array.
        """
        if channel in self.shared_channels():
            return self.get_shared(channel)
        return self.get_trace(channel, roi)

    # ------------------------------------------------------------------
    # /Parameters compound dataset helpers
    # ------------------------------------------------------------------

    def _read_params_array(self, f) -> np.ndarray:
        return f["Parameters"][:]

    def _write_params_array(self, f, arr: np.ndarray):
        del f["Parameters"]
        f.create_dataset("Parameters", data=arr)

    def _params_as_dict_list(self) -> list:
        """Return list of dicts (one per ROI), preserving all fields."""
        with h5py.File(self.filepath, "r") as f:
            arr = self._read_params_array(f)
        result = []
        for row in arr:
            result.append({name: row[name] for name in arr.dtype.names})
        return result

    # ------------------------------------------------------------------
    # Selection
    # ------------------------------------------------------------------

    def get_selection(self, roi: str) -> float:
        """Return this ROI's selection flag.

        :param roi: zero-padded ROI id.
        :return: 1.0 if selected (included in batch ops), 0.0 if rejected.
            Defaults to 1.0 when the ` selected` field has not been written yet.
        """
        roi_idx = int(roi)
        with h5py.File(self.filepath, "r") as f:
            arr = self._read_params_array(f)
            if SELECTED_FIELD in arr.dtype.names:
                return float(arr[roi_idx][SELECTED_FIELD])
        return 1.0

    def set_selection(self, roi: str, value: float):
        """Set this ROI's selection flag (1.0 = selected, 0.0 = rejected)."""
        self.set_parameter(SELECTED_FIELD, value, roi=roi)

    def get_selection_mask(self, roi: str) -> np.ndarray | None:
        """Return this ROI's per-point selection mask (1.0/0.0 float array).

        Stored at `/Traces/Selection/<roi>`, matching the layout of the other
        per-ROI channels. Returns None if no point selection has been made yet
        for this ROI (all points should then be treated as selected).
        """
        with h5py.File(self.filepath, "r") as f:
            group = f["Traces"].get(POINT_SELECTION_CHANNEL)
            if group is None or roi not in group:
                return None
            return group[roi][:].astype(np.float64)

    def set_selection_mask(self, roi: str, mask: np.ndarray):
        """Write this ROI's per-point selection mask, creating the dataset if needed.

        :param mask: 1-D array (1.0 = selected, 0.0 = rejected), length = n_samples.
        """
        data = np.asarray(mask, dtype=np.float32)
        with h5py.File(self.filepath, "r+") as f:
            group = f["Traces"].require_group(POINT_SELECTION_CHANNEL)
            if roi in group:
                del group[roi]
            group.create_dataset(roi, data=data)

    # ------------------------------------------------------------------
    # Parameters
    # ------------------------------------------------------------------

    def get_parameter(self, name: str, roi: str, default=np.nan) -> float:
        """Return a stored parameter value for one ROI.

        :param name: parameter field name, e.g. `"L (bp)"`.
        :param roi: zero-padded ROI id.
        :param default: value returned when the field is absent or stored as NaN.
        :return: the stored float, or `default`.
        """
        roi_idx = int(roi)
        with h5py.File(self.filepath, "r") as f:
            arr = self._read_params_array(f)
            if name in arr.dtype.names:
                val = float(arr[roi_idx][name])
                return val if not np.isnan(val) else default
        return default

    def set_parameter(self, name: str, value: float, roi: str):
        """Add or update a float parameter field for one ROI.

        Rewrites the `/Parameters` dataset. If `name` is a new field it is added
        as a float64 column, defaulting to NaN for every other ROI.

        :param name: parameter field name.
        :param value: value to store for this ROI.
        :param roi: zero-padded ROI id.
        """
        roi_idx = int(roi)
        with h5py.File(self.filepath, "r+") as f:
            arr = self._read_params_array(f)
            if name not in arr.dtype.names:
                # Extend dtype with new float64 field
                new_dtype = np.dtype(arr.dtype.descr + [(name, "<f8")])
                new_arr = np.zeros(arr.shape, dtype=new_dtype)
                for field in arr.dtype.names:
                    new_arr[field] = arr[field]
                # Default new field to NaN for all ROIs
                new_arr[name] = np.nan
                arr = new_arr
            arr[roi_idx][name] = value
            self._write_params_array(f, arr)

    def set_parameters_bulk(self, values_by_roi: dict):
        """Write many parameter fields for many ROIs in one open/read/write pass.

        Equivalent to calling `set_parameter` once per (roi, field) pair, but
        reads and rewrites the `/Parameters` array only once in total instead
        of once per call -- `set_parameter` is O(rows) per call (it rewrites
        the whole compound dataset each time), so writing N ROIs x M fields
        one field at a time is O(N*M*N); this is O(N*M).  Intended for batch
        fitting, where many ROIs' results need to land in one file at once.

        :param values_by_roi: `{roi: {field_name: value}}`. New field names
            are added as float64 columns, defaulting to NaN for every ROI not
            present in that field's `values_by_roi` entry.
        """
        if not values_by_roi:
            return
        with h5py.File(self.filepath, "r+") as f:
            arr = self._read_params_array(f)

            new_fields = sorted({
                name for fields in values_by_roi.values() for name in fields
                if name not in arr.dtype.names
            })
            if new_fields:
                new_dtype = np.dtype(arr.dtype.descr + [(name, "<f8") for name in new_fields])
                new_arr = np.zeros(arr.shape, dtype=new_dtype)
                for field in arr.dtype.names:
                    new_arr[field] = arr[field]
                for name in new_fields:
                    new_arr[name] = np.nan
                arr = new_arr

            for roi, fields in values_by_roi.items():
                roi_idx = int(roi)
                for name, value in fields.items():
                    arr[roi_idx][name] = value

            self._write_params_array(f, arr)

    def list_parameters(self, roi: str) -> dict:
        """Return {name: value} for all analysis parameters of this ROI.

        Excludes the instrument fields ROI, X (pix), Y (pix).
        """
        instrument_fields = {"ROI", "X (pix)", "Y (pix)"}
        roi_idx = int(roi)
        with h5py.File(self.filepath, "r") as f:
            arr = self._read_params_array(f)
        result = {}
        for name in arr.dtype.names:
            if name not in instrument_fields:
                val = float(arr[roi_idx][name])
                result[name] = val
        return result

    def delete_parameter(self, name: str, roi: str):
        self.set_parameter(name, np.nan, roi=roi)


# ------------------------------------------------------------------
# Fitting helpers (ported from TraceIO, no tkinter dependency)
# ------------------------------------------------------------------

def r2_label(model_name: str) -> str:
    """Parameter name a model's fit-quality R² is stored under, e.g. `"R2 zWLC"`."""
    return f"R2 {model_name}"


def r2_full_label(model_name: str) -> str:
    """Parameter name a model's whole-trace R² is stored under, e.g. `"R2 full zWLC"`.

    Unlike `r2_label`, evaluated against every point in the trace rather than
    just the selected/masked points used to fit — see `apply_fit`.
    """
    return f"R2 full {model_name}"


def prepare_fit(func_name: str, initial_parameters=None):
    """Return (lmfit.Model, params, names) for a function in ForceSpectroscopy."""
    try:
        func = getattr(fs, func_name)
    except AttributeError:
        raise ValueError(f"Function '{func_name}' not found in ForceSpectroscopy")

    model = Model(func)
    params = model.make_params()
    for p in params:
        params[p].vary = False

    doc = func.__doc__ or ""
    names = {}
    try:
        names["x"] = doc.split("x:")[1].split("#")[0].strip()
        names["y"] = doc.split("return:")[1].split("#")[0].strip()
    except IndexError:
        raise ValueError(f"'{func_name}' docstring missing :param x: or :return: entries")

    for key in params:
        try:
            names[key] = doc.split(f"{key}:")[1].split("#")[0].strip()
        except IndexError:
            names[key] = key
        try:
            range_str = (
                doc.split(f"{key}:")[-1]
                .split("\n")[0]
                .split("]")[0]
                .split("[")[-1]
                .strip()
                .split(", ")
            )
            range_vals = np.asarray(range_str).astype(float)
            params[key].set(min=range_vals[0], max=range_vals[1])
        except (ValueError, IndexError):
            pass

    if initial_parameters is not None:
        params.update(initial_parameters)

    return model, params, names


def apply_fit(func_name: str, x, y, selection=None, initial_parameters=None):
    """Fit func_name to (x, y) using selected points. Returns lmfit ModelResult.

    `result.r2` is R² over just the selected/masked points (can look
    artificially low or unstable for a small or low-variance selection, since
    its `ss_tot` denominator shrinks with it). `result.r2_full` is the same
    formula evaluated against the full, unmasked x/y — more reliable for
    judging whether the fitted model generalizes, for models (e.g. zWLC)
    whose parameters are physically expected to hold across the whole trace,
    not just the fit selection.
    """
    def r2(y_obs, y_fit):
        ss_res = np.sum((y_obs - y_fit) ** 2)
        ss_tot = np.sum((y_obs - np.mean(y_obs)) ** 2)
        return 1 - ss_res / ss_tot if ss_tot != 0 else np.nan

    model, params, names = prepare_fit(func_name, initial_parameters)

    x_full = np.asarray(x, dtype=float)
    y_full = np.asarray(y, dtype=float)
    x_fit = x_full
    y_fit_data = y_full
    if selection is not None:
        mask = np.asarray(selection) == 1.0
        x_fit = x_fit[mask]
        y_fit_data = y_fit_data[mask]

    all_fixed = all(not p.vary for p in params.values())
    if all_fixed:
        y_pred = model.eval(params=params, x=x_fit)
        result = lmfit.model.ModelResult(model, params)
        result.best_fit = y_pred
        result.r2 = r2(y_fit_data, y_pred)
        result.success = True
    else:
        result = model.fit(y_fit_data, params, x=x_fit)
        result.r2 = r2(y_fit_data, result.best_fit)

    y_pred_full = model.eval(params=result.params, x=x_full)
    result.r2_full = r2(y_full, y_pred_full)

    return result, names
