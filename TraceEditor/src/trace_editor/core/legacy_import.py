"""
Import legacy `.hdf` (PyTables) files into the new `.h5` (h5py) format.

Reads with `tables` -- see `oldhdf_versions/TraceIO.py` for the original read
logic this mirrors. Only fields with a clear equivalent in the new
`/Parameters` schema (see `docs/file-format.md`) are carried over; anything
else (force-calibration constants, pixel-space bead centers, step-fit scratch
columns, ...) is dropped rather than invented as a new field.
"""

import re
import warnings
from pathlib import Path

import numpy as np
import h5py
import tables

# Legacy field/group names ("Shift (mm)", " parameters", ...) aren't valid
# Python identifiers, so PyTables warns that natural-naming access won't work
# -- irrelevant here since this module only reads via get_node()/subscript.
warnings.filterwarnings("ignore", category=tables.NaturalNameWarning)


# Legacy /parameters field name -> new /Parameters field name, for fields that
# map directly to a documented app-written field (see docs/file-format.md).
# Legacy pixel-space centers (X0/Y0 (pix)) and force-calibration constants
# (Fmax/Fmin/L1/L2/alpha/offset (pts)) have no equivalent slot and are dropped.
_PARAM_FIELD_MAP = {
    "Z0 (um)": "Z0 (um)",
    "dZ_dt (nm_s)": "dZ_dt (nm_s)",
    "tau (s)": "tau (s)",
    "dZ (um)": "dZ (um)",
    "X0 (um)": "X0 (um)",
    "Y0 (um)": "Y0 (um)",
    "Zbead (um)": "Zbead (um)",
    "L (bp)": "L (bp)",
    "A (nm)": "A (nm)",
    "S (pN)": "S (pN)",
}


def _natural_key(text):
    return [int(t) if t.isdigit() else t for t in re.split(r"(\d+)", text)]


def import_legacy_hdf(hdf_path, h5_path=None) -> Path:
    """Convert a legacy `.hdf` (PyTables) file to a new-format `.h5` file.

    :param hdf_path: path to the legacy `.hdf` file.
    :param h5_path: output path; defaults to `hdf_path` with suffix `.h5`.
        Refuses to overwrite an existing file.
    :return: the path written.
    :raises FileExistsError: if the output path already exists.
    """
    hdf_path = Path(hdf_path)
    h5_path = Path(h5_path) if h5_path is not None else hdf_path.with_suffix(".h5")
    if h5_path.exists():
        raise FileExistsError(f"{h5_path} already exists")

    with tables.open_file(hdf_path, mode="r") as src:
        bead_labels = sorted(
            (
                l for l in src.root._v_children.keys()
                if l not in (" parameters", " shared", " settings")
            ),
            key=_natural_key,
        )

        try:
            shared = src.get_node("/ shared")[:]
            params_table = src.get_node("/ parameters")[:]
        except tables.NoSuchNodeError as e:
            raise ValueError(
                f"{hdf_path.name} is missing required node {e.args[0] if e.args else '/ shared or / parameters'} "
                "-- not a valid legacy TraceEditor .hdf file"
            ) from e
        n_samples = len(shared)
        labels_col = [v.decode("utf-8") for v in params_table[" label"]]

        with h5py.File(h5_path, "w") as dst:
            dst.attrs["Date"] = ""
            dst.attrs["User"] = ""
            dst.attrs["Notes"] = f"Imported from legacy {hdf_path.name}"
            dst.attrs["axis"] = "Shift_mm"
            dst.attrs["num_rois"] = len(bead_labels)
            dst.attrs["roi_size"] = 0
            dst.attrs["bin_roi_size"] = 0

            param_dtype = [
                ("ROI", np.int32),
                ("X (pix)", np.float32),
                ("Y (pix)", np.float32),
                (" selected", np.float64),
            ]
            new_field_names = set()
            for label in bead_labels:
                if label not in labels_col:
                    continue
                row = params_table[labels_col.index(label)]
                for legacy_name, new_name in _PARAM_FIELD_MAP.items():
                    if legacy_name in row.dtype.names and not np.isnan(row[legacy_name]):
                        new_field_names.add(new_name)
                        err_name = f"e:{legacy_name}"
                        if err_name in row.dtype.names and not np.isnan(row[err_name]):
                            new_field_names.add(f"{new_name} Error")

            param_dtype += [(name, np.float64) for name in sorted(new_field_names)]
            params_arr = np.zeros(len(bead_labels), dtype=np.dtype(param_dtype))
            for name in new_field_names:
                params_arr[name] = np.nan

            traces_grp = dst.create_group("Traces")
            traces_grp.create_dataset(
                "Time (s)", data=shared["Time (s)"].astype(np.float64)
            )
            if "Shift (mm)" in shared.dtype.names:
                traces_grp.create_dataset(
                    "Shift (mm)", data=shared["Shift (mm)"].astype(np.float64)
                )
            if "F (pN)" in shared.dtype.names:
                traces_grp.create_dataset(
                    "F (pN)", data=shared["F (pN)"].astype(np.float64)
                )

            channel_groups = {}

            for roi_idx, label in enumerate(bead_labels):
                params_arr[roi_idx]["ROI"] = roi_idx
                params_arr[roi_idx]["X (pix)"] = np.nan
                params_arr[roi_idx]["Y (pix)"] = np.nan

                if label in labels_col:
                    row = params_table[labels_col.index(label)]
                    selected = row[" selected"] if " selected" in row.dtype.names else 1.0
                    params_arr[roi_idx][" selected"] = 1.0 if np.isnan(selected) else selected
                    for legacy_name, new_name in _PARAM_FIELD_MAP.items():
                        if legacy_name not in row.dtype.names:
                            continue
                        val = row[legacy_name]
                        if not np.isnan(val):
                            params_arr[roi_idx][new_name] = val
                        err_name = f"e:{legacy_name}"
                        if err_name in row.dtype.names and not np.isnan(row[err_name]):
                            params_arr[roi_idx][f"{new_name} Error"] = row[err_name]
                else:
                    params_arr[roi_idx][" selected"] = 1.0

                roi_key = f"{roi_idx:03d}"
                bead_table = src.get_node(f"/{label}")[:]
                for legacy_channel in ("X (um)", "Y (um)", "Z (um)", "A (a.u.)", "Selection"):
                    if legacy_channel not in bead_table.dtype.names:
                        continue
                    dtype = np.float32
                    data = bead_table[legacy_channel].astype(dtype)
                    if len(data) != n_samples:
                        continue
                    group = channel_groups.get(legacy_channel)
                    if group is None:
                        channel_name = "Selection" if legacy_channel == "Selection" else legacy_channel
                        group = traces_grp.create_group(channel_name)
                        channel_groups[legacy_channel] = group
                    group.create_dataset(roi_key, data=data)

            dst.create_dataset("Parameters", data=params_arr)

    return h5_path
