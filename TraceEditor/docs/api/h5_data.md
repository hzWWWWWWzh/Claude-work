# Table of Contents

* [trace\_editor.core.h5\_data](#trace_editor.core.h5_data)
  * [H5Data](#trace_editor.core.h5_data.H5Data)
    * [\_\_init\_\_](#trace_editor.core.h5_data.H5Data.__init__)
    * [is\_lut\_file](#trace_editor.core.h5_data.H5Data.is_lut_file)
    * [num\_rois](#trace_editor.core.h5_data.H5Data.num_rois)
    * [get\_axis](#trace_editor.core.h5_data.H5Data.get_axis)
    * [get\_user](#trace_editor.core.h5_data.H5Data.get_user)
    * [set\_notes](#trace_editor.core.h5_data.H5Data.set_notes)
    * [get\_lut\_filename](#trace_editor.core.h5_data.H5Data.get_lut_filename)
    * [selected\_count](#trace_editor.core.h5_data.H5Data.selected_count)
    * [list\_rois](#trace_editor.core.h5_data.H5Data.list_rois)
    * [get\_notes](#trace_editor.core.h5_data.H5Data.get_notes)
    * [get\_date](#trace_editor.core.h5_data.H5Data.get_date)
    * [shared\_channels](#trace_editor.core.h5_data.H5Data.shared_channels)
    * [roi\_channels](#trace_editor.core.h5_data.H5Data.roi_channels)
    * [all\_channels](#trace_editor.core.h5_data.H5Data.all_channels)
    * [get\_shared](#trace_editor.core.h5_data.H5Data.get_shared)
    * [get\_trace](#trace_editor.core.h5_data.H5Data.get_trace)
    * [get\_channel](#trace_editor.core.h5_data.H5Data.get_channel)
    * [get\_selection](#trace_editor.core.h5_data.H5Data.get_selection)
    * [set\_selection](#trace_editor.core.h5_data.H5Data.set_selection)
    * [get\_selection\_mask](#trace_editor.core.h5_data.H5Data.get_selection_mask)
    * [set\_selection\_mask](#trace_editor.core.h5_data.H5Data.set_selection_mask)
    * [get\_parameter](#trace_editor.core.h5_data.H5Data.get_parameter)
    * [set\_parameter](#trace_editor.core.h5_data.H5Data.set_parameter)
    * [set\_parameters\_bulk](#trace_editor.core.h5_data.H5Data.set_parameters_bulk)
    * [list\_parameters](#trace_editor.core.h5_data.H5Data.list_parameters)
  * [r2\_label](#trace_editor.core.h5_data.r2_label)
  * [r2\_full\_label](#trace_editor.core.h5_data.r2_full_label)
  * [prepare\_fit](#trace_editor.core.h5_data.prepare_fit)
  * [apply\_fit](#trace_editor.core.h5_data.apply_fit)

<a id="trace_editor.core.h5_data"></a>

# trace\_editor.core.h5\_data

Read/write the new h5 instrument format (h5py).

Parameters and selection are stored by extending the /Parameters compound
dataset with additional float64 fields, one row per ROI.

<a id="trace_editor.core.h5_data.H5Data"></a>

## H5Data Objects

```python
class H5Data()
```

Read/write access to one new-format `.h5` instrument file.

Wraps an on-disk h5py file. Trace data is read on demand; fit parameters and
selection are stored by extending the `/Parameters` compound dataset with
float64 fields, one row per ROI. There is no in-memory cache — every getter
reads from disk and every setter writes to disk, so instances stay
consistent across processes.

See `docs/file-format.md` for the on-disk layout.

<a id="trace_editor.core.h5_data.H5Data.__init__"></a>

#### \_\_init\_\_

```python
def __init__(filepath)
```

**Arguments**:

- `filepath`: path to an existing `.h5` file.

**Raises**:

- `FileNotFoundError`: if the file does not exist.

<a id="trace_editor.core.h5_data.H5Data.is_lut_file"></a>

#### is\_lut\_file

```python
def is_lut_file() -> bool
```

LUT files use a .LUT extension, never .h5. This always returns False.

<a id="trace_editor.core.h5_data.H5Data.num_rois"></a>

#### num\_rois

```python
def num_rois() -> int
```

Return the number of tracked beads (ROIs) in the file.

<a id="trace_editor.core.h5_data.H5Data.get_axis"></a>

#### get\_axis

```python
def get_axis() -> str
```

Return the `axis` root attribute (`"Shift_mm"` or `"Focus_mm"`), or empty string.

<a id="trace_editor.core.h5_data.H5Data.get_user"></a>

#### get\_user

```python
def get_user() -> str
```

Return the `User` root attribute, or empty string.

<a id="trace_editor.core.h5_data.H5Data.set_notes"></a>

#### set\_notes

```python
def set_notes(notes: str)
```

Overwrite the `Notes` root attribute.

<a id="trace_editor.core.h5_data.H5Data.get_lut_filename"></a>

#### get\_lut\_filename

```python
def get_lut_filename() -> str
```

Return the LUT file used to calculate bead height, as a bare
filename with no path or extension (e.g. `"data_013"`).

Read from `/Settings` (`Imaging`/`lut_filename`). Returns empty string
if `/Settings` is absent (e.g. legacy-imported files) or the key isn't
set.

<a id="trace_editor.core.h5_data.H5Data.selected_count"></a>

#### selected\_count

```python
def selected_count() -> int
```

Return the number of ROIs with ` selected` = 1.0 (defaults to all if unset).

<a id="trace_editor.core.h5_data.H5Data.list_rois"></a>

#### list\_rois

```python
def list_rois() -> list
```

Return ROI ids as zero-padded 3-digit strings: `["000", "001", ...]`.

<a id="trace_editor.core.h5_data.H5Data.get_notes"></a>

#### get\_notes

```python
def get_notes() -> str
```

Return Notes attribute from root group, or empty string.

<a id="trace_editor.core.h5_data.H5Data.get_date"></a>

#### get\_date

```python
def get_date() -> str
```

Return Date attribute from root group, or empty string.

<a id="trace_editor.core.h5_data.H5Data.shared_channels"></a>

#### shared\_channels

```python
def shared_channels() -> list
```

Channels shared across all ROIs (e.g. `Time (s)`, `Shift (mm)`, `F (pN)`).

These are datasets directly under `/Traces`.

<a id="trace_editor.core.h5_data.H5Data.roi_channels"></a>

#### roi\_channels

```python
def roi_channels() -> list
```

Per-ROI channels (e.g. `Z (um)`, `X (um)`, `A (a.u.)`).

These are groups under `/Traces`, each holding one dataset per ROI.
Excludes `Selection`, which is a point-selection mask, not a plottable
or fittable physical channel.

<a id="trace_editor.core.h5_data.H5Data.all_channels"></a>

#### all\_channels

```python
def all_channels() -> list
```

All channel names: shared channels followed by per-ROI channels.

<a id="trace_editor.core.h5_data.H5Data.get_shared"></a>

#### get\_shared

```python
def get_shared(channel: str) -> np.ndarray
```

Return a shared channel's data as float64.

**Arguments**:

- `channel`: a name from `shared_channels()`.

**Returns**:

1-D array, length = number of samples.

<a id="trace_editor.core.h5_data.H5Data.get_trace"></a>

#### get\_trace

```python
def get_trace(channel: str, roi: str) -> np.ndarray
```

Return one ROI's trace for a per-ROI channel, as float64.

**Arguments**:

- `channel`: a name from `roi_channels()`.
- `roi`: zero-padded ROI id, e.g. `"003"`.

**Returns**:

1-D array, length = number of samples.

<a id="trace_editor.core.h5_data.H5Data.get_channel"></a>

#### get\_channel

```python
def get_channel(channel: str, roi: str) -> np.ndarray
```

Return data for any channel — shared or per-ROI — dispatching on type.

**Arguments**:

- `channel`: any name from `all_channels()`.
- `roi`: zero-padded ROI id; ignored for shared channels.

**Returns**:

1-D float64 array.

<a id="trace_editor.core.h5_data.H5Data.get_selection"></a>

#### get\_selection

```python
def get_selection(roi: str) -> float
```

Return this ROI's selection flag.

**Arguments**:

- `roi`: zero-padded ROI id.

**Returns**:

1.0 if selected (included in batch ops), 0.0 if rejected.
Defaults to 1.0 when the ` selected` field has not been written yet.

<a id="trace_editor.core.h5_data.H5Data.set_selection"></a>

#### set\_selection

```python
def set_selection(roi: str, value: float)
```

Set this ROI's selection flag (1.0 = selected, 0.0 = rejected).

<a id="trace_editor.core.h5_data.H5Data.get_selection_mask"></a>

#### get\_selection\_mask

```python
def get_selection_mask(roi: str) -> np.ndarray | None
```

Return this ROI's per-point selection mask (1.0/0.0 float array).

Stored at `/Traces/Selection/<roi>`, matching the layout of the other
per-ROI channels. Returns None if no point selection has been made yet
for this ROI (all points should then be treated as selected).

<a id="trace_editor.core.h5_data.H5Data.set_selection_mask"></a>

#### set\_selection\_mask

```python
def set_selection_mask(roi: str, mask: np.ndarray)
```

Write this ROI's per-point selection mask, creating the dataset if needed.

**Arguments**:

- `mask`: 1-D array (1.0 = selected, 0.0 = rejected), length = n_samples.

<a id="trace_editor.core.h5_data.H5Data.get_parameter"></a>

#### get\_parameter

```python
def get_parameter(name: str, roi: str, default=np.nan) -> float
```

Return a stored parameter value for one ROI.

**Arguments**:

- `name`: parameter field name, e.g. `"L (bp)"`.
- `roi`: zero-padded ROI id.
- `default`: value returned when the field is absent or stored as NaN.

**Returns**:

the stored float, or `default`.

<a id="trace_editor.core.h5_data.H5Data.set_parameter"></a>

#### set\_parameter

```python
def set_parameter(name: str, value: float, roi: str)
```

Add or update a float parameter field for one ROI.

Rewrites the `/Parameters` dataset. If `name` is a new field it is added
as a float64 column, defaulting to NaN for every other ROI.

**Arguments**:

- `name`: parameter field name.
- `value`: value to store for this ROI.
- `roi`: zero-padded ROI id.

<a id="trace_editor.core.h5_data.H5Data.set_parameters_bulk"></a>

#### set\_parameters\_bulk

```python
def set_parameters_bulk(values_by_roi: dict)
```

Write many parameter fields for many ROIs in one open/read/write pass.

Equivalent to calling `set_parameter` once per (roi, field) pair, but
reads and rewrites the `/Parameters` array only once in total instead
of once per call -- `set_parameter` is O(rows) per call (it rewrites
the whole compound dataset each time), so writing N ROIs x M fields
one field at a time is O(N*M*N); this is O(N*M).  Intended for batch
fitting, where many ROIs' results need to land in one file at once.

**Arguments**:

- `values_by_roi`: `{roi: {field_name: value}}`. New field names
are added as float64 columns, defaulting to NaN for every ROI not
present in that field's `values_by_roi` entry.

<a id="trace_editor.core.h5_data.H5Data.list_parameters"></a>

#### list\_parameters

```python
def list_parameters(roi: str) -> dict
```

Return {name: value} for all analysis parameters of this ROI.

Excludes the instrument fields ROI, X (pix), Y (pix).

<a id="trace_editor.core.h5_data.r2_label"></a>

#### r2\_label

```python
def r2_label(model_name: str) -> str
```

Parameter name a model's fit-quality R² is stored under, e.g. `"R2 zWLC"`.

<a id="trace_editor.core.h5_data.r2_full_label"></a>

#### r2\_full\_label

```python
def r2_full_label(model_name: str) -> str
```

Parameter name a model's whole-trace R² is stored under, e.g. `"R2 full zWLC"`.

Unlike `r2_label`, evaluated against every point in the trace rather than
just the selected/masked points used to fit — see `apply_fit`.

<a id="trace_editor.core.h5_data.prepare_fit"></a>

#### prepare\_fit

```python
def prepare_fit(func_name: str, initial_parameters=None)
```

Return (lmfit.Model, params, names) for a function in ForceSpectroscopy.

<a id="trace_editor.core.h5_data.apply_fit"></a>

#### apply\_fit

```python
def apply_fit(func_name: str, x, y, selection=None, initial_parameters=None)
```

Fit func_name to (x, y) using selected points. Returns lmfit ModelResult.

`result.r2` is R² over just the selected/masked points (can look
artificially low or unstable for a small or low-variance selection, since
its `ss_tot` denominator shrinks with it). `result.r2_full` is the same
formula evaluated against the full, unmasked x/y — more reliable for
judging whether the fitted model generalizes, for models (e.g. zWLC)
whose parameters are physically expected to hold across the whole trace,
not just the fit selection.

