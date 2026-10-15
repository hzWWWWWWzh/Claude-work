# HDF5 File Format Specification

## New instrument format (`.h5`, h5py)

Written by the new magnetic tweezers instrument. Read with `h5py`.

### Root attributes

| Attribute      | Type    | Example                    | Notes                                      |
|----------------|---------|----------------------------|--------------------------------------------|
| `Date`         | str     | `"2026-06-05 13:56:30"`    |                                            |
| `User`         | str     | `"wang"`                   |                                            |
| `Notes`        | str     | `"Experiment with Construct A"`                       | see Notes tagging convention below        |
| `axis`         | str     | `"Shift_mm"` / `"Focus_mm"`| which motor was moved during the experiment — `"Focus_mm"` is a focus-change experiment, not a LUT calibration file (see "LUT files" below) |
| `num_rois`     | int64   | `46`                       | number of tracked beads                   |
| `roi_size`     | int64   | `64`                       | ROI size in pixels                        |
| `bin_roi_size` | int64   | `64`                       |                                            |

#### Notes tagging convention

The Open File browser (`FileBrowserDialog`, see
[architecture.md](architecture.md)) extracts tags from `Notes` as free-form
`#hashtag` tokens (`extract_tags()` in `core/file_browser.py`) — any `#word`
anywhere in `Notes` becomes a tag, case-insensitive, with no fixed
vocabulary. E.g. `"Fail — #lowdensity #retry"` yields tags `lowdensity` and
`retry`. Query with `tag = lowdensity` in the browser's query box (see
architecture.md) to filter by tag. Notes without any `#tag` just show no
tags for that file.

### `/Parameters`  (compound dataset)

Shape `(num_rois,)`. Written by instrument with pixel coordinates; **extended by the app** with float64 fields as analysis proceeds. `H5Data.set_parameter()` rewrites the dataset adding new columns as needed (all new fields default to `NaN` for unprocessed ROIs).

#### Instrument fields (always present)

| Field      | dtype   | Notes              |
|------------|---------|--------------------|
| `ROI`      | int32   | zero-based index   |
| `X (pix)`  | float32 | camera x-position  |
| `Y (pix)`  | float32 | camera y-position  |

#### App-written fields (added as analysis runs)

| Field                | dtype   | Written by              | Notes |
|----------------------|---------|-------------------------|-------|
| ` selected`          | float64 | selection tools         | 1.0 = included in batch ops, 0.0 = rejected |
| `F (pN) calculated`  | float64 | Calculate Force         | 1.0 once force has been computed |
| `Z0 (um)`            | float64 | ExponentialDrift fit    | drift baseline |
| `dZ_dt (nm_s)`       | float64 | ExponentialDrift fit    | drift rate |
| `tau (s)`            | float64 | ExponentialDrift fit    | drift time constant |
| `dZ (um)`            | float64 | ExponentialDrift fit    | full tether extension |
| `X0 (um)`            | float64 | ExponentialDrift fit    | bead x-centre |
| `Y0 (um)`            | float64 | ExponentialDrift fit    | bead y-centre |
| `Zbead (um)`         | float64 | zWLC / geometry fit     | bead attachment offset (typically negative) |
| `L (bp)`             | float64 | zWLC / zFiber fit       | DNA contour length |
| `A (nm)`             | float64 | zWLC fit                | persistence length |
| `S (pN)`             | float64 | zWLC fit                | stretch modulus |
| `R2 <model>`         | float64 | any fit                 | e.g. `R2 zWLC` — goodness of fit |
| `<param> Error`      | float64 | any fit                 | stderr of that parameter, e.g. `L (bp) Error`. Named with a trailing ` Error` so it sorts next to its value column in h5 viewers. |

The Global/Local/Fit source of each parameter is **not** stored — it is inferred
at runtime from whether a value and a matching `<param> Error` exist. (Older
files may contain `_src_<param>` fields from a previous scheme; the app ignores
them.)

### `/Settings`  (compound dataset)

Shape `(40,)`. Instrument settings at acquisition time.

| Field   | dtype  |
|---------|--------|
| `Group` | object (bytes) |
| `Key`   | object (bytes) |
| `Value` | object (bytes) |

### `/Traces/`  (group)

#### Shared (one array per file)

| Dataset       | dtype   | Shape        | Present when          |
|---------------|---------|--------------|------------------------|
| `Time (s)`    | float64 | (n_samples,) | always                 |
| `Shift (mm)`  | float64 | (n_samples,) | `axis = "Shift_mm"`   |
| `Focus (mm)`  | float64 | (n_samples,) | `axis = "Focus_mm"` (focus-change experiments) |
| `Intensity`   | float64 | (n_samples,) | always                 |

#### Per-ROI (one dataset per bead, named `"000"`, `"001"`, ...)

| Group         | dtype   | Shape        | Notes                      |
|---------------|---------|--------------|----------------------------|
| `X (um)/NNN`  | float32 | (n_samples,) | bead x-position            |
| `Y (um)/NNN`  | float32 | (n_samples,) | bead y-position            |
| `Z (um)/NNN`  | float32 | (n_samples,) | bead z-position (extension)|
| `A (a.u.)/NNN`| float32 | (n_samples,) | diffraction ring amplitude |
| `Selection/NNN` | float32 | (n_samples,) | per-point selection mask (1.0 = selected, 0.0 = rejected), written by the Select Points dialog. Absent until the user makes a point selection for that ROI. Not a plottable/fittable channel — excluded from `roi_channels()`/`all_channels()`. |

ROI index `NNN` is zero-padded to 3 digits (`"000"` … `"045"`).

### Focus-change experiments (`axis = "Focus_mm"`)

Files where `axis = "Focus_mm"` recorded a focus-change trajectory (the
`Focus_mm` motor was moved instead of `Shift_mm`) and have `Focus (mm)`
instead of `Shift (mm)`. This is a legitimate experiment type, not a
calibration artifact — nothing in the app skips or filters these files.

### LUT files (`.LUT`, separate format)

Actual lookup-table calibration files used by the instrument firmware to
convert diffraction-ring diameter to bead height. These use a **`.LUT`
extension** and an entirely different on-disk format — they are not `.h5`
files and never appear in the Open File browser (`FileBrowserDialog`) or
elsewhere in this app. Each experiment `.h5` file records which LUT file was
used to calculate bead height, in `/Settings` (`Imaging`/`lut_filename`),
exposed via `H5Data.get_lut_filename()` (returns the bare filename with no
path or extension, e.g. `data_013`) and shown as the **LUT** column in the
Open File browser — useful for spotting which experiments share the same
bead-height calibration.

---

## Legacy format (`.hdf`, PyTables)

Written by the old instrument and by TraceEditor. Read with `tables` (PyTables).

### Layout

```
/ (root)
├── /0, /1, /2, ...    per-bead tables
├── / shared           shared time-series (note leading space in node name)
└── / parameters       fit parameters + metadata (note leading space)
```

Node names with leading spaces (` shared`, ` parameters`) are intentional — PyTables requires `_f_get_child()` to access them.

### `/ shared`  (PyTables Table)

One row per time point. Columns written by instrument; `F (pN)` added by TraceEditor after force calibration.

| Column       | dtype   | Notes                                     |
|--------------|---------|-------------------------------------------|
| `Time (s)`   | float64 | acquisition time                          |
| `Shift (mm)` | float32 | magnet position                           |
| `F (pN)`     | float32 | force — **absent until recalculate_force runs** |

### `/0`, `/1`, ...  (PyTables Tables)

One table per bead, named by integer index (no zero-padding). One row per time point.

| Column      | dtype   | Notes                          |
|-------------|---------|--------------------------------|
| `A (a.u.)`  | float64 | diffraction ring amplitude     |
| `Selection` | float64 | 1.0 = selected, 0.0 = rejected |
| `X (um)`    | float64 | bead x-position                |
| `Y (um)`    | float64 | bead y-position                |
| `Z (um)`    | float64 | bead z-position (extension)    |

### `/ parameters`  (PyTables Table)

One row per bead plus a `"global"` row at index 0. Fit parameters are added as columns by TraceEditor; only ` label` and ` selected` are present in a freshly converted file.

| Column       | dtype   | Notes                              |
|--------------|---------|------------------------------------|
| ` label`     | bytes   | `b"global"`, `b"0"`, `b"1"`, ...  |
| ` selected`  | float64 | 1.0 = included in batch jobs       |
| *(fit cols)* | float64 | added by TraceEditor after fitting |

---

## Comparison

| Property              | New `.h5`                          | Legacy `.hdf`                        |
|-----------------------|------------------------------------|--------------------------------------|
| Library               | h5py                               | PyTables (tables)                    |
| Bead indexing         | `"000"`, `"001"` (3-digit padded)  | `"0"`, `"1"` (no padding)            |
| Data layout           | per-channel groups, one ds/bead    | per-bead tables with all channels    |
| Shared data           | flat datasets in `Traces/`         | compound table `/ shared`            |
| dtype (per-bead)      | float32                            | float64                              |
| dtype (shared)        | float64                            | float64 / float32 mixed              |
| Metadata              | root attrs + `/Parameters` + `/Settings` | `/ parameters` table only       |
| Force `F (pN)`        | shared dataset in `Traces/` after Calculate Force | in `/ shared` after calibration |
| Fit parameters        | float64 columns in `/Parameters`   | columns in `/ parameters`            |
| Selection per bead (whole trace) | ` selected` column in `/Parameters`| `Selection` column in bead table (dual-purpose) |
| Selection per point   | `Traces/Selection/NNN` dataset     | `Selection` column in bead table (dual-purpose) |
| Focus-change experiments | identified by `axis = "Focus_mm"` | N/A                                  |
| LUT calibration file used | `/Settings` `Imaging`/`lut_filename` (separate `.LUT` file, not `.h5`) | not stored |
| Camera ROI positions  | `/Parameters` dataset              | not stored                           |
| Instrument settings   | `/Settings` dataset                | not stored                           |

## Conversion / import

`convert_h5_to_hdf.py` (repo root) converts `.h5` → `.hdf` for use with the **old** tkinter TraceEditor. This is not needed for the new PySide6 app.

The new app's **File > Import Legacy .hdf…** action (`core/legacy_import.py`,
`import_legacy_hdf()`) reads a legacy `.hdf` file with PyTables (`tables`,
a base dependency) and writes a new-format `.h5` file alongside it (refuses to
overwrite an existing file). Beads `"0", "1", ...` (unpadded) become ROIs
`"000", "001", ...` (zero-padded, sorted naturally); `/ shared` and per-bead
`X/Y/Z (um)`, `A (a.u.)`, `Selection` map directly onto the equivalent new
channels. Only `/ parameters` fields with a documented equivalent in the
app-written fields table above are carried over (`Z0/X0/Y0/Zbead (um)`,
`dZ_dt (nm_s)`, `tau (s)`, `dZ (um)`, `L (bp)`, `A (nm)`, `S (pN)`, plus their
`e:<field>` errors renamed to `<field> Error`); pixel-space bead centers
(`X0/Y0 (pix)`) and force-calibration constants (`Fmax`, `Fmin`, `L1`, `L2`,
`alpha`, `offset (pts)`) have no new-format slot and are dropped rather than
invented as new fields.
