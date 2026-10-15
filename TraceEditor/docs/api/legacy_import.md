# Table of Contents

* [trace\_editor.core.legacy\_import](#trace_editor.core.legacy_import)
  * [import\_legacy\_hdf](#trace_editor.core.legacy_import.import_legacy_hdf)

<a id="trace_editor.core.legacy_import"></a>

# trace\_editor.core.legacy\_import

Import legacy `.hdf` (PyTables) files into the new `.h5` (h5py) format.

Reads with `tables` -- see `oldhdf_versions/TraceIO.py` for the original read
logic this mirrors. Only fields with a clear equivalent in the new
`/Parameters` schema (see `docs/file-format.md`) are carried over; anything
else (force-calibration constants, pixel-space bead centers, step-fit scratch
columns, ...) is dropped rather than invented as a new field.

<a id="trace_editor.core.legacy_import.import_legacy_hdf"></a>

#### import\_legacy\_hdf

```python
def import_legacy_hdf(hdf_path, h5_path=None) -> Path
```

Convert a legacy `.hdf` (PyTables) file to a new-format `.h5` file.

**Arguments**:

- `hdf_path`: path to the legacy `.hdf` file.
- `h5_path`: output path; defaults to `hdf_path` with suffix `.h5`.
Refuses to overwrite an existing file.

**Raises**:

- `FileExistsError`: if the output path already exists.

**Returns**:

the path written.

