# Table of Contents

* [trace\_editor.core.file\_browser](#trace_editor.core.file_browser)
  * [extract\_tags](#trace_editor.core.file_browser.extract_tags)
  * [scan\_h5\_metadata](#trace_editor.core.file_browser.scan_h5_metadata)
  * [scan\_folder](#trace_editor.core.file_browser.scan_folder)

<a id="trace_editor.core.file_browser"></a>

# trace\_editor.core.file\_browser

Scan a folder of `.h5` files and extract cheap preview metadata for the
Open File browser dialog, without loading trace data.

<a id="trace_editor.core.file_browser.extract_tags"></a>

#### extract\_tags

```python
def extract_tags(notes: str) -> list
```

Return all ``tag`` tokens found in `notes`, lowercased, in order of

first appearance, with duplicates removed.

Free-form: any ``word`` in Notes becomes a tag — there is no fixed
vocabulary. E.g. `"Fail — `lowdensity` `retry`"` yields
`["lowdensity", "retry"]`.

**Arguments**:

- `notes`: free-text Notes attribute.

**Returns**:

list of tag strings, without the leading `#`.

<a id="trace_editor.core.file_browser.scan_h5_metadata"></a>

#### scan\_h5\_metadata

```python
def scan_h5_metadata(path) -> dict
```

Read root attributes and selection counts for one `.h5` file.

Cheap: only reads root attrs and the small `/Parameters` array, never
trace data.

**Arguments**:

- `path`: path to a `.h5` file.

**Returns**:

dict with keys `path`, `filename`, `axis`, `tags` (list),
`selected`, `num_rois`, `notes`, `date`, `user`, `lut_filename`, and
`error` (empty on success).

<a id="trace_editor.core.file_browser.scan_folder"></a>

#### scan\_folder

```python
def scan_folder(folder) -> list
```

Scan all `.h5` files directly inside `folder` (non-recursive).

**Arguments**:

- `folder`: directory to scan.

**Returns**:

list of metadata dicts from `scan_h5_metadata`, sorted by filename.

