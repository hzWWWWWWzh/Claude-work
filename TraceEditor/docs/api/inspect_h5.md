# Table of Contents

* [trace\_editor.core.inspect\_h5](#trace_editor.core.inspect_h5)
  * [InspectH5Error](#trace_editor.core.inspect_h5.InspectH5Error)
  * [open\_in\_vscode\_h5web](#trace_editor.core.inspect_h5.open_in_vscode_h5web)

<a id="trace_editor.core.inspect_h5"></a>

# trace\_editor.core.inspect\_h5

Open an `.h5` file in the H5Web VS Code extension, for inspecting the raw
HDF5 structure (all groups/datasets/attributes, not just what TraceEditor's
UI exposes).

Qt-free so it's testable without a QApplication; `gui/main_window.py` wraps
this with the dialog shown on failure.

<a id="trace_editor.core.inspect_h5.InspectH5Error"></a>

## InspectH5Error Objects

```python
class InspectH5Error(Exception)
```

Raised when the file can't be opened in VS Code's H5Web viewer.

:ivar reason: `"no_code_cli"` or `"no_extension"`, for callers that want
    to tailor their message/fallback UI.

<a id="trace_editor.core.inspect_h5.open_in_vscode_h5web"></a>

#### open\_in\_vscode\_h5web

```python
def open_in_vscode_h5web(filepath: str)
```

Open `filepath` in VS Code's H5Web viewer via the `code` CLI.

**Arguments**:

- `filepath`: path to the `.h5` file to inspect.

**Raises**:

- `InspectH5Error`: if the `code` CLI isn't on PATH, or the H5Web
extension isn't installed.

