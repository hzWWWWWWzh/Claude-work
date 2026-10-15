# Table of Contents

* [trace\_editor.core.force\_calc](#trace_editor.core.force_calc)
  * [calculate\_force\_for\_file](#trace_editor.core.force_calc.calculate_force_for_file)

<a id="trace_editor.core.force_calc"></a>

# trace\_editor.core.force\_calc

Convert Shift (mm) -> F (pN) for one file (Qt-free).

Shared by `gui.main_window.TraceEditor.calculate_force` (per-file body) and
`core.auto_fit_pipeline`'s `"calculate_force"` pipeline step, so the logic
exists once regardless of whether it's triggered from the GUI's Calculate
Force menu action or a headless auto-fit pipeline.

<a id="trace_editor.core.force_calc.calculate_force_for_file"></a>

#### calculate\_force\_for\_file

```python
def calculate_force_for_file(data: H5Data, force_defaults: dict)
```

Compute `F (pN)` for every ROI in `data` via `fs.MagnetForce`.

Uses each ROI's stored calibration parameters (`L1 (mm)`, `L2 (mm)`,
`Fmin (pN)`, `Fmax (pN)`, `alpha (-)`, `offset (pts)`), falling back to
`force_defaults` for any that aren't stored yet, and seeds the fallback
value into the ROI so it's stored for next time. Writes the result as
the shared `F (pN)` dataset and marks every ROI with
`F (pN) calculated = 1.0`.

**Arguments**:

- `data`: file to compute Force for; must have a `Shift (mm)` shared
channel.
- `force_defaults`: `[force_defaults]` from `app_config.toml`.

**Raises**:

- `KeyError`: if `data` has no `Shift (mm)` channel.

