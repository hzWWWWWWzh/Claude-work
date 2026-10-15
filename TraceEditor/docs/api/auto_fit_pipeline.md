# Table of Contents

* [trace\_editor.core.auto\_fit\_pipeline](#trace_editor.core.auto_fit_pipeline)
  * [load\_config](#trace_editor.core.auto_fit_pipeline.load_config)
  * [load\_pipeline](#trace_editor.core.auto_fit_pipeline.load_pipeline)
  * [list\_pipeline\_names](#trace_editor.core.auto_fit_pipeline.list_pipeline_names)
  * [find\_step\_for\_model](#trace_editor.core.auto_fit_pipeline.find_step_for_model)
  * [resolve\_params](#trace_editor.core.auto_fit_pipeline.resolve_params)
  * [resolve\_param\_sources](#trace_editor.core.auto_fit_pipeline.resolve_param_sources)
  * [step\_criteria](#trace_editor.core.auto_fit_pipeline.step_criteria)
  * [resolve\_criteria](#trace_editor.core.auto_fit_pipeline.resolve_criteria)
  * [resolve\_force\_defaults](#trace_editor.core.auto_fit_pipeline.resolve_force_defaults)
  * [run\_pipeline](#trace_editor.core.auto_fit_pipeline.run_pipeline)

<a id="trace_editor.core.auto_fit_pipeline"></a>

# trace\_editor.core.auto\_fit\_pipeline

Declarative fit pipelines: named, ordered sequences of fit steps read from
`app_config.toml` (Qt-free — see `batch_fit.py`'s module docstring for why).

A pipeline is a top-level TOML array-of-tables, each element one step:

    [[zwlc_pipeline]]
    model = "ExponentialDrift"          # a ForceSpectroscopy name, or "calculate_force"
    params."Z0 (um)" = {vary = true}    # per-param: any subset of vary/value/min/max
    [[zwlc_pipeline.criteria]]          # optional, AND-ed selection criteria
    channel = "Shift (mm)"
    min = "p90"

    [[zwlc_pipeline]]
    model = "calculate_force"
    params."L1 (mm)" = 1.4              # bare numbers -- calibration constants, no vary concept

See `app_config.default.toml` for a full worked example, and
`docs/architecture.md` for the fallback chain (`resolve_params`/
`resolve_criteria`) that replaced `[model_defaults]`/`[model_selection]`.

<a id="trace_editor.core.auto_fit_pipeline.load_config"></a>

#### load\_config

```python
def load_config() -> dict
```

Read the user's `app_config.toml` if present, else the checked-in default.

<a id="trace_editor.core.auto_fit_pipeline.load_pipeline"></a>

#### load\_pipeline

```python
def load_pipeline(config: dict, pipeline_name: str) -> list
```

Return `pipeline_name`'s ordered list of step dicts, or `[]` if absent.

<a id="trace_editor.core.auto_fit_pipeline.list_pipeline_names"></a>

#### list\_pipeline\_names

```python
def list_pipeline_names(config: dict) -> list
```

Return every top-level config key that looks like a pipeline.

A pipeline is a list of dicts each containing a `"model"` key -- this
distinguishes pipelines from unrelated top-level sections like
`[force_defaults]` (a single table, not a list).

<a id="trace_editor.core.auto_fit_pipeline.find_step_for_model"></a>

#### find\_step\_for\_model

```python
def find_step_for_model(steps: list, model_name: str)
```

Return the first step in `steps` whose `model` matches, or `None`.

A pipeline may list the same model more than once (e.g. a loose fit
followed by a tighter refit); `run_pipeline` runs every occurrence in
order regardless. This lookup only picks the first one, so anything that
resolves params/criteria by model name alone (the interactive Fit dialog,
`resolve_params`, `resolve_criteria`, `resolve_force_defaults`) only ever
sees the first occurrence's settings.

<a id="trace_editor.core.auto_fit_pipeline.resolve_params"></a>

#### resolve\_params

```python
def resolve_params(model_name: str, active_pipeline, config: dict)
```

Build `(model, params, names)` for `model_name`, seeded from

`prepare_fit`'s bare model-signature defaults, then `active_pipeline`'s
matching step (if any), which only overrides the specific
`vary`/`value`/`min`/`max` keys it actually sets. Pipelines are explicit:
there is no fallback to another pipeline for params a step doesn't list.

**Arguments**:

- `model_name`: a ForceSpectroscopy model name, e.g. `"zWLC"`.
- `active_pipeline`: pipeline name to read the matching step from, or `None`.
- `config`: parsed `app_config.toml`.

**Returns**:

`(model, params, names)`, same shape as `prepare_fit`.

<a id="trace_editor.core.auto_fit_pipeline.resolve_param_sources"></a>

#### resolve\_param\_sources

```python
def resolve_param_sources(model_name: str, active_pipeline,
                          config: dict) -> dict
```

Classify each of `model_name`'s parameters by where its dialog value

should come from, per `active_pipeline`'s matching step (mirrors
`resolve_params`'s layering; no fallback to another pipeline):

- `"Fit"`: the pipeline sets `vary = true` for this param (source: lmfit,
  starting from the pipeline's initial value if it also sets one).
- `"Global"`: the pipeline sets a fixed `value` (no `vary`) for this
  param (source: pipeline config, shared across every trace).
- `"Local"`: the pipeline doesn't mention this param at all (source:
  this ROI's own stored h5 value, falling back to the model's bare
  docstring default if nothing is stored yet).

**Arguments**:

- `model_name`: a ForceSpectroscopy model name, e.g. `"zWLC"`.
- `active_pipeline`: pipeline name to read the matching step from, or `None`.
- `config`: parsed `app_config.toml`.

**Returns**:

`{label: "Fit" | "Global" | "Local"}` for every non-x/y param.

<a id="trace_editor.core.auto_fit_pipeline.step_criteria"></a>

#### step\_criteria

```python
def step_criteria(step) -> list
```

Return one step's selection criteria list, or `[]` if it has none.

<a id="trace_editor.core.auto_fit_pipeline.resolve_criteria"></a>

#### resolve\_criteria

```python
def resolve_criteria(model_name: str, active_pipeline, config: dict) -> list
```

Return `active_pipeline`'s matching step's selection criteria, or `[]`
if the pipeline doesn't list criteria for this model (or no pipeline is
active). No fallback to another pipeline -- pipelines are explicit.

This is what `core.batch_fit._get_fit_selection`'s Default mode calls
instead of the old model-keyed `[[model_selection.<model>.criteria]]`.

<a id="trace_editor.core.auto_fit_pipeline.resolve_force_defaults"></a>

#### resolve\_force\_defaults

```python
def resolve_force_defaults(active_pipeline, config: dict) -> dict
```

Return calibration constants for `calculate_force`: `[force_defaults]`
merged under `active_pipeline`'s `calculate_force` step's `params`, if one
exists (the step's own values win). No fallback to another pipeline.

Shared by the interactive "Calculate Force" action and pipeline runs, so
both compute the same values from the same config.

<a id="trace_editor.core.auto_fit_pipeline.run_pipeline"></a>

#### run\_pipeline

```python
def run_pipeline(pipeline_name: str,
                 filepath: str,
                 config: dict,
                 progress_cb=None,
                 error_cb=None,
                 status_cb=None)
```

Run every step of `pipeline_name` against `filepath`, in order.

`"calculate_force"` steps call `calculate_force_for_file` directly (one
call, not per-ROI). Fit-model steps run every ROI via `run_batch`, which
writes results to disk before this function moves on to the next step, so
a later step can safely assume an earlier step's output (e.g. `F (pN)`)
is already persisted.

**Arguments**:

- `pipeline_name`: a pipeline name from `list_pipeline_names`.
- `filepath`: path to the `.h5` file to run the pipeline against.
- `config`: parsed `app_config.toml`.
- `progress_cb`: optional `(done: int, total: int) -> None`, forwarded
to `run_batch` for fit-model steps.
- `error_cb`: optional `(error) -> None`, forwarded to `run_batch`.
- `status_cb`: optional `(model_name: str, step: int, total_steps: int)
-> None`, called once as each step starts (before its work runs), so
a caller can show e.g. "file.h5: step 2/3 (zWLC)" instead of a bare
percentage that resets at every step/file boundary.

