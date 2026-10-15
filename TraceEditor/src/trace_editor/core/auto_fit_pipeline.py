"""Declarative fit pipelines: named, ordered sequences of fit steps read from
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
"""

import logging
import tomllib
from pathlib import Path

from .force_calc import calculate_force_for_file
from .h5_data import H5Data, prepare_fit

log = logging.getLogger(__name__)

_APP_CONFIG_FILE = Path(__file__).parent.parent / "gui" / "app_config.toml"
_APP_CONFIG_DEFAULT_FILE = Path(__file__).parent.parent / "gui" / "app_config.default.toml"

DEFAULT_PIPELINE_NAME = "default_pipeline"


def load_config() -> dict:
    """Read the user's `app_config.toml` if present, else the checked-in default."""
    config_file = _APP_CONFIG_FILE if _APP_CONFIG_FILE.exists() else _APP_CONFIG_DEFAULT_FILE
    with open(config_file, "rb") as f:
        return tomllib.load(f)


def load_pipeline(config: dict, pipeline_name: str) -> list:
    """Return `pipeline_name`'s ordered list of step dicts, or `[]` if absent."""
    return config.get(pipeline_name, [])


def list_pipeline_names(config: dict) -> list:
    """Return every top-level config key that looks like a pipeline.

    A pipeline is a list of dicts each containing a `"model"` key -- this
    distinguishes pipelines from unrelated top-level sections like
    `[force_defaults]` (a single table, not a list).
    """
    names = []
    for key, value in config.items():
        if isinstance(value, list) and all(isinstance(s, dict) and "model" in s for s in value):
            names.append(key)
    return names


def find_step_for_model(steps: list, model_name: str):
    """Return the first step in `steps` whose `model` matches, or `None`.

    A pipeline may list the same model more than once (e.g. a loose fit
    followed by a tighter refit); `run_pipeline` runs every occurrence in
    order regardless. This lookup only picks the first one, so anything that
    resolves params/criteria by model name alone (the interactive Fit dialog,
    `resolve_params`, `resolve_criteria`, `resolve_force_defaults`) only ever
    sees the first occurrence's settings.
    """
    for step in steps:
        if step.get("model") == model_name:
            return step
    return None


def _params_dict_for(pipeline_name, model_name: str, config: dict) -> dict:
    if pipeline_name is None:
        return {}
    steps = load_pipeline(config, pipeline_name)
    step = find_step_for_model(steps, model_name)
    return step.get("params", {}) if step is not None else {}


def resolve_params(model_name: str, active_pipeline, config: dict):
    """Build `(model, params, names)` for `model_name`, seeded from
    `prepare_fit`'s bare model-signature defaults, then `active_pipeline`'s
    matching step (if any), which only overrides the specific
    `vary`/`value`/`min`/`max` keys it actually sets. Pipelines are explicit:
    there is no fallback to another pipeline for params a step doesn't list.

    :param model_name: a ForceSpectroscopy model name, e.g. `"zWLC"`.
    :param active_pipeline: pipeline name to read the matching step from, or `None`.
    :param config: parsed `app_config.toml`.
    :return: `(model, params, names)`, same shape as `prepare_fit`.
    """
    model, params, names = prepare_fit(model_name)
    label_to_param = {label: p for p, label in names.items() if label not in ("x", "y")}

    for label, spec in _params_dict_for(active_pipeline, model_name, config).items():
        p = label_to_param.get(label)
        if p is None:
            log.warning("pipeline references unknown parameter %r for model %r", label, model_name)
            continue
        params[p].set(**spec)

    return model, params, names


def resolve_param_sources(model_name: str, active_pipeline, config: dict) -> dict:
    """Classify each of `model_name`'s parameters by where its dialog value
    should come from, per `active_pipeline`'s matching step (mirrors
    `resolve_params`'s layering; no fallback to another pipeline):

    - `"Fit"`: the pipeline sets `vary = true` for this param (source: lmfit,
      starting from the pipeline's initial value if it also sets one).
    - `"Global"`: the pipeline sets a fixed `value` (no `vary`) for this
      param (source: pipeline config, shared across every trace).
    - `"Local"`: the pipeline doesn't mention this param at all (source:
      this ROI's own stored h5 value, falling back to the model's bare
      docstring default if nothing is stored yet).

    :param model_name: a ForceSpectroscopy model name, e.g. `"zWLC"`.
    :param active_pipeline: pipeline name to read the matching step from, or `None`.
    :param config: parsed `app_config.toml`.
    :return: `{label: "Fit" | "Global" | "Local"}` for every non-x/y param.
    """
    _, _, names = prepare_fit(model_name)
    labels = [label for label in names.values() if label not in ("x", "y")]

    sources = {label: "Local" for label in labels}
    for label, spec in _params_dict_for(active_pipeline, model_name, config).items():
        if label not in sources:
            continue
        sources[label] = "Fit" if spec.get("vary") else "Global"
    return sources


def step_criteria(step) -> list:
    """Return one step's selection criteria list, or `[]` if it has none."""
    return step.get("criteria", []) if step is not None else []


def resolve_criteria(model_name: str, active_pipeline, config: dict) -> list:
    """Return `active_pipeline`'s matching step's selection criteria, or `[]`
    if the pipeline doesn't list criteria for this model (or no pipeline is
    active). No fallback to another pipeline -- pipelines are explicit.

    This is what `core.batch_fit._get_fit_selection`'s Default mode calls
    instead of the old model-keyed `[[model_selection.<model>.criteria]]`.
    """
    if active_pipeline is None:
        return []
    steps = load_pipeline(config, active_pipeline)
    return step_criteria(find_step_for_model(steps, model_name))


def resolve_force_defaults(active_pipeline, config: dict) -> dict:
    """Return calibration constants for `calculate_force`: `[force_defaults]`
    merged under `active_pipeline`'s `calculate_force` step's `params`, if one
    exists (the step's own values win). No fallback to another pipeline.

    Shared by the interactive "Calculate Force" action and pipeline runs, so
    both compute the same values from the same config.
    """
    defaults = dict(config.get("force_defaults", {}))
    if active_pipeline is not None:
        steps = load_pipeline(config, active_pipeline)
        step = find_step_for_model(steps, "calculate_force")
        if step is not None:
            defaults.update(step.get("params", {}))
    return defaults


def run_pipeline(pipeline_name: str, filepath: str, config: dict,
                  progress_cb=None, error_cb=None, status_cb=None):
    """Run every step of `pipeline_name` against `filepath`, in order.

    `"calculate_force"` steps call `calculate_force_for_file` directly (one
    call, not per-ROI). Fit-model steps run every ROI via `run_batch`, which
    writes results to disk before this function moves on to the next step, so
    a later step can safely assume an earlier step's output (e.g. `F (pN)`)
    is already persisted.

    :param pipeline_name: a pipeline name from `list_pipeline_names`.
    :param filepath: path to the `.h5` file to run the pipeline against.
    :param config: parsed `app_config.toml`.
    :param progress_cb: optional `(done: int, total: int) -> None`, forwarded
        to `run_batch` for fit-model steps.
    :param error_cb: optional `(error) -> None`, forwarded to `run_batch`.
    :param status_cb: optional `(model_name: str, step: int, total_steps: int)
        -> None`, called once as each step starts (before its work runs), so
        a caller can show e.g. "file.h5: step 2/3 (zWLC)" instead of a bare
        percentage that resets at every step/file boundary.
    """
    from .batch_fit import run_batch

    steps = load_pipeline(config, pipeline_name)
    for i, step in enumerate(steps):
        model_name = step["model"]
        if status_cb is not None:
            status_cb(model_name, i + 1, len(steps))

        if model_name == "calculate_force":
            data = H5Data(filepath)
            calculate_force_for_file(data, resolve_force_defaults(pipeline_name, config))
            continue

        _, params, names = resolve_params(model_name, pipeline_name, config)
        criteria = resolve_criteria(model_name, pipeline_name, config)
        tasks = [(filepath, roi) for roi in H5Data(filepath).list_rois()]
        run_batch(tasks, model_name, params, names,
                  progress_cb=progress_cb, error_cb=error_cb,
                  criteria=criteria)
