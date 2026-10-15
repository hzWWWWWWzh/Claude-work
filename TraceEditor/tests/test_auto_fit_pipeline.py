"""Unit tests for core.auto_fit_pipeline (Qt-free -- no QApplication/qtbot needed).

Exercises pipeline parsing, param/criteria resolution, and end-to-end
pipeline execution against real synthetic data.
"""

import numpy as np
import pytest

from trace_editor.core.auto_fit_pipeline import (
    DEFAULT_PIPELINE_NAME, find_step_for_model, list_pipeline_names,
    load_pipeline, resolve_criteria, resolve_params, run_pipeline, step_criteria,
)
from trace_editor.core.h5_data import H5Data

from conftest import write_workflow_h5, make_synthetic_workflow, TRUE_L_BP

pytestmark = pytest.mark.fitting


class TestLoadPipeline:
    def test_preserves_step_order(self):
        config = {
            "zwlc_pipeline": [
                {"model": "ExponentialDrift"},
                {"model": "calculate_force"},
                {"model": "zWLC"},
            ],
        }
        steps = load_pipeline(config, "zwlc_pipeline")
        assert [s["model"] for s in steps] == ["ExponentialDrift", "calculate_force", "zWLC"]

    def test_missing_pipeline_returns_empty(self):
        assert load_pipeline({}, "no_such_pipeline") == []


class TestListPipelineNames:
    def test_finds_pipelines_ignores_other_sections(self):
        config = {
            "force_defaults": {"L1 (mm)": 1.4},
            "zwlc_pipeline": [{"model": "zWLC"}],
            "default_pipeline": [{"model": "zWLC"}],
        }
        assert set(list_pipeline_names(config)) == {"zwlc_pipeline", "default_pipeline"}

    def test_empty_config_has_no_pipelines(self):
        assert list_pipeline_names({}) == []


class TestFindStepForModel:
    def test_finds_first_match(self):
        steps = [{"model": "ExponentialDrift"}, {"model": "zWLC"}]
        assert find_step_for_model(steps, "zWLC")["model"] == "zWLC"

    def test_no_match_returns_none(self):
        steps = [{"model": "ExponentialDrift"}]
        assert find_step_for_model(steps, "zWLC") is None


class TestStepCriteria:
    def test_returns_criteria_list(self):
        step = {"model": "ExponentialDrift", "criteria": [{"channel": "Shift (mm)", "min": "p90"}]}
        assert step_criteria(step) == [{"channel": "Shift (mm)", "min": "p90"}]

    def test_no_criteria_returns_empty(self):
        assert step_criteria({"model": "zWLC"}) == []

    def test_none_step_returns_empty(self):
        assert step_criteria(None) == []


class TestResolveParams:
    def test_active_pipeline_step_sets_params(self):
        config = {
            "custom_pipeline": [{"model": "zWLC", "params": {"L (bp)": {"value": 9000}}}],
        }
        _, params, names = resolve_params("zWLC", "custom_pipeline", config)
        label_to_param = {v: k for k, v in names.items() if v not in ("x", "y")}
        assert params[label_to_param["L (bp)"]].value == 9000

    def test_other_pipelines_dont_leak_in(self):
        """A step for `zWLC` in one pipeline must not affect resolving `zWLC`
        against a different pipeline that doesn't itself set that param --
        pipelines are explicit, with no fallback between them."""
        config = {
            "other_pipeline": [{"model": "zWLC", "params": {"L (bp)": {"value": 16500}}}],
            "custom_pipeline": [{"model": "ExponentialDrift"}],
        }
        _, params, names = resolve_params("zWLC", "custom_pipeline", config)
        label_to_param = {v: k for k, v in names.items() if v not in ("x", "y")}
        assert params[label_to_param["L (bp)"]].value == 1  # zWLC's signature default, untouched

    def test_falls_back_to_bare_model_default_when_no_pipeline_matches(self):
        _, params, names = resolve_params("zWLC", None, {})
        label_to_param = {v: k for k, v in names.items() if v not in ("x", "y")}
        assert params[label_to_param["L (bp)"]].value == 1  # zWLC's signature default
        assert params[label_to_param["L (bp)"]].vary is False

    def test_vary_true_sets_vary_flag(self):
        config = {"custom_pipeline": [{"model": "zWLC", "params": {"Zbead (um)": {"vary": True}}}]}
        _, params, names = resolve_params("zWLC", "custom_pipeline", config)
        label_to_param = {v: k for k, v in names.items() if v not in ("x", "y")}
        assert params[label_to_param["Zbead (um)"]].vary is True

    def test_unknown_param_label_is_ignored_not_raised(self):
        config = {"custom_pipeline": [{"model": "zWLC", "params": {"no such param": {"vary": True}}}]}
        resolve_params("zWLC", "custom_pipeline", config)  # must not raise


class TestResolveCriteria:
    def test_active_pipeline_criteria_used_when_present(self):
        config = {
            "custom_pipeline": [
                {"model": "ExponentialDrift", "criteria": [{"channel": "Shift (mm)", "min": "p90"}]},
            ],
        }
        assert resolve_criteria("ExponentialDrift", "custom_pipeline", config) == \
            [{"channel": "Shift (mm)", "min": "p90"}]

    def test_other_pipelines_dont_leak_in(self):
        config = {
            DEFAULT_PIPELINE_NAME: [
                {"model": "ExponentialDrift", "criteria": [{"channel": "Shift (mm)", "min": "p80"}]},
            ],
            "custom_pipeline": [{"model": "ExponentialDrift"}],  # no criteria of its own
        }
        assert resolve_criteria("ExponentialDrift", "custom_pipeline", config) == []

    def test_no_match_anywhere_returns_empty(self):
        assert resolve_criteria("zWLC", "custom_pipeline", {}) == []

    def test_none_active_pipeline_returns_empty(self):
        config = {DEFAULT_PIPELINE_NAME: [
            {"model": "zWLC", "criteria": [{"channel": "F (pN)", "min": 5}]},
        ]}
        assert resolve_criteria("zWLC", None, config) == []


class TestRunPipeline:
    def test_drift_then_force_then_zwlc_end_to_end(self, tmp_path):
        """A 3-step pipeline mirrors the app's real workflow: fit drift over
        low-force points, calculate Force, then fit zWLC with L fixed and
        Zbead free -- asserting each step's output is available to the next.
        """
        beads = [make_synthetic_workflow(noise_um=0.0, seed=0, Z0=12.0, dZ_dt=-1.5)]
        path = write_workflow_h5(tmp_path / "pipeline.h5", beads)

        config = {
            "zwlc_pipeline": [
                {
                    "model": "ExponentialDrift",
                    "params": {"Z0 (um)": {"vary": True}, "dZ_dt (nm_s)": {"vary": True}},
                    "criteria": [{"channel": "Shift (mm)", "min": "p90"}],
                },
                {"model": "calculate_force", "params": {
                    "L1 (mm)": 1.4, "L2 (mm)": 0.8, "Fmin (pN)": 0.02,
                    "Fmax (pN)": 85.0, "alpha (-)": 0.7, "offset (pts)": 0.0,
                }},
                {
                    "model": "zWLC",
                    "params": {"L (bp)": {"value": TRUE_L_BP}, "Zbead (um)": {"vary": True}},
                },
            ],
        }

        errors = []
        run_pipeline("zwlc_pipeline", str(path), config, error_cb=errors.append)
        assert errors == []

        data = H5Data(str(path))
        params = data.list_parameters("000")
        assert params["F (pN) calculated"] == 1.0
        assert not np.isnan(params["dZ_dt (nm_s)"])
        assert params["R2 zWLC"] > 0.9

    def test_unknown_model_raises(self, tmp_path):
        beads = [make_synthetic_workflow(noise_um=0.0, seed=0)]
        path = write_workflow_h5(tmp_path / "bad.h5", beads)
        config = {"bad_pipeline": [{"model": "NotARealModel"}]}
        with pytest.raises(ValueError):
            run_pipeline("bad_pipeline", str(path), config)

    def test_missing_pipeline_is_a_noop(self, tmp_path):
        beads = [make_synthetic_workflow(noise_um=0.0, seed=0)]
        path = write_workflow_h5(tmp_path / "noop.h5", beads)
        run_pipeline("no_such_pipeline", str(path), {})  # must not raise
