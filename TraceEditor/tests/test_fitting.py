"""Unit tests for the real fitting models and the app's fit machinery.

These drive ForceSpectroscopy models (ExponentialDrift, zWLC, MagnetForce)
and the app's `apply_fit`/`prepare_fit` directly. The full multi-stage
open->force->drift->WLC pipeline is covered in test_workflow.py.
"""

import numpy as np
import pytest
from lmfit import Parameters

from trace_editor.core import force_spectroscopy as fs
from trace_editor.core.h5_data import apply_fit, prepare_fit


pytestmark = pytest.mark.fitting


class TestExponentialDriftModel:
    """The real fs.ExponentialDrift model."""

    def test_zero_rate_is_constant(self):
        t = np.linspace(0, 50, 200)
        z = np.asarray(fs.ExponentialDrift(t, Z0=10.0, dZ_dt=0.0))
        assert np.allclose(z, 10.0)

    def test_linear_when_tau_infinite(self):
        """With tau=inf the drift is linear in time at rate dZ_dt (nm/s)."""
        t = np.linspace(0, 50, 200)
        z = np.asarray(fs.ExponentialDrift(t, Z0=10.0, dZ_dt=-2.0, tau=np.inf))
        # slope in um/s = dZ_dt/1000; over 50 s -> -0.1 um total.
        assert z[0] == pytest.approx(10.0, abs=1e-6)
        assert (z[-1] - z[0]) == pytest.approx(-2.0 / 1000.0 * 50.0, rel=1e-3)

    def test_recover_params_via_apply_fit(self):
        """apply_fit recovers Z0 and dZ_dt from clean ExponentialDrift data."""
        t = np.linspace(0, 50, 500)
        z = np.asarray(fs.ExponentialDrift(t, Z0=12.0, dZ_dt=-1.5, tau=np.inf))
        init = Parameters()
        init.add("Z0", value=0.0, vary=True)
        init.add("dZ_dt", value=0.0, vary=True)
        init.add("tau", value=1e9, vary=False)
        result, _ = apply_fit("ExponentialDrift", t, z, initial_parameters=init)
        assert result.params["Z0"].value == pytest.approx(12.0, abs=1e-3)
        assert result.params["dZ_dt"].value == pytest.approx(-1.5, abs=1e-2)


class TestMagnetForce:
    """The real fs.MagnetForce shift->force model."""

    def test_force_decreases_with_shift(self):
        """Higher Shift (magnet retracted) gives lower force."""
        shift = np.linspace(1.0, 10.8, 100)
        force = np.asarray(fs.MagnetForce(shift, fmax=100.0))
        assert np.all(np.diff(force) < 0)

    def test_force_positive(self):
        shift = np.linspace(1.0, 10.8, 100)
        force = np.asarray(fs.MagnetForce(shift, fmax=100.0))
        assert np.all(force > 0)


class TestZWLCModel:
    """The real fs.zWLC extensible worm-like-chain model."""

    def test_extension_increases_with_force(self):
        f = np.linspace(0.05, 40, 300)
        z = np.asarray(fs.zWLC(f, L=10000, z_bead=0.0))
        assert np.all(np.diff(z) > 0)

    def test_zbead_offsets_extension(self):
        """z_bead shifts the whole curve by a constant."""
        f = np.linspace(0.05, 40, 100)
        z0 = np.asarray(fs.zWLC(f, L=10000, z_bead=0.0))
        z1 = np.asarray(fs.zWLC(f, L=10000, z_bead=1.5))
        assert np.allclose(z1 - z0, 1.5, atol=1e-6)

    def test_recover_L_via_apply_fit(self):
        """apply_fit recovers contour length from clean zWLC data."""
        f = np.linspace(0.05, 40, 400)
        z = np.asarray(fs.zWLC(f, L=8000, A=50, S=1000, z_bead=0.0))
        init = Parameters()
        init.add("L", value=5000.0, min=1000, max=20000, vary=True)
        init.add("A", value=50.0, vary=False)
        init.add("S", value=1000.0, vary=False)
        init.add("z_bead", value=0.0, vary=True)
        result, _ = apply_fit("zWLC", f, z, initial_parameters=init)
        assert result.params["L"].value == pytest.approx(8000.0, rel=1e-3)


class TestPrepareFit:
    """prepare_fit metadata extraction from model docstrings."""

    def test_names_from_docstring(self):
        _, _, names = prepare_fit("zWLC")
        assert names["x"] == "F (pN)"
        assert names["y"] == "Z (um)"

    def test_params_start_fixed(self):
        _, params, _ = prepare_fit("zWLC")
        assert all(not p.vary for p in params.values())

    def test_unknown_model_raises(self):
        with pytest.raises(ValueError):
            prepare_fit("NoSuchModel")


class TestFitEvaluationMode:
    """apply_fit in all-fixed (evaluation) mode returns R2 against the data."""

    def test_all_fixed_eval_r2_high_on_match(self):
        f = np.linspace(0.05, 40, 200)
        z = np.asarray(fs.zWLC(f, L=10000, A=50, S=1000, z_bead=0.0))
        # All params fixed at the true values -> perfect eval, R2 == 1.
        init = Parameters()
        init.add("L", value=10000.0, vary=False)
        init.add("A", value=50.0, vary=False)
        init.add("S", value=1000.0, vary=False)
        init.add("z_bead", value=0.0, vary=False)
        result, _ = apply_fit("zWLC", f, z, initial_parameters=init)
        assert result.r2 == pytest.approx(1.0, abs=1e-9)


class TestWholeTraceR2:
    """apply_fit's result.r2_full — R2 evaluated against the full trace,
    not just the selected/masked fit points."""

    def test_matches_r2_when_no_selection(self):
        """With no selection, the fit and whole-trace domains are identical."""
        f = np.linspace(0.05, 40, 300)
        z = np.asarray(fs.zWLC(f, L=8000, A=50, S=1000, z_bead=0.0))
        init = Parameters()
        init.add("L", value=5000.0, min=1000, max=20000, vary=True)
        init.add("A", value=50.0, vary=False)
        init.add("S", value=1000.0, vary=False)
        init.add("z_bead", value=0.0, vary=True)
        result, _ = apply_fit("zWLC", f, z, initial_parameters=init)
        assert result.r2_full == pytest.approx(result.r2, abs=1e-9)

    def test_stays_high_on_narrow_selection_when_model_still_fits(self):
        """A model that genuinely holds outside its fit selection (here,
        the true zWLC everywhere) keeps a high whole-trace R2 even when the
        selected-domain R2 degenerates on a narrow, low-variance selection."""
        f = np.linspace(0.05, 40, 400)
        z = np.asarray(fs.zWLC(f, L=8000, A=50, S=1000, z_bead=0.0))
        # Select only the first 5% (a narrow, near-constant low-force plateau).
        selection = np.zeros_like(f)
        selection[:20] = 1.0
        init = Parameters()
        init.add("L", value=8000.0, min=1000, max=20000, vary=True)
        init.add("A", value=50.0, vary=False)
        init.add("S", value=1000.0, vary=False)
        init.add("z_bead", value=0.0, vary=True)
        result, _ = apply_fit("zWLC", f, z, selection=selection, initial_parameters=init)
        assert result.r2_full == pytest.approx(1.0, abs=1e-6)

    def test_drops_when_fit_does_not_generalize(self):
        """A fit constrained only within its (misspecified) selection can
        look perfect there while failing badly on the rest of the trace --
        whole-trace R2 should catch what selected-domain R2 misses."""
        f = np.linspace(0.05, 40, 400)
        z = np.asarray(fs.zWLC(f, L=8000, A=50, S=1000, z_bead=0.0))
        # Selected points are shifted upward, only within the selection --
        # a fit can match this subset perfectly with a wrong z_bead, which
        # will not describe the (unmodified) rest of the trace.
        selection = np.zeros_like(f)
        selection[:20] = 1.0
        z_bad = z.copy()
        z_bad[:20] += 5.0
        init = Parameters()
        init.add("L", value=8000.0, vary=False)
        init.add("A", value=50.0, vary=False)
        init.add("S", value=1000.0, vary=False)
        init.add("z_bead", value=0.0, min=-10, max=10, vary=True)
        result, _ = apply_fit("zWLC", f, z_bad, selection=selection, initial_parameters=init)
        assert result.r2 == pytest.approx(1.0, abs=1e-6)
        assert result.r2_full < 0.5
