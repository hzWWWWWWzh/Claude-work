"""Full analysis-workflow integration tests.

Ground-truth data is built with the real ForceSpectroscopy models (see
`conftest.make_synthetic_workflow`):

    Shift sweep -> MagnetForce -> zWLC(L=10000, z_bead=0) -> + ExponentialDrift

The workflow under test then tries to invert this:

    1. Fit ExponentialDrift over low-force points and subtract it. This also
       removes the (non-zero) WLC height at low force.
    2. Fit zWLC with free L and z_bead, which restores the WLC shape/height.

We test a noiseless version (tight tolerances — proving the logic is exact)
and a small-noise version (looser tolerances — proving robustness).

Note on the offset: a constant vertical shift is genuinely degenerate between
the drift Z0 and z_bead, so we assert on L, on the WLC *shape* (offset-removed
RMS), on the WLC *height* (peak-to-peak, offset-independent), and on the
z_bead prediction invariant — not on the absolute z_bead value.
"""

import numpy as np
import pytest

from trace_editor.core import force_spectroscopy as fs
from trace_editor.core.h5_data import H5Data
from workflow_helpers import run_workflow, low_force_selection, A_FIXED_NM, S_FIXED_PN
from conftest import make_synthetic_workflow, TRUE_L_BP


pytestmark = pytest.mark.workflow


class TestWorkflowNoiseless:
    """The workflow logic should be near-exact with no measurement noise."""

    @pytest.fixture(scope="class")
    def result(self):
        bead = make_synthetic_workflow(noise_um=0.0, seed=0)
        return bead, run_workflow(bead)

    def test_contour_length_recovered(self, result):
        """Fitted L returns to the ground-truth 10000 bp."""
        _, out = result
        assert abs(out["L_fit"] - TRUE_L_BP) < 50.0  # <0.5%

    def test_wlc_shape_restored(self, result):
        """Drift removal + zWLC refit restores the WLC shape (offset-removed)."""
        bead, out = result
        z_wlc_true = bead["z_wlc"]
        z_wlc_fit = out["z_wlc_fit"]
        # Remove the degenerate constant offset before comparing shape.
        offset = np.mean(z_wlc_fit - z_wlc_true)
        rms = np.sqrt(np.mean(((z_wlc_fit - offset) - z_wlc_true) ** 2))
        assert rms < 0.05  # um

    def test_wlc_height_restored(self, result):
        """The peak-to-peak WLC height is restored (offset-independent)."""
        bead, out = result
        h_true = float(np.ptp(bead["z_wlc"]))
        h_fit = float(np.ptp(out["z_wlc_fit"]))
        assert abs(h_fit - h_true) / h_true < 0.15

    def test_zbead_matches_lowforce_wlc_prediction(self, result):
        """Fitted z_bead == -mean(zWLC over selected low-force points).

        This pins down *why* the offset is what it is: the drift fit absorbs the
        mean WLC height at low force, and z_bead cancels exactly that.
        """
        _, out = result
        assert abs(out["z_bead_fit"] - out["z_bead_predicted"]) < 0.01

    def test_wlc_fit_quality(self, result):
        """zWLC fit is excellent on clean data."""
        _, out = result
        assert out["r2_wlc"] > 0.99

    def test_drift_removed_flattens_lowforce(self, result):
        """After drift removal the low-force plateau sits near the WLC low-F value."""
        bead, out = result
        sel = out["drift_selection"]
        # Corrected Z over selected points should be small (drift + Z0 removed),
        # on the order of the low-force WLC extension (~1-2 um), not ~12 um.
        assert np.mean(out["z_corrected"][sel]) < 3.0


class TestWorkflowNoisy:
    """The workflow should still recover parameters with small noise."""

    @pytest.fixture(scope="class")
    def result(self):
        bead = make_synthetic_workflow(noise_um=0.05, seed=1)
        return bead, run_workflow(bead)

    def test_contour_length_recovered(self, result):
        _, out = result
        assert abs(out["L_fit"] - TRUE_L_BP) < 300.0  # <3%

    def test_wlc_shape_restored(self, result):
        bead, out = result
        z_wlc_true = bead["z_wlc"]
        z_wlc_fit = out["z_wlc_fit"]
        offset = np.mean(z_wlc_fit - z_wlc_true)
        rms = np.sqrt(np.mean(((z_wlc_fit - offset) - z_wlc_true) ** 2))
        assert rms < 0.1  # um

    def test_wlc_height_restored(self, result):
        bead, out = result
        h_true = float(np.ptp(bead["z_wlc"]))
        h_fit = float(np.ptp(out["z_wlc_fit"]))
        assert abs(h_fit - h_true) / h_true < 0.2

    def test_zbead_matches_lowforce_wlc_prediction(self, result):
        _, out = result
        assert abs(out["z_bead_fit"] - out["z_bead_predicted"]) < 0.02

    def test_wlc_fit_quality(self, result):
        _, out = result
        assert out["r2_wlc"] > 0.95


class TestWorkflowReproducibility:
    """Running the workflow twice on identical data gives identical fits."""

    def test_deterministic(self):
        bead = make_synthetic_workflow(noise_um=0.05, seed=7)
        out1 = run_workflow(bead)
        out2 = run_workflow(bead)
        assert abs(out1["L_fit"] - out2["L_fit"]) < 1e-6
        assert abs(out1["z_bead_fit"] - out2["z_bead_fit"]) < 1e-9


class TestGeneratedFileIsRaw:
    """The permanent and regenerated files hold raw traces only — no fits."""

    def test_permanent_file_has_no_fits(self, permanent_workflow_h5):
        h5 = H5Data(permanent_workflow_h5)
        assert h5.num_rois() == 2
        # Only base instrument/selection fields; no fit columns yet.
        params = h5.list_parameters("000")
        assert "L (bp)" not in params
        assert "Z0 (um)" not in params
        # Raw traces are present.
        assert "Z (um)" in h5.roi_channels()
        assert "Shift (mm)" in h5.shared_channels()

    def test_generated_file_matches_permanent(self, generated_workflow_h5,
                                              permanent_workflow_h5):
        """Regenerated data reproduces the committed file's traces exactly."""
        gen = H5Data(generated_workflow_h5)
        perm = H5Data(permanent_workflow_h5)
        for roi in perm.list_rois():
            assert np.allclose(gen.get_channel("Z (um)", roi),
                               perm.get_channel("Z (um)", roi))


class TestWorkflowH5Roundtrip:
    """Run the workflow through a real .h5 file and write fits into a 3rd file.

    Reads raw traces from the generated file, runs the open->force->drift->WLC
    workflow, and writes fits + F(pN) into synthetic_workflow_fitted.h5. That
    file is left on disk (git-ignored) so it can be opened in the app to
    visually inspect the drift and zWLC fits.
    """

    def test_fit_and_write_back_to_h5(self, fitted_workflow_h5):
        import h5py
        h5 = H5Data(fitted_workflow_h5)
        assert h5.num_rois() == 2

        # Write the shared F (pN) dataset like calculate_force does.
        shift = h5.get_shared("Shift (mm)")
        force_shared = np.asarray(fs.MagnetForce(shift, fmax=100.0))
        with h5py.File(fitted_workflow_h5, "r+") as f:
            grp = f["Traces"]
            if "F (pN)" in grp:
                del grp["F (pN)"]
            grp.create_dataset("F (pN)", data=force_shared.astype(np.float64))

        for roi in h5.list_rois():
            # Read the traces back out of the file (exercises real file I/O).
            time = h5.get_channel("Time (s)", roi)
            shift = h5.get_channel("Shift (mm)", roi)
            z = h5.get_channel("Z (um)", roi)
            force = np.asarray(fs.MagnetForce(shift, fmax=100.0))

            bead = {"time": time, "shift": shift, "force": force, "z": z}
            out = run_workflow(bead)

            # Write fitted parameters back into /Parameters like the app does.
            h5.set_parameter("Z0 (um)", out["Z0_fit"], roi=roi)
            h5.set_parameter("dZ_dt (nm_s)", out["dZ_dt_fit"], roi=roi)
            h5.set_parameter("tau (s)", 1e9, roi=roi)
            h5.set_parameter("L (bp)", out["L_fit"], roi=roi)
            h5.set_parameter("A (nm)", A_FIXED_NM, roi=roi)
            h5.set_parameter("S (pN)", S_FIXED_PN, roi=roi)
            h5.set_parameter("Zbead (um)", out["z_bead_fit"], roi=roi)
            h5.set_parameter("R2 zWLC", out["r2_wlc"], roi=roi)
            h5.set_parameter("F (pN) calculated", 1.0, roi=roi)

        # Verify the written parameters persist and are sane.
        h5_reopened = H5Data(fitted_workflow_h5)
        assert "F (pN)" in h5_reopened.shared_channels()
        for roi in h5_reopened.list_rois():
            L = h5_reopened.get_parameter("L (bp)", roi)
            assert abs(L - TRUE_L_BP) < 500.0
            r2 = h5_reopened.get_parameter("R2 zWLC", roi)
            assert r2 > 0.9


class TestGroundTruthConstruction:
    """Sanity-check the synthetic ground truth itself (guards the fixtures)."""

    def test_force_spans_low_to_high(self):
        bead = make_synthetic_workflow(noise_um=0.0, seed=0)
        assert bead["force"].min() < 0.2      # near-zero force at high shift
        assert bead["force"].max() > 20.0     # high force at low shift

    def test_wlc_monotonic_in_force(self):
        """zWLC extension increases with force."""
        f = np.linspace(0.05, 40, 200)
        z = np.asarray(fs.zWLC(f, L=TRUE_L_BP, z_bead=0.0))
        assert np.all(np.diff(z) > 0)

    def test_drift_dominates_observed_low_force(self):
        """Observed Z at low force is ~Z0 (drift), far above the WLC height."""
        bead = make_synthetic_workflow(noise_um=0.0, seed=0)
        sel = bead["shift"] >= np.percentile(bead["shift"], 90)
        assert np.mean(bead["z"][sel]) > 10.0        # dominated by drift Z0
        assert np.mean(bead["z_wlc"][sel]) < 3.0     # true WLC height is small
