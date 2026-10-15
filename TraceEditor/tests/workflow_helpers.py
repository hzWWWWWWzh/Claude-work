"""Shared helpers driving the real analysis workflow via ForceSpectroscopy.

These wrap the app's `apply_fit`/`prepare_fit` so tests exercise the same
fitting code the GUI does. The two stages mirror the app:

  1. `fit_exponential_drift` — fit ExponentialDrift to Z(t) using only the
     low-force (high-Shift) points, matching `_get_fit_selection` in the GUI.
  2. `fit_zwlc` — fit zWLC to the drift-corrected Z vs force, with L and
     z_bead free (A, S fixed at their defaults).

`run_workflow` chains them and returns everything needed to assert on the
recovered parameters and reconstructed WLC curve.
"""

import numpy as np
from lmfit import Parameters

from trace_editor.core.h5_data import apply_fit
from trace_editor.core import force_spectroscopy as fs


# Fixed WLC nuisance parameters (kept at zWLC defaults, matching the app's
# typical zWLC fit where only L and z_bead vary).
A_FIXED_NM = 50.0
S_FIXED_PN = 1000.0


def low_force_selection(shift: np.ndarray, percentile: float = 90.0) -> np.ndarray:
    """Boolean mask of low-force points, mirroring GUI `_get_fit_selection`.

    High Shift = magnet retracted = LOW force, so we keep Shift >= 90th pct.
    (The GUI additionally applies an amplitude filter; our synthetic amplitude
    is high and flat, so that filter is a no-op here and is omitted.)
    """
    return shift >= np.percentile(shift, percentile)


def fit_exponential_drift(time, z, shift, tau_fixed=True):
    """Stage 1: fit ExponentialDrift to Z(t) over low-force points.

    Returns (drift_curve_full, result). `drift_curve_full` is the fitted drift
    evaluated over *all* time points (for subtraction), not just the selection.
    """
    sel = low_force_selection(shift)
    init = Parameters()
    init.add("Z0", value=float(np.mean(z[sel])), vary=True)
    init.add("dZ_dt", value=0.0, vary=True)
    # tau fixed at ~inf (linear drift) matches the common real-data case.
    init.add("tau", value=1e9, vary=not tau_fixed, min=1e-3)

    result, _ = apply_fit("ExponentialDrift", time, z,
                          selection=sel, initial_parameters=init)

    drift_full = np.asarray(fs.ExponentialDrift(
        time,
        Z0=result.params["Z0"].value,
        dZ_dt=result.params["dZ_dt"].value,
        tau=result.params["tau"].value,
    ))
    return drift_full, result, sel


def fit_zwlc(force, z_corrected, L_init=5000.0):
    """Stage 2: fit zWLC (free L, z_bead) to drift-corrected Z vs force."""
    init = Parameters()
    init.add("L", value=L_init, min=1000.0, max=20000.0, vary=True)
    init.add("A", value=A_FIXED_NM, vary=False)
    init.add("S", value=S_FIXED_PN, vary=False)
    init.add("z_bead", value=0.0, min=-5.0, max=5.0, vary=True)

    result, _ = apply_fit("zWLC", force, z_corrected,
                          selection=None, initial_parameters=init)
    return result


def run_workflow(bead: dict) -> dict:
    """Run the full open→force→drift→WLC workflow on one synthetic bead.

    `bead` is a dict from `conftest.make_synthetic_workflow`. Returns fitted
    parameters, the reconstructed WLC curve, and the drift selection so tests
    can assert on parameter recovery and on the z_bead prediction invariant.
    """
    time, shift, force, z = bead["time"], bead["shift"], bead["force"], bead["z"]

    drift_full, drift_res, sel = fit_exponential_drift(time, z, shift)
    z_corr = z - drift_full

    wlc_res = fit_zwlc(force, z_corr)
    L_fit = wlc_res.params["L"].value
    z_bead_fit = wlc_res.params["z_bead"].value

    z_wlc_fit = np.asarray(fs.zWLC(force, L=L_fit, A=A_FIXED_NM,
                                   S=S_FIXED_PN, z_bead=z_bead_fit))

    # Physical invariant: the drift fit's Z0 absorbs the mean WLC height over
    # the selected low-force points, so z_bead must cancel exactly that:
    #   z_bead ≈ -mean( zWLC(F_selected; L_fit, z_bead=0) )
    z_bead_predicted = -float(np.mean(
        fs.zWLC(force[sel], L=L_fit, A=A_FIXED_NM, S=S_FIXED_PN, z_bead=0.0)
    ))

    return {
        "drift_curve": drift_full,
        "z_corrected": z_corr,
        "z_wlc_fit": z_wlc_fit,
        "L_fit": L_fit,
        "z_bead_fit": z_bead_fit,
        "z_bead_predicted": z_bead_predicted,
        "drift_selection": sel,
        "drift_result": drift_res,
        "wlc_result": wlc_res,
        "Z0_fit": drift_res.params["Z0"].value,
        "dZ_dt_fit": drift_res.params["dZ_dt"].value,
        "r2_wlc": wlc_res.r2,
    }
