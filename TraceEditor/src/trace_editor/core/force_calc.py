"""Convert Shift (mm) -> F (pN) for one file (Qt-free).

Shared by `gui.main_window.TraceEditor.calculate_force` (per-file body) and
`core.auto_fit_pipeline`'s `"calculate_force"` pipeline step, so the logic
exists once regardless of whether it's triggered from the GUI's Calculate
Force menu action or a headless auto-fit pipeline.
"""

import math

import h5py
import numpy as np

from . import force_spectroscopy as fs
from .h5_data import H5Data


def calculate_force_for_file(data: H5Data, force_defaults: dict):
    """Compute `F (pN)` for every ROI in `data` via `fs.MagnetForce`.

    Uses each ROI's stored calibration parameters (`L1 (mm)`, `L2 (mm)`,
    `Fmin (pN)`, `Fmax (pN)`, `alpha (-)`, `offset (pts)`), falling back to
    `force_defaults` for any that aren't stored yet, and seeds the fallback
    value into the ROI so it's stored for next time. Writes the result as
    the shared `F (pN)` dataset and marks every ROI with
    `F (pN) calculated = 1.0`.

    :param data: file to compute Force for; must have a `Shift (mm)` shared
        channel.
    :param force_defaults: `[force_defaults]` from `app_config.toml`.
    :raises KeyError: if `data` has no `Shift (mm)` channel.
    """
    shift = data.get_shared("Shift (mm)")

    def gp(name):
        val = data.get_parameter(name, "000", default=float("nan"))
        return force_defaults.get(name, float("nan")) if math.isnan(val) else val

    force = fs.MagnetForce(
        shift,
        L1=gp("L1 (mm)"), L2=gp("L2 (mm)"),
        f0=gp("Fmin (pN)"), fmax=gp("Fmax (pN)"),
        alpha=gp("alpha (-)"), dt=int(gp("offset (pts)")),
    )

    for roi in data.list_rois():
        for name, default in force_defaults.items():
            if math.isnan(data.get_parameter(name, roi, default=float("nan"))):
                data.set_parameter(name, default, roi=roi)
        data.set_parameter("F (pN) calculated", 1.0, roi=roi)

    with h5py.File(data.filepath, "r+") as f:
        grp = f["Traces"] if "Traces" in f else f
        if "F (pN)" in grp:
            del grp["F (pN)"]
        grp.create_dataset("F (pN)", data=force.astype(np.float64))
