"""Tests for legacy .hdf -> new .h5 import (core/legacy_import.py).

Builds a small synthetic legacy PyTables file (mirroring oldhdf_versions/TraceIO.py's
on-disk layout, see docs/file-format.md) rather than depending on real data files
outside the test tree.
"""

import warnings

import numpy as np
import pytest
import tables
from tables import NaturalNameWarning

from trace_editor.core.legacy_import import import_legacy_hdf
from trace_editor.core.h5_data import H5Data

pytestmark = pytest.mark.io

warnings.simplefilter("ignore", NaturalNameWarning)


def write_legacy_hdf(filepath, n_beads=3, n_samples=50, with_fit=True):
    """Write a minimal legacy .hdf file with the on-disk layout TraceIO.hdf_data expects."""
    labels = [str(i) for i in range(n_beads)]
    time = np.linspace(0.0, 10.0, n_samples)
    shift = np.linspace(1.0, 9.0, n_samples).astype(np.float32)
    force = (10.0 / shift).astype(np.float32)

    with tables.open_file(filepath, mode="w") as f:
        shared_dtype = np.dtype([
            ("Time (s)", np.float64),
            ("Shift (mm)", np.float32),
            ("F (pN)", np.float32),
        ])
        shared = np.zeros(n_samples, dtype=shared_dtype)
        shared["Time (s)"] = time
        shared["Shift (mm)"] = shift
        shared["F (pN)"] = force
        f.create_table("/", " shared", obj=shared)

        bead_dtype = np.dtype([
            ("A (a.u.)", np.float64),
            ("Selection", np.float64),
            ("X (um)", np.float64),
            ("Y (um)", np.float64),
            ("Z (um)", np.float64),
        ])
        for i, label in enumerate(labels):
            bead = np.zeros(n_samples, dtype=bead_dtype)
            bead["A (a.u.)"] = 100.0 + i
            bead["Selection"] = 1.0
            bead["X (um)"] = 32.0 + i
            bead["Y (um)"] = 32.0
            bead["Z (um)"] = 10.0 + i + 0.01 * time
            f.create_table("/", label, obj=bead)

        param_fields = [(" label", "S16"), (" selected", np.float64)]
        if with_fit:
            param_fields += [
                ("Z0 (um)", np.float64), ("e:Z0 (um)", np.float64),
                ("dZ_dt (nm_s)", np.float64), ("tau (s)", np.float64),
                ("X0 (pix)", np.float64),  # no new-format equivalent -- must be dropped
            ]
        param_dtype = np.dtype(param_fields)
        rows = []
        for i, label in enumerate(labels):
            row = np.zeros(1, dtype=param_dtype)
            row[" label"] = label.encode()
            row[" selected"] = 1.0 if i != 1 else 0.0
            if with_fit:
                row["Z0 (um)"] = 10.0 + i
                row["e:Z0 (um)"] = 0.01
                row["dZ_dt (nm_s)"] = -0.5
                row["tau (s)"] = np.nan
                row["X0 (pix)"] = 640.0
            rows.append(row)
        global_row = np.zeros(1, dtype=param_dtype)
        global_row[" label"] = b"global"
        global_row[" selected"] = np.nan
        rows.append(global_row)
        params = np.concatenate(rows)
        f.create_table("/", " parameters", obj=params)

    return filepath


@pytest.fixture
def legacy_hdf(tmp_path):
    return write_legacy_hdf(tmp_path / "legacy.hdf")


class TestLegacyImport:
    def test_writes_h5_next_to_source(self, legacy_hdf):
        out = import_legacy_hdf(legacy_hdf)
        assert out.suffix == ".h5"
        assert out.exists()

    def test_refuses_to_overwrite(self, legacy_hdf):
        import_legacy_hdf(legacy_hdf)
        with pytest.raises(FileExistsError):
            import_legacy_hdf(legacy_hdf)

    def test_num_rois_matches_bead_count(self, legacy_hdf):
        out = import_legacy_hdf(legacy_hdf)
        data = H5Data(out)
        assert data.num_rois() == 3
        assert data.list_rois() == ["000", "001", "002"]

    def test_shared_channels_carried_over(self, legacy_hdf):
        out = import_legacy_hdf(legacy_hdf)
        data = H5Data(out)
        assert set(data.shared_channels()) >= {"Time (s)", "Shift (mm)", "F (pN)"}
        np.testing.assert_allclose(
            data.get_shared("Time (s)"), np.linspace(0.0, 10.0, 50)
        )

    def test_roi_channels_carried_over(self, legacy_hdf):
        out = import_legacy_hdf(legacy_hdf)
        data = H5Data(out)
        assert set(data.roi_channels()) >= {"X (um)", "Y (um)", "Z (um)", "A (a.u.)"}
        z = data.get_trace("Z (um)", "001")
        np.testing.assert_allclose(z[0], 11.0, atol=1e-3)

    def test_selection_flag_carried_over(self, legacy_hdf):
        out = import_legacy_hdf(legacy_hdf)
        data = H5Data(out)
        assert data.get_selection("000") == 1.0
        assert data.get_selection("001") == 0.0

    def test_fit_params_mapped_with_errors(self, legacy_hdf):
        out = import_legacy_hdf(legacy_hdf)
        data = H5Data(out)
        assert data.get_parameter("Z0 (um)", "000") == pytest.approx(10.0)
        assert data.get_parameter("Z0 (um) Error", "000") == pytest.approx(0.01)
        assert data.get_parameter("dZ_dt (nm_s)", "001") == pytest.approx(-0.5)

    def test_unmapped_legacy_fields_dropped(self, legacy_hdf):
        out = import_legacy_hdf(legacy_hdf)
        data = H5Data(out)
        params = data.list_parameters("000")
        assert "X0 (pix)" not in params
        assert "X0 (um)" not in params  # never set in this fixture, not carried over as NaN-only

    def test_no_fit_data_still_imports(self, tmp_path):
        legacy = write_legacy_hdf(tmp_path / "nofits.hdf", with_fit=False)
        out = import_legacy_hdf(legacy)
        data = H5Data(out)
        assert data.num_rois() == 3
        assert data.get_selection("000") == 1.0

    def test_invalid_file_raises_clean_error(self, tmp_path):
        bad = tmp_path / "bad.hdf"
        with tables.open_file(bad, mode="w") as f:
            f.create_array("/", "junk", obj=np.array([1, 2, 3]))
        with pytest.raises(ValueError):
            import_legacy_hdf(bad)

    def test_explicit_output_path(self, legacy_hdf, tmp_path):
        out_path = tmp_path / "custom_name.h5"
        out = import_legacy_hdf(legacy_hdf, h5_path=out_path)
        assert out == out_path
        assert out.exists()
