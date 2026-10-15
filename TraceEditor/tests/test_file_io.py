"""Tests for H5 file I/O operations.

Uses three data sources:
  * `temp_h5_file`         — small ad-hoc file, fast unit tests of the API.
  * `generated_workflow_h5`— realistic synthetic file (real fs models).
  * `real_h5_file`         — the committed real instrument file (skipped if absent).
"""

import pytest
import numpy as np
import h5py
from pathlib import Path
from trace_editor.core.h5_data import H5Data, SELECTED_FIELD
from trace_editor.core.file_browser import scan_folder, scan_h5_metadata, extract_tags

pytestmark = pytest.mark.io


class TestH5DataBasics:
    """Test basic H5 file reading operations."""

    def test_file_not_found(self):
        """Test that FileNotFoundError is raised for non-existent file."""
        with pytest.raises(FileNotFoundError):
            H5Data("/nonexistent/path/file.h5")

    def test_file_exists(self, temp_h5_file):
        """Test that H5Data can be instantiated with valid file."""
        h5data = H5Data(temp_h5_file)
        assert h5data.filepath == temp_h5_file

    def test_is_not_lut_file(self, temp_h5_file):
        """Test that .h5 files are not identified as LUT files."""
        h5data = H5Data(temp_h5_file)
        assert h5data.is_lut_file() is False

    def test_num_rois(self, temp_h5_file):
        """Test reading number of ROIs."""
        h5data = H5Data(temp_h5_file)
        n_rois = h5data.num_rois()
        assert n_rois == 5

    def test_list_rois(self, temp_h5_file):
        """Test listing ROI indices."""
        h5data = H5Data(temp_h5_file)
        rois = h5data.list_rois()
        assert len(rois) == 5
        assert rois == ["000", "001", "002", "003", "004"]

    def test_get_notes(self, temp_h5_file):
        """Test reading Notes attribute."""
        h5data = H5Data(temp_h5_file)
        notes = h5data.get_notes()
        assert notes == "Test experiment data"

    def test_get_date(self, temp_h5_file):
        """Test reading Date attribute."""
        h5data = H5Data(temp_h5_file)
        date = h5data.get_date()
        assert date == "2026-07-07 10:00:00"

    def test_get_axis(self, temp_h5_file):
        h5data = H5Data(temp_h5_file)
        assert h5data.get_axis() == "Shift_mm"

    def test_get_user(self, temp_h5_file):
        h5data = H5Data(temp_h5_file)
        assert h5data.get_user() == "test_user"

    def test_selected_count(self, temp_h5_file):
        h5data = H5Data(temp_h5_file)
        assert h5data.selected_count() == 5

    def test_get_lut_filename_absent_settings(self, temp_h5_file):
        # temp_h5_file has a /Settings dataset but no lut_filename key.
        h5data = H5Data(temp_h5_file)
        assert h5data.get_lut_filename() == ""

    def test_get_lut_filename_missing_settings_dataset(self, tmp_path):
        path = tmp_path / "no_settings.h5"
        with h5py.File(path, "w") as f:
            f.attrs["num_rois"] = 0
            f.create_dataset("Parameters", data=np.zeros(0, dtype=[("ROI", "i4")]))
        h5data = H5Data(path)
        assert h5data.get_lut_filename() == ""

    def test_get_lut_filename_reads_imaging_group(self, tmp_path):
        path = tmp_path / "with_lut.h5"
        with h5py.File(path, "w") as f:
            f.attrs["num_rois"] = 0
            f.create_dataset("Parameters", data=np.zeros(0, dtype=[("ROI", "i4")]))
            settings_dtype = np.dtype([("Group", "S50"), ("Key", "S50"), ("Value", "S200")])
            settings = np.array([
                (b"UI", b"dark_mode", b"True"),
                (b"Imaging", b"lut_filename", b"/data/muller/20260623/data_013.lut"),
            ], dtype=settings_dtype)
            f.create_dataset("Settings", data=settings)
        h5data = H5Data(path)
        assert h5data.get_lut_filename() == "data_013"

    def test_set_notes_overwrites_attribute(self, temp_h5_file):
        h5data = H5Data(temp_h5_file)
        h5data.set_notes("Success — updated notes")
        assert h5data.get_notes() == "Success — updated notes"

    def test_get_empty_notes(self):
        """Test reading Notes from file with empty/missing Notes."""
        import tempfile
        with tempfile.NamedTemporaryFile(suffix=".h5", delete=False) as f:
            filepath = Path(f.name)

        try:
            with h5py.File(filepath, "w") as f:
                f.attrs["num_rois"] = 1
                f.create_dataset("Parameters", data=np.array([(0, 0.0, 0.0)],
                    dtype=[("ROI", np.int32), ("X (pix)", np.float32), ("Y (pix)", np.float32)]))
                f.create_group("Traces")

            h5data = H5Data(filepath)
            notes = h5data.get_notes()
            assert notes == ""
        finally:
            filepath.unlink()


class TestChannelDiscovery:
    """Test channel discovery and listing."""

    def test_shared_channels(self, temp_h5_file):
        """Test identifying shared (non-ROI) channels."""
        h5data = H5Data(temp_h5_file)
        shared = h5data.shared_channels()
        assert "Time (s)" in shared
        assert "Shift (mm)" in shared
        assert "Intensity" in shared

    def test_roi_channels(self, temp_h5_file):
        """Test identifying per-ROI channels."""
        h5data = H5Data(temp_h5_file)
        roi_channels = h5data.roi_channels()
        assert "Z (um)" in roi_channels
        assert "X (um)" in roi_channels
        assert "Y (um)" in roi_channels
        assert "A (a.u.)" in roi_channels

    def test_all_channels(self, temp_h5_file):
        """Test listing all channels."""
        h5data = H5Data(temp_h5_file)
        all_ch = h5data.all_channels()
        shared = h5data.shared_channels()
        roi_ch = h5data.roi_channels()
        assert len(all_ch) == len(shared) + len(roi_ch)


class TestTraceDataReading:
    """Test reading trace data from file."""

    def test_get_shared_channel(self, temp_h5_file):
        """Test reading shared channel data."""
        h5data = H5Data(temp_h5_file)
        time = h5data.get_shared("Time (s)")
        assert isinstance(time, np.ndarray)
        assert len(time) == 1000
        assert time.dtype == np.float64

    def test_get_roi_trace(self, temp_h5_file):
        """Test reading per-ROI trace data."""
        h5data = H5Data(temp_h5_file)
        z = h5data.get_trace("Z (um)", "000")
        assert isinstance(z, np.ndarray)
        assert len(z) == 1000
        assert z.dtype == np.float64

    def test_get_channel_shared(self, temp_h5_file):
        """Test get_channel works for shared channels."""
        h5data = H5Data(temp_h5_file)
        shift = h5data.get_channel("Shift (mm)", "000")
        assert isinstance(shift, np.ndarray)
        assert len(shift) == 1000

    def test_get_channel_roi(self, temp_h5_file):
        """Test get_channel works for ROI channels."""
        h5data = H5Data(temp_h5_file)
        x = h5data.get_channel("X (um)", "001")
        assert isinstance(x, np.ndarray)
        assert len(x) == 1000

    def test_data_shape_consistency(self, temp_h5_file):
        """Test that all traces have consistent length."""
        h5data = H5Data(temp_h5_file)
        time = h5data.get_shared("Time (s)")
        for roi in h5data.list_rois():
            z = h5data.get_trace("Z (um)", roi)
            assert len(z) == len(time)


class TestParameterOperations:
    """Test parameter reading/writing."""

    def test_list_parameters_empty_roi(self, temp_h5_file):
        """Test listing parameters for fresh ROI (excludes instrument params)."""
        h5data = H5Data(temp_h5_file)
        params = h5data.list_parameters("000")
        # list_parameters excludes ROI, X (pix), Y (pix) but includes " selected"
        assert " selected" in params
        # Instrument fields should NOT be in list_parameters
        assert "ROI" not in params
        assert "X (pix)" not in params
        assert "Y (pix)" not in params

    def test_get_parameter_existing(self, temp_h5_file):
        """Test reading existing parameter."""
        h5data = H5Data(temp_h5_file)
        selected = h5data.get_parameter(" selected", "000")
        assert selected == 1.0

    def test_get_parameter_default(self, temp_h5_file):
        """Test reading non-existent parameter returns default."""
        h5data = H5Data(temp_h5_file)
        val = h5data.get_parameter("Z0 (um)", "000", default=99.0)
        # Should be NaN if not yet set, or default if we set it
        assert isinstance(val, (float, np.floating))

    def test_set_parameter_new_field(self, temp_h5_file):
        """Test adding new parameter to dataset."""
        h5data = H5Data(temp_h5_file)
        h5data.set_parameter("Z0 (um)", 5.5, roi="000")

        # Read it back
        val = h5data.get_parameter("Z0 (um)", "000")
        assert abs(val - 5.5) < 1e-6

    def test_set_parameter_scalar(self, temp_h5_file):
        """Test setting scalar parameter."""
        h5data = H5Data(temp_h5_file)
        h5data.set_parameter("L (bp)", 4000, roi="001")
        val = h5data.get_parameter("L (bp)", "001")
        assert val == 4000

    def test_set_selection(self, temp_h5_file):
        """Test setting ROI-level selection."""
        h5data = H5Data(temp_h5_file)
        h5data.set_selection("000", 0.0)
        selected = h5data.get_parameter(" selected", "000")
        assert selected == 0.0

    def test_set_parameters_bulk_writes_all_rois_and_fields(self, temp_h5_file):
        h5data = H5Data(temp_h5_file)
        h5data.set_parameters_bulk({
            "000": {"L (bp)": 8000.0, "R2 zWLC": 0.95},
            "001": {"L (bp)": 9000.0, "R2 zWLC": 0.90},
        })
        assert h5data.get_parameter("L (bp)", "000") == 8000.0
        assert h5data.get_parameter("R2 zWLC", "000") == 0.95
        assert h5data.get_parameter("L (bp)", "001") == 9000.0
        assert h5data.get_parameter("R2 zWLC", "001") == 0.90

    def test_set_parameters_bulk_leaves_other_rois_nan(self, temp_h5_file):
        """A field only written for some ROIs defaults to NaN for the rest,
        matching set_parameter's per-call behavior."""
        h5data = H5Data(temp_h5_file)
        h5data.set_parameters_bulk({"000": {"L (bp)": 8000.0}})
        assert np.isnan(h5data.get_parameter("L (bp)", "001", default=float("nan")))

    def test_set_parameters_bulk_matches_set_parameter(self, temp_h5_file):
        """Bulk and per-call writes produce identical stored values."""
        h5data = H5Data(temp_h5_file)
        h5data.set_parameters_bulk({
            "000": {"Z0 (um)": 12.0, "dZ_dt (nm_s)": -1.5},
            "002": {"Z0 (um)": 11.0, "dZ_dt (nm_s)": -1.0},
        })
        h5data.set_parameter("Z0 (um)", 13.0, roi="003")

        assert h5data.get_parameter("Z0 (um)", "000") == 12.0
        assert h5data.get_parameter("dZ_dt (nm_s)", "000") == -1.5
        assert h5data.get_parameter("Z0 (um)", "002") == 11.0
        assert h5data.get_parameter("Z0 (um)", "003") == 13.0

    def test_set_parameters_bulk_empty_is_noop(self, temp_h5_file):
        h5data = H5Data(temp_h5_file)
        before = h5data.list_parameters("000")
        h5data.set_parameters_bulk({})
        assert h5data.list_parameters("000") == before


class TestPointSelectionMask:
    """Per-point selection mask, stored at Traces/Selection/<roi>."""

    def test_no_mask_returns_none(self, temp_h5_file):
        """Before any point selection is made, the mask is absent."""
        h5data = H5Data(temp_h5_file)
        assert h5data.get_selection_mask("000") is None

    def test_set_then_get_roundtrips(self, temp_h5_file):
        h5data = H5Data(temp_h5_file)
        n = len(h5data.get_channel("Z (um)", "000"))
        mask = np.zeros(n, dtype=float)
        mask[10:20] = 1.0
        h5data.set_selection_mask("000", mask)

        readback = h5data.get_selection_mask("000")
        assert readback is not None
        assert np.array_equal(readback, mask)

    def test_mask_persists_across_instances(self, temp_h5_file):
        """Mirrors TestRoundTripPersistence: reopening the file sees the write."""
        n = len(H5Data(temp_h5_file).get_channel("Z (um)", "000"))
        mask = np.zeros(n, dtype=float)
        mask[:5] = 1.0
        H5Data(temp_h5_file).set_selection_mask("000", mask)

        readback = H5Data(temp_h5_file).get_selection_mask("000")
        assert np.array_equal(readback, mask)

    def test_overwrite_replaces_previous_mask(self, temp_h5_file):
        h5data = H5Data(temp_h5_file)
        n = len(h5data.get_channel("Z (um)", "000"))
        h5data.set_selection_mask("000", np.ones(n))
        h5data.set_selection_mask("000", np.zeros(n))
        assert np.array_equal(h5data.get_selection_mask("000"), np.zeros(n))

    def test_mask_is_per_roi(self, temp_h5_file):
        h5data = H5Data(temp_h5_file)
        n = len(h5data.get_channel("Z (um)", "000"))
        h5data.set_selection_mask("000", np.ones(n))
        assert h5data.get_selection_mask("001") is None

    def test_selection_channel_excluded_from_roi_channels(self, temp_h5_file):
        """Selection is a mask, not a plottable/fittable physical channel."""
        h5data = H5Data(temp_h5_file)
        n = len(h5data.get_channel("Z (um)", "000"))
        h5data.set_selection_mask("000", np.ones(n))
        assert "Selection" not in h5data.roi_channels()
        assert "Selection" not in h5data.all_channels()

    def test_stored_in_traces_group_as_float32(self, temp_h5_file):
        """Layout should match the other per-ROI channels for HDFView visibility."""
        h5data = H5Data(temp_h5_file)
        n = len(h5data.get_channel("Z (um)", "000"))
        h5data.set_selection_mask("000", np.ones(n))
        with h5py.File(temp_h5_file, "r") as f:
            ds = f["Traces"]["Selection"]["000"]
            assert ds.dtype == np.float32
            assert ds.shape == (n,)

    def test_set_multiple_parameters(self, temp_h5_file):
        """Test setting multiple parameters for same ROI."""
        h5data = H5Data(temp_h5_file)
        h5data.set_parameter("Z0 (um)", 10.0, roi="002")
        h5data.set_parameter("dZ_dt (nm_s)", 5.0, roi="002")
        h5data.set_parameter("tau (s)", 20.0, roi="002")

        z0 = h5data.get_parameter("Z0 (um)", "002")
        dz = h5data.get_parameter("dZ_dt (nm_s)", "002")
        tau = h5data.get_parameter("tau (s)", "002")

        assert abs(z0 - 10.0) < 1e-6
        assert abs(dz - 5.0) < 1e-6
        assert abs(tau - 20.0) < 1e-6


class TestDataIntegrity:
    """Test data integrity and edge cases."""

    def test_roi_index_padding(self, temp_h5_file):
        """Test that ROI indices are properly zero-padded."""
        h5data = H5Data(temp_h5_file)
        rois = h5data.list_rois()
        for roi in rois:
            assert len(roi) == 3
            assert roi.isdigit()

    def test_float_dtype_conversion(self, temp_h5_file):
        """Test that data is properly converted to float64."""
        h5data = H5Data(temp_h5_file)
        z = h5data.get_trace("Z (um)", "000")
        # Original data is float32, should be converted to float64
        assert z.dtype == np.float64

    def test_time_data_increasing(self, temp_h5_file):
        """Test that time data is monotonically increasing."""
        h5data = H5Data(temp_h5_file)
        time = h5data.get_shared("Time (s)")
        assert np.all(np.diff(time) > 0)

    def test_read_same_roi_twice(self, temp_h5_file):
        """Test that reading same ROI twice gives identical results."""
        h5data = H5Data(temp_h5_file)
        z1 = h5data.get_trace("Z (um)", "000")
        z2 = h5data.get_trace("Z (um)", "000")
        assert np.allclose(z1, z2)


# ---------------------------------------------------------------------------
# Realistic synthetic workflow file (real fs models, dwell-based Shift sweep)
# ---------------------------------------------------------------------------

class TestWorkflowFileReading:
    """Read the realistic synthetic file and check its channels/shape/values."""

    def test_channels_match_new_format(self, generated_workflow_h5):
        h5 = H5Data(generated_workflow_h5)
        assert set(h5.roi_channels()) == {"Z (um)", "X (um)", "Y (um)", "A (a.u.)"}
        assert "Time (s)" in h5.shared_channels()
        assert "Shift (mm)" in h5.shared_channels()

    def test_all_roi_traces_same_length_as_time(self, generated_workflow_h5):
        h5 = H5Data(generated_workflow_h5)
        n = len(h5.get_shared("Time (s)"))
        for roi in h5.list_rois():
            for ch in h5.roi_channels():
                assert len(h5.get_trace(ch, roi)) == n

    def test_shift_sweep_dwells_top_and_bottom(self, generated_workflow_h5):
        """The Shift trace starts and ends at the top and reaches the bottom."""
        h5 = H5Data(generated_workflow_h5)
        shift = h5.get_shared("Shift (mm)")
        assert shift[0] == pytest.approx(10.8, abs=0.2)   # start at top
        assert shift[-1] == pytest.approx(10.8, abs=0.2)  # end at top
        assert shift.min() == pytest.approx(1.0, abs=0.2)  # reaches bottom
        # Dwell at the top: several samples pinned near the maximum.
        assert np.sum(shift > 10.5) > 20

    def test_z_trace_dominated_by_drift(self, generated_workflow_h5):
        """Observed Z sits near the drift baseline (~11-12 um), not the WLC height."""
        h5 = H5Data(generated_workflow_h5)
        z = h5.get_trace("Z (um)", "000")
        assert z.mean() > 8.0

    def test_amplitude_channel_present_and_high(self, generated_workflow_h5):
        """A (a.u.) is high/flat so the drift-fit amplitude filter is a no-op."""
        h5 = H5Data(generated_workflow_h5)
        amp = h5.get_trace("A (a.u.)", "000")
        assert amp.mean() > 50.0


# ---------------------------------------------------------------------------
# Round-trip persistence: writes must survive reopening from disk
# ---------------------------------------------------------------------------

class TestRoundTripPersistence:
    """Parameters/selection written via H5Data persist on disk across reopen."""

    def test_parameter_survives_reopen(self, generated_workflow_h5):
        H5Data(generated_workflow_h5).set_parameter("L (bp)", 9999.0, roi="001")
        # Fresh instance -> reads from disk, not an in-memory cache.
        assert H5Data(generated_workflow_h5).get_parameter("L (bp)", "001") == 9999.0

    def test_selection_survives_reopen(self, generated_workflow_h5):
        H5Data(generated_workflow_h5).set_selection("001", 0.0)
        assert H5Data(generated_workflow_h5).get_selection("001") == 0.0

    def test_new_field_defaults_other_rois_to_nan(self, generated_workflow_h5):
        """Adding a field for one ROI leaves the others NaN, not zero."""
        h5 = H5Data(generated_workflow_h5)
        h5.set_parameter("A (nm)", 45.0, roi="000")
        reopened = H5Data(generated_workflow_h5)
        assert reopened.get_parameter("A (nm)", "000") == 45.0
        # ROI 001 was never written -> get_parameter returns the default (NaN sentinel).
        assert np.isnan(reopened.get_parameter("A (nm)", "001", default=np.nan))

    def test_overwrite_parameter(self, generated_workflow_h5):
        h5 = H5Data(generated_workflow_h5)
        h5.set_parameter("Z0 (um)", 1.0, roi="000")
        h5.set_parameter("Z0 (um)", 2.0, roi="000")
        assert H5Data(generated_workflow_h5).get_parameter("Z0 (um)", "000") == 2.0

    def test_traces_untouched_by_param_write(self, generated_workflow_h5):
        """Writing parameters must not alter the trace data (set_parameter
        rewrites only /Parameters)."""
        before = H5Data(generated_workflow_h5).get_trace("Z (um)", "000").copy()
        H5Data(generated_workflow_h5).set_parameter("tau (s)", 5.0, roi="000")
        after = H5Data(generated_workflow_h5).get_trace("Z (um)", "000")
        assert np.array_equal(before, after)


# ---------------------------------------------------------------------------
# Real instrument file (committed) — format parity and readability
# ---------------------------------------------------------------------------

class TestRealDataReading:
    """Smoke-read the real file. Skipped if real_data.h5 is not present."""

    @pytest.fixture(autouse=True)
    def _skip_if_absent(self, real_h5_file):
        if real_h5_file is None:
            pytest.skip("real_data.h5 not available")
        self.h5 = H5Data(real_h5_file)

    def test_metadata_readable(self):
        assert self.h5.num_rois() == 51
        assert self.h5.get_date() == "2026-06-23 16:31:26"
        assert "construct" in self.h5.get_notes()

    def test_channels_match_synthetic_format(self):
        """Real and synthetic files expose the same per-ROI channel set."""
        assert set(self.h5.roi_channels()) == {"Z (um)", "X (um)", "Y (um)", "A (a.u.)"}
        assert "Shift (mm)" in self.h5.shared_channels()
        assert "F (pN)" in self.h5.shared_channels()  # force already calculated

    def test_traces_readable_and_consistent(self):
        n = len(self.h5.get_shared("Time (s)"))
        for roi in ("000", "008", "050"):
            z = self.h5.get_trace("Z (um)", roi)
            assert len(z) == n
            assert z.dtype == np.float64

    def test_existing_fit_params_readable(self):
        """Bead 008 has stored ExponentialDrift + zWLC fits."""
        assert self.h5.get_parameter("L (bp)", "008") == pytest.approx(6000.0)
        r2 = self.h5.get_parameter("R2 zWLC", "008")
        assert 0.0 < r2 <= 1.0

    def test_roi_indices_zero_padded(self):
        rois = self.h5.list_rois()
        assert rois[0] == "000"
        assert rois[-1] == "050"
        assert all(len(r) == 3 for r in rois)


class TestExtractTags:
    """Free-form #hashtag tagging for the Open File browser."""

    @pytest.mark.parametrize("notes,expected", [
        ("Fail — bead drift too high #fail", ["fail"]),
        ("#Success great pulling curve", ["success"]),
        ("#lowdensity #retry same tag #retry again", ["lowdensity", "retry"]),
        ("Construct A, nothing notable", []),
        ("", []),
        ("no hashtag here, just a # alone", []),
    ])
    def test_tags(self, notes, expected):
        assert extract_tags(notes) == expected


class TestFileBrowserScan:
    """Metadata scanning for the Open File browser dialog."""

    def test_scan_single_file(self, temp_h5_file):
        row = scan_h5_metadata(temp_h5_file)
        assert row["axis"] == "Shift_mm"
        assert row["num_rois"] == 5
        assert row["selected"] == 5
        assert row["error"] == ""
        assert row["tags"] == []

    def test_scan_missing_attrs_reports_error_not_raise(self, tmp_path):
        bad_file = tmp_path / "broken.h5"
        with h5py.File(bad_file, "w") as f:
            f.attrs["Notes"] = "no num_rois here"
        row = scan_h5_metadata(bad_file)
        assert row["error"] != ""
        assert row["num_rois"] == 0

    def test_scan_folder_finds_h5_files(self, temp_h5_file, tmp_path):
        import shutil
        shutil.copy(temp_h5_file, tmp_path / "a.h5")
        shutil.copy(temp_h5_file, tmp_path / "b.h5")
        (tmp_path / "not_h5.txt").write_text("ignore me")
        rows = scan_folder(tmp_path)
        assert {r["filename"] for r in rows} == {"a.h5", "b.h5"}
