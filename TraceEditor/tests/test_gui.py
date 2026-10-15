"""GUI tests (pytest-qt), driving real TraceEditor behavior.

Qt runs offscreen (set in conftest via QT_QPA_PLATFORM=offscreen) so no window
appears and these are CI-safe. Files are loaded headlessly with
`conftest.load_file_headless`, bypassing the native file dialog, so we can drive
the real open -> select ROI -> plot -> fit flow and assert on outcomes (plotted
lines, populated combos, parameters written to the h5).
"""

from unittest.mock import patch

import numpy as np
import pytest
from PySide6.QtWidgets import QMessageBox

from trace_editor.core import force_spectroscopy as fs
from trace_editor.core.h5_data import prepare_fit
from trace_editor.gui.main_window import (
    _apply_point_selection, _get_fit_selection, FileBrowserDialog,
)
from conftest import load_file_headless, TRUE_L_BP

pytestmark = pytest.mark.gui


class TestConstruction:
    """The window builds with its real widgets and plot canvas wired up."""

    def test_title(self, editor):
        assert editor.windowTitle() == "Trace Editor"

    def test_plot_canvas_and_axes_exist(self, editor):
        # _setup_plot builds a real matplotlib Figure/axes/canvas.
        assert editor.canvas is not None
        assert editor.ax is not None
        assert editor.fig is not None

    def test_key_widgets_present(self, editor):
        # These are referenced throughout the app; their absence is a real bug.
        for name in ("fileComboBox", "roiComboBox", "xChannelComboBox",
                     "yChannelComboBox", "modelComboBox", "paramsListWidget"):
            assert hasattr(editor, name), f"missing widget {name}"

    def test_starts_empty(self, editor):
        assert editor._current_file is None
        assert editor._current_roi is None
        assert editor.fileComboBox.count() == 0


class TestFileLoading:
    """Loading a file populates the file/ROI/channel combos and a current ROI."""

    def test_file_registered_and_active(self, editor_with_file):
        win, _ = editor_with_file
        assert win._current_file is not None
        assert win.fileComboBox.count() == 1

    def test_roi_combo_populated(self, editor_with_file):
        win, data = editor_with_file
        assert win.roiComboBox.count() == data.num_rois()
        assert win._current_roi is not None

    def test_channel_combos_populated(self, editor_with_file):
        win, _ = editor_with_file
        x_items = [win.xChannelComboBox.itemText(i)
                   for i in range(win.xChannelComboBox.count())]
        y_items = [win.yChannelComboBox.itemText(i)
                   for i in range(win.yChannelComboBox.count())]
        assert "Time (s)" in x_items
        assert "Z (um)" in y_items
        # Sensible defaults chosen by _update_channels.
        assert win.xChannelComboBox.currentText() == "Time (s)"
        assert win.yChannelComboBox.currentText() == "Z (um)"

    def test_notes_panel_filled(self, editor_with_file):
        win, _ = editor_with_file
        # notesText is a QLabel.
        assert win.notesText.text().strip() != ""

    def test_file_meta_panel_shows_axis_lut_selected_tags(self, editor_with_file):
        win, data = editor_with_file
        meta = win.fileMetaText.text()
        assert "Date:" in meta
        assert "Axis:" in meta
        assert "LUT:" in meta
        assert "Selected:" in meta
        assert "Tags:" in meta

    def test_notes_panel_no_longer_duplicates_date(self, editor_with_file):
        win, data = editor_with_file
        assert data.get_date().strip() not in win.notesText.text()

    def test_edit_notes_updates_panel_and_file(self, editor_with_file, monkeypatch):
        win, data = editor_with_file
        from PySide6.QtWidgets import QInputDialog
        monkeypatch.setattr(QInputDialog, "getMultiLineText",
                             staticmethod(lambda *a, **k: ("#fail new notes from edit", True)))
        win._edit_notes()
        assert data.get_notes() == "#fail new notes from edit"
        assert "new notes from edit" in win.notesText.text()
        assert "#fail" in win.fileMetaText.text()

    def test_edit_notes_cancelled_leaves_file_unchanged(self, editor_with_file, monkeypatch):
        win, data = editor_with_file
        original = data.get_notes()
        from PySide6.QtWidgets import QInputDialog
        monkeypatch.setattr(QInputDialog, "getMultiLineText",
                             staticmethod(lambda *a, **k: ("ignored text", False)))
        win._edit_notes()
        assert data.get_notes() == original

    def test_edit_notes_without_open_file_is_safe(self, editor):
        editor._edit_notes()  # nothing loaded — must not raise


def _scatter_offsets(ax):
    """Y-values of the first scatter collection drawn by redraw()."""
    assert ax.collections, "redraw drew no scatter points"
    return ax.collections[0].get_offsets().data[:, 1]


class TestPlotting:
    """redraw() actually draws the selected trace onto the axes.

    Traces are rendered with ax.scatter, so points live in ax.collections
    (ax.lines only gets a model curve when a fit is displayed).
    """

    def test_redraw_scatters_points(self, editor_with_file):
        win, _ = editor_with_file
        win.redraw()
        assert len(win.ax.collections) >= 1

    def test_plotted_data_matches_trace(self, editor_with_file):
        win, data = editor_with_file
        win.xChannelComboBox.setCurrentText("Time (s)")
        win.yChannelComboBox.setCurrentText("A (a.u.)")  # not z-corrected
        win.redraw()
        expected = data.get_channel("A (a.u.)", win._current_roi)
        assert np.allclose(_scatter_offsets(win.ax), expected, rtol=1e-4)

    def test_changing_roi_updates_plot(self, editor_with_file):
        win, _ = editor_with_file
        win.yChannelComboBox.setCurrentText("A (a.u.)")
        win.roiComboBox.setCurrentIndex(0)
        win.redraw()
        y0 = _scatter_offsets(win.ax).copy()
        win.roiComboBox.setCurrentIndex(1)  # triggers _on_roi_changed -> redraw
        y1 = _scatter_offsets(win.ax)
        assert not np.array_equal(y0, y1)

    def test_next_roi_shortcut_advances(self, editor_with_file):
        win, _ = editor_with_file
        win.roiComboBox.setCurrentIndex(0)
        win._next_roi()
        assert win.roiComboBox.currentIndex() == 1

    def test_next_roi_stops_at_last_roi_does_not_change_file(self, editor_with_file):
        # Deliberately deviates from the legacy app: ROI nav must never
        # silently switch which file is open — F9/F8 do that explicitly.
        win, _ = editor_with_file
        last = win.roiComboBox.count() - 1
        win.roiComboBox.setCurrentIndex(last)
        file_before = win.fileComboBox.currentIndex()
        win._next_roi()
        assert win.roiComboBox.currentIndex() == last
        assert win.fileComboBox.currentIndex() == file_before

    def test_prev_roi_stops_at_first_roi_does_not_change_file(self, editor_with_file):
        win, _ = editor_with_file
        win.roiComboBox.setCurrentIndex(0)
        file_before = win.fileComboBox.currentIndex()
        win._prev_roi()
        assert win.roiComboBox.currentIndex() == 0
        assert win.fileComboBox.currentIndex() == file_before


class TestSelectionToggle:
    """Toggling ROI selection writes through to the h5 file."""

    def test_toggle_flips_stored_selection(self, editor_with_file):
        win, data = editor_with_file
        roi = win._current_roi
        before = data.get_selection(roi)
        win._toggle_current_roi_selected()
        after = data.get_selection(roi)
        assert after != before
        assert {before, after} == {0.0, 1.0}


class TestPointSelection:
    """Point-level criteria (Select Points > By Criteria) persist and drive both
    the live plot coloring and the fit point mask when Manual selection mode
    is active (see _get_fit_selection: "default" mode uses each model's
    configured criteria instead, and ignores the persisted mask).
    """

    def test_apply_point_selection_persists_mask(self, editor_with_file):
        win, data = editor_with_file
        roi = win._current_roi
        z = data.get_channel("Z (um)", roi)
        lo, hi = float(np.percentile(z, 30)), float(np.percentile(z, 70))

        criteria = [{"channel": "Z (um)", "min": str(lo), "max": str(hi),
                     "derivative": False, "invert": False}]
        sel = _apply_point_selection(criteria, data, roi)

        mask = data.get_selection_mask(roi)
        assert mask is not None
        assert np.array_equal(mask, sel)
        assert 0 < mask.sum() < len(mask)

    def test_fit_selection_uses_persisted_mask_for_zwlc_in_manual_mode(self, editor_with_file):
        win, data = editor_with_file
        roi = win._current_roi
        z = data.get_channel("Z (um)", roi)
        lo, hi = float(np.percentile(z, 30)), float(np.percentile(z, 70))

        criteria = [{"channel": "Z (um)", "min": str(lo), "max": str(hi),
                     "derivative": False, "invert": False}]
        sel = _apply_point_selection(criteria, data, roi)

        fit_sel = _get_fit_selection(data, roi, mode="manual")
        assert fit_sel is not None
        assert np.array_equal(fit_sel, sel == 1.0)

    def test_default_mode_ignores_persisted_mask_for_zwlc(self, editor_with_file):
        """zWLC has no [model_selection] config entry, so Default mode uses
        every point even when a manual mask happens to be persisted."""
        win, data = editor_with_file
        roi = win._current_roi
        z = data.get_channel("Z (um)", roi)
        lo, hi = float(np.percentile(z, 30)), float(np.percentile(z, 70))
        criteria = [{"channel": "Z (um)", "min": str(lo), "max": str(hi),
                     "derivative": False, "invert": False}]
        _apply_point_selection(criteria, data, roi)

        assert _get_fit_selection(data, roi, mode="default") is None

    def test_no_selection_means_fit_uses_all_points(self, editor_with_file):
        win, data = editor_with_file
        roi = win._current_roi
        assert data.get_selection_mask(roi) is None
        assert _get_fit_selection(data, roi, mode="manual") is None

    def test_fit_through_gui_respects_point_selection(self, editor_with_file, monkeypatch):
        """End-to-end: _run_fit passes exactly the persisted point mask through
        to apply_fit, restricting the point count actually used for the fit.

        Uses the noisy ROI so a restricted-range fit isn't numerically
        identical to the full-range fit (noiseless data fits exactly either
        way, which would make a result-based comparison meaningless).
        """
        import trace_editor.gui.main_window as mw

        win, data = editor_with_file
        win.calculate_force()
        win._selection_mode = "manual"
        roi = "001"  # small measurement noise

        calls = []
        real_apply_fit = mw.apply_fit

        def spying_apply_fit(func_name, x, y, selection=None, initial_parameters=None):
            calls.append(None if selection is None else int(np.sum(selection == 1.0)))
            return real_apply_fit(func_name, x, y, selection=selection,
                                   initial_parameters=initial_parameters)

        monkeypatch.setattr(mw, "apply_fit", spying_apply_fit)

        n_total = len(data.get_channel("Z (um)", roi))

        _, params, names = prepare_fit("zWLC")
        params["L"].set(value=5000.0, vary=True)
        params["A"].set(value=50.0, vary=False)
        params["S"].set(value=1000.0, vary=False)
        params["z_bead"].set(value=0.0, vary=True)
        win._run_fit(data, roi, "zWLC", params, names)
        assert calls[-1] is None  # no selection yet -> full trace

        # _run_fit stored a fitted Zbead, which _apply_point_selection applies
        # before thresholding (matching TraceIO), so build criteria against
        # that corrected Z, not the raw channel.
        z_bead = data.get_parameter("Zbead (um)", roi, default=0.0)
        z_corrected = data.get_channel("Z (um)", roi) - z_bead
        lo, hi = float(np.percentile(z_corrected, 0)), float(np.percentile(z_corrected, 50))
        criteria = [{"channel": "Z (um)", "min": str(lo), "max": str(hi),
                     "derivative": False, "invert": False}]
        sel = _apply_point_selection(criteria, data, roi)

        win._run_fit(data, roi, "zWLC", params, names)
        assert calls[-1] == int(sel.sum())
        assert 0 < calls[-1] < n_total

    def test_redraw_colors_by_persisted_point_mask(self, editor_with_file):
        """redraw() should color exactly the selected points, not the
        whole-trace ' selected' flag (which stays 1.0 here)."""
        win, data = editor_with_file
        roi = win._current_roi
        win._selection_mode = "manual"
        win.yChannelComboBox.setCurrentText("A (a.u.)")

        z = data.get_channel("Z (um)", roi)
        lo, hi = float(np.percentile(z, 30)), float(np.percentile(z, 70))
        criteria = [{"channel": "Z (um)", "min": str(lo), "max": str(hi),
                     "derivative": False, "invert": False}]
        sel = _apply_point_selection(criteria, data, roi)

        assert data.get_selection(roi) == 1.0  # whole-trace flag untouched
        win.redraw()

        # First collection = background (all points), second = highlighted selection.
        assert len(win.ax.collections) >= 2
        highlighted = win.ax.collections[1].get_offsets().data
        assert len(highlighted) == int(sel.sum())

    def test_redraw_highlights_all_points_when_no_default_criteria_configured(self, editor_with_file):
        """A model with no [model_selection.<model>] entry (e.g. zWLC) means
        the fit uses every point (_get_fit_selection returns None) -- the
        highlight must show every point too, regardless of the unrelated
        whole-trace ' selected' flag (regression: this used to fall back to
        that flag, so a deselected ROI showed nothing highlighted even
        though every point is what the fit actually uses)."""
        win, data = editor_with_file
        roi = win._current_roi
        win._selection_mode = "default"
        win.modelComboBox.setCurrentText("zWLC")
        data.set_selection(roi, 0.0)  # whole-trace flag off -- must not affect the highlight

        win.redraw()

        assert len(win.ax.collections) >= 2
        n_points = len(data.get_channel(win.xChannelComboBox.currentText(), roi))
        highlighted = win.ax.collections[1].get_offsets().data
        assert len(highlighted) == n_points


class TestCalculateForce:
    """The Calculate Force action writes F (pN) and marks ROIs."""

    def test_calculate_force_writes_channel(self, editor_with_file):
        win, data = editor_with_file
        win.calculate_force()
        assert "F (pN)" in data.shared_channels()
        # MagnetForce is positive and finite everywhere.
        force = data.get_shared("F (pN)")
        assert np.all(np.isfinite(force)) and np.all(force > 0)
        # Every ROI flagged as force-calculated.
        for roi in data.list_rois():
            assert data.get_parameter("F (pN) calculated", roi) == 1.0


class TestFitThroughGui:
    """Drive the real fit path (_run_fit) headlessly and check results land."""

    def _fit(self, win, data, roi, model_name, param_overrides):
        """Build params via prepare_fit, apply overrides, run the GUI fit."""
        _, params, names = prepare_fit(model_name)
        for key, (value, vary) in param_overrides.items():
            params[key].set(value=value, vary=vary)
        win._run_fit(data, roi, model_name, params, names)

    def test_drift_then_wlc_writes_params(self, editor_with_file):
        win, data = editor_with_file
        # Force must exist before a WLC (F vs Z) fit.
        win.calculate_force()
        roi = "000"  # noiseless bead

        # Stage 1: ExponentialDrift (free Z0, dZ_dt; tau fixed).
        self._fit(win, data, roi, "ExponentialDrift", {
            "Z0": (12.0, True), "dZ_dt": (0.0, True), "tau": (1e9, False),
        })
        assert np.isfinite(data.get_parameter("Z0 (um)", roi))
        assert np.isfinite(data.get_parameter("dZ_dt (nm_s)", roi))

        # Stage 2: zWLC (free L, z_bead; A, S fixed). _run_fit applies the
        # drift + Zbead correction to Z internally via _correct_z_for_fit.
        self._fit(win, data, roi, "zWLC", {
            "L": (5000.0, True), "A": (50.0, False),
            "S": (1000.0, False), "z_bead": (0.0, True),
        })
        L = data.get_parameter("L (bp)", roi)
        r2 = data.get_parameter("R2 zWLC", roi)
        assert abs(L - TRUE_L_BP) < 500.0     # recovered ~10000 bp
        assert r2 > 0.9

    def test_fit_populates_param_list_widget(self, editor_with_file):
        win, data = editor_with_file
        win.calculate_force()
        self._fit(win, data, "000", "ExponentialDrift", {
            "Z0": (12.0, True), "dZ_dt": (0.0, True), "tau": (1e9, False),
        })
        win.redraw()  # refreshes the params list for the current ROI
        labels = [win.paramsListWidget.item(i).text()
                  for i in range(win.paramsListWidget.count())]
        assert any("Z0" in t for t in labels)


class TestBatchFit:
    """Batch fit runs on a real QThread, which delegates to
    core.batch_fit.run_batch (a ProcessPoolExecutor) for the actual fitting.
    These guard against the thread (and the process pool it owns for the
    duration of the run) being torn down while still running -- crashes, not
    just warnings -- and are unaffected by the underlying worker being a
    process pool: monkeypatching functions in main_window has no effect
    inside worker *processes* (they don't share the patched module state), so
    these tests slow things down with enough ROIs/real fit work instead of a
    patched sleep.
    """

    def _params(self, model_name, overrides):
        _, params, names = prepare_fit(model_name)
        for key, (value, vary) in overrides.items():
            params[key].set(value=value, vary=vary)
        return params, names

    def test_batch_fit_runs_and_writes_params(self, editor_with_file, qtbot):
        win, data = editor_with_file
        win.calculate_force()
        for roi in data.list_rois():
            data.set_selection(roi, 1.0)
        tasks = win._batch_roi_tasks()
        assert tasks

        params, names = self._params("ExponentialDrift", {
            "Z0": (12.0, True), "dZ_dt": (0.0, True), "tau": (1e9, False),
        })
        win._start_batch(tasks, "ExponentialDrift", params, names)
        assert win._batch_thread is not None
        qtbot.waitUntil(lambda: win._batch_thread is None, timeout=10000)

        for _, roi in tasks:
            assert np.isfinite(data.get_parameter("Z0 (um)", roi))

    def test_start_batch_rejects_reentry_while_running(self, editor_with_file, qtbot):
        """A second _start_batch call while one is in flight must not replace
        the running thread/worker (regression: replacing them let the old
        QThread be garbage-collected mid-run -> hard crash on close/GC).

        Uses max_workers=1 (via a stop_cb-free single-task-at-a-time pool is
        not directly controllable from here, so instead we rely on process
        pool startup overhead -- spawning workers takes tens of ms -- to keep
        the thread "running" long enough to attempt the reentrant call).
        """
        win, data = editor_with_file
        win.calculate_force()
        for roi in data.list_rois():
            data.set_selection(roi, 1.0)
        tasks = win._batch_roi_tasks()

        params, names = self._params("ExponentialDrift", {
            "Z0": (12.0, True), "dZ_dt": (0.0, True), "tau": (1e9, False),
        })
        win._start_batch(tasks, "ExponentialDrift", params, names)
        first_thread = win._batch_thread
        assert first_thread is not None and first_thread.isRunning()

        win._start_batch(tasks, "ExponentialDrift", params, names)
        assert win._batch_thread is first_thread  # unchanged: re-entry rejected

        qtbot.waitUntil(lambda: win._batch_thread is None, timeout=10000)

    def test_close_while_batch_running_does_not_crash(self, qapp, generated_workflow_h5, qtbot):
        """Closing the window mid-batch must stop and join the worker thread
        (and the process pool it owns for the run's duration) before the
        window is torn down, instead of letting Qt destroy a still-running
        QThread or leaking worker processes.
        """
        import multiprocessing
        from trace_editor.gui.main_window import TraceEditor

        win = TraceEditor()
        data = load_file_headless(win, generated_workflow_h5)
        win.calculate_force()
        for roi in data.list_rois():
            data.set_selection(roi, 1.0)
        tasks = win._batch_roi_tasks()

        params, names = self._params("ExponentialDrift", {
            "Z0": (12.0, True), "dZ_dt": (0.0, True), "tau": (1e9, False),
        })
        win._start_batch(tasks, "ExponentialDrift", params, names)
        assert win._batch_thread is not None and win._batch_thread.isRunning()

        win.close()  # must not raise / must not leave a running QThread behind
        assert win._batch_thread is None or not win._batch_thread.isRunning()
        assert multiprocessing.active_children() == []


class TestResetAndClose:
    """Reset and close behave without leaving stale state."""

    def test_close_all_clears_state(self, editor_with_file):
        win, _ = editor_with_file
        win._close_all_files()
        assert win._current_file is None
        assert win.fileComboBox.count() == 0
        assert win.roiComboBox.count() == 0

    def test_close_without_file_is_safe(self, editor):
        editor._close_all_files()   # nothing loaded — must not raise
        assert editor._current_file is None

    def test_remove_all_processing_clears_selection_mask(self, editor_with_file):
        """Data reset wipes everything derived from raw instrument data --
        fit params, F (pN), and the manual point-selection mask alike."""
        win, data = editor_with_file
        roi = win._current_roi
        data.set_selection_mask(roi, np.array([1.0, 0.0, 1.0]))
        assert data.get_selection_mask(roi) is not None

        win.remove_all_processing()

        assert data.get_selection_mask(roi) is None


class TestInspectH5:
    """File > Inspect H5… dispatches to the H5Web VS Code extension."""

    def test_no_file_open_is_safe(self, editor):
        editor._inspect_h5()   # nothing loaded — must not raise
        assert "Open a file" in editor.statusbar.currentMessage()

    def test_success_dispatches_to_current_file(self, editor_with_file):
        win, data = editor_with_file
        with patch("trace_editor.gui.main_window.open_in_vscode_h5web") as mock_open:
            win._inspect_h5()
        mock_open.assert_called_once_with(win._current_file)
        assert "H5Web" in win.statusbar.currentMessage()

    def test_failure_shows_dialog_with_fallback_options(self, editor_with_file):
        win, data = editor_with_file
        from trace_editor.core.inspect_h5 import InspectH5Error

        with patch(
            "trace_editor.gui.main_window.open_in_vscode_h5web",
            side_effect=InspectH5Error("no code", "no_code_cli"),
        ), patch.object(QMessageBox, "exec", return_value=None) as mock_exec, \
           patch.object(QMessageBox, "clickedButton", return_value=None):
            win._inspect_h5()   # Cancel (clickedButton=None) — must not raise or hang
        assert mock_exec.called

    def test_try_again_retries_then_succeeds(self, editor_with_file):
        win, data = editor_with_file
        from trace_editor.core.inspect_h5 import InspectH5Error

        call_count = {"n": 0}

        def flaky_open(path):
            call_count["n"] += 1
            if call_count["n"] == 1:
                raise InspectH5Error("no code", "no_code_cli")

        with patch("trace_editor.gui.main_window.open_in_vscode_h5web", side_effect=flaky_open), \
             patch.object(QMessageBox, "exec", return_value=None), \
             patch.object(QMessageBox, "clickedButton") as mock_clicked:
            # First exec() call returns the "Try Again" button; loop re-invokes and succeeds.
            box_holder = {}
            orig_add = QMessageBox.addButton

            def capturing_add_button(self, *args, **kwargs):
                btn = orig_add(self, *args, **kwargs)
                if args and args[0] == "Try Again":
                    box_holder["retry"] = btn
                return btn

            with patch.object(QMessageBox, "addButton", capturing_add_button):
                mock_clicked.side_effect = lambda: box_holder.get("retry")
                win._inspect_h5()

        assert call_count["n"] == 2


class TestFileBrowserDialog:
    """Open File dialog: scans a folder into a filterable/sortable table."""

    def test_scans_folder_on_open(self, qapp, tmp_path):
        import h5py
        with h5py.File(tmp_path / "a.h5", "w") as f:
            f.attrs["axis"] = "Shift_mm"
            f.attrs["Notes"] = "#success good run"
            f.attrs["num_rois"] = 1
            f.create_dataset("Parameters", data=np.zeros(
                1, dtype=[("ROI", "i4"), (" selected", "f8")]))

        dlg = FileBrowserDialog(str(tmp_path))
        assert dlg._model.rowCount() == 1
        assert dlg._model.item(0, 0).text() == "a.h5"
        assert dlg._model.item(0, 2).text() == "#success"

    def test_focus_mm_rows_are_shown_not_hidden(self, qapp, tmp_path):
        # Focus_mm is a focus-change experiment, not a LUT file (LUT files
        # use a separate .LUT format and never appear in this browser) — so
        # both rows should always be visible, with no hide/show toggle.
        import h5py
        with h5py.File(tmp_path / "shift.h5", "w") as f:
            f.attrs["axis"] = "Shift_mm"
            f.attrs["num_rois"] = 0
            f.create_dataset("Parameters", data=np.zeros(0, dtype=[("ROI", "i4")]))
        with h5py.File(tmp_path / "focus.h5", "w") as f:
            f.attrs["axis"] = "Focus_mm"
            f.attrs["num_rois"] = 0
            f.create_dataset("Parameters", data=np.zeros(0, dtype=[("ROI", "i4")]))

        dlg = FileBrowserDialog(str(tmp_path))
        assert dlg._proxy.rowCount() == 2

    def test_lut_filename_column_populated(self, qapp, tmp_path):
        import h5py
        with h5py.File(tmp_path / "a.h5", "w") as f:
            f.attrs["axis"] = "Shift_mm"
            f.attrs["num_rois"] = 0
            f.create_dataset("Parameters", data=np.zeros(0, dtype=[("ROI", "i4")]))
            settings_dtype = np.dtype([("Group", "S50"), ("Key", "S50"), ("Value", "S200")])
            settings = np.array(
                [(b"Imaging", b"lut_filename", b"/data/muller/20260623/data_013.lut")],
                dtype=settings_dtype)
            f.create_dataset("Settings", data=settings)

        dlg = FileBrowserDialog(str(tmp_path))
        lut_col = dlg._model.horizontalHeaderItem(4)
        assert lut_col.text() == "LUT"
        assert dlg._model.item(0, 4).text() == "data_013"

    def test_notes_tooltip_shows_full_text(self, qapp, tmp_path):
        import h5py
        from trace_editor.gui.main_window import _NOTES_COL
        long_notes = "New Beads at lower density stuck to bottom of the flow cell, will retry"
        with h5py.File(tmp_path / "a.h5", "w") as f:
            f.attrs["axis"] = "Shift_mm"
            f.attrs["Notes"] = long_notes
            f.attrs["num_rois"] = 0
            f.create_dataset("Parameters", data=np.zeros(0, dtype=[("ROI", "i4")]))

        dlg = FileBrowserDialog(str(tmp_path))
        assert dlg._model.item(0, _NOTES_COL).toolTip() == long_notes

    def test_double_click_selects_and_accepts(self, qapp, tmp_path):
        import h5py
        with h5py.File(tmp_path / "a.h5", "w") as f:
            f.attrs["axis"] = "Shift_mm"
            f.attrs["num_rois"] = 0
            f.create_dataset("Parameters", data=np.zeros(0, dtype=[("ROI", "i4")]))

        dlg = FileBrowserDialog(str(tmp_path))
        dlg._table.selectRow(0)
        dlg._on_ok()
        assert dlg.selected_paths == [str(tmp_path / "a.h5")]

    def test_query_filters_by_axis_excluding_focus(self, qapp, tmp_path):
        # Regression: a plain substring search for "shift" used to also match
        # Focus_mm rows whose Notes happened to contain the word "shift".
        import h5py
        with h5py.File(tmp_path / "shift.h5", "w") as f:
            f.attrs["axis"] = "Shift_mm"
            f.attrs["Notes"] = "moved the stage"
            f.attrs["num_rois"] = 0
            f.create_dataset("Parameters", data=np.zeros(0, dtype=[("ROI", "i4")]))
        with h5py.File(tmp_path / "focus.h5", "w") as f:
            f.attrs["axis"] = "Focus_mm"
            f.attrs["Notes"] = "shift = 10.8 recorded before focus scan"
            f.attrs["num_rois"] = 0
            f.create_dataset("Parameters", data=np.zeros(0, dtype=[("ROI", "i4")]))

        dlg = FileBrowserDialog(str(tmp_path))
        dlg._query_edit.setText("axis = shift")

        visible = [r for r in range(dlg._proxy.rowCount()) if not dlg._table.isRowHidden(r)]
        visible_names = {dlg._proxy.index(r, 0).data() for r in visible}
        assert visible_names == {"shift.h5"}

    def test_query_combines_and_not(self, qapp, tmp_path):
        import h5py
        with h5py.File(tmp_path / "a.h5", "w") as f:
            f.attrs["axis"] = "Shift_mm"
            f.attrs["Notes"] = "#fail bad run"
            f.attrs["num_rois"] = 10
            f.create_dataset("Parameters", data=np.zeros(0, dtype=[("ROI", "i4")]))
        with h5py.File(tmp_path / "b.h5", "w") as f:
            f.attrs["axis"] = "Shift_mm"
            f.attrs["Notes"] = "#success good run"
            f.attrs["num_rois"] = 10
            f.create_dataset("Parameters", data=np.zeros(0, dtype=[("ROI", "i4")]))

        dlg = FileBrowserDialog(str(tmp_path))
        dlg._query_edit.setText("axis = shift and not tag = fail")

        visible = [r for r in range(dlg._proxy.rowCount()) if not dlg._table.isRowHidden(r)]
        visible_names = {dlg._proxy.index(r, 0).data() for r in visible}
        assert visible_names == {"b.h5"}

    def test_query_numeric_comparison(self, qapp, tmp_path):
        import h5py
        for name, n_selected in (("low.h5", 2), ("high.h5", 20)):
            with h5py.File(tmp_path / name, "w") as f:
                f.attrs["axis"] = "Shift_mm"
                f.attrs["num_rois"] = 30
                params = np.zeros(30, dtype=[("ROI", "i4"), (" selected", "f8")])
                params[" selected"][:n_selected] = 1.0
                f.create_dataset("Parameters", data=params)

        dlg = FileBrowserDialog(str(tmp_path))
        dlg._query_edit.setText("selected >= 10")

        visible = [r for r in range(dlg._proxy.rowCount()) if not dlg._table.isRowHidden(r)]
        visible_names = {dlg._proxy.index(r, 0).data() for r in visible}
        assert visible_names == {"high.h5"}

    def test_invalid_query_does_not_crash(self, qapp, tmp_path):
        import h5py
        with h5py.File(tmp_path / "a.h5", "w") as f:
            f.attrs["axis"] = "Shift_mm"
            f.attrs["num_rois"] = 0
            f.create_dataset("Parameters", data=np.zeros(0, dtype=[("ROI", "i4")]))

        dlg = FileBrowserDialog(str(tmp_path))
        dlg._query_edit.setText("axis = (")  # malformed
        # must not raise; row visibility is simply left unchanged on error


class TestFitDialogLivePreview:
    """FitDialog must preview a live fit immediately on open (for
    pipeline-driven "Fit" params) and re-fit on Default/Manual mode switch,
    not only in response to a manual param edit (regression: both used to
    require an explicit edit before _on_param_changed ever ran)."""

    _CONFIG = {
        "test_pipeline": [
            {"model": "ExponentialDrift",
             "params": {"Z0 (um)": {"vary": True}, "dZ_dt (nm_s)": {"vary": True}}},
        ],
    }

    def test_previews_fit_immediately_on_open(self, editor_with_file, monkeypatch):
        import math
        import trace_editor.gui.main_window as mw
        from trace_editor.gui.main_window import FitDialog

        win, data = editor_with_file
        monkeypatch.setattr(mw, "load_pipeline_config", lambda: self._CONFIG)
        win._active_pipeline = "test_pipeline"
        roi = win._current_roi
        assert math.isnan(data.get_parameter("dZ_dt (nm_s)", roi, default=float("nan")))

        dlg = FitDialog(win, model_name="ExponentialDrift", h5data=data,
                         roi=roi, live_preview_fn=win.redraw)
        try:
            fitted = data.get_parameter("dZ_dt (nm_s)", roi, default=float("nan"))
            assert not math.isnan(fitted), "dialog open should have run a live preview fit"
        finally:
            dlg.reject()

    def test_mode_switch_refits(self, editor_with_file, monkeypatch):
        import trace_editor.gui.main_window as mw
        from trace_editor.gui.main_window import FitDialog

        win, data = editor_with_file
        monkeypatch.setattr(mw, "load_pipeline_config", lambda: self._CONFIG)
        win._active_pipeline = "test_pipeline"
        roi = win._current_roi

        dlg = FitDialog(win, model_name="ExponentialDrift", h5data=data,
                         roi=roi, live_preview_fn=win.redraw)
        try:
            calls = []
            orig = dlg._on_param_changed
            monkeypatch.setattr(dlg, "_on_param_changed", lambda: (calls.append(None), orig())[-1])

            dlg._mode_tabs.setCurrentIndex(1)  # Default -> Manual
            assert calls, "switching mode should re-run the live preview fit"
        finally:
            dlg.reject()


class TestFitPipelines:
    """Pipeline selector, Run button, auto-fit-on-open toggle, and the
    in-app folder watcher (docs/architecture.md "Fit pipelines").
    """

    _CONFIG = {
        "pipeline_a": [
            {"model": "ExponentialDrift", "criteria": [{"channel": "Shift (mm)", "min": "p90"}]},
        ],
        "pipeline_b": [
            {"model": "ExponentialDrift", "criteria": [{"channel": "Shift (mm)", "min": "p10"}]},
        ],
    }

    def test_pipeline_combo_populated_from_config(self, qapp, editor, monkeypatch):
        import trace_editor.gui.main_window as mw
        monkeypatch.setattr(mw, "load_pipeline_config", lambda: self._CONFIG)
        # Population happens in __init__, so build a fresh window against the patch.
        win = mw.TraceEditor()
        try:
            names = {win.pipelineComboBox.itemText(i) for i in range(win.pipelineComboBox.count())}
            assert names == {"pipeline_a", "pipeline_b"}
        finally:
            win.close()

    def test_switching_pipeline_changes_default_selection_criteria(self, qapp, monkeypatch):
        import trace_editor.gui.main_window as mw
        monkeypatch.setattr(mw, "load_pipeline_config", lambda: self._CONFIG)
        win = mw.TraceEditor()
        try:
            win.pipelineComboBox.setCurrentText("pipeline_a")
            assert win._active_pipeline == "pipeline_a"
        finally:
            win.close()

    def test_run_button_fits_open_files_with_active_pipeline(self, qapp, qtbot, editor_with_file, monkeypatch):
        import trace_editor.gui.main_window as mw
        win, data = editor_with_file
        config = {
            "test_pipeline": [
                {"model": "ExponentialDrift",
                 "params": {"Z0 (um)": {"vary": True}, "dZ_dt (nm_s)": {"vary": True}}},
            ],
        }
        monkeypatch.setattr(mw, "load_pipeline_config", lambda: config)
        win._active_pipeline = "test_pipeline"

        win._on_run_pipeline_clicked()

        def done():
            assert not np.isnan(data.get_parameter("dZ_dt (nm_s)", win._current_roi, default=float("nan")))

        qtbot.waitUntil(done, timeout=5000)

    def test_auto_fit_on_open_runs_pipeline_on_newly_opened_file(self, qapp, qtbot, editor, generated_workflow_h5, monkeypatch):
        import trace_editor.gui.main_window as mw
        config = {
            "test_pipeline": [
                {"model": "ExponentialDrift",
                 "params": {"Z0 (um)": {"vary": True}, "dZ_dt (nm_s)": {"vary": True}}},
            ],
        }
        monkeypatch.setattr(mw, "load_pipeline_config", lambda: config)
        editor._active_pipeline = "test_pipeline"
        editor._auto_fit_on_open = True

        assert editor.open_single_file(str(generated_workflow_h5))
        editor._activate_opened_file(str(generated_workflow_h5))

        from trace_editor.core.h5_data import H5Data
        data = H5Data(str(generated_workflow_h5))

        def done():
            assert not np.isnan(data.get_parameter("dZ_dt (nm_s)", "000", default=float("nan")))

        qtbot.waitUntil(done, timeout=5000)

    def test_auto_fit_on_open_disabled_by_default(self, editor_with_file):
        win, data = editor_with_file
        assert win._auto_fit_on_open is False

    def test_redraw_survives_criteria_referencing_missing_channel(self, editor_with_file, monkeypatch):
        """A pipeline step's criteria can reference a channel the currently
        viewed file doesn't have (e.g. switched axis, imported file) --
        redraw() must degrade gracefully (status message, no highlight),
        not crash (regression: this used to propagate out of a QThread
        finished-callback and corrupt Qt's thread bookkeeping)."""
        import trace_editor.gui.main_window as mw
        win, data = editor_with_file
        config = {
            "test_pipeline": [
                {"model": "ExponentialDrift",
                 "criteria": [{"channel": "No Such Channel (mm)", "min": "p90"}]},
            ],
        }
        monkeypatch.setattr(mw, "load_pipeline_config", lambda: config)
        win._active_pipeline = "test_pipeline"
        win.modelComboBox.setCurrentText("ExponentialDrift")
        assert win.modelComboBox.currentText() == "ExponentialDrift"

        win.redraw()   # must not raise

        assert "selection error" in win.statusbar.currentMessage()

    def test_watcher_opens_new_file_dropped_into_watched_folder(self, qapp, qtbot, editor, generated_workflow_h5, tmp_path):
        import shutil
        watch_dir = tmp_path / "watched"
        watch_dir.mkdir()
        editor._watch_folder = str(watch_dir)
        editor._start_watching()

        dest = watch_dir / "new_trace.h5"
        shutil.copy(generated_workflow_h5, dest)

        def opened():
            assert str(dest) in editor._open_files

        qtbot.waitUntil(opened, timeout=6000)
        editor._watch_timer.stop()

    def test_watcher_ignores_pre_existing_files_at_start(self, qapp, editor, generated_workflow_h5, tmp_path):
        import shutil
        watch_dir = tmp_path / "watched2"
        watch_dir.mkdir()
        pre_existing = watch_dir / "already_here.h5"
        shutil.copy(generated_workflow_h5, pre_existing)

        editor._watch_folder = str(watch_dir)
        editor._start_watching()

        assert str(pre_existing) not in editor._open_files
        editor._watch_timer.stop()
