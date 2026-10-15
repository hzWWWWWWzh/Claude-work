import sys
import math
import re
import tomllib
from pathlib import Path

from PySide6.QtWidgets import (
    QApplication, QMainWindow, QVBoxLayout, QHBoxLayout, QFileDialog,
    QListWidgetItem, QDialog, QLabel, QLineEdit, QComboBox,
    QDialogButtonBox, QPushButton, QTableWidget, QTableWidgetItem,
    QHeaderView, QAbstractItemView, QSizePolicy, QMessageBox, QTabBar,
    QTableView,
)
from PySide6.QtCore import (
    Qt, QSettings, QThread, QObject, QTimer, Signal, QUrl, QSortFilterProxyModel,
    QItemSelectionModel,
)
from PySide6.QtGui import QKeySequence, QShortcut, QImage, QDesktopServices, QIcon, QStandardItemModel, QStandardItem, QActionGroup
import io

import matplotlib
matplotlib.use("QtAgg")
import matplotlib.pyplot as plt
from matplotlib.figure import Figure
from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg as FigureCanvas

import numpy as np

import warnings
warnings.filterwarnings("ignore", message=".*pkg_resources.*")

from trace_editor.core.h5_data import H5Data, prepare_fit, apply_fit, r2_label, r2_full_label
from trace_editor.core.batch_fit import (
    _exponential_drift, _correct_z_for_fit, _get_fit_selection,
    _store_fit_result, run_batch, _apply_point_selection,
)
from trace_editor.core.auto_fit_pipeline import (
    DEFAULT_PIPELINE_NAME, list_pipeline_names, resolve_params,
    resolve_param_sources, resolve_criteria, resolve_force_defaults, run_pipeline,
)
from trace_editor.core.auto_fit_pipeline import load_config as load_pipeline_config
from trace_editor.core.force_calc import calculate_force_for_file
from trace_editor.core.inspect_h5 import (
    open_in_vscode_h5web, InspectH5Error, MYHDF5_FALLBACK_URL,
)
from trace_editor.core.file_browser import scan_folder, extract_tags
from trace_editor.core.file_query import parse_query, evaluate_query, QueryError
from trace_editor.gui.ui_main_window import Ui_MainWindow


_APP_CONFIG_FILE = Path(__file__).parent / "app_config.toml"
_APP_CONFIG_DEFAULT_FILE = Path(__file__).parent / "app_config.default.toml"
plt.style.use(Path(__file__).parent / "traceeditor.mplstyle")

_THEMES_FILE = Path(__file__).parent / "themes.toml"
with open(_THEMES_FILE, "rb") as _f:
    _THEMES = tomllib.load(_f)

# Selected-point scatter color palette, in tab10 order (cyan/green swapped).
_PALETTE_COLORS = {
    "theme":     "#226b67",  # App primary theme color (default)
    "blue":      "#4069E0",  # Royal Blue
    "orange":    "#FD8C00",  # Orange
    "green":     "#4BAB6E",  # Go Green
    "pink":      "#FF58B3",  # Hot Pink
    "lightblue": "#5DCBF3",  # Sky Blue
    "purple":    "#D16EFF",  # Parmaviolet
    "fuchia":    "#EA2353",  # Crimson Pink
    "yellow":    "#FFE066",  # Soft Yellow
    "teal":      "#35CBC4",  # Bright Teal
    "red":       "#F14124",  # Red-Orange
}


def _apply_trace_selection(criteria: list, h5data: H5Data, roi: str) -> bool:
    """Apply trace-level selection criteria. Returns True if ROI passes all criteria.

    Mirrors TraceIO.select_trace(). Each criterion: {Parameter, min, max, invert}.
    """
    stored = h5data.list_parameters(roi)
    total = 1
    for c in criteria:
        name = c.get("Parameter", "")
        val = stored.get(name, float("nan"))
        if math.isnan(val):
            continue
        selected = 1
        lo = c.get("min", "")
        hi = c.get("max", "")
        if lo != "":
            try:
                selected *= int(val > float(lo))
            except ValueError:
                pass
        if hi != "":
            try:
                selected *= int(val < float(hi))
            except ValueError:
                pass
        if c.get("invert", False):
            selected = 1 - selected
        total *= selected
    h5data.set_selection(roi, float(total))
    return bool(total)


# ---------------------------------------------------------------------------
# Open File browser dialog
# ---------------------------------------------------------------------------

_BROWSER_COLUMNS = [
    "Filename", "Axis", "Tags", "Selected/Total", "LUT", "Notes", "Date",
]
_NOTES_COL = _BROWSER_COLUMNS.index("Notes")

_FILE_NUM_RE = re.compile(r"(\d+)")


def _file_sort_key(path: str):
    """Sort key grouping files by parent folder name, then by the numeric
    index in the filename (e.g. `data_042.h5` -> 42), so files from the same
    folder stay together and in `data_XXX` order."""
    p = Path(path)
    m = _FILE_NUM_RE.search(p.stem)
    num = int(m.group(1)) if m else -1
    return (p.parent.name, num, p.name)


def _display_label(path: str) -> str:
    """`"Folder/filename"` label for the file combo box, e.g.
    `20260710/data_003.h5`."""
    p = Path(path)
    return f"{p.parent.name}/{p.name}"


class FileBrowserDialog(QDialog):
    """Open File dialog: a sortable/queryable table of `.h5` metadata for one
    folder, so files can be triaged (axis, tags, selected/total beads, LUT
    file used) before opening them.

    Populated via `core.file_browser.scan_folder`, which reads only root
    attributes, `/Settings`, and the small `/Parameters` array — no trace data
    is loaded.
    """

    def __init__(self, folder: str, parent=None, open_paths: set[str] | None = None):
        super().__init__(parent)
        self.setWindowTitle("Open h5 files")
        self.resize(820, 420)
        self.selected_paths = []
        self.closed_paths = []
        self._open_paths = open_paths or set()
        self._rows = []  # raw scan_folder() rows, keyed by source model row

        self._folder_edit = QLineEdit(folder)
        browse_btn = QPushButton("Browse…")
        browse_btn.clicked.connect(self._browse)
        folder_row = QHBoxLayout()
        folder_row.addWidget(QLabel("Folder:"))
        folder_row.addWidget(self._folder_edit)
        folder_row.addWidget(browse_btn)

        self._query_edit = QLineEdit()
        self._query_edit.setPlaceholderText(
            "Query, e.g. axis = shift and not tag = fail and selected >= 10")
        self._query_edit.textChanged.connect(self._on_query_changed)
        filter_row = QHBoxLayout()
        filter_row.addWidget(self._query_edit)

        self._model = QStandardItemModel(0, len(_BROWSER_COLUMNS), self)
        self._model.setHorizontalHeaderLabels(_BROWSER_COLUMNS)
        self._proxy = QSortFilterProxyModel(self)
        self._proxy.setSourceModel(self._model)

        self._table = QTableView(self)
        self._table.setModel(self._proxy)
        self._table.setSortingEnabled(True)
        self._table.setSelectionBehavior(QAbstractItemView.SelectRows)
        self._table.setSelectionMode(QAbstractItemView.ExtendedSelection)
        self._table.setEditTriggers(QAbstractItemView.NoEditTriggers)
        self._table.horizontalHeader().setSectionResizeMode(_NOTES_COL, QHeaderView.Stretch)
        self._table.doubleClicked.connect(self._on_ok)

        btns = QDialogButtonBox(QDialogButtonBox.Open | QDialogButtonBox.Cancel)
        btns.button(QDialogButtonBox.Open).clicked.connect(self._on_ok)
        btns.rejected.connect(self.reject)

        lay = QVBoxLayout(self)
        lay.addLayout(folder_row)
        lay.addLayout(filter_row)
        lay.addWidget(self._table)
        lay.addWidget(btns)

        self._rescan()
        self._folder_edit.editingFinished.connect(self._rescan)

    def _browse(self):
        folder = QFileDialog.getExistingDirectory(self, "Choose folder", self._folder_edit.text())
        if folder:
            self._folder_edit.setText(folder)
            self._rescan()

    def _rescan(self):
        folder = self._folder_edit.text()
        self._model.setRowCount(0)
        self._rows = []
        if not folder or not Path(folder).is_dir():
            return
        for row in scan_folder(folder):
            self._rows.append(row)
            notes_item = QStandardItem(row["error"] or row["notes"])
            notes_item.setToolTip(row["error"] or row["notes"])
            items = [
                QStandardItem(row["filename"]),
                QStandardItem(row["axis"]),
                QStandardItem(" ".join(f"#{t}" for t in row["tags"])),
                QStandardItem(f"{row['selected']}/{row['num_rois']}"),
                QStandardItem(row["lut_filename"]),
                notes_item,
                QStandardItem(row["date"]),
            ]
            items[0].setData(row["path"], Qt.UserRole)
            for item in items:
                item.setEditable(False)
            self._model.appendRow(items)
        self._table.resizeColumnsToContents()
        self._table.horizontalHeader().setSectionResizeMode(_NOTES_COL, QHeaderView.Stretch)
        self._select_open_rows()

    def _select_open_rows(self):
        if not self._open_paths:
            return
        selection = self._table.selectionModel()
        for src_row, row in enumerate(self._rows):
            if row["path"] in self._open_paths:
                proxy_idx = self._proxy.mapFromSource(self._model.index(src_row, 0))
                selection.select(
                    proxy_idx,
                    QItemSelectionModel.Select | QItemSelectionModel.Rows)

    def _on_query_changed(self, text):
        try:
            node = parse_query(text)
        except QueryError:
            self._query_edit.setStyleSheet("QLineEdit { background: #5a2a2a; }")
            return
        self._query_edit.setStyleSheet("")
        for src_row, row in enumerate(self._rows):
            try:
                match = evaluate_query(node, row)
            except QueryError:
                match = True
            self._table.setRowHidden(
                self._proxy.mapFromSource(self._model.index(src_row, 0)).row(), not match)

    def _on_ok(self):
        rows = {i.row() for i in self._table.selectionModel().selectedIndexes()}
        paths = []
        for row in rows:
            src_row = self._proxy.mapToSource(self._proxy.index(row, 0)).row()
            paths.append(self._model.item(src_row, 0).data(Qt.UserRole))
        self.selected_paths = paths
        # Already-open files in this folder that were deselected -> close them.
        folder_open_paths = {row["path"] for row in self._rows} & self._open_paths
        self.closed_paths = list(folder_open_paths - set(paths))
        if paths or self.closed_paths:
            self.accept()


# ---------------------------------------------------------------------------
# Select Points dialog
# ---------------------------------------------------------------------------

class SelectPointsDialog(QDialog):
    """Channel min/max/derivative/invert table with live preview on current ROI."""

    def __init__(self, h5data: H5Data, roi: str, live_fn, parent=None):
        super().__init__(parent)
        self.setWindowTitle("Select points …")
        self._h5data = h5data
        self._roi = roi
        self._live_fn = live_fn
        self.result_token = "Cancel"
        self.criteria = []

        channels = h5data.all_channels()
        self._channels = channels

        # Table: channel | min | max | derivative | invert
        self._table = QTableWidget(len(channels), 5, self)
        self._table.setHorizontalHeaderLabels(
            ["Channel", "Min", "Max", "Derivative", "Invert"])
        self._table.horizontalHeader().setSectionResizeMode(0, QHeaderView.ResizeToContents)
        self._table.horizontalHeader().setSectionResizeMode(QHeaderView.Stretch)
        self._table.setSelectionMode(QAbstractItemView.NoSelection)

        for i, ch in enumerate(channels):
            self._table.setItem(i, 0, QTableWidgetItem(ch))
            self._table.item(i, 0).setFlags(Qt.ItemIsEnabled)
            self._table.setItem(i, 1, QTableWidgetItem(""))
            self._table.setItem(i, 2, QTableWidgetItem(""))
            deriv = QComboBox()
            deriv.addItems(["False", "True"])
            deriv.setCurrentText("True" if ch == "Shift (mm)" else "False")
            self._table.setCellWidget(i, 3, deriv)
            inv = QComboBox()
            inv.addItems(["False", "True"])
            self._table.setCellWidget(i, 4, inv)

        self._table.cellChanged.connect(self._on_changed)

        btns = QDialogButtonBox()
        self._btn_cancel = btns.addButton("Cancel", QDialogButtonBox.RejectRole)
        self._btn_batch  = btns.addButton("Batch",  QDialogButtonBox.ActionRole)
        self._btn_ok     = btns.addButton("OK",     QDialogButtonBox.AcceptRole)
        self._btn_cancel.clicked.connect(self._on_cancel)
        self._btn_batch.clicked.connect(self._on_batch)
        self._btn_ok.clicked.connect(self._on_ok)

        lay = QVBoxLayout(self)
        lay.addWidget(self._table)
        lay.addWidget(btns)
        self.resize(560, 400)

    def _build_criteria(self):
        rows = []
        for i, ch in enumerate(self._channels):
            lo = (self._table.item(i, 1).text() if self._table.item(i, 1) else "").strip()
            hi = (self._table.item(i, 2).text() if self._table.item(i, 2) else "").strip()
            deriv = self._table.cellWidget(i, 3).currentText() == "True"
            inv   = self._table.cellWidget(i, 4).currentText() == "True"
            rows.append({"channel": ch, "min": lo, "max": hi,
                         "derivative": deriv, "invert": inv})
        return rows

    def _on_changed(self):
        self.criteria = self._build_criteria()
        _apply_point_selection(self.criteria, self._h5data, self._roi)
        if self._live_fn:
            self._live_fn()

    def _on_cancel(self):
        self.result_token = "Cancel"
        self.criteria = []
        self.reject()

    def _on_batch(self):
        self.criteria = [c for c in self._build_criteria()
                         if c["min"] != "" or c["max"] != ""]
        self.result_token = "Batch"
        self.accept()

    def _on_ok(self):
        self.result_token = "OK"
        self.criteria = []
        self.accept()


# ---------------------------------------------------------------------------
# Select Traces dialog
# ---------------------------------------------------------------------------

class SelectTracesDialog(QDialog):
    """Parameter min/max/invert table; OK applies to all selected ROIs."""

    def __init__(self, h5data: H5Data, roi: str, parent=None):
        super().__init__(parent)
        self.setWindowTitle("Select traces …")
        self.criteria = []
        self.result_token = "Cancel"

        stored = h5data.list_parameters(roi)
        # Only show params with non-NaN values and no internal prefix
        self._params = [k for k, v in stored.items()
                        if not math.isnan(v)
                        and not k.startswith("_")
                        and not k.startswith(" ")]

        self._table = QTableWidget(len(self._params), 4, self)
        self._table.setHorizontalHeaderLabels(["Parameter", "Min", "Max", "Invert"])
        self._table.horizontalHeader().setSectionResizeMode(0, QHeaderView.ResizeToContents)
        self._table.horizontalHeader().setSectionResizeMode(QHeaderView.Stretch)
        self._table.setSelectionMode(QAbstractItemView.NoSelection)

        for i, name in enumerate(self._params):
            self._table.setItem(i, 0, QTableWidgetItem(name))
            self._table.item(i, 0).setFlags(Qt.ItemIsEnabled)
            self._table.setItem(i, 1, QTableWidgetItem(""))
            self._table.setItem(i, 2, QTableWidgetItem(""))
            inv = QComboBox()
            inv.addItems(["False", "True"])
            self._table.setCellWidget(i, 3, inv)

        btns = QDialogButtonBox(QDialogButtonBox.Ok | QDialogButtonBox.Cancel)
        btns.accepted.connect(self._on_ok)
        btns.rejected.connect(self.reject)

        lay = QVBoxLayout(self)
        lay.addWidget(self._table)
        lay.addWidget(btns)
        self.resize(480, 380)

    def _on_ok(self):
        self.criteria = []
        for i, name in enumerate(self._params):
            lo = (self._table.item(i, 1).text() if self._table.item(i, 1) else "").strip()
            hi = (self._table.item(i, 2).text() if self._table.item(i, 2) else "").strip()
            inv = self._table.cellWidget(i, 3).currentText() == "True"
            if lo != "" or hi != "":
                self.criteria.append({"Parameter": name, "min": lo, "max": hi, "invert": inv})
        self.result_token = "OK"
        self.accept()


TYPE_OPTIONS = ["Global", "Local", "Fit"]


# ---------------------------------------------------------------------------
# Fit Parameters Dialog  (mirrors original SetFitParametersWindow)
# ---------------------------------------------------------------------------

class FitDialog(QDialog):
    """Modal dialog for setting fit parameters.

    Columns: Name | Value | Type (Global/Local/Fit) | Min | Max
    Buttons: Cancel | Fit Current | Fit All

    On Fit Current → caller fits current ROI with these params.
    On Fit All → caller fits all selected ROIs with these params.
    On Cancel (or closing via the window's X) → live-preview writes made
    while the dialog was open are reverted, restoring the H5 parameters
    for this ROI to their pre-dialog state.

    `live_preview_fn` is called whenever any value/type changes so the
    main plot updates in real time (matching original apply_model behaviour).
    """

    # Result token set before close
    result_token: str = "Cancel"   # "OK" | "Batch" | "Cancel"

    def __init__(self, parent, model_name: str, h5data: H5Data, roi: str,
                 live_preview_fn=None):
        super().__init__(parent)
        self.setWindowTitle("Set Fit Parameters …")
        self.setModal(True)
        self._live_preview_fn = live_preview_fn
        # `parent` is the TraceEditor main window; selection mode is shared
        # session state (see TraceEditor.__init__), read/written here so the
        # Default/Manual toggle and the main plot's highlight always agree.
        self._main_window = parent

        # Build lmfit params with initial values seeded from the active
        # pipeline's own step for this model, falling back to the model's
        # bare docstring default -- see core.auto_fit_pipeline.resolve_params.
        config = load_pipeline_config()
        try:
            self._model, self.params, self.names = resolve_params(
                model_name, self._main_window._active_pipeline, config)
        except Exception as e:
            raise ValueError(f"Cannot prepare model '{model_name}': {e}") from e

        # Point-selection criteria for "Default" mode, resolved once here
        # (same pattern as params above) so the fit-execution layer
        # (_get_fit_selection/run_batch/_fit_one_roi) only ever sees a plain
        # criteria list, never the active pipeline itself.
        self._criteria = resolve_criteria(model_name, self._main_window._active_pipeline, config)

        self._model_name = model_name
        self._h5data = h5data
        self._roi = roi

        # Seed type + value from stored parameters (mirrors original logic)
        stored = h5data.list_parameters(roi)
        self._pre_dialog_params = dict(stored)  # snapshot for revert on Cancel/close

        # Zbead baseline for live-preview correction/accumulation, frozen at
        # dialog-open time. _correct_z_for_fit/_store_fit_result normally
        # read+accumulate "Zbead (um)" straight from h5, which is safe only
        # when that slot isn't touched between reads. Live preview edits the
        # z_bead value field, and if that edit were written to h5 before the
        # next preview fit, the "stored baseline" the correction reads would
        # actually be the user's typed number, not the true accumulated
        # baseline -- causing the offset to compound on every edit instead of
        # settling like the first fit does. Freezing the baseline here and
        # threading it through every live-preview call keeps preview fits
        # anchored to the pre-dialog value no matter how many edits happen;
        # only committing (Fit Current/Fit All) advances the real baseline.
        zbead_baseline = stored.get("Zbead (um)", 0.0)
        self._zbead_baseline = 0.0 if math.isnan(zbead_baseline) else zbead_baseline

        # Type = where this param's value comes from, per the active
        # pipeline's step for this model (see resolve_param_sources): "Fit"
        # (pipeline says vary=true), "Global" (pipeline sets a fixed value,
        # shared across every trace), or "Local" (the pipeline doesn't
        # mention this param -- use this ROI's own stored value, or the
        # model's bare default if nothing is stored yet).
        sources = resolve_param_sources(model_name, self._main_window._active_pipeline, config)
        for p in self.params:
            label = self.names[p]
            if label in ("x", "y"):
                continue
            source = sources.get(label, "Local")
            # Stored in lmfit's `user_data` slot (not a bare custom attribute)
            # so it survives `copy.deepcopy`/pickling across the batch-fit
            # process boundary -- see core.batch_fit._fit_one_roi.
            self.params[p].user_data = source
            if source == "Fit":
                err = stored.get(f"{label} Error", float("nan"))
                if not math.isnan(err) and err > 0:
                    self.params[p].stderr = err
            elif source == "Local":
                val = stored.get(label, float("nan"))
                if not math.isnan(val):
                    self.params[p].value = val
                self.params[p].vary = False

        self._build_ui()
        self._on_param_changed()   # initial live preview -- see _on_param_changed

    def _build_ui(self):
        layout = QVBoxLayout(self)

        # Default/Manual point-selection toggle (see TraceEditor._selection_mode).
        # QTabBar rather than checkable QPushButtons: macOS's native style
        # doesn't render checkable QPushButton state reliably (the "checked"
        # look is indistinguishable from ordinary focus/default-button
        # styling and doesn't move when the selection changes), while tab
        # selection has unambiguous, platform-consistent active-state styling.
        self._mode_tabs = QTabBar()
        self._mode_tabs.addTab("Default selection")
        self._mode_tabs.addTab("Manual selection")
        self._mode_tabs.setCurrentIndex(1 if self._main_window._selection_mode == "manual" else 0)
        self._mode_tabs.currentChanged.connect(self._on_mode_toggled)
        layout.addWidget(self._mode_tabs)

        # Table: Name | Value | Type | Min | Max
        param_keys = [p for p in self.params if self.names[p] not in ("x", "y")]
        self._param_keys = param_keys

        tbl = QTableWidget(len(param_keys), 5, self)
        tbl.setHorizontalHeaderLabels(["Parameter", "Value", "Type", "Min", "Max"])
        tbl.horizontalHeader().setSectionResizeMode(0, QHeaderView.ResizeMode.Stretch)
        for col in (1, 2, 3, 4):
            tbl.horizontalHeader().setSectionResizeMode(col, QHeaderView.ResizeMode.ResizeToContents)
        tbl.setSelectionMode(QAbstractItemView.SelectionMode.NoSelection)
        tbl.verticalHeader().setVisible(False)
        self._tbl = tbl

        self._value_edits: list[QLineEdit] = []
        self._type_combos: list[QComboBox] = []
        self._min_edits: list[QLineEdit] = []
        self._max_edits: list[QLineEdit] = []

        for row, p in enumerate(param_keys):
            label = self.names[p]
            prm = self.params[p]

            # Name (read-only)
            name_item = QTableWidgetItem(label)
            name_item.setFlags(name_item.flags() & ~Qt.ItemFlag.ItemIsEditable)
            tbl.setItem(row, 0, name_item)

            # Value
            val_str = f"{prm.value:.6g}" if not math.isnan(float(prm.value)) else ""
            val_edit = QLineEdit(val_str)
            val_edit.editingFinished.connect(self._on_param_changed)
            tbl.setCellWidget(row, 1, val_edit)
            self._value_edits.append(val_edit)

            # Type
            type_cb = QComboBox()
            type_cb.addItems(TYPE_OPTIONS)
            type_cb.setCurrentText(prm.user_data or "Global")
            type_cb.currentIndexChanged.connect(self._on_param_changed)
            tbl.setCellWidget(row, 2, type_cb)
            self._type_combos.append(type_cb)

            # Min
            min_str = f"{prm.min:.6g}" if math.isfinite(prm.min) else ""
            min_edit = QLineEdit(min_str)
            min_edit.editingFinished.connect(self._on_param_changed)
            tbl.setCellWidget(row, 3, min_edit)
            self._min_edits.append(min_edit)

            # Max
            max_str = f"{prm.max:.6g}" if math.isfinite(prm.max) else ""
            max_edit = QLineEdit(max_str)
            max_edit.editingFinished.connect(self._on_param_changed)
            tbl.setCellWidget(row, 4, max_edit)
            self._max_edits.append(max_edit)

        layout.addWidget(tbl)

        # Buttons: Cancel | Fit Current | Fit All
        btn_layout = QHBoxLayout()
        cancel_btn = QPushButton("Cancel")
        ok_btn     = QPushButton("Fit Current")
        batch_btn  = QPushButton("Fit All")
        for btn in (cancel_btn, ok_btn, batch_btn):
            btn.setSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Fixed)
            btn.setAutoDefault(False)
            btn_layout.addWidget(btn)
        batch_btn.setAutoDefault(True)
        batch_btn.setDefault(True)
        cancel_btn.clicked.connect(self._on_cancel)
        batch_btn.clicked.connect(self._on_batch)
        ok_btn.clicked.connect(self._on_ok)
        layout.addLayout(btn_layout)

        # Fit every row with no scrolling. Rather than compute a pixel
        # height by hand (which drifted out of sync once before, and is
        # fragile across platforms/DPI/font metrics -- exactly what
        # happened here), cap the table's height at its real content height
        # and let Qt's own layout system size the dialog around it via
        # adjustSize(). QTableWidget defaults to an expanding size policy
        # (it assumes it lives in a scroll area), so without this the
        # layout has no correct height to report and just clips the table
        # to whatever resize() guessed.
        tbl.resizeRowsToContents()
        content_height = tbl.horizontalHeader().height() + tbl.frameWidth() * 2
        for row in range(tbl.rowCount()):
            content_height += tbl.rowHeight(row)
        tbl.setMaximumHeight(content_height)
        tbl.setMinimumHeight(content_height)

        self.setMinimumWidth(600)
        self.adjustSize()

    def _on_mode_toggled(self, index: int):
        """Default/Manual toggle: update shared session state, re-run the
        live preview fit under the new selection, and redraw.

        Switching modes only changes which mask is highlighted/used — it
        never writes anything (the "manual" mask is whatever's already
        persisted at Traces/Selection/<roi>, set independently via
        SelectPointsDialog/"Select points in view"). Goes through
        `_on_param_changed` (not just `_live_preview_fn`) so any "Fit" params
        actually re-fit against the newly selected points, not just redraw
        with stale values from before the switch.
        """
        self._main_window._selection_mode = "manual" if index == 1 else "default"
        self._on_param_changed()

    def _read_ui_into_params(self):
        """Pull current widget values into self.params."""
        for row, p in enumerate(self._param_keys):
            prm = self.params[p]
            try:
                prm.value = float(self._value_edits[row].text())
            except ValueError:
                pass
            try:
                prm.min = float(self._min_edits[row].text())
            except ValueError:
                prm.min = -np.inf
            try:
                prm.max = float(self._max_edits[row].text())
            except ValueError:
                prm.max = np.inf

            typ = self._type_combos[row].currentText()
            prm.user_data = typ
            prm.vary = (typ == "Fit")

    def _on_param_changed(self):
        """Live preview: read params, run fit if any vary=True, update plot."""
        self._read_ui_into_params()

        # Write all current values to h5 so the live preview uses them, except
        # z_bead: its h5 slot ("Zbead (um)") is the accumulated baseline that
        # _correct_z_for_fit/_store_fit_result read and add to, not a plain
        # live value. Writing the dialog's in-progress edit there would
        # corrupt that baseline for every subsequent preview fit (see
        # `_zbead_baseline` comment in __init__). z_bead is exempt from the
        # displayed model curve anyway (`_draw_model_curve` forces it to 0),
        # so skipping this write only holds the data's display shift at its
        # pre-dialog value until the edit is committed via Fit Current/All.
        for p in self._param_keys:
            if self.names[p] == "z_bead":
                continue
            self._h5data.set_parameter(self.names[p], float(self.params[p].value),
                                        roi=self._roi)

        # If any params are set to Fit, run the fit immediately for live preview
        any_vary = any(self.params[p].vary for p in self._param_keys)
        if any_vary:
            try:
                x = self._h5data.get_channel(self.names["x"], self._roi)
                y = self._h5data.get_channel(self.names["y"], self._roi)
                y = _correct_z_for_fit(y, self._model_name, self.names,
                                        self._h5data, self._roi,
                                        zbead_baseline=self._zbead_baseline)
                sel = _get_fit_selection(self._h5data, self._roi,
                                          mode=self._main_window._selection_mode,
                                          criteria=self._criteria)
                result, out_names = apply_fit(
                    self._model_name, x, y, selection=sel, initial_parameters=self.params)
                if result.success:
                    for row, p in enumerate(self._param_keys):
                        if self.params[p].vary:
                            fitted_val = float(result.params[p].value)
                            stderr = result.params[p].stderr or 0.0
                            self.params[p].value = fitted_val
                            self.params[p].stderr = stderr
                            # Update value edit to show fitted result
                            self._value_edits[row].blockSignals(True)
                            self._value_edits[row].setText(f"{fitted_val:.6g}")
                            self._value_edits[row].blockSignals(False)
                            # Write fit result + error to h5 (z_bead is accumulated
                            # onto the frozen pre-dialog baseline, not onto
                            # whatever h5 currently holds)
                            _store_fit_result(p, self.names, result,
                                              self._model_name, self._h5data, self._roi,
                                              zbead_baseline=self._zbead_baseline)
            except Exception:
                pass

        if self._live_preview_fn:
            self._live_preview_fn()

    def _revert_live_preview(self):
        """Restore H5 parameters for this ROI to their pre-dialog snapshot,
        undoing any writes made by live preview while the dialog was open."""
        for name, val in self._pre_dialog_params.items():
            if math.isnan(val):
                self._h5data.delete_parameter(name, self._roi)
            else:
                self._h5data.set_parameter(name, val, roi=self._roi)

    def reject(self):
        self._revert_live_preview()
        self.result_token = "Cancel"
        super().reject()

    def _on_cancel(self):
        self.reject()

    def _on_batch(self):
        self._read_ui_into_params()
        self.result_token = "Batch"
        self.accept()

    def _on_ok(self):
        self._read_ui_into_params()
        self.result_token = "OK"
        self.accept()

    def get_params(self):
        return self.params

    def get_names(self):
        return self.names

    def get_criteria(self):
        return self._criteria


# ---------------------------------------------------------------------------
# Batch fitting worker (runs in a QThread)
# ---------------------------------------------------------------------------

class FitWorker(QObject):
    """Runs a batch fit on a background QThread by delegating to
    `core.batch_fit.run_batch`, which parallelizes the per-ROI fits across a
    `ProcessPoolExecutor`.

    This class is a thin Qt adapter: it owns no fitting logic itself, only
    translates `run_batch`'s plain `progress_cb`/`stop_cb`/`error_cb`
    callables into Qt signals and the `stop_requested` flag, so the QThread
    plumbing in `TraceEditor._start_batch` (and the crash-safety fixes in
    `_clear_batch_thread`/`closeEvent`) doesn't need to change.
    """
    progress = Signal(int, int)
    finished = Signal()
    error = Signal(str)

    def __init__(self, tasks, model_name, open_files, lmfit_params, names,
                 selection_mode="default", criteria=None):
        super().__init__()
        self.tasks = tasks
        self.model_name = model_name
        self.open_files = open_files   # unused by run_batch (kept for API compat)
        self.lmfit_params = lmfit_params   # already-configured lmfit Parameters
        self.names = names
        self.selection_mode = selection_mode
        self.criteria = criteria   # already resolved by the caller (FitDialog)
        self.stop_requested = False

    def run(self):
        run_batch(
            self.tasks, self.model_name, self.lmfit_params, self.names,
            progress_cb=lambda done, total: self.progress.emit(done, total),
            stop_cb=lambda: self.stop_requested,
            error_cb=lambda err: self.error.emit(str(err)),
            selection_mode=self.selection_mode,
            criteria=self.criteria,
        )
        self.finished.emit()


class PipelineWorker(QObject):
    """Runs a `core.auto_fit_pipeline.run_pipeline` call on a background
    `QThread` for one or more files, one after another.

    Thin Qt adapter, same shape as `FitWorker` — `run_pipeline`'s
    `progress_cb`/`error_cb`/`status_cb` are plain callables, translated here
    into Qt signals so `TraceEditor._run_pipeline_on_files` can reuse the
    same thread-lifecycle handling `_start_batch` already has.
    """
    progress = Signal(int, int)
    status = Signal(str, int, int, int, int)   # model_name, step, total_steps, file, total_files
    finished = Signal()
    error = Signal(str)

    def __init__(self, pipeline_name, filepaths, config):
        super().__init__()
        self.pipeline_name = pipeline_name
        self.filepaths = filepaths
        self.config = config

    def run(self):
        for i, path in enumerate(self.filepaths):
            try:
                run_pipeline(
                    self.pipeline_name, path, self.config,
                    progress_cb=lambda done, total: self.progress.emit(done, total),
                    error_cb=lambda err: self.error.emit(str(err)),
                    status_cb=lambda model_name, step, total_steps, i=i: self.status.emit(
                        model_name, step, total_steps, i + 1, len(self.filepaths)),
                )
            except Exception as e:
                self.error.emit(str(e))
        self.progress.emit(len(self.filepaths), len(self.filepaths))
        self.finished.emit()


# ---------------------------------------------------------------------------
# Main window
# ---------------------------------------------------------------------------

class TraceEditor(QMainWindow, Ui_MainWindow):
    """Main application window.

    A 3-column PySide6 window (model/params | plot | axis controls) over one or
    more open `H5Data` files. Holds the open files, the current file/ROI
    selection, and the matplotlib canvas. Single-ROI fitting delegates to
    `core.h5_data.apply_fit` directly (`_run_fit`); batch fitting delegates to
    `core.batch_fit.run_batch` via the background `FitWorker`/`QThread`
    (`_start_batch`), which parallelizes per-ROI fits across a process pool.
    See `docs/architecture.md` for the layout, the Z-correction logic, and the
    batch/selection flow.
    """

    def __init__(self):
        super().__init__()
        self.setupUi(self)
        self.menubar.setNativeMenuBar(False)

        self._open_files: dict[str, H5Data] = {}
        self._current_file: str | None = None
        self._current_roi: str | None = None
        self._zoom_stack: list[tuple] = []
        self._drag_start = None
        self._batch_thread: QThread | None = None
        self._batch_worker: FitWorker | None = None
        self._pipeline_thread: QThread | None = None
        self._pipeline_worker: PipelineWorker | None = None
        self._theme = _THEMES["light"]

        # Point-selection mode shared by the plot highlight and FitDialog's
        # Default/Manual toggle (see docs/architecture.md "Selection modes").
        # Global/session, not per-ROI: whichever is active applies uniformly
        # across Fit Current and Fit All. "manual" reads each ROI's own
        # persisted Traces/Selection/<roi> mask (set via SelectPointsDialog
        # or "Select points in view"); "default" resolves criteria from the
        # active pipeline (see core.auto_fit_pipeline.resolve_criteria).
        self._selection_mode = "default"   # "default" | "manual"

        # Active fit pipeline (core.auto_fit_pipeline) driving FitDialog's
        # param/criteria seeding, the Run button, and auto-fit-on-open.
        # Remembers the last-selected pipeline across restarts (QSettings);
        # falls back to default_pipeline, then the first available, only
        # the first time the app runs or if the remembered name no longer
        # exists in config.
        pipeline_names = list_pipeline_names(load_pipeline_config())
        remembered = QSettings("TraceEditor", "TraceEditor").value("active_pipeline", "")
        if remembered and remembered in pipeline_names:
            self._active_pipeline: str | None = remembered
        else:
            self._active_pipeline: str | None = (
                DEFAULT_PIPELINE_NAME if DEFAULT_PIPELINE_NAME in pipeline_names
                else (pipeline_names[0] if pipeline_names else None)
            )
        self._auto_fit_on_open = False

        # In-app folder watcher (File > Watch Folder...): polls for new .h5
        # files and opens (+ optionally auto-fits) them as they appear.
        self._watch_folder: str | None = None
        self._watch_known_files: set[str] = set()
        self._watch_pending: dict = {}   # path -> (stable_count, (size, mtime))
        self._watch_timer = QTimer(self)
        self._watch_timer.setInterval(1500)
        self._watch_timer.timeout.connect(self._check_watched_folder)

        self._setup_plot()
        self._connect_signals()
        self._setup_shortcuts()
        self.progressBar.setVisible(False)

        # Fix the width to the checked-state label so toggling doesn't
        # resize the button (both states share the same base text). Measured
        # via sizeHint() (not fontMetrics + a guessed pad) so it reflects
        # whatever the active theme's button padding/weight/case-transform
        # actually renders at, rather than a native-style assumption.
        self.selectedButton.setText("✓ ROI selected")
        width = self.selectedButton.sizeHint().width()
        self.selectedButton.setText("ROI selected")
        self.selectedButton.setFixedWidth(width)

        self.pipelineComboBox.blockSignals(True)
        self.pipelineComboBox.addItems(pipeline_names)
        if self._active_pipeline is not None:
            self.pipelineComboBox.setCurrentText(self._active_pipeline)
        self.pipelineComboBox.blockSignals(False)

        settings = QSettings("TraceEditor", "TraceEditor")
        watch_folder = settings.value("watch_folder", "")
        if watch_folder:
            self._watch_folder = str(watch_folder)
        if settings.value("watch_enabled", False, type=bool) and self._watch_folder:
            self._start_watching()

    # ------------------------------------------------------------------
    # Plot setup
    # ------------------------------------------------------------------

    def _setup_plot(self):
        self.fig = Figure(tight_layout=True)
        self.ax = self.fig.add_subplot(111)
        self.canvas = FigureCanvas(self.fig)

        layout = QVBoxLayout(self.plotWidget)
        layout.setContentsMargins(0, 0, 0, 0)
        layout.addWidget(self.canvas)

        self.canvas.mpl_connect("button_press_event", self._on_canvas_press)
        self.canvas.mpl_connect("button_release_event", self._on_canvas_release)

    # ------------------------------------------------------------------
    # Signals
    # ------------------------------------------------------------------

    def _connect_signals(self):
        self.actionOpen.triggered.connect(self.open_files)
        self.actionImportLegacy.triggered.connect(self._import_legacy_hdf)
        self.actionInspectH5.triggered.connect(self._inspect_h5)
        self.editNotesButton.clicked.connect(self._edit_notes)
        self.actionCloseCurrent.triggered.connect(self._close_current_file)
        self.actionCloseAll.triggered.connect(self._close_all_files)
        self.actionWatchFolder.triggered.connect(self._set_watch_folder)
        self.actionWatchingEnabled.toggled.connect(self._on_watching_toggled)
        self.actionExport.triggered.connect(self._export_data)
        self.actionExit.triggered.connect(self.close)

        self.actionAll.triggered.connect(self.select_all_traces)
        self.actionNone.triggered.connect(self.select_none_traces)
        self.actionToggle.triggered.connect(self._toggle_current_roi_selected)
        self.actionCriteria.triggered.connect(self._select_by_criteria)
        self.actionPointsInRange.triggered.connect(self._select_points_in_range)
        self.actionPointsInView.triggered.connect(self._select_points_in_view)

        self.actionFit.triggered.connect(self.open_fit_dialog)
        self.actionRemoveModelViewed.triggered.connect(self.remove_model_viewed_trace)
        self.actionRemoveModelAll.triggered.connect(self.remove_model_all_traces)
        self.actionRemoveAllProcessing.triggered.connect(self.remove_all_processing)
        self.actionPlotHistogram.triggered.connect(self._plot_histogram)
        self.actionAutoFitOnOpen.toggled.connect(self._on_auto_fit_on_open_toggled)

        self.pipelineComboBox.currentTextChanged.connect(self._on_pipeline_changed)
        self.runPipelineButton.clicked.connect(self._on_run_pipeline_clicked)

        self.actionCalculateForce.triggered.connect(self.calculate_force)
        self.actionFreeDNA.triggered.connect(self._free_dna)
        self.actionRuptures.triggered.connect(self._ruptures)

        self.actionTypeZt.triggered.connect(lambda: self._set_plot_preset("Z(t)"))
        self.actionTypeFZ.triggered.connect(lambda: self._set_plot_preset("F(Z)"))
        self.actionTypeYX.triggered.connect(lambda: self._set_plot_preset("Y(X)"))
        self.actionScaleAuto.triggered.connect(lambda: self._set_scale_preset("Auto"))
        self.actionScaleFix.triggered.connect(lambda: self._set_scale_preset("Fix"))
        self.actionScaleContour.triggered.connect(lambda: self._set_scale_preset("Contour"))
        self.actionCopyImage.triggered.connect(self._copy_image)
        self.actionSaveImage.triggered.connect(self._save_image)
        self.actionBulkSaveImages.triggered.connect(self._bulk_save_images)

        self.fileComboBox.currentIndexChanged.connect(self._on_file_changed)
        self.roiComboBox.currentIndexChanged.connect(self._on_roi_changed)
        self.selectedButton.toggled.connect(self._on_trace_selected_toggled)
        self.selectedOnlyCheckBox.toggled.connect(self._on_selected_only_toggled)

        self.xChannelComboBox.currentIndexChanged.connect(self.redraw)
        self.yChannelComboBox.currentIndexChanged.connect(self.redraw)
        self.xMinEdit.editingFinished.connect(self.redraw)
        self.xMaxEdit.editingFinished.connect(self.redraw)
        self.yMinEdit.editingFinished.connect(self.redraw)
        self.yMaxEdit.editingFinished.connect(self.redraw)
        self.xLogCheckBox.toggled.connect(self.redraw)
        self.yLogCheckBox.toggled.connect(self.redraw)
        self.colorComboBox.currentIndexChanged.connect(self.redraw)
        self.modelComboBox.currentIndexChanged.connect(self.redraw)

        self._color_actions = [
            self.actionColorTheme, self.actionColorBlue, self.actionColorOrange,
            self.actionColorGreen, self.actionColorPink, self.actionColorLightblue,
            self.actionColorPurple, self.actionColorFuchia, self.actionColorYellow,
            self.actionColorTeal, self.actionColorRed,
        ]
        self._color_action_group = QActionGroup(self)
        self._color_action_group.setExclusive(True)
        for i, action in enumerate(self._color_actions):
            self._color_action_group.addAction(action)
            action.triggered.connect(lambda checked, i=i: self.colorComboBox.setCurrentIndex(i))
        self.colorComboBox.currentIndexChanged.connect(self._sync_color_menu)
        self._sync_color_menu(self.colorComboBox.currentIndex())

    def _setup_shortcuts(self):
        QShortcut(QKeySequence("Shift+F12"), self).activated.connect(self._next_roi)
        QShortcut(QKeySequence("Shift+F11"), self).activated.connect(self._prev_roi)
        QShortcut(QKeySequence("Shift+F10"), self).activated.connect(self._toggle_current_roi_selected)
        QShortcut(QKeySequence("Shift+F1"),  self).activated.connect(self._select_by_criteria)
        QShortcut(QKeySequence("Escape"),    self).activated.connect(self._stop_batch)
        QShortcut(QKeySequence("F9"),        self).activated.connect(self._next_file)
        QShortcut(QKeySequence("F8"),        self).activated.connect(self._prev_file)

    # ------------------------------------------------------------------
    # File management
    # ------------------------------------------------------------------

    def open_files(self):
        """Prompt for `.h5` files and load them.

        Opens the `FileBrowserDialog` (a table of file metadata — axis, Notes
        tag, selected/total beads — for the chosen folder), pre-highlighting
        rows already open elsewhere in the app. Registers each newly selected
        file as an `H5Data`, closes any deselected already-open file, adds
        newly opened files to the file combo, and activates the last newly
        added file. LUT files are skipped. Errors per file are shown in the
        status bar, not raised. If auto-fit-on-open is enabled, the active
        pipeline is run once across every newly added file.
        """
        settings = QSettings("TraceEditor", "TraceEditor")
        last_dir = str(settings.value("last_dir", ""))
        dlg = FileBrowserDialog(last_dir, self, open_paths=set(self._open_files.keys()))
        if dlg.exec() != QDialog.Accepted:
            return
        paths = dlg.selected_paths
        if paths:
            settings.setValue("last_dir", str(Path(paths[0]).parent))
        for path in dlg.closed_paths:
            self._close_file(path)
        added = [path for path in paths if self.open_single_file(path)]

        if added:
            self._select_opened_file(added[-1])
            if self._auto_fit_on_open and self._active_pipeline:
                self._run_pipeline_on_files(added)

    def open_single_file(self, path: str) -> bool:
        """Register `path` as an open `H5Data` file and add it to the file combo.

        Shared by `open_files()` (manual, interactive) and the in-app folder
        watcher (`_check_watched_folder`). Does not activate the file — call
        `_activate_opened_file(path)` after, if `path` was newly added.
        Returns whether `path` is now open (`True` if newly added or already
        open; `False` on error or if it's a LUT file).
        """
        if path in self._open_files:
            return True
        try:
            data = H5Data(path)
            if data.is_lut_file():
                self.statusbar.showMessage(f"Skipped LUT file: {Path(path).name}", 3000)
                return False
            self._open_files[path] = data
            self._rebuild_file_combo()
            return True
        except Exception as e:
            self.statusbar.showMessage(f"Error opening {Path(path).name}: {e}", 5000)
            return False

    def _rebuild_file_combo(self):
        """Repopulate `fileComboBox` from `_open_files`, sorted by folder name
        then by the numeric index in the filename, showing "Folder/filename"
        labels. Preserves the current selection."""
        current = self._current_file
        self.fileComboBox.blockSignals(True)
        self.fileComboBox.clear()
        for path in sorted(self._open_files, key=_file_sort_key):
            self.fileComboBox.addItem(_display_label(path), userData=path)
        if current:
            idx = self.fileComboBox.findData(current)
            if idx >= 0:
                self.fileComboBox.setCurrentIndex(idx)
        self.fileComboBox.blockSignals(False)

    def _select_opened_file(self, path: str):
        """Make `path` the current file in the combo box (and activate it)."""
        idx = self.fileComboBox.findData(path)
        if idx >= 0 and self.fileComboBox.currentIndex() != idx:
            self.fileComboBox.setCurrentIndex(idx)
        else:
            # Already at the right index (or not found) — activate directly
            self._activate_file(path)

    def _activate_opened_file(self, path: str):
        self._select_opened_file(path)
        if self._auto_fit_on_open and self._active_pipeline:
            self._run_pipeline_on_files([path])

    def _import_legacy_hdf(self):
        """Import a legacy `.hdf` (PyTables) file, converting it to a new `.h5`
        alongside it and opening the result.
        """
        settings = QSettings("TraceEditor", "TraceEditor")
        last_dir = str(settings.value("last_dir", ""))
        path, _ = QFileDialog.getOpenFileName(
            self, "Import legacy .hdf file", last_dir, "Legacy HDF files (*.hdf);;All files (*)"
        )
        if not path:
            return
        settings.setValue("last_dir", str(Path(path).parent))

        from trace_editor.core.legacy_import import import_legacy_hdf

        try:
            h5_path = import_legacy_hdf(path)
        except FileExistsError:
            self.statusbar.showMessage(
                f"{Path(path).with_suffix('.h5').name} already exists -- delete it or rename to re-import", 5000
            )
            return
        except Exception as e:
            self.statusbar.showMessage(f"Import failed: {e}", 5000)
            return

        h5_path = str(h5_path)
        try:
            data = H5Data(h5_path)
        except Exception as e:
            self.statusbar.showMessage(f"Imported {Path(h5_path).name} but failed to open it: {e}", 5000)
            return

        self._open_files[h5_path] = data
        self._rebuild_file_combo()
        idx = self.fileComboBox.findData(h5_path)
        if idx >= 0 and self.fileComboBox.currentIndex() != idx:
            self.fileComboBox.setCurrentIndex(idx)
        else:
            self._activate_file(h5_path)
        self.statusbar.showMessage(f"Imported {Path(path).name} -> {Path(h5_path).name}", 5000)

    def _inspect_h5(self):
        """Open the current file's raw HDF5 structure in the H5Web VS Code
        extension. Falls back to a dialog pointing at setup instructions and
        the myhdf5.hdfgroup.org viewer if VS Code / the extension aren't found.
        """
        if self._current_file is None:
            self.statusbar.showMessage("Open a file first.", 3000)
            return

        reason_messages = {
            "no_code_cli": "The 'code' command was not found on PATH.",
            "no_extension": "The H5Web VS Code extension is not installed.",
        }
        while True:
            try:
                open_in_vscode_h5web(self._current_file)
                self.statusbar.showMessage(
                    f"Opened {Path(self._current_file).name} in VS Code's H5Web viewer", 3000
                )
                return
            except InspectH5Error as e:
                box = QMessageBox(self)
                box.setWindowTitle("Can't open in VS Code")
                box.setText(
                    f"{reason_messages.get(e.reason, str(e))}\n\n"
                    "See the README for setup instructions, or use the "
                    "myhdf5.hdfgroup.org viewer instead."
                )
                guide_btn = box.addButton("Open Guide", QMessageBox.ActionRole)
                web_btn = box.addButton("Open myhdf5.hdfgroup.org", QMessageBox.ActionRole)
                retry_btn = box.addButton("Try Again", QMessageBox.ActionRole)
                box.addButton(QMessageBox.Cancel)
                box.exec()
                clicked = box.clickedButton()
                if clicked == guide_btn:
                    readme = Path(__file__).resolve().parents[3] / "README.md"
                    QDesktopServices.openUrl(QUrl.fromLocalFile(str(readme)))
                elif clicked == web_btn:
                    QDesktopServices.openUrl(QUrl(MYHDF5_FALLBACK_URL))
                elif clicked == retry_btn:
                    continue
                return

    def _activate_file(self, path: str):
        self._current_file = path
        self._update_channels()
        self._refresh_roi_list()
        self._update_notes_panel()

    def _update_notes_panel(self):
        if not self._current_file:
            self.fileMetaText.setText("")
            self.notesText.setText("")
            return
        data = self._open_files[self._current_file]
        notes = data.get_notes().replace("\\n", "\n").strip()
        self.notesText.setText(notes)

        date = data.get_date() or "—"
        axis = data.get_axis() or "—"
        lut = data.get_lut_filename() or "—"
        tags = extract_tags(data.get_notes())
        tags_text = " ".join(f"#{t}" for t in tags) if tags else "—"
        selected = data.selected_count()
        total = data.num_rois()
        self.fileMetaText.setText(
            f"Date: {date}\nAxis: {axis}\nLUT: {lut}\n"
            f"Selected: {selected}/{total}\nTags: {tags_text}")

    def _edit_notes(self):
        if not self._current_file:
            return
        data = self._open_files[self._current_file]
        from PySide6.QtWidgets import QInputDialog
        text, ok = QInputDialog.getMultiLineText(
            self, f"Edit Notes — {Path(self._current_file).name}", "Notes:", data.get_notes())
        if not ok:
            return
        data.set_notes(text)
        self._update_notes_panel()

    def _on_file_changed(self, index: int):
        path = self.fileComboBox.itemData(index)
        if path and path != self._current_file:
            self._activate_file(path)

    def _close_current_file(self):
        if self._current_file is None:
            return
        self._close_file(self._current_file)

    def _close_file(self, path: str):
        """Close `path`, whether or not it's the currently active file."""
        if path not in self._open_files:
            return
        was_current = path == self._current_file
        self._open_files.pop(path, None)
        self._rebuild_file_combo()
        if not was_current:
            return
        self._current_file = None
        self._current_roi = None
        if self._open_files:
            new_path = self.fileComboBox.itemData(self.fileComboBox.currentIndex())
            if new_path:
                self._activate_file(new_path)
        else:
            self.roiComboBox.clear()
            self.ax.clear()
            self.canvas.draw()

    def _close_all_files(self):
        self._open_files.clear()
        self._current_file = None
        self._current_roi = None
        self.fileComboBox.blockSignals(True)
        self.fileComboBox.clear()
        self.fileComboBox.blockSignals(False)
        self.roiComboBox.clear()
        self.ax.clear()
        self.canvas.draw()

    def _all_file_paths(self) -> list[str]:
        return list(self._open_files.keys())

    # ------------------------------------------------------------------
    # Channel / ROI population
    # ------------------------------------------------------------------

    def _update_channels(self):
        if self._current_file is None:
            return
        data = self._open_files[self._current_file]
        channels = data.all_channels()

        for cb, default in [(self.xChannelComboBox, "Time (s)"),
                            (self.yChannelComboBox, "Z (um)")]:
            cb.blockSignals(True)
            old = cb.currentText()
            cb.clear()
            cb.addItems(channels)
            chosen = old if old in channels else (default if default in channels else channels[0])
            cb.setCurrentText(chosen)
            cb.blockSignals(False)

    def _refresh_roi_list(self):
        if self._current_file is None:
            return
        data = self._open_files[self._current_file]
        all_rois = data.list_rois()
        rois = ([r for r in all_rois if data.get_selection(r) == 1.0]
                if self.selectedOnlyCheckBox.isChecked() else all_rois)

        self.roiComboBox.blockSignals(True)
        current = self.roiComboBox.currentText()
        self.roiComboBox.clear()
        self.roiComboBox.addItems(rois)
        if current in rois:
            self.roiComboBox.setCurrentText(current)
        self.roiComboBox.blockSignals(False)

        self._current_roi = self.roiComboBox.currentText() or (rois[0] if rois else None)
        self._update_trace_selected_checkbox()
        self.redraw()

    def _on_roi_changed(self):
        self._current_roi = self.roiComboBox.currentText() or None
        self._update_trace_selected_checkbox()
        self.redraw()

    def _on_selected_only_toggled(self):
        self._refresh_roi_list()

    def _update_trace_selected_checkbox(self):
        if self._current_file and self._current_roi:
            sel = self._open_files[self._current_file].get_selection(self._current_roi)
            self.selectedButton.blockSignals(True)
            self.selectedButton.setChecked(sel == 1.0)
            self.selectedButton.blockSignals(False)
            self._style_selected_button(sel == 1.0)

    def _on_trace_selected_toggled(self, checked: bool):
        self._style_selected_button(checked)
        if self._current_file and self._current_roi:
            self._open_files[self._current_file].set_selection(
                self._current_roi, 1.0 if checked else 0.0)
            self.redraw()

    def _style_selected_button(self, checked: bool):
        self.selectedButton.setText("✓ ROI selected" if checked else "ROI selected")

    # ------------------------------------------------------------------
    # ROI navigation
    # ------------------------------------------------------------------

    def _next_roi(self):
        # Deliberately does not roll over into the next file (unlike the
        # legacy app's F12 handler) — F9 is the explicit way to change files,
        # so ROI navigation never silently switches which file is open.
        cb = self.roiComboBox
        if cb.currentIndex() < cb.count() - 1:
            cb.setCurrentIndex(cb.currentIndex() + 1)

    def _prev_roi(self):
        cb = self.roiComboBox
        if cb.currentIndex() > 0:
            cb.setCurrentIndex(cb.currentIndex() - 1)

    def _next_file(self):
        # fileComboBox is kept in global (folder, then data_XXX index) order
        # by _rebuild_file_combo, so stepping the index and wrapping around
        # walks files within a folder before moving to the next folder.
        fcb = self.fileComboBox
        if fcb.count() > 0:
            fcb.setCurrentIndex((fcb.currentIndex() + 1) % fcb.count())

    def _prev_file(self):
        fcb = self.fileComboBox
        if fcb.count() > 0:
            fcb.setCurrentIndex((fcb.currentIndex() - 1) % fcb.count())

    def _toggle_current_roi_selected(self):
        if self._current_file and self._current_roi:
            data = self._open_files[self._current_file]
            current = data.get_selection(self._current_roi)
            data.set_selection(self._current_roi, 0.0 if current == 1.0 else 1.0)
            self._update_trace_selected_checkbox()
            self.redraw()

    # ------------------------------------------------------------------
    # Ctrl+drag zoom
    # ------------------------------------------------------------------

    def _on_canvas_press(self, event):
        if event.inaxes != self.ax:
            return
        mods = QApplication.keyboardModifiers()
        if mods & Qt.KeyboardModifier.ControlModifier:
            if event.button == 1:
                self._drag_start = (event.xdata, event.ydata)
            elif event.button == 3:
                if self._zoom_stack:
                    xlim, ylim = self._zoom_stack.pop()
                    self.ax.set_xlim(xlim)
                    self.ax.set_ylim(ylim)
                    self.canvas.draw_idle()

    def _on_canvas_release(self, event):
        if self._drag_start is None:
            return
        mods = QApplication.keyboardModifiers()
        if (mods & Qt.KeyboardModifier.ControlModifier) and event.button == 1 and event.inaxes == self.ax:
            x0, y0 = self._drag_start
            x1, y1 = event.xdata, event.ydata
            if x0 != x1 and y0 != y1:
                self._zoom_stack.append((self.ax.get_xlim(), self.ax.get_ylim()))
                self.ax.set_xlim(sorted([x0, x1]))
                self.ax.set_ylim(sorted([y0, y1]))
                self.canvas.draw_idle()
        self._drag_start = None

    # ------------------------------------------------------------------
    # Z-axis corrections (drift + zbead)
    # ------------------------------------------------------------------

    def _apply_z_corrections(self, axis_name: str, data_array: np.ndarray,
                              time_array: np.ndarray, h5data: H5Data, roi: str) -> np.ndarray:
        if axis_name == "Z (um)" and self.modelComboBox.currentText() != "ExponentialDrift":
            Z0     = h5data.get_parameter("Z0 (um)",      roi, default=0.0)
            dZ_dt  = h5data.get_parameter("dZ_dt (nm_s)", roi, default=0.0)
            tau    = h5data.get_parameter("tau (s)",       roi, default=np.inf)
            z_bead = h5data.get_parameter("Zbead (um)",   roi, default=0.0)
            for v in (Z0, dZ_dt, z_bead):
                pass  # assign below with nan-check
            if math.isnan(Z0):    Z0 = 0.0
            if math.isnan(dZ_dt): dZ_dt = 0.0
            if math.isnan(tau):   tau = np.inf
            if math.isnan(z_bead): z_bead = 0.0
            drift = _exponential_drift(time_array, Z0, dZ_dt, tau)
            return data_array - drift - z_bead
        if axis_name in ("X (um)", "Y (um)"):
            x0 = h5data.get_parameter(axis_name[0] + "0 (um)", roi, default=0.0)
            if math.isnan(x0): x0 = 0.0
            return data_array - x0
        return data_array

    def _apply_axis_limits(self):
        def _parse(text):
            try:
                return float(text)
            except (ValueError, TypeError):
                return None

        xmin = _parse(self.xMinEdit.text())
        xmax = _parse(self.xMaxEdit.text())
        ymin = _parse(self.yMinEdit.text())
        ymax = _parse(self.yMaxEdit.text())

        # Each side is independent: None means autoscale that side
        cur_xmin, cur_xmax = self.ax.get_xlim()
        cur_ymin, cur_ymax = self.ax.get_ylim()
        new_xmin = xmin if xmin is not None else cur_xmin
        new_xmax = xmax if xmax is not None else cur_xmax
        new_ymin = ymin if ymin is not None else cur_ymin
        new_ymax = ymax if ymax is not None else cur_ymax

        if new_xmin != new_xmax:
            self.ax.set_xlim(new_xmin, new_xmax)
        if new_ymin != new_ymax:
            self.ax.set_ylim(new_ymin, new_ymax)

    # ------------------------------------------------------------------
    # Plot / redraw
    # ------------------------------------------------------------------

    def redraw(self, *_):
        """Redraw the plot for the current file/ROI from current control values.

        Scatters the selected X vs Y channels (with drift/Zbead correction
        applied to Z/X/Y), overlays the current model curve if one is chosen, and
        refreshes the parameter list. A no-op with no file/ROI selected. Called
        on every control change and after each fit.
        """
        self.ax.clear()
        t = self._theme
        self.fig.set_facecolor(t["fig_face"])
        self.ax.set_facecolor(t["axes_face"])
        for spine in self.ax.spines.values():
            spine.set_edgecolor(t["axes_edge"])
        self.ax.tick_params(colors=t["axes_edge"], labelcolor=t["label"])

        if not self._current_file or not self._current_roi:
            self.canvas.draw()
            return

        data = self._open_files[self._current_file]
        x_ch = self.xChannelComboBox.currentText()
        y_ch = self.yChannelComboBox.currentText()
        if not x_ch or not y_ch:
            self.canvas.draw()
            return

        try:
            x_raw = data.get_channel(x_ch, self._current_roi)
            y_raw = data.get_channel(y_ch, self._current_roi)
        except Exception as e:
            self.statusbar.showMessage(f"Plot error: {e}", 3000)
            self.canvas.draw()
            return

        try:
            time = data.get_channel("Time (s)", self._current_roi)
        except Exception:
            time = np.arange(len(x_raw), dtype=float)

        x = self._apply_z_corrections(x_ch, x_raw, time, data, self._current_roi)
        y = self._apply_z_corrections(y_ch, y_raw, time, data, self._current_roi)

        model_for_selection = self.modelComboBox.currentText()
        selection_error = None
        try:
            criteria = resolve_criteria(model_for_selection, self._active_pipeline, load_pipeline_config())
            point_mask = _get_fit_selection(
                data, self._current_roi, mode=self._selection_mode, criteria=criteria)
        except Exception as e:
            selection_error = str(e)
            point_mask = None
        # None means "use every point" (see _get_fit_selection) -- highlight
        # all of them, matching what the fit will actually use. This is
        # unrelated to the whole-trace ' selected' flag (include/exclude
        # this ROI from batch ops), which used to be the fallback here but
        # would show nothing highlighted for a deselected ROI even though
        # every point is still what a fit would use.
        mask = np.ones(len(x), dtype=bool) if point_mask is None else point_mask

        scatter_color = _PALETTE_COLORS.get(self.colorComboBox.currentText(), t["scatter_sel"])
        self.ax.scatter(x, y, s=8, color=t["scatter_bg"], alpha=0.5, linewidths=0)
        if mask.any():
            self.ax.scatter(x[mask], y[mask], s=8, color=scatter_color, alpha=0.8, linewidths=0)

        model_name = self.modelComboBox.currentText()
        if model_name and model_name != "none":
            self._draw_model_curve(data, x_ch, y_ch, time)

        self.ax.set_xlabel(x_ch, color=t["label"])
        self.ax.set_ylabel(y_ch, color=t["label"])
        self._apply_axis_limits()
        if self.xLogCheckBox.isChecked():
            self.ax.set_xscale("log")
        if self.yLogCheckBox.isChecked():
            self.ax.set_yscale("log")

        self.canvas.draw()
        self._refresh_params_listbox()
        status = f"{Path(self._current_file).name}  ROI {self._current_roi}"
        if selection_error is not None:
            status += f"  (selection error: {selection_error})"
        self.statusbar.showMessage(status)

    # ------------------------------------------------------------------
    # Image export
    # ------------------------------------------------------------------

    def _copy_image(self):
        """Copy the current plot to the clipboard as a PNG."""
        buf = io.BytesIO()
        self.fig.savefig(buf, format="png", dpi=150, facecolor=self.fig.get_facecolor())
        image = QImage.fromData(buf.getvalue(), "PNG")
        QApplication.clipboard().setImage(image)
        self.statusbar.showMessage("Plot copied to clipboard", 3000)

    def _save_image(self):
        """Prompt for a path and save the current plot to it."""
        settings = QSettings("TraceEditor", "TraceEditor")
        last_dir = str(settings.value("last_save_dir", ""))
        path, _ = QFileDialog.getSaveFileName(
            self, "Save plot image", last_dir,
            "PNG image (*.png);;SVG image (*.svg);;PDF document (*.pdf)")
        if not path:
            return
        settings.setValue("last_save_dir", str(Path(path).parent))
        try:
            self.fig.savefig(path, dpi=150, facecolor=self.fig.get_facecolor())
        except Exception as e:
            self.statusbar.showMessage(f"Save failed: {e}", 5000)
            return
        self.statusbar.showMessage(f"Saved {Path(path).name}", 3000)

    def _bulk_save_images(self):
        """Save one PNG per selected ROI across all open files.

        Redraws each selected ROI into the shared figure (temporarily changing
        the current file/ROI selection) and restores the original selection
        afterwards.
        """
        targets = [(path, roi) for path, data in self._open_files.items()
                   for roi in data.list_rois() if data.get_selection(roi) == 1.0]
        if not targets:
            self.statusbar.showMessage("No selected traces to save", 3000)
            return

        settings = QSettings("TraceEditor", "TraceEditor")
        last_dir = str(settings.value("last_save_dir", ""))
        out_dir = QFileDialog.getExistingDirectory(self, "Bulk save images to folder", last_dir)
        if not out_dir:
            return
        settings.setValue("last_save_dir", out_dir)

        orig_file, orig_roi = self._current_file, self._current_roi
        saved, failed = 0, 0
        for path, roi in targets:
            self._current_file = path
            self._current_roi = roi
            self.redraw()
            out_path = Path(out_dir) / f"{Path(path).stem}_{roi}.png"
            try:
                self.fig.savefig(out_path, dpi=150, facecolor=self.fig.get_facecolor())
                saved += 1
            except Exception:
                failed += 1

        self._current_file, self._current_roi = orig_file, orig_roi
        self.redraw()
        msg = f"Bulk save: {saved} image(s) written to {out_dir}"
        if failed:
            msg += f" ({failed} failed)"
        self.statusbar.showMessage(msg, 5000)

    def _draw_model_curve(self, data: H5Data, x_ch: str, y_ch: str, time: np.ndarray):
        model_name = self.modelComboBox.currentText()
        try:
            model, params, names = prepare_fit(model_name)
        except Exception:
            return

        # Seed params from stored h5 values
        stored = data.list_parameters(self._current_roi)
        for p in params:
            label = names.get(p, p)
            if label in ("x", "y"):
                continue
            val = stored.get(label, float("nan"))
            if not math.isnan(val):
                params[p].value = val
            # z_bead is forced to 0 during display, matching original TraceEditor.py line 507-508.
            # The data is shifted by -Zbead in _apply_z_corrections; the model doesn't need to shift too.
            if p == "z_bead":
                params[p].value = 0.0

        # Evaluate model on its canonical x channel (raw, uncorrected) — matches original line 510-511
        try:
            x_canonical_raw = data.get_channel(names["x"], self._current_roi)
            y_model = model.eval(params=params, x=x_canonical_raw)
        except Exception:
            return

        # The model output (y_model) is already in the same units as the displayed corrected data,
        # because z_bead=0 means the model output starts at 0 just like the zbead-subtracted data.
        # For ExponentialDrift, y_model IS the drift curve in Z(um) and plots directly against time.
        # No z-correction is applied to y_model — it's plotted as-is (matching original).

        # x-axis for plotting: the currently displayed x channel (corrected)
        try:
            x_display_raw = data.get_channel(x_ch, self._current_roi)
            x_display = self._apply_z_corrections(x_ch, x_display_raw, time, data, self._current_roi)
        except Exception:
            return

        if names["y"] == y_ch:
            # Model output matches displayed y-axis → plot (x_display, y_model)
            # Works for both exact-match (ExponentialDrift Z(t)) and
            # indirect-match (zWLC F→Z displayed as Z(t)): arrays share the same index
            self.ax.plot(x_display, y_model, color=self._theme["fit_line"], lw=2)

        elif names["y"] == x_ch:
            # Model output matches the displayed x-axis (e.g. viewing F(Z) with axes swapped)
            try:
                y_display_raw = data.get_channel(y_ch, self._current_roi)
                y_display = self._apply_z_corrections(y_ch, y_display_raw, time, data, self._current_roi)
            except Exception:
                return
            self.ax.plot(y_model, y_display, color=self._theme["fit_line"], lw=2)

    # ------------------------------------------------------------------
    # Params listbox  (Name = Value  [Type / ± error])
    # ------------------------------------------------------------------

    def _refresh_params_listbox(self):
        self.paramsListWidget.blockSignals(True)
        self.paramsListWidget.clear()

        if not self._current_file or not self._current_roi:
            self.paramsListWidget.blockSignals(False)
            return

        data = self._open_files[self._current_file]
        stored = data.list_parameters(self._current_roi)
        model_name = self.modelComboBox.currentText()

        if model_name and model_name != "none":
            try:
                _, params, names = prepare_fit(model_name)
                param_labels = [names[p] for p in params if names[p] not in ("x", "y")]
                param_labels.append(r2_label(model_name))
                param_labels.append(r2_full_label(model_name))
            except Exception:
                param_labels = []
        else:
            param_labels = [name for name, val in stored.items()
                            if name != " selected"
                            and not name.startswith("_src_")
                            and not name.endswith(" Error")
                            and not math.isnan(val)]

        for label in param_labels:
            val = stored.get(label, float("nan"))
            err = stored.get(f"{label} Error", float("nan"))

            if math.isnan(val):
                val_str = "—"
            elif not math.isnan(err) and err > 0:
                val_str = f"{val:.6g} ± {err:.2g}"
            else:
                val_str = f"{val:.6g}"

            item = QListWidgetItem(f"{label} = {val_str}")
            item.setData(Qt.ItemDataRole.UserRole, label)
            self.paramsListWidget.addItem(item)

        self.paramsListWidget.blockSignals(False)

    # ------------------------------------------------------------------
    # Fit dialog
    # ------------------------------------------------------------------

    def open_fit_dialog(self):
        if not self._current_file or not self._current_roi:
            return
        model_name = self.modelComboBox.currentText()
        if not model_name or model_name == "none":
            box = QMessageBox(QMessageBox.Icon.NoIcon, "No model selected",
                               "Select a model before fitting.", QMessageBox.StandardButton.Ok, self)
            box.setOption(QMessageBox.Option.DontUseNativeDialog)
            box.exec()
            return

        data = self._open_files[self._current_file]

        try:
            _, _, names = prepare_fit(model_name)
        except ValueError as e:
            self.statusbar.showMessage(str(e), 5000)
            return

        if names["x"] == "F (pN)":
            calculated = data.get_parameter(
                "F (pN) calculated", self._current_roi, default=0.0)
            if calculated != 1.0:
                box = QMessageBox(
                    QMessageBox.Icon.NoIcon, "Force required",
                    f"'{model_name}' fits against Force, which hasn't been "
                    "calculated for this file yet.\n\nCalculate force for all "
                    "open files now?",
                    QMessageBox.StandardButton.Yes | QMessageBox.StandardButton.Cancel,
                    self,
                )
                box.setOption(QMessageBox.Option.DontUseNativeDialog)
                reply = box.exec()
                if reply != QMessageBox.StandardButton.Yes:
                    return
                self.calculate_force()

        try:
            dlg = FitDialog(
                self,
                model_name=model_name,
                h5data=data,
                roi=self._current_roi,
                live_preview_fn=self.redraw,
            )
        except ValueError as e:
            self.statusbar.showMessage(str(e), 5000)
            return

        dlg.exec()

        if dlg.result_token == "Cancel":
            self.redraw()
            return

        params = dlg.get_params()
        names = dlg.get_names()
        criteria = dlg.get_criteria()

        if dlg.result_token == "OK":
            self._run_fit(data, self._current_roi, model_name, params, names, criteria)
            self.redraw()

        elif dlg.result_token == "Batch":
            tasks = self._batch_roi_tasks()
            if not tasks:
                self.statusbar.showMessage("No selected ROIs for batch fit", 3000)
                return
            self._start_batch(tasks, model_name, params, names, criteria)

    def _run_fit(self, data: H5Data, roi: str, model_name: str, params, names, criteria=None):
        """Fit one model to one ROI and write the results to the h5 file.

        Reads the model's x/y channels, applies drift+Zbead correction to Z via
        `_correct_z_for_fit`, restricts to the model's fit selection (per the
        current Default/Manual toggle, `self._selection_mode`), calls
        `apply_fit`, then stores each parameter (accumulating Zbead),
        `R2 <model_name>`, and `R2 full <model_name>`. Failures are reported
        to the status bar, not raised.

        :param data: the file containing `roi`.
        :param roi: zero-padded ROI id.
        :param model_name: a ForceSpectroscopy model name, e.g. `"zWLC"`.
        :param params: lmfit `Parameters` with vary flags set by the caller.
        :param names: channel/label map from `prepare_fit`.
        :param criteria: selection criteria for "Default" mode, already
            resolved by the caller (`FitDialog.get_criteria`) from the
            active pipeline — see `_get_fit_selection`.
        """
        import copy
        p = copy.deepcopy(params)
        try:
            x = data.get_channel(names["x"], roi)
            y = data.get_channel(names["y"], roi)
            y = _correct_z_for_fit(y, model_name, names, data, roi)
            sel = _get_fit_selection(data, roi, mode=self._selection_mode, criteria=criteria)
            result, out_names = apply_fit(model_name, x, y, selection=sel, initial_parameters=p)
        except Exception as e:
            self.statusbar.showMessage(f"Fit failed: {e}", 5000)
            return

        for pk in result.params:
            _store_fit_result(pk, out_names, result, model_name, data, roi)

        r2 = getattr(result, "r2", np.nan)
        r2_full = getattr(result, "r2_full", np.nan)
        data.set_parameter(r2_label(model_name), r2, roi=roi)
        data.set_parameter(r2_full_label(model_name), r2_full, roi=roi)
        self.statusbar.showMessage(f"Fit done  R²={r2:.4f}  R² (full)={r2_full:.4f}", 5000)

    def _batch_roi_tasks(self) -> list:
        """Return [(filepath, roi)] for every selected ROI across all open files."""
        tasks = []
        for path in self._all_file_paths():
            d = self._open_files[path]
            for roi in d.list_rois():
                if d.get_selection(roi) == 1.0:
                    tasks.append((path, roi))
        return tasks

    def _start_batch(self, tasks, model_name, params, names, criteria=None):
        import copy
        if self._batch_thread is not None and self._batch_thread.isRunning():
            self.statusbar.showMessage("A batch fit is already running", 3000)
            return
        self._batch_thread = QThread(self)
        self._batch_worker = FitWorker(
            tasks, model_name, self._open_files,
            copy.deepcopy(params), names, selection_mode=self._selection_mode,
            criteria=criteria)
        self._batch_worker.moveToThread(self._batch_thread)
        self._batch_thread.started.connect(self._batch_worker.run)
        self._batch_worker.progress.connect(self._on_batch_progress)
        self._batch_worker.finished.connect(self._on_batch_done)
        self._batch_worker.error.connect(self._on_batch_error)
        self._batch_worker.finished.connect(self._batch_thread.quit)
        self._batch_worker.finished.connect(self._batch_worker.deleteLater)
        # Clear our references only once the QThread has actually stopped
        # running (finished fires after the event loop exits), not from the
        # worker's finished handler — nulling self._batch_thread there drops
        # the last Python reference to the QThread while it's still shutting
        # down, which can destroy the C++ object mid-run (crash).
        self._batch_thread.finished.connect(self._clear_batch_thread)
        self.progressBar.setValue(0)
        self.progressBar.setVisible(True)
        self._set_controls_enabled(False)
        self.statusbar.showMessage(f"Fitting {model_name} on {len(tasks)} ROI(s)…")
        self._batch_thread.start()

    def _on_batch_error(self, msg: str):
        self.statusbar.showMessage(f"Fit error: {msg}", 5000)

    def _on_batch_progress(self, done: int, total: int):
        self.progressBar.setValue(100 * done // total)
        self.statusbar.showMessage(f"Fitting… {done}/{total} ROIs")

    def _on_batch_done(self):
        self.progressBar.setVisible(False)
        self._set_controls_enabled(True)
        self.statusbar.showMessage("Batch fit done", 5000)
        self.redraw()

    def _clear_batch_thread(self):
        self._batch_worker = None
        thread, self._batch_thread = self._batch_thread, None
        if thread is not None:
            thread.deleteLater()

    def _stop_batch(self):
        if self._batch_worker:
            self._batch_worker.stop_requested = True

    # ------------------------------------------------------------------
    # Fit pipelines (core.auto_fit_pipeline)
    # ------------------------------------------------------------------

    def _on_pipeline_changed(self, name: str):
        self._active_pipeline = name or None
        QSettings("TraceEditor", "TraceEditor").setValue("active_pipeline", self._active_pipeline)

    def _on_run_pipeline_clicked(self):
        if not self._open_files:
            self.statusbar.showMessage("No files open", 3000)
            return
        self._run_pipeline_on_files(list(self._open_files.keys()))

    def _on_auto_fit_on_open_toggled(self, checked: bool):
        self._auto_fit_on_open = checked

    def _run_pipeline_on_files(self, paths: list):
        """Run the active pipeline against `paths` on a background thread.

        Shared by the Run button (all open files) and auto-fit-on-open (the
        file(s) just opened, whether from a bulk Open or the folder watcher).
        A no-op if no pipeline is selected or `paths` is empty; if a pipeline
        run is already in progress, queues nothing further — the caller can
        retry later (mirrors `_start_batch`'s "already running" guard).
        """
        if not paths or not self._active_pipeline:
            return
        if self._pipeline_thread is not None and self._pipeline_thread.isRunning():
            self.statusbar.showMessage("A pipeline is already running", 3000)
            return

        self._pipeline_thread = QThread(self)
        self._pipeline_worker = PipelineWorker(self._active_pipeline, paths, load_pipeline_config())
        self._pipeline_worker.moveToThread(self._pipeline_thread)
        self._pipeline_thread.started.connect(self._pipeline_worker.run)
        self._pipeline_worker.progress.connect(self._on_pipeline_progress)
        self._pipeline_worker.status.connect(self._on_pipeline_status)
        self._pipeline_worker.finished.connect(self._on_pipeline_done)
        self._pipeline_worker.error.connect(self._on_pipeline_error)
        self._pipeline_worker.finished.connect(self._pipeline_thread.quit)
        self._pipeline_worker.finished.connect(self._pipeline_worker.deleteLater)
        self._pipeline_thread.finished.connect(self._clear_pipeline_thread)
        self.progressBar.setValue(0)
        self.progressBar.setVisible(True)
        self.statusbar.showMessage(
            f"Running pipeline '{self._active_pipeline}' on {len(paths)} file(s)…")
        self._pipeline_thread.start()

    def _on_pipeline_error(self, msg: str):
        self.statusbar.showMessage(f"Pipeline error: {msg}", 5000)

    def _on_pipeline_progress(self, done: int, total: int):
        self.progressBar.setValue(100 * done // total)

    def _on_pipeline_status(self, model_name: str, step: int, total_steps: int,
                             file_index: int, total_files: int):
        """Show which file/step the pipeline is currently on, since the
        progress bar alone resets at every step/file boundary and can't
        convey that on its own."""
        name = Path(self._pipeline_worker.filepaths[file_index - 1]).name if self._pipeline_worker else ""
        self.statusbar.showMessage(
            f"Pipeline: {name} ({file_index}/{total_files})"
            f"  step {step}/{total_steps} ({model_name})")

    def _on_pipeline_done(self):
        self.progressBar.setVisible(False)
        self.statusbar.showMessage("Pipeline done", 5000)
        if self._current_file:
            self._update_channels()
        self.redraw()

    def _clear_pipeline_thread(self):
        self._pipeline_worker = None
        thread, self._pipeline_thread = self._pipeline_thread, None
        if thread is not None:
            thread.deleteLater()

    # ------------------------------------------------------------------
    # Folder watcher (File > Watch Folder...)
    # ------------------------------------------------------------------

    def _set_watch_folder(self):
        settings = QSettings("TraceEditor", "TraceEditor")
        last_dir = self._watch_folder or str(settings.value("last_dir", ""))
        folder = QFileDialog.getExistingDirectory(self, "Watch Folder", last_dir)
        if not folder:
            return
        self._watch_folder = folder
        settings.setValue("watch_folder", folder)
        self.actionWatchingEnabled.setChecked(True)   # triggers _start_watching via toggled

    def _on_watching_toggled(self, checked: bool):
        settings = QSettings("TraceEditor", "TraceEditor")
        settings.setValue("watch_enabled", checked)
        if checked:
            self._start_watching()
        else:
            self._watch_timer.stop()

    def _start_watching(self):
        if not self._watch_folder or not Path(self._watch_folder).is_dir():
            self.statusbar.showMessage("No valid watch folder set", 3000)
            self.actionWatchingEnabled.setChecked(False)
            return
        # Files already present when watching starts are seeded as known,
        # not auto-opened -- only files that appear/finish writing afterward
        # trigger _activate_opened_file.
        self._watch_known_files = {str(p) for p in Path(self._watch_folder).glob("*.h5")}
        self._watch_pending = {}
        self.actionWatchingEnabled.blockSignals(True)
        self.actionWatchingEnabled.setChecked(True)
        self.actionWatchingEnabled.blockSignals(False)
        self._watch_timer.start()

    def _check_watched_folder(self):
        """Poll `self._watch_folder` for new `.h5` files.

        A candidate is opened once its (size, mtime) is unchanged across two
        consecutive polls and it can be opened via `h5py.File(path, "r")`
        (the tweezer software has finished writing and closed it).
        """
        import h5py

        folder = Path(self._watch_folder)
        try:
            candidates = [p for p in folder.glob("*.h5") if str(p) not in self._watch_known_files]
        except OSError as e:
            self.statusbar.showMessage(f"Watch folder error: {e}", 5000)
            self._watch_timer.stop()
            return

        for path in candidates:
            path_str = str(path)
            try:
                stat = path.stat()
            except FileNotFoundError:
                self._watch_pending.pop(path_str, None)
                continue
            fingerprint = (stat.st_size, stat.st_mtime)
            count, prev = self._watch_pending.get(path_str, (0, None))
            count = count + 1 if fingerprint == prev else 1
            self._watch_pending[path_str] = (count, fingerprint)
            if count < 2:
                continue
            try:
                with h5py.File(path, "r"):
                    pass
            except OSError:
                continue

            self._watch_known_files.add(path_str)
            self._watch_pending.pop(path_str, None)
            if self.open_single_file(path_str):
                self._activate_opened_file(path_str)

    def closeEvent(self, event):
        if self._batch_thread is not None and self._batch_thread.isRunning():
            self._stop_batch()
            self._batch_thread.quit()
            self._batch_thread.wait()
        if self._pipeline_thread is not None and self._pipeline_thread.isRunning():
            self._pipeline_thread.quit()
            self._pipeline_thread.wait()
        self._watch_timer.stop()
        super().closeEvent(event)

    def _set_controls_enabled(self, enabled: bool):
        for w in (self.fileComboBox, self.roiComboBox, self.modelComboBox):
            w.setEnabled(enabled)
        self.actionFit.setEnabled(enabled)

    def remove_model_viewed_trace(self):
        """Remove the selected model's fit parameters from the currently viewed trace only."""
        if not self._current_file or not self._current_roi:
            return
        model_name = self.modelComboBox.currentText()
        if not model_name or model_name == "none":
            return
        try:
            _, _, names = prepare_fit(model_name)
        except Exception:
            return
        data = self._open_files[self._current_file]
        for p, label in names.items():
            if label not in ("x", "y"):
                data.delete_parameter(label, self._current_roi)
                data.delete_parameter(f"{label} Error", self._current_roi)
        self.redraw()

    def remove_model_all_traces(self):
        """Remove the selected model's fit parameters from every ROI in every open file."""
        model_name = self.modelComboBox.currentText()
        if not model_name or model_name == "none":
            return
        try:
            _, _, names = prepare_fit(model_name)
        except Exception:
            return
        labels = [label for p, label in names.items() if label not in ("x", "y")]
        for path in self._all_file_paths():
            data = self._open_files[path]
            for roi in data.list_rois():
                for label in labels:
                    data.delete_parameter(label, roi)
                    data.delete_parameter(f"{label} Error", roi)
        self.redraw()

    def remove_all_processing(self):
        """Remove all calculated parameters, F (pN), and manual point
        selection masks from all open files -- everything derived from the
        raw instrument data, which is left untouched."""
        import h5py
        instrument_fields = {"ROI", "X (pix)", "Y (pix)", " selected"}
        paths = self._all_file_paths()
        self.progressBar.setValue(0)
        self.progressBar.setVisible(True)
        self._set_controls_enabled(False)
        self.statusbar.showMessage(f"Resetting {len(paths)} file(s)…")
        QApplication.processEvents()
        for i, path in enumerate(paths):
            data = self._open_files[path]
            self.statusbar.showMessage(
                f"Resetting {Path(path).name} ({i + 1}/{len(paths)})…")
            for roi in data.list_rois():
                stored = data.list_parameters(roi)
                for name in list(stored.keys()):
                    if name not in instrument_fields:
                        data.delete_parameter(name, roi)
            with h5py.File(data.filepath, "r+") as f:
                grp = f["Traces"] if "Traces" in f else f
                if "F (pN)" in grp:
                    del grp["F (pN)"]
                if "Selection" in grp:
                    del grp["Selection"]
            self.progressBar.setValue(100 * (i + 1) // len(paths))
            QApplication.processEvents()
        self.progressBar.setVisible(False)
        self._set_controls_enabled(True)
        if self._current_file:
            self._update_channels()
        self.redraw()
        self.statusbar.showMessage("Data reset", 5000)

    # ------------------------------------------------------------------
    # Force calculation
    # ------------------------------------------------------------------

    def calculate_force(self):
        """Convert Shift (mm) → F (pN) for every open file via `MagnetForce`.

        Uses the active pipeline's `calculate_force` step params (falling back
        to `[force_defaults]`, via `resolve_force_defaults` — the same
        resolution `run_pipeline` uses for a `calculate_force` step), writes
        the result as the shared `F (pN)` dataset, and marks every ROI with
        `F (pN) calculated = 1.0`. Runs across all open files with a progress
        bar.
        """
        config = load_pipeline_config()
        force_defaults = resolve_force_defaults(self._active_pipeline, config)
        paths = self._all_file_paths()
        self.progressBar.setValue(0)
        self.progressBar.setVisible(True)
        self._set_controls_enabled(False)
        self.statusbar.showMessage(f"Calculating force for {len(paths)} file(s)…")
        QApplication.processEvents()
        count = 0
        for i, path in enumerate(paths):
            data = self._open_files[path]
            self.statusbar.showMessage(
                f"Calculating force: {Path(path).name} ({i + 1}/{len(paths)})…")
            try:
                calculate_force_for_file(data, force_defaults)
            except Exception as e:
                self.statusbar.showMessage(f"No Shift (mm) in {Path(path).name}: {e}", 5000)
                self.progressBar.setValue(100 * (i + 1) // len(paths))
                QApplication.processEvents()
                continue
            count += len(data.list_rois())

            self.progressBar.setValue(100 * (i + 1) // len(paths))
            QApplication.processEvents()

        self.progressBar.setVisible(False)
        self._set_controls_enabled(True)
        self.statusbar.showMessage(f"Force calculated for {count} ROIs", 5000)
        if self._current_file:
            self._update_channels()
        self.redraw()

    # ------------------------------------------------------------------
    # Trace selection
    # ------------------------------------------------------------------

    def select_all_traces(self):
        for path in self._all_file_paths():
            data = self._open_files[path]
            for roi in data.list_rois():
                data.set_selection(roi, 1.0)
        self._update_trace_selected_checkbox()
        self.redraw()

    def select_none_traces(self):
        for path in self._all_file_paths():
            data = self._open_files[path]
            for roi in data.list_rois():
                data.set_selection(roi, 0.0)
        self._update_trace_selected_checkbox()
        self.redraw()

    def _select_points_in_range(self):
        if not self._current_file or not self._current_roi:
            return
        data = self._open_files[self._current_file]
        roi  = self._current_roi

        dlg = SelectPointsDialog(data, roi, live_fn=self.redraw, parent=self)
        dlg.exec()

        if dlg.result_token != "Cancel":
            # The dialog's live preview already persisted a mask for this ROI
            # (see SelectPointsDialog._on_changed) -- switch to Manual mode so
            # the Fit dialog picks it up.
            self._selection_mode = "manual"

        if dlg.result_token == "Batch" and dlg.criteria:
            tasks = self._batch_roi_tasks()
            for fpath, r in tasks:
                _apply_point_selection(dlg.criteria, self._open_files[fpath], r)
            self.redraw()
            self.statusbar.showMessage(
                f"Point selection applied to {len(tasks)} ROIs", 4000)
        else:
            self.redraw()

    def _select_points_in_view(self):
        """Select points in view: use the current axis limits as range criteria.

        Legacy shortcut (`select_points_in_view` in `oldhdf_versions/TraceEditor.py`)
        reimagined as an alternate way to build criteria for `_apply_point_selection` —
        the channel and min/max come from whatever's currently on the x/y axis and the
        plot's current view limits, instead of being typed into `SelectPointsDialog`.
        Persists the resulting mask the same way `SelectPointsDialog` does, and
        switches to Manual selection mode so the Fit dialog picks it up.
        """
        if not self._current_file or not self._current_roi:
            return
        data = self._open_files[self._current_file]
        roi = self._current_roi
        x_ch = self.xChannelComboBox.currentText()
        y_ch = self.yChannelComboBox.currentText()
        if not x_ch or not y_ch:
            return

        xmin, xmax = self.ax.get_xlim()
        ymin, ymax = self.ax.get_ylim()
        criteria = [
            {"channel": x_ch, "min": xmin, "max": xmax, "derivative": False, "invert": False},
            {"channel": y_ch, "min": ymin, "max": ymax, "derivative": False, "invert": False},
        ]
        _apply_point_selection(criteria, data, roi)
        self._selection_mode = "manual"
        self.redraw()
        self.statusbar.showMessage("Selected points in view", 3000)

    def _select_by_criteria(self):
        if not self._current_file or not self._current_roi:
            return
        data = self._open_files[self._current_file]

        dlg = SelectTracesDialog(data, self._current_roi, parent=self)
        if dlg.exec() != QDialog.Accepted or not dlg.criteria:
            return

        tasks = self._batch_roi_tasks()
        n_selected = 0
        for fpath, roi in tasks:
            if _apply_trace_selection(dlg.criteria, self._open_files[fpath], roi):
                n_selected += 1
        self._refresh_roi_list()
        self.redraw()
        self.statusbar.showMessage(
            f"{n_selected}/{len(tasks)} ROIs selected by criteria", 4000)

    # ------------------------------------------------------------------
    # Plot presets
    # ------------------------------------------------------------------

    def _sync_color_menu(self, index: int):
        if 0 <= index < len(self._color_actions):
            self._color_actions[index].setChecked(True)

    def _set_plot_preset(self, preset: str):
        presets = {
            "Z(t)": ("Time (s)", "Z (um)"),
            "F(Z)": ("Z (um)",   "F (pN)"),
            "Y(X)": ("X (um)",   "Y (um)"),
        }
        x_ch, y_ch = presets.get(preset, (None, None))
        if not x_ch:
            return
        channels = (self._open_files[self._current_file].all_channels()
                    if self._current_file else [])
        if x_ch in channels:
            self.xChannelComboBox.setCurrentText(x_ch)
        if y_ch in channels:
            self.yChannelComboBox.setCurrentText(y_ch)
        for w in (self.xMinEdit, self.xMaxEdit, self.yMinEdit, self.yMaxEdit):
            w.clear()

    def _set_scale_preset(self, preset: str):
        if preset == "Auto":
            for w in (self.xMinEdit, self.xMaxEdit, self.yMinEdit, self.yMaxEdit):
                w.clear()
            self.redraw()
        elif preset == "Fix":
            self.xMinEdit.setText(f"{self.ax.get_xlim()[0]:.4g}")
            self.xMaxEdit.setText(f"{self.ax.get_xlim()[1]:.4g}")
            self.yMinEdit.setText(f"{self.ax.get_ylim()[0]:.4g}")
            self.yMaxEdit.setText(f"{self.ax.get_ylim()[1]:.4g}")
        elif preset == "Contour":
            if self._current_file and self._current_roi:
                data = self._open_files[self._current_file]
                L = data.get_parameter("L (bp)", self._current_roi, default=float("nan"))
                if not math.isnan(L) and L > 1:
                    z_max = L * 0.34 / 1000
                    self.yMinEdit.setText(f"{-0.1 * z_max:.4g}")
                    self.yMaxEdit.setText(f"{1.1 * z_max:.4g}")
                    self.redraw()

    # ------------------------------------------------------------------
    # Stubs
    # ------------------------------------------------------------------

    def _export_data(self):
        self.statusbar.showMessage("Export: not yet implemented", 3000)

    def _plot_histogram(self):
        self.statusbar.showMessage("Parameter histogram: not yet implemented", 3000)

    def _free_dna(self):
        self.statusbar.showMessage("Free DNA: not yet implemented", 3000)

    def _ruptures(self):
        self.statusbar.showMessage("Ruptures: not yet implemented", 3000)


def main():
    import qt_material

    app = QApplication(sys.argv)
    # resources/theme.xml's colors are matched to resources/icon.svg's
    # palette. A real theme file, not qt_material's `extra` param: `extra`
    # only patches the generated QSS, but qt_material renders the actual
    # checkbox/radio/arrow icon assets straight from the theme file *before*
    # `extra` is merged in, so icons would still use the stock theme's
    # colors (this is why the checkbox tick looked like leftover neon teal).
    # primaryColor and primaryLightColor are both the icon's darkest teal:
    # primaryColor fills those solid icons, where the icon's own bright mint
    # accent reads as too vivid, and primaryLightColor is the menu/combobox
    # hover background drawn under white text, which needs the same darkness.
    theme_file = Path(__file__).parent / "resources" / "theme.xml"
    qt_material.apply_stylesheet(app, theme=str(theme_file))
    app.setWindowIcon(QIcon(str(Path(__file__).parent / "resources" / "icon.png")))
    window = TraceEditor()
    window.show()
    sys.exit(app.exec())


if __name__ == "__main__":
    main()
