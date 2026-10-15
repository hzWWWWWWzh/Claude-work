# -*- coding: utf-8 -*-

################################################################################
## Form generated from reading UI file 'main_window.ui'
##
## Created by: Qt User Interface Compiler version 6.11.1
##
## WARNING! All changes made in this file will be lost when recompiling UI file!
################################################################################

from PySide6.QtCore import (QCoreApplication, QDate, QDateTime, QLocale,
    QMetaObject, QObject, QPoint, QRect,
    QSize, QTime, QUrl, Qt)
from PySide6.QtGui import (QAction, QBrush, QColor, QConicalGradient,
    QCursor, QFont, QFontDatabase, QGradient,
    QIcon, QImage, QKeySequence, QLinearGradient,
    QPainter, QPalette, QPixmap, QRadialGradient,
    QTransform)
from PySide6.QtWidgets import (QApplication, QCheckBox, QComboBox, QFrame,
    QHBoxLayout, QLabel, QLineEdit, QListWidget,
    QListWidgetItem, QMainWindow, QMenu, QMenuBar,
    QProgressBar, QPushButton, QScrollArea, QSizePolicy,
    QSpacerItem, QStatusBar, QVBoxLayout, QWidget)

class Ui_MainWindow(object):
    def setupUi(self, MainWindow):
        if not MainWindow.objectName():
            MainWindow.setObjectName(u"MainWindow")
        MainWindow.resize(1280, 900)
        self.actionOpen = QAction(MainWindow)
        self.actionOpen.setObjectName(u"actionOpen")
        self.actionImportLegacy = QAction(MainWindow)
        self.actionImportLegacy.setObjectName(u"actionImportLegacy")
        self.actionInspectH5 = QAction(MainWindow)
        self.actionInspectH5.setObjectName(u"actionInspectH5")
        self.actionCloseCurrent = QAction(MainWindow)
        self.actionCloseCurrent.setObjectName(u"actionCloseCurrent")
        self.actionCloseAll = QAction(MainWindow)
        self.actionCloseAll.setObjectName(u"actionCloseAll")
        self.actionExport = QAction(MainWindow)
        self.actionExport.setObjectName(u"actionExport")
        self.actionExit = QAction(MainWindow)
        self.actionExit.setObjectName(u"actionExit")
        self.actionWatchFolder = QAction(MainWindow)
        self.actionWatchFolder.setObjectName(u"actionWatchFolder")
        self.actionWatchingEnabled = QAction(MainWindow)
        self.actionWatchingEnabled.setObjectName(u"actionWatchingEnabled")
        self.actionWatchingEnabled.setCheckable(True)
        self.actionPointsInRange = QAction(MainWindow)
        self.actionPointsInRange.setObjectName(u"actionPointsInRange")
        self.actionPointsInView = QAction(MainWindow)
        self.actionPointsInView.setObjectName(u"actionPointsInView")
        self.actionAll = QAction(MainWindow)
        self.actionAll.setObjectName(u"actionAll")
        self.actionNone = QAction(MainWindow)
        self.actionNone.setObjectName(u"actionNone")
        self.actionToggle = QAction(MainWindow)
        self.actionToggle.setObjectName(u"actionToggle")
        self.actionCriteria = QAction(MainWindow)
        self.actionCriteria.setObjectName(u"actionCriteria")
        self.actionFit = QAction(MainWindow)
        self.actionFit.setObjectName(u"actionFit")
        self.actionPlotHistogram = QAction(MainWindow)
        self.actionPlotHistogram.setObjectName(u"actionPlotHistogram")
        self.actionAutoFitOnOpen = QAction(MainWindow)
        self.actionAutoFitOnOpen.setObjectName(u"actionAutoFitOnOpen")
        self.actionAutoFitOnOpen.setCheckable(True)
        self.actionRemoveModelViewed = QAction(MainWindow)
        self.actionRemoveModelViewed.setObjectName(u"actionRemoveModelViewed")
        self.actionRemoveModelAll = QAction(MainWindow)
        self.actionRemoveModelAll.setObjectName(u"actionRemoveModelAll")
        self.actionRemoveAllProcessing = QAction(MainWindow)
        self.actionRemoveAllProcessing.setObjectName(u"actionRemoveAllProcessing")
        self.actionCalculateForce = QAction(MainWindow)
        self.actionCalculateForce.setObjectName(u"actionCalculateForce")
        self.actionFreeDNA = QAction(MainWindow)
        self.actionFreeDNA.setObjectName(u"actionFreeDNA")
        self.actionRuptures = QAction(MainWindow)
        self.actionRuptures.setObjectName(u"actionRuptures")
        self.actionTypeZt = QAction(MainWindow)
        self.actionTypeZt.setObjectName(u"actionTypeZt")
        self.actionTypeFZ = QAction(MainWindow)
        self.actionTypeFZ.setObjectName(u"actionTypeFZ")
        self.actionTypeYX = QAction(MainWindow)
        self.actionTypeYX.setObjectName(u"actionTypeYX")
        self.actionScaleAuto = QAction(MainWindow)
        self.actionScaleAuto.setObjectName(u"actionScaleAuto")
        self.actionScaleFix = QAction(MainWindow)
        self.actionScaleFix.setObjectName(u"actionScaleFix")
        self.actionColorTheme = QAction(MainWindow)
        self.actionColorTheme.setObjectName(u"actionColorTheme")
        self.actionColorTheme.setCheckable(True)
        self.actionColorBlue = QAction(MainWindow)
        self.actionColorBlue.setObjectName(u"actionColorBlue")
        self.actionColorBlue.setCheckable(True)
        self.actionColorOrange = QAction(MainWindow)
        self.actionColorOrange.setObjectName(u"actionColorOrange")
        self.actionColorOrange.setCheckable(True)
        self.actionColorGreen = QAction(MainWindow)
        self.actionColorGreen.setObjectName(u"actionColorGreen")
        self.actionColorGreen.setCheckable(True)
        self.actionColorPink = QAction(MainWindow)
        self.actionColorPink.setObjectName(u"actionColorPink")
        self.actionColorPink.setCheckable(True)
        self.actionColorLightblue = QAction(MainWindow)
        self.actionColorLightblue.setObjectName(u"actionColorLightblue")
        self.actionColorLightblue.setCheckable(True)
        self.actionColorPurple = QAction(MainWindow)
        self.actionColorPurple.setObjectName(u"actionColorPurple")
        self.actionColorPurple.setCheckable(True)
        self.actionColorFuchia = QAction(MainWindow)
        self.actionColorFuchia.setObjectName(u"actionColorFuchia")
        self.actionColorFuchia.setCheckable(True)
        self.actionColorYellow = QAction(MainWindow)
        self.actionColorYellow.setObjectName(u"actionColorYellow")
        self.actionColorYellow.setCheckable(True)
        self.actionColorTeal = QAction(MainWindow)
        self.actionColorTeal.setObjectName(u"actionColorTeal")
        self.actionColorTeal.setCheckable(True)
        self.actionColorRed = QAction(MainWindow)
        self.actionColorRed.setObjectName(u"actionColorRed")
        self.actionColorRed.setCheckable(True)
        self.actionScaleContour = QAction(MainWindow)
        self.actionScaleContour.setObjectName(u"actionScaleContour")
        self.actionCopyImage = QAction(MainWindow)
        self.actionCopyImage.setObjectName(u"actionCopyImage")
        self.actionSaveImage = QAction(MainWindow)
        self.actionSaveImage.setObjectName(u"actionSaveImage")
        self.actionBulkSaveImages = QAction(MainWindow)
        self.actionBulkSaveImages.setObjectName(u"actionBulkSaveImages")
        self.centralwidget = QWidget(MainWindow)
        self.centralwidget.setObjectName(u"centralwidget")
        self.mainLayout = QHBoxLayout(self.centralwidget)
        self.mainLayout.setSpacing(4)
        self.mainLayout.setObjectName(u"mainLayout")
        self.mainLayout.setContentsMargins(4, 4, 4, 4)
        self.leftWidget = QWidget(self.centralwidget)
        self.leftWidget.setObjectName(u"leftWidget")
        self.leftWidget.setMinimumWidth(200)
        self.leftWidget.setMaximumWidth(260)
        self.leftLayout = QVBoxLayout(self.leftWidget)
        self.leftLayout.setSpacing(4)
        self.leftLayout.setObjectName(u"leftLayout")
        self.leftLayout.setContentsMargins(0, 0, 0, 0)
        self.modelLabel = QLabel(self.leftWidget)
        self.modelLabel.setObjectName(u"modelLabel")

        self.leftLayout.addWidget(self.modelLabel)

        self.modelComboBox = QComboBox(self.leftWidget)
        self.modelComboBox.addItem("")
        self.modelComboBox.addItem("")
        self.modelComboBox.addItem("")
        self.modelComboBox.addItem("")
        self.modelComboBox.addItem("")
        self.modelComboBox.addItem("")
        self.modelComboBox.addItem("")
        self.modelComboBox.setObjectName(u"modelComboBox")

        self.leftLayout.addWidget(self.modelComboBox)

        self.parametersLabel = QLabel(self.leftWidget)
        self.parametersLabel.setObjectName(u"parametersLabel")

        self.leftLayout.addWidget(self.parametersLabel)

        self.paramsListWidget = QListWidget(self.leftWidget)
        self.paramsListWidget.setObjectName(u"paramsListWidget")
        sizePolicy = QSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Expanding)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(1)
        sizePolicy.setHeightForWidth(self.paramsListWidget.sizePolicy().hasHeightForWidth())
        self.paramsListWidget.setSizePolicy(sizePolicy)

        self.leftLayout.addWidget(self.paramsListWidget)

        self.pipelineLabel = QLabel(self.leftWidget)
        self.pipelineLabel.setObjectName(u"pipelineLabel")

        self.leftLayout.addWidget(self.pipelineLabel)

        self.pipelineLayout = QHBoxLayout()
        self.pipelineLayout.setObjectName(u"pipelineLayout")
        self.pipelineComboBox = QComboBox(self.leftWidget)
        self.pipelineComboBox.setObjectName(u"pipelineComboBox")

        self.pipelineLayout.addWidget(self.pipelineComboBox)

        self.runPipelineButton = QPushButton(self.leftWidget)
        self.runPipelineButton.setObjectName(u"runPipelineButton")

        self.pipelineLayout.addWidget(self.runPipelineButton)


        self.leftLayout.addLayout(self.pipelineLayout)


        self.mainLayout.addWidget(self.leftWidget)

        self.centerWidget = QWidget(self.centralwidget)
        self.centerWidget.setObjectName(u"centerWidget")
        sizePolicy1 = QSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Preferred)
        sizePolicy1.setHorizontalStretch(3)
        sizePolicy1.setVerticalStretch(0)
        sizePolicy1.setHeightForWidth(self.centerWidget.sizePolicy().hasHeightForWidth())
        self.centerWidget.setSizePolicy(sizePolicy1)
        self.centerLayout = QVBoxLayout(self.centerWidget)
        self.centerLayout.setSpacing(4)
        self.centerLayout.setObjectName(u"centerLayout")
        self.centerLayout.setContentsMargins(0, 0, 0, 0)
        self.fileSelectorLayout = QHBoxLayout()
        self.fileSelectorLayout.setObjectName(u"fileSelectorLayout")
        self.filenameLabel = QLabel(self.centerWidget)
        self.filenameLabel.setObjectName(u"filenameLabel")

        self.fileSelectorLayout.addWidget(self.filenameLabel)

        self.fileComboBox = QComboBox(self.centerWidget)
        self.fileComboBox.setObjectName(u"fileComboBox")
        sizePolicy2 = QSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Fixed)
        sizePolicy2.setHorizontalStretch(0)
        sizePolicy2.setVerticalStretch(0)
        sizePolicy2.setHeightForWidth(self.fileComboBox.sizePolicy().hasHeightForWidth())
        self.fileComboBox.setSizePolicy(sizePolicy2)

        self.fileSelectorLayout.addWidget(self.fileComboBox)

        self.roiLabel = QLabel(self.centerWidget)
        self.roiLabel.setObjectName(u"roiLabel")

        self.fileSelectorLayout.addWidget(self.roiLabel)

        self.roiComboBox = QComboBox(self.centerWidget)
        self.roiComboBox.setObjectName(u"roiComboBox")

        self.fileSelectorLayout.addWidget(self.roiComboBox)

        self.selectedButton = QPushButton(self.centerWidget)
        self.selectedButton.setObjectName(u"selectedButton")
        self.selectedButton.setCheckable(True)
        sizePolicy3 = QSizePolicy(QSizePolicy.Policy.Fixed, QSizePolicy.Policy.Fixed)
        sizePolicy3.setHorizontalStretch(0)
        sizePolicy3.setVerticalStretch(0)
        sizePolicy3.setHeightForWidth(self.selectedButton.sizePolicy().hasHeightForWidth())
        self.selectedButton.setSizePolicy(sizePolicy3)

        self.fileSelectorLayout.addWidget(self.selectedButton)


        self.centerLayout.addLayout(self.fileSelectorLayout)

        self.plotWidget = QWidget(self.centerWidget)
        self.plotWidget.setObjectName(u"plotWidget")
        sizePolicy.setHeightForWidth(self.plotWidget.sizePolicy().hasHeightForWidth())
        self.plotWidget.setSizePolicy(sizePolicy)

        self.centerLayout.addWidget(self.plotWidget)

        self.progressBar = QProgressBar(self.centerWidget)
        self.progressBar.setObjectName(u"progressBar")
        self.progressBar.setMaximum(100)
        self.progressBar.setValue(0)
        self.progressBar.setVisible(False)

        self.centerLayout.addWidget(self.progressBar)


        self.mainLayout.addWidget(self.centerWidget)

        self.rightWidget = QWidget(self.centralwidget)
        self.rightWidget.setObjectName(u"rightWidget")
        self.rightWidget.setMinimumWidth(180)
        self.rightWidget.setMaximumWidth(240)
        self.rightLayout = QVBoxLayout(self.rightWidget)
        self.rightLayout.setSpacing(8)
        self.rightLayout.setObjectName(u"rightLayout")
        self.rightLayout.setContentsMargins(4, 4, 4, 4)
        self.selectedOnlyCheckBox = QPushButton(self.rightWidget)
        self.selectedOnlyCheckBox.setObjectName(u"selectedOnlyCheckBox")
        self.selectedOnlyCheckBox.setCheckable(True)

        self.rightLayout.addWidget(self.selectedOnlyCheckBox)

        self.sp1 = QSpacerItem(0, 6, QSizePolicy.Policy.Minimum, QSizePolicy.Policy.Fixed)

        self.rightLayout.addItem(self.sp1)

        self.hline1 = QFrame(self.rightWidget)
        self.hline1.setObjectName(u"hline1")
        self.hline1.setFrameShape(QFrame.Shape.HLine)
        self.hline1.setFrameShadow(QFrame.Shadow.Sunken)

        self.rightLayout.addWidget(self.hline1)

        self.sp2 = QSpacerItem(0, 6, QSizePolicy.Policy.Minimum, QSizePolicy.Policy.Fixed)

        self.rightLayout.addItem(self.sp2)

        self.xAxisHeaderLayout = QHBoxLayout()
        self.xAxisHeaderLayout.setObjectName(u"xAxisHeaderLayout")
        self.xAxisLabel = QLabel(self.rightWidget)
        self.xAxisLabel.setObjectName(u"xAxisLabel")

        self.xAxisHeaderLayout.addWidget(self.xAxisLabel)

        self.xChannelComboBox = QComboBox(self.rightWidget)
        self.xChannelComboBox.setObjectName(u"xChannelComboBox")
        sizePolicy4 = QSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Fixed)
        sizePolicy4.setHorizontalStretch(1)
        sizePolicy4.setVerticalStretch(0)
        sizePolicy4.setHeightForWidth(self.xChannelComboBox.sizePolicy().hasHeightForWidth())
        self.xChannelComboBox.setSizePolicy(sizePolicy4)

        self.xAxisHeaderLayout.addWidget(self.xChannelComboBox)


        self.rightLayout.addLayout(self.xAxisHeaderLayout)

        self.xMinMaxLayout = QHBoxLayout()
        self.xMinMaxLayout.setObjectName(u"xMinMaxLayout")
        self.xMinLabel = QLabel(self.rightWidget)
        self.xMinLabel.setObjectName(u"xMinLabel")

        self.xMinMaxLayout.addWidget(self.xMinLabel)

        self.xMinEdit = QLineEdit(self.rightWidget)
        self.xMinEdit.setObjectName(u"xMinEdit")
        sizePolicy4.setHeightForWidth(self.xMinEdit.sizePolicy().hasHeightForWidth())
        self.xMinEdit.setSizePolicy(sizePolicy4)

        self.xMinMaxLayout.addWidget(self.xMinEdit)

        self.xMaxLabel = QLabel(self.rightWidget)
        self.xMaxLabel.setObjectName(u"xMaxLabel")

        self.xMinMaxLayout.addWidget(self.xMaxLabel)

        self.xMaxEdit = QLineEdit(self.rightWidget)
        self.xMaxEdit.setObjectName(u"xMaxEdit")
        sizePolicy4.setHeightForWidth(self.xMaxEdit.sizePolicy().hasHeightForWidth())
        self.xMaxEdit.setSizePolicy(sizePolicy4)

        self.xMinMaxLayout.addWidget(self.xMaxEdit)

        self.xLogCheckBox = QCheckBox(self.rightWidget)
        self.xLogCheckBox.setObjectName(u"xLogCheckBox")

        self.xMinMaxLayout.addWidget(self.xLogCheckBox)


        self.rightLayout.addLayout(self.xMinMaxLayout)

        self.sp3 = QSpacerItem(0, 6, QSizePolicy.Policy.Minimum, QSizePolicy.Policy.Fixed)

        self.rightLayout.addItem(self.sp3)

        self.hline2 = QFrame(self.rightWidget)
        self.hline2.setObjectName(u"hline2")
        self.hline2.setFrameShape(QFrame.Shape.HLine)
        self.hline2.setFrameShadow(QFrame.Shadow.Sunken)

        self.rightLayout.addWidget(self.hline2)

        self.sp4 = QSpacerItem(0, 6, QSizePolicy.Policy.Minimum, QSizePolicy.Policy.Fixed)

        self.rightLayout.addItem(self.sp4)

        self.yAxisHeaderLayout = QHBoxLayout()
        self.yAxisHeaderLayout.setObjectName(u"yAxisHeaderLayout")
        self.yAxisLabel = QLabel(self.rightWidget)
        self.yAxisLabel.setObjectName(u"yAxisLabel")

        self.yAxisHeaderLayout.addWidget(self.yAxisLabel)

        self.yChannelComboBox = QComboBox(self.rightWidget)
        self.yChannelComboBox.setObjectName(u"yChannelComboBox")
        sizePolicy4.setHeightForWidth(self.yChannelComboBox.sizePolicy().hasHeightForWidth())
        self.yChannelComboBox.setSizePolicy(sizePolicy4)

        self.yAxisHeaderLayout.addWidget(self.yChannelComboBox)


        self.rightLayout.addLayout(self.yAxisHeaderLayout)

        self.yMinMaxLayout = QHBoxLayout()
        self.yMinMaxLayout.setObjectName(u"yMinMaxLayout")
        self.yMinLabel = QLabel(self.rightWidget)
        self.yMinLabel.setObjectName(u"yMinLabel")

        self.yMinMaxLayout.addWidget(self.yMinLabel)

        self.yMinEdit = QLineEdit(self.rightWidget)
        self.yMinEdit.setObjectName(u"yMinEdit")
        sizePolicy4.setHeightForWidth(self.yMinEdit.sizePolicy().hasHeightForWidth())
        self.yMinEdit.setSizePolicy(sizePolicy4)

        self.yMinMaxLayout.addWidget(self.yMinEdit)

        self.yMaxLabel = QLabel(self.rightWidget)
        self.yMaxLabel.setObjectName(u"yMaxLabel")

        self.yMinMaxLayout.addWidget(self.yMaxLabel)

        self.yMaxEdit = QLineEdit(self.rightWidget)
        self.yMaxEdit.setObjectName(u"yMaxEdit")
        sizePolicy4.setHeightForWidth(self.yMaxEdit.sizePolicy().hasHeightForWidth())
        self.yMaxEdit.setSizePolicy(sizePolicy4)

        self.yMinMaxLayout.addWidget(self.yMaxEdit)

        self.yLogCheckBox = QCheckBox(self.rightWidget)
        self.yLogCheckBox.setObjectName(u"yLogCheckBox")

        self.yMinMaxLayout.addWidget(self.yLogCheckBox)


        self.rightLayout.addLayout(self.yMinMaxLayout)

        self.sp5 = QSpacerItem(0, 6, QSizePolicy.Policy.Minimum, QSizePolicy.Policy.Fixed)

        self.rightLayout.addItem(self.sp5)

        self.hline3 = QFrame(self.rightWidget)
        self.hline3.setObjectName(u"hline3")
        self.hline3.setFrameShape(QFrame.Shape.HLine)
        self.hline3.setFrameShadow(QFrame.Shadow.Sunken)

        self.rightLayout.addWidget(self.hline3)

        self.sp6 = QSpacerItem(0, 6, QSizePolicy.Policy.Minimum, QSizePolicy.Policy.Fixed)

        self.rightLayout.addItem(self.sp6)

        self.colorComboBox = QComboBox(self.rightWidget)
        self.colorComboBox.addItem("")
        self.colorComboBox.addItem("")
        self.colorComboBox.addItem("")
        self.colorComboBox.addItem("")
        self.colorComboBox.addItem("")
        self.colorComboBox.addItem("")
        self.colorComboBox.addItem("")
        self.colorComboBox.addItem("")
        self.colorComboBox.addItem("")
        self.colorComboBox.addItem("")
        self.colorComboBox.addItem("")
        self.colorComboBox.setObjectName(u"colorComboBox")
        self.colorComboBox.setVisible(False)

        self.rightLayout.addWidget(self.colorComboBox)

        self.fileMetaLabel = QLabel(self.rightWidget)
        self.fileMetaLabel.setObjectName(u"fileMetaLabel")

        self.rightLayout.addWidget(self.fileMetaLabel)

        self.fileMetaText = QLabel(self.rightWidget)
        self.fileMetaText.setObjectName(u"fileMetaText")
        self.fileMetaText.setWordWrap(True)
        self.fileMetaText.setAlignment(Qt.AlignmentFlag.AlignTop|Qt.AlignmentFlag.AlignLeft)

        self.rightLayout.addWidget(self.fileMetaText)

        self.sp9 = QSpacerItem(0, 6, QSizePolicy.Policy.Minimum, QSizePolicy.Policy.Fixed)

        self.rightLayout.addItem(self.sp9)

        self.hline5 = QFrame(self.rightWidget)
        self.hline5.setObjectName(u"hline5")
        self.hline5.setFrameShape(QFrame.Shape.HLine)
        self.hline5.setFrameShadow(QFrame.Shadow.Sunken)

        self.rightLayout.addWidget(self.hline5)

        self.sp10 = QSpacerItem(0, 6, QSizePolicy.Policy.Minimum, QSizePolicy.Policy.Fixed)

        self.rightLayout.addItem(self.sp10)

        self.notesHeaderLayout = QHBoxLayout()
        self.notesHeaderLayout.setObjectName(u"notesHeaderLayout")
        self.notesLabel = QLabel(self.rightWidget)
        self.notesLabel.setObjectName(u"notesLabel")

        self.notesHeaderLayout.addWidget(self.notesLabel)

        self.notesHeaderSpacer = QSpacerItem(0, 0, QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Minimum)

        self.notesHeaderLayout.addItem(self.notesHeaderSpacer)

        self.editNotesButton = QPushButton(self.rightWidget)
        self.editNotesButton.setObjectName(u"editNotesButton")

        self.notesHeaderLayout.addWidget(self.editNotesButton)


        self.rightLayout.addLayout(self.notesHeaderLayout)

        self.notesScrollArea = QScrollArea(self.rightWidget)
        self.notesScrollArea.setObjectName(u"notesScrollArea")
        self.notesScrollArea.setWidgetResizable(True)
        self.notesScrollArea.setFrameShape(QFrame.Shape.NoFrame)
        sizePolicy.setHeightForWidth(self.notesScrollArea.sizePolicy().hasHeightForWidth())
        self.notesScrollArea.setSizePolicy(sizePolicy)
        self.notesScrollAreaContents = QWidget()
        self.notesScrollAreaContents.setObjectName(u"notesScrollAreaContents")
        self.notesScrollAreaLayout = QVBoxLayout(self.notesScrollAreaContents)
        self.notesScrollAreaLayout.setObjectName(u"notesScrollAreaLayout")
        self.notesScrollAreaLayout.setContentsMargins(0, 0, 0, 0)
        self.notesText = QLabel(self.notesScrollAreaContents)
        self.notesText.setObjectName(u"notesText")
        self.notesText.setWordWrap(True)
        self.notesText.setAlignment(Qt.AlignmentFlag.AlignTop|Qt.AlignmentFlag.AlignLeft)

        self.notesScrollAreaLayout.addWidget(self.notesText)

        self.notesScrollArea.setWidget(self.notesScrollAreaContents)

        self.rightLayout.addWidget(self.notesScrollArea)


        self.mainLayout.addWidget(self.rightWidget)

        MainWindow.setCentralWidget(self.centralwidget)
        self.menubar = QMenuBar(MainWindow)
        self.menubar.setObjectName(u"menubar")
        self.menubar.setGeometry(QRect(0, 0, 1280, 24))
        self.menuFile = QMenu(self.menubar)
        self.menuFile.setObjectName(u"menuFile")
        self.menuSelect = QMenu(self.menubar)
        self.menuSelect.setObjectName(u"menuSelect")
        self.menuPoints = QMenu(self.menuSelect)
        self.menuPoints.setObjectName(u"menuPoints")
        self.menuTraces = QMenu(self.menuSelect)
        self.menuTraces.setObjectName(u"menuTraces")
        self.menuParameters = QMenu(self.menubar)
        self.menuParameters.setObjectName(u"menuParameters")
        self.menuRemove = QMenu(self.menuParameters)
        self.menuRemove.setObjectName(u"menuRemove")
        self.menuProcess = QMenu(self.menubar)
        self.menuProcess.setObjectName(u"menuProcess")
        self.menuPlot = QMenu(self.menubar)
        self.menuPlot.setObjectName(u"menuPlot")
        self.menuType = QMenu(self.menuPlot)
        self.menuType.setObjectName(u"menuType")
        self.menuScale = QMenu(self.menuPlot)
        self.menuScale.setObjectName(u"menuScale")
        self.menuColor = QMenu(self.menuPlot)
        self.menuColor.setObjectName(u"menuColor")
        MainWindow.setMenuBar(self.menubar)
        self.statusbar = QStatusBar(MainWindow)
        self.statusbar.setObjectName(u"statusbar")
        MainWindow.setStatusBar(self.statusbar)

        self.menubar.addAction(self.menuFile.menuAction())
        self.menubar.addAction(self.menuSelect.menuAction())
        self.menubar.addAction(self.menuParameters.menuAction())
        self.menubar.addAction(self.menuProcess.menuAction())
        self.menubar.addAction(self.menuPlot.menuAction())
        self.menuFile.addAction(self.actionOpen)
        self.menuFile.addAction(self.actionImportLegacy)
        self.menuFile.addSeparator()
        self.menuFile.addAction(self.actionCloseCurrent)
        self.menuFile.addAction(self.actionCloseAll)
        self.menuFile.addSeparator()
        self.menuFile.addAction(self.actionWatchFolder)
        self.menuFile.addAction(self.actionWatchingEnabled)
        self.menuFile.addSeparator()
        self.menuFile.addAction(self.actionExport)
        self.menuFile.addSeparator()
        self.menuFile.addAction(self.actionInspectH5)
        self.menuFile.addSeparator()
        self.menuFile.addAction(self.actionExit)
        self.menuSelect.addAction(self.menuPoints.menuAction())
        self.menuSelect.addAction(self.menuTraces.menuAction())
        self.menuPoints.addAction(self.actionPointsInRange)
        self.menuPoints.addAction(self.actionPointsInView)
        self.menuTraces.addAction(self.actionAll)
        self.menuTraces.addAction(self.actionNone)
        self.menuTraces.addAction(self.actionToggle)
        self.menuTraces.addAction(self.actionCriteria)
        self.menuParameters.addAction(self.actionFit)
        self.menuParameters.addAction(self.actionPlotHistogram)
        self.menuParameters.addSeparator()
        self.menuParameters.addAction(self.actionAutoFitOnOpen)
        self.menuParameters.addSeparator()
        self.menuParameters.addAction(self.menuRemove.menuAction())
        self.menuRemove.addAction(self.actionRemoveModelViewed)
        self.menuRemove.addAction(self.actionRemoveModelAll)
        self.menuRemove.addSeparator()
        self.menuRemove.addAction(self.actionRemoveAllProcessing)
        self.menuProcess.addAction(self.actionCalculateForce)
        self.menuProcess.addAction(self.actionFreeDNA)
        self.menuProcess.addAction(self.actionRuptures)
        self.menuPlot.addAction(self.menuType.menuAction())
        self.menuPlot.addAction(self.menuScale.menuAction())
        self.menuPlot.addAction(self.menuColor.menuAction())
        self.menuPlot.addSeparator()
        self.menuPlot.addAction(self.actionCopyImage)
        self.menuPlot.addAction(self.actionSaveImage)
        self.menuPlot.addAction(self.actionBulkSaveImages)
        self.menuType.addAction(self.actionTypeZt)
        self.menuType.addAction(self.actionTypeFZ)
        self.menuType.addAction(self.actionTypeYX)
        self.menuScale.addAction(self.actionScaleAuto)
        self.menuScale.addAction(self.actionScaleFix)
        self.menuScale.addAction(self.actionScaleContour)
        self.menuColor.addAction(self.actionColorTheme)
        self.menuColor.addAction(self.actionColorBlue)
        self.menuColor.addAction(self.actionColorOrange)
        self.menuColor.addAction(self.actionColorGreen)
        self.menuColor.addAction(self.actionColorPink)
        self.menuColor.addAction(self.actionColorLightblue)
        self.menuColor.addAction(self.actionColorPurple)
        self.menuColor.addAction(self.actionColorFuchia)
        self.menuColor.addAction(self.actionColorYellow)
        self.menuColor.addAction(self.actionColorTeal)
        self.menuColor.addAction(self.actionColorRed)

        self.retranslateUi(MainWindow)

        QMetaObject.connectSlotsByName(MainWindow)
    # setupUi

    def retranslateUi(self, MainWindow):
        MainWindow.setWindowTitle(QCoreApplication.translate("MainWindow", u"Trace Editor", None))
        self.actionOpen.setText(QCoreApplication.translate("MainWindow", u"Open\u2026", None))
#if QT_CONFIG(shortcut)
        self.actionOpen.setShortcut(QCoreApplication.translate("MainWindow", u"Ctrl+O", None))
#endif // QT_CONFIG(shortcut)
        self.actionImportLegacy.setText(QCoreApplication.translate("MainWindow", u"Import Legacy .hdf\u2026", None))
        self.actionInspectH5.setText(QCoreApplication.translate("MainWindow", u"Inspect H5\u2026", None))
        self.actionCloseCurrent.setText(QCoreApplication.translate("MainWindow", u"Close Current", None))
        self.actionCloseAll.setText(QCoreApplication.translate("MainWindow", u"Close All", None))
        self.actionExport.setText(QCoreApplication.translate("MainWindow", u"Export\u2026", None))
#if QT_CONFIG(shortcut)
        self.actionExport.setShortcut(QCoreApplication.translate("MainWindow", u"Ctrl+E", None))
#endif // QT_CONFIG(shortcut)
        self.actionExit.setText(QCoreApplication.translate("MainWindow", u"Exit", None))
#if QT_CONFIG(shortcut)
        self.actionExit.setShortcut(QCoreApplication.translate("MainWindow", u"Ctrl+W", None))
#endif // QT_CONFIG(shortcut)
        self.actionWatchFolder.setText(QCoreApplication.translate("MainWindow", u"Watch Folder\u2026", None))
        self.actionWatchingEnabled.setText(QCoreApplication.translate("MainWindow", u"Watching", None))
        self.actionPointsInRange.setText(QCoreApplication.translate("MainWindow", u"By Criteria\u2026", None))
        self.actionPointsInView.setText(QCoreApplication.translate("MainWindow", u"In view", None))
        self.actionAll.setText(QCoreApplication.translate("MainWindow", u"All", None))
        self.actionNone.setText(QCoreApplication.translate("MainWindow", u"None", None))
        self.actionToggle.setText(QCoreApplication.translate("MainWindow", u"Toggle", None))
        self.actionCriteria.setText(QCoreApplication.translate("MainWindow", u"Criteria\u2026", None))
        self.actionFit.setText(QCoreApplication.translate("MainWindow", u"Fit\u2026", None))
        self.actionPlotHistogram.setText(QCoreApplication.translate("MainWindow", u"Plot Histogram", None))
        self.actionAutoFitOnOpen.setText(QCoreApplication.translate("MainWindow", u"Auto-fit on File Open", None))
        self.actionRemoveModelViewed.setText(QCoreApplication.translate("MainWindow", u"Model from Viewed Trace", None))
        self.actionRemoveModelAll.setText(QCoreApplication.translate("MainWindow", u"Model from All Traces", None))
        self.actionRemoveAllProcessing.setText(QCoreApplication.translate("MainWindow", u"All Processing (Params + Force)", None))
        self.actionCalculateForce.setText(QCoreApplication.translate("MainWindow", u"Calculate Force", None))
        self.actionFreeDNA.setText(QCoreApplication.translate("MainWindow", u"Free DNA", None))
        self.actionRuptures.setText(QCoreApplication.translate("MainWindow", u"Ruptures", None))
        self.actionTypeZt.setText(QCoreApplication.translate("MainWindow", u"Z(t)", None))
        self.actionTypeFZ.setText(QCoreApplication.translate("MainWindow", u"F(Z)", None))
        self.actionTypeYX.setText(QCoreApplication.translate("MainWindow", u"Y(X)", None))
        self.actionScaleAuto.setText(QCoreApplication.translate("MainWindow", u"Auto", None))
        self.actionScaleFix.setText(QCoreApplication.translate("MainWindow", u"Fix", None))
        self.actionColorTheme.setText(QCoreApplication.translate("MainWindow", u"theme", None))
        self.actionColorBlue.setText(QCoreApplication.translate("MainWindow", u"blue", None))
        self.actionColorOrange.setText(QCoreApplication.translate("MainWindow", u"orange", None))
        self.actionColorGreen.setText(QCoreApplication.translate("MainWindow", u"green", None))
        self.actionColorPink.setText(QCoreApplication.translate("MainWindow", u"pink", None))
        self.actionColorLightblue.setText(QCoreApplication.translate("MainWindow", u"lightblue", None))
        self.actionColorPurple.setText(QCoreApplication.translate("MainWindow", u"purple", None))
        self.actionColorFuchia.setText(QCoreApplication.translate("MainWindow", u"fuchia", None))
        self.actionColorYellow.setText(QCoreApplication.translate("MainWindow", u"yellow", None))
        self.actionColorTeal.setText(QCoreApplication.translate("MainWindow", u"teal", None))
        self.actionColorRed.setText(QCoreApplication.translate("MainWindow", u"red", None))
        self.actionScaleContour.setText(QCoreApplication.translate("MainWindow", u"Contour length", None))
        self.actionCopyImage.setText(QCoreApplication.translate("MainWindow", u"Copy Image", None))
#if QT_CONFIG(shortcut)
        self.actionCopyImage.setShortcut(QCoreApplication.translate("MainWindow", u"Ctrl+Shift+C", None))
#endif // QT_CONFIG(shortcut)
        self.actionSaveImage.setText(QCoreApplication.translate("MainWindow", u"Save Image\u2026", None))
#if QT_CONFIG(shortcut)
        self.actionSaveImage.setShortcut(QCoreApplication.translate("MainWindow", u"Ctrl+Shift+S", None))
#endif // QT_CONFIG(shortcut)
        self.actionBulkSaveImages.setText(QCoreApplication.translate("MainWindow", u"Bulk Save Images\u2026", None))
        self.modelLabel.setText(QCoreApplication.translate("MainWindow", u"Model", None))
        self.modelComboBox.setItemText(0, QCoreApplication.translate("MainWindow", u"none", None))
        self.modelComboBox.setItemText(1, QCoreApplication.translate("MainWindow", u"ExponentialDrift", None))
        self.modelComboBox.setItemText(2, QCoreApplication.translate("MainWindow", u"zMagnetForce", None))
        self.modelComboBox.setItemText(3, QCoreApplication.translate("MainWindow", u"zWLC", None))
        self.modelComboBox.setItemText(4, QCoreApplication.translate("MainWindow", u"zFiber", None))
        self.modelComboBox.setItemText(5, QCoreApplication.translate("MainWindow", u"zHMf", None))
        self.modelComboBox.setItemText(6, QCoreApplication.translate("MainWindow", u"zHMf2", None))

        self.parametersLabel.setText(QCoreApplication.translate("MainWindow", u"Parameters", None))
        self.paramsListWidget.setStyleSheet(QCoreApplication.translate("MainWindow", u"QListWidget::item { padding: 2px; }", None))
        self.pipelineLabel.setText(QCoreApplication.translate("MainWindow", u"Pipeline", None))
        self.runPipelineButton.setText(QCoreApplication.translate("MainWindow", u"Run", None))
        self.filenameLabel.setText(QCoreApplication.translate("MainWindow", u"Filename:", None))
        self.roiLabel.setText(QCoreApplication.translate("MainWindow", u"ROI:", None))
        self.selectedButton.setText(QCoreApplication.translate("MainWindow", u"ROI selected", None))
        self.selectedOnlyCheckBox.setText(QCoreApplication.translate("MainWindow", u"Show Selected Only", None))
        self.xAxisLabel.setText(QCoreApplication.translate("MainWindow", u"X axis", None))
        self.xMinLabel.setText(QCoreApplication.translate("MainWindow", u"min", None))
        self.xMinEdit.setStyleSheet(QCoreApplication.translate("MainWindow", u"QLineEdit { padding: 2px 3px; }", None))
        self.xMaxLabel.setText(QCoreApplication.translate("MainWindow", u"max", None))
        self.xMaxEdit.setStyleSheet(QCoreApplication.translate("MainWindow", u"QLineEdit { padding: 2px 3px; }", None))
        self.xLogCheckBox.setText(QCoreApplication.translate("MainWindow", u"Log", None))
        self.yAxisLabel.setText(QCoreApplication.translate("MainWindow", u"Y axis", None))
        self.yMinLabel.setText(QCoreApplication.translate("MainWindow", u"min", None))
        self.yMinEdit.setStyleSheet(QCoreApplication.translate("MainWindow", u"QLineEdit { padding: 2px 3px; }", None))
        self.yMaxLabel.setText(QCoreApplication.translate("MainWindow", u"max", None))
        self.yMaxEdit.setStyleSheet(QCoreApplication.translate("MainWindow", u"QLineEdit { padding: 2px 3px; }", None))
        self.yLogCheckBox.setText(QCoreApplication.translate("MainWindow", u"Log", None))
        self.colorComboBox.setItemText(0, QCoreApplication.translate("MainWindow", u"theme", None))
        self.colorComboBox.setItemText(1, QCoreApplication.translate("MainWindow", u"blue", None))
        self.colorComboBox.setItemText(2, QCoreApplication.translate("MainWindow", u"orange", None))
        self.colorComboBox.setItemText(3, QCoreApplication.translate("MainWindow", u"green", None))
        self.colorComboBox.setItemText(4, QCoreApplication.translate("MainWindow", u"pink", None))
        self.colorComboBox.setItemText(5, QCoreApplication.translate("MainWindow", u"lightblue", None))
        self.colorComboBox.setItemText(6, QCoreApplication.translate("MainWindow", u"purple", None))
        self.colorComboBox.setItemText(7, QCoreApplication.translate("MainWindow", u"fuchia", None))
        self.colorComboBox.setItemText(8, QCoreApplication.translate("MainWindow", u"yellow", None))
        self.colorComboBox.setItemText(9, QCoreApplication.translate("MainWindow", u"teal", None))
        self.colorComboBox.setItemText(10, QCoreApplication.translate("MainWindow", u"red", None))

        self.fileMetaLabel.setText(QCoreApplication.translate("MainWindow", u"Properties", None))
        self.fileMetaText.setText("")
        self.notesLabel.setText(QCoreApplication.translate("MainWindow", u"Notes", None))
        self.editNotesButton.setText(QCoreApplication.translate("MainWindow", u"Edit\u2026", None))
        self.notesText.setText("")
        self.menuFile.setTitle(QCoreApplication.translate("MainWindow", u"File", None))
        self.menuSelect.setTitle(QCoreApplication.translate("MainWindow", u"Select", None))
        self.menuPoints.setTitle(QCoreApplication.translate("MainWindow", u"Points", None))
        self.menuTraces.setTitle(QCoreApplication.translate("MainWindow", u"Traces", None))
        self.menuParameters.setTitle(QCoreApplication.translate("MainWindow", u"Parameters", None))
        self.menuRemove.setTitle(QCoreApplication.translate("MainWindow", u"Remove", None))
        self.menuProcess.setTitle(QCoreApplication.translate("MainWindow", u"Process", None))
        self.menuPlot.setTitle(QCoreApplication.translate("MainWindow", u"Plot", None))
        self.menuType.setTitle(QCoreApplication.translate("MainWindow", u"Type", None))
        self.menuScale.setTitle(QCoreApplication.translate("MainWindow", u"Scale", None))
        self.menuColor.setTitle(QCoreApplication.translate("MainWindow", u"Color", None))
    # retranslateUi

