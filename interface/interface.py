# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'interface.ui'
#
# Created by: PyQt5 UI code generator 5.9.2
#
# WARNING! All changes made in this file will be lost!

from PyQt5 import QtCore, QtGui, QtWidgets

class Ui_MainWindow(object):
    def setupUi(self, MainWindow):
        MainWindow.setObjectName("MainWindow")
        MainWindow.resize(1024, 680)
        self.centralwidget = QtWidgets.QWidget(MainWindow)
        self.centralwidget.setObjectName("centralwidget")
        self.gridLayout_28 = QtWidgets.QGridLayout(self.centralwidget)
        self.gridLayout_28.setObjectName("gridLayout_28")
        self.tabWidget = QtWidgets.QTabWidget(self.centralwidget)
        self.tabWidget.setMinimumSize(QtCore.QSize(500, 0))
        self.tabWidget.setObjectName("tabWidget")
        self.tab_6 = QtWidgets.QWidget()
        self.tab_6.setObjectName("tab_6")
        self.gridLayout_9 = QtWidgets.QGridLayout(self.tab_6)
        self.gridLayout_9.setObjectName("gridLayout_9")
        self.gridLayout_2 = QtWidgets.QGridLayout()
        self.gridLayout_2.setObjectName("gridLayout_2")
        self.spin_list = QtWidgets.QTextEdit(self.tab_6)
        self.spin_list.setObjectName("spin_list")
        self.gridLayout_2.addWidget(self.spin_list, 4, 0, 1, 1)
        self.label = QtWidgets.QLabel(self.tab_6)
        self.label.setObjectName("label")
        self.gridLayout_2.addWidget(self.label, 0, 0, 1, 1)
        self.spin_combination = QtWidgets.QComboBox(self.tab_6)
        self.spin_combination.setObjectName("spin_combination")
        self.spin_combination.addItem("")
        self.spin_combination.addItem("")
        self.spin_combination.addItem("")
        self.spin_combination.addItem("")
        self.spin_combination.addItem("")
        self.gridLayout_2.addWidget(self.spin_combination, 2, 1, 1, 1)
        self.nsites = QtWidgets.QLineEdit(self.tab_6)
        self.nsites.setObjectName("nsites")
        self.gridLayout_2.addWidget(self.nsites, 0, 1, 1, 1)
        self.label_5 = QtWidgets.QLabel(self.tab_6)
        self.label_5.setObjectName("label_5")
        self.gridLayout_2.addWidget(self.label_5, 3, 0, 1, 1)
        self.initialize_spins = QtWidgets.QPushButton(self.tab_6)
        self.initialize_spins.setObjectName("initialize_spins")
        self.gridLayout_2.addWidget(self.initialize_spins, 2, 0, 1, 1)
        self.gridLayout_9.addLayout(self.gridLayout_2, 0, 0, 1, 1)
        self.tabWidget.addTab(self.tab_6, "")
        self.tab_5 = QtWidgets.QWidget()
        self.tab_5.setObjectName("tab_5")
        self.gridLayout_4 = QtWidgets.QGridLayout(self.tab_5)
        self.gridLayout_4.setObjectName("gridLayout_4")
        self.tabWidget_5 = QtWidgets.QTabWidget(self.tab_5)
        self.tabWidget_5.setObjectName("tabWidget_5")
        self.tab = QtWidgets.QWidget()
        self.tab.setObjectName("tab")
        self.gridLayout_20 = QtWidgets.QGridLayout(self.tab)
        self.gridLayout_20.setObjectName("gridLayout_20")
        self.gridLayout_19 = QtWidgets.QGridLayout()
        self.gridLayout_19.setObjectName("gridLayout_19")
        self.label_19 = QtWidgets.QLabel(self.tab)
        self.label_19.setObjectName("label_19")
        self.gridLayout_19.addWidget(self.label_19, 0, 2, 1, 1)
        self.label_18 = QtWidgets.QLabel(self.tab)
        self.label_18.setObjectName("label_18")
        self.gridLayout_19.addWidget(self.label_18, 0, 1, 1, 1)
        self.jzz_table = QtWidgets.QTextEdit(self.tab)
        self.jzz_table.setObjectName("jzz_table")
        self.gridLayout_19.addWidget(self.jzz_table, 1, 2, 1, 1)
        self.jyy_table = QtWidgets.QTextEdit(self.tab)
        self.jyy_table.setObjectName("jyy_table")
        self.gridLayout_19.addWidget(self.jyy_table, 1, 1, 1, 1)
        self.label_17 = QtWidgets.QLabel(self.tab)
        self.label_17.setObjectName("label_17")
        self.gridLayout_19.addWidget(self.label_17, 0, 0, 1, 1)
        self.jxx_table = QtWidgets.QTextEdit(self.tab)
        self.jxx_table.setObjectName("jxx_table")
        self.gridLayout_19.addWidget(self.jxx_table, 1, 0, 1, 1)
        self.gridLayout_20.addLayout(self.gridLayout_19, 0, 0, 1, 1)
        self.tabWidget_5.addTab(self.tab, "")
        self.tab_14 = QtWidgets.QWidget()
        self.tab_14.setObjectName("tab_14")
        self.gridLayout_23 = QtWidgets.QGridLayout(self.tab_14)
        self.gridLayout_23.setObjectName("gridLayout_23")
        self.gridLayout_3 = QtWidgets.QGridLayout()
        self.gridLayout_3.setObjectName("gridLayout_3")
        self.label_2 = QtWidgets.QLabel(self.tab_14)
        self.label_2.setObjectName("label_2")
        self.gridLayout_3.addWidget(self.label_2, 0, 0, 1, 1)
        self.label_8 = QtWidgets.QLabel(self.tab_14)
        self.label_8.setObjectName("label_8")
        self.gridLayout_3.addWidget(self.label_8, 1, 0, 1, 1)
        self.gridLayout_8 = QtWidgets.QGridLayout()
        self.gridLayout_8.setObjectName("gridLayout_8")
        self.Jxx = QtWidgets.QLineEdit(self.tab_14)
        self.Jxx.setObjectName("Jxx")
        self.gridLayout_8.addWidget(self.Jxx, 0, 1, 1, 1)
        self.Jyy = QtWidgets.QLineEdit(self.tab_14)
        self.Jyy.setObjectName("Jyy")
        self.gridLayout_8.addWidget(self.Jyy, 1, 1, 1, 1)
        self.Jzz = QtWidgets.QLineEdit(self.tab_14)
        self.Jzz.setObjectName("Jzz")
        self.gridLayout_8.addWidget(self.Jzz, 2, 1, 1, 1)
        self.label_3 = QtWidgets.QLabel(self.tab_14)
        self.label_3.setObjectName("label_3")
        self.gridLayout_8.addWidget(self.label_3, 1, 0, 1, 1)
        self.label_9 = QtWidgets.QLabel(self.tab_14)
        self.label_9.setMinimumSize(QtCore.QSize(40, 0))
        self.label_9.setObjectName("label_9")
        self.gridLayout_8.addWidget(self.label_9, 0, 0, 1, 1)
        self.label_21 = QtWidgets.QLabel(self.tab_14)
        self.label_21.setObjectName("label_21")
        self.gridLayout_8.addWidget(self.label_21, 2, 0, 1, 1)
        self.gridLayout_3.addLayout(self.gridLayout_8, 0, 1, 1, 1)
        self.dimerization_exchange = QtWidgets.QLineEdit(self.tab_14)
        self.dimerization_exchange.setObjectName("dimerization_exchange")
        self.gridLayout_3.addWidget(self.dimerization_exchange, 1, 1, 1, 1)
        self.gridLayout_23.addLayout(self.gridLayout_3, 0, 0, 1, 1)
        self.tabWidget_5.addTab(self.tab_14, "")
        self.gridLayout_4.addWidget(self.tabWidget_5, 0, 0, 1, 1)
        self.tabWidget.addTab(self.tab_5, "")
        self.tab_9 = QtWidgets.QWidget()
        self.tab_9.setObjectName("tab_9")
        self.gridLayout_25 = QtWidgets.QGridLayout(self.tab_9)
        self.gridLayout_25.setObjectName("gridLayout_25")
        self.tabWidget_3 = QtWidgets.QTabWidget(self.tab_9)
        self.tabWidget_3.setObjectName("tabWidget_3")
        self.tab_7 = QtWidgets.QWidget()
        self.tab_7.setObjectName("tab_7")
        self.gridLayout_24 = QtWidgets.QGridLayout(self.tab_7)
        self.gridLayout_24.setObjectName("gridLayout_24")
        self.gridLayout_7 = QtWidgets.QGridLayout()
        self.gridLayout_7.setObjectName("gridLayout_7")
        self.bz_table = QtWidgets.QTextEdit(self.tab_7)
        self.bz_table.setObjectName("bz_table")
        self.gridLayout_7.addWidget(self.bz_table, 1, 2, 1, 1)
        self.bx_table = QtWidgets.QTextEdit(self.tab_7)
        self.bx_table.setObjectName("bx_table")
        self.gridLayout_7.addWidget(self.bx_table, 1, 0, 1, 1)
        self.by_table = QtWidgets.QTextEdit(self.tab_7)
        self.by_table.setObjectName("by_table")
        self.gridLayout_7.addWidget(self.by_table, 1, 1, 1, 1)
        self.label_22 = QtWidgets.QLabel(self.tab_7)
        self.label_22.setObjectName("label_22")
        self.gridLayout_7.addWidget(self.label_22, 0, 0, 1, 1)
        self.label_23 = QtWidgets.QLabel(self.tab_7)
        self.label_23.setObjectName("label_23")
        self.gridLayout_7.addWidget(self.label_23, 0, 1, 1, 1)
        self.label_24 = QtWidgets.QLabel(self.tab_7)
        self.label_24.setObjectName("label_24")
        self.gridLayout_7.addWidget(self.label_24, 0, 2, 1, 1)
        self.gridLayout_24.addLayout(self.gridLayout_7, 0, 0, 1, 1)
        self.tabWidget_3.addTab(self.tab_7, "")
        self.tab_8 = QtWidgets.QWidget()
        self.tab_8.setObjectName("tab_8")
        self.gridLayout_27 = QtWidgets.QGridLayout(self.tab_8)
        self.gridLayout_27.setObjectName("gridLayout_27")
        self.gridLayout_6 = QtWidgets.QGridLayout()
        self.gridLayout_6.setObjectName("gridLayout_6")
        self.label_4 = QtWidgets.QLabel(self.tab_8)
        self.label_4.setObjectName("label_4")
        self.gridLayout_6.addWidget(self.label_4, 2, 0, 1, 1)
        self.modulation_field = QtWidgets.QComboBox(self.tab_8)
        self.modulation_field.setObjectName("modulation_field")
        self.modulation_field.addItem("")
        self.modulation_field.addItem("")
        self.gridLayout_6.addWidget(self.modulation_field, 2, 1, 1, 1)
        self.label_25 = QtWidgets.QLabel(self.tab_8)
        self.label_25.setObjectName("label_25")
        self.gridLayout_6.addWidget(self.label_25, 0, 0, 1, 1)
        self.gridLayout_26 = QtWidgets.QGridLayout()
        self.gridLayout_26.setObjectName("gridLayout_26")
        self.label_26 = QtWidgets.QLabel(self.tab_8)
        self.label_26.setMinimumSize(QtCore.QSize(40, 0))
        self.label_26.setObjectName("label_26")
        self.gridLayout_26.addWidget(self.label_26, 0, 0, 1, 1)
        self.label_27 = QtWidgets.QLabel(self.tab_8)
        self.label_27.setObjectName("label_27")
        self.gridLayout_26.addWidget(self.label_27, 1, 0, 1, 1)
        self.label_28 = QtWidgets.QLabel(self.tab_8)
        self.label_28.setObjectName("label_28")
        self.gridLayout_26.addWidget(self.label_28, 2, 0, 1, 1)
        self.Bx = QtWidgets.QLineEdit(self.tab_8)
        self.Bx.setObjectName("Bx")
        self.gridLayout_26.addWidget(self.Bx, 0, 1, 1, 1)
        self.By = QtWidgets.QLineEdit(self.tab_8)
        self.By.setObjectName("By")
        self.gridLayout_26.addWidget(self.By, 1, 1, 1, 1)
        self.Bz = QtWidgets.QLineEdit(self.tab_8)
        self.Bz.setObjectName("Bz")
        self.gridLayout_26.addWidget(self.Bz, 2, 1, 1, 1)
        self.gridLayout_6.addLayout(self.gridLayout_26, 0, 1, 1, 1)
        self.label_29 = QtWidgets.QLabel(self.tab_8)
        self.label_29.setObjectName("label_29")
        self.gridLayout_6.addWidget(self.label_29, 1, 0, 1, 1)
        self.location_field = QtWidgets.QComboBox(self.tab_8)
        self.location_field.setObjectName("location_field")
        self.location_field.addItem("")
        self.location_field.addItem("")
        self.location_field.addItem("")
        self.gridLayout_6.addWidget(self.location_field, 1, 1, 1, 1)
        self.gridLayout_27.addLayout(self.gridLayout_6, 0, 0, 1, 1)
        self.tabWidget_3.addTab(self.tab_8, "")
        self.gridLayout_25.addWidget(self.tabWidget_3, 0, 0, 1, 1)
        self.tabWidget.addTab(self.tab_9, "")
        self.tab_10 = QtWidgets.QWidget()
        self.tab_10.setObjectName("tab_10")
        self.gridLayout_29 = QtWidgets.QGridLayout(self.tab_10)
        self.gridLayout_29.setObjectName("gridLayout_29")
        self.gridLayout_11 = QtWidgets.QGridLayout()
        self.gridLayout_11.setObjectName("gridLayout_11")
        self.maxm = QtWidgets.QLineEdit(self.tab_10)
        self.maxm.setObjectName("maxm")
        self.gridLayout_11.addWidget(self.maxm, 0, 1, 1, 1)
        self.nsweeps = QtWidgets.QLineEdit(self.tab_10)
        self.nsweeps.setObjectName("nsweeps")
        self.gridLayout_11.addWidget(self.nsweeps, 1, 1, 1, 1)
        self.label_10 = QtWidgets.QLabel(self.tab_10)
        self.label_10.setObjectName("label_10")
        self.gridLayout_11.addWidget(self.label_10, 0, 0, 1, 1)
        self.label_11 = QtWidgets.QLabel(self.tab_10)
        self.label_11.setObjectName("label_11")
        self.gridLayout_11.addWidget(self.label_11, 1, 0, 1, 1)
        self.gridLayout_29.addLayout(self.gridLayout_11, 0, 0, 1, 1)
        self.tabWidget.addTab(self.tab_10, "")
        self.gridLayout_28.addWidget(self.tabWidget, 0, 0, 1, 1)
        self.tabWidget_2 = QtWidgets.QTabWidget(self.centralwidget)
        self.tabWidget_2.setMinimumSize(QtCore.QSize(500, 0))
        self.tabWidget_2.setObjectName("tabWidget_2")
        self.tab_2 = QtWidgets.QWidget()
        self.tab_2.setObjectName("tab_2")
        self.gridLayout_22 = QtWidgets.QGridLayout(self.tab_2)
        self.gridLayout_22.setObjectName("gridLayout_22")
        self.gridLayout_21 = QtWidgets.QGridLayout()
        self.gridLayout_21.setObjectName("gridLayout_21")
        self.label_20 = QtWidgets.QLabel(self.tab_2)
        self.label_20.setObjectName("label_20")
        self.gridLayout_21.addWidget(self.label_20, 0, 0, 1, 1)
        self.magnetization_type = QtWidgets.QComboBox(self.tab_2)
        self.magnetization_type.setObjectName("magnetization_type")
        self.magnetization_type.addItem("")
        self.magnetization_type.addItem("")
        self.magnetization_type.addItem("")
        self.gridLayout_21.addWidget(self.magnetization_type, 0, 1, 1, 1)
        self.get_magnetization = QtWidgets.QPushButton(self.tab_2)
        self.get_magnetization.setObjectName("get_magnetization")
        self.gridLayout_21.addWidget(self.get_magnetization, 1, 1, 1, 1)
        self.gridLayout_22.addLayout(self.gridLayout_21, 0, 0, 1, 1)
        self.tabWidget_2.addTab(self.tab_2, "")
        self.tab_3 = QtWidgets.QWidget()
        self.tab_3.setObjectName("tab_3")
        self.gridLayout = QtWidgets.QGridLayout(self.tab_3)
        self.gridLayout.setObjectName("gridLayout")
        self.gridLayout_10 = QtWidgets.QGridLayout()
        self.gridLayout_10.setObjectName("gridLayout_10")
        self.label_6 = QtWidgets.QLabel(self.tab_3)
        self.label_6.setObjectName("label_6")
        self.gridLayout_10.addWidget(self.label_6, 0, 0, 1, 1)
        self.static_arrangement = QtWidgets.QComboBox(self.tab_3)
        self.static_arrangement.setObjectName("static_arrangement")
        self.static_arrangement.addItem("")
        self.static_arrangement.addItem("")
        self.static_arrangement.addItem("")
        self.static_arrangement.addItem("")
        self.gridLayout_10.addWidget(self.static_arrangement, 0, 1, 1, 1)
        self.label_7 = QtWidgets.QLabel(self.tab_3)
        self.label_7.setObjectName("label_7")
        self.gridLayout_10.addWidget(self.label_7, 1, 0, 1, 1)
        self.static_type = QtWidgets.QComboBox(self.tab_3)
        self.static_type.setObjectName("static_type")
        self.static_type.addItem("")
        self.static_type.addItem("")
        self.static_type.addItem("")
        self.static_type.addItem("")
        self.gridLayout_10.addWidget(self.static_type, 1, 1, 1, 1)
        self.get_static_correlator = QtWidgets.QPushButton(self.tab_3)
        self.get_static_correlator.setObjectName("get_static_correlator")
        self.gridLayout_10.addWidget(self.get_static_correlator, 2, 0, 1, 1)
        self.gridLayout.addLayout(self.gridLayout_10, 0, 0, 1, 1)
        self.tabWidget_2.addTab(self.tab_3, "")
        self.tab_4 = QtWidgets.QWidget()
        self.tab_4.setObjectName("tab_4")
        self.gridLayout_12 = QtWidgets.QGridLayout(self.tab_4)
        self.gridLayout_12.setObjectName("gridLayout_12")
        self.tabWidget_4 = QtWidgets.QTabWidget(self.tab_4)
        self.tabWidget_4.setObjectName("tabWidget_4")
        self.tab_11 = QtWidgets.QWidget()
        self.tab_11.setObjectName("tab_11")
        self.gridLayout_18 = QtWidgets.QGridLayout(self.tab_11)
        self.gridLayout_18.setObjectName("gridLayout_18")
        self.gridLayout_15 = QtWidgets.QGridLayout()
        self.gridLayout_15.setObjectName("gridLayout_15")
        self.dynamic_site_i_single = QtWidgets.QLineEdit(self.tab_11)
        self.dynamic_site_i_single.setObjectName("dynamic_site_i_single")
        self.gridLayout_15.addWidget(self.dynamic_site_i_single, 0, 1, 1, 1)
        self.label_15 = QtWidgets.QLabel(self.tab_11)
        self.label_15.setObjectName("label_15")
        self.gridLayout_15.addWidget(self.label_15, 0, 0, 1, 1)
        self.get_dynamical_correlator_single = QtWidgets.QPushButton(self.tab_11)
        self.get_dynamical_correlator_single.setObjectName("get_dynamical_correlator_single")
        self.gridLayout_15.addWidget(self.get_dynamical_correlator_single, 2, 1, 1, 1)
        self.dynamic_type_single = QtWidgets.QComboBox(self.tab_11)
        self.dynamic_type_single.setObjectName("dynamic_type_single")
        self.dynamic_type_single.addItem("")
        self.dynamic_type_single.addItem("")
        self.dynamic_type_single.addItem("")
        self.gridLayout_15.addWidget(self.dynamic_type_single, 1, 1, 1, 1)
        self.label_16 = QtWidgets.QLabel(self.tab_11)
        self.label_16.setObjectName("label_16")
        self.gridLayout_15.addWidget(self.label_16, 1, 0, 1, 1)
        self.gridLayout_18.addLayout(self.gridLayout_15, 0, 0, 1, 1)
        self.tabWidget_4.addTab(self.tab_11, "")
        self.tab_12 = QtWidgets.QWidget()
        self.tab_12.setObjectName("tab_12")
        self.gridLayout_16 = QtWidgets.QGridLayout(self.tab_12)
        self.gridLayout_16.setObjectName("gridLayout_16")
        self.gridLayout_13 = QtWidgets.QGridLayout()
        self.gridLayout_13.setObjectName("gridLayout_13")
        self.comboBox_2 = QtWidgets.QComboBox(self.tab_12)
        self.comboBox_2.setObjectName("comboBox_2")
        self.comboBox_2.addItem("")
        self.comboBox_2.addItem("")
        self.comboBox_2.addItem("")
        self.gridLayout_13.addWidget(self.comboBox_2, 0, 1, 1, 1)
        self.label_12 = QtWidgets.QLabel(self.tab_12)
        self.label_12.setObjectName("label_12")
        self.gridLayout_13.addWidget(self.label_12, 0, 0, 1, 1)
        self.gridLayout_16.addLayout(self.gridLayout_13, 0, 0, 1, 1)
        self.tabWidget_4.addTab(self.tab_12, "")
        self.tab_13 = QtWidgets.QWidget()
        self.tab_13.setObjectName("tab_13")
        self.gridLayout_17 = QtWidgets.QGridLayout(self.tab_13)
        self.gridLayout_17.setObjectName("gridLayout_17")
        self.gridLayout_14 = QtWidgets.QGridLayout()
        self.gridLayout_14.setObjectName("gridLayout_14")
        self.comboBox_6 = QtWidgets.QComboBox(self.tab_13)
        self.comboBox_6.setObjectName("comboBox_6")
        self.comboBox_6.addItem("")
        self.comboBox_6.addItem("")
        self.gridLayout_14.addWidget(self.comboBox_6, 1, 1, 1, 1)
        self.smearing_dynamic = QtWidgets.QLineEdit(self.tab_13)
        self.smearing_dynamic.setObjectName("smearing_dynamic")
        self.gridLayout_14.addWidget(self.smearing_dynamic, 0, 1, 1, 1)
        self.label_13 = QtWidgets.QLabel(self.tab_13)
        self.label_13.setObjectName("label_13")
        self.gridLayout_14.addWidget(self.label_13, 0, 0, 1, 1)
        self.label_14 = QtWidgets.QLabel(self.tab_13)
        self.label_14.setObjectName("label_14")
        self.gridLayout_14.addWidget(self.label_14, 1, 0, 1, 1)
        self.gridLayout_17.addLayout(self.gridLayout_14, 0, 0, 1, 1)
        self.tabWidget_4.addTab(self.tab_13, "")
        self.gridLayout_12.addWidget(self.tabWidget_4, 0, 0, 1, 1)
        self.tabWidget_2.addTab(self.tab_4, "")
        self.gridLayout_28.addWidget(self.tabWidget_2, 0, 1, 1, 1)
        MainWindow.setCentralWidget(self.centralwidget)
        self.menubar = QtWidgets.QMenuBar(MainWindow)
        self.menubar.setGeometry(QtCore.QRect(0, 0, 1024, 25))
        self.menubar.setObjectName("menubar")
        MainWindow.setMenuBar(self.menubar)
        self.statusbar = QtWidgets.QStatusBar(MainWindow)
        self.statusbar.setObjectName("statusbar")
        MainWindow.setStatusBar(self.statusbar)

        self.retranslateUi(MainWindow)
        self.tabWidget.setCurrentIndex(0)
        self.tabWidget_5.setCurrentIndex(0)
        self.tabWidget_3.setCurrentIndex(0)
        self.tabWidget_2.setCurrentIndex(0)
        self.tabWidget_4.setCurrentIndex(0)
        QtCore.QMetaObject.connectSlotsByName(MainWindow)

    def retranslateUi(self, MainWindow):
        _translate = QtCore.QCoreApplication.translate
        MainWindow.setWindowTitle(_translate("MainWindow", "SpinFlare"))
        self.label.setText(_translate("MainWindow", "Number of spin sites"))
        self.spin_combination.setItemText(0, _translate("MainWindow", "1/2"))
        self.spin_combination.setItemText(1, _translate("MainWindow", "1"))
        self.spin_combination.setItemText(2, _translate("MainWindow", "3/2"))
        self.spin_combination.setItemText(3, _translate("MainWindow", "2"))
        self.spin_combination.setItemText(4, _translate("MainWindow", "5/2"))
        self.nsites.setText(_translate("MainWindow", "6"))
        self.label_5.setText(_translate("MainWindow", "Table with the spins"))
        self.initialize_spins.setText(_translate("MainWindow", "Initialize spins"))
        self.tabWidget.setTabText(self.tabWidget.indexOf(self.tab_6), _translate("MainWindow", "Spins"))
        self.label_19.setText(_translate("MainWindow", "i     j    Jzz"))
        self.label_18.setText(_translate("MainWindow", "i     j    Jyy"))
        self.label_17.setText(_translate("MainWindow", "i     j    Jxx"))
        self.tabWidget_5.setTabText(self.tabWidget_5.indexOf(self.tab), _translate("MainWindow", "Exchange"))
        self.label_2.setText(_translate("MainWindow", "Couplings"))
        self.label_8.setText(_translate("MainWindow", "Dimerization"))
        self.Jxx.setText(_translate("MainWindow", "1.0"))
        self.Jyy.setText(_translate("MainWindow", "1.0"))
        self.Jzz.setText(_translate("MainWindow", "1.0"))
        self.label_3.setText(_translate("MainWindow", "Jyy"))
        self.label_9.setText(_translate("MainWindow", "Jxx"))
        self.label_21.setText(_translate("MainWindow", "Jzz"))
        self.dimerization_exchange.setText(_translate("MainWindow", "0.0"))
        self.tabWidget_5.setTabText(self.tabWidget_5.indexOf(self.tab_14), _translate("MainWindow", "Generators"))
        self.tabWidget.setTabText(self.tabWidget.indexOf(self.tab_5), _translate("MainWindow", "Exchange"))
        self.label_22.setText(_translate("MainWindow", "i     Bx"))
        self.label_23.setText(_translate("MainWindow", "i     By"))
        self.label_24.setText(_translate("MainWindow", "i     Bz"))
        self.tabWidget_3.setTabText(self.tabWidget_3.indexOf(self.tab_7), _translate("MainWindow", "Field table"))
        self.label_4.setText(_translate("MainWindow", "Modulation"))
        self.modulation_field.setItemText(0, _translate("MainWindow", "Ferromagnetic"))
        self.modulation_field.setItemText(1, _translate("MainWindow", "Antiferromagnetic"))
        self.label_25.setText(_translate("MainWindow", "Field strength"))
        self.label_26.setText(_translate("MainWindow", "Bx"))
        self.label_27.setText(_translate("MainWindow", "By"))
        self.label_28.setText(_translate("MainWindow", "Bz"))
        self.Bx.setText(_translate("MainWindow", "0.0"))
        self.By.setText(_translate("MainWindow", "0.0"))
        self.Bz.setText(_translate("MainWindow", "0.0"))
        self.label_29.setText(_translate("MainWindow", "Location"))
        self.location_field.setItemText(0, _translate("MainWindow", "Everywhere"))
        self.location_field.setItemText(1, _translate("MainWindow", "First"))
        self.location_field.setItemText(2, _translate("MainWindow", "Last"))
        self.tabWidget_3.setTabText(self.tabWidget_3.indexOf(self.tab_8), _translate("MainWindow", "Generator"))
        self.tabWidget.setTabText(self.tabWidget.indexOf(self.tab_9), _translate("MainWindow", "Fields"))
        self.maxm.setToolTip(_translate("MainWindow", "Bond dimension, this controls the accuracy of the MPS antsaz. Increase it for more accurate results"))
        self.maxm.setText(_translate("MainWindow", "10"))
        self.nsweeps.setToolTip(_translate("MainWindow", "Number of sweeps of the DMRG algorithm"))
        self.nsweeps.setText(_translate("MainWindow", "4"))
        self.label_10.setText(_translate("MainWindow", "Bond dimension"))
        self.label_11.setText(_translate("MainWindow", "Number of sweeps"))
        self.tabWidget.setTabText(self.tabWidget.indexOf(self.tab_10), _translate("MainWindow", "Params"))
        self.label_20.setText(_translate("MainWindow", "Component"))
        self.magnetization_type.setItemText(0, _translate("MainWindow", "X"))
        self.magnetization_type.setItemText(1, _translate("MainWindow", "Y"))
        self.magnetization_type.setItemText(2, _translate("MainWindow", "Z"))
        self.get_magnetization.setText(_translate("MainWindow", "Plot magnetization"))
        self.tabWidget_2.setTabText(self.tabWidget_2.indexOf(self.tab_2), _translate("MainWindow", "Ground state"))
        self.label_6.setText(_translate("MainWindow", "Sweep type"))
        self.static_arrangement.setToolTip(_translate("MainWindow", "Kind of sequence of static correlators to compute"))
        self.static_arrangement.setItemText(0, _translate("MainWindow", "From edge"))
        self.static_arrangement.setItemText(1, _translate("MainWindow", "From bulk"))
        self.static_arrangement.setItemText(2, _translate("MainWindow", "First neighbor right"))
        self.static_arrangement.setItemText(3, _translate("MainWindow", "First neighbor left"))
        self.label_7.setText(_translate("MainWindow", "Correlator"))
        self.static_type.setToolTip(_translate("MainWindow", "Type of the static correlator to compute"))
        self.static_type.setItemText(0, _translate("MainWindow", "SS"))
        self.static_type.setItemText(1, _translate("MainWindow", "XX"))
        self.static_type.setItemText(2, _translate("MainWindow", "YY"))
        self.static_type.setItemText(3, _translate("MainWindow", "ZZ"))
        self.get_static_correlator.setToolTip(_translate("MainWindow", "Compute and show the static correlator"))
        self.get_static_correlator.setText(_translate("MainWindow", "Show"))
        self.tabWidget_2.setTabText(self.tabWidget_2.indexOf(self.tab_3), _translate("MainWindow", "Static correlator"))
        self.dynamic_site_i_single.setToolTip(_translate("MainWindow", "Site on which to compute the correlator"))
        self.dynamic_site_i_single.setText(_translate("MainWindow", "0"))
        self.label_15.setText(_translate("MainWindow", "SIte"))
        self.get_dynamical_correlator_single.setToolTip(_translate("MainWindow", "Compute and show the dynamical correlator"))
        self.get_dynamical_correlator_single.setText(_translate("MainWindow", "Show correlator"))
        self.dynamic_type_single.setToolTip(_translate("MainWindow", "Type of correlator to compute"))
        self.dynamic_type_single.setItemText(0, _translate("MainWindow", "XX"))
        self.dynamic_type_single.setItemText(1, _translate("MainWindow", "YY"))
        self.dynamic_type_single.setItemText(2, _translate("MainWindow", "ZZ"))
        self.label_16.setText(_translate("MainWindow", "Correlator"))
        self.tabWidget_4.setTabText(self.tabWidget_4.indexOf(self.tab_11), _translate("MainWindow", "Single"))
        self.comboBox_2.setItemText(0, _translate("MainWindow", "Onsite"))
        self.comboBox_2.setItemText(1, _translate("MainWindow", "From edge"))
        self.comboBox_2.setItemText(2, _translate("MainWindow", "From bulk"))
        self.label_12.setText(_translate("MainWindow", "Type"))
        self.tabWidget_4.setTabText(self.tabWidget_4.indexOf(self.tab_12), _translate("MainWindow", "Map"))
        self.comboBox_6.setToolTip(_translate("MainWindow", "Method used to compute the dynamical correlator"))
        self.comboBox_6.setItemText(0, _translate("MainWindow", "KPM"))
        self.comboBox_6.setItemText(1, _translate("MainWindow", "TD"))
        self.smearing_dynamic.setToolTip(_translate("MainWindow", "Frequency smearing of the dynamical correlator"))
        self.smearing_dynamic.setText(_translate("MainWindow", "0.1"))
        self.label_13.setText(_translate("MainWindow", "Smearing"))
        self.label_14.setText(_translate("MainWindow", "Technique"))
        self.tabWidget_4.setTabText(self.tabWidget_4.indexOf(self.tab_13), _translate("MainWindow", "Params"))
        self.tabWidget_2.setTabText(self.tabWidget_2.indexOf(self.tab_4), _translate("MainWindow", "Dynamical correlator"))

