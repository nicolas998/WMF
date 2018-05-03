# -*- coding: utf-8 -*-
"""
/***************************************************************************
 HydroSEDPluginDockWidget
                                 A QGIS plugin
 Este plugin facilita la ejecucion de modelos hidrologicos implementados desde QGIS.
                             -------------------
        begin                : 2018-05-02
        git sha              : $Format:%H$
        copyright            : (C) 2018 by Universidad Nacional de Colombia - Sede Medellín
        email                : jctrujil@unal.edu.co
 ***************************************************************************/

/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/
"""

import os

from PyQt4 import QtGui, uic
from PyQt4.QtCore import pyqtSignal

from qgis.gui import QgsMessageBar
from PyQt4.QtGui import QFileDialog

from wmf import wmf

FORM_CLASS, _ = uic.loadUiType(os.path.join(
    os.path.dirname(__file__), 'HydroSEDPlugin_dockwidget_base.ui'))


class HydroSEDPluginDockWidget(QtGui.QDockWidget, FORM_CLASS):

    closingPlugin = pyqtSignal()

    def __init__(self, parent=None):

        """Constructor."""
        super(HydroSEDPluginDockWidget, self).__init__(parent)
        # Set up the user interface from Designer.
        # After setupUI you can access any designer object by doing
        # self.<objectname>, and you can use autoconnect slots - see
        # http://qt-project.org/doc/qt-4.8/designer-using-a-ui-file.html
        # #widgets-and-dialogs-with-auto-connect
        self.setupUi(self)
        self.setupUIInputsOutputs ()
        self.setupUIButtonEvents ()

    def closeEvent(self, event):

        self.closingPlugin.emit()
        event.accept()

    def handleClickEventEjecutarTrazadorCorrientes (self):

        print wmf
        print self.spinBoxLatitudTrazadorCorrientes.value ()
        print self.spinBoxLongitudTrazadorCorrientes.value ()
        print self.lineEditSelectorOutputCorrienteShapefileTrazadorCorrientes.text ()


    def handleClickEventEjecutarTrazadorCuencas (self):

        print wmf
        print self.spinBoxLatitudTrazadorCuencas.value ()
        print self.spinBoxLongitudTrazadorCuencas.value ()
        print self.lineEditInputCorrienteShapefileTrazadorCuencas.text ()
        print self.lineEditOutputCuencaShapefileTrazadorCuencas.text ()
        print self.lineEditOutputCuencaNCTrazadorCuencas.text ()


#        from PyQt4.QtGui import QFileDialog
#        filename1 = QFileDialog.getOpenFileName()
#        filename2 = QFileDialog.getSaveFileName()
#        print filename1
#        print filename2

    def setupUIInputsOutputs (self):

        def setupLineEditButtonOpenFileDialog (lineEditHolder, fileDialogHolder):

            lineEditHolder.setText (fileDialogHolder.getOpenFileName ())

        def setupLineEditButtonSaveFileDialog (lineEditHolder, fileDialogHolder):

            lineEditHolder.setText (fileDialogHolder.getSaveFileName ())

        def clickEventSelectorMapaDEM ():

            setupLineEditButtonOpenFileDialog (self.lineEditMapaDEM, QFileDialog)

        def clickEventSelectorMapaDIR ():

            setupLineEditButtonOpenFileDialog (self.lineEditMapaDIR, QFileDialog)

        def clickEventSelectorBinarioNC ():

            setupLineEditButtonOpenFileDialog (self.lineEditSelectorBinarioNC, QFileDialog)

        def clickEventSelectorOutputCorrienteShapefileTrazadorCorrientes ():

            setupLineEditButtonSaveFileDialog (self.lineEditSelectorOutputCorrienteShapefileTrazadorCorrientes, QFileDialog)

        def clickEventSelectorInputCorrienteShapefileTrazadorCuencas ():

            setupLineEditButtonOpenFileDialog (self.lineEditInputCorrienteShapefileTrazadorCuencas, QFileDialog)

        def clickEventSelectorOutputCuencaShapefileTrazadorCuencas ():

            setupLineEditButtonSaveFileDialog (self.lineEditOutputCuencaShapefileTrazadorCuencas, QFileDialog)

        def clickEventSelectorOutputCuencaNCTrazadorCuencas ():

            setupLineEditButtonSaveFileDialog (self.lineEditOutputCuencaNCTrazadorCuencas, QFileDialog)

        self.botonSelectorMapaDEM.clicked.connect (clickEventSelectorMapaDEM)
        self.botonSelectorMapaDIR.clicked.connect (clickEventSelectorMapaDIR)
        self.botonSelectorBinarioNC.clicked.connect (clickEventSelectorBinarioNC)
        self.botonSelectorOutputCorrienteShapefileTrazadorCorrientes.clicked.connect (clickEventSelectorOutputCorrienteShapefileTrazadorCorrientes)
        self.botonInputCorrienteShapefileTrazadorCuencas.clicked.connect (clickEventSelectorInputCorrienteShapefileTrazadorCuencas)
        self.botonOutputCuencaShapefileTrazadorCuencas.clicked.connect (clickEventSelectorOutputCuencaShapefileTrazadorCuencas)
        self.botonOutputCuencaNCTrazadorCuencas.clicked.connect (clickEventSelectorOutputCuencaNCTrazadorCuencas)

    def setupUIButtonEvents (self):

        self.botonEjecutarTrazadorCorrientes.clicked.connect (self.handleClickEventEjecutarTrazadorCorrientes)
        self.botonEjecutarTrazadorCuencas.clicked.connect (self.handleClickEventEjecutarTrazadorCuencas)





