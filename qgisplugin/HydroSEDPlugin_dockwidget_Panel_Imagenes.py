# -*- coding: utf-8 -*-
"""
/***************************************************************************
 HydroSEDPlugin_dockwidget_Panel_Imagenes
                                 A QGIS plugin
 Este plugin facilita la ejecucion de modelos hidrologicos implementados desde QGIS.
                             -------------------
        begin                : 2018-05-02
        git sha              : $Format:%H$
        copyright            : (C) 2018 by Universidad Nacional de Colombia - Sede Medellín
        email                : 
 ***************************************************************************/

/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***
"""

import os

from PyQt4 import QtGui, uic
from PyQt4.QtCore import pyqtSignal, QUrl
from PyQt4.QtWebKit import QWebView;

import numpy as np
import os

from plotly.offline import plot
import plotly.graph_objs as go

FORM_CLASS, _ = uic.loadUiType(os.path.join(
    os.path.dirname(__file__), 'HydroSEDPlugin_dockwidget_Panel_Imagenes_base.ui'))
#FORM_CLASS_I, _ = uic.loadUiType(os.path.join(
#    os.path.dirname(__file__), 'HydroSEDPlugin_dockwidget_base.ui'))
#FORM_CLASS_I, _ = uic.loadUiType(os.path.join(
#    os.path.dirname(__file__), 'HydroSEDPlugin_dockwidget_base.ui'))

#print FORM_CLASS_I
#print FORM_CLASS_I.labelImagen1

class HydroSEDPluginDockWidgetPanelImagenes(QtGui.QDockWidget, FORM_CLASS):

    closingPlugin = pyqtSignal()

    def __init__(self, iface = None, parent=None):

        """Constructor."""
        super(HydroSEDPluginDockWidgetPanelImagenes, self).__init__(parent)
        self.setupUi(self)

        if not (iface is None):
            self.iface = iface

        self.loadImage (os.path.dirname (os.path.realpath (__file__)) + "/" + "imagenPrueba.jpg")

    def closeEvent(self, event):

        self.closingPlugin.emit()
        event.accept()

    def loadImage (self, pathImage):

        qtPixmap = QtGui.QPixmap (pathImage)

        wLabel = self.labelContenidoImagen1.width ();
        hLabel = self.labelContenidoImagen1.height ();

        self.labelContenidoImagen1.setPixmap (qtPixmap.scaled (wLabel, hLabel, 1));

        N = 1000
        random_x = np.random.randn(N)
        random_y = np.random.randn(N)

        # Create a trace
        trace = go.Scatter(
                             x = random_x,
                             y = random_y,
                             mode = 'markers'
                          )

        data = [trace]

        plot(data, filename='/tmp/testGraficaPPlotly.html', auto_open=False)

        self.myWV = QWebView(None); self.myWV.load(QUrl.fromLocalFile('/tmp/testGraficaPPlotly.html')); self.myWV.show()








