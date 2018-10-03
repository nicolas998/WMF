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
 comando para cargar todos los iconos: pyrcc4 -o resources.py resources.qrc
"""

import os
import numpy as np
import glob 
import datetime as dt 

from PyQt4 import QtGui, uic, QtCore, Qt
from PyQt4.QtCore import pyqtSignal, QUrl
from PyQt4.QtWebKit import QWebView

from qgis.gui import QgsMessageBar
from PyQt4.QtGui import QFileDialog, QTableWidgetItem, QAbstractItemView

import os.path

import HydroSEDPluginUtils as HSutils
import HydroGetCoordinates as HSCoord
import HydroSEDPlots as HSplots

import GdalTools_utils as GdalTools_utils

FORM_CLASS, _ = uic.loadUiType(os.path.join(
    os.path.dirname(__file__), 'HydroSEDPlugin_dockwidget_base.ui'))


class HydroSEDPluginDockWidget(QtGui.QDockWidget, FORM_CLASS):

    closingPlugin = pyqtSignal()

    def __init__(self, iface = None, parent=None):

        """Constructor."""
        super(HydroSEDPluginDockWidget, self).__init__(parent)
        # Set up the user interface from Designer.
        # After setupUI you can access any designer object by doing
        # self.<objectname>, and you can use autoconnect slots - see
        # http://qt-project.org/doc/qt-4.8/designer-using-a-ui-file.html
        # #widgets-and-dialogs-with-auto-connect
        self.setupUi(self)
        self.setupUIInputsOutputs ()
        self.setupHidro_Balance()
        self.setupBasinManager()
        self.setupGeomorfologia()
        #self.setupUIButtonEvents ()
        self.setupRainfall()
        self.setupSimulation()
        self.setupNcVariables()
        
        self.TablaFila_WMF = 0
        self.TablaFila_NC = 0
        self.botonClicked_NC = False
        self.botonClicked_WMF = False
        self.segundaCarga = False
        self.Segunda_carga_Qobs = False 
        self.Segunda_carga_Qobs_Sed = False 
        self.Path2Radar = ''
        
        self.GeoPlots = HSplots.PlotGeomorphology()
        self.SimStoragePlots = HSplots.PlotStorage()
        
        #Inicia el comboBox de la seleccion de categoria para transformar raster a WMF
        for k in ['base','Geomorfo','SimHidro','Hidro']:
            self.ComboBoxRaster2WMF.addItem(k)        
        self.setupRaster2WMF()
        
        if not (iface is None):
            self.iface = iface
        self.HSutils = HSutils.controlHS()        
        self.GetCoordsCorriente = HSCoord.PointTool(self.iface.mapCanvas(), self.spinBoxLatitudTrazadorCorrientes,
            self.spinBoxLongitudTrazadorCorrientes)
        self.GetCoordsCuenca = HSCoord.PointTool(self.iface.mapCanvas(), self.spinBoxLatitudTrazadorCuencas,
            self.spinBoxLongitudTrazadorCuencas)
        #self.iface.mapCanvas().setMapTool(GetCoords)

    def closeEvent(self, event):

        self.closingPlugin.emit()
        event.accept()

    def handleClickCoordCorrientes(self):
        self.iface.mapCanvas().setMapTool(self.GetCoordsCorriente)

    def handleClickCoordCuencas(self):
        self.iface.mapCanvas().setMapTool(self.GetCoordsCuenca)


    def handleClickEventEjecutarTrazadorCorrientes (self):
        #Traza la corriente
        y = self.spinBoxLatitudTrazadorCorrientes.value ()
        x = self.spinBoxLongitudTrazadorCorrientes.value ()
        #Camino a los temporales
        if len(self.PathCorriente_out.text()) == 0:
            self.PathCorriente_out.setText('/tmp/HydroSED/vector/Corriente.shp')
        OutPath = self.PathCorriente_out.text()
        try:
            self.HSutils.trazador_corriente(x,y, OutPath)
            ret, layer = self.HSutils.cargar_mapa_vector(OutPath, self.HSutils.TIPO_STYLE_POLILINEA)
            
            self.iface.mapCanvas().refresh() 
            self.iface.legendInterface().refreshLayerSymbology(layer)
            
            self.iface.messageBar().pushInfo(u'HidroSIG',u'Se ha trazado la corriente de forma exitosa')
        except:
            pass

    def handleClickEventEjecutarTrazadorCuencas (self):
        #Obtiene coordenadas
        y = self.spinBoxLatitudTrazadorCuencas.value ()
        x = self.spinBoxLongitudTrazadorCuencas.value ()
        #Camino a los temporales
        if len(self.PathOutputDivisoria.text()) == 0:
            self.PathOutputDivisoria.setText('/tmp/HydroSED/vector/Cuenca.shp')
        if len(self.PathOutputRed.text()) == 0:
            self.PathOutputRed.setText('/tmp/HydroSED/vector/Cuenca_Red.shp')
        #Paths para guardar la cuenca 
        OutPathDivisoria = self.PathOutputDivisoria.text()
        OutPathRed = self.PathOutputRed.text()
        OutPathNC = self.PathOutputNETCDF.text()
        #Traza la cuenca
        self.HSutils.trazador_cuenca(x,y, self.spinBox_dxPlano.value(),
            self.spinBoxUmbralRed.value(),OutPathDivisoria, OutPathRed, OutPathNC, self.lineEditMapaDEM,
            self.lineEditMapaDIR)
        #Establece a la cuenca trazada como el proyecto actual.
        if len(OutPathNC)>2:
            self.lineEditRutaCuenca.setText(OutPathNC)
            self.ButtonLoadBasinProyect.setEnabled(True)
                
        #Carga la divisoria
        ret, layer = self.HSutils.cargar_mapa_vector(OutPathDivisoria, self.HSutils.TIPO_STYLE_POLIGONO, color = (255,0,0), width = 0.6)            
        self.iface.mapCanvas().refresh() 
        self.iface.legendInterface().refreshLayerSymbology(layer)            

        #Carga la red 
        ret, layer = self.HSutils.cargar_mapa_vector(OutPathRed, self.HSutils.TIPO_STYLE_POLILINEA, width = 0.4)
        self.iface.mapCanvas().refresh() 
        self.iface.legendInterface().refreshLayerSymbology(layer)         

        #mensaje de exito
        self.iface.messageBar().pushInfo(u'HidroSIG',u'Se ha trazado la cuenca de forma exitosa')        

    def setupBasinManager(self):
        '''Conjunto de herramientas y variables que permiten cargar y actualizar archivos .nc
        de cuencas'''
        
        def setupLineEditButtonOpenBasinFileDialog (lineEditHolder, fileDialogHolder):
            '''Busca un proyecto de cuenca ya guardado anteriormente'''
            lineEditHolder.setText (fileDialogHolder.getOpenFileName (QtGui.QDialog (), "Cargador de cuencas", "*","Cuenca en NetCDF4(*.nc);;", 
                QtGui.QFileDialog.DontUseNativeDialog))
        
        def clickEventSelectShape2SaveNc():
            self.ruta2network = SelectNetworkshapeDialog (QFileDialog)
            self.Tabla_Prop_NC.setSelectionMode(QAbstractItemView.SingleSelection)
            self.ruta2network
        
        def SelectNetworkshapeDialog (fileDialogHolder):
            '''Hace que cuando se busquen shapes solo se encuetren formatos vectoriales'''            
            lineEditHolder = fileDialogHolder.getSaveFileName (QtGui.QDialog (), "Guardar en capa vectorial...", "*", "Shapefiles (*.shp);;")
            return lineEditHolder
            #if ((os.path.exists (lineEditHolder.text ().strip ())) and (not (self.iface is None))):
             #   self.iface.addVectorLayer (lineEditHolder.text ().strip (), os.path.basename (lineEditHolder.text ()).strip (), "ogr")
        
        #Funciones para Cargar variables
        def clickEventSelectorBasin():
            '''click para seleccionar un proyecto de cuenca'''
            #Pone el texto de la ruta 
            setupLineEditButtonOpenBasinFileDialog(self.lineEditRutaCuenca, QFileDialog)
            #Habilita visualizarlo 
            if len(self.lineEditRutaCuenca.text())>2:
                self.ButtonLoadBasinProyect.setEnabled(True)
        def clickEventBasin2WMF():
            '''Agrega el proyecto de cuenca a WMF'''
            #Checks para simulacion de sedimentos e hidrologica
            Simhidro = self.checkBox_simBasin.isChecked()
            SimSed = self.checkBox_simSed.isChecked()
            #Si hay una tabla vieja le trata de borrar todas sus entradas
            if self.segundaCarga:
                self.TabNC.EmptyTable(self.Tabla_Prop_NC)
                self.TabWMF.EmptyTable(self.Tabla_Prop_WMF)
            self.segundaCarga = True
            #Cargado de la cuenca
            Area, self.EPSG, dxp, self.noData, self.umbral = self.HSutils.Basin_LoadBasin(self.lineEditRutaCuenca.text().strip(),
                Simhidro, SimSed)
            #self.HSutils.Basin_LoadBasin(self.lineEditRutaCuenca.text().strip(), Simhidro, SimSed)
            self.TableStart()
            #Actualiza tabla de Nc y comboBox 
            self.VarFromNC.clear()
            for k in self.HSutils.DicBasinNc.keys():
                self.TabNC.NewEntry(self.HSutils.DicBasinNc[k],k, self.Tabla_Prop_NC)
                self.VarFromNC.addItem(k)
            #Area, self.EPSG, dxp, self.noData, self.umbral = self.HSutils.Basin_LoadBasin(self.lineEditRutaCuenca.text().strip(),
        #       Simhidro, SimSed)
            #print self.umbral
            #Habilita los botones de visualizacion de red hidrica y divisoria 
            self.Boton_verDivisoria.setEnabled(True)
            self.Boton_verRedHidrica.setEnabled(True)
            #Pone el area de la cuenca 
            texto = '%.1f'%Area
            self.LabelBasinArea.setText(texto)
            self.LineEditEPSG.setText(self.EPSG)
            texto = '%.1f' % self.noData
            self.LineEditNoData.setText(texto)
            self.spinBox_dxPlano.setValue(dxp)
            self.SpinGeoUmbralCanal.setValue(self.umbral)
            #Actualiza comboBox de goemorfo
            for k in self.HSutils.DicBasinWMF.keys():
                self.ComboGeoMaskVar.addItem(k)
            #Actualiza Param de claibracion 
            if Simhidro:
                self.ParamNamesCombo.clear()
                for k in self.HSutils.DicParameters.keys():
                    self.ParamNamesCombo.addItem(k)
            self.iface.messageBar().pushMessage (u'Hydro-SIG:', u'El proyecto de cuenca ha sido cargado con exito',
                level=QgsMessageBar.INFO, duration=3)
            
        def clickEventBasinLoadDivisory():
            '''Carga la divisoria de la cuenca cargada a WMF'''
            OutPathDivisoria = '/tmp/HydroSED/vector/CuencaCargada.shp'
            self.HSutils.Basin_LoadBasinDivisory(OutPathDivisoria)
            #Carga la divisoria
            ret, layer = self.HSutils.cargar_mapa_vector(OutPathDivisoria, self.HSutils.TIPO_STYLE_POLIGONO, color = (255,0,0), width = 0.6)            
            self.iface.mapCanvas().refresh() 
            self.iface.legendInterface().refreshLayerSymbology(layer)      
        def clickEventBasinLoadNetwork():
            '''Carga la divisoria de la cuenca cargada a WMF'''
            OutPathNetwork = '/tmp/HydroSED/vector/RedCargada.shp'
            self.HSutils.Basin_LoadBasinNetwork(OutPathNetwork)
            #Carga la red 
            ret, layer = self.HSutils.cargar_mapa_vector(OutPathNetwork, self.HSutils.TIPO_STYLE_POLILINEA, width = 0.4)
            self.iface.mapCanvas().refresh() 
            self.iface.legendInterface().refreshLayerSymbology(layer)   
        
        def clickEventBasinVarNC2Network():
            '''Convierte un conjunto de variables a una red hidrica de la cuenca'''
            if self.botonClicked_NC is False:
                #Hace que la seleccion sea multiple en la columna de NC
                self.Tabla_Prop_NC.setSelectionMode(QAbstractItemView.MultiSelection)
                self.ruta2Network = ''
                #imprime las variables seleccionadas
                self.Vars2Network = []
                self.botonClicked_NC = True
                return
            elif self.botonClicked_NC is True:
                for i in self.Tabla_Prop_NC.selectedItems()[::4]:
                    self.Vars2Network.append(i.text())
                #En el segundo click selecciona el archivo y ejecuta 
                self.ruta2Network = SelectNetworkshapeDialog(QFileDialog)
                self.Tabla_Prop_NC.setSelectionMode(QAbstractItemView.SingleSelection)
                self.botonClicked_NC = False
                #Guardado de las variables en formato red hidrica
                self.HSutils.BasinNc2Network(self.ruta2Network,self.Vars2Network)
                self.iface.messageBar().pushMessage (u'Hydro-SIG:', u'Las variables seleccionadas han sido convertidas a la red hidrica',
                    level=QgsMessageBar.INFO, duration=3)
        
        def clickEventBasinVarWMF2Network():
            '''Convierte un conjunto de variables a una red hidrica de la cuenca'''
            if self.botonClicked_WMF is False:
                #Hace que la seleccion sea multiple en la columna de NC
                self.Tabla_Prop_WMF.setSelectionMode(QAbstractItemView.MultiSelection)
                self.ruta2Network = ''
                #imprime las variables seleccionadas
                self.Vars2Network = []
                self.botonClicked_WMF = True
                return
            elif self.botonClicked_WMF is True:
                for i in self.Tabla_Prop_WMF.selectedItems()[::4]:
                    self.Vars2Network.append(i.text())
                #En el segundo click selecciona el archivo y ejecuta 
                self.ruta2Network = SelectNetworkshapeDialog(QFileDialog)
                self.Tabla_Prop_WMF.setSelectionMode(QAbstractItemView.SingleSelection)
                self.botonClicked_WMF = False
                #Guardado de las variables en formato red hidrica
                self.HSutils.BasinWMF2Network(self.ruta2Network,self.Vars2Network)
                self.iface.messageBar().pushMessage (u'Hydro-SIG:', u'Las variables seleccionadas han sido convertidas a la red hidrica',
                    level=QgsMessageBar.INFO, duration=3)
        
        #Botones para variables de entrada 
        self.botonSelectorProyectBasin.clicked.connect(clickEventSelectorBasin)
        self.ButtonLoadBasinProyect.clicked.connect(clickEventBasin2WMF)
        #Botones para visualizar polilineas y poligonos 
        self.Boton_verDivisoria.clicked.connect(clickEventBasinLoadDivisory)
        self.Boton_verRedHidrica.clicked.connect(clickEventBasinLoadNetwork)
        #boton para convertir variables a red hidrica
        self.ButtonVar2Net_NC.clicked.connect(clickEventBasinVarNC2Network)
        self.ButtonVar2Net_WMF.clicked.connect(clickEventBasinVarWMF2Network)
        

    def setupGeomorfologia(self):
        '''Conjunto de herramientas para manejar parametros geomorfologicos de la cuenca analizada'''
        
        

        def clickEventActivateGeoCheckBoxes():
            '''Selecciona y des-selecciona todas las opciones de calculo de una ves'''
            #Si esta checked todas se seleccionan
            if self.checkBoxTodos.isChecked():
                self.checkBoxArea.setChecked(True)
                self.checkBoxOrder.setChecked(True)
                self.checkBoxChannels.setChecked(True)
                self.checkBoxDist2Out.setChecked(True)
                self.checkBoxHAND.setChecked(True)
                self.checkBoxIT.setChecked(True)
                self.checkBoxOCG.setChecked(True)
                self.checkBoxKubota.setChecked(True)
                self.checkBoxRunoff.setChecked(True)
            #Si no, des-selecciona a todas
            else:
                self.checkBoxArea.setChecked(False)
                self.checkBoxOrder.setChecked(False)
                self.checkBoxChannels.setChecked(False)
                self.checkBoxDist2Out.setChecked(False)
                self.checkBoxHAND.setChecked(False)
                self.checkBoxIT.setChecked(False)
                self.checkBoxOCG.setChecked(False)
                self.checkBoxKubota.setChecked(False)
                self.checkBoxRunoff.setChecked(False)
            
        def GeoTableStart():
            '''Inicia la tabla donde monta los parametros geomorfologicos de la cuenca'''
            self.GeoTableNumRows = 23
            self.GeoTableNumItems = 0
            self.GeoTableHeader = ["Parametro", "Valor", "Unidades"]
            self.TableGeoParameters.setRowCount(self.GeoTableNumRows)
            self.TableGeoParameters.setColumnCount(len(self.GeoTableHeader))
            self.TableGeoParameters.setHorizontalHeaderLabels(self.GeoTableHeader)
        
        def clickEventGeoProperties():
            '''Calcula los parametros geomorfologicos de la cuenca.'''
            #Reinicia la talba para que no se llene de cosas
            self.GeoTableNumItems = 0
            self.TableGeoParameters.clear()
            self.TableGeoParameters.clearContents()
            #Calcula los parametros.
            Param,Tc = self.HSutils.Basin_GeoGetParameters()
            print Param
            #Inicia la tabla 
            GeoTableStart()
            #le pone los parametros
            for d in Param.keys():
                #Obtiene nombre y unidad
                cadena = d.split('_')
                unidad = cadena[-1].split('[')[-1].split(']')[0]
                nombre = ' '.join(cadena[:-1])
                valor = '%.3f' % Param[d]
                #Actualiza la tabla
                self.TableGeoParameters.setItem (self.GeoTableNumItems, 0, QTableWidgetItem(nombre))
                self.TableGeoParameters.setItem (self.GeoTableNumItems, 1, QTableWidgetItem(valor))
                self.TableGeoParameters.setItem (self.GeoTableNumItems, 2, QTableWidgetItem(unidad))
                self.GeoTableNumItems += 1
            for d in Tc.keys():
                nombre = d
                valor = '%.3f' % Tc[d]
                unidad = 'Hrs'
                #Actualiza la tabla
                self.TableGeoParameters.setItem (self.GeoTableNumItems, 0, QTableWidgetItem(nombre))
                self.TableGeoParameters.setItem (self.GeoTableNumItems, 1, QTableWidgetItem(valor))
                self.TableGeoParameters.setItem (self.GeoTableNumItems, 2, QTableWidgetItem(unidad))
                self.GeoTableNumItems += 1
            #Habilita los botones de geomrofologia 
            self.ButtonGeoParameters2Excel.setEnabled(True)
            self.Button_GeomorfoPerfil.setEnabled(True)
            self.Button_GeomorfoTc.setEnabled(True)
        
        def setupLineEditButtonSaveFileDialog (fileDialogHolder):
            '''Pone la ruta elegida en el dialogo de texto para guardado'''
            return fileDialogHolder.getSaveFileName (QtGui.QDialog (), "Guardar parametros en un archivo de excel", "*", "Excel (*.xlsx);;")
        
        def clickEventExportParam2Excel():
            '''Exporta los parametros geomorfologicos a excel'''
            #Set de la ruta donde se guarda el archivo 
            PathExcel = setupLineEditButtonSaveFileDialog(QFileDialog)
            #Funcion que pone param en excel 
            Retorno = self.HSutils.Basin_Geo2Excel(PathExcel)
            #Mensage de exito o error
            if Retorno == 0:
                self.iface.messageBar().pushInfo(u'HidroSIG',u'Los parametros geomorfologicos se han exportado correctamente')
            else:
                self.iface.messageBar().pushMessage (u'Hydro-SIG:', u'No ha sido posible exportar los parametros geomorfologicos de la cuenca',
                    level=QgsMessageBar.WARNING, duration=5)
                
        def clickEventGeoRasterProp():
            '''calcula los parametros geomorfologicos de la cuenca por raster'''
            #Lista de variables a calcular
            ListaVar = []
            #Revisa cada checkbox
            if self.checkBoxArea.isChecked():
                self.HSutils.Basin_GeoGetAcumSlope()
                ListaVar.extend(['Area','Pendiente'])
            if self.checkBoxOrder.isChecked():
                self.HSutils.Basin_GeoGetOrder()
                ListaVar.extend(['Order_hills','Order_channels'])
            if self.checkBoxHAND.isChecked():
                self.HSutils.Basin_GeoGetHAND()
                ListaVar.extend(['HAND','HDND','HAND_class'])
            if self.checkBoxIT.isChecked():
                self.HSutils.Basin_GeoGetIT()
                ListaVar.extend(['Topo_index'])
            if self.checkBoxDist2Out.isChecked():
                self.HSutils.Basin_GeoGetDist2Out()
                ListaVar.extend(['Dist2Out'])
            if self.checkBoxChannels.isChecked():
                self.HSutils.Basin_GeoGetChannels()
                ListaVar.extend(['Channels'])
            if self.checkBoxOCG.isChecked():
                self.HSutils.Basin_GeoGetOCG()
                ListaVar.extend(['OCG_coef','OCG_exp'])
            if self.checkBoxKubota.isChecked():
                self.HSutils.Basin_GeoGetKubota()                
                ListaVar.extend(['Kubota_coef','Kubota_exp'])
                #mensajes de advertencia por si no han sido cargadas las variables. 
                if 'h1_max' not in self.HSutils.DicBasinNc.keys(): 
                    self.iface.messageBar().pushMessage (u'Hydro-SIG:',
                    'Como no se ha cargado h1_max al NC, se estima con h1_max = 100 mm ',
                    level=QgsMessageBar.WARNING, duration=3)    
                if 'v_coef' not in self.HSutils.DicBasinNc.keys(): 
                    self.iface.messageBar().pushMessage (u'Hydro-SIG:',
                    'Como no se ha cargado Ks al NC, se estima con Ks = 0.003 mm/s ',
                    level=QgsMessageBar.WARNING, duration=3)    
            if self.checkBoxRunoff.isChecked():
                E1 = self.RunoffE1.value()
                Epsilon1 = self.RunoffEpsi.value()
                self.HSutils.Basin_GeoGetRunoff(e1=E1,Epsilon=Epsilon1)
                ListaVar.extend(['Runoff_coef'])
                ListaVar.extend(['Runoff_exp'])
                #Mensaje de advertencia por si no ha sido cargada la variable. 
                if 'Manning' not in self.HSutils.DicBasinWMF.keys() and 'Manning' not in self.HSutils.DicBasinNc.keys():
                    self.iface.messageBar().pushMessage (u'Hydro-SIG:',
                    'Como no se ha cargado Manning al NC, se estima con n_manning = 0.05',
                    level=QgsMessageBar.WARNING, duration=3)
            #mensaje de caso de exito
            self.iface.messageBar().pushMessage(u'HidroSIG:',
            u'Calculo de geomorfologia distribuida realizado, revisar la tabla Variables WMF.',
            level=QgsMessageBar.INFO, duration=2)
            #Actualiza la tabla de variables temporales 
            for k in ListaVar:
                self.TabWMF.NewEntry(self.HSutils.DicBasinWMF[k],k, self.Tabla_Prop_WMF)
                #self.ComboGeoMaskVar.addItem(k)
                #self.ComboGeoVar2Acum.addItem(k)
        
        def PlotTiempoConcentracion():
            '''Plot con el tiempo de concentracion estimado por diferentes metodologias'''
            #Obtiene la variable e invoca la funciond e grafica
            PathFigure = '/tmp/HydroSED/Plots_Geomorfo/TiemposConcentracion.html'
            self.GeoPlots.TiempoConcentracion(self.HSutils.cuenca.Tc, PathFigure)
            #Set de la ventana que contiene la figura.
            self.VistaRainWeb = QWebView(None)
            self.VistaRainWeb.load(QUrl.fromLocalFile(PathFigure))
            self.VistaRainWeb.setWindowTitle('Tiempos de concentracion')
            self.VistaRainWeb.setMinimumWidth(100)
            self.VistaRainWeb.setMaximumWidth(600)
            self.VistaRainWeb.setMinimumHeight(100)
            self.VistaRainWeb.setMaximumHeight(500)
            self.VistaRainWeb.show()
         
        def PlotCurvaHipsometricaPerfil():
            '''Hace el plot de la curva hipsometrica de la cuenca y del perfil del canal ppal'''
            #Obtiene la variable e invoca la funciond e grafica
            PathFigure = '/tmp/HydroSED/Plots_Geomorfo/PerfilYCurvaHipso.html'
            self.HSutils.cuenca.GetGeo_Ppal_Hipsometric()
            #Variables
            Hipso = np.copy(self.HSutils.cuenca.hipso_basin)
            Stream = np.vstack([self.HSutils.cuenca.ppal_stream[1]/1000.,
                self.HSutils.cuenca.ppal_stream[0]])
            #Plot
            self.GeoPlots.StreamProfileHipsometric(PathFigure, Hipso, Stream)
            #Set de la ventana que contiene la figura.
            self.VistaRainWeb = QWebView(None)
            self.VistaRainWeb.load(QUrl.fromLocalFile(PathFigure))
            self.VistaRainWeb.setWindowTitle('Perfil longitudinal y curva hipsométrica')
            self.VistaRainWeb.setMinimumWidth(100)
            self.VistaRainWeb.setMaximumWidth(700)
            self.VistaRainWeb.setMinimumHeight(100)
            self.VistaRainWeb.setMaximumHeight(400)
            self.VistaRainWeb.show()
         
        def clickEventUpdateAcumVarMask():
            '''Actuaqliza la lista de variables existentes en las tablas para usar como mascara'''
            self.ComboGeoMaskVar.clear()
            for i in self.HSutils.DicBasinWMF.keys():
                self.ComboGeoMaskVar.addItem(i+'-WMF')
            for i in self.HSutils.DicBasinNc.keys():
                self.ComboGeoMaskVar.addItem(i+'-NC')
            self.ComboGeoMaskVar.addItem('sin mascara')
            self.ComboGeoMaskVar.setCurrentIndex(self.ComboGeoMaskVar.findText('sin mascara'))
        
        def clickEventUpdateAcumVar():
            '''Actuaqliza la lista de variables existentes en las tablas para acumular'''
            self.ComboGeoVar2Acum.clear()
            for i in self.HSutils.DicBasinWMF.keys():
                self.ComboGeoVar2Acum.addItem(i+'-WMF')
            for i in self.HSutils.DicBasinNc.keys():
                self.ComboGeoVar2Acum.addItem(i+'-NC')        
        
        def clickEventAcumVar():
            '''Acumula la variable seleccionada.'''
            #Toma el nombre de la variable 
            Nombre = self.NameGeoAcumVar.text().strip()
            umbral = self.UmbralAcum.value()
            #Variable mascara
            MascaraText = self.ComboGeoMaskVar.currentText()
            if MascaraText == 'sin mascara':
                VarMaskName = None
                WhatMaskDic = None
            else:
                WhatMaskDic = MascaraText.split('-')[-1]
                VarMaskName = MascaraText.split('-')[0]
            #Variable a acumular
            VarAcumText = self.ComboGeoVar2Acum.currentText()
            WhatAcumDic = VarAcumText.split('-')[-1]
            VarAcumName = VarAcumText.split('-')[0]
            #Invoca funcion para acumular variable 
            Retorno = self.HSutils.Basin_GeoAcumVar(Nombre, 
                VarAcumName, 
                WhatAcumDic,
                VarMaskName,
                WhatMaskDic,
                umbral)
            #Actualiza la tabla de WMF 
            self.TabWMF.NewEntry(self.HSutils.DicBasinWMF[Nombre],Nombre, self.Tabla_Prop_WMF)
            #Mensaje de exito o error 
            if Retorno == 0:
                self.iface.messageBar().pushMessage(u'HidroSIG:',u'La variable '+VarAcumName+' se ha acumulado como '+Nombre+' en la tabla WMF.',
                    level=QgsMessageBar.INFO,
                    duration = 5)
            else:
                self.iface.messageBar().pushMessage (u'Hydro-SIG:', 
                    u'No se ha logrado acumular la variable '+ VarAcumName,
                    level=QgsMessageBar.WARNING, duration=3)
            
        
        #Botones de ejecucion
        self.ButtonGeomorfoRasterVars.clicked.connect(clickEventGeoRasterProp)
        self.checkBoxTodos.clicked.connect(clickEventActivateGeoCheckBoxes)
        self.ButtonGeoParameters.clicked.connect(clickEventGeoProperties)
        self.ButtonUpdateMask.clicked.connect(clickEventUpdateAcumVarMask)
        self.ButtonUpdateVarList.clicked.connect(clickEventUpdateAcumVar)
        self.ButtonGeomorfoAcumVar.clicked.connect(clickEventAcumVar)
        #self.ComboGeoMaskVar.activated.connect(clickEventUpdateComboBoxMask)
        #Botones de figuras
        self.Button_GeomorfoTc.clicked.connect(PlotTiempoConcentracion)
        self.Button_GeomorfoPerfil.clicked.connect(PlotCurvaHipsometricaPerfil)
        #Boton para exporar datos a excel 
        self.ButtonGeoParameters2Excel.clicked.connect(clickEventExportParam2Excel)
        
        
    
    def setupHidro_Balance(self):
        
        #Inicia Combo box de funciones de distribucion
        self.ComboQmaxPDF.addItem('gumbel')
        self.ComboQmaxPDF.addItem('lognorm')
        self.ComboQmaxPDF.setCurrentIndex(self.ComboQmaxPDF.findData('Gumbel'))
        
        
        def setupLineEditButtonOpenRasterFileDialog (lineEditHolder, fileDialogHolder):
            '''Hace que solo se busquen formatos aceptados por GDAL'''
            lineEditHolder.setText (fileDialogHolder.getOpenFileName (QtGui.QDialog (), 'Open File',"", GdalTools_utils.FileFilter.allRastersFilter (),
                QtGui.QFileDialog.DontUseNativeDialog))
        def setupLineEditButtonSaveFileDialog (lineEditHolder, fileDialogHolder):
            '''Pone la ruta elegida en el dialogo de texto para guardado'''
            lineEditHolder.setText (fileDialogHolder.getSaveFileName ())
        
        #Funciones para Cargar variables
        def clickEventSelectorRaster():
            '''click para seleccionar el raster de lluvia'''
            #Pone el texto de la ruta 
            setupLineEditButtonOpenRasterFileDialog(self.PathInHydro_Rain, QFileDialog)
            #Habilita visualizarlo 
            if len(self.PathInHydro_Rain.text())>2:
                self.Button_HidroViewRain.setEnabled(True)
        
        #Funciones para visualizar variables 
        def clickEventViewRainfall():
            pathMapaRain = self.PathInHydro_Rain.text().strip()
            flagCargaMapaRain = self.HSutils.cargar_mapa_raster(pathMapaRain)
        def clickEventViewQmedNetwork():
            OutPathRed = self.PathOutHydro_Qmed.text().strip()
            ret, layer = self.HSutils.cargar_mapa_vector(OutPathRed, self.HSutils.TIPO_STYLE_POLILINEA, width = 0.4)
            self.iface.mapCanvas().refresh() 
            self.iface.legendInterface().refreshLayerSymbology(layer) 
        
        #Funciones para decir donde se van a guardar las variables.
        def clickEventOutQmed():
            '''Habilita la ruta para decir donde se va a guardar la variable de caudal medio'''
            setupLineEditButtonSaveFileDialog(self.PathOutHydro_Qmed, QFileDialog)
            
        #Funciones para el calculo de balance    
        def hadleClickEventEjecutarBalance():
            '''Hace el balance hidrologico una ves que se da click en el boton: Butto_Ejec_HidroBalance'''
            #Selecciona el tipo de etr
            if self.RadioBalance_ETR_Choudry.isChecked():
                TipoETR = 3
            if self.RadioBalance_ETR_Cenicafe.isChecked():
                TipoETR = 2
            if self.RadioBalance_ETR_Turc.isChecked():
                TipoETR = 1
            #Invoca la funcion
            Retorno, QSalida = self.HSutils.hidologia_balance(self.spinBox_dxPlano.value(),
                self.spinBoxUmbralRed.value(), 
                self.PathInHydro_Rain.text(), 
                TipoETR)
            #Pone el valor de cadual medio en el cuadro
            textoCaudal = '%.3f' % QSalida
            self.ShowResultQmed.setText(textoCaudal)
            #Mensaje de exito 
            if Retorno == 0:
                self.iface.messageBar().pushMessage(u'HidroSIG:',
                    u'Balance realizado: variables Caudal, ETR y Runoff cargadas a Tabla de propiedades WMF.',
                    level=QgsMessageBar.INFO, duration=3)
            else:
                self.iface.messageBar().pushMessage (u'Hydro-SIG:', 
                    u'No se ha logrado realizar el balance hidrológico en la cuenca',
                    level=QgsMessageBar.WARNING, duration=3)
            #Actualiza la tabla de variables temporales y actualiza comboBox de gomorfo
            for k in ['Caudal','ETR','Runoff']:
                self.TabWMF.NewEntry(self.HSutils.DicBasinWMF[k],k, self.Tabla_Prop_WMF)
                
        def handleClickEventExtremeByRegionalization():
            '''Calcula caudales extremos a partir de regionalizacion'''
            #Obtiene el caudal medio 
            NombreQmed = self.PathInHydro_Qmed4Qmax.text().strip()
            try:
                Qmed = np.copy(self.HSutils.DicBasinWMF[NombreQmed]['var'])
            except:
                try:
                    Qmed = np.copy(self.HSutils.DicBasinNc[NombreQmed]['var'])
                except:
                    self.iface.messageBar().pushMessage (u'Hydro-SIG:', 
                        u'No se encuentra la variable '+NombreQmed+' en WMF o NC, no es posible estimar Qmax o Qmin',
                        level=QgsMessageBar.ERROR, duration=3)
                    return
            #Si es maximo establece valores por defecto de coeifciente y exponentes
            if self.RadButton_maximos.isChecked():
                CoefExpoList = [6.71, 0.82, 3.29, 0.64]
                StrVal = ['%.4f' % i for i in CoefExpoList]
                isMaxorMin = 'QMax'
            elif self.RadButton_minimos.isChecked():
                CoefExpoList = [0.4168, 1.058, 0.2, 0.98]
                StrVal = ['%.4f' % i for i in CoefExpoList]
                isMaxorMin = 'QMin'
            #Valores por defecto
            for i in range(1,5):
                val = getattr(self, 'SpinBox_C'+str(i)).value()
                if val == 0.0:
                    comando = 'self.SpinBox_C'+str(i)+'.setValue('+StrVal[i-1]+')'
                    eval(comando)
            #Obtiene la funcion de distribucion 
            Pdf2Use = self.ComboQmaxPDF.currentText().encode()
            #Ejecuta el calculo de caudal a largo plazo 
            retorno = self.HSutils.hidrologia_extremos_regional(Qmed, CoefExpoList, Pdf2Use, isMaxorMin)
            #Mensajes de exito o error
            if retorno == 0:
                #Actualiza la tabla de WMF
                for Tr in [2.33, 5, 10, 25, 50, 100]:
                    nombre = isMaxorMin+'_'+str(Tr)
                    self.TabWMF.NewEntry(self.HSutils.DicBasinWMF[nombre],nombre, self.Tabla_Prop_WMF)
                #Mensaje de exito 
                self.iface.messageBar().pushMessage (u'Hydro-SIG:', 
                    u'Caudales extremos calculados con exito',
                    level=QgsMessageBar.INFO, duration=3)
            else:
                self.iface.messageBar().pushMessage (u'Hydro-SIG:', 
                    u'No ha sido posible determinar caudales extremos para esta cuenca',
                    level=QgsMessageBar.WARNING, duration=3)
                    
            
        
        #Botones para ejecutar
        self.Butto_Ejec_HidroBalance.clicked.connect(hadleClickEventEjecutarBalance)
        self.Butto_Ejec_HidroExtremos.clicked.connect(handleClickEventExtremeByRegionalization)
        
    def TableStart (self):
        '''Arranca las tablas de NC y WMF'''
        self.TabNC = Tabla(50,self.Tabla_Prop_NC)
        self.TabWMF = Tabla(50, self.Tabla_Prop_WMF)
    
    def setupRaster2WMF(self):
               
        def setupLineEditButtonOpenRasterFileDialog (lineEditHolder, fileDialogHolder):
            '''Hace que solo0 se busquen formatos aceptados por GDAL'''
            lineEditHolder.setText (fileDialogHolder.getOpenFileName (QtGui.QDialog (), 'Open File',"", GdalTools_utils.FileFilter.allRastersFilter (),
                QtGui.QFileDialog.DontUseNativeDialog))
        def clickEventSelectorMapaRaster():
            '''Evento de click: selecciona mapa Raster'''
            setupLineEditButtonOpenRasterFileDialog (self.PathRaster2WMF, QFileDialog)
            
        #Funciones para Cargar variables
        def clickEventSelectorRaster():
            '''click para seleccionar un proyecto de cuenca'''
            #Pone el texto de la ruta 
            setupLineEditButtonOpenRasterFileDialog(self.PathRaster2WMF, QFileDialog)
        
        #Funcion para pasar variable raster a WMF           
        def handleClickConnectRaster2WMF():
            #Chequeos de variables
            Valor = -9
            if len(self.NameRaster2WMF.text())<2:
                self.iface.messageBar().pushMessage (u'Hydro-SIG:', u'Debe ingresar un nombre para el mapa raster a convertir.',
                    level=QgsMessageBar.WARNING, duration=3)
                return 1            
            try:
                Valor = float(self.PathRaster2WMF.text())
                EsNumero = True
            except:
                EsNumero = False
            if EsNumero:
                try:
                    Valor = float(self.PathRaster2WMF.text())
                    PathRaster = 'noMapa'
                    PathOrValue = 'Value'
                    print Valor
                except:
                    self.iface.messageBar().pushMessage (u'Hydro-SIG:', u'Debe seleccionar un mapa raster para ser convertido a la cuenca.',
                        level=QgsMessageBar.WARNING, duration=3)
                    return 1
            else:
                PathRaster = self.PathRaster2WMF.text()
                PathOrValue = 'Path'
            #Parametros para la conversion
            Nombre = self.NameRaster2WMF.text()
            Grupo = self.ComboBoxRaster2WMF.currentText()
            #Conversion, convierte la variable y actualiza el diccionario.
            Retorno = self.HSutils.Basin_Raster2WMF(Nombre, PathRaster, Grupo, PathOrValue, Valor)
            self.TabWMF.NewEntry(self.HSutils.DicBasinWMF[Nombre],Nombre, self.Tabla_Prop_WMF)
            if Retorno == 0:
                self.iface.messageBar().pushInfo (u'Hydro-SIG:', u'El mapa raster ha ingresado a WMF.')
            else:
                self.iface.messageBar().pushMessage (u'Hydro-SIG:', 
                    u'El mapa raster no ha ingresado a WMF.',
                    level=QgsMessageBar.WARNING, duration=3)
        
        #Habilita botones.
        self.ButtonPathRaster2WMF.clicked.connect(clickEventSelectorMapaRaster)
        self.Button_Raster2WMF.clicked.connect(handleClickConnectRaster2WMF)
        
    def setupNcVariables(self):
        '''Conjunto de herramientas para gestionar las variables del NC'''
        
        #Llena de datos los combobox 
        ListaMetodos = ['No metodo','moda','media','min','P10','P25','P50','P75','P90','max']
        map(self.ComboMethod4Conversion.addItem,ListaMetodos)
        #Lista de grupos posibles para una variable 
        #ListaGrupos = ['base','Geomorfo','SimHidro','Hidro']
        ListaGrupos = ['Geomorfo','Hidro']
        map(self.ComboBoxNewWMFVarGroup.addItem, ListaGrupos)
        #Lista de unidades de conversion 
        ListaUnidades = ['Celdas','Laderas','Canales']
        map(self.ComboConversionUnits.addItem, ListaUnidades)
        
        def clickEventEvalString():
            '''Funcion que eavlua una expresion escrita'''
            #toma la expresion y la evalua
            exp = self.LineaComando.text().strip()
            Var = self.HSutils.ExpressionParser(exp)
            EsCuenca = False
            if Var.size == self.HSutils.cuenca.ncells: EsCuenca = True
            #conversion si solo si la nueva variable tiene el ncells de la cuenca
            Agregado = self.ComboConversionUnits.currentText().strip()
            if EsCuenca and Agregado <> 'Celdas':
                #Revisar si es laderas o canales y metodo
                Metodo = self.ComboMethod4Conversion.currentText().strip()
                #Agrega
                Var = self.HSutils.BasinConvert2HillsOrChannels(Var, Metodo, Agregado)
                SiAgrego = True
            else:
                SiAgrego = False
            #Si es al nc sobre-escribe una entrada del diccionario
            if self.Radio2NcVar.isChecked():
                VarDestinoName = self.VarFromNC.currentText().encode()
                Capa = self.ObjectiveLayer.value()
                #Solo pasa al NC si el tamano es el de la cuenca
                if EsCuenca:
                    #Busca si la variable tiene una o mas dimensiones
                    try:
                        self.HSutils.DicBasinNc[VarDestinoName]['var'][Capa-1] = np.copy(Var)
                    except:
                        self.HSutils.DicBasinNc[VarDestinoName]['var'] = np.copy(Var)
                    #Establece a la variable nueva como no guardarda
                    self.HSutils.DicBasinNc[VarDestinoName]['saved'] = False
                    self.TabNC.EditedEntry(VarDestinoName, self.Tabla_Prop_NC)
                    #Mensaje de exito 
                    self.iface.messageBar().pushMessage (u'Hydro-SIG:', 
                        u'La variable '+VarDestinoName+' ha sido actualizada en NC',
                        level=QgsMessageBar.INFO, duration=3) 
                #Si no es cuenca da un mensaje de alerta
                else:
                    #Mensaje de exito o error 
                    self.iface.messageBar().pushMessage (u'Hydro-SIG:', 
                        u'El cálculo arroja una variable que no contiene '+str(self.HSutils.cuenca.ncells)+' celdas. Favor enviar a WMF.',
                        level=QgsMessageBar.WARNING, duration=3) 
            #si es una variable nueva la tira al WMF
            elif self.Radio2WMF.isChecked():
                #Variable al diccionario
                VarDestinoName = self.lineEditNewVarName.text().strip()
                Grupo = self.ComboBoxNewWMFVarGroup.currentText().strip().encode()
                self.HSutils.DicBasinWMF.update({VarDestinoName:
                    {'nombre':VarDestinoName,
                    'tipo':Var.dtype.name,
                    'shape':Var.shape,
                    'raster':EsCuenca,
                    'basica': False,
                    'categoria': Grupo,
                    'var': np.copy(Var),
                    'saved':False}})
                #Actualiza la tabla 
                self.TabWMF.NewEntry(self.HSutils.DicBasinWMF[VarDestinoName],
                    VarDestinoName, self.Tabla_Prop_WMF)
                #Mensaje de exito o error 
                self.iface.messageBar().pushMessage (u'Hydro-SIG:', 
                    u'La variable '+VarDestinoName+' ha sido cargado a la tabla WMF',
                    level=QgsMessageBar.INFO, duration=3) 
            #Pone de nuevo la capa destino igual a 0
            self.ObjectiveLayer.setValue(0)
            
            
        
        def clickEventConvertVariable2NC():
            '''Convierte una variable clickeada a la tabla de NC con una transformacion'''
            selectedItems = self.Tabla_Prop_WMF.currentRow ()
            VarName1 = self.Tabla_Prop_WMF.item(selectedItems,0).text()
            print VarName1
            
            selectedItems = self.Tabla_Prop_NC.currentRow ()
            VarName1 = self.Tabla_Prop_NC.item(selectedItems,0).text()
            print VarName1
            
        
        self.Button_EditNcVariable.clicked.connect(clickEventEvalString)
    
    def setupRainfall(self):
        '''Conjunto de herramientas dispuestas para interpolar campos de precipitacion'''
      
        #ListaRadarDates = []
        #FechasRadar = []
      
        def setupLineEditButtonOpenShapeFileDialog (lineEditHolder, fileDialogHolder):
            '''Hace que cuando se busquen shapes solo se encuetren formatos vectoriales'''
            lineEditHolder.setText (fileDialogHolder.getOpenFileName (QtGui.QDialog (), "", "*", "Shapefiles (*.shp);;"))
            if ((os.path.exists (lineEditHolder.text ().strip ())) and (not (self.iface is None))):
                self.iface.addVectorLayer (lineEditHolder.text ().strip (), os.path.basename (lineEditHolder.text ()).strip (), "ogr")
        
        def setupLineEditButtonSaveFileDialog (lineEditHolder, fileDialogHolder):
            '''Pone la ruta elegida en el dialogo de texto para guardado'''
            lineEditHolder.setText (fileDialogHolder.getSaveFileName (QtGui.QDialog (), "Guardar (cargar) binario de lluvia", "*", "RainBin (*.bin);;"))
        
        def setupLineEditButtonOpenExcelFileDialog (lineEditHolder, fileDialogHolder):
            '''Hace que cuando se busquen shapes solo se encuetren formatos vectoriales'''
            lineEditHolder.setText (fileDialogHolder.getOpenFileName (QtGui.QDialog (), "", "*", "Excel (*.xlsx);;"))
        
        def setupLineEditButtonOpenRadarFileDialog (lineEditHolder, fileDialogHolder):
            '''Para cambiar la carpeta por defecto donde se buscan las imagenes de radar'''
            #Obtiene la ruta donde esta lo de radar
            OutputFolder = fileDialogHolder.getExistingDirectory(QtGui.QDialog(), "Cargador de Cuencas", "/tmp/", QFileDialog.ShowDirsOnly)
            lineEditHolder.setText (OutputFolder)
        
        def clickEventSelectorMapaPuntosPluvio():
            '''Evento de click: selecciona el shp con las estaciones'''
            #Reinicia el vector con los nombres de las variables y abre el archivo de puntos
            self.HSutils.Interpol_Columnas = []
            setupLineEditButtonOpenShapeFileDialog (self.PathInHydro_Interpol_Pluvios, QFileDialog)
            #Esculca los nombres de las columnas del mismo.
            Path2Vect = self.PathInHydro_Interpol_Pluvios.text().strip()
            self.HSutils.Interpol_GetFields(Path2Vect)
            #Llena las opciones del ComboBox de Ids para elegir
            for l in self.HSutils.Interpol_Columnas:
                self.comboBox_Interpol.addItem(l)
            
        def clickEventSelectorArchivoExcel():
            '''Evento de click: selecciona el archivo de excel con los datos de pracipitacion a interpolar'''
            #Busca el archivo
            setupLineEditButtonOpenExcelFileDialog (self.PathInHydro_Interpol_Serie, QFileDialog)
            #Encuentra la fecha inicio, fin y paso de tiempo, sugiere esos datos al usuario
            Path2Excel = self.PathInHydro_Interpol_Serie.text().strip()
            self.HSutils.Interpol_GetDateTimeParams(Path2Excel)
            #Pone las fechas
            Date = QtCore.QDateTime(self.HSutils.Interpol_fi)
            self.Interpol_DateTimeStart.setDateTime(Date)
            Date = QtCore.QDateTime(self.HSutils.Interpol_ff)
            self.Interpol_DateTimeEnd.setDateTime(Date)
            #Pone el intervalo de tiempo de interpolacion
            self.Interpol_SpinBox_delta.setValue(self.HSutils.Interpol_fd.total_seconds())
        
        def clickEventSelectorRadarFiles():
            '''Evento de seleccion de la carpeta contenedora de los archivos de radar'''
            #Obtiene la ruta donde estan los archivos de radar
            setupLineEditButtonOpenRadarFileDialog(self.PathInHydro_Radar, QFileDialog)
            self.Path2Radar = self.PathInHydro_Radar
            #Actualiza lista con variables del radar
            self.ListaRadarDates = glob.glob(self.Path2Radar.text().strip()+'/*.nc')
            self.ListaRadarDates.sort()
            self.FechasRadar = []
            for L in self.ListaRadarDates:
                try:
                    self.FechasRadar.append(dt.datetime.strptime(L[-23:-11],'%Y%m%d%H%M'))
                except:
                    pass
            #print len(self.ListaRadarDates)
            #print self.ListaRadarDates[0]
            #print self.ListaRadarDates[-1]
        
        def clickEventSelectorArchivoBinarioLluvia():
            '''Selecciona la ruta en donde se guardara el binario de salida.'''
            #pone el camino del archivo con la lluvia de la cuenca
            setupLineEditButtonSaveFileDialog(self.PathOutHydro_Interpol,QFileDialog)
            #Trata de leer datos de lluvia en caso de que ya existan
            try:
                PathData = self.PathOutHydro_Interpol.text().strip()
                self.HSplots = HSplots.PlotRainfall(PathData)
            except:
                pass
        
        def clickEventSelectorArchivoBinarioLluviaRadar():
            '''Selecciona la ruta en donde se guardara el binario de salida.'''
            #pone el camino del archivo con la lluvia de la cuenca
            setupLineEditButtonSaveFileDialog(self.PathOutHydro_Radar,QFileDialog)
            #Trata de leer datos de lluvia en caso de que ya existan
            try:
                PathData = self.PathOutHydro_Radar.text().strip()
                self.HSplots = HSplots.PlotRainfall(PathData)
            except:
                pass
        
        def clickEventEjecutarInterpolacion():
            '''Interpola los campos de precipitacion con los parametros ingresados.'''
            #Toma los parametros para la interpolacion
            Path2Shp = self.PathInHydro_Interpol_Pluvios.text().strip()
            Campo2Read = self.comboBox_Interpol.currentText()
            fi = self.Interpol_DateTimeStart.dateTime().toPyDateTime()
            ff = self.Interpol_DateTimeEnd.dateTime().toPyDateTime()
            fd = self.Interpol_SpinBox_delta.value()
            expo = self.Interpol_SpinBox_expIDW.value()
            PathOut = self.PathOutHydro_Interpol.text().strip()
            #Interpola para la cuenca seleccionada
            self.HSutils.Interpol_GetInterpolation(Path2Shp,Campo2Read,fi,ff,fd,expo, PathOut)
            #Trata de leer datos de lluvia en caso de que ya existan
            try:
                PathData = self.PathOutHydro_Interpol.text().strip()
                self.HSplots = HSplots.PlotRainfall(PathData)
            except:
                pass
            #Aviso de existo
            self.iface.messageBar().pushInfo(u'HidroSIG:',u'Interpolacion de campos de precipitacion realizada con exito')
        
        def clickEventEjecutarConversionRadar():
            '''Convierte campos de radar en la cuenca usando el archivo netCDF que se encuentran en la ruta especificada'''
            #Toma los parametros para la interpolacion
            PathRadar = self.PathInHydro_Radar.text().strip()
            fi = self.Interpol_DateTimeStart_Radar.dateTime().toPyDateTime()
            ff = self.Interpol_DateTimeEnd_Radar.dateTime().toPyDateTime()
            Step = '%dS' % self.Interpol_SpinBox_delta_Radar.value()
            PathOut = self.PathOutHydro_Radar.text().strip()
            #Obtiene las fechas para conversion y la lista de valores
            self.HSutils.Radar_FechasProcess(fi, ff, Step, self.FechasRadar, self.ListaRadarDates)
            #Interpola para la cuenca seleccionada
            self.HSutils.Radar_Conver2Basin(PathOut, Step)
            #Trata de leer datos de lluvia en caso de que ya existan
            try:
                PathData = self.PathOutHydro_Radar.text().strip()
                self.HSplots = HSplots.PlotRainfall(PathData)
            except:
                pass
            #Aviso de existo
            self.iface.messageBar().pushInfo(u'HidroSIG:',u'Conversión de campos de radar realizada con exito')
        
        def clickEventViewSerieRainfall():
            '''Genera y visualiza la grafica de lluvia interpolada para la cuenca'''
            #Hace la figura
            PathFigure = '/tmp/HydroSED/Plots_Rainfall/RainfallPlot.html'
            self.HSplots.Plot_Rainfall(PathFigure)
            #Set de la ventana que contiene la figura.
            self.VistaRainWeb = QWebView(None)
            self.VistaRainWeb.load(QUrl.fromLocalFile(PathFigure))
            self.VistaRainWeb.setWindowTitle('Precipitacion media en la cuenca')
            self.VistaRainWeb.setMinimumWidth(1100)
            self.VistaRainWeb.setMaximumWidth(3000)
            self.VistaRainWeb.setMinimumHeight(100)
            self.VistaRainWeb.setMaximumHeight(400)
            self.VistaRainWeb.show()
        
        def clickEventViewHistogramRainfall():
            '''Genera y visualiza la grafica de lluvia interpolada para la cuenca'''
            #Hace la figura
            PathFigure = '/tmp/HydroSED/Plots_Rainfall/RainfallHistogram.html'
            self.HSplots.Plot_Histogram(PathFigure)
            #Set de la ventana que contiene la figura.
            self.VistaRainWeb = QWebView(None)
            self.VistaRainWeb.load(QUrl.fromLocalFile(PathFigure))
            self.VistaRainWeb.setWindowTitle('Histograma precipitacion')
            self.VistaRainWeb.setMinimumWidth(200)
            self.VistaRainWeb.setMaximumWidth(400)
            self.VistaRainWeb.setMinimumHeight(400)
            self.VistaRainWeb.setMaximumHeight(400)
            self.VistaRainWeb.show()
        
        def clickEventViewMediaMensualRainfall():
            '''Genera y visualiza la grafica de la lluvia media mensual en la cuenca'''
            #Hace la figura
            PathFigure = '/tmp/HydroSED/Plots_Rainfall/RainfallMediaMensual.html'
            self.HSplots.Plot_MediaMensual(PathFigure)
            #Set de la ventana que contiene la figura.
            self.VistaRainWeb = QWebView(None)
            self.VistaRainWeb.load(QUrl.fromLocalFile(PathFigure))
            self.VistaRainWeb.setWindowTitle('Media mensual multi-anual de precipitacion')
            self.VistaRainWeb.setMinimumWidth(200)
            self.VistaRainWeb.setMaximumWidth(800)
            self.VistaRainWeb.setMinimumHeight(100)
            self.VistaRainWeb.setMaximumHeight(400)
            self.VistaRainWeb.show()
        
        def clickEventGetAcumRainfall():
            '''Obtiene el acumulado de lluvia en el periodo especifico'''
            #Punto inicial y final 
            Path = self.PathOutHydro_Interpol.text().strip()
            inicio = int(self.spinBoxCampoInicio.value())
            fin = int(self.spinBoxCampoFin.value())
            #Obtiene el campo acumulado para el periodo seleccionado
            self.HSutils.Interpol_GetRainfallAcum(Path, inicio, fin)
            #Actualiza la tabla WMF
            k = 'Lluvia'
            self.TabWMF.NewEntry(self.HSutils.DicBasinWMF[k],k, self.Tabla_Prop_WMF)
            #Aviso de existo
            self.iface.messageBar().pushInfo(u'HidroSIG:',u'Los campos de han acumulado en la variable Lluvia de la tabla WMF')
        
        def clickEventGetAcumRainfallRadar():
            '''Obtiene el acumulado de lluvia en el periodo especifico'''
            #Punto inicial y final 
            Path = self.PathOutHydro_Radar.text().strip()
            inicio = int(self.spinBoxCampoInicio_Radar.value())
            fin = int(self.spinBoxCampoFin_Radar.value())
            #Obtiene el campo acumulado para el periodo seleccionado
            self.HSutils.Interpol_GetRainfallAcum(Path, inicio, fin)
            #Actualiza la tabla WMF
            k = 'Lluvia'
            self.TabWMF.NewEntry(self.HSutils.DicBasinWMF[k],k, self.Tabla_Prop_WMF)
            #Aviso de existo
            self.iface.messageBar().pushInfo(u'HidroSIG:',u'Los campos de han acumulado en la variable Lluvia de la tabla WMF')
        
        
        #Botones de set de interpolacion
        self.Boton_HidroLoad_Pluvios.clicked.connect(clickEventSelectorMapaPuntosPluvio)
        self.Boton_HidroLoad_Serie.clicked.connect(clickEventSelectorArchivoExcel)
        self.Button_HidroSaveInterpol.clicked.connect(clickEventSelectorArchivoBinarioLluvia)
        #Botones de ejecucion
        self.Button_Ejec_HidroInterpol.clicked.connect(clickEventEjecutarInterpolacion)
        self.Button_InterpolViewField.clicked.connect(clickEventGetAcumRainfall)
        #Botones de visualizacion
        self.Button_InterpolSerieView.clicked.connect(clickEventViewSerieRainfall)
        self.Button_InterpolHistogram.clicked.connect(clickEventViewHistogramRainfall)
        self.Button_InterpolCiclo.clicked.connect(clickEventViewMediaMensualRainfall)
        
        #botones set de radar
        self.Boton_HidroLoad_RadarData.clicked.connect(clickEventSelectorRadarFiles)
        self.Button_HidroSaveRadar.clicked.connect(clickEventSelectorArchivoBinarioLluviaRadar)
        #Boton de ejecucion de radar
        self.Button_Ejec_HidroRadar.clicked.connect(clickEventEjecutarConversionRadar)
        #Botones de visualizacion
        self.Button_InterpolCiclo_Radar.clicked.connect(clickEventViewMediaMensualRainfall)
        self.Button_InterpolHistogram_Radar.clicked.connect(clickEventViewHistogramRainfall)
        self.Button_InterpolSerieView_Radar.clicked.connect(clickEventViewSerieRainfall)
        self.Button_InterpolViewRadar.clicked.connect(clickEventGetAcumRainfallRadar)
    
    def setupSimulation(self):
        '''Herramientas para gestionar la simulacion hidrologica con la cuenca cargada'''
        
        def setupLineEditButtonOpenExcelCaudalFileDialog (lineEditHolder, fileDialogHolder):
            '''Que solo busque los archivos con esa extension .xlsx'''
            lineEditHolder.setText (fileDialogHolder.getOpenFileName (QtGui.QDialog (), "", "*", "Excel (*.xlsx);;")) 
            
        def setupLineEditButtonOpenBinFileDialog (lineEditHolder, fileDialogHolder):
            '''Que solo busque los archivos con esa extension .bin'''
            lineEditHolder.setText (fileDialogHolder.getOpenFileName (QtGui.QDialog (), "", "*", "Binarios (*.bin);;")) 
        
        def setupLineEditButtonSaveBinStoFileDialog (lineEditHolder, fileDialogHolder):
            '''Que solo busque los archivos con esa extension .bin'''
            lineEditHolder.setText (fileDialogHolder.getSaveFileName (QtGui.QDialog (), "Guardar estados almacenamiento", "*", "Binario (*.StObin);;")) 
            
        def setupLineEditButtonOpenBinStoFileDialog (lineEditHolder, fileDialogHolder):
            '''Que solo busque los archivos con esa extension .bin'''
            lineEditHolder.setText (fileDialogHolder.getOpenFileName (QtGui.QDialog (), "Buscar estados de almacenamiento", "*", "Binarios (*.StObin);;")) 
        
        def changeEventUpdateScalarParameters():
            '''Actualiza los parametros escalares de las tablas de acuerdo al set seleccionado'''
            #Obtiene el nombre de la param seleccionada
            key = self.ParamNamesCombo.currentText().strip().encode()
            #Itera en el diccionario de param de la cuenca 
            print key
            if key <> '':
                for c,values in enumerate(self.HSutils.DicParameters[key]['var'][:11]):
                    codigo = 'self.Param'+str(c+1)+'.setValue('+str(values)+')'
                    eval(codigo)
                #for c,values in enumerate(self.HSutils.DicParameters[key]['var'][11:]):
                 #   codigo = 'self.ParamExp'+str(c+1)+'.setValue('+str(values)+')'
                  #  eval(codigo)
        
        def clickEventSetPath2States():
            '''Eventod e click para establecer el path de archivo binario con estados'''
            setupLineEditButtonOpenBinStoFileDialog (self.StateFilePath, QFileDialog)
        
        def clickEventUpdateParamMapValues():
            '''Muestra en los campos de simulacion el valor medio de los mapas de simulacion'''
            VarNames = ['h1_max','h3_max', 'v_coef','v_coef','v_coef','v_coef','h_coef',
                'h_coef','h_coef','h_coef', 'Krus','Crus','Prus',
                'PArLiAc','PArLiAc','PArLiAc','h_exp','h_exp','h_exp','h_exp']
            Ejes = [0,0,0,1,2,3,0,1,2,3,0,0,0,0,1,2,0,1,2,3]
            for name, i, eje in zip(VarNames, range(1,21), Ejes):
                Campo = getattr(self, 'ParamVal'+str(i))
                Campo.setMinimum(-99999)
                Campo.setValue(0)
                try:
                    Value = self.HSutils.DicBasinNc[name]['var'].mean(axis = 1)[eje]
                except:
                    try:
                        Value = self.HSutils.DicBasinNc[name]['var'].mean()
                    except:
                        Value = -9999
                Campo.setValue(Value)
        
        def clickEventAddNewParamSet():
            '''Agrega un nuevo conjunto de parametros en el proyecto de cuenca'''
            #Obtiene los parametros
            PathNC = self.lineEditRutaCuenca.text().strip()
            ParamName = self.ParamName.text().strip().encode()
            #Itera para los parametros escalares y de sedimentos
            ListaParam = []
            for i in range(1,12):
                ListaParam.append(getattr(self, 'Param'+str(i)).value())
            #Itera en los exponentes
            #for i in range(1,5):
             #   ListaParam.append(getattr(self, 'ParamExp'+str(i)).value())    
            #Mete el set nuevo de calibracion
            ListaParam.extend([0,0,0,0])
            self.HSutils.Sim_SaveParameters(PathNC, ParamName, ListaParam)
            #Actualiza la lista de parametros 
            self.ParamNamesCombo.clear()
            for k in self.HSutils.DicParameters.keys():
                self.ParamNamesCombo.addItem(k)
            #mensaje de exito
            self.iface.messageBar().pushMessage (u'Hydro-SIG:', u'La parametrización: ' + ParamName+' se ha guardado en el proyecto',
                level=QgsMessageBar.INFO, duration=3)
                                  
        def clickEventSelectorArchivoExcelCaudales():
            '''Evento de click: selecciona el archivo de excel con los datos caudal observado a utilizar'''
            #Busca el archivo
            setupLineEditButtonOpenExcelCaudalFileDialog (self.PathinSimu_Qobs, QFileDialog) 
            #vacia la lista desplegable 
            if self.Segunda_carga_Qobs == True:
                self.comboBox_Selec_Qobs.clear()
                
            for l in self.HSutils.Sim_GetQobsInfo(self.PathinSimu_Qobs.text().strip())[0]:
                self.comboBox_Selec_Qobs.addItem(str(l))
                
            self.Segunda_carga_Qobs = True
        def clickEventSelectorArchivoExcelCaudalesSed():
            '''Evento de click: selecciona el archivo de excel con los datos de Caudal solido a usar'''
            #Busca el archivo
            setupLineEditButtonOpenExcelCaudalFileDialog (self.PathinSimu_Qobs_Sed, QFileDialog) 
            
            #vacia la lista desplegable    
            if self.Segunda_carga_Qobs_Sed == True:
                self.comboBox_Selec_Qobs_Sed.clear()
            
            for l in self.HSutils.Sim_GetQobsInfo(self.PathinSimu_Qobs_Sed.text().strip())[0]:
                self.comboBox_Selec_Qobs_Sed.addItem(str(l))
                
            self.Segunda_carga_Qobs_Sed = True
            
        def clickEventSelectorArchivoBinarioSimulacion():
            '''Evento de click: selecciona el archivo binario con la precipitacion para correr el modelo'''
            #Busca el archivo
            setupLineEditButtonOpenBinFileDialog (self.PathinSimu_Precipitacion,QFileDialog)
            #Obtiene el Dt de modelacion a partir de la lluvia
            PathRBin, PathRHdr = HSutils.wmf.__Add_hdr_bin_2route__(self.PathinSimu_Precipitacion.text().strip())
            RainStruct = HSutils.wmf.read_rain_struct(PathRHdr)
            #Obtiene el Dt
            Dt = RainStruct.index[1] - RainStruct.index[0]
            Texto = '%d' % Dt.seconds
            self.DtShowFrame.setText(Texto)
            #SEt para la figura de lluvia
            try:
                PathData = self.PathinSimu_Precipitacion.text().strip()
                self.HSplots = HSplots.PlotRainfall(PathData)
            except:
                pass
                
        def clickEventViewSerieSimuRainfall():
            '''Genera y visualiza la grafica de lluvia interpolada para la cuenca con la que se va a simular'''
            #Hace la figura
            PathFigure = '/tmp/HydroSED/Plots_Rainfall/RainfallPlotSimu.html'
            self.HSplots.Plot_Rainfall(PathFigure)
            #Set de la ventana que contiene la figura.
            self.VistaRainWeb = QWebView(None)
            self.VistaRainWeb.load(QUrl.fromLocalFile(PathFigure))
            self.VistaRainWeb.setWindowTitle('Precipitacion media en la cuenca')
            self.VistaRainWeb.setMinimumWidth(1100)
            self.VistaRainWeb.setMaximumWidth(3000)
            self.VistaRainWeb.setMinimumHeight(100)
            self.VistaRainWeb.setMaximumHeight(400)
            self.VistaRainWeb.show()
            
        def clickEventViewSerieQobs():
            PathFigure =  '/tmp/HydroSED/Plots_Rainfall/QobsPlotSimu.html'
            self.PathQobs = self.PathinSimu_Qobs.text().strip()
            #Obtiene el id de la estación 
            id_est = int(self.comboBox_Selec_Qobs.currentText().encode())
            #Obtiene la serie de caudales 
            self.DataQ = self.HSutils.Sim_GetQobsInfo(self.PathQobs)[1]
            self.HSplotsCaudal = HSplots.PlotCaudal()
            self.HSplotsCaudal.Plot_Caudal(PathFigure,self.DataQ,id_est,'blue')
            #Set de la ventana que contiene la figura.
            self.VistaQobsWeb = QWebView(None)
            self.VistaQobsWeb.load(QUrl.fromLocalFile(PathFigure))
            self.VistaQobsWeb.setWindowTitle('Serie de Caudales medios')
            self.VistaQobsWeb.setMinimumWidth(1100)
            self.VistaQobsWeb.setMaximumWidth(3000)
            self.VistaQobsWeb.setMinimumHeight(100)
            self.VistaQobsWeb.setMaximumHeight(400)
            self.VistaQobsWeb.show()
            
            
        def clickEventViewSerieQobsSed():
            PathFigure =  '/tmp/HydroSED/Plots_Rainfall/QobsSedPlotSimu.html'
            PathQobs = self.PathinSimu_Qobs_Sed.text().strip()
            #Obtiene el id de la estación 
            id_est = int(self.comboBox_Selec_Qobs_Sed.currentText().encode())
            #Obtiene la serie de caudales 
            self.DataQsed = self.HSutils.Sim_GetQobsInfo(PathQobs)[1]
            self.HSplotsCaudal = HSplots.PlotCaudal()
            self.HSplotsCaudal.Plot_Caudal(PathFigure,self.DataQsed,id_est,'goldenrod')
            #Set de la ventana que contiene la figura.
            self.VistaQobsWeb = QWebView(None)
            self.VistaQobsWeb.load(QUrl.fromLocalFile(PathFigure))
            self.VistaQobsWeb.setWindowTitle('Serie de Caudales medios')
            self.VistaQobsWeb.setMinimumWidth(1100)
            self.VistaQobsWeb.setMaximumWidth(3000)
            self.VistaQobsWeb.setMinimumHeight(100)
            self.VistaQobsWeb.setMaximumHeight(400)
            self.VistaQobsWeb.show()
              
        def clickEventViewSimStorageSeries():
            '''Hace plot de las condiciones medias de almacenamiento de la cuenca'''
            #Path de la figura
            PathFigure =  '/tmp/HydroSED/Plots_Sim/SimStoredConditions.html'
            #Obtiene el path de los datos
            Path = self.Where2SaveStates.text().strip()
            PathBin, PathHdr = HSutils.wmf.__Add_hdr_bin_2route__(Path, storage = True)
            #Invoca la figura 
            self.SimStoragePlots.Plot_Storages(PathHdr, PathFigure)
            #Set de la ventana que contiene la figura.
            self.VistaQobsWeb = QWebView(None)
            self.VistaQobsWeb.load(QUrl.fromLocalFile(PathFigure))
            self.VistaQobsWeb.setWindowTitle('Series de almacenamiento medio')
            self.VistaQobsWeb.setMinimumWidth(600)
            self.VistaQobsWeb.setMaximumWidth(600)
            self.VistaQobsWeb.setMinimumHeight(500)
            self.VistaQobsWeb.setMaximumHeight(700)
            self.VistaQobsWeb.show()
        
        def __ParseCalibValues__():
            '''Obtiene una lista de los param de calibracion a partir de los elem que estan en al interfaz'''
            Calibracion = []
            for i in range(1,12):
                Calibracion.append(getattr(self, 'Param'+str(i)).value())
            return Calibracion
        
        def __Dates2Start_Nsteps__(fi,ff,PathRain):
            '''A partir de la fecha inicio y fin obtien el paso de inicio y la cantidad de pasos'''
            #Obtiene el path para header y para el binario, lee la estructura de la lluvia
            PathRBin, PathRHdr = HSutils.wmf.__Add_hdr_bin_2route__(PathRain)
            RainStruct = HSutils.wmf.read_rain_struct(PathRHdr)
            #Obtiene el Dt
            Dt = RainStruct.index[1] - RainStruct.index[0]
            Dt = Dt.seconds
            #Obtiene punto de inicio 
            Start = RainStruct.index.get_loc(fi.strftime('%Y-%m-%d %H:%M'))
            End = RainStruct.index.get_loc(ff.strftime('%Y-%m-%d %H:%M'))
            Nsteps = End - Start
            #Retorna
            return Start, Nsteps, PathRBin, Dt
        
        def __UpdateVar2WMFTable__(varName, varValue):
            '''Agrega una variable a la tabla de WMF'''
            self.HSutils.DicBasinWMF.update({varName:
                    {'nombre':varName,
                    'tipo':varValue.dtype.name,
                    'shape':varValue.shape,
                    'raster':True,
                    'basica': False,
                    'categoria': 'Hidro',
                    'var': np.copy(varValue),
                    'saved':False}})
            self.TabWMF.NewEntry(self.HSutils.DicBasinWMF[varName],
                varName, self.Tabla_Prop_WMF)
        
        def clickEventSaveAlmacenamientos():
            '''Establece la ruta de guardado de los almacenamientos'''
            #Establece el punto de guardado
            setupLineEditButtonSaveBinStoFileDialog(self.Where2SaveStates, QFileDialog)
            Path = self.Where2SaveStates.text().strip()
            #Hace check del boton
            self.checkBox_Simu_MeanSpeed.setChecked(True)
            self.checkBox_Simu_MeanStorage.setChecked(True)
            HSutils.wmf.models.show_storage = 1
        
        def clickEventSetAlmacenamientos():
            '''Establece que clase de almacenamientos se van a usar para la simulacion'''
            #Itera por opciones y almacenamientos
            self.DictStatesOptions = {}
            for state in range(1,7):
                for option in range(1,4):
                    #Determina en que opcion esta 
                    boton = 'RadioO'+str(option)+'S'+str(state)
                    Estado = getattr(self, boton).isChecked()
                    #Si encontro el estado lo habilita
                    if Estado:
                        #Valor dado por el susuario
                        if option == 1:
                            Value = getattr(self, 'Storage'+str(state)+'Val').value()
                            self.HSutils.Sim_setStates_ConstValue(state, Value)
                            Valor = np.ones(self.HSutils.cuenca.ncells)*Value
                        #Valor tomado de un archivo de estados de almacenamiento
                        if option == 2:
                            FilePath = self.StateFilePath.text().strip()
                            Date = self.StateDate.dateTime().toPyDateTime()
                            Return, Valor = self.HSutils.Sim_setStates_FileValue(state, FilePath, Date)
                            Texto = '%.2f' % Return
                            Comando = 'self.Storage'+str(state)+'Val.setValue('+Texto+')'
                            eval(Comando)
                        #Valor tomado de la variable storage del NC
                        if option == 3:
                            Return, Valor = self.HSutils.Sim_setStates_NcValue(state)
                            Texto = '%.2f' % Return
                            Comando = 'self.Storage'+str(state)+'Val.setValue('+Texto+')'
                            eval(Comando)
                if state<6:
                    #Actualiza el diccionario NC de storage para visualizar variables
                    self.HSutils.DicBasinNc['storage']['var'][state-1] = Valor
            #Establece las condiciones de esa entrada en la tabla de NC
            self.HSutils.DicBasinNc['storage']['saved'] = False
            self.TabNC.EditedEntry('storage', self.Tabla_Prop_NC)
            #Mensaje de exito y cambio
            self.iface.messageBar().pushMessage (u'Hydro-SIG:', 
                u'Se han actualizado las condiciones de estado, mirar variable storage en Tabla NC',
                level=QgsMessageBar.INFO, duration=3) 
        
        def clickEventSimulationDeCuenca():
            '''Hace la simulacion hidrologica con el set de param seleccionados y los mapas propios de la cuenca'''
            #Obtiene lo que se necesita para ejecutar 
            PathRain = self.PathinSimu_Precipitacion.text().strip()
            Calib = __ParseCalibValues__()
            #Obtiene pubnto de inicio y cantidad de pasos 
            Fi = self.Simulacion_DateTimeStart.dateTime().toPyDateTime()
            Ff = self.Simulacion_DateTimeEnd.dateTime().toPyDateTime()
            Inicio, Npasos, PathBin, TimeDelta = __Dates2Start_Nsteps__(Fi,Ff,PathRain)
            #Exponenetes de funciones lineales o no lineales
            Exponenetes = []
            for i in [17,18,19]:
                Campo = getattr(self, 'ParamVal'+str(i))
                Exponenetes.append(Campo)
            #Banderas de ejecucion
            if self.checkBox_Simu_Sedimentos.isChecked():
                HSutils.wmf.models.sim_sediments = 1
            if self.checkBox_Simu_MeanSpeed.isChecked():
                HSutils.wmf.models.show_mean_speed = 1
            if self.checkBox_Simu_MeanStorage.isChecked():
                HSutils.wmf.models.show_storage = 1
            if self.checkBox_Simu_Retorno.isChecked():
                HSutils.wmf.models.retorno = 1
                HSutils.wmf.models.show_mean_retorno = 1
            #Guarda o no estados de almacenamiento
            pathStorage = None
            if self.checkBox_Simu_MeanSpeed.isChecked():
                pathStorage = self.Where2SaveStates.text().strip()
            #Simulacion hidrologica
            self.HSutils.Sim_Basin(Inicio, Npasos, Calib, PathBin, TimeDelta, Exponenetes, pathStorage)
            #Pone estados de almacenamiento finales en el NC
            self.HSutils.DicBasinNc['storage']['var'] = np.copy(self.HSutils.Sim_Storage)
            self.HSutils.DicBasinNc['storage']['saved'] = False
            self.TabNC.EditedEntry('storage', self.Tabla_Prop_NC)
            #Pone el campo acumulado de lluvia usado en el WMF
            __UpdateVar2WMFTable__('Sim_Rain',self.HSutils.Sim_RainfallField)
            #Campos de erosion en caso de que este se simule 
            if self.checkBox_Simu_Sedimentos.isChecked():
                __UpdateVar2WMFTable__('Sim_Erosion',self.HSutils.Sim_ErosionMap)
                __UpdateVar2WMFTable__('Sim_Deposit',self.HSutils.Sim_DepositionMap)
            #Mensaje de exito 
            self.iface.messageBar().pushMessage (u'Hydro-SIG:', 
                u'El modelo se ha ejecutado con exito',
                level=QgsMessageBar.INFO, duration=3) 
                
                
        def clickEventViewSerieQobsQsim():
            #Fechas de inicio y fin de simulacion
            self.f_ini = self.Simulacion_DateTimeStart.dateTime().toPyDateTime()
            self.f_fin = self.Simulacion_DateTimeEnd.dateTime().toPyDateTime()
            #llama el df de caudal observado 
            dfDataQsim = self.HSutils.Sim_Streamflow
            #identifica el tramo para graficar 
            id_tramo = str(self.spinBox.value())
            #Llama el df en el id de tramo que puso el usuario 
            self.dfDataQsim = dfDataQsim[id_tramo]
            #Llama la clase de plot caudal 
            self.HSplotsCaudal = HSplots.PlotCaudal()
            PathFigure =  '/tmp/HydroSED/Plots_Rainfall/QobsQsimPlotSimu.html'
            
            if self.checkBox_Simu_Qs.isChecked():
                #Selecciona el id que haya puesto el usuario
                id_est = int(self.comboBox_Selec_Qobs.currentText().encode())
                #Obtiene la serie de caudales 
                self.PathQobs = self.PathinSimu_Qobs.text().strip()
                #llama el df de caudal observado 
                self.DataQ = self.HSutils.Sim_GetQobsInfo(self.PathQobs)[1]
                #Crea el dataframe 
                self.dfDataQobs = self.DataQ[id_est][self.f_ini:self.f_fin]
            else:
                self.dfDataQobs = self.dfDataQsim*np.nan
            print self.HSutils.Sim_Rainfall
            self.HSplotsCaudal.Plot_Caudal_Simu(PathFigure,self.dfDataQobs,self.dfDataQsim,self.HSutils.Sim_Rainfall)
            #Set de la ventana que contiene la figura
            self.VistaQobsWeb = QWebView(None)
            self.VistaQobsWeb.load(QUrl.fromLocalFile(PathFigure))
            self.VistaQobsWeb.setWindowTitle('Series de Caudales')
            self.VistaQobsWeb.setMinimumWidth(1100)
            self.VistaQobsWeb.setMaximumWidth(3000)
            self.VistaQobsWeb.setMinimumHeight(200)
            self.VistaQobsWeb.setMaximumHeight(500)
            self.VistaQobsWeb.show()
            
        def clickEventViewCDCQobsQsim():
            self.f_ini = self.Simulacion_DateTimeStart.dateTime().toPyDateTime()
            self.f_fin = self.Simulacion_DateTimeEnd.dateTime().toPyDateTime()
            self.HSplotsCaudal = HSplots.PlotCaudal()
            PathFigure1 =  '/tmp/HydroSED/Plots_Rainfall/CDCQobsQsimuPlotSimu.html'
            #llama el df de caudal observado 
            dfDataQsim = self.HSutils.Sim_Streamflow
            #identifica el tramo para graficar 
            id_tramo = str(self.spinBox_2.value())
            #Llama el df en el id de tramo que puso el usuario 
            self.dfDataQsim = dfDataQsim[id_tramo]
            
            if self.checkBox_Simu_Qs_2.isChecked():
                #Selecciona el id que haya puesto el usuario
                id_est = int(self.comboBox_Selec_Qobs.currentText().encode())
                #Obtiene la serie de caudales 
                self.PathQobs = self.PathinSimu_Qobs.text().strip()
                #llama el df de caudal observado 
                self.DataQ = self.HSutils.Sim_GetQobsInfo(self.PathQobs)[1]
                #Crea el dataframe 
                self.dfDataQobs = self.DataQ[id_est][self.f_ini:self.f_fin]
            else:
                self.dfDataQobs = self.dfDataQsim*np.nan
            
            #Obtiene la serie de caudales 
            self.HSplotsCaudal.Plot_CDC_caudal(PathFigure1,self.dfDataQobs,self.dfDataQsim)
            #Set de la ventana que contiene la figura
            self.VistaQobsWeb = QWebView(None)
            self.VistaQobsWeb.load(QUrl.fromLocalFile(PathFigure1))
            self.VistaQobsWeb.setWindowTitle('Curvas de duracion de Caudales')
            self.VistaQobsWeb.setMinimumWidth(400)
            self.VistaQobsWeb.setMaximumWidth(500)
            self.VistaQobsWeb.setMinimumHeight(400)
            self.VistaQobsWeb.setMaximumHeight(400)
            self.VistaQobsWeb.show()
            
        def clickEventViewSerieSedimentos():
            #Fechas de inicio y fin de simulacion
            self.f_ini = self.Simulacion_DateTimeStart.dateTime().toPyDateTime()
            self.f_fin = self.Simulacion_DateTimeEnd.dateTime().toPyDateTime()
            #llama el df de caudal observado 
            dfDataQsim = self.HSutils.Sim_Sediments
            #identifica el tramo para graficar 
            id_tramo = str(self.spinBox_3.value())
            #Llama el df en el id de tramo que puso el usuario 
            dfDataQsim = dfDataQsim[id_tramo]
            #Llama la clase de plot caudal 
            self.HSplotsCaudal = HSplots.PlotCaudal()
            self.Arenas=dfDataQsim['Sand']
            self.Limos=dfDataQsim['Lime']
            self.Arcillas= dfDataQsim['Clay']
            self.Sed_total = self.Arenas + self.Limos + self.Arcillas 
            
            PathFigure =  '/tmp/HydroSED/Plots_Rainfall/QSedPlotSimu.html'
            
            if self.checkBox_Simu_Qs_3.isChecked():
                #Selecciona el id que haya puesto el usuario
                id_est = int(self.comboBox_Selec_Qobs_Sed.currentText().encode())
                #Obtiene la serie de caudales 
                self.PathQobs = self.PathinSimu_Qobs_Sed.text().strip()
                #llama el df de caudal observado 
                self.DataQ = self.HSutils.Sim_GetQobsInfo(self.PathQobs)[1]
                #Recorta el DF a las fechas 
                self.dfDataQobs = self.DataQ[id_est][self.f_ini:self.f_fin]
            else:
                self.dfDataQobs = self.Arenas*np.nan
            #print self.HSutils.Sim_Rainfall
            self.HSplotsCaudal.Plot_Series_Sedimentos(PathFigure,self.Arenas,self.Limos,self.Arcillas,self.Sed_total,self.dfDataQobs)
            #Set de la ventana que contiene la figura
            self.VistaQobsWeb = QWebView(None)
            self.VistaQobsWeb.load(QUrl.fromLocalFile(PathFigure))
            self.VistaQobsWeb.setWindowTitle('Series de Caudales solidos')
            self.VistaQobsWeb.setMinimumWidth(1100)
            self.VistaQobsWeb.setMaximumWidth(3000)
            self.VistaQobsWeb.setMinimumHeight(200)
            self.VistaQobsWeb.setMaximumHeight(500)
            self.VistaQobsWeb.show()                                  
                                  
                                  
         
        self.ButtonSimCalib2Nc.clicked.connect(clickEventAddNewParamSet)    
        self.tabPanelDockOpciones.currentChanged.connect(clickEventUpdateParamMapValues)
        self.ParamNamesCombo.currentIndexChanged.connect(changeEventUpdateScalarParameters) 
             
        self.Boton_SimuSelec_Qobs.clicked.connect(clickEventSelectorArchivoExcelCaudales)
        self.Boton_SimuSelec_Qobs_Sed.clicked.connect(clickEventSelectorArchivoExcelCaudalesSed)
        self.Boton_SimuSelec_Binario.clicked.connect(clickEventSelectorArchivoBinarioSimulacion)
        
        self.Boton_Visualizar_Binario.clicked.connect(clickEventViewSerieSimuRainfall)
        self.Boton_Visualizar_Qobs.clicked.connect(clickEventViewSerieQobs)
        self.Boton_Visualizar_Qobs_Sed.clicked.connect(clickEventViewSerieQobsSed)
        self.ButtonSim_ViewStorage.clicked.connect(clickEventViewSimStorageSeries)
        
        self.ButtonSelectPath2States.clicked.connect(clickEventSetPath2States)
        self.BotonSelectSaveStateRute.clicked.connect(clickEventSaveAlmacenamientos)
        self.ButtonSimSetStates.clicked.connect(clickEventSetAlmacenamientos)
        self.ButtonSim_RunSimulacion.clicked.connect(clickEventSimulationDeCuenca)
        self.ButtonSim_ViewStreamflow.clicked.connect(clickEventViewSerieQobsQsim)
        self.ButtonSim_ViewCDC.clicked.connect(clickEventViewCDCQobsQsim)
        self.ButtonSim_ViewSediments.clicked.connect(clickEventViewSerieSedimentos)
        
  
    def setupUIInputsOutputs (self):
        
        def DrawClickEvent_histogram_WMF():
            '''Hace un histograma de la variable seleccionada'''
            #Selecciona el item y su nombre
            selectedItems = self.Tabla_Prop_WMF.currentRow ()
            ItemName =  self.Tabla_Prop_WMF.item(selectedItems,0).text()
            if ItemName[-1] == '*': ItemName = ItemName[:-1]
            #Obtiene la variable e invoca la funciond e grafica
            Var = self.HSutils.DicBasinWMF[ItemName]['var']
            PathFigure = '/tmp/HydroSED/Plots_Geomorfo/VarHistogram.html'
            self.GeoPlots.VarHistogram(Var, PathFigure, self.HSutils.cuenca.ncells)
            #Set de la ventana que contiene la figura.
            self.VistaRainWeb = QWebView(None)
            self.VistaRainWeb.load(QUrl.fromLocalFile(PathFigure))
            self.VistaRainWeb.setWindowTitle('Histograma variable')
            self.VistaRainWeb.setMinimumWidth(100)
            self.VistaRainWeb.setMaximumWidth(500)
            self.VistaRainWeb.setMinimumHeight(100)
            self.VistaRainWeb.setMaximumHeight(500)
            self.VistaRainWeb.show()
            
        def DrawClickEvent_histogram_NC():
            '''Hace un histograma de la variable seleccionada y su capa seleccionada'''
            #Selecciona el item y su nombre
            selectedItems = self.Tabla_Prop_NC.currentRow ()
            ItemName =  self.Tabla_Prop_NC.item(selectedItems,0).text()
            if ItemName[-1] == '*': ItemName = ItemName[:-1]
            #Selecciona la capa 
            if self.SpinBoxNCLayer.value() == 0:
                Var = self.HSutils.DicBasinNc[ItemName]['var']
            else:
                Capa = self.SpinBoxNCLayer.value() - 1
                Var = self.HSutils.DicBasinNc[ItemName]['var'][Capa]
            #Obtiene la variable e invoca la funciond e grafica
            PathFigure = '/tmp/HydroSED/Plots_Geomorfo/VarHistogram.html'
            self.GeoPlots.VarHistogram(Var, PathFigure, self.HSutils.cuenca.ncells)
            #Set de la ventana que contiene la figura.
            self.VistaRainWeb = QWebView(None)
            self.VistaRainWeb.load(QUrl.fromLocalFile(PathFigure))
            self.VistaRainWeb.setWindowTitle('Histograma variable')
            self.VistaRainWeb.setMinimumWidth(100)
            self.VistaRainWeb.setMaximumWidth(500)
            self.VistaRainWeb.setMinimumHeight(100)
            self.VistaRainWeb.setMaximumHeight(500)
            self.VistaRainWeb.show()
            #Vuelve el spinBox a su set original 
            self.SpinBoxNCLayer.setValue(0)
        
        
        def handleClickEventButton_Eliminar_Desde_WMF ():
            #Selecciona el item y su nombre
            selectedItems = self.Tabla_Prop_WMF.currentRow ()
            ItemName =  self.Tabla_Prop_WMF.item(selectedItems,0).text()
            #Remueve de la tabla visible y de los demas elementos.
            self.Tabla_Prop_WMF.removeRow (selectedItems)
            self.TabWMF.DelEntry(ItemName)
            self.HSutils.DicBasinWMF.pop(ItemName)
            
        def handleClickEventButton_Eliminar_Desde_NC ():
            #Selecciona el item y su nombre
            selectedItems = self.Tabla_Prop_NC.currentRow ()
            ItemName = str(self.Tabla_Prop_NC.item(selectedItems,0).text())
            UnsavedInTable = False
            if ItemName[-1] == '*': 
                ItemName = ItemName[:-1]
                UnsavedInTable = True
            #Revisa que todavia no este guardado 
            if UnsavedInTable and self.HSutils.DicBasinNc[ItemName]['saved'] is False or self.HSutils.DicBasinNc[ItemName]['basica'] is False:
                #ItemName = ItemName[:-1]
                #Remueve de la tabla visible y de los demas elementos.
                self.Tabla_Prop_NC.removeRow (selectedItems)
                self.TabNC.DelEntry(ItemName)
                self.HSutils.DicBasinNc.pop(ItemName)
                print self.HSutils.Nc2Save
                self.HSutils.Nc2Save.remove(ItemName)
                print self.HSutils.Nc2Save
                #self.HSutils.Nc2Erase.append(ItemName)
                #Mensaje de exito 
                self.iface.messageBar().pushInfo (u'Hydro-SIG:', u'La variable '+ItemName + ' ha sido borrada')
            else:
                #Mensaje de no exito 
                self.iface.messageBar().pushMessage (u'Hydro-SIG:', u'La variable '+ItemName+' no puede ser borrada del Nc (ya esta guardada)',
                    level=QgsMessageBar.WARNING, duration=3)

        def handleClickEventButton_Actualizar_WMF_Desde_NC ():
            rows = sorted (set (index.row () for index in self.Tabla_Prop_NC.selectedIndexes ()))
            for row in rows:
                print('Row %d is selected' % row)
    
        def handleClickEventButton_Ver_Desde_NC():
            '''Visualiza una de las variables de la cuenca en Qgis'''
            #Nombre de la variable a observar
            selectedItems = self.Tabla_Prop_NC.currentRow ()
            VarName = self.Tabla_Prop_NC.item(selectedItems,0).text()
            if VarName[-1] == '*': VarName = VarName[:-1]
            #Capa de la variable 
            if self.SpinBoxNCLayer.value() == 0:
                Capa = None
            else:
                Capa = self.SpinBoxNCLayer.value() - 1
            #Ejecuta la conversion a mapa raster
            pathMapa = self.HSutils.Basin_LoadVariableFromDicNC(VarName, capa = Capa)
            #Visualiza 
            flagCargaMapa = self.HSutils.cargar_mapa_raster(pathMapa)
            if flagCargaMapa:
                self.iface.messageBar().pushInfo (u'Hydro-SIG:', u'Se cargó la variable de forma exitosa')
            else:
                self.iface.messageBar().pushMessage (u'Hydro-SIG:', u'No fue posible cargar la variable',
                    level=QgsMessageBar.WARNING, duration=3)
            #Hace que el spinbox vuelva a lo normal 
            self.SpinBoxNCLayer.setValue(0)
        
        def handleClickEventButton_Ver_Desde_WMF():
            '''Visualiza una de las variables de la cuenca en Qgis'''
            #Ejecucion de la transformacion de la variable cuenca a raster
            selectedItems = self.Tabla_Prop_WMF.currentRow ()
            VarName = self.Tabla_Prop_WMF.item(selectedItems,0).text()
            pathMapa = self.HSutils.Basin_LoadVariableFromDicWMF(VarName)
            #Visualiza 
            flagCargaMapa = self.HSutils.cargar_mapa_raster(pathMapa)
            if flagCargaMapa:
                self.iface.messageBar().pushInfo (u'Hydro-SIG:', u'Se cargó la variable de forma exitosa')
            else:
                self.iface.messageBar().pushMessage (u'Hydro-SIG:', u'No fue posible cargar la variable',
                    level=QgsMessageBar.WARNING, duration=3)
        
        def handleClickEventButton_NC2WMF():
            '''Mueve variables de NC a WMF en la tabla.'''
            # Elemento seleccionado
            selectedItems = self.Tabla_Prop_NC.currentRow ()
            VarName = self.Tabla_Prop_NC.item(selectedItems,0).text()
            #Verifica que la variable ahun no este guardada en el Nc             
            print VarName.strip()
            if VarName[-1] == '*':
                # Copia la entrada a WMF y la saca de NC
                VarName = VarName[:-1]
                self.HSutils.DicBasinWMF.update({VarName:self.HSutils.DicBasinNc[VarName]})
                # Mete la entrada en la tabla de WMF y la saca de la tabla de NC
                self.TabWMF.NewEntry(self.HSutils.DicBasinNc[VarName], VarName, self.Tabla_Prop_WMF)
                self.Tabla_Prop_NC.removeRow (selectedItems)
                self.TabNC.DelEntry(VarName)
                self.HSutils.DicBasinNc.pop(VarName)
                self.HSutils.Nc2Save.remove(VarName)
                #Mensaje de exito 
                self.iface.messageBar().pushMessage (u'Hydro-SIG:', u'La variable '+VarName+' ha sido movida de NC a WMF')
            else:
                #Mensaje de no exito
                self.iface.messageBar().pushMessage (u'Hydro-SIG:', 
                    u'La variable '+VarName+' no ha podido ser movida a WMF, ya se encuentra guardada en NC', 
                    level=QgsMessageBar.WARNING, duration=3)
        
        def handleClickEventButton_WMF2NC():
            '''Mueve variables de NC a WMF en la tabla.'''
            # Elemento seleccionado
            selectedItems = self.Tabla_Prop_WMF.currentRow ()
            VarName = self.Tabla_Prop_WMF.item(selectedItems,0).text()
            # Copia la entrada a WMF y la saca de NC
            self.HSutils.DicBasinNc.update({VarName:self.HSutils.DicBasinWMF[VarName]})
            self.HSutils.DicBasinNc[VarName]['var'] = np.copy(self.HSutils.DicBasinWMF[VarName]['var'])
            # Mete la entrada en la tabla de WMF y la saca de la tabla de NC
            self.TabNC.NewEntry(self.HSutils.DicBasinWMF[VarName], VarName, self.Tabla_Prop_NC, New = True)
            self.Tabla_Prop_WMF.removeRow (selectedItems)
            self.TabWMF.DelEntry(VarName)
            self.HSutils.DicBasinWMF.pop(VarName)
            self.HSutils.Nc2Save.append(VarName)
            # Mensaje de exito
            self.iface.messageBar().pushInfo (u'Hydro-SIG:', u'La variable '+VarName+' ha sido movida de la Tabla WMF a NC') 
        
        def clickEventBasinUpdateNC():
            '''Actualiza el archivo .nc de la cuenca con las variables cargadas en la TablaNC'''
            RutaNC = self.lineEditRutaCuenca.text().strip()
            self.HSutils.Basin_Update(RutaNC)
            self.TabNC.SavedEntry(self.Tabla_Prop_NC)
        
        def setupLineEditButtonOpenShapeFileDialog (lineEditHolder, fileDialogHolder):
            '''Hace que cuando se busquen shapes solo se encuetren formatos vectoriales'''
            lineEditHolder.setText (fileDialogHolder.getOpenFileName (QtGui.QDialog (), "", "*", "Shapefiles (*.shp);;"))
            if ((os.path.exists (lineEditHolder.text ().strip ())) and (not (self.iface is None))):
                self.iface.addVectorLayer (lineEditHolder.text ().strip (), os.path.basename (lineEditHolder.text ()).strip (), "ogr")

        def setupLineEditButtonOpenRasterFileDialog (lineEditHolder, fileDialogHolder):
            '''Hace que solo0 se busquen formatos aceptados por GDAL'''
            lineEditHolder.setText (fileDialogHolder.getOpenFileName (QtGui.QDialog (), 'Open File',"", GdalTools_utils.FileFilter.allRastersFilter (),
                QtGui.QFileDialog.DontUseNativeDialog))

        def setupLineEditButtonOpenFileDialog (lineEditHolder, fileDialogHolder):
            '''Pone la ruta elegida en el dialogo de texto para cargar'''
            lineEditHolder.setText (fileDialogHolder.getOpenFileName ())

        def setupLineEditButtonSaveFileDialog (lineEditHolder, fileDialogHolder):
            '''Pone la ruta elegida en el dialogo de texto para guardado'''
            lineEditHolder.setText (fileDialogHolder.getSaveFileName ())

        def clickEventSelectorMapaDEM ():
            '''Evento de click: selecciona mapa DEM'''
            setupLineEditButtonOpenRasterFileDialog (self.lineEditMapaDEM, QFileDialog)

        def clickEventSelectorMapaDIR ():

#            setupLineEditButtonOpenFileDialog (self.lineEditMapaDIR, QFileDialog)
            setupLineEditButtonOpenRasterFileDialog (self.lineEditMapaDIR, QFileDialog)

        def clickEventVisualizarMapaDEM ():
            pathMapaDEM = self.lineEditMapaDEM.text ().strip ()
            flagCargaMapaDEM = self.HSutils.cargar_mapa_raster (pathMapaDEM)
            if flagCargaMapaDEM:
                self.iface.messageBar().pushInfo (u'Hydro-SED', u'Se cargó el mapa MDE de forma exitosa')
            else:
                self.iface.messageBar().pushMessage (u'Hydro-SED', 
                    u'No fue posible cargar el mapa MDE. Verifique su ruta. Verifique su formato. Y por favor intente de nuevo.',
                    level=QgsMessageBar.WARNING, duration=3)
            
        def clickEventVisualizarMapaDIR ():
            pathMapaDIR = self.lineEditMapaDIR.text ().strip ()
            flagCargaMapaDIR = self.HSutils.cargar_mapa_raster (pathMapaDIR)
            if flagCargaMapaDIR:
                self.iface.messageBar().pushInfo (u'Hydro-SED', u'Se cargó el mapa DIR de forma exitosa')
            else:
                self.iface.messageBar().pushMessage (u'Hydro-SED', 
                    u'No fue posible cargar el mapa DIR. Verifique su ruta. Verifique su formato. Y por favor intente de nuevo.',
                    level=QgsMessageBar.WARNING, duration=3)

        def clickEventCargarWMFMapaDEM ():
            '''Carga el mapa dDM base para WMF'''
            pathMapaDEM = self.lineEditMapaDEM.text ().strip ()
            dxpMapaDEM  = self.spinBox_dxPlano.value()
            flagCargaMapaDEM_WMF, self.EPSG, self.noData = self.HSutils.cargar_mapa_dem_wmf (pathMapaDEM, dxpMapaDEM)
            if flagCargaMapaDEM_WMF:
                self.iface.messageBar().pushInfo (u'Hydro-SED', u'Se cargó el mapa MDE al WMF de forma exitosa')
            else:
                self.iface.messageBar().pushMessage (u'Hydro-SED', 
                    u'No fue posible cargar el mapa MDE al WMF. Verifique su ruta. Verifique su formato. Y por favor intente de nuevo.',
                    level=QgsMessageBar.WARNING, duration=3)
            #Pone el nombre del codigo EPSG en el dialogo.
            self.LineEditEPSG.setText(self.EPSG)
            t = '%.1f' % self.noData
            self.LineEditNoData.setText(t)
            
        def clickEventCargarWMFMapaDIR():
            '''Carga el mapa DIR base para WMF'''
            pathMapaDIR = self.lineEditMapaDIR.text ().strip ()
            dxpMapaDIR  = self.spinBox_dxPlano.value ()
            flagCargaMapaDIR_WMF, self.EPSG = self.HSutils.cargar_mapa_dir_wmf (pathMapaDIR, dxpMapaDIR)
            if flagCargaMapaDIR_WMF:
                self.iface.messageBar().pushInfo (u'Hydro-SED', u'Se cargó el mapa DIR al WMF de forma exitosa')
            else:
                self.iface.messageBar().pushMessage (u'Hydro-SED', 
                    u'No fue posible cargar el mapa DIR al WMF. Verifique su ruta. Verifique su formato. Y por favor intente de nuevo.',
                    level=QgsMessageBar.WARNING, duration=3)

        def clickEventSelectorBinarioNC ():
            setupLineEditButtonOpenFileDialog (self.lineEditSelectorBinarioNC, QFileDialog)

        def clickEventSelectorOutputCorrienteShapefileTrazadorCorrientes():
            setupLineEditButtonSaveFileDialog(self.PathCorriente_out, QFileDialog)

        def clickEventSelectorInputCuencaShapefileTrazadorCuencas():
            #Hace set del path para la divisoria
            setupLineEditButtonSaveFileDialog(self.PathOutputDivisoria, QFileDialog)
            #Obtiene solo el path y el nombre
            PathOnly = os.path.dirname(self.PathOutputDivisoria.text().strip())
            #PathOnly = self.PathOutputDivisoria
            NameOnly = os.path.basename(self.PathOutputDivisoria.text().strip())
            NameOnly = os.path.splitext(NameOnly)[0]
            #Genera las rutas de la red y del nc
            PathRed = PathOnly +'/' +NameOnly + '_Red.shp'
            PathNC = PathOnly + '/'+NameOnly + 'Topo.nc'
            #Hace un set del path de la red y del nc.
            self.PathOutputRed.setText(PathRed)
            self.PathOutputNETCDF.setText(PathNC)
        
        def clickEventSelectorInputRedShapefileTrazadorCuencas():
            #Hace set del path para la red
            setupLineEditButtonSaveFileDialog(self.PathOutputRed, QFileDialog)
        
        def clickEventSelectorInputNC_ShapefileTrazadorCuencas():
            #Hace set del path para el nc
            setupLineEditButtonSaveFileDialog(self.PathOutputNETCDF, QFileDialog)

        def clickEventSelectorOutputCuencaShapefileTrazadorCuencas ():
            setupLineEditButtonSaveFileDialog (self.lineEditOutputCuencaShapefileTrazadorCuencas, QFileDialog)

        def clickEventSelectorOutputCuencaNCTrazadorCuencas ():
            setupLineEditButtonSaveFileDialog (self.lineEditOutputCuencaNCTrazadorCuencas, QFileDialog)

        self.botonSelectorMapaDEM.clicked.connect (clickEventSelectorMapaDEM)
        self.botonSelectorMapaDIR.clicked.connect (clickEventSelectorMapaDIR)
        
        #Botones para cargar y visualizar mapas DEM y DIR
        self.Boton_verMDE.clicked.connect (clickEventVisualizarMapaDEM)
        self.Boton_verDIR.clicked.connect (clickEventVisualizarMapaDIR)
        self.Boton_MDE2WMF.clicked.connect (clickEventCargarWMFMapaDEM)
        self.Boton_DIR2WMF.clicked.connect (clickEventCargarWMFMapaDIR)

        #self.botonSelectorBinarioNC.clicked.connect (clickEventSelectorBinarioNC)
        
        #Botones para establecer la ruta de guardado del terazador de corrientes yd e cuencas
        self.botonSelectorOutputCorrienteShapefileTrazadorCorrientes.clicked.connect (clickEventSelectorOutputCorrienteShapefileTrazadorCorrientes)
        self.BotonPathDivisoria.clicked.connect(clickEventSelectorInputCuencaShapefileTrazadorCuencas)
        #Botones para ruta de red y de nc.
        self.BotonPathRed.clicked.connect(clickEventSelectorInputRedShapefileTrazadorCuencas)
        self.BotonPathNC.clicked.connect(clickEventSelectorInputNC_ShapefileTrazadorCuencas)

        #Botones para agarrar coordenadqas de trazado
        self.BotonCoord_corriente.clicked.connect(self.handleClickCoordCorrientes)
        self.BotonCoord_cuenca.clicked.connect(self.handleClickCoordCuencas)
        
        #Botones ejecucion trazadores
        self.botonEjecutarTrazadorCorrientes.clicked.connect (self.handleClickEventEjecutarTrazadorCorrientes)
        self.botonEjecutarTrazadorCuencas.clicked.connect(self.handleClickEventEjecutarTrazadorCuencas)
        #self.botonOutputCuencaShapefileTrazadorCuencas.clicked.connect (clickEventSelectorOutputCuencaShapefileTrazadorCuencas)
        #self.botonOutputCuencaNCTrazadorCuencas.clicked.connect (clickEventSelectorOutputCuencaNCTrazadorCuencas)
        #self.spinBox_dxPlano.setValue(self.spinBox_dxPlano.valueFromText())
        def set_dxplano():
            #self.spinBox_dxPlano.valueFromText()
            self.spinBox_dxPlano.value()
        self.spinBox_dxPlano.valueChanged.connect(set_dxplano)
        
        #Botones de borrado de variables de NC y WMF 
        self.Tabla_Prop_WMF.setSelectionMode(QtGui.QAbstractItemView.SingleSelection)
        self.Tabla_Prop_WMF.setSelectionBehavior(QtGui.QAbstractItemView.SelectRows)
        self.Button_Eliminar_Desde_WMF.clicked.connect(handleClickEventButton_Eliminar_Desde_WMF)
        self.Tabla_Prop_NC.setSelectionMode(QtGui.QAbstractItemView.SingleSelection)
        self.Tabla_Prop_NC.setSelectionBehavior(QtGui.QAbstractItemView.SelectRows)        
        self.Button_Eliminar_Desde_NC.clicked.connect(handleClickEventButton_Eliminar_Desde_NC)
        #Botones de visualizacion de variables de NC
        self.Tabla_Prop_NC.setSelectionMode(QtGui.QAbstractItemView.SingleSelection)
        self.Tabla_Prop_NC.setSelectionBehavior(QtGui.QAbstractItemView.SelectRows)        
        self.Button_Visualizar_Desde_NC.clicked.connect(handleClickEventButton_Ver_Desde_NC)
        #Botones de visualizacion de variables de WMF
        self.Tabla_Prop_WMF.setSelectionMode(QtGui.QAbstractItemView.SingleSelection)
        self.Tabla_Prop_WMF.setSelectionBehavior(QtGui.QAbstractItemView.SelectRows)        
        self.Button_Visualizar_Desde_WMF.clicked.connect(handleClickEventButton_Ver_Desde_WMF)
        #Botones movimiento variables NC a WMF y de WMF a NC
        self.Button_NC2WMF.clicked.connect(handleClickEventButton_NC2WMF)
        self.Button_WMF2NC.clicked.connect(handleClickEventButton_WMF2NC)
        #Boton para actualizar los archivos que se encuentran guardados en un netCDF
        self.Button_Update_NC.clicked.connect(clickEventBasinUpdateNC)
        #Botones para graficas de variables de la cuenca
        self.ButtonDraw_histogram_WMF.clicked.connect(DrawClickEvent_histogram_WMF)
        self.ButtonDraw_histogram_NC.clicked.connect(DrawClickEvent_histogram_NC)
 
class Tabla():
    
    def __init__(self, NumRows, TabElement):
        self.NumRows = 0
        self.TabNames = []
        Header = ["Nombre", "Tipo", "Forma", "Categoria"]
        TabElement.setRowCount(NumRows)
        TabElement.setColumnCount(len(Header))
        TabElement.setHorizontalHeaderLabels(Header)
    
    def EmptyTable(self, TabElement):
        '''Vacia la tabla, quita todos los elementos'''
        TabElement.setRowCount(0)
        self.TabNames = []
        self.NumRows = 0 
        
    def DelEntry(self, KeyToDel):
        '''Borra una entrada en el diccionario de datos'''
        pos = self.TabNames.index(KeyToDel)
        self.TabNames.pop(pos)
        self.NumRows -= 1

    def EditedEntry(self, KeyEdited, TabElement, New_or_Edited = 'Edited'):
        '''Establece en la tabla visual si una entrada ha sido editada'''
        #Encuentra la pos
        pos = self.TabNames.index(KeyEdited)
        #Edita el nombre en la tabla 
        TabElement.setItem(pos, 0, QTableWidgetItem(KeyEdited+'*'))
    
    def SavedEntry(self, TabElement):
        '''Busca los elementos de la tabla que terminen con * y se los quita, solo para
        que el usuario sepa que han sido guardados'''
        #Busca en cada entrada
        for i in range(self.NumRows):
            Nombre = TabElement.takeItem(i,0).text()
            if Nombre[-1] == '*':
                TabElement.setItem(i,0, QTableWidgetItem(Nombre[:-1]))
            else:
                TabElement.setItem(i,0, QTableWidgetItem(Nombre))
    
    def NewEntry(self, Dic, DicKey,TabElement, New = False):
        '''Actualiza la lista de las variables en una tabla'''
        #Busca si ese nombre ya se encuentra en la tabla
        try:
            #Si esta, remplaza esa posicion
            pos = self.TabNames.index(DicKey)
            suma = 0
        except:
            #Si no esta, lo pone al final.
            pos = self.NumRows
            self.TabNames.append(DicKey)
            suma = 1
        #for keyParam in Dic:
        if New:
            TabElement.setItem (pos, 0, QTableWidgetItem (Dic["nombre"]+'*'))
        else:
            TabElement.setItem (pos, 0, QTableWidgetItem (Dic["nombre"]))
        TabElement.setItem (pos, 1, QTableWidgetItem (Dic["tipo"]))
        TabElement.setItem (pos, 2, QTableWidgetItem (str (Dic["shape"])))
        TabElement.setItem (pos, 3, QTableWidgetItem (Dic["categoria"]))
        self.NumRows += suma
            
    
    

        
