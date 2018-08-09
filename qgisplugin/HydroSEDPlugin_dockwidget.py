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
import numpy as np

from PyQt4 import QtGui, uic, QtCore
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
        self.setupRainfallInterpolation()
        self.setupSimulation()
        
        self.TablaFila_WMF = 0
        self.TablaFila_NC = 0
        self.botonClicked = False
        
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
            self.PathCorriente_out.setText('/tmp/HydroSED/Corriente.shp')
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
            self.PathOutputDivisoria.setText('/tmp/HydroSED/Cuenca.shp')
        if len(self.PathOutputRed.text()) == 0:
            self.PathOutputRed.setText('/tmp/HydroSED/Cuenca_Red.shp')
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
            lineEditHolder = fileDialogHolder.getOpenFileName (QtGui.QDialog (), "Guardar en capa vectorial...", "*", "Shapefiles (*.shp);;")
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
            #Cargado de la cuenca
            self.HSutils.Basin_LoadBasin(self.lineEditRutaCuenca.text().strip(), Simhidro, SimSed)
            self.TableStart()
            for k in self.HSutils.DicBasinNc.keys():
                self.TabNC.NewEntry(self.HSutils.DicBasinNc[k],k, self.Tabla_Prop_NC)
            Area, self.EPSG, dxp, self.noData, self.umbral = self.HSutils.Basin_LoadBasin(self.lineEditRutaCuenca.text().strip())
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
            
        def clickEventBasinLoadDivisory():
            '''Carga la divisoria de la cuenca cargada a WMF'''
            OutPathDivisoria = '/tmp/HydroSED/CuencaCargada.shp'
            self.HSutils.Basin_LoadBasinDivisory(OutPathDivisoria)
            #Carga la divisoria
            ret, layer = self.HSutils.cargar_mapa_vector(OutPathDivisoria, self.HSutils.TIPO_STYLE_POLIGONO, color = (255,0,0), width = 0.6)            
            self.iface.mapCanvas().refresh() 
            self.iface.legendInterface().refreshLayerSymbology(layer)      
        def clickEventBasinLoadNetwork():
            '''Carga la divisoria de la cuenca cargada a WMF'''
            OutPathNetwork = '/tmp/HydroSED/RedCargada.shp'
            self.HSutils.Basin_LoadBasinNetwork(OutPathNetwork)
            #Carga la red 
            ret, layer = self.HSutils.cargar_mapa_vector(OutPathNetwork, self.HSutils.TIPO_STYLE_POLILINEA, width = 0.4)
            self.iface.mapCanvas().refresh() 
            self.iface.legendInterface().refreshLayerSymbology(layer)   
        
        def clickEventBasinVarNC2Network():
            '''Convierte un conjunto de variables a una red hidrica de la cuenca'''
            if self.botonClicked is False:
                #Hace que la seleccion sea multiple en la columna de NC
                self.Tabla_Prop_NC.setSelectionMode(QAbstractItemView.MultiSelection)
                self.ruta2Network = ''
                #imprime las variables seleccionadas
                self.Vars2Network = []
                self.botonClicked = True
                self.ButtonVar2Net_NC.setText('Variables')
                return
            elif self.botonClicked is True:
                for i in self.Tabla_Prop_NC.selectedItems()[::4]:
                    self.Vars2Network.append(i.text())
                #En el segundo click selecciona el archivo y ejecuta 
                self.ruta2Network = SelectNetworkshapeDialog(QFileDialog)
                self.Tabla_Prop_NC.setSelectionMode(QAbstractItemView.SingleSelection)
                self.botonClicked = False
                #Guardado de las variables en formato red hidrica
                self.HSutils.BasinNc2Network(self.ruta2Network,self.Vars2Network)
                self.ButtonVar2Net_NC.setText('Var2Net')
        
        #Botones para variables de entrada 
        self.botonSelectorProyectBasin.clicked.connect(clickEventSelectorBasin)
        self.ButtonLoadBasinProyect.clicked.connect(clickEventBasin2WMF)
        #Botones para visualizar polilineas y poligonos 
        self.Boton_verDivisoria.clicked.connect(clickEventBasinLoadDivisory)
        self.Boton_verRedHidrica.clicked.connect(clickEventBasinLoadNetwork)
        #boton para convertir variables a red hidrica
        self.ButtonVar2Net_NC.clicked.connect(clickEventBasinVarNC2Network)
        

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
                ListaVar.extend(['OCG_coef'])
            if self.checkBoxKubota.isChecked():
                self.HSutils.Basin_GeoGetKubota()
                ListaVar.extend(['kubota_coef'])
            if self.checkBoxRunoff.isChecked():
                self.HSutils.Basin_GeoGetRunoff()
                ListaVar.extend(['Runoff_coef'])
                
            #mensaje de caso de exito
            self.iface.messageBar().pushInfo(u'HidroSIG:',u'Calculo de geomorfologia distribuida realizado, revisar la tabla Variables WMF.')
            #Actualiza la tabla de variables temporales 
            for k in ListaVar:
                self.TabWMF.NewEntry(self.HSutils.DicBasinWMF[k],k, self.Tabla_Prop_WMF)
                self.ComboGeoMaskVar.addItem(k)
                self.ComboGeoVar2Acum.addItem(k)
        
        #Botones de ejecucion
        self.ButtonGeomorfoRasterVars.clicked.connect(clickEventGeoRasterProp)
        self.checkBoxTodos.clicked.connect(clickEventActivateGeoCheckBoxes)
        self.ButtonGeoParameters.clicked.connect(clickEventGeoProperties)
        #self.ComboGeoMaskVar.activated.connect(clickEventUpdateComboBoxMask)
    
    def setupHidro_Balance(self):
        
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
            QSalida = self.HSutils.hidologia_balance(self.spinBox_dxPlano.value(),
                self.spinBoxUmbralRed.value(), 
                self.PathInHydro_Rain.text(), 
                TipoETR, 
                self.PathOutHydro_Qmed.text())
            #Pone el valor de cadual medio en el cuadro
            textoCaudal = '%.3f' % QSalida
            self.ShowResultQmed.setText(textoCaudal)
            #Mensaje de exito 
            self.iface.messageBar().pushInfo(u'HidroSIG:',u'Balance realizado: variables Caudal, ETR y Runoff cargadas a Tabla de propiedades WMF.')
            #Habilita botones de visualizacion de variables 
            if len(self.PathOutHydro_Qmed.text()) > 2:
                self.Button_HidroViewQmed.setEnabled(True)
            #Actualiza la tabla de variables temporales y actualiza comboBox de gomorfo
            for k in ['Caudal','ETR','Runoff']:
                self.TabWMF.NewEntry(self.HSutils.DicBasinWMF[k],k, self.Tabla_Prop_WMF)
                
        #Botones para variables de entrada 
        self.Boton_HidroLoadRain.clicked.connect(clickEventSelectorRaster)
        #Botones para variables de salida
        self.Button_HidroSaveQmed.clicked.connect(clickEventOutQmed)
        #Botones para visualizar variables 
        self.Button_HidroViewRain.clicked.connect(clickEventViewRainfall)
        self.Button_HidroViewQmed.clicked.connect(clickEventViewQmedNetwork)
        #Botones para ejecutar
        self.Butto_Ejec_HidroBalance.clicked.connect(hadleClickEventEjecutarBalance)
        
    def TableStart (self):
        '''Arranca las tablas de NC y WMF'''
        self.TabNC = Tabla(self.HSutils.NumDicBasinNcVariables,self.Tabla_Prop_NC)
        self.TabWMF = Tabla(self.HSutils.NumDicBasinWMFVariables, self.Tabla_Prop_WMF)
    
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
            #if len(self.NameRaster2WMF.text())<2:
             #   self.iface.messageBar().pushError (u'Hydro-SIG:', u'Debe ingresar un nombre para el mapa raster a convertir.')
              #  return 1
            #if len(self.PathRaster2WMF.text())<2:
             #   self.iface.messageBar().pushError (u'Hydro-SIG:', u'Debe seleccionar un mapa raster para ser convertido a la cuenca.')
             #   return 1
            #Parametros para la conversion
            Nombre = self.NameRaster2WMF.text()
            PathRaster = self.PathRaster2WMF.text()
            Grupo = self.ComboBoxRaster2WMF.currentText()
            #Conversion, convierte la variable y actualiza el diccionario.
            Retorno = self.HSutils.Basin_Raster2WMF(Nombre, PathRaster, Grupo)
            self.TabWMF.NewEntry(self.HSutils.DicBasinWMF[Nombre],Nombre, self.Tabla_Prop_WMF)
            if Retorno == 0:
                self.iface.messageBar().pushInfo (u'Hydro-SIG:', u'El mapa raster ha ingresado a WMF.')
        
        #Habilita botones.
        self.ButtonPathRaster2WMF.clicked.connect(clickEventSelectorMapaRaster)
        self.Button_Raster2WMF.clicked.connect(handleClickConnectRaster2WMF)
        
    
    def setupRainfallInterpolation(self):
        '''Conjunto de herramientas dispuestas para interpolar campos de precipitacion'''
        
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
        
        def clickEventViewSerieRainfall():
            '''Genera y visualiza la grafica de lluvia interpolada para la cuenca'''
            #Hace la figura
            PathFigure = '/tmp/HydroSED/RainfallPlot.html'
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
            PathFigure = '/tmp/HydroSED/RainfallHistogram.html'
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
            PathFigure = '/tmp/HydroSED/RainfallMediaMensual.html'
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
        
    
    def setupSimulation(self):
        '''Herramientas para gestionar la simulacion hidrologica con la cuenca cargada'''
        
        def clickEventUpdateParamMapValues():
            '''Muestra en los campos de simulacion el valor medio de los mapas de simulacion'''
            VarNames = ['h1_max','h3_max', 'v_coef','v_coef','v_coef','v_coef','h_coef',
                'h_coef','h_coef','h_coef', 'Krusle','Crusle','Prusle',
                'Arenas','Limos','Arcillas']
            Ejes = [0,0,0,1,2,3,0,1,2,3,0,0,0,0,0,0]
            for name, i, eje in zip(VarNames, range(1,17), Ejes):
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
            ParamName = self.ParamName.text().strip()
            #Itera para los parametros escalares y de sedimentos
            ListaParam = []
            for i in range(1,12):
                ListaParam.append(getattr(self, 'Param'+str(i)).value())
            #Itera en los exponentes
            for i in range(1,5):
                ListaParam.append(getattr(self, 'ParamExp'+str(i)).value())    
            #Mete el set nuevo de calibracion
            self.HSutils.Sim_SaveParameters(PathNC, ParamName, ListaParam)
            
            
        self.ButtonSimCalib2Nc.clicked.connect(clickEventAddNewParamSet)    
        self.tabPanelDockOpciones.currentChanged.connect(clickEventUpdateParamMapValues)
        
    def setupUIInputsOutputs (self):
        
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
            #Remueve de la tabla visible y de los demas elementos.
            self.Tabla_Prop_NC.removeRow (selectedItems)
            self.TabNC.DelEntry(ItemName)
            self.HSutils.DicBasinNc.pop(ItemName)
            self.HSutils.Nc2Erase.append(ItemName)

        def handleClickEventButton_Actualizar_WMF_Desde_NC ():
            rows = sorted (set (index.row () for index in self.Tabla_Prop_NC.selectedIndexes ()))
            for row in rows:
                print('Row %d is selected' % row)
    
        def handleClickEventButton_Ver_Desde_NC():
            '''Visualiza una de las variables de la cuenca en Qgis'''
            #Ejecucion de la transformacion de la variable cuenca a raster
            selectedItems = self.Tabla_Prop_NC.currentRow ()
            VarName = self.Tabla_Prop_NC.item(selectedItems,0).text()
            try:
                pathMapa = self.HSutils.Basin_LoadVariableFromDicNC(VarName)
            except:
                pathMapa = self.HSutils.Basin_LoadVariableFromDicNC(VarName[:-1])
            #Visualiza 
            flagCargaMapa = self.HSutils.cargar_mapa_raster(pathMapa)
            if flagCargaMapa:
                self.iface.messageBar().pushInfo (u'Hydro-SIG:', u'Se cargó la variable de forma exitosa')
            else:
                self.iface.messageBar().pushError (u'Hydro-SIG:', u'No fue posible cargar la variable')
        
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
                self.iface.messageBar().pushError (u'Hydro-SIG:', u'No fue posible cargar la variable')
        
        def handleClickEventButton_NC2WMF():
            '''Mueve variables de NC a WMF en la tabla.'''
            # Elemento seleccionado
            selectedItems = self.Tabla_Prop_NC.currentRow ()
            VarName = self.Tabla_Prop_NC.item(selectedItems,0).text()
            # Copia la entrada a WMF y la saca de NC
            self.HSutils.DicBasinWMF.update({VarName:self.HSutils.DicBasinNc[VarName]})
            # Mete la entrada en la tabla de WMF y la saca de la tabla de NC
            self.TabWMF.NewEntry(self.HSutils.DicBasinNc[VarName], VarName, self.Tabla_Prop_WMF)
            self.Tabla_Prop_NC.removeRow (selectedItems)
            self.TabNC.DelEntry(VarName)
            self.HSutils.DicBasinNc.pop(VarName)            
        
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
        
        def clickEventBasinUpdateNC():
            '''Actualiza el archivo .nc de la cuenca con las variables cargadas en la TablaNC'''
            RutaNC = self.lineEditRutaCuenca.text().strip()
            self.HSutils.Basin_Update(RutaNC)
        
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
                self.iface.messageBar().pushError (u'Hydro-SED', u'No fue posible cargar el mapa MDE. Verifique su ruta. Verifique su formato. Y por favor intente de nuevo.')
            
        def clickEventVisualizarMapaDIR ():
            pathMapaDIR = self.lineEditMapaDIR.text ().strip ()
            flagCargaMapaDIR = self.HSutils.cargar_mapa_raster (pathMapaDIR)
            if flagCargaMapaDIR:
                self.iface.messageBar().pushInfo (u'Hydro-SED', u'Se cargó el mapa DIR de forma exitosa')
            else:
                self.iface.messageBar().pushError (u'Hydro-SED', u'No fue posible cargar el mapa DIR. Verifique su ruta. Verifique su formato. Y por favor intente de nuevo.')

        def clickEventCargarWMFMapaDEM ():
            '''Carga el mapa dDM base para WMF'''
            pathMapaDEM = self.lineEditMapaDEM.text ().strip ()
            dxpMapaDEM  = self.spinBox_dxPlano.value()
            flagCargaMapaDEM_WMF, self.EPSG, self.noData = self.HSutils.cargar_mapa_dem_wmf (pathMapaDEM, dxpMapaDEM)
            if flagCargaMapaDEM_WMF:
                self.iface.messageBar().pushInfo (u'Hydro-SED', u'Se cargó el mapa MDE al WMF de forma exitosa')
            else:
                self.iface.messageBar().pushError (u'Hydro-SED', u'No fue posible cargar el mapa MDE al WMF. Verifique su ruta. Verifique su formato. Y por favor intente de nuevo.')
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
                self.iface.messageBar().pushError (u'Hydro-SED', u'No fue posible cargar el mapa DIR al WMF. Verifique su ruta. Verifique su formato. Y por favor intente de nuevo.')

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
        #Botones de visualizacion de variables de NC
        self.Tabla_Prop_NC.setSelectionMode(QtGui.QAbstractItemView.SingleSelection)
        self.Tabla_Prop_NC.setSelectionBehavior(QtGui.QAbstractItemView.SelectRows)        
        self.Button_Visualizar_Desde_NC.clicked.connect(handleClickEventButton_Ver_Desde_NC)
        #Botones de visualizacion de variables de WMF
        self.Tabla_Prop_WMF.setSelectionMode(QtGui.QAbstractItemView.SingleSelection)
        self.Tabla_Prop_WMF.setSelectionBehavior(QtGui.QAbstractItemView.SelectRows)        
        self.Button_Visualizar_Desde_WMF.clicked.connect(handleClickEventButton_Ver_Desde_WMF)
        #Botones movimiento variables NC a WMF y de WMF a NC
        #self.Button_NC2WMF.clicked.connect(handleClickEventButton_NC2WMF)
        self.Button_WMF2NC.clicked.connect(handleClickEventButton_WMF2NC)
        #Boton para actualizar los archivos que se encuentran guardados en un netCDF
        self.Button_Update_NC.clicked.connect(clickEventBasinUpdateNC)
 
class Tabla():
    
    def __init__(self, NumRows, TabElement):
        self.NumRows = 0
        self.TabNames = []
        Header = ["Nombre", "Tipo", "Forma", "Categoria"]
        TabElement.setRowCount(NumRows)
        TabElement.setColumnCount(len(Header))
        TabElement.setHorizontalHeaderLabels(Header)
    
    def DelEntry(self, KeyToDel):
        '''Borra una entrada en el diccionario de datos'''
        pos = self.TabNames.index(KeyToDel)
        self.TabNames.pop(pos)
        self.NumRows -= 1
    
    def SavedEntry(self, TabElement):
		'''Busca los elementos de la tabla que terminen con * y se los quita, solo para
		que el usuario sepa que han sido guardados'''
		#Busca en cada entrada
		for i in range(self.NumRows):
			Nombre = TabElement.ItemAt(i,0).text()
			if Nombre[-1] == '*':
				TabElement.setItem(i,0, QTableWidgetItem(Nombre[:-1]))
    
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
            
    
    

        
