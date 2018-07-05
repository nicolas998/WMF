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
from PyQt4.QtGui import QFileDialog, QTableWidgetItem

import os.path

import HydroSEDPluginUtils as HSutils
import HydroGetCoordinates as HSCoord

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
        
        self.TablaFila_WMF = 0
        self.TablaFila_NC = 0
        
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
            self.HSutils.Basin_LoadBasin(self.lineEditRutaCuenca.text().strip())
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
        
        #Botones para variables de entrada 
        self.botonSelectorProyectBasin.clicked.connect(clickEventSelectorBasin)
        self.ButtonLoadBasinProyect.clicked.connect(clickEventBasin2WMF)
        #Botones para visualizar polilineas y poligonos 
        self.Boton_verDivisoria.clicked.connect(clickEventBasinLoadDivisory)
        self.Boton_verRedHidrica.clicked.connect(clickEventBasinLoadNetwork)

    def setupGeomorfologia(self):
        '''Conjunto de herramientas para manejar parametros geomorfologicos de la cuenca analizada'''

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
            #mensaje de caso de exito
            self.iface.messageBar().pushInfo(u'HidroSIG:',u'Calculo de geomorfologia distribuida realizado, revisar la tabla Variables WMF.')
            #Actualiza la tabla de variables temporales 
            for k in ListaVar:
                self.TabWMF.NewEntry(self.HSutils.DicBasinWMF[k],k, self.Tabla_Prop_WMF)
        #Botones de ejecucion
        self.ButtonGeomorfoRasterVars.clicked.connect(clickEventGeoRasterProp)
    
    
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
            #Actualiza la tabla de variables temporales 
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

        def handleClickEventButton_Actualizar_WMF_Desde_NC ():
            rows = sorted (set (index.row () for index in self.Tabla_Prop_NC.selectedIndexes ()))
            for row in rows:
                print('Row %d is selected' % row)
    
        def handleClickEventButton_Ver_Desde_NC():
            '''Visualiza una de las variables de la cuenca en Qgis'''
            #Ejecucion de la transformacion de la variable cuenca a raster
            selectedItems = self.Tabla_Prop_NC.currentRow ()
            VarName = self.Tabla_Prop_NC.item(selectedItems,0).text()
            pathMapa = self.HSutils.Basin_LoadVariableFromDicNC(VarName)
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
            # Mete la entrada en la tabla de WMF y la saca de la tabla de NC
            self.TabNC.NewEntry(self.HSutils.DicBasinWMF[VarName], VarName, self.Tabla_Prop_NC)
            self.Tabla_Prop_WMF.removeRow (selectedItems)
            self.TabWMF.DelEntry(VarName)
            self.HSutils.DicBasinWMF.pop(VarName)            
                
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
    
    def NewEntry(self, Dic, DicKey,TabElement):
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
        TabElement.setItem (pos, 0, QTableWidgetItem (Dic["nombre"]))
        TabElement.setItem (pos, 1, QTableWidgetItem (Dic["tipo"]))
        TabElement.setItem (pos, 2, QTableWidgetItem (str (Dic["shape"])))
        TabElement.setItem (pos, 3, QTableWidgetItem (Dic["categoria"]))
        self.NumRows += suma
            
    
    

        
