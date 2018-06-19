import os.path

from qgis.core import QgsRasterLayer, QgsMapLayerRegistry, QgsVectorLayer, QgsFillSymbolV2
from PyQt4 import QtGui, uic
import netCDF4
from wmf import wmf

class controlHS:
    
    TIPO_STYLE_POLIGONO  = 1
    TIPO_STYLE_POLILINEA = 2
    
    def __init__(self):
        self.DEM = 0
        self.DIR = 0
        self.xll = 0
        self.yll = 0
        self.nodata = 0
        self.BasinsCount = 0
        self.StreamsCount = 0
        self.DicBasinNc = {}
        self.DicBasinWMF = {}


    def cargar_mapa_raster (self,pathMapaRaster):
    
        retornoCargaLayerMapaRaster = False
    
        pathMapaRaster = pathMapaRaster.strip ()
    
        if (os.path.exists (pathMapaRaster)):
    
            baseNameMapaRaster   = os.path.basename (pathMapaRaster)
            baseNameMapaRaster = os.path.splitext(baseNameMapaRaster)[0]
            layerMapaRaster = QgsRasterLayer (pathMapaRaster, baseNameMapaRaster)
            QgsMapLayerRegistry.instance ().addMapLayer (layerMapaRaster)
    
            retornoCargaLayerMapaRaster = layerMapaRaster.isValid ()
    
        return retornoCargaLayerMapaRaster
    
    def cargar_mapa_vector(self, pathMapaVector, tipo_style, color = (50,50,250), width = 0.5):
        #Inicia vandera de cargado y ruta del vector
        retornoCargarMapaVector = False
        pathMapaVector = pathMapaVector.strip()
        #verifica existencia y dado el caso carga

        if os.path.exists(pathMapaVector):

            baseNameMapaVector = os.path.basename(pathMapaVector)
            baseNameMapaVector = os.path.splitext(baseNameMapaVector)[0]
            layerMapaVector = QgsVectorLayer(pathMapaVector, baseNameMapaVector, 'ogr')
            QgsMapLayerRegistry.instance().addMapLayer(layerMapaVector)

            if tipo_style == self.TIPO_STYLE_POLILINEA:

                symbols = layerMapaVector.rendererV2().symbols()
                symbol = symbols[0]
                symbol.setColor(QtGui.QColor.fromRgb(color[0],color[1],color[2]))
                symbol.setWidth(width)

            #try:
            #    symbol.setWidth(width)
            #except:
            #    symbol.setBorderWidth(width)
            #if layerMapVector.geometryType() == QGis.Polygon:

            elif tipo_style == self.TIPO_STYLE_POLIGONO:

                Render = layerMapaVector.rendererV2()
                mySymbol1 = QgsFillSymbolV2.createSimple({'color':'blue', 
                                                          'color_border':'#%02x%02x%02x' % color,
                                                          'width_border':str(width),
                                                          'style':'no',
                                                          'style_border':'solid'})

                Render.setSymbol(mySymbol1)
                layerMapaVector.triggerRepaint()

            retornoCargarMapaVector = layerMapaVector.isValid()
        return retornoCargarMapaVector, layerMapaVector
    
    def cargar_mapa_dem_wmf (self,pathMapaDEM, dxp):
        retornoCargaLayerMapaRaster = False
        pathMapaDEM = pathMapaDEM.strip ()
        EPSG_code = -999
        try:
            self.DEM, EPSG_code = wmf.read_map_raster (pathMapaDEM, isDEMorDIR = True, dxp = dxp, noDataP = -9999)
            retornoCargaLayerMapaRaster = True
        except:
            retornoCargaLayerMapaRaster = False
        return retornoCargaLayerMapaRaster, EPSG_code, wmf.cu.nodata
    
    def cargar_mapa_dir_wmf (self,pathMapaDIR, dxp):
        retornoCargaLayerMapaRaster = False
        pathMapaDIR = pathMapaDIR.strip ()
        EPSG_code = -999
        try:
            self.DIR, EPSG_code = wmf.read_map_raster (pathMapaDIR, isDEMorDIR = True, isDIR = True, dxp = dxp, noDataP = -9999)
            retornoCargaLayerMapaRaster = True    
        except:
            self.DIR = 1
        return retornoCargaLayerMapaRaster, EPSG_code
    
    def trazador_corriente(self,x,y, path = None):
        #Traza la corriente en las coordenadas especificadas
        self.stream = wmf.Stream(x, y, self.DEM, self.DIR)
        #Guarda la corriente.
        if path is not None:
            self.stream.Save_Stream2Map(path)
        
    def trazador_cuenca(self,x,y, dxp,umbral,PathDiv, PathRed, PathNC,PathDEM, PathDIR,TopoNodes = False, LastStream = True):
        # Traza la cuenca con y sin la ultima corriente.
        if LastStream:
            self.cuenca = wmf.SimuBasin(x, y, self.DEM, self.DIR, stream=self.stream)
        else:
            self.cuenca = wmf.SimuBasin(x, y, self.DEM, self.DIR)
        # Guarda los shapes de divisoria y de red hidrica.
        self.cuenca.Save_Basin2Map(PathDiv, dxp)
        self.cuenca.Save_Net2Map(PathRed, dxp, umbral)
        # Guarda el nc de la cuenca 
        if len(PathNC)>2:
            self.cuenca.Save_SimuBasin(PathNC)
        
    def hidologia_balance(self, dxp, umbral, PathRain, PathETR, PathQmed):
        #Se fija si la lluvia es un path o un valor 
        try:
            Rain = float(PathRain)
        except:
            Rain, prop = wmf.read_map_raster(PathRain) 
            Rain = self.cuenca.Transform_Map2Basin(Rain, prop)
        #Realiza el balance 
        self.cuenca.GetQ_Balance(Rain, Tipo_ETR = PathETR)
        #Actualiza el diccionario de WMF 
        self.DicBasinWMF.update({'Caudal':
            {'nombre':'Caudal',
            'tipo':'float32',
            'shape':self.cuenca.CellQmed.shape,
            'raster':True,
            'basica': False,
            'categoria': 'Hidro',
            'var': self.cuenca.CellQmed}})
        self.DicBasinWMF.update({'ETR':
            {'nombre':'ETR',
            'tipo':'float32',
            'shape':self.cuenca.CellETR.shape,
            'raster':True,
            'basica': False,
            'categoria': 'Hidro',
            'var': self.cuenca.CellETR}})
        Runoff = Rain - self.cuenca.CellETR
        self.DicBasinWMF.update({'Runoff':
            {'nombre':'Runoff',
            'tipo':'float32',
            'shape':self.cuenca.CellETR.shape,
            'raster':True,
            'basica': False,
            'categoria': 'Hidro',
            'var': Runoff}})
        # Guarda el resultado 
        if len(PathQmed)>2:
            self.cuenca.Save_Net2Map(PathQmed, dxp, umbral, qmed = self.cuenca.CellQmed)
        #Retorna el resultado a la salida 
        return self.cuenca.CellQmed[-1]

    def Basin_LoadBasin(self, PathNC):
        # Numero Total de Variables
        self.NumDicBasinNcVariables = 0
        # Numero Total de Variables Basicas
        self.NumDicBasinNcVariablesBasicas = 0
        #Cargar la cuenca y sus variables base a WMF 
        self.cuenca = wmf.SimuBasin(rute = PathNC)
        #Cargar las variables de la cuenca a un diccionario.
        g = netCDF4.Dataset(PathNC)
        for k in g.variables.keys():
            #Evalua si tiene la misma cantidad de celdas y puede ser un mapa
            shape = g.variables[k].shape
            MapaRaster = False
            for s in shape:
                if s == self.cuenca.ncells:
                    MapaRaster = True
                #Actualiza el diccionario
                self.DicBasinNc.update({k:
                    {'nombre':k,
                    'tipo':g.variables[k].dtype.name,
                    'shape':g.variables[k].shape,
                    'raster':MapaRaster,
                    'basica': True,
                    'categoria': 'Base'}})
                self.NumDicBasinNcVariables = self.NumDicBasinNcVariables + 1
                self.NumDicBasinNcVariablesBasicas = self.NumDicBasinNcVariablesBasicas + 1
        g.close()
        #Cargar la cuenca y sus variables base a WMF 
        self.cuenca = wmf.SimuBasin(rute = PathNC)
        #Area de la cuenca y codigo EPSG  
        return self.cuenca.ncells*wmf.cu.dxp**2./1e6, self.cuenca.epsg, wmf.models.dxp, wmf.cu.nodata
    
    def Basin_LoadBasinDivisory(self, PathDivisory):
        # Guarda los shapes de divisoria y de red hidrica.
        self.cuenca.Save_Basin2Map(PathDivisory, wmf.cu.dxp)
    
    def Basin_LoadBasinNetwork(self, PathNetwork):
        # Guarda los shapes de divisoria y de red hidrica.
        self.cuenca.Save_Net2Map(PathNetwork, wmf.cu.dxp, self.cuenca.umbral)
    
    def Basin_LoadBasicVariable(self,PathNC, VarName):
        '''Toma una variable del diccionario de variables basicas y la carga a Qgis'''
        #Lee los datos del nc de la cuenca de la variable indicada 
        g = netCDF4.Dataset(PathNC)
        Data = g.variables[VarName][:]
        g.close()
        #Transforma a un raster 
        rutaSalida = '/tmp/HydroSED/Raster_'+VarName+'.tiff'
        self.cuenca.Transform_Basin2Map(Data,
            ruta = rutaSalida,
            EPSG = self.cuenca.epsg)
        return rutaSalida
        
    def Basin_GeoGetHAND(self, umbral):
        self.cuenca.GetGeo_HAND(umbral)
        self.DicBasinWMF.update({'HAND':
            {'nombre':'HAND',
            'tipo':'float32',
            'shape':self.cuenca.CellHAND.shape,
            'raster':True,
            'basica': False,
            'categoria': 'Geo',
            'var': self.cuenca.CellHAND}})
        self.DicBasinWMF.update({'HDND':
            {'nombre':'HDND',
            'tipo':'float32',
            'shape':self.cuenca.CellHDND.shape,
            'raster':True,
            'basica': False,
            'categoria': 'Geo',
            'var': self.cuenca.CellHDND}})
        self.DicBasinWMF.update({'HAND_class':
            {'nombre':'HAND_class',
            'tipo':'float32',
            'shape':self.cuenca.CellHAND_class.shape,
            'raster':True,
            'basica': False,
            'categoria': 'Geo',
            'var': self.cuenca.CellHAND_class}})
