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
    
        try:
    
            self.DEM = wmf.read_map_raster (pathMapaDEM, isDEMorDIR = True, dxp = dxp, noDataP = -9999)
            retornoCargaLayerMapaRaster = True
    
        except:
    
            retornoCargaLayerMapaRaster = False
    
        return retornoCargaLayerMapaRaster
    
    
    def cargar_mapa_dir_wmf (self,pathMapaDIR, dxp):
        retornoCargaLayerMapaRaster = False
        pathMapaDIR = pathMapaDIR.strip ()
        try:
            self.DIR = wmf.read_map_raster (pathMapaDIR, isDEMorDIR = True, isDIR = True, dxp = dxp, noDataP = -9999)
            retornoCargaLayerMapaRaster = True    
        except:
            self.DIR = 1
        return retornoCargaLayerMapaRaster
    
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
        
    def hidologia_balance(self, dxp, umbral, PathRain, PathETR, PathQmed, PathETROUT, PathRunoff):
        #Se fija si la lluvia es un path o un valor 
        try:
            Rain = float(PathRain)
        except:
            Rain, prop = wmf.read_map_raster(PathRain) 
            Rain = self.cuenca.Transform_Map2Basin(Rain, prop)
        #Realiza el balance 
        self.cuenca.GetQ_Balance(Rain, Tipo_ETR = PathETR)
        # Guarda el resultado 
        if len(PathQmed)>2:
            self.cuenca.Save_Net2Map(PathQmed, dxp, umbral, qmed = self.cuenca.CellQmed)
        if len(PathETROUT)>2:
            self.cuenca.Transform_Basin2Map(self.cuenca.CellETR, PathETROUT)
        if len(PathRunoff)>2:
            Runoff = Rain - self.cuenca.CellETR
            self.cuenca.Transform_Basin2Map(Runoff, PathRunoff)
                
        #Retorna el resultado a la salida 
        return self.cuenca.CellQmed[-1]

    def Basin_LoadBasin(self, PathNC):
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
				'raster':MapaRaster}})
		g.close()
		print self.DicBasinNc.keys()
	
