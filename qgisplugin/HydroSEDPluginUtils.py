import os.path

from qgis.core import QgsRasterLayer, QgsMapLayerRegistry

from wmf import wmf

class controlHS:
	
	def __init__(self):
		self.DEM = 0
		self.DIR = 0
		self.xll = 0
		self.yll = 0
		self.nodata = 0

	def cargar_mapa_raster (self,pathMapaRaster):
	
	    retornoCargaLayerMapaRaster = False
	
	    pathMapaRaster = pathMapaRaster.strip ()
	
	    if (os.path.exists (pathMapaRaster)):
	
	        baseNameMapaRaster   = os.path.basename (pathMapaRaster)
	        layerMapaRaster = QgsRasterLayer (pathMapaRaster, baseNameMapaRaster)
	        QgsMapLayerRegistry.instance ().addMapLayer (layerMapaRaster)
	
	        retornoCargaLayerMapaRaster = layerMapaRaster.isValid ()
	
	    return retornoCargaLayerMapaRaster
	
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
	
