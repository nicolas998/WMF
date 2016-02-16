from cuencas import *
import numpy as np
import glob
from mpl_toolkits.basemap import Basemap, addcyclic, shiftgrid, cm
import pylab as pl
import scipy.ndimage as nd
from netCDF4 import Dataset  
import osgeo.ogr, osgeo.osr
import gdal
from scipy.spatial import Delaunay

#-----------------------------------------------------------------------
#Lectura de informacion y mapas 
#-----------------------------------------------------------------------
def read_map_raster(ruta_map,isDEMorDIR=False):
	'Funcion: read_map\n'\
	'Descripcion: Lee un mapa raster soportado por GDAL.\n'\
	'Parametros Obligatorios:.\n'\
	'	-ruta_map: Ruta donde se encuentra el mpaa.\n'\
	'Parametros Opcionales:.\n'\
	'	-isDEMorDIR: Pasa las propiedades de los mapas al modulo cuencas \n'\
	'		escrito en fortran \n'\
	'Retorno:.\n'\
	'	Si no es DEM o DIR retorna todas las propieades del elemento en un vector.\n'\
	'		En el siguiente orden: ncols,nrows,xll,yll,dx,nodata.\n'\
	'	Si es DEM o DIR le pasa las propieades a cuencas para el posterior trazado.\n'\
	'		de cuencas y tramos.\n' \
    #Abre el mapa
	direction=gdal.Open(ruta_map)
	#lee la informacion del mapa
	ncols=direction.RasterXSize
	nrows=direction.RasterYSize
	banda=direction.GetRasterBand(1)
	noData=banda.GetNoDataValue()
	geoT=direction.GetGeoTransform()
	dx=geoT[1]
	xll=geoT[0]; yll=geoT[3]-nrows*dx
	#lee el mapa
	Mapa=direction.ReadAsArray()
	direction.FlushCache()
	del direction
	if isDEMorDIR==True:
		cuencas.ncols=ncols
		cuencas.nrows=nrows
		cuencas.nodata=noData
		cuencas.dx=dx
		cuencas.xll=xll
		cuencas.yll=yll
		return Mapa.T
	else:
		return Mapa.T,[ncols,nrows,xll,yll,dx,noData]


#-----------------------------------------------------------------------
#Clase de cuencas
#-----------------------------------------------------------------------

class basin:
	#------------------------------------------------------
	# Subrutinas de lectura y prerpocesamiento inicial
	#------------------------------------------------------
	#Inicia la clase
	def __init__(self,lat,lon,DEM,DIR):
		'Descripcion: Inicia la variable de la cuenca, y la traza \n'\
		'	obtiene las propiedades basicas de la cuenca. \n'\
		'\n'\
		'Parametros\n'\
		'----------\n'\
		'self : Inicia las variables vacias.\n'\
		'lat : Coordenada en X de la salida de la cuenca.\n'\
		'lon : Coordenada en Y de la salida de la cuenca.\n'\
		'\n'\
		'Retornos\n'\
		'----------\n'\
		'self : Con las variables iniciadas.\n'\
		#Traza la cuenca 
		basin.nceldas = cuencas.basin_find(lat,lon,DIR,
			cuencas.ncols,cuencas.nrows)
		basin.cuenca = cuencas.basin_cut(basin.nceldas)
		#Obtiene parametros geomorfologicos
		basin.Param,basin.Tc=utils.basin_Tc(basin,
			DEM,DIR,cu.dxp,nceldas,cu.ncols,cu.nrows)
