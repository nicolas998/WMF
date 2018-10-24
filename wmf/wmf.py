#!Modelos: acoplado con cuencas presenta una serie de modelos hidrologicos distribuidos
#!Copyright (C) <2016>  <Nicolas Velasquez Giron>

#!This program is free software: you can redistribute it and/or modify
#!it under the terms of the GNU General Public License as published by
#!the Free Software Foundation, either version 3 of the License, or
#!(at your option) any later version.

#!This program is distributed in the hope that it will be useful,
#!but WITHOUT ANY WARRANTY; without even the implied warranty of
#!MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#!GNU General Public License for more details.

#!You should have received a copy of the GNU General Public License
#!along with this program.  If not, see <http://www.gnu.org/licenses/>.
#Algo
import matplotlib
matplotlib.use('Agg')
from cu import *
from models import *
import numpy as np
import pylab as pl
import osgeo.ogr, osgeo.osr
import gdal
from scipy.spatial import Delaunay
from scipy.stats import norm
import os
import pandas as pd
import datetime as datetime
from multiprocessing import Pool
import matplotlib.path as mplPath
try:
    import netcdf as netcdf
except:
    try:
        import netCDF4 as netcdf
    except:
        print('No netcdf en esta maquina, se desabilita la funcion SimuBasin.save_SimuBasin')
        pass
try:
    from mpl_toolkits.basemap import Basemap, addcyclic, shiftgrid, cm
    from matplotlib.patches import Polygon
    from matplotlib.collections import PatchCollection
except:
    print('No se logra importar basemap, por lo tanto no funciona Plot_basin')
    pass
try:
    from deap import base, creator
    from deap import tools
    FlagCalib_NSGAII = True
except:
    print('No se logra importar deap tools, por lo tanto se deshabilita SimuBasin.Calib_NSGAII')
    FlagCalib_NSGAII = False
try:
    from rasterio import features as __fea__
    FlagBasinPolygon = True
except:
    print('No se logra importar rasterio, se deshabilita obtencion de poligono de cuenca')
    FlagBasinPolygon = False

import random
#Variable codigo EPSG
Global_EPSG = -9999

#-----------------------------------------------------------------------
#Ploteo de variables
#-----------------------------------------------------------------------
def plot_sim_single(Qs,Qo=None,mrain=None,Dates=None,ruta=None,
    figsize=(8,4.5),ids=None,legend=True,ax1 = None,**kwargs):
    '''ENTRADAS:
    Qs = Caudal simulado, variable tipo lista
    Qo = Caudal observado, variable tipo lista
    mrain = Lluvia media observada,variable tipo lista
    Dates = Indice de tiempo (Datetime), variable tipo lista
    ids = Identificador de puntos de control
    ax1 = axis o entorno para graficar
    **kwargs = Argumentos por defecto, se pueden cambiar.
    Estos argumentos son:

    ARGUMENTOS POR DEFECTO
    rain_alpha: 0.4, Transparencia de la lluvia
    rain_color: blue, Color de la lluvia
    rain_lw : 0, Ancho de linea de lluvia
    rain_ylabel: Precipitation [$mm$], Etiqueta del eje y de la lluvia
    label_size: 14, Tamano de fuente
    rain_ylim : Se ajusta, Limite del eje y
    ColorSim : ['r','g','k','c','y'], Colores para simulacion de caudal
    Qs_lw : 1.5, Ancho de linea de caudal simulado
    Qo_lw : 2.0, Ancho de linea de caudal observado
    Qs_color : red, Color de linea de caudal simulado
    Qo_color : blue, Color de linea de caudal observado
    Qo_label : Observed, Leyenda de caudal observado
    Qs_label : Simulated, Leyenda de caudal simulado
    xlabel : Time [$min$], Etiqueta del eje x
    ylabel : Streamflow $[m^3/seg]$, Etiqueta del eje y
    legend_loc : upper center, Ubicacion de la leyenda
    bbox_to_anchor : (0.5,-0.12), Ajuste de la caja de la leyenda
    legend_ncol : 4, Numero de columnas para la leyenda
    '''
    #Argumentos kw
    show = kwargs.get('show',True)
    if ax1 == None:
        fig=pl.figure(facecolor='w',edgecolor='w',figsize=figsize)
        ax1=fig.add_subplot(111)
    else:
        show = False
    #Fechas
    if Dates==None:
        if len(Qs.shape)>1:
            ejeX=range(len(Qs[0]))
        else:
            ejeX = range(len(Qs))
    else:
        ejeX=Dates
    #Grafica la lluvia
    if mrain is not None:
        ax1AX=pl.gca()
        ax2=ax1.twinx()
        ax2AX=pl.gca()
        alpha = kwargs.get('rain_alpha',0.4)
        color = kwargs.get('rain_color','blue')
        lw = kwargs.get('rain_lw',0)
        ax2.fill_between(ejeX,0,mrain,alpha=alpha,color=color,lw=lw)
        ylabel = kwargs.get('rain_ylabel','Precipitation [$mm$]')
        label_size = kwargs.get('label_size',14)
        ax2.set_ylabel(ylabel,size=label_size)
        ylim = kwargs.get('rain_ylim',ax2AX.get_ylim() [::-1])
        ax2AX.set_ylim(ylim)
    else:
        ax2 = None
    #grafica las hidrografas
    ColorSim=kwargs.get('ColorSim',['r','g','k','c','y'])
    Qs_lw = kwargs.get('Qs_lw',1.5)
    Qo_lw = kwargs.get('Qo_lw',2.0)
    Qs_color = kwargs.get('Qs_color','r')
    Qo_color = kwargs.get('Qo_color','b')
    Qo_label = kwargs.get('Qo_label','Observed')
    Qs_label = kwargs.get('Qs_label','Simulated')
    if len(Qs.shape)>1:
        if ids is None:
            ids = np.arange(1,Qs.shape[0]+1)
        for i,c,d in zip(Qs,ColorSim,ids):
            ax1.plot(ejeX,i,c,lw=Qs_lw,label=str(d))
    else:
        ax1.plot(ejeX,Qs,Qs_color,lw=Qs_lw,label=Qs_label)
    if Qo is not None: ax1.plot(ejeX,Qo,Qo_color,lw=Qo_lw,label=Qo_label)
    #Pone elementos en la figura
    xlabel = kwargs.get('xlabel','Time [$min$]')
    ylabel = kwargs.get('ylabel','Streamflow $[m^3/seg]$')
    label_size = kwargs.get('label_size', 16)
    ax1.set_xlabel(xlabel,size=label_size)
    ax1.set_ylabel(ylabel,size=label_size)
    ax1.grid(True)
    loc = kwargs.get('legend_loc','upper center')
    bbox_to_anchor = kwargs.get('bbox_to_anchor',(0.5,-0.12))
    ncol = kwargs.get('legend_ncol',4)
    if legend == True:
        lgn1=ax1.legend(loc=loc, bbox_to_anchor=bbox_to_anchor,
            fancybox=True, shadow=True, ncol=ncol)
    if ruta is not None:
        pl.savefig(ruta, bbox_inches='tight')
    if show == True:
        pl.show()
    return ax1, ax2

def plot_mean_storage(Mean_Storage, Dates = None, mrain = None,
    rute = None, **kwargs):
    'Funcion: plot_mean_storage\n'\
    'Descripcion: Plotea como se encuentran los almacenamientos medios en la cuenca.\n'\
    'Parametros Obligatorios:.\n'\
    '   -Mean_Storage: almacenamiento medio en cada uno de los tanques (5,Npasos).\n'\
    '       Esta es la variables obtenida por wmf.SimuBasin.run_shia como Mean_Storage.\n'\
    'Parametros Opcionales:.\n'\
    '   -Dates: Fechas para el plot\n'\
    '   -marin: Lluvia promedio sobre la cuenca\n'\
    '   -rute: Ruta donde se va a guardar\n'\
    '   -kwargs: algunos argumentos de set del plot\n'\
    '       - figsize: tamano de la figura.\n'\
    '       - color: Color de los plot.\n'\
    '       - lw: ancho de la linea del plot.\n'\
    '       - labelsize: tamano de los label.\n'\
    '       - ysize: tamano de la fuente en el eje Y.\n'\
    'Retorno:.\n'\
    '   Plot de los almacenamientos medios.\n'\
    #Argumentos kw
    figsize = kwargs.get('figsize',(13,9))
    color = kwargs.get('color','k')
    lw = kwargs.get('lw',4)
    labelsize = kwargs.get('labelsize', 14)
    ysize = kwargs.get('ysize', 16)
    show = kwargs.get('show',True)
    alpha = kwargs.get('alpha',0.5)
    colorRain = kwargs.get('colorRain','b')
    lwRain = kwargs.get('lwRain',0.1)
    #Inicio de la funcion
    if Dates==None:
        ejeX=range(Mean_Storage.shape[1])
    else:
        ejeX=Dates
    #figura
    fig = pl.figure(figsize = figsize)
    nombres = ['Hu','Runoff','Hg','Sub','Stream']
    for c,i in enumerate(Mean_Storage):
        ax = fig.add_subplot(5,1,c+1)
        ax.plot(ejeX, i, color, lw = lw)
        ax.grid()
        # Deja el eje X solo en el ultimo plot
        if c<4:
            ax.set_xticklabels([])
        ax.tick_params(labelsize = labelsize)
        #Pinta la lluvia en el primer plot
        if c == 0 and mrain  is not  None:
            ax2=ax.twinx()
            ax2AX=pl.gca()
            ax2.fill_between(ejeX,0,mrain,alpha=alpha,color=colorRain,lw=lwRain)
            ylim = ax2AX.get_ylim()[::-1]; ylim = list(ylim); ylim[1] = 0
            ax2AX.set_ylim(ylim)
        #Nombre de cada tanque
        ax.set_ylabel(nombres[c], size = ysize)
    if rute is not None:
        pl.savefig(rute, bbox_inches='tight')
    if show == True:
        pl.show()
    if c == 0 and mrain  is not  None:
        return ax, ax2
    else:
        return ax

#-----------------------------------------------------------------------
#Lectura de informacion y mapas
#-----------------------------------------------------------------------
def read_map_raster(ruta_map,isDEMorDIR=False,dxp=None, noDataP = None,isDIR = False,DIRformat = 'r.watershed'):
    'Funcion: read_map\n'\
    'Descripcion: Lee un mapa raster soportado por GDAL.\n'\
    'Parametros Obligatorios:.\n'\
    '   -ruta_map: Ruta donde se encuentra el mapa.\n'\
    'Parametros Opcionales:.\n'\
    '   -isDEMorDIR: Pasa las propiedades de los mapas al modulo cuencas \n'\
    '       escrito en fortran \n'\
    '   -dxp: tamano plano del mapa\n'\
    '   -noDataP: Valor para datos nulos en el mapa (-9999)\n'\
    '   -DIRformat: donde se ha conseguido el mapa dir (r.watershed) \n'\
    '       - r.watershed: mapa de direcciones obtenido por la funcion de GRASS\n'\
    '       - opentopo: mapa de direcciones de http://www.opentopography.org/\n'\
    '   -isDIR: (FALSE) es este un mapa de direcciones\n'\
    'Retorno:.\n'\
    '   Si no es DEM o DIR retorna todas las propieades del elemento en un vector.\n'\
    '       En el siguiente orden: ncols,nrows,xll,yll,dx,nodata.\n'\
    '   Si es DEM o DIR le pasa las propieades a cuencas para el posterior trazado.\n'\
    '       de cuencas y tramos.\n' \
    #Abre el mapa
    direction=gdal.Open(ruta_map)
    #Projection
    proj = osgeo.osr.SpatialReference(wkt=direction.GetProjection())
    EPSG_code = proj.GetAttrValue('AUTHORITY',1)
    #lee la informacion del mapa
    ncols=direction.RasterXSize
    nrows=direction.RasterYSize
    banda=direction.GetRasterBand(1)
    noData=banda.GetNoDataValue()
    geoT=direction.GetGeoTransform()
    dx=geoT[1]
    dy = np.abs(geoT[-1])
    xll=geoT[0]; yll=geoT[3]-nrows*dy
    #lee el mapa
    Mapa=direction.ReadAsArray()
    direction.FlushCache()
    del direction
    if isDEMorDIR==True:
        cu.ncols=ncols
        cu.nrows=nrows
        if noDataP is not None:
            cu.nodata = noDataP
            Mapa[Mapa<0] = cu.nodata
        else:
            cu.nodata=noData
        cu.dx=dx
        cu.dy = dy
        cu.xll=xll
        cu.yll=yll
        if dxp==None:
            cu.dxp=30.0
        else:
            cu.dxp=dxp
        #Guarda la variable para el proyecto
        global Global_EPSG
        Global_EPSG = EPSG_code
        # si es un dir se fija si es de r.watershed
        if isDIR:
            if DIRformat == 'r.watershed':
                Mapa[Mapa<=0] = cu.nodata.astype(int)
                Mapa = cu.dir_reclass_rwatershed(Mapa.T,cu.ncols,cu.nrows)
                return Mapa, EPSG_code
            if DIRformat == 'opentopo':
                Mapa[Mapa<=0] = cu.nodata.astype(int)
                Mapa = cu.dir_reclass_opentopo(Mapa.T,cu.ncols,cu.nrows)
                return Mapa, EPSG_code
        #retorna el mapa
        return Mapa.T.astype(float),EPSG_code
    else:
        return Mapa.T.astype(float),[ncols,nrows,xll,yll,dx,dy,noData],EPSG_code

def read_map_points(ruta_map, ListAtr = None):
    'Funcion: read_map_points\n'\
    'Descripcion: Lee un mapa vectorial de puntos soportado por GDAL.\n'\
    'Parametros Obligatorios:.\n'\
    '   -ruta_map: Ruta donde se encuentra el mapa.\n'\
    'Parametros Opcionales:.\n'\
    '   -ListAtr: Lista con los nombres de los atributos de las columnas\n'\
    '       que se quieren leer dentro de la variable Dict\n'\
    'Retorno:.\n'\
    '   Si ListAtr == None: Retorna unicamente las coordenadas.\n'\
    '   Si ListAtr == [Nombre1, Nombre2, ...]: Retorna: Coord y diccionario con variables.\n'\
    #Obtiene el acceso
    dr = osgeo.ogr.Open(ruta_map)
    l = dr.GetLayer()
    #Lee las coordenadas
    Cord = []
    for i in range(l.GetFeatureCount()):
        f = l.GetFeature(i)
        g = f.GetGeometryRef()
        pt = [g.GetX(),g.GetY()]
        Cord.append(pt)
    Cord = np.array(Cord).T
    #Si hay atributos para buscar los lee y los m
    if ListAtr is not None:
        Dict = {}
        for j in ListAtr:
            #Busca si el atributo esta
            pos = f.GetFieldIndex(j)
            #Si esta lee la info del atributo
            if pos is not -1:
                vals = []
                for i in range(l.GetFeatureCount()):
                    f = l.GetFeature(i)
                    vals.append(f.GetField(pos))
                Dict.update({j:np.array(vals)})
        #Cierra el mapa
        dr.Destroy()
        return Cord, Dict
    else:
        #Cierra el mapa
        dr.Destroy()
        return Cord

def Save_Array2Raster(Array, ArrayProp, ruta, EPSG = 4326, Format = 'GTiff'):
    dst_filename = ruta
    #Formato de condiciones del mapa
    x_pixels = Array.shape[0]  # number of pixels in x
    y_pixels = Array.shape[1]  # number of pixels in y
    PIXEL_SIZE = ArrayProp[4]  # size of the pixel...
    x_min = ArrayProp[2]
    y_max = ArrayProp[3] + ArrayProp[4] * ArrayProp[1] # x_min & y_max are like the "top left" corner.
    driver = gdal.GetDriverByName(Format)
    #Para encontrar el formato de GDAL
    NP2GDAL_CONVERSION = {
      "uint8": 1,
      "int8": 1,
      "uint16": 2,
      "int16": 3,
      "uint32": 4,
      "int32": 5,
      "float32": 6,
      "float64": 7,
      "complex64": 10,
      "complex128": 11,
    }
    gdaltype = NP2GDAL_CONVERSION[Array.dtype.name]
    # Crea el driver
    dataset = driver.Create(
        dst_filename,
        x_pixels,
        y_pixels,
        1,
        gdaltype,)
    #coloca la referencia espacial
    dataset.SetGeoTransform((
        x_min,    # 0
        PIXEL_SIZE,  # 1
        0,                      # 2
        y_max,    # 3
        0,                      # 4
        -PIXEL_SIZE))
    #coloca la proyeccion a partir de un EPSG
    proj = osgeo.osr.SpatialReference()
    texto = 'EPSG:' + str(EPSG)
    proj.SetWellKnownGeogCS( texto )
    dataset.SetProjection(proj.ExportToWkt())
    #Coloca el nodata
    band = dataset.GetRasterBand(1)
    if ArrayProp[-1] == None:
        band.SetNoDataValue(wmf.cu.nodata.astype(int).max())
    else:
        band.SetNoDataValue(-9999)
    #Guarda el mapa
    dataset.GetRasterBand(1).WriteArray(Array.T)
    dataset.FlushCache()

def Save_Points2Map(XY,ids,ruta,EPSG = 4326, Dict = None,
    DriverFormat='ESRI Shapefile'):
    'Funcion: Save_Points2Map\n'\
    'Descripcion: Guarda coordenadas X,Y en un mapa tipo puntos.\n'\
    'Parametros Opcionales:.\n'\
    'XY: Coordenadas de los puntos X,Y.\n'\
    'ids: Numero que representa a cada punto.\n'\
    'ruta: Ruta de escritura de los puntos.\n'\
    'Opcionales:.\n'\
    'EPSG : Codigo del tipo de proyeccion usada, defecto 4326 de WGS84.\n'\
    'Dict : Diccionario con prop de los puntos, Defecto: None.\n'\
    'DriverFormat : Formato del archivo vectorial, Defecto: ESRI Shapefile.\n'\
    'Retorno:.\n'\
    '   Escribe el mapa en la ruta especificada.\n'\
    #Genera el shapefile
    spatialReference = osgeo.osr.SpatialReference()
    spatialReference.ImportFromEPSG(EPSG)
    driver = osgeo.ogr.GetDriverByName(DriverFormat)
    if os.path.exists(ruta):
         driver.DeleteDataSource(ruta)
    shapeData = driver.CreateDataSource(ruta)
    layer = shapeData.CreateLayer('layer1', spatialReference, osgeo.ogr.wkbPoint)
    layerDefinition = layer.GetLayerDefn()
    new_field=osgeo.ogr.FieldDefn('Estacion',osgeo.ogr.OFTReal)
    layer.CreateField(new_field)
    if Dict is not None:
        for p in Dict.keys():
            tipo = type(Dict[p][0])
            if tipo is np.float64:
                new_field=osgeo.ogr.FieldDefn(p[:10],osgeo.ogr.OFTReal)
            elif tipo is np.int64:
                new_field=osgeo.ogr.FieldDefn(p[:10],osgeo.ogr.OFTInteger)
            elif tipo is str:
                new_field=osgeo.ogr.FieldDefn(p[:10],osgeo.ogr.OFTString)
            layer.CreateField(new_field)
    #Mete todos los puntos
    featureIndex=0
    contador=0
    #Calcula el tamano de la muestra
    N=np.size(XY,axis=1)
    if N>1:
        for i in XY.T:
            if i[0]!=-9999.0 and i[1]!=-9999.0:
                #inserta el punto
                point = osgeo.ogr.Geometry(osgeo.ogr.wkbPoint)
                point.SetPoint(0, float(i[0]), float(i[1]))
                feature = osgeo.ogr.Feature(layerDefinition)
                feature.SetGeometry(point)
                feature.SetFID(featureIndex)
                feature.SetField('Estacion',ids[contador])
                #Le coloca lo del diccionario
                if Dict is not None:
                    for p in Dict.keys():
                        tipo = type(Dict[p][0])
                        if tipo is np.float64:
                            feature.SetField(p[:10],float("%.2f" % Dict[p][contador]))
                        elif tipo is np.int64:
                            feature.SetField(p[:10],int("%d" % Dict[p][contador]))
                        elif tipo is str:
                            feature.SetField(p[:10],Dict[p][contador])
                #Actualiza contadores
                contador+=1
                layer.CreateFeature(feature)
                point.Destroy()
                feature.Destroy()
    else:
        point = osgeo.ogr.Geometry(osgeo.ogr.wkbPoint)
        point.SetPoint(0, float(XY[0,0]), float(XY[1,0]))
        feature = osgeo.ogr.Feature(layerDefinition)
        feature.SetGeometry(point)
        feature.SetFID(featureIndex)
        feature.SetField('Estacion',ids)
        for p in Dict.keys():
            feature.SetField(p[:p.index('[')].strip()[:10],float("%.2f" % Param[p]))
        contador+=1
        layer.CreateFeature(feature)
        point.Destroy()
        feature.Destroy()
    shapeData.Destroy()

def __ListaRadarNames__(ruta,FechaI,FechaF,fmt,exten,string,dt):
    'Funcion: OCG_param\n'\
    'Descripcion: Obtiene una lista con los nombres para leer datos de radar.\n'\
    'Parametros:.\n'\
    '   -FechaI, FechaF: Fechas inicial y final en formato datetime.\n'\
    '   -fmt : Formato en el que se pasa FechaI a texto, debe couincidir con.\n'\
    '       el del texto que se tenga en los archivos.\n'\
    '   -exten : Extension de los archivos .asc, .nc, .bin ...\n'\
    '   -string : texto antes de la fecha.\n'\
    '   -dt : Intervalos de tiempo entre los eventos.\n'\
    'Retorno:.\n'\
    '   Lista : la lista de python con los nombres de los binarios.\n'\
    #Crea lista de fechas y mira que archivos estan en esas fechas
    Dates=[FechaI]; date=FechaI
    while date<FechaF:
        date+=datetime.timedelta(minutes=dt)
        Dates.append(date)
    #Mira que archivos estan en esas fechas
    Lista=[]; L=os.listdir(ruta)
    DatesFin = []
    for i in Dates:
        try:
            stringB=string+i.strftime(fmt)+exten
            Pos=L.index(stringB)
            Lista.append(L[Pos])
            DatesFin.append(i)
        except:
            pass
    return Lista,DatesFin

def __Add_hdr_bin_2route__(rute,storage=False):
    if storage == False:
        if rute.endswith('.bin') == False and rute.endswith('.hdr') == True:
            ruteBin = rute[:-3] + 'bin'
            ruteHdr = rute
        elif rute.endswith('.bin') == True and rute.endswith('.hdr') == False:
            ruteBin = rute
            ruteHdr = rute[:-3] + 'hdr'
        elif  rute.endswith('.bin') == False and rute.endswith('.hdr') == False:
            ruteBin = rute + '.bin'
            ruteHdr = rute + '.hdr'
    else:
        if rute.endswith('.StObin') == False and rute.endswith('.StOhdr') == True:
            ruteBin = rute[:-6] + 'StObin'
            ruteHdr = rute[:-6] + 'StOhdr'
        elif rute.endswith('.StObin') == True and rute.endswith('.StOhdr') == False:
            ruteBin = rute[:-6] + 'StObin'
            ruteHdr = rute[:-6] + 'StOhdr'
        elif  rute.endswith('.StObin') == False and rute.endswith('.StOhdr') == False:
            ruteBin = rute + '.StObin'
            ruteHdr = rute + '.StOhdr'
        elif rute.endswith('.bin') == False and rute.endswith('.hdr') == True:
            ruteBin = rute[:-3] + 'StObin'
            ruteHdr = rute[:-3] + 'StOhdr'
        elif rute.endswith('.bin') == True and rute.endswith('.hdr') == False:
            ruteBin = rute[:-3] + 'StObin'
            ruteHdr = rute[:-3] + 'StOhdr'
        elif  rute.endswith('.bin') == False and rute.endswith('.hdr') == False:
            ruteBin = rute + '.StObin'
            ruteHdr = rute + '.StOhdr'
    return ruteBin,ruteHdr

def read_mean_rain(ruta,Nintervals=None,FirstInt=0):
    #Abrey cierra el archivo plano
    Data = np.loadtxt(ruta,skiprows=6,usecols=(2,3),delimiter=',',dtype='str')
    Rain = np.array([float(i[0]) for i in Data])
    Dates = [datetime.datetime.strptime(i[1],' %Y-%m-%d-%H:%M  ') for i in Data]
    #Corrige pedazo para capturar
    if Nintervals == None: Nintervals = Rain.size
    #Obtiene el pedazo
    Dates = Dates[FirstInt:FirstInt+Nintervals]
    Rain = Rain[FirstInt:FirstInt+Nintervals]
    #Regresa el resultado de la funcion
    return pd.Series(Rain,index = pd.to_datetime(Dates))

def read_rain_struct(ruta):
    D = pd.read_csv(ruta,skiprows=5,
    index_col=2, parse_dates=True,
    infer_datetime_format=True,
    usecols = (1,2,3))
    return D

def read_storage_struct(ruta):
    '''Lee la estructura del archivo encabezado de almacenamiento'''
    #Obtiene rutaHdr
    PathBin, PathHdr = __Add_hdr_bin_2route__(ruta,storage=True)
    #Lee el archivo
    Data = pd.read_csv(PathHdr,
        skiprows=4,
        index_col=6,
        parse_dates=True)
    return Data

def __Save_storage_hdr__(rute,rute_rain,Nintervals,FirstInt,cuenca,
    Mean_Storage):
    #Lee fechas para el intervalo de tiempo
    S = read_mean_rain(rute_rain,Nintervals,FirstInt)
    #Escribe el encabezado del archivo
    f=open(rute,'w')
    f.write('Numero de celdas: %d \n' % cuenca.ncells)
    f.write('Numero de laderas: %d \n' % cuenca.nhills)
    f.write('Numero de registros: %d \n' % Nintervals)
    f.write('Tipo Modelo: %s \n' % cuenca.modelType)
    f.write('IDfecha, Tanque 1, Tanque 2, Tanque 3, Tanque 4, Tanque 5, Fecha \n')
    c = 1
    #Si no hay almacenamiento medio lo coloca en -9999
    #Escribe registros medios y fechas de los almacenamientos
    for d,sto in zip(S.index.to_pydatetime(),Mean_Storage.T):
        f.write('%d, \t %.2f, \t %.4f, \t %.4f, \t %.2f, \t %.2f, %s \n' %
            (c,sto[0],sto[1],sto[2],sto[3],sto[4],d.strftime('%Y-%m-%d-%H:%M')))
        c+=1
    f.close()

def __Save_speed_hdr__(rute,rute_rain,Nintervals,FirstInt,cuenca,
    Mean_Speed = None):
    #Lee fechas para el intervalo de tiempo
    S = read_mean_rain(rute_rain,Nintervals,FirstInt)
    #Escribe el encabezado del archivo
    f=open(rute,'w')
    f.write('Numero de celdas: %d \n' % cuenca.ncells)
    f.write('Numero de laderas: %d \n' % cuenca.nhills)
    f.write('Numero de registros: %d \n' % Nintervals)
    f.write('Tipo Modelo: %s \n' % cuenca.modelType)
    f.write('IDfecha, Tanque 2, Tanque 3, Tanque 4, Tanque 5, Fecha \n')
    c = 1
    #Si no hay almacenamiento medio lo coloca en -9999
    if Mean_Speed == None:
        Mean_Speed = np.ones((5,Nintervals))*-9999
    #Escribe registros medios y fechas de los almacenamientos
    for d,sto in zip(S.index.to_pydatetime(),Mean_Speed.T):
        f.write('%d, \t %.2f, \t %.4f, \t %.2f, \t %.2f, %s \n' %
            (c,sto[0],sto[1],sto[2],sto[3],d.strftime('%Y-%m-%d-%H:%M')))
        c+=1
    f.close()

def __Save_retorno_hdr__(rute,rute_rain,Nintervals,FirstInt,cuenca,
    Mean_retorno = None):
    #Lee fechas para el intervalo de tiempo
    S = read_mean_rain(rute_rain,Nintervals,FirstInt)
    #Escribe el encabezado del archivo
    f=open(rute,'w')
    f.write('Numero de celdas: %d \n' % cuenca.ncells)
    f.write('Numero de laderas: %d \n' % cuenca.nhills)
    f.write('Numero de registros: %d \n' % Nintervals)
    f.write('Tipo Modelo: %s \n' % cuenca.modelType)
    f.write('IDfecha, Retorno[mm], Fecha \n')
    c = 1
    #Si no hay almacenamiento medio lo coloca en -9999
    if Mean_retorno == None:
        Mean_retorno = np.ones(Nintervals)*-9999
    #Escribe registros medios y fechas de los almacenamientos
    for d,sto in zip(S.index.to_pydatetime(),Mean_retorno):
        f.write('%d, \t %.2f, %s \n' % (c, sto,d.strftime('%Y-%m-%d-%H:%M')))
        c+=1
    f.close()

def map_acum_to_stream(ACUM,umbral):
    'Funcion: map_acum_to_stream\n'\
    'Descripcion: Calcula red hidrica a partir del area acumulada y un.\n'\
    '   umbral.\n'\
    'Parametros :.\n'\
    '   -ACUM: Mapa de celdas acumuladas.\n'\
    '   -umbral: Umbral para la generacion de cauce.\n'\
    'Retorno:.\n'\
    '   CAUCE: Mapa binario con los cauces: 1 cauce, 0 ladera.\n'\
    #Invoca funcion de fortran
    CAUCE = cu.geo_acum_to_cauce(ACUM,umbral,cu.ncols,cu.nrows)
    return CAUCE

def SimuBains_Update_DEM_DIR(ruta_basin, rute_dem, rute_dir):
    'Funcion: map_acum_to_stream\n'\
    'Descripcion: Actualiza la ruta al DEM y al DIR de un proyecto de simulacion.\n'\
    'Parametros :.\n'\
    '   -ruta_basin: ruta del proyecto de la cuenca.\n'\
    '   -ruta_dem: Ruta al mapa dem.\n'\
    '   -ruta_dir: Ruta al mapa dir .\n'\
    'Retorno:.\n'\
    '   actualiza las rutas en el proyecto.\n'\
    #Lee el nc y le actualiza ambas rutas
    g = netcdf.Dataset(ruta_basin,'a')
    g.DEM = rute_dem
    g.DIR = rute_dir
    g.close()

#-----------------------------------------------------------------------
#Ecuaciones Que son de utilidad
#-----------------------------------------------------------------------
def OCG_param(alfa=[0.75,0.2],sigma=[0.0,0.225,0.225],    c1=5.54,k=0.5,fhi=0.95,Omega=0.13,pend=None,area=None,):
    'Funcion: OCG_param\n'\
    'Descripcion: Calcula los parametros de la onda cinematica.\n'\
    'geomorfologica (Velez, 2001).\n'\
    'Parametros Opcionales:.\n'\
    '   -isDEMorDIR: Pasa las propiedades de los mapas al modulo cuencas \n'\
    '       escrito en fortran \n'\
    'Retorno:.\n'\
    '   Parametros: B, w1, w2 y w3.\n'\
    '       si se entregan los mapas de pend, aacum entrega h_coef(4,:) .\n' \
    '       se asume que w1 corresponde a h_exp(4,:) .\n' \
    #Calcula los parametros de la ecuacion de onda cinematica
    B = Omega*(c1*k**(alfa[0]-alfa[1]))**((2.0/3.0)-alfa[1])
    eB=1.0/(1+alfa[1]*((2/3.0)-sigma[1]))
    w1=((2/3.0)-sigma[1])*(1.0-alfa[1])*eB
    w2=(1+alfa[1]*((2/3.0)-sigma[1]))/(fhi*((2/3.0)-sigma[1])*(alfa[0]-alfa[1])+sigma[0])
    w2=(fhi*(0.667-sigma[1])*(alfa[1]-alfa[0])+sigma[0])*eB
    w3=(0.5-sigma[2])*eB
    B=B**(-eB)
    if pend is not None and area is not None:
        var = B*(pend**w2)*(area**w3)
        return var,w1
    else:
        return B,w1,w2,w3

def PotCritica(S,D,te = 0.056):
    ti = te * (D*1600*9.8)
    return ti *(8.2* (((ti/(1000*9.8*S))/D)**(1.0/6.0)) * np.sqrt(ti/1000.0))/9800.0

def __eval_nash__(s_o,s_s):
    med_s_o=np.nansum(s_o)/np.sum(np.isfinite(s_o))
    Dif_sim=(s_o-s_s[:len(s_o)])**2
    Dif_med=(s_o-med_s_o)**2
    E=1-(np.nansum(Dif_sim)/np.nansum(Dif_med))
    return E

def __eval_rmse__(So,Ss):
    Dif=(So-Ss)**2
    Dif=np.ma.array(Dif,mask=np.isnan(Dif))
    return np.sqrt(Dif.sum()/Dif.shape[0])

def __eval_rmse_log__(So,Ss):
    So=np.ma.array(So,mask=np.isnan(So))
    Ss=np.ma.array(Ss,mask=np.isnan(Ss))
    So=np.log(So); Ss=np.log(Ss)
    Dif=(So-Ss)**2
    Dif=np.ma.array(Dif,mask=np.isnan(Dif))
    return np.sqrt(Dif.sum()/Dif.shape[0])

def __eval_t_pico__(s_o,s_s,dt):
    max_o=np.argmax(np.ma.array(s_o,mask=np.isnan(s_o)))
    max_s=np.argmax(np.ma.array(s_s,mask=np.isnan(s_s)))
    dif_tpico=(max_o-max_s)*dt
    return dif_tpico

def __eval_q_pico__(s_o,s_s):
    max_o=np.argmax(np.ma.array(s_o,mask=np.isnan(s_o)))
    max_s=np.argmax(np.ma.array(s_s,mask=np.isnan(s_s)))
    Qo_max=s_o[max_o]
    Qs_max=s_s[max_s]
    dif_qpico=((Qo_max-Qs_max)/Qo_max)*100
    return dif_qpico


#Funciones para mirar como es un netCDf por dentro
def netCDf_varSumary2DataFrame(ruta, print_netCDF = False):
    'Funcion: netCDf_var_view\n'\
    'Descripcion: Muestra las variables que estan cargadas en el netCDF.\n'\
    'Parametros Opcionales:.\n'\
    '   -ruta: Ruta donde se encuentra el netCDF\n'\
    '   -print_netCDF: Imprime la info generica del netCDF.\n'\
    'Retorno:.\n'\
    '   DataFram de pandas con las variables del netCDF.\n'\
    # lectura
    g = netcdf.Dataset(ruta)
    Dict = {}
    for k in g.variables.keys():
        D = {'type': g.variables[k].datatype, 'dimensions': g.variables[k].dimensions}
        Dict.update({k:D})
    # if print
    if print_netCDF:
        print(g)
    g.close()
    return pd.DataFrame.from_dict(Dict, orient='index')

#Funciones de ejecucion en paralelo del modelo
def __multiprocess_Warper__(Lista):
    Res = Lista[4].run_shia(Lista[0],Lista[1],Lista[2],Lista[3])
    return Res

def __ejec_parallel__(ListEjecs, nproc, nodo):
    P = Pool(processes=nproc)
    Res = P.map(__multiprocess_Warper__, ListEjecs)
    Lista = [i['Qsim'][nodo][0] for i in Res]
    P.close()
    P.join()
    return Lista

#-----------------------------------------------------------------------
#Transformacion de datos
#-----------------------------------------------------------------------
def __ModifyElevErode__(X,slope=0.01,d2 = 0.03, window = 25):
    # Obtiene la corriente quitando los puntos en donde se eleva
    Pos = []
    Y = np.copy(X)
    c = 1
    for i,j in zip(X[:-1],X[1:]):
        if j>i:
            Flag = True
            c2 = c
            c3 = 0
            while Flag:
                if i<Y[c2]:
                    Y[c2] = i-slope*c3
                    c2+=1
                    c3+=1
                    if c2 >= Y.size: Flag = False
                else:
                    Flag=False
                    Pos.append(c2)
        c+=1
    # Obtiene la segunda derivada de la corriente corregida

    Y = pd.Series(Y)
    Y = Y.rolling(window).mean()
    #Y = pd.rolling_mean(Y,window)
    l = Y[np.isnan(Y)].shape[0]
    Y[:l]=Y[l+1]
    slope = np.diff(Y,1)
    slope = np.hstack([slope,slope[-1]])
    return Y,np.abs(slope)

#-----------------------------------------------------------------------
#Clase de cuencas
#-----------------------------------------------------------------------

class Basin:
    #------------------------------------------------------
    # Subrutinas de trazado de cuenca y obtencion de parametros
    #------------------------------------------------------
    #Inicia la cuenca
    def __init__(self,lat=0,lon=0,DEM=None,DIR=None,name='NaN',stream=None,
        umbral=1000, ruta = None, useCauceMap = None):
        'Descripcion: Inicia la variable de la cuenca, y la traza \n'\
        '   obtiene las propiedades basicas de la cuenca. \n'\
        '\n'\
        'Parametros\n'\
        '----------\n'\
        'self : Inicia las variables vacias.\n'\
        'lat : Coordenada en X de la salida de la cuenca.\n'\
        'lon : Coordenada en Y de la salida de la cuenca.\n'\
        'name : Nombre con el que se va a conocer la cuenca.\n'\
        'stream : Opcional, si se coloca, las coordenadas no tienen.\n'\
        '   que ser exactas, estas se van a corregir para ubicarse.\n'\
        '   en el punto mas cercano dentro de la corriente, este.\n'\
        '   debe ser un objeto del tipo stream.\n'\
        'umbral : umbral minimo para la creacion de cauce (defecto =1000).\n'\
        'ruta : Ruta donde se encuentra un archivo binario de la cuenca.\n'\
        'useCauceMap : Si se asigna un mapa binario donde 1 es cauce y 0 es ladera.\n'\
        '   el trazador corregira las coordenadas para que estas lleguen a una celda.\n'\
        '   tipo cauce y se trace la cuenca de forma correcta (esta opcion deshabilita \n'\
        '   a stream)\n'\
        '\n'\
        'Retornos\n'\
        '----------\n'\
        'self : Con las variables iniciadas.\n'\
        #Relaciona con el DEM y el DIR
        self.DEM=DEM
        self.DIR=DIR
        #Si se da la opcion de que use el useCauceMap deshabilita stream
        if useCauceMap is not None:
            stream = None
        #Si se entrega cauce corrige coordenadas
        if ruta == None:
            # Si se da el stream corrige por corriente
            if stream is not None:
                error=[]
                for i in stream.structure.T:
                    error.append( np.sqrt((lat-i[0])**2+(lon-i[1])**2) )
                loc=np.argmin(error)
                lat=stream.structure[0,loc]
                lon=stream.structure[1,loc]
            if useCauceMap is not None and useCauceMap.shape == DEM.shape:
                lat,lon = cu.stream_find_to_corr(lat,lon,DEM,DIR,useCauceMap,
                    cu.ncols,cu.nrows)
            #copia la direccion de los mapas de DEM y DIR, para no llamarlos mas
            self.name=name
            #Traza la cuenca
            self.ncells = cu.basin_find(lat,lon,DIR,
                cu.ncols,cu.nrows)
            self.structure = cu.basin_cut(self.ncells)
            self.umbral = umbral
            self.DEMvec = self.Transform_Map2Basin(DEM,[cu.ncols, cu.nrows, cu.xll, cu.yll, cu.dx, cu.dy])
            self.DIRvec = self.Transform_Map2Basin(DIR,[cu.ncols, cu.nrows, cu.xll, cu.yll, cu.dx, cu.dy])
        else:
            self.__Load_BasinNc(ruta)
            #Genera el poligono de la cuenca
            
        self.__GetBasinPolygon__()
    #Cargador de cuenca
    def __Load_BasinNc(self,ruta,Var2Search=None):
        'Descripcion: Lee una cuenca posteriormente guardada\n'\
        '   La cuenca debio ser guardada con Basin.Save_Basin2nc\n'\
        '\n'\
        'Parametros\n'\
        '----------\n'\
        'self : Inicia las variables vacias.\n'\
        'ruta : ruta donde se encuentra ubicada la cuenca guardada\n'\
        'Retornos\n'\
        '----------\n'\
        'self : La cuenca con sus parametros ya cargada.\n'\
        #Abre el archivo binario de la cuenca
        self.rutaNC = ruta
        gr = netcdf.Dataset(ruta,'a')
        #obtiene las prop de la cuenca
        self.name = gr.nombre
        self.ncells = gr.ncells
        self.umbral = gr.umbral
        #Obtiene las prop de los mapas
        cu.ncols=gr.ncols
        cu.nrows=gr.nrows
        cu.nodata=gr.noData
        cu.dx=gr.dx
        cu.xll=gr.xll
        cu.yll=gr.yll
        cu.dxp=gr.dxp
        #Obtiene las variables vectoriales
        self.structure = gr.variables['structure'][:]
        #Cierra el archivo
        gr.close()

    def Load_BasinVar(self, varName):
        'Descripcion: Lee una variable especifica del netCDf de la cuenca\n'\
        '\n'\
        'Parametros\n'\
        '----------\n'\
        'varName: nombre de la variable guardada\n'\
        'Retornos\n'\
        '----------\n'\
        'var : Retorna la variable como un numpy array.\n'\
        #Lectura del netCDf y lectura de variable
        gr = netcdf.Dataset(self.rutaNC,'a')
        Var = gr.variables[varName][:]
        gr.close()
        return Var

    #Guardado de de la cuenca en nc
    def Save_Basin2nc(self,ruta,qmed=None,q233=None,q5=None,
        ExtraVar=None):
        'Descripcion: guarda una cuenca previamente ejecutada\n'\
        '\n'\
        'Parametros\n'\
        '----------\n'\
        'ruta : Ruta donde la cuenca sera guardada.\n'\
        'qmed : Matriz con caudal medio estimado.\n'\
        'q233 : Matriz con caudal minimo de 2.33.\n'\
        'q5 : Matriz con caudal minimo de 5.\n'\
        'ExtraVar: Diccionario con variables extras que se quieran guardar,.\n'\
        '   se guardan como flotantes.\n'\
        '\n'\
        'Retornos\n'\
        '----------\n'\
        'self : Con las variables iniciadas.\n'\
        #Diccionario con las propiedades
        Dict = {'nombre':self.name,
            'noData':cu.nodata,
            'ncells':self.ncells,
            'umbral':self.umbral,
            'dxp':cu.dxp,
            'dx':cu.dx,
            'xll':cu.xll,
            'yll':cu.yll,
            'ncols':cu.ncols,
            'nrows':cu.nrows}
        #abre el archivo
        gr = netcdf.Dataset(ruta,'w',format='NETCDF4')
        #Establece tamano de las variables
        DimNcell = gr.createDimension('ncell',self.ncells)
        DimCol3 = gr.createDimension('col3',3)
        #Crea variables
        VarStruc = gr.createVariable('structure','i4',('col3','ncell'),zlib=True)
        VarQmed = gr.createVariable('q_med','f4',('ncell',),zlib=True)
        VarQ233 = gr.createVariable('q_233','f4',('ncell',),zlib=True)
        VarQ5 = gr.createVariable('q_5','f4',('ncell',),zlib=True)
        #Variables opcionales
        if type(ExtraVar) is dict:
            for k in ExtraVar.keys():
                Var = gr.createVariable(k,'f4',('ncell',),zlib=True)
                Var[:] = ExtraVar[k]
        #Asigna valores a las variables
        VarStruc[:] = self.structure
        if qmed is not None:
            VarQmed[:] = qmed
        if q233 is not None:
            VarQ233[:] = q233
        if q5 is not None:
            VarQ5[:] = q5
        #asignlas prop a la cuenca
        gr.setncatts(Dict)
        #Cierra el archivo
        gr.close()
        #Sale del programa
        return
    #Parametros Geomorfologicos
    def GetGeo_Parameters(self,rutaParamASC=None,plotTc=False,
        rutaTcPlot = None, figsize=(8,5), GetPerim=True):
        'Descripcion: Obtiene los parametros geomorfologicos de la cuenca \n'\
        '   y los tiempos de concentracion calculados por diferentes metodologias. \n'\
        '\n'\
        'Parametros\n'\
        '   rutaParamASC: ruta del ascii donde se escriben los param.\n'\
        '   plotTc: Plotea o no los tiempos de concentracion.\n'\
        '   rutaTcPlot: Si se da se guarda la figura de tiempos de concentracion.\n'\
        '----------\n'\
        '\n'\
        'Retornos\n'\
        '----------\n'\
        'GeoParameters : Parametros de la cuenca calculados.\n'\
        'Tc :  Tiempo de concentracion calculado para la cuenca.\n'\
        #Calcula lo que se necesita para sacar los parametros
        acum,longCeld,slope,Elev=cu.basin_basics(self.structure,
            self.DEMvec,self.DIRvec,self.ncells)
        Lpma,puntto=cu.basin_findlong(self.structure,self.ncells)
        cauce,nodos,trazado,n_nodos,n_cauce = cu.basin_stream_nod(self.structure,
            acum,self.umbral,self.ncells)
        ppal_nceldas,punto = cu.basin_ppalstream_find(self.structure,
            nodos,longCeld,Elev,self.ncells)
        ppal = cu.basin_ppalstream_cut(ppal_nceldas,self.ncells)
        self.hipso_main,self.hipso_basin=cu.basin_ppal_hipsometric(
            self.structure,Elev,punto,30,ppal_nceldas,self.ncells)
        self.main_stream=ppal
        #Obtiene los parametros
        Perim = self.Polygon.shape[1]*cu.dxp/1000.
        Area=(self.ncells*cu.dxp**2)/1e6
        TotalCauces = self.CellCauce*self.CellLong
        TotalCauces = TotalCauces.sum() / 1000. #[km]
        Densidad = TotalCauces / Area #[km/km2]
        Lcau=ppal[1,-1]/1000.0
        Scau=np.polyfit(ppal[1,::-1],ppal[0],1)[0]*100
        Scue=slope.mean()*100
        Hmin=Elev[-1]; Hmax=Elev[puntto]; Hmean=Elev.mean()
        HCmax=Elev[punto]
        x,y = cu.basin_coordxy(self.structure,self.ncells)
        CentXY = [np.median(x),np.median(y)]
        #Genera un diccionario con las propiedades de la cuenca
        self.GeoParameters={'Area_[km2]': Area,
            'Pend_Cauce_[%]':Scau,
            'Long_Cauce_[km]': Lcau,
            'Pend_Cuenca_[%]': Scue,
            'Long_Cuenca_[km]': Lpma,
            'Hmax_[m]':Hmax,
            'Hmin_[m]':Hmin,
            'Hmean_[m]':Hmean,
            'H_Cauce_Max_[m]':HCmax,
            'Centro_[X]': CentXY[0],
            'Centro_[Y]': CentXY[1],
            'Long_tot_cauces_[km]': TotalCauces,
            'Densidad_drenaje_[km/km2]': Densidad}
        if GetPerim:
            self.GeoParameters.update({'Perimetro_[km]':Perim})
        #Calcula los tiempos de concentracion
        Tiempos={}
        Tc=0.3*(Lcau/(Scue**0.25))**0.75
        Tiempos.update({'US Army': Tc})
        Tc=0.3*(Lcau/((Hmax-Hmin)/Lcau)**0.25)**0.75
        Tiempos.update({'Direccion Carreteras Espana': Tc})
        Tc=(0.02*(Lpma*1000.0)**0.77)/((Scau/100.0)**0.385)/60.0
        Tiempos.update({'Kiprich': Tc})
        Tc=8.157*((Area**0.316)/(((Scau*100)**0.17)*Scue**0.565))
        Tiempos.update({'Campo y Munera': Tc})
        Tc=(4*np.sqrt(Area)+1.5*Lcau)/(0.8*np.sqrt(Hmean))
        Tiempos.update({'Giandotti':Tc})
        Tc=0.4623*(Lcau**0.5)*((Scue/100.0)**(-0.25))
        Tiempos.update({'John Stone': Tc})
        Tc=(Lcau/Area)*(np.sqrt(Area)/Scau)
        Tiempos.update({'Ventura': Tc})
        Tc=0.3*(Lcau/(((HCmax-Hmin)/Lcau)*100)**0.25)**0.75
        Tiempos.update({'Temez': Tc})
        self.Tc=Tiempos
        #Si se habilita la funcion para guardar el ascii de param lo hace
        if rutaParamASC is not None:
            self.__WriteGeoParam__(rutaParamASC)
        # Grafica Tc si se habilita
        if plotTc is True:
            self.Plot_Tc(ruta = rutaTcPlot, figsize=figsize)
    #Funcion para escribir los parametros de la cuenca en un ascii
    def __WriteGeoParam__(self,ruta):
        f=open(ruta,'w')
        f.write('------------------------------------------------------------ \n')
        f.write('Parametros Cuenca \n')
        f.write('------------------------------------------------------------ \n')
        k=self.GeoParameters.keys()
        for i in k:
            v=self.GeoParameters[i]
            if type(v) is not list:
                f.write('%s : %.4f \n' % (i,v))
        f.write('------------------------------------------------------------ \n')
        f.write('Tiempos de concentracion \n')
        f.write('------------------------------------------------------------ \n')
        k=self.Tc.keys()
        for i in k:
            v=self.Tc[i]
            f.write('%s : %.4f \n' % (i,v))
        f.close()

    #Obtiene la envolvente de la cuenca
    def __GetBasinPolygon__(self):
            'Descripcion: obtiene la envolvente de la cuenca, en coordenadas \n'\
            '   x,y, esta informacion luego sirve para plot y para escribir el\n'\
            '   shpfile de la cuenca\n'\
            #Evalua si se cuenta con rasterio en el sistema o no.
            if FlagBasinPolygon:
                #Obtiene mapa raster, propiedades geo y formato para escribir
                Map, Prop = self.Transform_Basin2Map(np.ones(self.ncells))
                tt = [Prop[2], Prop[4].tolist(), 0.0,Prop[1]*Prop[-2] + Prop[3], 0.0, -1*Prop[5].tolist()]
                #Obtiene los shps con la forma de la o las envolventes de cuenca
                Map = Map.T
                mask = Map != -9999.
                c,a,b,f,d,e = tt;tt = __fea__.rasterio.Affine(a,b,c,d,e,f) if (float('.'.join(__fea__.rasterio.__version__.split('.')[:2]))>=0.8) else tt
                shapes = __fea__.shapes(Map, mask=mask, transform=tt)
                #Obtiene el poligono de la cuenca completo
                Shtemp = []
                flag = True
                
                while flag:
                        try:
                                Shtemp.append(next(shapes))
                        except:
                                flag = False
                DicPoly = {}
                for Sh in Shtemp:
                    Coord = Sh[0]['coordinates']
                    Value = int(Sh[1])
                    DicPoly.update({str(Value):{}})
                    for cont,co in enumerate(Coord):
                        DicPoly[str(Value)].update({str(cont):np.array(co).T})
                #Por si hay mas de un poligono al interior
                Tam = 0
                #return DicPoly
                for k in DicPoly['1'].keys():
                    if DicPoly['1'][k].size > Tam:
                        Tam = DicPoly['1'][k].size
                        goodKey = k
                self.Polygon = DicPoly['1'][goodKey]
                return 0
            else:
                return 1

    #Parametros por mapas (distribuidos)
    def GetGeo_Cell_Basics(self):
        'Descripcion: Obtiene: area acumulada, long de celdas, Pendiente \n'\
        '   y Elevacion en cada elemento de la cuenca. \n'\
        '\n'\
        'Parametros\n'\
        '----------\n'\
        'self : no necesita nada es autocontenido.\n'\
        '\n'\
        'Retornos\n'\
        '----------\n'\
        'CellAcum : Cantidad de celdas acumuladas.\n'\
        'CellLong : Longitud de cada una de las celdas [mts].\n'\
        'CellSlope : Pendiente de cada una de las celdas [y/x].\n'\
        'CellHeight : Elevacion de cada una de las celdas [m.s.n.m].\n'\
        #obtiene los parametros basicos por celdas
        acum,longCeld,S0,Elev=cu.basin_basics(self.structure,
            self.DEMvec,self.DIRvec,self.ncells)
        self.CellAcum=acum; self.CellLong=longCeld
        self.CellSlope=S0; self.CellHeight=Elev
        #Obtiene el canal en la cuenca
        self.CellCauce = np.zeros(self.ncells)
        self.CellCauce[self.CellAcum>self.umbral]=1
    def GetGeo_StreamOrder(self, MajorBasins = False, umbral = 100, verbose = False, FirtsOrder = 1):
        'Descripcion: Obtiene el orden de horton para cada celda de \n'\
        '   cada ladera y para las celdas de cada cauce.\n'\
        '\n'\
        'Parametros\n'\
        '----------\n'\
        'self : no necesita nada es autocontenido.\n'\
        'MajorBasins : Obtiene binarios con las sub-cuencas drenando.\n'\
        '   unicamente a cuencas de orden mayor (ej: todas las orden 2 que \n'\
        '   drenan a orden 3 o major).\n'\
        'umbral: Cantidad minima de celdas para considerar canal (aplica cuando\n'\
        '   se obtienen las sub-cuencas mayores.\n'\
        'verbose: Muestra el paso de calculo de cuencas mayores.\n'\
        'FirtsOrder: Primer orden a partir dle cual se analizan ordenes mayores.\n'\
        '\n'\
        'Retornos\n'\
        '----------\n'\
        'CellHorton_Hill : Orden de horton de cada ladera.\n'\
        'CellHorton_Stream : Orden de horton de cada elemento de cauce.\n'\
        #obtiene los parametros basicos por celdas
        cauce,nodos_fin,n_nodos = cu.basin_subbasin_nod(self.structure,self.CellAcum,self.umbral,self.ncells)
        sub_pert,sub_basin = cu.basin_subbasin_find(self.structure,nodos_fin,n_nodos,self.ncells)
        sub_basins = cu.basin_subbasin_cut(n_nodos)
        sub_horton,nod_horton = cu.basin_subbasin_horton(sub_basins,self.ncells,n_nodos)
        self.CellHorton_Hill,sub_basin = cu.basin_subbasin_find(self.structure,nod_horton,n_nodos,self.ncells)
        #Obtiene el canal en la cuenca
        self.CellHorton_Stream = self.CellCauce * self.CellHorton_Hill
        #Obtiene las cuencas mayores
        if MajorBasins:
            pos = np.where(models.control>0)[1]
            X,Y = cu.basin_coordxy(self.structure, self.ncells)
            DictBasins = {}
            for Orden in range(FirtsOrder,self.CellHorton_Hill[-1]):
                #Encuentra cuencas de un orden que drenen a un orden mayor
                pos2 = np.where(self.CellHorton_Stream[pos] == Orden)[0]
                drena = self.ncells - self.structure[0]
                pos3 = np.where(self.CellHorton_Stream[drena[pos[pos2]]] > Orden)[0]
                #Las ubica en un mapa dentro d ela cuenca
                SubCuencas = np.zeros(self.ncells)
                cont = 1
                for x,y in zip(X[pos[pos2[pos3]]], Y[pos[pos2[pos3]]]):
                    #Traza la cuenca
                    cuTemp = SimuBasin(x,y, self.DEM, self.DIR, umbral=umbral)
                    #La pega en una mascara con las cub-cuencas
                    Map, prop = cuTemp.Transform_Basin2Map(np.ones(cuTemp.ncells),)
                    Map[Map == -9999] = 0
                    Var = self.Transform_Map2Basin(Map, prop)
                    Var[Var == 0] = -9999
                    Var[Var == 1] = cont
                    ptemp = np.where(Var == cont)[0]
                    #Lo pega en la mascara de sub-uencas
                    SubCuencas[ptemp] = SubCuencas[ptemp] + Var[ptemp]
                    cont+=1
                #Agrega al diccionario
                DictBasins.update({str(Orden):SubCuencas})
                #Si es verbose muestra en que paso va
                if verbose:
                    print('Sub-cuencas orden '+str(Orden)+' calculadas')
            #Traza la cuenca original para no danar la estructura de guardado
            cuTemp = SimuBasin(X[-1], Y[-1], self.DEM, self.DIR, umbral = umbral)
            #Retorna el diccionario con las sub-cuencas mayore
            return DictBasins
    def GetGeo_IsoChrones(self,Tc,Niter=4):
        'Descripcion: Obtiene el tiempo de viaje aproximado de cada  \n'\
        '   celda a la salida de la cuenca, para eso usa el tiempo de . \n'\
        '   concentracion obtenido por la funcion GetGeo_Parameters . \n'\
        '\n'\
        'Parametros\n'\
        '----------\n'\
        '   self : no necesita nada es autocontenido.\n'\
        '   Tc : Valor escalar de tiempo de concentracion.\n'\
        '   Niter: Cantidad de iteraciones para aproximar vel, defecto 4.\n'\
        '\n'\
        'Retornos\n'\
        '----------\n'\
        'isochrones : Mapa de viaje de cada celda a la salida [hrs].\n'\
        #Calcula la velocidad adecuada para que el tiempo coincida
        acum,longCeld,S0,Elev=cu.basin_basics(self.structure,
            self.DEM,self.DIR,cu.ncols,cu.nrows,self.ncells)
        rangos=[50,25,1]
        for i in range(Niter):
            times=[]
            for r in rangos:
                speed = r*S0**(0.5)
                time = cu.basin_time_to_out(self.structure,
                    longCeld,speed,self.ncells)/3600.0
                times.append(time[np.isfinite(time)].mean())
            for j in range(2):
                if Tc>times[j] and Tc<times[j+1]:
                    rangos=[rangos[j],(rangos[j]+rangos[j+1])/2.0,
                        rangos[j+1]]
        #Calcula los intervalos
        intervalos=np.arange(0,np.ceil(time.max())+1,
            np.ceil(time.max())/10.0)
        timeC=np.zeros(self.ncells)
        tamano=[]
        for i,j in zip(intervalos[:-1],intervalos[1:]):
            timeC[(time>=i) & (time<j)]=(i+j)/2.0
            tamano.append(timeC[timeC==(i+j)/2.0].shape[0])
        tamano=np.array(tamano)
        aportes=(tamano/float(self.ncells))*((self.ncells*cu.dxp**2)/1e6)
        self.CellTravelTime=time

    def GetGeo_WidthFunction(self, binsC = 50, binsN = 50,
        ruta = None, Npos = 10000, **kwargs):
        'Descripcion: Obtiene la funcion de ancho de la cuenca  \n'\
        '   Entrega como resultado la distancia de cada elemento a la salida. \n'\
        '   y grafica la funcion de ancho. \n'\
        '\n'\
        'Parametros\n'\
        '----------\n'\
        '   self : no necesita nada es autocontenido.\n'\
        '   Previo a esta funcion se debe usar: cu.Get_Cell_Basics.\n'\
        '\n'\
        'Param Opcionales\n'\
        '   binsC: intervalos de clase para la cuenca.\n'\
        '   binsN: intervalos de clase para la red.\n'\
        '   ruta: Ruta donde se guarda la figura.\n'\
        '   Npos: Cantidad de puntos muestreados.\n'\
        'Param Kwargs.\n'\
        '   show: Muestra o no la figura (True).\n'\
        '   ccolor: Color del histograma de la cuenca.\n'\
        '   ncolor: Color del histograma de la red hidrica.\n'\
        '   figsize: Tamano del plot (8,5).\n'\
        'Retornos\n'\
        '----------\n'\
        'self.CellDist2Out : Distancia de cada elemento a la salida guardado en objeto\n'\
        'Figura de la funcion de ancho.\n'\
        #Parametros kwargs
        show = kwargs.get("show", True)
        ccolor = kwargs.get("ccolor",'b')
        ncolor = kwargs.get("ncolor",'r')
        figsize = kwargs.get("figsize", (8,5))
        #Calcula la velocidad adecuada para que el tiempo coincida
        self.CellDist2Out = cu.basin_time_to_out(self.structure,
            self.CellLong, np.ones(self.ncells), self.ncells)
        #Cantidad de puntos
        pos = np.random.choice(self.ncells, Npos)
        #Calcula la funcion de ancho para la cuenca y para la red hidrica
        hc,bc = np.histogram(self.CellDist2Out[pos], bins = binsC)
        ViajeCauce = self.CellDist2Out*self.CellCauce
        hn,bn = np.histogram(ViajeCauce[ViajeCauce>0], bins = binsN)
        #Variables del elemento
        self.width_hits = hn
        self.width_distances = bn[:-1]/1000.
        #Estandariza en terminos de probabilidad
        hc = hc.astype(float); hc = hc / hc.sum()
        hn = hn.astype(float); hn = hn / hn.sum()
        #Grafica
        fig = pl.figure(figsize=figsize)
        ax = fig.add_subplot(111)
        ax.plot(bc[:-1]/1000.0, hc, ccolor, lw = 3, label = 'Basin Width')
        ax.plot(bn[:-1]/1000.0, hn, ncolor,lw = 3, label = 'Network Width')
        ax.grid(True)
        ax.set_xlabel('Distance to Outlet [$km$]', size = 16)
        ax.set_ylabel('PDF', size = 16)
        ax.fill_between(bc[:-1]/1000.0, hc, color = ccolor, alpha = 0.3)
        ax.fill_between(bn[:-1]/1000.0, hn, color = ncolor, alpha = 0.3)
        ax.set_ylim(0, np.max([hc.max(), hn.max()])+0.01)
        ax.tick_params(labelsize = 15)
        ax.legend(loc = 0)
        ax2 = ax.twinx()
        ax2.plot(bc[:-1]/1000.0, hc.cumsum(), ccolor, lw = 3, ls = '--', label = 'Basin Width')
        ax2.plot(bn[:-1]/1000.0, hn.cumsum(), ncolor, lw = 3, ls = '--',label = 'Network Width')
        ax2.tick_params(labelsize = 15)
        ax2.set_ylabel('CDF', size = 16)
        #Guarda la figura
        if ruta is not None:
            pl.savefig(ruta, bbox_inches='tight',pad_inches = 0.25)
        if show:
            pl.show()
        #Retorna ejes de manipulacion
        return ax,ax2

    def GetGeo_Ppal_Hipsometric(self,umbral=1000,
        intervals = 30):
        'Descripcion: Calcula y grafica la curva hipsometrica de\n'\
        '   la cuenca.\n'\
        '\n'\
        'Parametros\n'\
        '----------\n'\
        'umbral : cantidad minima de celdas para el trazado.\n'\
        'intervals: Cantidad de intervalos en los cuales se haran .\n'\
        '   los muestreos de la curva hipsometrica.\n'\
        '\n'\
        'Retornos\n'\
        '----------\n'\
        'Hipso : Curva hipsometrica en formato vectorial, contiene.\n'\
        '   - Curva sobre cauce ppal.\n'\
        '   - Curva como histograma de la cuenca.\n'\
        'Figura: Figura de ambas curvas.\n'\
        # Obtiene los nodos de la red hidrica
        cauce,nodos,trazado,n_nodos,n_cauce = cu.basin_stream_nod(
            self.structure,
            self.CellAcum,
            umbral,
            self.ncells)
        # Obtiene el cauce ppal
        ppal_nceldas,punto = cu.basin_ppalstream_find(
            self.structure,
            nodos,
            self.CellLong,
            self.CellHeight,
            self.ncells)
        self.ppal_stream = cu.basin_ppalstream_cut(ppal_nceldas,
            self.ncells)
        #Corrige el perfil del cauce ppal
        self.ppal_stream[0],self.ppal_slope = __ModifyElevErode__(self.ppal_stream[0])
        # Obtiene la curva hipsometrica
        self.hipso_ppal, self.hipso_basin = cu.basin_ppal_hipsometric(
            self.structure,
            self.CellHeight,
            punto,
            intervals,
            ppal_nceldas)
        self.hipso_ppal[1],self.hipso_ppal_slope = __ModifyElevErode__(self.hipso_ppal[1])
    def GetGeo_IT(self):
        'Descripcion: Calcula el indice topografico para cada celda (Beven)\n'\
        '   Internamente calcula el area para cada elemento y la pendiente en radianes.\n'\
        '\n'\
        'Parametros\n'\
        '----------\n'\
        'self : no necesita nada es autocontenido.\n'\
        '\n'\
        'Retornos\n'\
        '----------\n'\
        'IT : Indice topografico adimensional, a mayor valor se supone un suelo mas humedo.\n'\
        #Obtiene el area
        acum = cu.basin_acum(self.structure, self.ncells)
        #Obtiene la pendiente en radianes
        acum,longCeld,slope,Elev=cu.basin_basics(self.structure,
            self.DEMvec,self.DIRvec,self.ncells)
        slope = np.arctan(slope)
        slope[slope == 0] = 0.0001
        return np.log((acum*cu.dxp) / np.tan(slope))

    def GetGeo_HAND(self,umbral=1000):
        'Descripcion: Calcula Height Above the Nearest Drainage (HAND) \n'\
        '   y Horizontal Distance to the Nearest Drainage (HDND) (Renno, 2008). \n'\
        '\n'\
        'Parametros\n'\
        '----------\n'\
        'self : no necesita nada es autocontenido.\n'\
        'umbral : cantidad minima de celdas para el trazado.\n'\
        '\n'\
        'Retornos\n'\
        '----------\n'\
        'HAND : Elevacion sobre la red de drenaje cercana [mts].\n'\
        'HDND : Distancia horizontal a la red de drenaje cercana [mts].\n'\
        #obtiene los parametros basicos por celdas
        acum,longCeld,S0,Elev=cu.basin_basics(self.structure,
            self.DEMvec,self.DIRvec,self.ncells)
        cauce,nodos,trazado,n_nodos,n_cauce = cu.basin_stream_nod(
            self.structure,acum,umbral,self.ncells)
        hand,hdnd,hand_destiny = cu.geo_hand(self.structure,Elev,longCeld,cauce,self.ncells)
        handC=np.zeros(self.ncells)
        handC[hand<5.3]=1
        handC[(hand>=5.3) & (hand<=15.0)]=2
        handC[(hand>15.0) & (S0<0.076)]=4
        handC[(hand>15.0) & (S0>=0.076)]=3
        self.CellHAND=hand
        self.CellHAND_class=handC
        self.CellHDND=hdnd
        self.CellHAND_drainCell=hand_destiny

    def GetGeo_Sections(self, NumCeldas = 6):
        'Descripcion: Obtiene secciones transversales a traves de todos.\n'\
        '   los elementos de la red de drenaje, las secciones se obtienen\n'\
        '   en la direccion perpendicular al flujo, es decir si el mapa de\n'\
        '   direcciones indica en una celda la direccion norte, las secciones\n'\
        '   se obtienen en lsa direcciones oriente, occidente, con NumCeldas\n'\
        '   a cada lado de la celda tipo cauce.\n'\
        '\n'\
        'Parametros\n'\
        '----------\n'\
        'NumCeldas: Cantidad de celdas para elaborar secciones a ambos lados.\n'\
        '\n'\
        'Retornos\n'\
        '----------\n'\
        'self.Sections : Secciones a traves de los elementos del cauce.\n'\
        '   su tamano es [NumCeldas*2 + 1, self.ncells]\n'\
        #Obtiene mapa de cauces
        self.GetGeo_Cell_Basics()
        #Obtiene vector de direcciones
        directions = self.Transform_Map2Basin(self.DIR,
            [cu.ncols, cu.nrows, cu.xll, cu.yll, cu.dx, cu.dy ,0.0])
        #Obtiene las secciones
        self.Sections, self.Sections_Cells = cu.basin_stream_sections(self.structure,
            self.CellCauce, directions, self.DEM, NumCeldas,
            self.ncells, cu.ncols, cu.nrows)

    #------------------------------------------------------
    # Subrutinas para el calculo de extremos mediante hidrografa unitaria sintetica
    #------------------------------------------------------
    def GetHU_Snyder(self,Area,Tc,Cp=0.8,Fc=2.9,pequena='si'):
        'Descripcion: Obtiene HU segun Snyder.\n'\
        '\n'\
        'Parametros\n'\
        '----------\n'\
        'Area: area de la cuenca en km2.\n'\
        'Tc: Tiempo de concentracion en Horas.\n'\
        'Cp : Factor de escalamiento [0.8].\n'\
        'Fc : Factor de escalamiento [2.9].\n'\
        '\n'\
        'Retornos\n'\
        '----------\n'\
        'Tiempo : Tiempo total de duracion [M].\n'\
        'Q : Caudal unitario.\n'\
        'HU : Reglas del hidrograma unitario.\n'\
        #Calcula parametros para la hidrografa
        Tr=0.6*Tc # [h]
        T=0.1*Tc # [h]
        Ts=Tr/5.5 # [h]
        up=Cp*640/(Tr+(T-Ts)/4) #pie3/s/mi2/Pulg
        Up=up*0.386*Area
        Up=Up*(0.3048)**3.0/25.4
        Tp=T/2+Tr # horas
        Tb=3+3*Tr/24 # dias
        W50=770/up**1.08
        W75=440/up**1.08
        # Factor de multiplicacion
        tb=Tp*Fc #horas
        #obtiene las coordenadas del hidrograma
        HU=[[0,0]]
        At=(Tp-W50/3)*60; Aq=0.5*Up; HU.append([At,Aq])
        Bt=(Tp-W75/3)*60; Bq=0.75*Up; HU.append([Bt,Bq])
        Tpt=Tp*60; Tpq=Up; HU.append([Tpt,Tpq])
        Dt=(Tp+2*W75/3)*60; Dq=0.75*Up; HU.append([Dt,Dq])
        Et=(Tp+2*W50/3)*60; Eq=0.5*Up; HU.append([Et,Eq])
        if pequena=='si':
            Tbt=tb*60
        else:
            Tbt=Tb*60
        HU.append([Tbt,0])
        HU=np.array(HU).T
        #Obtiene los intervalod de tiempo
        Dt=Tc*0.1*60 #min
        Tiempo=np.arange(0,Tbt+Dt,Dt)
        #Obtiene la HU para los intervalos de tiempo necesarios
        Q=np.zeros(Tiempo.size)
        for cont,t in enumerate(Tiempo):
            for h1,h2 in zip(HU.T[:-1],HU.T[1:]):
                if t>h1[0] and t<h2[0]:
                    p=np.polyfit([h1[0],h2[0]],[h1[1],h2[1]],1)
                    Q[cont]=t*p[0]+p[1]
        AreaBajoCurva=Q.sum()*Dt*60/1000.0
        Diferencia=100*((Area-AreaBajoCurva)/Area)
        #Devuelve el hu y la diferencia
        return Tiempo,Q,HU,Diferencia
    def GetHU_Williams(self,Area,LongCuenca,PendCauce,Tc):
        'Descripcion: Obtiene HU segun Williams.\n'\
        '\n'\
        'Parametros\n'\
        '----------\n'\
        'Area: area de la cuenca en km2.\n'\
        'LongCuenca: Longitu de la cuenca en km.\n'\
        'PendCauce: Pendiente del cauce en m/m.\n'\
        'Tc: Tiempo de concentracion en Horas.\n'\
        '\n'\
        'Retornos\n'\
        '----------\n'\
        'Tiempo : Tiempo total de duracion [M].\n'\
        'Q : Caudal unitario.\n'\
        'HU : Reglas del hidrograma unitario.\n'\
        #Cambio de unidades en area y longitu de cuenca
        Area=Area*0.386; LongCuenca=LongCuenca*0.62
        #Calculo parametros hidrograma
        w=(Area)/(LongCuenca) #millas
        Lw=(LongCuenca)/w #adim
        tp=4.63*(Area**0.422)*(PendCauce**-0.46)*(Lw**0.133)
        k=27*(Area**0.231)*(PendCauce**-0.777)*(Lw**0.124)
        ktp=k/tp
        n=1+((0.5/ktp) + ((0.25/ktp**2) + (1.0/ktp))**0.5)**2.0
        B=262.9*np.log(n)+35.723
        fc=0.3048**3
        to=tp*(1+1/(n-1)**0.5); t1=to+2*k
        HU=[]; HU.append([tp,B*0.0394*(Area/tp)])
        HU.append([to,HU[0][1]*(to/tp)**(n-1)*np.exp((1-n)*((to/tp)-1))])
        HU.append([t1,HU[1][1]*np.exp((to-t1)/k)])
        HU=np.array(HU).T; HU[1,:]=HU[1,:]*fc
        #Obtiene los intervalos de tiempo
        Dt=Tc*0.1*60 #min
        Tiempo=np.arange(0,3*t1*60+Dt,Dt)/60.0
        ttp=(Tiempo)/tp
        Q=np.zeros(Tiempo.size)
        for cont,t in enumerate(zip(ttp,Tiempo)):
            if t[1]<to:
                Q[cont]=HU[1,0]*(t[0]**(n-1))*np.exp((1-n)*(t[0]-1))
            if t[1]>t1:
                Q[cont]=HU[1,2]*np.exp((t1-t[1])/(3*k))
            if t[1]<t1 and t[1]>to:
                Q[cont]=HU[1,1]*np.exp((to-t[1])/k)
        return Tiempo*60,Q,HU
    def GetHU_SCS(self,AreaCuenca,Tc,N=25):
        'Descripcion: Obtiene HU segun SCS.\n'\
        '\n'\
        'Parametros\n'\
        '----------\n'\
        'AreaCuenca: area de la cuenca en km2.\n'\
        'Tc: Tiempo de concentracion en Horas.\n'\
        '\n'\
        'Retornos\n'\
        '----------\n'\
        'Tiempo : Tiempo total de duracion [M].\n'\
        'Q : Caudal unitario.\n'\
        'HU : Reglas del hidrograma unitario.\n'\
        #Parametros fijos
        ttp=np.array([0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.1,
        1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2.0,2.2,2.4,2.6,2.8,3.0,3.2,
        3.4,3.6,3.8,4.0,4.5,5.0])
        UUp=np.array([0.0,0.03,0.1,0.19,0.31,0.47,0.66,0.82,0.93,0.99,1.0,0.99,0.93,0.86,
        0.78,0.68,0.56,0.46,0.39,0.33,0.28,0.21,0.15,0.11,0.08,0.06,0.04,0.03,0.02,0.02,0.01,
        0.01,0.0])
        #Parametros de la hidrografa
        Trezago=(3.0/5.0)*float(Tc)
        Tlluvia=0.133*Tc
        Tpico=Tlluvia/2.0+Trezago
        Upa=484.0*AreaCuenca*0.386/Tpico
        Upm=Upa*0.3048**3/25.4
        #Obtiene el HU
        HU=[]
        for t,u in zip(ttp,UUp):
            uTemp=u*Upa
            HU.append([t*Tpico,uTemp*0.3048**3/25.4])
        HU=np.array(HU)
        #Obtiene los intervalod de tiempo
        Dt=Tc*0.1*60 #min
        #Obtiene la HU para los intervalos de tiempo necesarios
        Q=[]; flag=True
        Tiempo=[];t=0.0
        for i in range(N):
            for h1,h2 in zip(HU[:-1],HU[1:]):
                if t>=h1[0] and t<h2[0]:
                    p=np.polyfit([h1[0],h2[0]],[h1[1],h2[1]],1)
                    Q.append(t*p[0]+p[1])
                    Tiempo.append(t*60.0)
                    t+=Dt/60.0
        return np.array(Tiempo),np.array(Q),HU.T
    # Funcion Tormenta de diseno
    def GetHU_DesingStorm(self,IntTr,Dur,CN=70,plot='no',ruta=None,Tr=None,
        CurvaHuff=np.array([0.18,0.47,0.65,0.74,0.80,0.85,0.89,0.92,0.94,1.0])):
        'Descripcion: Obtiene la tormenta de diseno a partir de una tormenta\n'\
        '\n'\
        'Parametros\n'\
        '----------\n'\
        'IntTr : Intensidades de precipitacion para diferentes Tr [Ntr].\n'\
        'Dur : Duracion del evento de tormenta.\n'\
        'CN: Numero de curva para el calculo de abstracciones.\n'\
        '\n'\
        'Retornos\n'\
        '----------\n'\
        'Qtr : Hidrografa para cada periodo de retorno [M, Tr].\n'\
        'QmaxTr : Caudal maximo para cada periodo de retorno [Tr].\n'\
        'Tiempo : Tiempo total de duracion [M].\n'\
        #Calcula la lamina para cada periodo de retorno
        lluviaTrAcum=[]
        for P in IntTr*Dur:
            lluviaTrAcum.append([c*P for c in CurvaHuff])
        lluviaTrAcum=np.array(lluviaTrAcum)
        lluviaTr=[]
        for l in lluviaTrAcum:
            a=np.copy(l)
            a[1:]=l[1:]-l[:-1]
            lluviaTr.append(a)
        lluviaTr=np.array(lluviaTr)
        #si se da el valor de CN lo calcula
        if CN is not None:
            S=(1000.0/float(CN)-10.0)*25.4; Ia=0.2*S
            lluviaTrTemp=[]
            for l in lluviaTrAcum:
                L=[]
                for i in l:
                    if i>0:
                        L.append(((i-Ia)**2)/(i-Ia+S))
                    else:
                        L.append(0)
                lluviaTrTemp.append(L)
            lluviaTrTemp=np.array(lluviaTrTemp)
            lluviaTrEfect=lluviaTrTemp[:,1:]-lluviaTrTemp[:,:-1]
            lluviaTrEfect[lluviaTrEfect<0]=0
            lluviaTrEfect=np.insert(lluviaTrEfect,0,lluviaTrTemp[:,0],axis=1)
        if plot=='si':
            #if Tr==None or len(Tr)!=IntTr.size():
            #   Lista=[2.33,5,10,25,50,100,500,1000]
            #   Tr=[l for l in Lista[:IntTr.size]]
            fig=pl.figure(edgecolor='w',facecolor='w')
            ax=fig.add_subplot(111)
            X=np.linspace(0,Dur,lluviaTrEfect.shape[1])
            Grosor=np.arange(0.5,4,0.2)
            for l,le,t,g in zip(lluviaTr,lluviaTrEfect,Tr,Grosor):
                ax.plot(X,l,c='b',lw=g,label=str(t))
            if CN is not None:
                for l,le,t,g in zip(lluviaTr,lluviaTrEfect,Tr,Grosor):
                    ax.plot(X,le,c='r',lw=g,label=str(t))
            ax.set_xlabel('Tiempo $[h]$',size=16)
            ax.set_ylabel('Precipitacion $[mm]$',size=16)
            ax.grid(True)
            ax.tick_params(labelsize=15)
            pl.legend(loc=0,ncol=2)
            if ruta is not None:
                pl.savefig(ruta,bbox_inches='tight')
            pl.show()
            if CN is not None:
                return lluviaTr,np.array(lluviaTrEfect),S
            else:
                return lluviaTr
    #Convolucion de la tormenta de diseno
    def GetHU_Convolution(self,Tiempo,Qhu,lluvEfect):
        'Descripcion: Realia la convolucion de los hidrogramas unitarios \n'\
        '\n'\
        'Parametros\n'\
        '----------\n'\
        'Tiempo : Vector de la duracion de la hidrografa unitaria [Nt].\n'\
        'Qhu : Valores del caudal unitario de la hidrografa unitaria [Nt].\n'\
        'lluvEfect : Lluvia efectiva en un evento de tormenta cualquiera [Nt, Tr].\n'\
        '\n'\
        'Retornos\n'\
        '----------\n'\
        'Qtr : Hidrografa para cada periodo de retorno [M, Tr].\n'\
        'QmaxTr : Caudal maximo para cada periodo de retorno [Tr].\n'\
        'Tiempo : Tiempo total de duracion [M].\n'\
        #Calcula el caudal convulucionado
        Qtr=[]; QmaxTr=[]
        for le in lluvEfect:
            Qi=np.array([i*Qhu for i in le])
            Qconv=np.zeros(Qhu.size+lluvEfect.shape[1]-1)
            for cont,q in enumerate(Qi):
                Qconv[cont:cont+Qhu.size]=Qconv[cont:cont+Qhu.size]+q
            Qtr.append(Qconv)
            QmaxTr.append(Qconv.max())
        #Calcula el tiempo para toda la hidrografa
        Dt=Tiempo[1]-Tiempo[0]; Tlast=Tiempo[-1]
        T=list(Tiempo)
        for i in range(1,lluvEfect.shape[1]):
            T.append(Tlast+i*Dt)
        return np.array(Qtr), np.array(QmaxTr),np.array(T)
    #Grafica los hidrogramas sinteticos
    def PlotHU_Synthetic(self,DictHU,ruta = None):
        fig=pl.figure(edgecolor='w',facecolor='w')
        ax=fig.add_subplot(111)
        colors = ['b','r','k','g','m']
        for co,k in enumerate(DictHU.keys()):
            ax.plot(DictHU[k]['time'],
                DictHU[k]['HU'],
                c=colors[co],lw=1.5,label=k)
        ax.grid(True)
        ax.set_xlabel('Tiempo $[min]$',size=14)
        ax.set_ylabel('HU $[m^3/seg/mm]$',size=14)
        ax.legend(loc=0)
        if ruta is not None:
            pl.savefig(ruta,bbox_inches='tight')
        pl.show()

    #------------------------------------------------------
    # Trabajo con mapas externos y variables fisicas
    #------------------------------------------------------
    def Transform_Map2Basin(self,Map,MapProp):
        'Descripcion: A partir de un mapa leido obtiene un vector \n'\
        '   con la forma de la cuenca, el cual luego puede ser agregado a esta. \n'\
        '\n'\
        'Parametros\n'\
        '----------\n'\
        'self : Inicia las variables vacias.\n'\
        'Map : Matriz con la informacion del mapa.\n'\
        'MapProp : Propiedades del mapa.\n'\
        '   1. Ncols Mapa.\n'\
        '   2. Nrows Mapa.\n'\
        '   3. Xll Mapa.\n'\
        '   4. Yll Mapa.\n'\
        '   5. dx Mapa.\n'\
        '\n'\
        'Retornos\n'\
        '----------\n'\
        'vecMap : Vector conla informacion del mapa al interio de la cuenca.\n'\
        #Comienza le codifgo
        vec = cu.basin_map2basin(self.structure,
            Map,MapProp[2],MapProp[3],MapProp[4],MapProp[5],
            cu.nodata,
            self.ncells,
            MapProp[0],MapProp[1])
        return vec
    def Transform_Basin2Map(self, BasinVar, ruta = None, DriverFormat='GTiff',
        EPSG=4326):
        'Descripcion: A partir de un vector con propiedades de la cuenca en celdas\n'\
        '   obtiene un mapa (matriz) con las propiedades del DEM, este puede ser escrito \n'\
        '   en algun formato legible por GDAL. \n'\
        '\n'\
        'Parametros\n'\
        '----------\n'\
        'self : la cuenca misma.\n'\
        'BasinVar : Vector con la variable de la cuenca a escribir [ncells].\n'\
        '\n'\
        'Retornos\n'\
        '----------\n'\
        'Map : Matriz con la variable de la cuenca, donde no hay cuenca es wmf.cu.nodata .\n'\
        '   el tamano de la matriz es igual a DEM.shape.\n'\
        # Convierte la variable a mapa
        map_ncols,map_nrows = cu.basin_2map_find(self.structure,self.ncells)
        M,mxll,myll = cu.basin_2map(self.structure, BasinVar,
            map_ncols, map_nrows, self.ncells)
        # Si exporta el mapa lo guarda si no simplemente devuelve la matriz
        if ruta is not None:
            Save_Array2Raster(M, [map_ncols,map_nrows,mxll,myll,cu.dx,cu.nodata],
                ruta = ruta, EPSG = EPSG, Format = DriverFormat)
        return M, [map_ncols,map_nrows,mxll,myll,cu.dx,cu.dy,cu.nodata]

    def Transform_Hills2Basin(self,HillsMap):
        'Descripcion: A partir de un vector con propiedades de las laderas\n'\
        '   obtiene un vector con las propiedades por celda, ojo estas \n'\
        '   quedan con las formas de las laderas y la variable queda agregada. \n'\
        '\n'\
        'Parametros\n'\
        '----------\n'\
        'self : la cuenca misma.\n'\
        'MapHills : Vector con las variables por laderas [nhills].\n'\
        '\n'\
        'Retornos\n'\
        '----------\n'\
        'CellMap : Vector con la variable agregada por laderas, pero .\n'\
        '   pasada a celdas.\n'\
        #Genera el mapa de basin vacio
        CellMap = np.ones(self.ncells)
        #itera por la cantidad de elementos y les va asignando
        for i,k in enumerate(HillsMap[::-1]):
            CellMap[self.hills_own==i+1] = k
        return CellMap
    def Transform_Basin2Hills(self,CellMap,mask=None,SumMeanMax=0):
        'Descripcion: A partir de un vector tipo Basin obtiene un\n'\
        '   vector del tipo laderas, en donde las propiedades se \n'\
        '   agregan para cada ladera. \n'\
        '\n'\
        'Parametros\n'\
        '----------\n'\
        'self : la cuenca misma.\n'\
        'CellMap : Vector con las propiedades por celdas [ncells].\n'\
        'mask : Celdas sobre las cuales se agrega la variable (1), y\n'\
        '   sobre las que no (0).\n'\
        'SumMeanMax : si la variable sera agregada como un promedio (0)\n'\
        '   o como una suma (1), o como el maximo valor (2).\n'\
        '\n'\
        'Retornos\n'\
        '----------\n'\
        'HillsMap : Vector con las prop agregadas a laderas .\n'\
        #Si hay mascara la tiene en cuenta
        if mask is not None:
            Ma = np.zeros(self.ncells)
            if type(mask) is float or type(mask) is int:
                Ma[CellMap==mask] = 1
            elif type(mask) is np.ndarray:
                Ma = np.copy(mask)
        else:
            Ma = np.ones(self.ncells)
        #Pasa el mapa de celdas a mapa de laderas
        HillsMap = cu.basin_subbasin_map2subbasin(self.hills_own,
            CellMap, self.nhills, Ma, SumMeanMax, self.ncells)
        return HillsMap


    def Transform_Basin2Polygon(self, Vector,):
        'Descripcion: convierte una variable de topologia de la cuenca en varios poligonos\n'\
        '   cada poligono corresponde al numero de esa variable\n'\
        'parametros\n'\
        '----------\n'\
        'Vector: Vector con la topologia de la cuenca que contiene los datos a transformar'
        'Retorna\n'\
        '----------\n'\
        'DicPoly: Diccionario donde el key indica el poligono encontrado y contiene las coord'\
        #Obtiene mapa raster, propiedades geo y formato para escribir
        Map, Prop = self.Transform_Basin2Map(Vector)
        tt = [Prop[2], Prop[4].tolist(), 0.0,
            Prop[1]*Prop[-2] + Prop[3], 0.0, -1*Prop[5].tolist()]
        #Obtiene los shps con la forma de la o las envolventes de cuenca
        Map = Map.T
        mask = Map != -9999.
        shapes = __fea__.shapes(Map, mask=mask, transform=tt)
        #Obtiene el poligono de la cuenca completo
        DicPoly = {}
        for Sh in shapes:
            Coord = Sh[0]['coordinates']
            Value = int(Sh[1])
            DicPoly.update({str(Value):{}})
            for cont,co in enumerate(Coord):
                DicPoly[str(Value)].update({str(cont):np.array(co).T})
        return DicPoly

    #------------------------------------------------------
    # Trabajo con datos puntuales puntos
    #------------------------------------------------------
    def Points_Points2Stream(self,coordXY,ids):
        'Descripcion: toma las coordenadas de puntos XY y las mueve\n'\
        '   hacia los cauces de la cuenca\n'\
        '\n'\
        'Parametros\n'\
        '----------\n'\
        'coordXY : coordenadas de los puntos [2,Ncoord].\n'\
        'ids : Identificacion de los puntos que se van a mover.\n'\
        '\n'\
        'Retornos\n'\
        '----------\n'\
        'coordXYNew : Coordenadas transportadas a los cauces.\n'\
        'basin_pts : Vector con la forma de la cuenca con los puntos.\n'\
        '   donde quedaron localizadas las coordenadas.\n'\
        'idsOrdered : El orden en que quedaron los Ids ingresados.\n'\
        'Ver Tambien\n'\
        '----------\n'\
        'Points_Point2Basin\n'\
        #Obtiene el cauce
        self.GetGeo_Cell_Basics()
        #modifica los puntos
        res_coord,basin_pts,xy_new = cu.basin_stream_point2stream(
            self.structure,
            self.CellCauce,
            ids,
            coordXY,
            coordXY.shape[1],
            self.ncells)
        return xy_new, basin_pts, basin_pts[basin_pts!=0]
    def Points_Points2Basin(self,coordXY,ids):
        'Descripcion: toma las coordenadas de puntos XY y las pone\n'\
        '   en la cuenca, no las mueve hacia los cuaces\n'\
        '\n'\
        'Parametros\n'\
        '----------\n'\
        'coordXY : coordenadas de los puntos [2,Ncoord].\n'\
        'ids : Identificacion de los puntos que se van a mover.\n'\
        '\n'\
        'Retornos\n'\
        '----------\n'\
        'basin_pts : Vector con la forma de la cuenca con los puntos.\n'\
        '   donde quedaron localizadas las coordenadas.\n'\
        'idsOrdered : El orden en que quedaron los Ids ingresados.\n'\
        'Ver Tambien\n'\
        '----------\n'\
        'Points_Point2Stream\n'\
        #Obtiene el cauce
        self.GetGeo_Cell_Basics()
        #modifica los puntos
        res_coord,basin_pts = cu.basin_point2var(
            self.structure,
            ids,
            coordXY,
            coordXY.shape[1],
            self.ncells)
        return basin_pts,basin_pts[basin_pts!=0]
    #------------------------------------------------------
    # Caudales de largo plazo y regionalizacion
    #------------------------------------------------------
    #Caudal de largo plazo
    def GetQ_Balance(self,Precipitation, Tipo_ETR = 1, mu_choud = 1.37):
        'Descripcion: Calcula el caudal medio por balance de largo plazo\n'\
        '   para ello requiere conocer la precipitacion y el metodo de\n'\
        '   estimacion de la evaporacion.\n'\
        '\n'\
        'Parametros\n'\
        '----------\n'\
        'self : no necesita nada es autocontenido.\n'\
        'Precipitation : Cantidad anual de lluvia, escalar o vector.\n'\
        'Elevacion : Elevacion en cada punto de la cuenca.\n'\
        'Tipo_ETR : Tipo de ecuacion para calcular la evaporacion.\n'\
        '   -1. Turc.\n'\
        '   -2. Cenicafe Budyko.\n'\
        '   -3. Choundry.\n'\
        '   Defecto: 1.\n'\
        '\n'\
        'Retornos\n'\
        '----------\n'\
        'self.CellQmed : Caudal medio calculado para toda la cuenca.\n'\
        'self.CellETR : ETR calculada para toda la cuenca.\n'\
        #Calcula las propiedades de la cuenca
        self.GetGeo_Cell_Basics()
        #Determina si la precipitacion es un vector o un escalar
        if type(Precipitation) is int or type(Precipitation) is float:
            precip = np.ones(self.ncells)*Precipitation
        elif type(Precipitation) is np.ndarray:
            precip = Precipitation
        #Calcula el qmed
        self.CellQmed,self.CellETR = cu.basin_qmed(
            self.structure,
            self.CellHeight,
            precip,
            Tipo_ETR,
            mu_choud,
            self.ncells,)

    #Caudales extremos
    def GetQ_Max(self,Qmed,Coef=[6.71, 3.29], Expo= [0.82, 0.64], Cv = 0.5,
        Tr = [2.33, 5, 10, 25, 50, 100], Dist = 'gumbel', metodo = 'poveda',
        Expo2 = [0.7745, 0.4608]):
        'Descripcion: Calcula el caudal maximo para diferentes\n'\
        '   periodos de retorno. Recibe como entrada el caudal medio \n'\
        '   calculado con GetQ_Balance.\n'\
        '   el calculo se hace a partir de la ecuacion:.\n'\
        '   MedMax = Coef[0] * Qmed ** Exp[0] .\n'\
        '   DesMax = Coef[1] * Qmed ** Exp[1] .\n'\
        '   Qmax = MedMax + K(Tr) * DesMax.\n'\
        '   Donde K es un coeficiente que depende de Tr y f(x).\n'\
        '\n'\
        'Parametros\n'\
        '----------\n'\
        'Qmed : Caudal medio calculado con GetQ_Balance.\n'\
        'Coef : Coeficientes para media y desviacion [6.71, 3.29].\n'\
        'Expo : Exponentes para media y desviacion [0.82, 0.64].\n'\
        'Cv : Coeficiente de escalamiento [0.5].\n'\
        'Tr : Periodo de retorno [2.33, 5, 10, 25, 50, 100].\n'\
        'Dist : Distrbucion [gumbel/lognorm].\n'\
        'metodo: Poveda (u = cQA^a) o atlas (u = C(P-E)^a A^b), en este segundo caso Qmed = P-E.\n'\
        '\n'\
        'Retornos\n'\
        '----------\n'\
        'CellQmax : Caudal maximo para los diferentes periodos.\n'\
        # Calcula la media y desviacion
        if metodo == 'poveda':
            MedMax = Coef[0] * Qmed ** Expo[0]
            DesMax = Coef[1] * Qmed ** Expo[1]
        elif metodo == 'atlas':
            Acum = cu.basin_acum(self.structure, self.ncells)
            Acum = Acum * (cu.dxp**2)/1e6
            MedMax = Coef[0] * (Qmed ** Expo[0]) * (Acum**Expo2[0])
            DesMax = Coef[1] * (Qmed ** Expo[1]) * (Acum**Expo2[1])
        #Itera para todos los periodos de retorno
        Qmax=[]
        for t in Tr:
            #Calcula k
            if Dist is 'gumbel':
                k=-1*(0.45+0.78*np.log(-1*np.log(1-1/float(t))))
            elif Dist is 'lognorm':
                Ztr=norm.ppf(1-1/float(t))
                k=(np.exp(Ztr*np.sqrt(np.log(1+Cv**2))-0.5*np.log(1-Cv**2))-1)/Cv
            #Calcula el caudal maximo
            Qmax.append(list(MedMax+k*DesMax))
        return np.array(Qmax)

    def GetQ_Min(self,Qmed,Coef=[0.4168, 0.2], Expo= [1.058, 0.98], Cv = 0.5,
        Tr = [2.33, 5, 10, 25, 50, 100], Dist = 'gumbel'):
        'Descripcion: Calcula el caudal minimo para diferentes\n'\
        '   periodos de retorno. Recibe como entrada el caudal medio \n'\
        '   calculado con GetQ_Balance.\n'\
        '   el calculo se hace a partir de la ecuacion:.\n'\
        '   MedMin = Coef[0] * Qmed ** Exp[0] .\n'\
        '   DesMin = Coef[1] * Qmed ** Exp[1] .\n'\
        '   Qmin = MedMin + K(Tr) * DesMin.\n'\
        '   Donde K es un coeficiente que depende de Tr y f(x).\n'\
        '\n'\
        'Parametros\n'\
        '----------\n'\
        'Qmed : Caudal medio calculado con GetQ_Balance.\n'\
        'Coef : Coeficientes para media y desviacion [6.71, 3.29].\n'\
        'Expo : Exponentes para media y desviacion [0.82, 0.64].\n'\
        'Cv : Coeficiente de escalamiento [0.5].\n'\
        'Tr : Periodo de retorno [2.33, 5, 10, 25, 50, 100].\n'\
        'Dist : Distrbucion [gumbel/lognorm].\n'\
        '\n'\
        'Retornos\n'\
        '----------\n'\
        'CellQmax : Caudal maximo para los diferentes periodos.\n'\
        # Calcula la media y desviacion
        MedMin = Coef[0] * Qmed ** Expo[0]
        DesMin = Coef[1] * Qmed ** Expo[1]
        #Itera para todos los periodos de retorno
        Qmin=[]
        for t in Tr:
            #Calcula k
            if Dist is 'gumbel':
                k = (-1*np.sqrt(6)/np.pi)*(0.5772+np.log(-1*np.log(1/float(t))))
            elif Dist is 'lognorm':
                Ztr=norm.ppf(1/float(t))
                k = 1*(np.exp(Ztr*np.sqrt(np.log(1+Cv**2))-0.5*np.log(1-Cv**2))-1)/Cv
            #Calcula el caudal maximo
            Qmin.append(list(MedMin+k*DesMin))
        return np.array(Qmin)

    #------------------------------------------------------
    # Guardado shp de cuencas y redes hidricas
    #------------------------------------------------------
    def Save_Net2Map(self,ruta,dx=cu.dxp,umbral=None,
        qmed=None,Dict=None,DriverFormat='ESRI Shapefile',
        EPSG=4326, NumTramo = True, formato = '%.2f'):
        'Descripcion: Guarda la red hidrica simulada de la cuenca en .shp \n'\
        '   Puede contener un diccionario con propiedades de la red hidrica. \n'\
        '\n'\
        'Parametros\n'\
        '----------\n'\
        'self : no necesita nada es autocontenido.\n'\
        'ruta : Lugar y nombre donde se va a guardar la red hidrica.\n'\
        'dx : Longitud de las celdas planas (Valor de Dx plano asignado a wmf.cu.dxp).\n'\
        'umbral : cantidad de celdas necesarias para corriente (Valor del umbral asignado a self.umbral).\n'\
        'qmed : caudal medio calculado por alguna metodologia.\n'\
        'Dict : Diccionario con parametros de la red hidrica que se quieren imprimir.\n'\
        'DriverFormat : nombre del tipo de archivo vectorial de salida (ver OsGeo).\n'\
        'EPSG : Codigo de proyeccion utilizada para los datos, defecto WGS84.\n'\
        'NumTramo: Poner o no el numero de tramo en cada elemento de la red.\n'\
        '\n'\
        'Retornos\n'\
        '----------\n'\
        'Escribe un archivo vectorial con la estructura de la red hidrica y sus propiedades.\n'\
        #varia el umbral en funcion de self
        if umbral == None:
            umbral = self.umbral
        #division de la cuenca
        acum=cu.basin_acum(self.structure,self.ncells)
        cauce,nod_f,n_nodos=cu.basin_subbasin_nod(self.structure,acum,umbral,self.ncells)
        sub_pert,sub_basin=cu.basin_subbasin_find(self.structure,nod_f,n_nodos,self.ncells)
        sub_basins=cu.basin_subbasin_cut(n_nodos)
        sub_horton,nod_hort=cu.basin_subbasin_horton(sub_basins,self.ncells,n_nodos)
        sub_hort=cu.basin_subbasin_find(self.structure,nod_hort,n_nodos,self.ncells)[0]
        cauceHorton=sub_hort*cauce
        #Obtiene la red en manera vectorial
        nodos = cu.basin_stream_nod(self.structure,acum,umbral,self.ncells)[1]
        netsize = cu.basin_netxy_find(self.structure,nodos,cauceHorton,self.ncells)
        net=cu.basin_netxy_cut(netsize,self.ncells)
        #Para net con caudal medio
        if qmed is not None:
            netsize = cu.basin_netxy_find(self.structure,nodos,cauce*qmed,self.ncells)
            netQmed=cu.basin_netxy_cut(netsize,self.ncells)
        #Para net con tramos
        if NumTramo:
            netsize2 = cu.basin_netxy_find(self.structure,nodos,sub_pert*cauce,self.ncells)
            netTramo = cu.basin_netxy_cut(netsize2,self.ncells)
        #Cortes
        cortes=np.where(net[0,:]==-999)
        cortes=cortes[0].tolist()
        cortes.insert(0,0)
        #Escribe el shp de la red hidrica
        spatialReference = osgeo.osr.SpatialReference()
        spatialReference.ImportFromEPSG(EPSG)
        driver = osgeo.ogr.GetDriverByName(DriverFormat)
        if os.path.exists(ruta):
             driver.DeleteDataSource(ruta)
        shapeData = driver.CreateDataSource(ruta)
        layer = shapeData.CreateLayer('layer1', spatialReference, osgeo.ogr.wkbLineString)
        layerDefinition = layer.GetLayerDefn()
        new_field=osgeo.ogr.FieldDefn('Long[km]',osgeo.ogr.OFTReal)
        layer.CreateField(new_field)
        new_field=osgeo.ogr.FieldDefn('Horton',osgeo.ogr.OFTInteger)
        layer.CreateField(new_field)
        #coloca los nodos
        if NumTramo:
            new_field=osgeo.ogr.FieldDefn('Tramo',osgeo.ogr.OFTInteger)
            layer.CreateField(new_field)
        if qmed is not None:
            new_field=osgeo.ogr.FieldDefn('Qmed[m3s]',osgeo.ogr.OFTReal)
            layer.CreateField(new_field)
        if Dict is not None:
            if type(Dict==dict):
                netDict=[]
                for k in Dict.keys():
                    print (k[:10])
                    print (osgeo.ogr.OFTReal)
                    new_field=osgeo.ogr.FieldDefn(k[:10],osgeo.ogr.OFTReal)
                    layer.CreateField(new_field)
                    netsizeT = cu.basin_netxy_find(self.structure,nodos,cauce*Dict[k],self.ncells)
                    netDict.append(cu.basin_netxy_cut(netsizeT,self.ncells))
        #Para cada tramo
        featureFID=0
        for i,j in zip(cortes[:-1],cortes[1:]):
            line = osgeo.ogr.Geometry(osgeo.ogr.wkbLineString)
            for x,y in zip(net[1,i+1:j],net[2,i+1:j]):
                line.AddPoint_2D(float(x),float(y))
            feature = osgeo.ogr.Feature(layerDefinition)
            feature.SetGeometry(line)
            feature.SetFID(0)
            feature.SetField('Long[km]',(net[1,i+1:j].size*dx)/1000.0)
            feature.SetField('Horton',int(net[0,i+1]))
            if qmed is not None:
                feature.SetField('Qmed[m3s]',float(netQmed[0,i+1]))
            if NumTramo:
                feature.SetField('Tramo',int(netTramo[0,i+1]))
            if Dict is not None:
                if type(Dict==dict):
                    for n,k in zip(netDict,Dict.keys()):
                        feature.SetField(k[:10],float(formato % n[0,i+1]))
            #featureFID+=1
            layer.CreateFeature(feature)
            line.Destroy()
            feature.Destroy()
        shapeData.Destroy()
    def Save_Basin2Map(self,ruta,dx=30.0,Param={},
        DriverFormat='ESRI Shapefile',EPSG=4326, GeoParam = False):
        'Descripcion: Guarda un archivo vectorial de la cuenca en .shp \n'\
        '   Puede contener un diccionario con propiedades. \n'\
        '\n'\
        'Parametros\n'\
        '----------\n'\
        'self : no necesita nada es autocontenido.\n'\
        'ruta : Lugar y nombre donde se va a guardar la cuenca.\n'\
        'dx : Longitud de las celdas planas.\n'\
        'DriverFormat : nombre del tipo de archivo vectorial de salida (ver OsGeo).\n'\
        'EPSG : Codigo de proyeccion utilizada para los datos, defecto WGS84.\n'\
        'GeoParam: (False) determina si calcular de una los parametros geomorfo o no.\n'\
        '\n'\
        'Retornos\n'\
        '----------\n'\
        'Escribe un archivo vectorial de la cuenca.\n'\
        #Obtiene el perimetro de la cuenca
        #nperim = cu.basin_perim_find(self.structure,self.ncells)
        #basinPerim=cu.basin_perim_cut(nperim)

        #Parametros geomorfo
        if GeoParam:
            self.GetGeo_Parameters()
            DictParam = {}
            for k in self.GeoParameters.keys():
                DictParam.update({k[:8]: self.GeoParameters[k]})
        #Genera el shapefile
        spatialReference = osgeo.osr.SpatialReference()
        spatialReference.ImportFromEPSG(EPSG)
        driver = osgeo.ogr.GetDriverByName(DriverFormat)
        if os.path.exists(ruta):
             driver.DeleteDataSource(ruta)
        shapeData = driver.CreateDataSource(ruta)
        layer = shapeData.CreateLayer('layer1', spatialReference, osgeo.ogr.wkbPolygon)
        layerDefinition = layer.GetLayerDefn()
        for p in Param.keys():
            #new_field=osgeo.ogr.FieldDefn(p[:p.index('[')].strip()[:10],osgeo.ogr.OFTReal)
            new_field=osgeo.ogr.FieldDefn(p,osgeo.ogr.OFTReal)
            layer.CreateField(new_field)
        if GeoParam:
            for p in DictParam.keys():
                new_field=osgeo.ogr.FieldDefn(p,osgeo.ogr.OFTReal)
                layer.CreateField(new_field)
        #Calcula el tamano de la muestra
        ring = osgeo.ogr.Geometry(osgeo.ogr.wkbLinearRing)
        for i in self.Polygon.T:
            ring.AddPoint(x=float(i[0]),y=float(i[1]))
        poly=osgeo.ogr.Geometry(osgeo.ogr.wkbPolygon)
        poly.AddGeometry(ring)
        feature = osgeo.ogr.Feature(layerDefinition)
        feature.SetGeometry(poly)
        feature.SetFID(0)
        for p in Param.keys():
            #feature.SetField(p[:p.index('[')].strip()[:10],float("%.2f" % Param[p]))
            feature.SetField(p,float("%.2f" % Param[p]))
        #Si calcula parametros geomorfo
        if GeoParam:
            for p in DictParam.keys():
                feature.SetField(p,float("%.2f" % DictParam[p]))
        layer.CreateFeature(feature)
        poly.Destroy()
        ring.Destroy()
        feature.Destroy()
        shapeData.Destroy()

    #------------------------------------------------------
    # Graficas de la cuenca
    #------------------------------------------------------
    def Plot_basin(self,vec=None,Min=None,
            Max=None,ruta=None,figsize=(10,7),
            ZeroAsNaN ='no',extra_lat=0.0,extra_long=0.0,lines_spaces='Default',
            xy=None,xycolor='b',colorTable=None,alpha=1.0,vmin=None,vmax=None,
            colorbar=True, colorbarLabel = None,axis=None,rutaShp=None,shpWidth = 0.7,
            shpColor = 'r',axloc = 111, fig = None, EPSG = 4326,backMap = False,
            **kwargs):
            #Plotea en la terminal como mapa un vector de la cuenca
            'Funcion: Plot_basin\n'\
            'Descripcion: Genera un plot del mapa entregado.\n'\
            'del mismo en forma de mapa \n'\
            'Parametros Obligatorios:.\n'\
            '   -basin: Vector con la forma de la cuenca.\n'\
            '   -vec: Vector con los valores a plotear.\n'\
            'Parametros Opcionales:.\n'\
            '   -Min: Valor minimo del plot, determina el rango de colores.\n'\
            '   -Max: Valor maximo del plot, determina el rango de colores.\n'\
            '   -ruta: Ruta en la cual se guarda la grafica.\n'\
            '   -colorbar: Muestra o no la barra de colores, defecto: True.\n'\
            '   -figsize: tamano de la ventana donde se muestra la cuenca.\n'\
            '   -ZeroAsNaN: Convierte los valores de cero en NaN.\n'\
            '   -rutaShp: Ruta a un vectorial que se quiera mostrar en el mapa.\n'\
            '   -shpWidth: Ancho de las lineas del shp cargado.\n'\
            '   -shpColor: Color de las lineas del shp cargado.\n'\
            '   -backMap: Pone de fondo un mapa tipo arcGIS.\n'\
            'Otros argumentos:.\n'\
            '   -axis = Entorno de grafica que contiene elementos de las figuras.\n'\
            '   -parallels = Grafica Paralelos, list-like.\n'\
            '   -parallels_labels = Etiquetas de los paralelos, list-like.\n'\
            '       labels = [left,right,top,bottom].\n'\
            '       Ejemplo:\n'\
            '       m.drawparallels(parallels,labels=[False,True,True,False]).\n'\
            '   -parallels_offset: desplazamiento de las coordenadas en X.\n'\
            '   -meridians = Grafica Meridianos, list-like.\n'\
            '   -meridians_labels = Etiquetas de los meridianos, list-like.\n'\
            '   -meridians_offset: desplazamiento de las coordeandas en Y.\n'\
            '   -per_color = Color del perimetro.\n'\
            '   -per_lw = Ancho de linea del perimetro.\n'\
            '   -colorbarLabel = Titulo del colorbar.\n'\
            '   -xy_edgecolor = Color de la linea exterior del scatter .\n'\
            '   -xy_lw = Ancho de linea del Scatter .\n'\
            '   -xy_s = Tamano del scatter.\n'\
            '   -xy_colorbar: Coloca o no barra de colores del scatter enviado (False).\n'\
            '   -show = boolean, si es True muestra la grafica.\n'\
            '   -cbar_ticks: (None) ubicacion de los ticks del cbar.\n'\
            '   -cbar_ticklabels: (None) Labels a poner sobre los ticks.\n'\
            '   -cbar_ticksize: (14) Tamano de los ticks.\n'\
            '   -ShpIsPolygon: Indica si ese shp cargado es poligono (True) o no (False).\n'\
            '   -shpAlpha: transparencia del shp (0.5).\n'\
            'Retorno:.\n'\
            '   -Actualizacion del binario .int\n'\
            '   -m = Para continuar graficando encima del entorno creado en Basemap\n'\
            '       Ejemplo:.\n'\
            '       m = Plot_basin(**args).\n'\
            '       m.scatter(coordenada_x,coordenada_y).\n'
            #Prop de la barra de colores
            cbar_ticklabels = kwargs.get('cbar_ticklabels', None)
            cbar_ticks = kwargs.get('cbar_ticks', None)
            cbar_ticksize = kwargs.get('cbar_ticksize', 14)
            show = kwargs.get('show', True)
            ShpIsPolygon = kwargs.get('ShpIsPolygon',None)
            shpAlpha = kwargs.get('shpAlpha',0.5)
            xy_colorbar = kwargs.get('xy_colorbar', False)
            if lines_spaces == 'Default':
                lines_spaces = cu.dx*cu.ncols*0.05
            #El mapa
            Mcols,Mrows=cu.basin_2map_find(self.structure,self.ncells)
            Map,mxll,myll=cu.basin_2map(self.structure,self.structure[0]
                ,Mcols,Mrows,self.ncells)
            longs=np.array([mxll+0.5*cu.dx+i*cu.dx for i in range(Mcols)])
            lats=np.array([myll+0.5*cu.dy+i*cu.dy for i in range(Mrows)])
            X,Y=np.meshgrid(longs,lats)
            Y=Y[::-1]
            show = kwargs.get('show',True)
            if fig is None:
                fig = pl.figure(figsize = figsize)
            if axis == None:
                ax = fig.add_subplot(axloc)
            else:
                show = False
            m = Basemap(projection='merc',
                llcrnrlat=lats.min()-extra_lat,
                urcrnrlat=lats.max()+extra_lat,
                llcrnrlon=longs.min()-extra_long,
                urcrnrlon=longs.max()+extra_long,
                resolution='c',
                epsg = EPSG)
            parallels = kwargs.get('parallels',np.arange(lats.min(),
                lats.max(),lines_spaces))
            parallels_labels = kwargs.get('parallels_labels',[1,0,0,0])
            parallel_offset = kwargs.get('parallels_offset', 0.001)
            m.drawparallels(parallels,
                labels = parallels_labels,
                fmt="%.2f",
                rotation='vertical',
                xoffset=parallel_offset)
            meridians = kwargs.get('meridians', np.arange(longs.min(),
                longs.max(),lines_spaces))
            meridians_labels = kwargs.get('meridians_labels', [0,0,1,0])
            meridians_offset = kwargs.get('meridians_offset', 0.001)
            m.drawmeridians(meridians,
                labels=meridians_labels,
                fmt="%.2f",
                yoffset=meridians_offset)
            Xm,Ym=m(X,Y)
            #plotea el mapa de fondo de arcGIS
            if backMap:
                m.arcgisimage(server='http://server.arcgisonline.com/ArcGIS', service='World_Topo_Map', xpixels = 1500, verbose = True)
            #Plotea el contorno de la cuenca y la red
            xp,yp = m(self.Polygon[0], self.Polygon[1])
            per_color = kwargs.get('per_color','r')
            per_lw = kwargs.get('per_lw',2)
            m.plot(xp, yp, color=per_color,lw=per_lw)
            #hay una variable la monta
            if vec is not None:
                if vmin is None:
                    vmin = vec.min()
                if vmax is None:
                    vmax = vec.max()
                MapVec,mxll,myll=cu.basin_2map(self.structure,vec,Mcols,Mrows,
                    self.ncells)
                MapVec[MapVec==cu.nodata]=np.nan
                if ZeroAsNaN is 'si':
                    MapVec[MapVec == 0] = np.nan
                if colorTable is not None:
                    cs = m.contourf(Xm, Ym, MapVec.T, 25, alpha=alpha,cmap=colorTable,
                        vmin=vmin,vmax=vmax)
                else:
                    cs = m.contourf(Xm, Ym, MapVec.T, 25, alpha=alpha,
                        vmin=vmin,vmax=vmax)
                cbar_label_size = kwargs.get('cbar_label_size',16)
                if colorbar:
                    cbar = pl.colorbar(cs,orientation='horizontal',pad=0.05)
                    if colorbarLabel is not None:
                        cbar.set_label(colorbarLabel, size = cbar_label_size)
                    if cbar_ticks is not None:
                        cbar.set_ticks(cbar_ticks)
                    if cbar_ticklabels is not None:
                        cbar.ax.set_yticklabels(cbar_ticklabels, size = cbar_ticksize,)
                        cbar.ax.set_xticklabels(cbar_ticklabels, size = cbar_ticksize,)
            #Si hay coordenadas de algo las plotea
            xy_edgecolor = kwargs.get('xy_edgecolor','black')
            xy_lw = kwargs.get('xy_lw',30)
            xy_s = kwargs.get('xy_s',0.5)
            if xy is not None:
                xc,yc=m(xy[0],xy[1])
                sx = m.scatter(xc,yc,c=xycolor,
                    s=xy_s,
                    linewidth=xy_lw,
                    edgecolor=xy_edgecolor)
                if xy_colorbar:
                    pl.colorbar(sx)
            #Si hay una ruta a un shp lo plotea
            if rutaShp is not None:
                if type(rutaShp) == str:
                    m.readshapefile(rutaShp, 'mapashp', linewidth = shpWidth, color = shpColor)
                elif type(rutaShp) == list:
                    for c,shape in enumerate(rutaShp):
                        m.readshapefile(shape, 'mapashp', linewidth = shpWidth[c], color = shpColor[c])
                        if ShpIsPolygon[c]:
                            patches   = []
                            for shape in m.mapashp:
                                patches.append(Polygon(np.array(shape), True))
                            ax.add_collection(PatchCollection(patches, facecolor= shpColor[c],
                                edgecolor=shpColor[c], linewidths=shpWidth[c], zorder=2, alpha = shpAlpha[c]))
            #Guarda
            if ruta is not None:
                pl.savefig(ruta, bbox_inches='tight',pad_inches = 0.25)
            if show is True:
                pl.show()
            if xy is None:
                return m,ax
            else:
                return m, ax, sx
    #Grafica de plot para montar en paginas web o presentaciones
    def Plot_basinClean(self, vec, ruta = None, umbral = 0.0,
        vmin = 0.0, vmax = None, show_cbar = False, **kwargs):
        'Funcion: Plot_basinClean\n'\
        'Descripcion: Genera un plot del mapa entregado en un lienzo limpio.\n'\
        'Parametros Obligatorios:.\n'\
        '   -vec: Vector con los valores a plotear.\n'\
        'Parametros Opcionales:.\n'\
        '   -ruta: ruta donde se guarda el png.\n'\
        '   -umbral: Umbral a partir del cual se plotea variable.\n'\
        '   -vmin: Valor minimo de la variable.\n'\
        '   -vmax: valor maximo de la variable.\n'\
        '   -show_cbar: muestra o no el Cbar del plot.\n'\
        'Otros argumentos:.\n'\
        '   -cmap: Esquema de colores.\n'\
        '   -figsize = Tamano de la figura.\n'\
        '   -cbar_aspect: (20) relacion largo ancho del cbar.\n'\
        '   -cbar_ticks: (None) ubicacion de los ticks del cbar.\n'\
        '   -cbar_ticklabels: (None) Labels a poner sobre los ticks.\n'\
        '   -cbar_ticksize: (14) Tamano de los ticks.\n'\
        '   -show: SE muestra por defecto la figura (True).\n'\
        '   -interpolation: Tipo de interpolacion utilizada por la funcion imshow (ver opciones en matplotlib).\n'\
        'Retorno:.\n'\
        '   -Figura se muestra y se guarda.\n'\
        '   -Coordenadas de los bordes del mapa.\n'\
        #Argumentos kw
        cmap = kwargs.get('cmap','Spectral')
        figsize = kwargs.get('figsize', (10,8))
        cbar_aspect = kwargs.get('cbar_aspect', 20)
        cbar_ticklabels = kwargs.get('cbar_ticklabels', None)
        cbar_ticks = kwargs.get('cbar_ticks', None)
        cbar_ticksize = kwargs.get('cbar_ticksize', 14)
        interpolation = kwargs.get('None')
        show = kwargs.get('show', True)
        #Obtiene la matriz
        M,p = self.Transform_Basin2Map(vec)
        M[(M == -9999) | (M<umbral)] = np.nan
        #Calcula esquinas: izquierda, derecha, abajo, arriba.
        Corners = [p[2],
        p[2]+p[0]*p[4],
        p[3],
        p[3]+p[1]*p[5]]
        #Crea la figura
        fig = pl.figure(figsize=figsize)
        ax = fig.add_subplot(111)
        #plot
        if vmax == None:
            im = pl.imshow(M.T,
                interpolation=interpolation,
                cmap= pl.get_cmap(cmap),
                vmin = vmin)
        else:
            im = pl.imshow(M.T,
                interpolation=interpolation,
                cmap= pl.get_cmap(cmap),
                vmin = vmin,
                vmax = vmax)
        #Quita ejes
        ax.set_xticklabels([])
        ax.set_yticklabels([])
        ax.axis('off')
        #colorca colorbar
        if show_cbar:
            cbar = pl.colorbar(im, aspect = cbar_aspect, )
            if cbar_ticks is not None:
                cbar.set_ticks(cbar_ticks)
            if cbar_ticklabels is not None:
                cbar.ax.set_yticklabels(cbar_ticklabels, size = cbar_ticksize,)
        #Guarda transparente y ajustando bordes
        if ruta is not None:
            pl.savefig(ruta,
                bbox_inches = 'tight',
                pad_inches = 0,
                transparent = True)
        else:
            if show:
                pl.show()
        return Corners, ax
    #Grafica de variables sobre la red
    def Plot_Net(self, vec, vec_c = None,ruta = None,
        q_compare = None, show = True,
        vmin = 0, vmax = 1, umbral = 0.1,**kwargs):
        'Descripcion: Hace un plot de lass variables sobre el cauce \n'\
        '\n'\
        'Parametros\n'\
        '----------\n'\
        'vec : Variable a ser ploteada sobre la red.\n'\
        'vec_c : Variable utilizada para darle color al plot.\n'\
        'ruta : Lugar y nombre donde se va a guardar la cuenca.\n'\
        'q_compare : Caudal o variable de comparacion.\n'\
        'show : Mostrar o no la figura.\n'\
        'vmin : Maximo valor para plotear.\n'\
        'vmax : Minimo valor para plotear.\n'\
        'umbral : Umbral para determinar cauce (constante o la var self.CellCauce).\n'\
        '\n'\
        'kwargs\n'\
        '----------\n'\
        'cmap : mapa de colores "Spectral".\n'\
        'figsize : tamano de la figura (10,8).\n'\
        'escala : Factor para escalar el tamano de los puntos.\n'\
        'show_cbar : Pinta o no barra de colores (True).\n'\
        'clean : Pinta la figura limpia sin ejes (False).\n'\
        'transparent: Guarda la figura sin fondo (False).\n'\
        'grid: Pinta la grilla (True).\n'\
        'size: Tamano de texto en los ejes (16).\n'\
        'ticksize: Tamano de texto en los valores de los ejes (16).\n'\
        'norm: normalizacion de la escala de colores para el cmap.\n'\
        '-cbar_aspect: (20) relacion largo ancho del cbar.\n'\
        '-cbar_ticks: (None) ubicacion de los ticks del cbar.\n'\
        '-cbar_ticklabels: (None) Labels a poner sobre los ticks.\n'\
        '-cbar_ticksize: (14) Tamano de los ticks.\n'\
        '.\n'\
        'Retornos\n'\
        '----------\n'\
        'Hace un plot de una variable sobre la red de drenaje.\n'\
        #Argumentos kw
        cmap = kwargs.get('cmap','Spectral')
        figsize = kwargs.get('figsize', (10,8))
        escala = kwargs.get('escala', 1)
        show_cbar = kwargs.get('show_cbar',True)
        clean = kwargs.get('clean',False)
        transparent = kwargs.get('transparent', False)
        grid = kwargs.get('grid', True)
        size = kwargs.get('size', 14)
        ticksize = kwargs.get('ticksize', 14)
        norm = kwargs.get('norm', None)
        cbar_aspect = kwargs.get('cbar_aspect', 20)
        cbar_ticklabels = kwargs.get('cbar_ticklabels', None)
        cbar_ticks = kwargs.get('cbar_ticks', None)
        cbar_ticksize = kwargs.get('cbar_ticksize', 14)
        #Donde plotea
        if type(umbral) == float or type(umbral) == int :
            pos = np.where(vec>umbral)[0]
        elif type(umbral) == np.ndarray and umbral.shape[0] == self.ncells:
            pos = np.where(umbral == 1)[0]
        x,y = cu.basin_coordxy(self.structure, self.ncells)
        #Vector para pintar si no tiene el vec_c usa vec
        if vec_c is None:
            vec_c = np.copy(vec)
        #Compara o no
        if q_compare is not None:
            vec = vec/q_compare.astype(float)
        #Figura
        fig = pl.figure(figsize=figsize)
        ax = fig.add_subplot(111)
        sca = pl.scatter(x[pos],y[pos],
            s = vec[pos]*escala,
            c = vec_c[pos],
            lw = 0,
            vmin = vmin,
            vmax = vmax,
            cmap = cmap,
            norm = norm)
        ax.patch.set_facecolor('w')
        ax.patch.set_alpha(0.0)
        if grid:
            pl.grid(True)
        #colorca colorbar
        if show_cbar:
            cbar = pl.colorbar(sca, aspect = cbar_aspect, )
            if cbar_ticks  is not  None:
                cbar.set_ticks(cbar_ticks)
            if cbar_ticklabels  is not  None:
                cbar.ax.set_yticklabels(cbar_ticklabels, size = cbar_ticksize,)

        ax.set_xlim(x[pos].min(),x[pos].max())
        ax.set_ylim(y[pos].min(),y[pos].max())
        #Quita ejes
        if clean:
            ax.set_xticklabels([])
            ax.set_yticklabels([])
            ax.axis('off')
        if clean == False:
            ax.set_xlabel('Latitud', size = size)
            ax.set_ylabel('Longitud', size = size)
            ax.tick_params(labelsize = ticksize)
        #Guarda transparente y ajustando bordes
        if ruta is not None:
            pl.savefig(ruta,
                bbox_inches = 'tight',
                pad_inches = 0,
                transparent = transparent,
                edgecolor = 'none',
                facecolor = 'none')
        if show:
            pl.show()
        pl.close(fig)
        #Retorna las 4 coordenadas de las esquinas
        return [x[pos].min(), x[pos].max(), y[pos].min(), y[pos].max()]


    # Grafica barras de tiempos de concentracion
    def Plot_Tc(self,ruta=None,figsize=(8,6),**kwargs):
        keys=self.Tc.keys()
        keys[2]=u'Carr Espana'
        Media=np.array(self.Tc.values()).mean()
        Desv=np.array(self.Tc.values()).std()
        Mediana=np.percentile(self.Tc.values(),50)
        rango=[Media-Desv,Media+Desv]
        color1 = kwargs.get('color1','b')
        color2 = kwargs.get('color2','r')
        colores=[]
        for t in self.Tc.values():
            if t>rango[0] and t<rango[1]:
                colores.append(color1)
            else:
                colores.append(color2)
        axis = kwargs.get('axis',None)
        show = kwargs.get('show',True)
        if axis == None:
            fig=pl.figure(edgecolor='w',facecolor='w',figsize=figsize)
            ax=fig.add_subplot(111)
        else:
            show = False
            ax = axis
        box = ax.get_position()
        ax.set_position([box.x0, box.y0 + box.height * 0.18,
            box.width, box.height * 0.9])
        ax.set_xlim(-0.4,len(keys)+1-0.8)
        ax.bar(range(len(keys)),self.Tc.values(),color=colores)
        Texto='%.2f' % Media
        ax.hlines(Media,-0.4,len(keys)+1-0.8,'k',lw=2,label='$\\mu='+Texto+'$')
        Texto='%.2f' % Mediana
        ax.hlines(Mediana,-0.4,len(keys)+1-0.8,'r',lw=2,label = '$P_{50}='+Texto+'$')
        Texto='%.2f' % (Media+Desv)
        ax.hlines(Media+Desv,-0.4,len(keys)+1-0.8,'b',lw=2,label = u'$\mu+\sigma='+Texto+'$')
        ax.set_xticks(list(np.arange(1,len(keys)+1)-0.8))
        ax.set_xticklabels(keys,rotation=60)
        ylabel = kwargs.get('ylabel',u'Tiempo de concentracion $T_c[hrs]$')
        ax.set_ylabel(ylabel,size=14)
        ax.grid(True)
        ax.legend(loc='upper_right',ncol='3',fontsize='medium')
        if ruta is not None:
            pl.savefig(ruta,bbox_inches='tight')
        if show == True:
            pl.show()
    #Plot de cuace ppal
    def PlotPpalStream(self,ruta = None, figsize = (8,6),axis=None,**kwargs):
        '''GRAFICA EL PERFIL DEL CAUCE PRINCIPAL
        ===================================================================================
        ARGUMENTO  - DEFAULT               - DESCRIPCIoN                  - TIPO
        ruta       - None                  - Ruta para guardar grafica    - str
        figsize    - (8,6)                 - Tamano de la figura          - tuple
        axis       - None                  - Axis o entorno para graficar - matplotlib.axes
        show       - True                  - Muestra grafica              - bool
        fontsize   - 16                    - Tamano de letra              - int,float
        ylabel     - Elevacion $[m.s.n.m]$ - Etiqueta de la ordenada      - str
        cbar_label - None                  - Etiqueta del colorbar        - str
        =================================================================================== '''
        show = kwargs.get('show',True)
        if axis == None:
            fig = pl.figure(figsize = figsize, edgecolor = 'w',facecolor = 'w')
            ax = fig.add_subplot(111)
        else:
            show=False
            ax = axis
        pl.scatter(self.ppal_stream[1]/1000.0,
            self.ppal_stream[0],
            c = self.ppal_slope,
            linewidth = 0)
        ax.grid(True)
        fontsize = kwargs.get('fontsize',16)
        ylabel = kwargs.get('ylabel','Elevacion $[m.s.n.m]$')
        ax.set_xlabel('Distancia $[km]$',size = fontsize)
        ax.set_ylabel(ylabel,size = fontsize)
        ax.tick_params(labelsize = 14)
        cb = pl.colorbar()
        cbar_label = kwargs.get('cbar_label',None)
        if cbar_label  is not  None:
            cb.set_label(cbar_label,fontsize=fontsize)
        if ruta is not None:
            pl.savefig(ruta,bbox_inches='tight')
        if show == True:
            pl.show()
        return ax,cb
    #Plot de histograma de pendientes
    def PlotSlopeHist(self,ruta=None,bins=[0,2,0.2],
        Nsize=1000, figsize = (8,6), fig = None, show = True,**kwargs):
        '''Hace un plot del histograma de distribucion de pendientes en la cuenca.
                Requiere:
                    - ruta: ruta de guaradado de la imagen,
                    - bins: rango inferior, superior y paso para intervalos.
                    - Nsize: Cantidad de datos a usar para realizar el histograma.
                Retorna:
                    - Plot.
                    - Eje del plot si se quiere continuar editando
                **kwargs:
                    - lw: ancho de la linea.
                    - axissize: tamano de los ejes
                    - labelsize: tamano de los nombres'''
        lw = kwargs.get('lw',3)
        labelsize = kwargs.get('labelsize',16)
        axissize = kwargs.get('axissize',15)
        if Nsize>self.ncells:
            Nsize = self.ncells
        pos = np.random.choice(self.ncells,Nsize)
        h,b=np.histogram(self.CellSlope[pos],bins=np.arange(bins[0],bins[1],bins[2]))
        b=(b[:-1]+b[1:])/2.0
        h=h.astype(float)/h.astype(float).sum()
        if fig is None:
            fig=pl.figure(figsize = figsize, edgecolor='w',facecolor='w')
        ax=fig.add_subplot(111)
        ax.plot(b,h,lw=lw)
        ax.grid(True)
        ax.tick_params(labelsize = axissize)
        ax.set_xlabel('Pendiente',size=labelsize)
        ax.set_ylabel('$pdf [\%]$',size=labelsize)
        if ruta is not None:
            pl.savefig(ruta,bbox_inches='tight')
        if show:
            pl.show()
        return ax

    #Plot de histograma de tiempos de viajes en la cuenca
    def Plot_Travell_Hist(self,ruta=None,Nint=10.0):
        #comparacion histogramas de tiempos de respuestas
        bins=np.arange(0,np.ceil(self.CellTravelTime.max()),
            np.ceil(self.CellTravelTime.max())/Nint)
        h_lib,b_lib=np.histogram(self.CellTravelTime,bins=bins)
        h_lib=h_lib.astype(float)/h_lib.sum()
        b_lib=(b_lib[:-1]+b_lib[1:])/2.0
        hc_lib=np.cumsum(h_lib)
        fig=pl.figure(facecolor='w',edgecolor='w')
        ax=fig.add_subplot(111)
        ax.plot(b_lib,h_lib,'b',lw=2,label='Tiempos')
        ax2=ax.twinx()
        ax2.plot(b_lib,hc_lib,'r',lw=2)
        ax2.set_ylim(0,1.1)
        ax.set_xlim(0,np.ceil(self.CellTravelTime.max()))
        ax.grid(True)
        ax.set_xlabel('Tiempo $t [hrs]$',size=14)
        ax.set_ylabel('$pdf[\%]$',size=14)
        ax2.set_ylabel('$cdf[\%]$',size=14)
        ax.set_xticks(b_lib)
        ax.legend(loc=4)
        if ruta is not None:
            pl.savefig(ruta,bbox_inches='tight')
        pl.show()
    #Plot de curva hipsometrica
    def Plot_Hipsometric(self,ruta=None,ventana=10,normed=False,
        figsize = (8,6)):
        #Suaviza la elevacion en el cuace ppal
        elevPpal=pd.Series(self.hipso_ppal[1])
        elevBasin=self.hipso_basin[1]
        if normed==True:
            elevPpal=elevPpal-elevPpal.min()
            elevPpal=(elevPpal/elevPpal.max())*100.0
            elevBasin=elevBasin-elevBasin.min()
            elevBasin=(elevBasin/elevBasin.max())*100.0
        elevPpal=pd.rolling_mean(elevPpal,ventana)
        ppal_acum=(self.hipso_ppal[0]/self.hipso_ppal[0,-1])*100
        basin_acum=(self.hipso_basin[0]/self.hipso_basin[0,0])*100
        #Genera el plot
        fig=pl.figure(edgecolor='w',facecolor='w',figsize = figsize)
        ax=fig.add_subplot(111)
        #box = ax.get_position()
        #ax.set_position([box.x0, box.y0 + box.height * 0.1,
        #   box.width, box.height * 0.9])
        ax.plot(ppal_acum,elevPpal,c='b',lw=3,label='Cacuce Principal')
        ax.plot(basin_acum,elevBasin,c='r',lw=3,label='Cuenca')
        ax.tick_params(labelsize = 14)
        ax.grid()
        ax.set_xlabel('Porcentaje Area Acumulada $[\%]$',size=16)
        if normed==False:
            ax.set_ylabel('Elevacion $[m.s.n.m]$',size=16)
        elif normed==True:
            ax.set_ylabel('Elevacion $[\%]$',size=16)
        lgn1=ax.legend(loc=0)
        if ruta is not None:
            pl.savefig(ruta, bbox_inches='tight')
            pl.close('all')
        else:
            pl.show()
class SimuBasin(Basin):

    def __init__(self,lat=None,lon=None,DEM=None,DIR=None,rute = None, name='NaN',stream=None,
        umbral=500,useCauceMap = None,
        noData=-999,modelType='cells',SimSed=False,SimSlides=False,dt=60,
        SaveStorage='no',SaveSpeed='no',retorno = 0,
        SeparateFluxes = 'no',SeparateRain='no',ShowStorage='no', SimFloods = 'no',
        controlNodos = True, storageConstant = 0.001):
        'Descripcion: Inicia un objeto para simulacion \n'\
        '   el objeto tiene las propieades de una cuenca con. \n'\
        '   la diferencia de que inicia las variables requeridas. \n'\
        '   para simular. \n'\
        '\n'\
        'Parametros\n'\
        '----------\n'\
        'self : Inicia las variables vacias.\n'\
        'lat : Coordenada en X de la salida de la cuenca.\n'\
        'lon : Coordenada en Y de la salida de la cuenca.\n'\
        'name : Nombre con el que se va a conocer la cuenca.\n'\
        '   (defecto = NaN).\n'\
        'stream : Opcional, si se coloca, las coordenadas no tienen.\n'\
        '   que ser exactas, estas se van a corregir para ubicarse.\n'\
        '   en el punto mas cercano dentro de la corriente, este.\n'\
        '   debe ser un objeto del tipo stream.\n'\
        'useCauceMap: Utiliza un mapa de 1 y 0 representando las celdas que.\n'\
        '   son canal, con este mapa corrige el punto de trazado de la cuenca.\n'\
        'umbral : Cantidad minima de celdas para la creacion de cauces.\n'\
        '   (defecto = 500 ).\n'\
        'noData : Valor correspondiente a valores sin dato (defecto = -999).\n'\
        'modelType : Tipo de modelo, por celdas o por laderas (defecto = cells).\n'\
        '   opciones: .\n'\
        '       cells => modela por celdas.\n'\
        '       hills => modela por laderas.\n'\
        'SimSed : Simula si, o no simula sedimentos no.\n'\
        'SimSlides : Simula si, o no simula deslizamientos no.\n'\
        'dt : Tamano del intervlao de tiempo en que trabaj el modelo (defecto=60seg) [secs].\n'\
        'SaveStorage : Guarda o no el almacenamiento.\n'\
        'SaveSpeed : Guarda o no mapas de velocidad.\n'\
        'rute : por defecto es None, si se coloca la ruta, el programa no.\n'\
        '   hace una cuenca, si no que trata de leer una en el lugar especificado.\n'\
        '   Ojo: la cuenca debio ser guardada con la funcion: Save_SimuBasin().\n'\
        'retorno : (defecto = 0), si es cero no se considera alm maximo en .\n'\
        '   el tanque 3, si es 1, si se considera.\n'\
        'SeparateFluxes : Separa el flujo en base, sub-superficial y escorrentia.\n'\
        'SeparateRain : Separa el flujo proveniente de convectivas y de estratiformes.\n'\
        'ShowStorage : Muestra en la salida del modelo el alm promedio en cada uno de los tanques.\n'\
        'controlNodos: Coloca por defecto puntos de control en todos los nodos (True) o no (False).\n'\
        '\n'\
        'Retornos\n'\
        '----------\n'\
        'self : Con las variables iniciadas.\n'\
        #Esta variable es para controlar cuando se reinician las variables de sedimentos 
        self.segunda_cuenca = False
        #Variables de radar
        self.radarDates = []
        self.radarPos = []
        self.radarMeanRain = []
        self.radarCont = 1
        #Si no hay ruta y el global del codigo EPSG existe, traza la cuenca
        if rute is None and int(Global_EPSG) > 0:
            #Si se entrega cauce corrige coordenadas
            if stream is not None:
                error=[]
                for i in stream.structure.T:
                    error.append( np.sqrt((lat-i[0])**2+(lon-i[1])**2) )
                loc=np.argmin(error)
                lat=stream.structure[0,loc]
                lon=stream.structure[1,loc]
            if useCauceMap is not None and useCauceMap.shape == DEM.shape:
                lat,lon = cu.stream_find_to_corr(lat,lon,DEM,DIR,useCauceMap,
                    cu.ncols,cu.nrows)
            #copia la direccion de los mapas de DEM y DIR, para no llamarlos mas
            self.name=name
            self.DEM=DEM
            self.DIR=DIR
            self.modelType=modelType
            self.nodata=noData
            self.umbral = umbral
            self.epsg = Global_EPSG
            #Traza la cuenca
            self.ncells = cu.basin_find(lat,lon,DIR,
                cu.ncols,cu.nrows)
            self.structure = cu.basin_cut(self.ncells)
            #Obtiene las propiedades para el tamano de la cuenca.
            self.DEMvec = self.Transform_Map2Basin(DEM,
                [cu.ncols, cu.nrows, cu.xll, cu.yll, cu.dx, cu.dy])
            self.DIRvec = self.Transform_Map2Basin(DIR,
                [cu.ncols, cu.nrows, cu.xll, cu.yll, cu.dx, cu.dy])
            acum=cu.basin_acum(self.structure,self.ncells)
            cauce,nodos,self.nhills = cu.basin_subbasin_nod(self.structure
                ,acum,umbral,self.ncells)
            self.hills_own,sub_basin = cu.basin_subbasin_find(self.structure,
                nodos,self.nhills,self.ncells)
            self.hills = cu.basin_subbasin_cut(self.nhills)
            models.drena=self.structure
            print(1)
            #Determina la cantidad de celdas para alojar
            if modelType=='cells':
                N=self.ncells
            elif modelType=='hills':
                N=self.nhills
            #aloja variables
            models.v_coef = np.ones((4,N))
            models.h_coef = np.ones((4,N))
            models.v_exp = np.ones((4,N))
            models.h_exp = np.ones((4,N))
            models.max_capilar = np.ones((1,N))
            models.max_gravita = np.ones((1,N))
            models.storage = np.zeros((5,N))
            models.dt = dt
            models.calc_niter = 5
            models.retorno = 0
            models.verbose = 0
            #Define los puntos de control
            models.control = np.zeros((1,N))
            #Si se da la opcion de puntos de control en toda la red lo hace
            if controlNodos:
                if modelType == 'cells':
                    self.GetGeo_Cell_Basics()
                    cauce,nodos,n_nodos = cu.basin_subbasin_nod(
                        self.structure,
                        self.CellAcum,
                        umbral,
                        self.ncells)
                    pos = np.where(nodos!=0)[0]
                    x,y = cu.basin_coordxy(self.structure,self.ncells)
                    idsOrd,xy = self.set_Control(np.vstack([x[pos],y[pos]]),nodos[pos])
                elif modelType == 'hills':
                    models.control = np.ones((1,self.nhills)) * nodos[nodos!=0]
            #Puntos de control de humedad sin control por defecto
            models.control_h = np.zeros((1,N))
            #Define las simulaciones que se van a hacer
            models.sim_sediments=0
            if SimSed is 'si':
                models.sim_sediments=1
            models.sim_slides=0
            if SimSlides:
                models.sim_slides=1
            models.save_storage=0
            if SaveStorage is 'si':
                models.save_storage=1
            models.save_speed=0
            if SaveSpeed is 'si':
                models.save_speed=1
            models.separate_fluxes = 0
            if SeparateFluxes is 'si':
                models.separate_fluxes = 1
            models.separate_rain = 0
            if SeparateRain is 'si':
                models.separate_rain = 1
            models.show_storage = 0
            if ShowStorage is 'si':
                models.show_storage = 1
            if SimFloods == 'si':
                models.sim_floods = 1
            models.storage_constant = storageConstant
            #Determina que la geomorfologia no se ha estimado
            self.isSetGeo = False
        # si hay tura lee todo lo de la cuenca
        elif rute is not None:
            self.__Load_SimuBasin(rute, SimSlides)
        # Obtiene la envolvente de la cuenca
        self.__GetBasinPolygon__()

    def __Load_SimuBasin(self,ruta, sim_slides = False):
        'Descripcion: Lee una cuenca posteriormente guardada\n'\
        '   La cuenca debio ser guardada con SimuBasin.save_SimuBasin\n'\
        '\n'\
        'Parametros\n'\
        '----------\n'\
        'self : Inicia las variables vacias.\n'\
        'ruta : ruta donde se encuentra ubicada la cuenca guardada\n'\
        'Opcionales\n'\
        '----------\n'\
        'sim_slides : (False, True) Carga variables de modelos de deslizamientos previmaente guardadas.\n'\
        'sim_sed : (False, True) Carga variables de modelo de sedimentos previamente guardadas\n'\
        'Retornos\n'\
        '----------\n'\
        'self : La cuenca con sus parametros ya cargada.\n'\
        #Abre el archivo binario de la cuenca
        self.rutaNC = ruta
        gr = netcdf.Dataset(ruta,'a')
        #obtiene las prop de la cuenca
        self.name = gr.nombre
        self.modelType = gr.modelType#.encode()
        self.nodata = gr.noData
        self.umbral = gr.umbral
        self.ncells = gr.ncells
        self.nhills = gr.nhills
        self.epsg = gr.epsg
        models.dt = gr.dt
        models.dxp = gr.dxp
        models.retorno = gr.retorno
        try:
            models.storage_constant = gr.storageConst
        except:
            models.storage_constant = 0.0
        #Si carga deslizamientos
        if sim_slides:
            models.sim_slides = 1
            models.sl_fs = gr.sl_fs
            models.gullienogullie = gr.sl_gullie
            models.sl_gammaw = gr.sl_gammaw
        #Nueva metodologia de geoespacial de la cuenca
        cu.ncols = gr.ncols
        cu.nrows = gr.nrows
        cu.xll = gr.xll
        cu.yll = gr.yll
        cu.dx = gr.dx
        cu.dy = gr.dy
        cu.dxp = gr.dxp
        cu.nodata = gr.noData
        #de acuerdo al tipo de modeloe stablece numero de elem
        if self.modelType[0] is 'c':
            N = self.ncells
        elif self.modelType[0] is 'h':
            N = self.nhills
        #Obtiene las variables base
        GrupoBase = gr.groups['base']
        self.structure = GrupoBase.variables['structure'][:]
        self.hills = GrupoBase.variables['hills'][:]
        self.hills_own = GrupoBase.variables['hills_own'][:]
        #Obtiene las geomorfologicas
        GrupoGeo = gr.groups['Geomorfo']
        self.DEMvec = GrupoGeo.variables['DEM'][:]
        self.DIRvec = GrupoGeo.variables['DIR'][:]
        self.DEM = self.Transform_Basin2Map(self.DEMvec)
        self.DIR = self.Transform_Basin2Map(self.DIRvec)
        #obtiene las propieades del modelo
        GrupoSimHid = gr.groups['SimHidro']
        models.h_coef = np.ones((4,N)) * GrupoSimHid.variables['h_coef'][:]
        models.v_coef = np.ones((4,N)) * GrupoSimHid.variables['v_coef'][:]
        models.h_exp = np.ones((4,N)) * GrupoSimHid.variables['h_exp'][:]
        models.v_exp = np.ones((4,N)) * GrupoSimHid.variables['v_exp'][:]
        models.max_capilar = np.ones((1,N)) * GrupoSimHid.variables['h1_max'][:]
        models.max_gravita = np.ones((1,N)) * GrupoSimHid.variables['h3_max'][:]

        #Propiedades de Sedimentos
        GrupoSimSed = gr.groups['SimSediments']
        models.krus = np.ones((1,N))*GrupoSimSed.variables['Krus'][:]
        models.prus = np.ones((1,N))*GrupoSimSed.variables['Prus'][:]
        models.crus = np.ones((1,N))*GrupoSimSed.variables['Crus'][:]
        models.parliac = np.ones((3,N))*GrupoSimSed.variables['PArLiAc'][:]

        #Variable de drena de acuerdo al tipo de modelo
        if self.modelType[0] is 'c':
            models.drena = np.ones((3,N)) *GrupoSimHid.variables['drena'][:]
        elif self.modelType[0] is 'h':
            models.drena = np.ones((1,N)) * GrupoSimHid.variables['drena'][:]
        models.unit_type = np.ones((1,N)) * GrupoSimHid.variables['unit_type'][:]
        models.hill_long = np.ones((1,N)) * GrupoSimHid.variables['hill_long'][:]
        models.hill_slope = np.ones((1,N)) * GrupoSimHid.variables['hill_slope'][:]
        models.stream_long = np.ones((1,N)) * GrupoSimHid.variables['stream_long'][:]
        models.stream_slope = np.ones((1,N)) * GrupoSimHid.variables['stream_slope'][:]
        models.stream_width = np.ones((1,N)) * GrupoSimHid.variables['stream_width'][:]
        models.elem_area = np.ones((1,N)) * GrupoSimHid.variables['elem_area'][:]
        models.speed_type = np.ones((3)) * GrupoSimHid.variables['speed_type'][:]
        models.storage = np.ones((5,N)) * GrupoSimHid.variables['storage'][:]
        models.control = np.ones((1,N)) * GrupoSimHid.variables['control'][:]
        models.control_h = np.ones((1,N)) * GrupoSimHid.variables['control_h'][:]

        #Propiedades de deslizamientos
        if sim_slides:
            GrupoSlides = gr.groups['SimSlides']
            models.sl_gammas = np.ones((1,N)) * GrupoSlides.variables['gamma_soil'][:]
            models.sl_cohesion = np.ones((1,N)) * GrupoSlides.variables['cohesion'][:]
            models.sl_frictionangle = np.ones((1,N)) * GrupoSlides.variables['friction_angle'][:]
            models.sl_radslope = np.ones((1,N)) * GrupoSlides.variables['rad_slope'][:]
            models.sl_zs = np.ones((1,N)) * GrupoSlides.variables['z_soil'][:]
        #Cierra el archivo
        gr.close()
        #Determina que por defecto debe estar set la geomorfologia
        self.isSetGeo = True

    #------------------------------------------------------
    # Subrutinas de lluvia, interpolacion, lectura, escritura, agrega funcion para variar evp en
    # el tiempo
    #------------------------------------------------------
    def __GetEVP_Serie__(self, index):
        '''Descripcion: Genera una serie que pondera la evp '''
        rng=index
        rad=np.zeros(rng.size)
        for pos,time in enumerate(rng):
            Hora=time
            # Dia del Ano
            dn = Hora.timetuple().tm_yday
            Theta_d = (2 * np.pi * (dn-1))/ 365.
            # (d/d)2
            an = [1.000110, 0.034221, 0.000719]
            bn = [0,        0.001280, 0.000077]
            #
            d   = 0
            tmp = 0
            for i in range(3):
                tmp = (an[i] * np.cos(i*Theta_d)) + (bn[i] * np.sin(i*Theta_d))
                d = d + tmp
            # Delta
            a_n = [0.006918, -0.399912, -0.006758, -0.002697]
            b_n = [0,         0.070257,  0.000907,  0.001480]
            #
            Delta = 0
            tmp   = 0
            for i in range(4):
                tmp = (a_n[i] * np.cos(i*Theta_d)) + (b_n[i] * np.sin(i*Theta_d))
                Delta = Delta + tmp
            #Angulo horario (cada minuto)
            Minutos = (Hora.hour * 60) + Hora.minute
            Horario = 180 - (0.25 * Minutos)
            Horario = (Horario * np.pi)/180.
            # Coseno de Theta
            Latitud = (6.2593 * np.pi)/180.
            Cos_Theta = (np.sin(Latitud)*np.sin(Delta)) + (np.cos(Latitud)*np.cos(Delta)*np.cos(Horario))
            # Radiacion Teorica
            So = 1367 #w/m2
            Q = So * d * Cos_Theta
            # Escala entre 0 y 1
            rad_max=1369.8721876806876
            Q=1*Q/rad_max
            #Guarda
            rad[pos]=Q
        #Se vuelven cero los valores negativos.
        rad[rad<0]=0
        #Serie
        rad=pd.Series(rad,index=rng)
        models.evpserie = np.copy(rad.values)
        return rad

    def rain_interpolate_mit(self,coord,registers,ruta, umbral = 0.01):
        'Descripcion: Interpola la lluvia mediante una malla\n'\
        '   irregular de triangulos, genera campos que son. \n'\
        '   guardados en un binario para luego ser leido por el. \n'\
        '   modelo en el momento de simular. \n'\
        '\n'\
        'Parametros\n'\
        '----------\n'\
        'self : .\n'\
        'coord : Array (2,Ncoord) con las coordenadas de estaciones.\n'\
        'registers : Array (Nest,Nregisters) con los registros de lluvia.\n'\
        'ruta : Ruta con nombre en donde se guardara el binario con.\n'\
        '   la informacion de lluvia.\n'\
        'umbral: Umbral a partir del cual se considera que un campo contiene lluvia\n'\
        '\n'\
        'Retornos\n'\
        '----------\n'\
        'Guarda el binario, no hay retorno\n'\
        'meanRain :  La serie de lluvia promedio interpolada para la cuenca\n'\
        '\n'\
        'Mirar Tambien\n'\
        '----------\n'\
        'rain_interpolate_idw: interpola campos mediante la metodologia idw.\n'\
        'rain_read_bin: Lee binario de registros para ver que tiene.\n'\
        'rain_ncf2bin: Convierte ncf a binario en formato del modelo (para imagenes de radar).\n'\
        #Mira si los registros son un data frame de pandas
        isPandas=False
        if type(registers)==pd.core.frame.DataFrame:
            reg=registers.values.T
            isPandas=True
        else:
            reg=registers
        #inventa estaciones en las esquinas del DEM
        for i,j in zip([0,0,1,1],[0,1,1,0]):
            x=cu.xll+cu.ncols*cu.dx*i
            y=cu.yll+cu.nrows*cu.dx*j
            d=np.sqrt((coord[0,:]-x)**2+(coord[1,:]-y)**2)
            pos=np.argmin(d)
            #Actualiza las coordenadas
            coord=np.vstack((coord.T,np.array([x,y]))).T
            #pone lluvia en ese registro
            reg=np.vstack((reg,reg[pos]))
        #Obtiene las coordenadas de cada celda de la cuenca
        x,y = cu.basin_coordxy(self.structure,self.ncells)
        xy_basin=np.vstack((x,y))
        #Obtiene la malla irregular
        TIN_mesh=Delaunay(coord.T)
        TIN_mesh=TIN_mesh.vertices.T+1
        #Obtiene las pertenencias en la cuenca a la malla
        #TIN_perte = models.rain_pre_mit(xy_basin,TIN_mesh,coord,self.ncells,
        #   TIN_mesh.shape[1],coord.shape[1])
        c = 1
        TIN_perte = np.zeros((1,self.ncells))
        for t in TIN_mesh.T:
            t = t-1
            t = np.append(t, t[0])
            Poly = np.vstack([coord[0][t], coord[1][t]])
            bbPath = mplPath.Path(Poly.T, closed=True)
            Contiene = bbPath.contains_points(xy_basin.T)
            TIN_perte[0][Contiene == True] = c
            c+=1
        #Revisa si todas las celdas quedaron asignadas
        if len(TIN_perte[TIN_perte == 0]) == 0:
            #Selecciona si es por laderas o por celdas
            if self.modelType[0] is 'h':
                maskVector = np.copy(self.hills_own)
            elif self.modelType[0] is 'c':
                maskVector = np.ones(self.ncells)
            #Interpola con tin
            meanRain,posIds = models.rain_mit(xy_basin,
                coord,
                reg,
                TIN_mesh,
                TIN_perte,
                self.nhills,
                ruta,
                umbral,
                maskVector,
                self.ncells,
                coord.shape[1],
                TIN_mesh.shape[1],
                reg.shape[1])
            #Guarda un archivo con informacion de la lluvia
            f=open(ruta[:-3]+'hdr','w')
            f.write('Numero de celdas: %d \n' % self.ncells)
            f.write('Numero de laderas: %d \n' % self.nhills)
            f.write('Numero de registros: %d \n' % reg.shape[1])
            f.write('Numero de campos no cero: %d \n' % posIds.max())
            f.write('Tipo de interpolacion: TIN\n')
            f.write('IDfecha, Record, Lluvia, Fecha \n')
            if isPandas:
                dates=registers.index.to_pydatetime()
                c = 1
                for d,pos,m in zip(dates,posIds,meanRain):
                    f.write('%d, \t %d, \t %.2f, %s \n' % (c,pos,m,d.strftime('%Y-%m-%d-%H:%M')))
                    c+=1
            f.close()
            return meanRain, posIds
        else:
            print('Error: Existen celdas de la cuenca sin asignacion de triangulos, se retornan las coordenadas no asignadas')
            pos = np.where(TIN_perte == 0)[1]
            return xy_basin[0,pos], xy_basin[1,pos]

    def rain_interpolate_idw(self,coord,registers,ruta,p=1,umbral=0.0):
        'Descripcion: Interpola la lluvia mediante la metodologia\n'\
        '   del inverso de la distancia ponderado. \n'\
        '\n'\
        'Parametros\n'\
        '----------\n'\
        'self : .\n'\
        'coord : Array (2,Ncoord) con las coordenadas de estaciones.\n'\
        'registers : DataFrame de pandas (Nest,Nregisters) con los registros de lluvia.\n'\
        'p :  exponente para la interpolacion de lluvia.\n'\
        'ruta : Ruta con nombre en donde se guardara el binario con.\n'\
        '   la informacion de lluvia.\n'\
        'umbral : Umbral de suma total de lluvia bajo el cual se considera\n'\
        '   que un intervalo tiene suficiente agua como para generar reaccion\n'\
        '   (umbral = 0.0) a medida que incremente se generaran archivos mas\n'\
        '   livianos, igualmente existe la posibilidad de borrar informacion.\n'\
        '\n'\
        'Retornos\n'\
        '----------\n'\
        'Guarda el binario, no hay retorno\n'\
        'meanRain :  La serie de lluvia promedio interpolada para la cuenca\n'\
        '\n'\
        'Mirar Tambien\n'\
        '----------\n'\
        'rain_interpolate_mit: interpola campos mediante la metodologia idw.\n'\
        'rain_read_bin: Lee binario de registros para ver que tiene.\n'\
        'rain_ncf2bin: Convierte ncf a binario en formato del modelo (para imagenes de radar).\n'\
        #Mira si los registros son un data frame de pandas
        isPandas=False
        if type(registers)==pd.core.frame.DataFrame:
            reg=registers.values.T
            isPandas=True
        else:
            reg=registers
        #Obtiene las coordenadas de cada celda de la cuenca
        x,y = cu.basin_coordxy(self.structure,self.ncells)
        xy_basin=np.vstack((x,y))
        #Interpola con idw
        if self.modelType[0] is 'h':
            meanRain,posIds = models.rain_idw(xy_basin, coord, reg, p, self.nhills,
                ruta, umbral, self.hills_own, self.ncells, coord.shape[1],reg.shape[1])
        elif self.modelType[0] is 'c':
            meanRain,posIds = models.rain_idw(xy_basin, coord, reg, p, self.nhills,
                ruta, umbral, np.ones(self.ncells), self.ncells, coord.shape[1],reg.shape[1])
        #Guarda un archivo con informacion de la lluvia
        f=open(ruta[:-3]+'hdr','w')
        f.write('Numero de celdas: %d \n' % self.ncells)
        f.write('Numero de laderas: %d \n' % self.nhills)
        f.write('Numero de registros: %d \n' % reg.shape[1])
        f.write('Numero de campos no cero: %d \n' % posIds.max())
        f.write('Tipo de interpolacion: IDW, p= %.2f \n' % p)
        f.write('IDfecha, Record, Lluvia, Fecha \n')
        if isPandas:
            dates=registers.index.to_pydatetime()
            c = 1
            for d,pos,m in zip(dates,posIds,meanRain):
                f.write('%d, \t %d, \t %.2f, %s \n' % (c,pos,m,d.strftime('%Y-%m-%d-%H:%M')))
                c+=1
        f.close()
        return meanRain,posIds

    def rain_radar2basin_from_asc(self,ruta_in,ruta_out,fechaI,fechaF,dt,
        pre_string,post_string,fmt = '%Y%m%d%H%M',conv_factor=1.0/12.0,
        umbral = 0.0):
        'Descripcion: Genera campos de lluvia a partir de archivos asc. \n'\
        '\n'\
        'Parametros\n'\
        '----------\n'\
        'self : .\n'\
        'ruta_in: Ruta donde se encuentran loas .asc.\n'\
        'ruta_out: Ruta donde escribe el binario con la lluvia.\n'\
        'fechaI: Fecha de inicio de registros.\n'\
        'fechaF: Fecha de finalizacion de registros.\n'\
        'dt: Intervalo de tiempo entre registros.\n'\
        'Retornos\n'\
        '----------\n'\
        'Guarda el binario, no hay retorno\n'\
        'meanRain :  La serie de lluvia promedio.\n'\
        'meanRain :  La serie de lluvia promedio.\n'\
        '\n'\
        'Mirar Tambien\n'\
        '----------\n'\
        'rain_interpolate_idw: interpola campos mediante la metodologia idw.\n'\
        'rain_radar2basin_from_array: Mete campos de lluvia mediante multiples arrays.\n'\
        #Edita la ruta de salida
        if ruta_out.endswith('.hdr') or ruta_out.endswith('.bin'):
            ruta_bin = ruta_out[:-3]+'.bin'
            ruta_hdr = ruta_out[:-3]+'.hdr'
        else:
            ruta_bin = ruta_out+'.bin'
            ruta_hdr = ruta_out+'.hdr'
        #Establece la cantidad de elementos de acuerdo al tipo de cuenca
        if self.modelType[0] is 'c':
            N = self.ncells
        elif self.modelType[0] is 'h':
            N = self.nhills
        #Guarda la primera entrada como un mapa de ceros
        models.write_int_basin(ruta_bin,np.zeros((1,N)),1,N,1)
        #Genera la lista de las fechas.
        ListDates,dates = __ListaRadarNames__(ruta_in,
            fechaI,fechaF,
            fmt,post_string,pre_string,dt)
        #Lee los mapas y los transforma
        cont = 1
        meanRain = []
        posIds = []
        for l in ListDates:
            Map,p = read_map_raster(ruta_in + l)
            vec = self.Transform_Map2Basin(Map,p) * conv_factor
            #Si el mapa tiene mas agua de un umbral
            if vec.sum() > umbral:
                #Actualiza contador, lluvia media y pocisiones
                cont +=1
                meanRain.append(vec.mean())
                posIds.append(cont)
                #Guarda el vector
                vec = vec*1000; vec = vec.astype(int)
                models.write_int_basin(ruta_bin,np.zeros((1,N))+vec,cont,N,1)
            else:
                #lluvia media y pocisiones
                meanRain.append(0.0)
                posIds.append(1)
        posIds = np.array(posIds)
        meanRain = np.array(meanRain)
        #Guarda un archivo con informacion de la lluvia
        f=open(ruta_hdr[:-3]+'hdr','w')
        f.write('Numero de celdas: %d \n' % self.ncells)
        f.write('Numero de laderas: %d \n' % self.nhills)
        f.write('Numero de registros: %d \n' % meanRain.shape[0])
        f.write('Numero de campos no cero: %d \n' % posIds.max())
        f.write('Tipo de interpolacion: radar \n')
        f.write('IDfecha, Record, Lluvia, Fecha \n')
        c = 1
        for d,pos,m in zip(dates,posIds,meanRain):
            f.write('%d, \t %d, \t %.2f, %s \n' % (c,pos,m,d.strftime('%Y-%m-%d-%H:%M')))
            c+=1
        f.close()
        return np.array(meanRain),np.array(posIds)

    def rain_radar2basin_from_array(self,vec=None,ruta_out=None,fecha=None,dt=None,
        status='update',umbral = 0.01, doit = False):
        'Descripcion: Genera campos de lluvia a partir de archivos array\n'\
        '\n'\
        'Parametros\n'\
        '----------\n'\
        'self : .\n'\
        'vec: Array en forma de la cuenca con la informacion.\n'\
        'ruta_out: Ruta donde escribe el binario con la lluvia.\n'\
        'fecha: Fecha del registro actual.\n'\
        'dt: Intervalo de tiempo entre registros.\n'\
        'status: Estado en el cual se encuentra el binario que se va a pasar a campo.\n'\
        '   update: (Defecto) Con este estado se genera un binario y se agregan nuevos campos.\n'\
        '   old: Estado para abrir y tomar las propiedades de self.radar.. para la generacion de un binario.\n'\
        '   close: Cierra un binario que se ha generado mediante update.\n'\
        '   reset: Reinicia las condiciones de self.radar... para la creacion de un campo nuevo.\n'\
        'doit: Independiente del umbral escribe el binario en la siguiente entrada.\n'\
        'Retornos\n'\
        '----------\n'\
        'Guarda el binario, no hay retorno\n'\
        'meanRain :  La serie de lluvia promedio.\n'\
        '\n'\
        'Mirar Tambien\n'\
        '----------\n'\
        'rain_interpolate_idw: interpola campos mediante la metodologia idw.\n'\
        'rain_radar2basin_from_asc: Mete campos de lluvia mediante multiples arrays.\n'\
        #Edita la ruta de salida
        if ruta_out  is not  None:
            if ruta_out.endswith('.hdr') or ruta_out.endswith('.bin'):
                ruta_bin = ruta_out[:-4]+'.bin'
                ruta_hdr = ruta_out[:-4]+'.hdr'
            else:
                ruta_bin = ruta_out+'.bin'
                ruta_hdr = ruta_out+'.hdr'
        #Establece la cantidad de elementos de acuerdo al tipo de cuenca
        if self.modelType[0] is 'c':
            N = self.ncells
        elif self.modelType[0] is 'h':
            N = self.nhills
            try:
                if vec.shape[0]  == self.ncells:
                    vec = self.Transform_Basin2Hills(vec,sumORmean=1)
            except:
                pass
        # De acerudo al estado actualiza las variables o guarda el
        # binario final
        actualizo = 1
        if status == 'update':
            #Entrada 1 es la entrada de campos sin lluvia
            if len(self.radarDates) == 0:
                models.write_int_basin(ruta_bin,np.zeros((1,N)),1,N,1)
            if vec.mean() > umbral or doit:
                #Actualiza contador, lluvia media y pocisiones
                self.radarCont +=1
                self.radarMeanRain.append(vec.mean())
                self.radarPos.append(self.radarCont)
                #Guarda el vector
                vec = vec*1000; vec = vec.astype(int)
                models.write_int_basin(ruta_bin,np.zeros((1,N))+vec,
                    self.radarCont,N,1)
                actualizo = 0
            else:
                #lluvia media y pocisiones
                self.radarMeanRain.append(0.0)
                self.radarPos.append(1)
            self.radarDates.append(fecha)
        #Si ya no va a agregar nada, no agrega mas campos y genera el .hdr
        elif status == 'close':
            self.radarMeanRain = np.array(self.radarMeanRain)
            self.radarPos = np.array(self.radarPos)
            #Guarda un archivo con informacion de la lluvia
            f=open(ruta_hdr[:-3]+'hdr','w')
            f.write('Numero de celdas: %d \n' % self.ncells)
            f.write('Numero de laderas: %d \n' % self.nhills)
            f.write('Numero de registros: %d \n' % self.radarMeanRain.shape[0])
            f.write('Numero de campos no cero: %d \n' % self.radarPos.max())
            f.write('Tipo de interpolacion: radar \n')
            f.write('IDfecha, Record, Lluvia, Fecha \n')
            c = 1
            for d,pos,m in zip(self.radarDates,
                self.radarPos,self.radarMeanRain):
                f.write('%d, \t %d, \t %.2f, %s \n' % (c,pos,m,d.strftime('%Y-%m-%d-%H:%M')))
                c+=1
            f.close()
            #Vuelve las variables listas de nuevo
            self.radarMeanRain = self.radarMeanRain.tolist()
            self.radarPos = self.radarPos.tolist()
        elif status == 'reset':
            #Variables de radar
            self.radarDates = []
            self.radarPos = []
            self.radarMeanRain = []
            self.radarCont = 1
        elif status == 'old':
            #si es un archivo viejo, lo abre para tomar las variables y continuar en ese punto
            f=open(ruta_hdr[:-3]+'hdr','r')
            Lista = f.readlines()
            self.radarCont = int(Lista[3].split()[-1])
            cantidadIds = int(Lista[2].split()[-1])
            f.close()
            #Abre con numpy para simplificar las cosas
            a = np.loadtxt(ruta_hdr,skiprows=6,dtype='str').T
            if self.radarCont >= 1 and cantidadIds > 1:
                self.radarPos = [int(i.split(',')[0]) for i in a[1]]
                self.radarMeanRain = [float(i.split(',')[0]) for i in a[2]]
                for i in a[3]:
                    d = datetime.datetime.strptime(i,'%Y-%m-%d-%H:%M')
                    self.radarDates.append(d)
            else:
                self.radarPos = [int(a[1].split(',')[0])]
                self.radarMeanRain = [float(a[2].split(',')[0])]
                self.radarDates = [datetime.datetime.strptime(a[-1], '%Y-%m-%d-%H:%M')]
        return actualizo

    #------------------------------------------------------
    # Subrutinas para preparar modelo
    #------------------------------------------------------
    def set_Geomorphology(self,umbrales=[30,500],stream_width=None):
        'Descripcion: calcula las propiedades geomorfologicas necesarias \n'\
        '   para la simulacion. \n'\
        '\n'\
        'Parametros\n'\
        '----------\n'\
        'self : Inicia las variables vacias.\n'\
        'umbrales : Lista con la cantidad de celdas necesarias .\n'\
        '   para que una celda sea: ladera, carcava o cauce .\n'\
        'stream_width = Ancho del canal en cada tramo (opcional).\n'\
        '\n'\
        'Retornos\n'\
        '----------\n'\
        'self : Con las variables geomorfologicas de simulacion iniciadas.\n'\
        '   models.drena : Numero de celda o ladera destino. \n'\
        '   models.nceldas : Numero de celdas o laderas. \n'\
        '   models.unit_type : tipo de celda, en el caso de ladera no aplica.\n'\
        '       1: Celda tipo ladera.\n'\
        '       2: Celda tipo transitorio.\n'\
        '       3: Celda tipo cauce.\n'\
        '   models.hill_long : Longitud promedio de la ladera (o celda). \n'\
        '       - Celdas: Longitud de cada elemento (cu.dxp) para ortogonal y 1.43*cu.dxp para diagonal. \n'\
        '       - Laderas: Promedio de las longitudes recorridas entre una celda y la corriente de la ladera. \n'\
        '   models.hill_slope : Pendiente de cada ladera (promedio) o celda.\n'\
        '   models.stream_long : Longitud de cada tramo de cuace. \n'\
        '       - Celdas: Es cu.dxp: ortogonal, 1.42*cu.dxp: Diagonal. \n'\
        '       - Laderas: Es la Longitud de los elementos del cauce que componen la ladera. \n'\
        '   models.stream_slope : Pendiente de cada tramo de cauce. \n'\
        '       - Celdas: Es la pendiente estimada por self.GetGeo_Cell_Basics(). \n'\
        '       - Laderas: Pendiente promedio de las laderas. \n'\
        '   models.stream_width : Ancho de cada tramo de cauce. \n'\
        '       - Celdas: No aplica y se hacen todos iguales a la unidad. \n'\
        '       - Laderas: Es el ancho del canal calculado a partir de geomrofologia o entregado. \n'\
        '   models.elem_area : Area de cada celda (cu.dxp**2) o ladera (nceldasLadera * cu.dxp**2). \n'\
        #Obtiene lo basico para luego pasar argumentos
        acum,hill_long,pend,elev = cu.basin_basics(self.structure,
            self.DEMvec,self.DIRvec,self.ncells)
        #Obtiene la pendiente y la longitud de las corrientes
        cauce,nodos,trazado,n_nodos,n_cauce = cu.basin_stream_nod(
            self.structure,acum,umbrales[1],self.ncells)
        stream_s,stream_l = cu.basin_stream_slope(
            self.structure,elev,hill_long,nodos,n_cauce,self.ncells)
        stream_s[np.isnan(stream_s)]=self.nodata
        #Obtiene para metros por subn cuencas
        sub_horton,nod_horton = cu.basin_subbasin_horton(self.hills,self.ncells,
            self.hills.shape[1])
        sub_basin_long,max_long,nodo_max_long = cu.basin_subbasin_long(
            self.hills_own,cauce,hill_long,self.hills,
            sub_horton,self.hills.shape[1],self.ncells)
        #Obtiene las propiedades por laderas de los cauces
        stream_slope,stream_long = cu.basin_subbasin_stream_prop(
            self.hills_own,cauce,hill_long,
            pend,self.hills.shape[1],self.ncells)
        #opbtiene el ancho si noe s dado lo asume igual a uno
        if stream_width is None:
            stream_width=np.ones(self.ncells)
        #De acuerdo a si el modelo es por laderas o por celdas agrega lass varaibeles
        if self.modelType[0]=='c':
            #Obtiene el tipo de celdas
            unit_type = cu.basin_stream_type(self.structure,
                acum,umbrales,len(umbrales),self.ncells)
            #Asigna variables
            models.drena = np.ones((1,self.ncells))*self.structure
            models.nceldas = self.ncells
            models.unit_type = np.ones((1,self.ncells))*unit_type
            models.hill_long = np.ones((1,self.ncells))*hill_long
            models.hill_slope = np.ones((1,self.ncells))*pend
            models.stream_long = np.ones((1,self.ncells))*np.percentile(hill_long, 40)
            models.stream_slope = np.ones((1,self.ncells))*pend
            models.stream_width = np.ones((1,self.ncells))*stream_width
            models.elem_area = np.ones((1,self.ncells))*cu.dxp**2.0
        elif self.modelType[0]=='h':
            N=self.hills.shape[1]
            models.drena = np.ones((1,N))*self.hills[1]
            models.nceldas = self.hills.shape[1]
            models.unit_type = np.ones((1,N))*np.ones(N)*3
            models.hill_long = np.ones((1,N))*sub_basin_long
            models.hill_slope = np.ones((1,N))*self.Transform_Basin2Hills(pend)
            models.stream_long = np.ones((1,N))*stream_long
            models.stream_slope = np.ones((1,N))*stream_slope
            models.stream_width = np.ones((1,N))*cu.basin_subbasin_map2subbasin(self.hills_own,stream_width,self.nhills,cauce,0,self.ncells)
            no0min = models.stream_width[models.stream_width!=0].min()
            models.stream_width[models.stream_width==0] = no0min
            models.elem_area = np.ones((1,N))*np.array([self.hills_own[self.hills_own==i].shape[0] for i in range(1,self.hills.shape[1]+1)])*cu.dxp**2.0
        #Ajusta variable de que la geomorfologia esta calculada
        self.isSetGeo = True

    def set_Speed_type(self,types=np.ones(3)):
        'Descripcion: Especifica el tipo de velocidad a usar en cada \n'\
        '   nivel del modelo. \n'\
        '\n'\
        'Parametros\n'\
        '----------\n'\
        'self : Inicia las variables vacias.\n'\
        'types : tipos de velocidad .\n'\
        '   1. Velocidad tipo embalse lineal, no se especifica h_exp.\n'\
        '   2. Velocidad onda cinematica, se debe especificar h_exp.\n'\
        '\n'\
        'Retornos\n'\
        '----------\n'\
        'self : Con la variable models.speed_type especificada.\n'\
        #Especifica la ecuacion de velocidad a usar en cada nivel del modelo
        for c,i in enumerate(types):
            if i==1 or i==2:
                models.speed_type[c]=i
            else:
                models.speed_type[c]=1

    def set_Floods(self,var,VarName, umbral = 1000, NumCeldas = 6, Default = False):
        'Descripcion: Aloja las variables del sub modelo de inundaciones\n'\
        '\n'\
        'Parametros\n'\
        'Mierda PUTA\n'\
        '----------\n'\
        'var : Variable que describe la propiedad (constante, ruta, mapa o vector).\n'\
        'varName: Nombre de la variable a ingresar en el modelo.\n'\
        '   GammaWater : Densidad del agua (Defecto: 1000).\n'\
        '   GammaSoil : Densidad del sedimento (Defecto: 2600).\n'\
        '   VelArea: Factor conversion Velocidad - Area (Defecto: 1/200).\n'\
        '   Cmax: Maxima concentracion de sedimentos (Defecto 0.75).\n'\
        '   MaxIter: Cantidad maxima de iteraciones para obtener la mancha de inunudacion.\n'\
        '   VelUmbral: Velocidad minima para que se de el flujo de escombros ( Defecto: 3 m/s).\n'\
        '   Stream_W: Ancho del canal en cada celda (ncells).\n'\
        '   Stream_D50: Tamano de particula percentil 50 (ncells).\n'\
        '   HAND: Modelo de elevacion relativa de celda ladera a celda cauce (ncells), no se ponde nada.\n'\
        '   umbral: Umbral para determinar corriente segun HAND debe coincidir con el umbral del modelo.\n'\
        '   Slope: Pendiente del canal, (no se ponde variable)\n'\
        '   Default: Poner o no parametros por defecto (False)\n'\
        'Retornos\n'\
        '----------\n'\
        'self : variables iniciadas en el modelo bajo los nombres de:.\n'\
        '   wmf.models.flood_dw.\n'\
        '   wmf.models.flood_dsed.\n'\
        '   wmf.models.flood_av.\n'\
        '   wmf.models.flood_cmax.\n'\
        '   wmf.models.flood_d50.\n'\
        '   wmf.models.flood_hand.\n'\
        '   wmf.models.flood_aquien.\n'\
        #Si el modelo es tipo ladera agrega la variable
        if self.modelType[0] == 'h':
            return 'El modelo por laderas no simula inundaciones.'
        #Pone el gamma del agua por defecto
        if Default:
            models.flood_dw = 1000
            models.flood_dsed = 2600
            models.flood_av = 1./200.0
            models.flood_cmax = 0.75
            models.flood_umbral = 3.0
            models.flood_max_iter = 10
            models.flood_step = 1.0
        #Obtiene el vector que va a alojar en el modelo
        if VarName != 'GammaWater' and VarName != 'GammaSoil' and VarName != 'VelArea' and VarName != 'Cmax' and VarName != 'VelUmbral':
            isVec=False
            if type(var) is str:
                #Si es un string lee el mapa alojado en esa ruta
                Map,Pp = read_map_raster(var)
                Vec = self.Transform_Map2Basin(Map,Pp)
                isVec=True
            elif type(var) is int or float:
                Vec = np.ones((1,self.ncells))*var
                isVec=True
            elif type(var) is np.ndarray and var.shape[0] == self.ncells:
                Vec = var
                isVec=True
            #finalmente mete la variable en el modelo
            N = self.ncells
            if VarName is 'Stream_W' :
                models.flood_w = np.ones((1,N))*Vec
            elif VarName is 'Stream_D50':
                models.flood_d50 = np.ones((1,N))*Vec
            elif VarName is 'HAND':
                self.GetGeo_HAND(umbral = umbral)
                models.flood_hand = np.ones((1,N))*np.copy(self.CellHAND)
                models.flood_aquien = np.ones((1,N))*np.copy(self.CellHAND_drainCell)
            elif VarName is 'Slope':
                self.GetGeo_Cell_Basics()
                models.flood_slope = np.ones((1,N))*np.sin(np.arctan(self.CellSlope))
            elif VarName is 'Sections':
                self.GetGeo_Sections(NumCeldas = NumCeldas)
                models.flood_sections = np.ones((NumCeldas*2+1,N)) * self.Sections
                models.flood_sec_cells = np.ones((NumCeldas*2+1,N)) * self.Sections_Cells
        elif VarName == 'GammaWater':
            models.flood_dw = var
        elif VarName == 'GammaSoil':
            models.flood_dsed = var
        elif VarName == 'VelArea':
            models.flood_av = var
        elif VarName == 'Cmax':
            models.flood_cmax = var
        elif VarName == 'VelUmbral':
            models.flood_umbral = var
        elif VarName == 'MaxIter':
            models.flood_max_iter = var

    def set_PhysicVariables(self,modelVarName,var,pos,mask=None):
        'Descripcion: Coloca las variables fisicas en el modelo \n'\
        '   Se debe assignarel nombre del tipo de variable, la variable\n'\
        '   y la posicion en que esta va a ser insertada\n'\
        '\n'\
        'Parametros\n'\
        '----------\n'\
        'self : Objeto de la cuenca 2que luego se va a simular.\n'\
        'modelVarName : Nombre de la variable del modelo que se va a incertar.\n'\
        '   - h_coef: coeficientes horizontales\n'\
        '       [0]: Flujo de Escorrentia.\n'\
        '       [1]: Flujo Sub-superficial.\n'\
        '       [2]: Flujo subterraneo.\n'\
        '       [3]: Flujo en cauces.\n'\
        '   - h_exp: exponentes del tipo v = h_coef*(A**h_exp)\n'\
        '       [0]: Escorrentia.\n'\
        '       [1]: sub-superficial.\n'\
        '       [2]: subterraneo.\n'\
        '       [3]: cauce.\n'\
        '   - v_coef.\n'\
        '       [0]: Tasa evaporacion.\n'\
        '       [1]: Infiltracion.\n'\
        '       [2]: Percolacion.\n'\
        '       [3]: Perdidas (0).\n'\
        '   - v_exp: procesos verticales no lineales\n'\
        '       (no implementado dentro del modelo)\n'
        '   - capilar.\n'\
        '   - gravit.\n'\
        'var : variable que ingresa en el modelo, esta puede ser:.\n'\
        '   - Ruta: una ruta del tipo string.\n'\
        '   - Escalar : Un valor escalar que se asignara a toda la cuenca.\n'\
        '   - Vector : Un vector con la informacion leida (1,ncells).\n'\
        'pos : Posicion de insercion, aplica para : h_coef, v_coef,.\n'\
        '   h_exp, v_exp.\n'\
        '\n'\
        'Retornos\n'\
        '----------\n'\
        'Guarda el binario, no hay retorno\n'\
        '\n'\
        'Mirar Tambien\n'\
        '----------\n'\
        'rain_interpolate_idw: interpola campos mediante la metodologia idw.\n'\
        'rain_read_bin: Lee binario de registros para ver que tiene.\n'\
        'rain_ncf2bin: Convierte ncf a binario en formato del modelo (para imagenes de radar).\n'\

        #Obtiene el vector que va a alojar en el modelo
        isVec=False
        if type(var) is str:
            #Si es un string lee el mapa alojado en esa ruta
            Map,Pp = read_map_raster(var)
            Vec = self.Transform_Map2Basin(Map,Pp)
            isVec=True
        elif type(var) is int or float:
            Vec = np.ones((1,self.ncells))*var
            isVec=True
        elif type(var) is np.ndarray and var.shape[0] == self.ncells:
            Vec = var
            isVec=True
        #Si el modelo es tipo ladera agrega la variable
        if self.modelType[0] is 'h':
            Vec = self.Transform_Basin2Hills(Vec,mask=mask)
        #finalmente mete la variable en el modelo
        if modelVarName is 'h_coef':
            models.h_coef[pos] = Vec
        elif modelVarName is 'h_exp':
            models.h_exp[pos] = Vec
        elif modelVarName is 'v_coef':
            models.v_coef[pos] = Vec
        elif modelVarName is 'v_exp':
            models.v_exp[pos] = Vec
        elif modelVarName is 'capilar':
            models.max_capilar[0] = Vec
        elif modelVarName is 'gravit':
            models.max_gravita[0] = Vec

    def set_Storage(self,var,pos,hour_scale=False):
        'Descripcion: \n'\
        '   Establece el almacenamiento inicial del modelo\n'\
        '   la variable puede ser un valor, una ruta o un vector.\n'\
        '\n'\
        'Parametros\n'\
        '----------\n'\
        'var : Variable con la cual se va a iniciar el almancenamiento.\n'\
        '   - ruta : es una ruta a un archivo binario de almacenamiento.\n'\
        '   - escalar : valor de almacenamiento constate para toda la cuenca.\n'\
        '   - vector : Vector con valores para cadda unidad de la cuenca.\n'\
        'pos : Posicion de insercion,.\n'\
        '   - 0 :  alm cpailar.\n'\
        '   - 1 :  alm superficial.\n'\
        '   - 2 :  alm sub-superficial.\n'\
        '   - 3 :  alm subterraneo.\n'\
        '   - 4 :  alm cauce.\n'\
        '   - fecha: en caso de que var sea la ruta a StOhdr:\n'\
        '       - valor: puede ser un entero con la posicion.\n'\
        '       - str fecha: puede ser un srint con la fecha: YYYY-MM-DD-HH:MM :\n'\
        'hour_scale: Buscar fechas a escala horaria o minutos?.\n'\
        '\n'\
        'Retornos\n'\
        '----------\n'\
        'Guarda el binario, no hay retorno\n'\
        '\n'\
        'Mirar Tambien\n'\
        '----------\n'\
        'save_storage(slef,storage).\n'\
        #Determina el tipo de unidades del modelo
        if self.modelType[0] is 'c':
            N = self.ncells
        elif self.modelType[0] is 'h':
            N = self.nhills
        #Obtiene el vector que va a alojar en el modelo
        isVec=False
        if type(var) is str:
            var_bin,var_hdr = __Add_hdr_bin_2route__(var,storage = True)
            if type(pos) is not str:
                #Si es un string lee el binario de almacenamiento alojado en esa ruta
                Vec,res = models.read_float_basin_ncol(var_bin,pos+1,N,5)
            if type(pos) == str:
                # Lee las fechas
                L = np.loadtxt(var_hdr,skiprows = 5, usecols = (0,6), dtype = str)
                if hour_scale == False:
                    Fechas = [datetime.datetime.strptime(i[1], '%Y-%m-%d-%H:%M') for i in L]
                else:
                    Fechas = [datetime.datetime.strptime(i[1][:-3], '%Y-%m-%d-%H') for i in L]
                try:
                    if hour_scale == False:
                        posFecha = Fechas.index(datetime.datetime.strptime(pos, '%Y-%m-%d-%H:%M'))
                    else:
                        posFecha = Fechas.index(datetime.datetime.strptime(pos[:-3], '%Y-%m-%d-%H'))
                except:
                    print('Error: no se encuentra la fecha especificada en el archivo '+var)
                Vec,res = models.read_float_basin_ncol(var_bin,posFecha+1,N,5)
            isVec=True
            for p in range(5):
                models.storage[p] = Vec[p]
        elif type(var) is int or float:
            Vec = np.ones((1,N))*var
            isVec=True
            models.storage[pos] = Vec
        elif type(var) is np.ndarray and var.shape[0] == N:
            Vec = var
            isVec=True
            models.storage[pos] = Vec

    def set_Control(self,coordXY,ids,tipo = 'Q'):
        'Descripcion: \n'\
        '   Establece los puntos deonde se va a realizar control del caudal\n'\
        '   y de la humedad simulada.\n'\
        '\n'\
        'Parametros\n'\
        '----------\n'\
        'coordXY : Coordenadas de los puntos de control [2,Ncontrol].\n'\
        'ids : Identificacion de los puntos de control.\n'\
        'tipo: Control para caudal [Q] o para humedad [H].\n'\
        '\n'\
        'Retornos\n'\
        '----------\n'\
        'Define los puntos de control en las variables:\n'\
        '   - models.control : control de caudal y sedimentos.\n'\
        '   - models.control_h : control de la humedad del suelo.\n'\
        'IdsControl : El orden en que quedaron los puntos de control al interior.\n'\
        '\n'\
        'Mirar Tambien\n'\
        '----------\n'\
        'set_record, set_storage.\n'\
        #Obtiene los puntos donde hay coordenadas
        if tipo is 'Q':
            xyNew, basinPts, order = self.Points_Points2Stream(coordXY,ids)
            if self.modelType[0] is 'c':
                models.control[0] = basinPts
                IdsConvert = basinPts[basinPts!=0]
            elif self.modelType[0] is 'h':
                unitario = basinPts / basinPts
                pos = self.hills_own * self.CellCauce * unitario
                posGrande = self.hills_own * self.CellCauce * basinPts
                IdsConvert = posGrande[posGrande!=0] / pos[pos!=0]
                models.control[0][pos[pos!=0].astype(int).tolist()] = IdsConvert
        elif tipo is 'H':
            xyNew = coordXY
            basinPts, order = self.Points_Points2Basin(coordXY,ids)
            if self.modelType[0] == 'c':
                models.control_h[0] = basinPts
                IdsConvert = basinPts[basinPts!=0]
            elif self.modelType[0] == 'h':
                unitario = basinPts / basinPts
                pos = self.hills_own * unitario
                posGrande = self.hills_own * basinPts
                IdsConvert = posGrande[posGrande!=0] / pos[pos!=0]
                models.control_h[0][pos[pos!=0].astype(int).tolist()] = IdsConvert
        return IdsConvert,xyNew


    def set_sediments(self,var,VarName, wi = [0.036, 2.2e-4, 8.6e-7],
        diametro = [0.35, 0.016, 0.001], G = 9.8):
        'Descripcion: Alojas las variables requeridas para la ejecucion\n'\
        '   del modelo de sedimentos.\n'\
        '\n'\
        'Parametros\n'\
        '----------\n'\
        'var : Variable que describe la propiedad (constante, ruta, mapa o vector).\n'\
        'varName: Nombre de la variable a ingresar en el modelo.\n'\
        '   Krus : Erosividad del suelo (RUSLE).\n'\
        '   Crus : Cobertura del suelo (RUSLE).\n'\
        '   Prus : Practicas proteccion (RUSLE).\n'\
        '   PArLiAc : Porcentaje Arenas, Limos y Arcillas.\n'\
        '   wi : velocidad de caida de cada particula de sed [m/s].\n'\
        '   diametro : diametro de cada tipo de sedimento [mm].\n'\
        '\n'\
        'Retornos\n'\
        '----------\n'\
        'self : variables iniciadas en el modelo bajo los nombres de:.\n'\
        '   wmf.models.krus.\n'\
        '   wmf.models.crus.\n'\
        '   wmf.models.prus.\n'\
        '   wmf.models.parliac.\n'\
        '   wmf.models.wi.\n'\
        '   wmf.models.diametro.\n'\
                #Determina el tipo de unidades del modelo
        if self.modelType[0] is 'c':
            N = self.ncells
        elif self.modelType[0] is 'h':
            N = self.nhills
                #Se fija que tipo de variable es
        isVec=False
        if type(var) is str:
            #Si es un string lee el mapa alojado en esa ruta
            Map,Pp = read_map_raster(var)
            Vec = self.Transform_Map2Basin(Map,Pp)
            isVec=True
        elif type(var) is int or float:
            Vec = np.ones((1,self.ncells))*var
            isVec=True
        elif type(var) is np.ndarray and var.shape[0] == self.ncells:
            Vec = var
            isVec=True
        #Si el modelo es tipo ladera agrega la variable
        if self.modelType[0] is 'h':
            Vec = self.Transform_Basin2Hills(Vec,mask=mask)
        #Inicia las variables
        if VarName == 'Krus':
            models.krus = np.ones((1,N))*Vec
        if VarName == 'Prus':
            models.prus = np.ones((1,N))*Vec
        if VarName == 'Crus':
            models.crus = np.ones((1,N))*Vec
        if VarName == 'PArLiAc':
            models.parliac = np.ones((3,N))*Vec
        #Variables de diametro y velocidad de caida
        models.wi = wi
        models.diametro = diametro
        models.g = G

    def set_Slides(self,var,VarName):
        'Descripcion: Alojas las variables requeridas para la ejecucion\n'\
        '   del modelo de deslizamientos.\n'\
        '\n'\
        'Parametros\n'\
        '----------\n'\
        'var : Variable que describe la propiedad (constante, ruta, mapa o vector).\n'\
        'varName: Nombre de la variable a ingresar en el modelo.\n'\
        '   Zs : Profundidad del suelo.\n'\
        '   GammaSoil: Densidad del suelo.\n'\
        '   Cohesion: Cohesion del suelo.\n'\
        '   FrictionAngle : Angulo de friccion del suelo.\n'\
        '   FS: Factor de Seguridad, en este caso se envia una constante.\n'\
        '   RadSlope : .\n'\
        '\n'\
        'Retornos\n'\
        '----------\n'\
        'self : variables iniciadas en el modelo bajo los nombres de:.\n'\
        '   wmf.models.sl_zs.\n'\
        '   wmf.models.sl_gammas.\n'\
        '   wmf.models.sl_cohesion.\n'\
        '   wmf.models.sl_frictionangle.\n'\
        '   wmf.models.sl_radslope.\n'\
        '   wmf.models.sl_fs.\n'\
        #Pone el gamma del agua por defecto
        models.sl_gammaw = 9.8
        #Obtiene el vector que va a alojar en el modelo
        if VarName != 'FS':
            isVec=False
            if type(var) is str:
                #Si es un string lee el mapa alojado en esa ruta
                Map,Pp = read_map_raster(var)
                Vec = self.Transform_Map2Basin(Map,Pp)
                isVec=True
            elif type(var) is int or float:
                Vec = np.ones((1,self.ncells))*var
                isVec=True
            elif type(var) is np.ndarray and var.shape[0] == self.ncells:
                Vec = var
                isVec=True
        #Si el modelo es tipo ladera agrega la variable
        if self.modelType[0] == 'h':
            return 'El modelo por laderas no simula deslizamientos.'
        #finalmente mete la variable en el modelo
        N = self.ncells
        if VarName is 'GammaSoil' :
            models.sl_gammas = np.ones((1,N))*Vec
        elif VarName is 'Cohesion':
            models.sl_cohesion = np.ones((1,N))*Vec
        elif VarName is 'FrictionAngle':
            models.sl_frictionangle = np.ones((1,N))*np.deg2rad(Vec)
        elif VarName is 'Zs':
            models.sl_zs = np.ones((1,N))*Vec
        elif VarName is 'FS':
            models.sl_fs = var
        elif VarName is 'Slope':
            models.sl_radslope = np.ones((1,N))*np.arctan(Vec)
            models.sl_radslope[models.sl_radslope == 0] = 0.01
    #------------------------------------------------------
    # Guardado y Cargado de modelos de cuencas preparados
    #------------------------------------------------------
    def Save_SimuBasin(self,ruta,SimSlides = False,
        ExtraVar = None):
        'Descripcion: guarda una cuenca previamente ejecutada\n'\
        '\n'\
        'Parametros\n'\
        '----------\n'\
        'ruta : Ruta donde la cuenca sera guardada.\n'\
        'ruta_dem : direccion donde se aloja el DEM (se recomienda absoluta).\n'\
        'ruta_dir : direccion donde se aloja el DIR (se recomienda absoluta).\n'\
        'SimSlides: indica a la funcion si va a guardar o no informacion para la simulacion.\n'\
        '   de deslizamientos.\n'\
        'ExtraVar: Variables extras de simulacion deben ir en un diccionario.\n'\
        '   Forma del diccionario Dict = {"varName": {"Data": vector[ncells], "type": "tipo"}}.\n'\
        '   Los tipos de variables son: flotante: "f4", entero "i4".\n'\
        '\n'\
        'Retornos\n'\
        '----------\n'\
        'self : Con las variables iniciadas.\n'\

        #Si esta o no set el Geomorphology, de acuerdo a eso lo estima por defecto  
        if self.isSetGeo is False:
            self.set_Geomorphology()
            print('Aviso: SE ha estimado la geomorfologia con los umbrales por defecto umbral = [30, 500]')
        #Guarda la cuenca
        if self.modelType[0] is 'c':
            N = self.ncells
        elif self.modelType[0] is 'h':
            N = self.nhills

        Dict = {'nombre':self.name,
                    'modelType':self.modelType,'noData':cu.nodata,'umbral':self.umbral,
                    'ncells':self.ncells,'nhills':self.nhills,
                    'dt':models.dt,'Nelem':N,'dxp':cu.dxp,'retorno':models.retorno,
                    'storageConst' :models.storage_constant,
                    'ncols':cu.ncols,
                    'nrows':cu.nrows,
                    'xll':cu.xll,
                    'yll':cu.yll,
                    'dx':cu.dx,
                    'dy':cu.dy,
                    'epsg': self.epsg}
        if SimSlides:
            Dict.update({'sl_fs':models.sl_fs, 'sl_gullie':models.sl_gullienogullie, 'sl_gammaw':models.sl_gammaw})
        #abre el archivo
        gr = netcdf.Dataset(ruta,'w',format='NETCDF4')
        #Grupo base
        GrupoBase = gr.createGroup('base')
        GrupoSimHid = gr.createGroup('SimHidro')
        GrupoSimSed = gr.createGroup('SimSediments')
        GrupoSimSli = gr.createGroup('SimSlides')
        GrupoHidro = gr.createGroup('Hidro')
        GrupoGeo = gr.createGroup('Geomorfo')
        GrupoCalib = gr.createGroup('Parametros')
        #Variables grupo base
        DimNcell = GrupoBase.createDimension('ncell',self.ncells)
        DimNhill = GrupoBase.createDimension('nhills',self.nhills)
        DimCol3 = GrupoBase.createDimension('col3',3)
        DimCol2 = GrupoBase.createDimension('col2',2)
        VarStruc = GrupoBase.createVariable('structure','i4',('col3','ncell'),zlib=True)
        VarHills = GrupoBase.createVariable('hills','i4',('col2','nhills'),zlib=True)
        VarHills_own = GrupoBase.createVariable('hills_own','i4',('ncell',),zlib=True)
        VarStruc[:] = self.structure
        VarHills[:] = self.hills
        VarHills_own[:] = self.hills_own
        #Variables de Simulacion hidrologica
        DimNelem = GrupoSimHid.createDimension('Nelem',N)
        DimCol3 = GrupoSimHid.createDimension('col3',3)
        DimCol2 = GrupoSimHid.createDimension('col2',2)
        DimCol4 = GrupoSimHid.createDimension('col4',4)
        DimCol5 = GrupoSimHid.createDimension('col5',5)
        VarH_coef = GrupoSimHid.createVariable('h_coef','f4',('col4','Nelem'),zlib=True)
        VarV_coef = GrupoSimHid.createVariable('v_coef','f4',('col4','Nelem'),zlib=True)
        VarH_exp = GrupoSimHid.createVariable('h_exp','f4',('col4','Nelem'),zlib=True)
        VarV_exp = GrupoSimHid.createVariable('v_exp','f4',('col4','Nelem'),zlib=True)
        Var_H1max = GrupoSimHid.createVariable('h1_max','f4',('Nelem',),zlib = True)
        Var_H3max = GrupoSimHid.createVariable('h3_max','f4',('Nelem',),zlib = True)
        Control = GrupoSimHid.createVariable('control','i4',('Nelem',),zlib = True)
        ControlH = GrupoSimHid.createVariable('control_h','i4',('Nelem',),zlib = True)
        if self.modelType[0] is 'c':
            drena = GrupoSimHid.createVariable('drena','i4',('col3','Nelem'),zlib = True)
        elif self.modelType[0] is 'h':
            drena = GrupoSimHid.createVariable('drena','i4',('Nelem'),zlib = True)
        unitType = GrupoSimHid.createVariable('unit_type','i4',('Nelem',),zlib = True)
        hill_long = GrupoSimHid.createVariable('hill_long','f4',('Nelem',),zlib = True)
        hill_slope = GrupoSimHid.createVariable('hill_slope','f4',('Nelem',),zlib = True)
        stream_long = GrupoSimHid.createVariable('stream_long','f4',('Nelem',),zlib = True)
        stream_slope = GrupoSimHid.createVariable('stream_slope','f4',('Nelem',),zlib = True)
        stream_width = GrupoSimHid.createVariable('stream_width','f4',('Nelem',),zlib = True)
        elem_area = GrupoSimHid.createVariable('elem_area','f4',('Nelem',),zlib = True)
        speed_type = GrupoSimHid.createVariable('speed_type','i4',('col3',),zlib = True)
        storage = GrupoSimHid.createVariable('storage','i4',('col5','Nelem'),zlib = True)
        VarH_coef[:] = models.h_coef
        VarV_coef[:] = models.v_coef
        VarH_exp[:] = models.h_exp
        VarV_exp[:] = models.v_exp
        Var_H1max[:] = models.max_capilar
        Var_H3max[:] = models.max_gravita
        Control[:] = models.control
        ControlH[:] = models.control_h
        drena[:] = models.drena
        unitType[:] = models.unit_type
        hill_long[:] = models.hill_long
        hill_slope[:] = models.hill_slope
        stream_long[:] = models.stream_long
        stream_slope[:] = models.stream_slope
        stream_width[:] = models.stream_width
        elem_area[:] = models.elem_area
        speed_type[:] = models.speed_type
        storage[:] = models.storage
        #Variables grupo GEomorfologia
        DimNcell = GrupoGeo.createDimension('ncell',self.ncells)
        VarDEM = GrupoGeo.createVariable('DEM','f4',('ncell',),zlib = True)
        VarDIR = GrupoGeo.createVariable('DIR','i4',('ncell',),zlib = True)
        VarDEM[:] = self.DEMvec
        VarDIR[:] = self.DIRvec
        #Variables de sedimentos
        DimNelem = GrupoSimSed.createDimension('Nelem',N)
        DimCol3 = GrupoSimSed.createDimension('col3',3)
        VarKrus = GrupoSimSed.createVariable('Krus','f4',('Nelem'),zlib=True)
        VarPrus = GrupoSimSed.createVariable('Prus','f4',('Nelem'),zlib=True)
        VarCrus = GrupoSimSed.createVariable('Crus','f4',('Nelem'),zlib=True)
        VarParliac = GrupoSimSed.createVariable('PArLiAc','f4',('col3','Nelem'),zlib=True)
        def __reshape (Variable,dim):
            try:
                if Variable==None:
                    new_var = np.ones((dim,N))
            except:
                if len(Variable) != N:
                    new_var = np.ones((dim,N))
                else:
                    new_var = Variable
            return new_var
                    
        VarKrus[:]=__reshape(models.krus,1)
        VarPrus[:]=__reshape(models.prus,1)
        VarCrus[:]=__reshape(models.crus,1)
        VarParliac[:]=__reshape(models.parliac,3)
        #Variables de deslizamientos
        if SimSlides:
            DimNelem = GrupoSimSli.createDimension('Nelem',N)
            frictionAngle = GrupoSimSli.createVariable('friction_angle','f4',('Nelem',),zlib = True)
            Cohesion = GrupoSimSli.createVariable('cohesion','f4',('Nelem',),zlib = True)
            GammaSoil = GrupoSimSli.createVariable('gamma_soil','f4',('Nelem',),zlib = True)
            ZSoil = GrupoSimSli.createVariable('z_soil','f4',('Nelem',),zlib = True)
            RadSlope = GrupoSimSli.createVariable('rad_slope','f4',('Nelem',),zlib = True)
            frictionAngle[:] = models.sl_frictionangle
            Cohesion[:] = models.sl_cohesion
            GammaSoil[:] = models.sl_gammas
            ZSoil[:] = models.sl_zs
            RadSlope[:] = models.sl_radslope
        #Introduce variables extras en caso de que el usuario las incluyera
        if type(ExtraVar) is dict:
            for k in ExtraVar.keys():
                Var = gr.createVariable(k,ExtraVar[k]['type'],('ncell',),zlib=True)
                Var[:] = ExtraVar[k]['Data']
        #asigna las prop a la cuenca
        gr.setncatts(Dict)
        #Cierra el archivo
        gr.close()
        #Sale del programa
        self.segunda_cuenca = True 
        return

        #------------------------------------------------------
        # Ejecucion del modelo
        #------------------------------------------------------
    def run_shia(self,Calibracion,
        rain_rute, N_intervals, start_point = 1, StorageLoc = None, HspeedLoc = None,ruta_storage = None, ruta_speed = None,
        ruta_conv = None, ruta_stra = None, ruta_retorno = None,kinematicN = 5, QsimDataFrame = True, EvpVariable = False):
        'Descripcion: Ejecuta el modelo una ves este es preparado\n'\
        '   Antes de su ejecucion se deben tener listas todas las . \n'\
        '   variables requeridas . \n'\
        '\n'\
        'Parametros\n'\
        '----------\n'\
        'self : Cuenca a ejecutar con todo listo para ser ejecutada.\n'\
        'Calibracion : Parametros de calibracion del modelo, orden:.\n'\
        '   - Evaporacion.\n'\
        '   - Infiltracion.\n'\
        '   - Percolacion.\n'\
        '   - Perdidas.\n'\
        '   - Vel Superficial .\n'\
        '   - Vel Sub-superficial.\n'\
        '   - Vel Subterranea.\n'\
        '   - Vel Cauce.\n'\
        '   - Max Capilar.\n'\
        '   - Max Gravitacional.\n'\
        'rain_rute : Ruta donde se encuentra el archivo binario de lluvia:.\n'\
        '   generado por rain_interpolate_* o por rain_radar2basin.\n'\
        'N_intervals : Numero de intervalos de tiempo.\n'\
        'start_point : Punto donde comienza a usar registros de lluvia.\n'\
        '   los binarios generados por rain_* generan un archivo de texto.\n'\
        '   que contiene fechas par aayudar a ubicar el punto de inicio deseado.\n'\
        'StorageLoc: Variable local de almacenamiento, esto es en el caso de que no se.\n'\
        '   desee usar la configuracion global de almacenamiento del modelo (5, N).\n'\
        'HspeedLoc: Variable local para setear la velocidad horizontal inicial de ejecucion.\n'\
        'ruta_storage : Ruta donde se guardan los estados del modelo en cada intervalo.\n'\
        '   de tiempo, esta es opcional, solo se guardan si esta variable es asignada.\n'\
        'ruta_conv : Ruta al binario y hdr indicando las nubes que son convectivas.\n'\
        'ruta_stra : Ruta al binario y hdr indicando las nubes que son estratiformes.\n'\
        'ruta_retorno : Ruta al binario y hdr en donde escribe la serie con los milimetros retornados al tanque runoff.\n'\
        'kinematicN: Cantidad de iteraciones para la solucion de la onda cinematica.\n'\
        '   De forma continua: 5 iteraciones, recomendado para cuando el modelo se\n'\
        '       ejecuta en forma continua Ej: cu.run_shia(Calib, rain_rute, 100)\n'\
        '   Por intervalos: 10 iteraciones, recomendado cuando el modelo se ejecuta\n'\
        '       por intervalos de un paso ej:\n'\
        '       for i in range(1,N)\n'\
        '           Results = cu.run_shia(Calib, ruta_rain, 1, i)\n'\
        '           for c,j in enumerate(Results[''Storage'']):\n'\
        '               cu.set_Storage(j,c)\n'\
        'QsimDataFrame: Retorna un data frame con los caudales simulados indicando su id de acuerdo con el\n'\
        '   que guarda la funcion Save_Net2Map con la opcion NumTramo = True. \n'\
                'EvpVariable: (False) Asume que la evp del modelo cambia en funcion o no de la radiacion\n'\
        '\n'\
        'Retornos\n'\
        '----------\n'\
        'Qsim : Caudal simulado en los puntos de control.\n'\
        'Hsim : Humedad simulada en los puntos de control.\n'\
                #genera las rutas
        rain_ruteBin,rain_ruteHdr = __Add_hdr_bin_2route__(rain_rute)
        #Obtiene las fechas
        Rain = read_mean_rain(rain_ruteHdr, N_intervals, start_point)
        # De acuerdo al tipo de modelo determina la cantidad de elementos
        if self.modelType[0] is 'c':
            N = self.ncells
        elif self.modelType[0] is 'h':
            N = self.nhills
        #prepara variables globales
        models.rain_first_point = start_point
        models.calc_niter = kinematicN
        #Prepara terminos para control
        if np.count_nonzero(models.control) == 0 :
            NcontrolQ = 1
        else:
            NcontrolQ = np.count_nonzero(models.control)+1
        if np.count_nonzero(models.control_h) == 0 :
            NcontrolH = 1
        else:
            NcontrolH = np.count_nonzero(models.control_h)
        #Prepara variables de guardado de almacenamiento
        if ruta_storage is not None:
            models.save_storage = 1
            ruta_sto_bin, ruta_sto_hdr = __Add_hdr_bin_2route__(ruta_storage,
                storage = True)
        else:
            models.save_storage = 0
            ruta_sto_bin = 'no_guardo_nada.StObin'
            ruta_sto_hdr = 'no_guardo_nada.StOhdr'
        #prepara variable para guardado de velocidad
        if ruta_speed  is not  None:
            models.save_speed = 1
            ruta_speed_bin, ruta_speed_hdr = __Add_hdr_bin_2route__(ruta_speed,
                storage = True)
        else:
            models.save_speed = 0
            ruta_speed_bin = 'no_guardo_nada.bin'
            ruta_speed_hdr = 'no_guardo_nada.hdr'
        #Prepara variables para el almacenamiento de retorno
        if ruta_retorno is not None:
            models.save_retorno = 1
            ruta_ret_bin, ruta_ret_hdr = __Add_hdr_bin_2route__(ruta_retorno)
        else:
            models.save_retorno = 0
            ruta_ret_bin = 'no_guardo_nada.bin'
            ruta_ret_hdr = 'no_guardo_nada.hdr'
        #Variables de separacion de flujo por tipo de lluvia
        if ruta_conv  is not  None and ruta_stra  is not  None:
            ruta_binConv,ruta_hdrConv = __Add_hdr_bin_2route__(ruta_conv)
            ruta_binStra,ruta_hdrStra = __Add_hdr_bin_2route__(ruta_stra)
            models.separate_rain = 1
        else:
            models.separate_rain = 0
            ruta_binConv = 'none'
            ruta_hdrConv = 'none'
            ruta_binStra = 'none'
            ruta_hdrStra = 'none'
        #Prepara la variable de almacenamiento
        if StorageLoc is not None:
            if StorageLoc.shape != (5, N):
                print('Error: almacenamiento debe ser: (5,'+str(N)+'), y es: ('+str(StorageLoc.shape[0])+','+str(StorageLoc.shape[1])+')')
                StorageLoc = np.zeros((5,N))*-9999.0
        else:
            StorageLoc = np.zeros((5,N))*-9999.0
        #Prepara la variable de velocidad horizontal inicial
        if HspeedLoc is not None:
            if HspeedLoc.shape != (4, N):
                print('Error: velocidad debe ser: (4,'+str(N)+'), y es: ('+str(HspeedLoc.shape[0])+','+str(HspeedLoc.shape[1])+')')
                HspeedLoc = np.zeros((4,N))*-9999.0
        else:
            HspeedLoc = np.zeros((4,N))*-9999.0
        #Implementa o no la EVP variable en funcion de la radiacion
        if EvpVariable:
            Rad = self.__GetEVP_Serie__(Rain.index)
        else:
            models.evpserie = np.ones(N_intervals)
        # Ejecuta el modelo
        Qsim,Qsed,Qseparated,Humedad,St1,St3,Balance,Speed,Area,Alm,Qsep_byrain = models.shia_v1(
            rain_ruteBin,
            rain_ruteHdr,
            Calibracion,
            #N,
            NcontrolQ,
            NcontrolH,
            N_intervals,
            StorageLoc,
            HspeedLoc,
            N,
            ruta_sto_bin,
            ruta_speed_bin,
            ruta_binConv,
            ruta_binStra,
            ruta_hdrConv,
            ruta_hdrStra,
            ruta_ret_bin)
        #Retorno de variables de acuerdo a lo simulado
        Retornos={'Qsim' : Qsim}
        Retornos.update({'Balance' : Balance})
        Retornos.update({'Storage' : Alm})
        if np.count_nonzero(models.control_h)>0:
            Retornos.update({'Humedad' : Humedad})
            Retornos.update({'Humedad_t1' : St1})
            Retornos.update({'Humedad_t2' : St3})
        if models.sim_sediments == 1:
            Retornos.update({'Sediments' : Qsed})
            #print (Qsed)
        if models.separate_fluxes == 1:
            Retornos.update({'Fluxes' : Qseparated})
        if models.separate_rain == 1:
            Retornos.update({'Rain_sep' : Qsep_byrain})
        if models.show_storage == 1:
            Retornos.update({'Mean_Storage' : np.copy(models.mean_storage)})
        if models.save_storage == 1:
            rutaStorageHdr = __Add_hdr_bin_2route__(ruta_storage)
            #Caso en el que se registra el alm medio
            if models.show_storage == 1:
                __Save_storage_hdr__(ruta_sto_hdr,rain_ruteHdr,N_intervals,
                    start_point,self,Mean_Storage = np.copy(models.mean_storage))
            #Caso en el que no hay alm medio para cada uno de los
            else:
                __Save_storage_hdr__(ruta_sto_hdr,rain_ruteHdr,N_intervals,
                    start_point,self,Mean_Storage=np.zeros((5,N))*-9999)
        #Area de la seccion
        if models.show_area == 1:
            Retornos.update({'Sec_Area': Area})
        #Velocidades
        if models.show_speed == 1:
            Retornos.update({'Speed' : Speed})
        if models.show_mean_speed == 1:
            Retornos.update({'Mean_Speed' : np.copy(models.mean_speed)})
        if models.save_speed == 1:
            rutaSpeedHdr = __Add_hdr_bin_2route__(ruta_speed)
            #Caso en el que hay velocidad media para todos los tanques
            if models.show_mean_speed == 1:
                __Save_speed_hdr__(ruta_speed_hdr,rain_ruteHdr,N_intervals,
                    start_point,self,Mean_Speed = np.copy(models.mean_speed))
            #Caso en el que no hay alm medio para cada uno de los
            else:
                __Save_speed_hdr__(ruta_speed_hdr,rain_ruteHdr,N_intervals,
                    start_point,self)

        if models.save_retorno == 1:
            if models.show_mean_retorno == 1:
                __Save_retorno_hdr__(ruta_ret_hdr, rain_ruteHdr, N_intervals,
                    start_point, self, Mean_retorno = models.mean_retorno)
            else:
                __Save_retorno_hdr__(ruta_ret_hdr, rain_ruteHdr, N_intervals,
                    start_point, self)

        #Campo de lluvia acumulado para el evento
        Retornos.update({'Rain_Acum': models.acum_rain})
        Retornos.update({'Rain_hietogram': models.mean_rain})
        #Retornos en caso de simular deslizamientos
        if models.sim_slides == 1:
            Retornos.update({'Slides_Map': np.copy(models.sl_slideocurrence)})
            Retornos.update({'Slides_NCell_Serie': np.copy(models.sl_slidencelltime)})
            Retornos.update({'Slides_Acum':np.copy(models.sl_slideacumulate)})
        #Caudal simulado en un dataframe
        if QsimDataFrame:
            #Obtiene ids
            ids = models.control[models.control!=0]
            Qdict = {}
            for i,j in zip(Retornos['Qsim'][1:], ids):
                Qdict.update({str(j): i})
            Qdict = pd.DataFrame(Qdict, index=Rain.index)
            #Si separa flujos, entrega el caudal tambien en data frame por flujos
            if models.separate_fluxes == 1:
                Qsep = []
                tupla = []
                for z,i,j in zip(Retornos['Qsim'][1:], Retornos['Fluxes'][1:], ids):
                    tupla.append((str(j),'run'))
                    tupla.append((str(j),'sub'))
                    tupla.append((str(j),'sup'))
                    Qsep.extend([i[0],i[1],z-i[0]-i[1]])
                index = pd.MultiIndex.from_tuples(tupla, names=['reach','flux'])
                Qsep = np.array(Qsep)
                QsepDict = pd.DataFrame(Qsep.T, index=Rain.index, columns=index)
                #Si simula sedimentos, hace el dataframe
            if models.sim_sediments == 1:
                Qsedi = []
                tupla = []
                for z,i,j in zip(Retornos['Qsim'][1:], Retornos['Sediments'][1:], ids):
                    tupla.append((str(j),'Sand'))
                    tupla.append((str(j),'Lime'))
                    tupla.append((str(j),'Clay'))
                    Qsedi.extend([i[0],i[1],i[2]])#  [i[0],i[1],z-i[0]-i[1]])
                index = pd.MultiIndex.from_tuples(tupla, names=['reach','Sediments'])
                Qsedi = np.array(Qsedi)
                QsediDict = pd.DataFrame(Qsedi.T, index=Rain.index, columns=index)
            if models.separate_fluxes == 1 and models.sim_sediments == 0:
                return Retornos, Qdict, QsepDict
            if models.separate_fluxes == 1 and models.sim_sediments == 1:
                return Retornos, Qdict, QsepDict, QsediDict
            if models.separate_fluxes == 0 and models.sim_sediments == 1:
                return Retornos, Qdict, QsediDict
            return Retornos, Qdict
        return Retornos

    def efficiencia(self, Qobs, Qsim):
        'Descripcion: Calcula diferentes indices de desempeno del modelo\n'\
        '   nash, qpico, rmse, rmseLog, t_pico\n'\
        '\n'\
        'Parametros\n'\
        '----------\n'\
        'Qobs : Caudal observado.\n'\
        'Qsim : Caudal simulado.\n'\
        'Retornos\n'\
        '----------\n'\
        'DicEff: diccionario con la eficiencia obtenida en cada param.\n'\
        #Comienza a evaluar cada parametro
        DictEff = {'Nash': __eval_nash__(Qobs, Qsim)}
        DictEff.update({'Qpico': __eval_q_pico__(Qobs, Qsim)})
        DictEff.update({'Tpico': __eval_t_pico__(Qobs, Qsim, models.dt)})
        DictEff.update({'RmseLog': __eval_rmse_log__(Qobs, Qsim)})
        DictEff.update({'Rmse': __eval_rmse__(Qobs, Qsim)})
        return DictEff

    def Calib_NSGAII(self, nsga_el, nodo_eval, pop_size = 40, process = 4,
        NGEN = 6, MUTPB = 0.5, CXPB = 0.5):
        'Descripcion: Algoritmo para calibrar de forma automatica el modelo\n'\
        '   hidrologico utilizando DEAP y su funcion de seleccion NSGAII.\n'\
        '\n'\
        'Parametros\n'\
        '----------\n'\
        '   - nsga_el : un objeto de la clase nsgaii_element, el cual determinara las reglas.\n'\
        '       para la generacion de calibraciones.\n'\
        '   - nodo_eval: nodo donde se evalua el modelo.\n'\
        '   - pop_size: tamano de la poblacion (cantidad de calibraciones a evaluar, multiplo de 4).\n'\
        '   - process: Cantidad de procesadores que se utilizaran en la generacion.\n'\
        '   - NGEN: Cantidad de generaciones para obtener la poblacion final.\n'\
        '   - MUTPB: Probabilidad generica de que un gen mute.\n'\
        '   - CXPB: Probabilidad generica de que dos genes se crucen.\n'\
        '\n'\
        'Retornos\n'\
        '----------\n'\
        '   - pop: poblacion final.\n'\
        '   - Qsim: Caudales simulados en el punto evaluado.\n'\
        '   - fitness: desempeno de la poblacion obtenida.\n'\
        #Check de que la poblacion sea multiplo de 4
        Flag = True
        while Flag:
            if np.mod(pop_size,4) == 0:
                Flag = False
            else:
                pop_size += 1
        #Inicia el elemento nsga con los parametros de el
        nsga_el.set_nsgaII()
        #Crea las poblaciones y las ejecuciones
        pop = nsga_el.toolbox.population(pop_size)
        Ejecs = map(nsga_el.__crea_ejec__, pop)
        #Ejecuta a la poblacion
        QsimPar = __ejec_parallel__(Ejecs, process, nodo_eval)
        fitnesses = map(nsga_el.toolbox.evaluate, QsimPar)
        for ind, fit in zip(pop, fitnesses):
            ind.fitness.values = fit
        #Itera la poblacion hasta encontrar a la mejor
        #toolbox = base.Toolbox()
        for g in range(NGEN):
            offspring = tools.selTournamentDCD(pop, len(pop))
            offspring = map(nsga_el.toolbox.clone, offspring)
            # Apply crossover and mutation on the offspring
            for child1, child2 in zip(offspring[::2], offspring[1::2]):
                if random.random() < CXPB:
                    nsga_el.toolbox.mate(child1[0], child2[0])
                    del child1.fitness.values
                    del child2.fitness.values
            #Muta a algunos de los genomas
            for mutant in offspring:
                if random.random() < MUTPB:
                    nsga_el.toolbox.mutate(mutant[0])
                    del mutant.fitness.values
            #Ejecuta a la nueva generacion
            Ejecs = map(nsga_el.__crea_ejec__, offspring)
            QsimPar = __ejec_parallel__(Ejecs, process, nodo_eval)
            #Identifica a la gente que no cumple de la generacion
            Qsim_invalid = []
            invalid_ind = []
            for i,q in zip(offspring, QsimPar):
                if not i.fitness.valid:
                    Qsim_invalid.append(q)
                    invalid_ind.append(i)
            #Evaluate the individuals with an invalid fitness
            fitnesses = map(nsga_el.toolbox.evaluate, Qsim_invalid)
            for ind, fit in zip(invalid_ind, fitnesses):
                ind.fitness.values = fit
            #Toma la siguiente generacion
            pop = nsga_el.toolbox.select(pop + offspring, pop_size)
        #Retorno
        return pop, QsimPar, np.array(fitnesses).T


class nsgaii_element:
    def __init__(self, rutaLluvia, Qobs, npasos, inicio, SimuBasinElem ,evp =[0,1], infil = [1,200], perco = [1, 40],
        losses = [0,1],velRun = [0.1, 1], velSub = [0.1, 1], velSup =[0.1, 1],
        velStream = [0.1, 1], Hu = [0.1, 1], Hg = [0.1, 1],
        probCruce = np.ones(10)*0.5, probMutacion = np.ones(10)*0.5,
        rangosMutacion = [[0,1], [1,200], [1,40], [0,1], [0.1,1], [0.1, 1], [0.1,1], [0.1,1], [0.1, 1], [0.1,1]],
        MaxMinOptima = (1.0, -1.0), CrowDist = 0.5):
        'Descripcion: Inicia el objeto de calibracion genetica tipo NSGAII\n'\
        '   este objeto contiene las reglas principales para la implementacion\n'\
        '   de todo el algoritmo de calibracion genetico.\n'\
        '\n'\
        'Parametros\n'\
        '----------\n'\
        '   -rutaLluvia : Ruta donde se encuentra el archivo de lluvia binario.\n'\
        '   -Qobs: Array numpy con el caudal observado [npasos].\n'\
        '   -npasos: Cantidad de pasos en simulacion.\n'\
        '   -inicio: Punto de inicio en la simulacion.\n'\
        '   -SimuBasinElem: Objeto de simulacion.\n'\
        'Retornos\n'\
        '----------\n'\
        'self : Con las variables iniciadas.\n'\
        #Rangos parametros
        self.evp_range = evp
        self.infil_range = infil
        self.perco_range = perco
        self.losses_range = losses
        self.velRun_range = velRun
        self.velSub_range = velSub
        self.velSup_range = velSup
        self.velStream_range = velStream
        self.hu_range = Hu
        self.hg_range = Hg
        #Establece propiedades para ejecucion
        self.npasos = npasos
        self.inicio = inicio
        self.ruta_lluvia = rutaLluvia
        self.Qobs = Qobs
        self.simelem = SimuBasinElem
        #Propiedades de cruce y mutacion
        self.prob_cruce = probCruce
        self.prob_mutacion = probMutacion
        self.rangos_mutacion = rangosMutacion
        self.optimiza = MaxMinOptima
        self.crowdist = CrowDist

    def __crea_calibracion__(self):
        #Evp
        if self.evp_range[0] != self.evp_range[1]:
            evp = np.random.uniform(self.evp_range[0], self.evp_range[1],1)[0]
        else:
            evp = self.evp_range[0]
        #Infiltracion
        if self.infil_range[0] != self.infil_range[1]:
            infil = np.random.uniform(self.infil_range[0], self.infil_range[1],1)[0]
        else:
            infil = self.infil_range[0]
        #Percolacion
        if self.perco_range[0] != self.perco_range[1]:
            perco = np.random.uniform(self.perco_range[0], self.perco_range[1],1)[0]
        else:
            perco = self.perco_range[0]
        #Perdidas
        if self.losses_range[0] != self.losses_range[1]:
            losses = np.random.uniform(self.losses_range[0], self.losses_range[1],1)[0]
        else:
            losses = self.losses_range[0]
        #runoff
        if self.velRun_range[0] != self.velRun_range[1]:
            velRun = np.random.uniform(self.velRun_range[0], self.velRun_range[1],1)[0]
        else:
            velRun = self.velRun_range[0]
        #Vel subsuperficial
        if self.velSub_range[0] != self.velSub_range[1]:
            velSub = np.random.uniform(self.velSub_range[0], self.velSub_range[1],1)[0]
        else:
            velSub = self.velSub_range[0]
        #Vel acuifero
        if self.velSup_range[0] != self.velSup_range[1]:
            velSup = np.random.uniform(self.velSup_range[0], self.velSup_range[1],1)[0]
        else:
            velSup = self.velSup_range[0]
        #Vel cauce
        if self.velStream_range[0] != self.velStream_range[1]:
            velStream = np.random.uniform(self.velStream_range[0], self.velStream_range[1],1)[0]
        else:
            velStream = self.velStream_range[0]
        #almacenamiento hu
        if self.hu_range[0] != self.hu_range[1]:
            hu = np.random.uniform(self.hu_range[0], self.hu_range[1],1)[0]
        else:
            hu = self.hu_range[0]
        #almacenamiento hg
        if self.hg_range[0] != self.hg_range[1]:
            hg = np.random.uniform(self.hg_range[0], self.hg_range[1],1)[0]
        else:
            hg = self.hg_range[0]
        #Retorna una calibracion aleatoria
        return [evp, infil, perco, losses, velRun, velSub, velSup, velStream, hu, hg]

    def __crea_ejec__(self, calibracion):
        return [calibracion, self.ruta_lluvia, self.npasos, self.inicio, self.simelem]

    def __evalfunc__(self, Qsim, f1 = __eval_nash__, f2 = __eval_q_pico__):
        E1 = f1(self.Qobs, Qsim)
        E2 = f2(self.Qobs, Qsim)
        return E1, E2

    def __cruce__(self, indi1, indi2):
        for i,u in zip(range(10), self.prob_cruce):
            p = np.random.uniform(0,1,1)
            if p>u:
                a = indi1[i]; b = indi2[i]
                indi1[i] = b
                indi2[i] = a
        return indi1, indi2

    def __mutacion__(self, indi):
        c = 0
        for i,u in zip(range(10), self.prob_mutacion):
            p = np.random.uniform(0,1,1)
            if p>u:
                indi[i] = np.random.uniform(self.rangos_mutacion[c][0],self.rangos_mutacion[c][0],1)[0]
            c+=1
        return indi

    def set_nsgaII(self):
        self.toolbox = base.Toolbox()
        creator.create("FitnessMin", base.Fitness,
            weights = self.optimiza,
            crowding_dist = self.crowdist)
        creator.create("Individual", list, fitness=creator.FitnessMin)
        self.toolbox.register("attr1", self.__crea_calibracion__)
        self.toolbox.register("individual", tools.initRepeat, creator.Individual,
            self.toolbox.attr1, n=1)
        self.toolbox.register("population", tools.initRepeat, list, self.toolbox.individual)
        self.toolbox.register("evaluate", self.__evalfunc__)
        self.toolbox.register("mate", self.__cruce__)
        self.toolbox.register("mutate", self.__mutacion__)
        self.toolbox.register("select", tools.selNSGA2)

class Stream:
    #------------------------------------------------------
    # Subrutinas de trazado de corriente y obtencion de parametros
    #------------------------------------------------------
    #Inicia la cuenca
    def __init__(self,lat,lon,DEM,DIR,name='NaN'):
        'Descripcion: Traza un cauce e inicia la variable de este \n'\
        '\n'\
        'Parametros\n'\
        '----------\n'\
        'self : Inicia las variables vacias.\n'\
        'lat : Coordenada en X del punto mas alto del cauce.\n'\
        'lon : Coordenada en Y del punto mas alto del cauce.\n'\
        '\n'\
        'Retornos\n'\
        '----------\n'\
        'self : Con las variables iniciadas y estructura del cauce.\n'\
        #Realiza copia de los mapas y obtiene el cauce
        self.DEM = DEM
        self.DIR = DIR
        self.ncells = cu.stream_find(lat,lon,self.DEM,
            self.DIR,cu.ncols,cu.nrows)
        self.structure = cu.stream_cut(self.ncells)
    #------------------------------------------------------
    # Guardado shp de cauce
    #------------------------------------------------------
    def Save_Stream2Map(self,ruta,DriverFormat='ESRI Shapefile',
        EPSG=4326):
        'Descripcion: Guarda el cauce trazado en un mapa \n'\
        '\n'\
        'Parametros\n'\
        '----------\n'\
        'self : Inicia las variables vacias.\n'\
        'ruta : Nombre del lugar donde se va a guardar el cauce.\n'\
        'DriverFormat : Tipo de mapa vectorial.\n'\
        'EPSG : Codigo de tipo de proyeccion usada (defecto 4326).\n'\
        '\n'\
        'Retornos\n'\
        '----------\n'\
        'Guarda el cauce en el formato especificado.\n'\
        #Escribe el shp del cauce
        if ruta.endswith('.shp')==False:
            ruta=ruta+'.shp'
        spatialReference = osgeo.osr.SpatialReference()
        spatialReference.ImportFromEPSG(EPSG)
        driver = osgeo.ogr.GetDriverByName(DriverFormat)
        if os.path.exists(ruta):
             driver.DeleteDataSource(ruta)
        shapeData = driver.CreateDataSource(ruta)
        layer = shapeData.CreateLayer('layer1',
            spatialReference, osgeo.ogr.wkbLineString)
        layerDefinition = layer.GetLayerDefn()
        line = osgeo.ogr.Geometry(osgeo.ogr.wkbLineString)
        for x,y in zip(self.structure[0],self.structure[1]):
            line.AddPoint_2D(float(x),float(y))
        feature = osgeo.ogr.Feature(layerDefinition)
        feature.SetGeometry(line)
        feature.SetFID(0)
        layer.CreateFeature(feature)
        line.Destroy()
        feature.Destroy()
        shapeData.Destroy()
    #------------------------------------------------------
    # Plot de variables
    #------------------------------------------------------
    def Plot_Profile(self,ruta=None):
        'Descripcion: Grafica el perfil del cauce trazado \n'\
        '\n'\
        'Parametros\n'\
        '----------\n'\
        'self : Inicia las variables vacias.\n'\
        'ruta : Nombre de la imagen si se va a guardar la imagen del cauce.\n'\
        '\n'\
        'Retornos\n'\
        '----------\n'\
        'Guarda el cauce en el formato especificado.\n'\
        #Escribe el shp del cauce
        fig=pl.figure(facecolor='w',edgecolor='w')
        ax=fig.add_subplot(111)
        ax.plot(self.structure[3],self.structure[2],lw=2)
        ax.set_xlabel('Distancia $[mts]$',size=14)
        ax.set_ylabel('Elevacion $[m.s.n.m]$',size=14)
        ax.grid(True)
        if ruta is not None:
            pl.savefig(ruta,bbox_inches='tight')
        pl.show()
    #def Plot_Map(self,ruta=None
