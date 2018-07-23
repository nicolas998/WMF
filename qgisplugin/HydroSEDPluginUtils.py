import os.path

from qgis.core import QgsRasterLayer, QgsMapLayerRegistry, QgsVectorLayer, QgsFillSymbolV2
from PyQt4 import QtGui, uic
import netCDF4
from wmf import wmf
import numpy as np 
import pandas as pd
import osgeo.ogr as ogr

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
        self.Nc2Save = []
        self.Interpol_Columnas = []

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
            self.cuenca = wmf.SimuBasin(x, y, self.DEM, self.DIR, stream=self.stream,  umbral = umbral)
        else:
            self.cuenca = wmf.SimuBasin(x, y, self.DEM, self.DIR,  umbral = umbral)
        # Guarda los shapes de divisoria y de red hidrica.
        self.cuenca.Save_Basin2Map(PathDiv, dxp)
        self.cuenca.Save_Net2Map(PathRed, dxp, umbral)
        # Guarda el nc de la cuenca 
        if len(PathNC)>2:
            self.cuenca.Save_SimuBasin(PathNC)
    
    def BasinNc2Network(self,pathNetwork, names):
        '''Guarda una red hidrica con las variables seleccionadas del diccionario NC'''
        # genera el diccionario de las variables a guardar
        DicVar = {}
        for n in names:
            DicVar.update({n.encode():self.DicBasinNc[n]['var'].data})
        print DicVar
        #Guarda la red hidrica 
        self.cuenca.Save_Net2Map(pathNetwork, wmf.cu.dxp, 
            umbral = self.cuenca.umbral,
            EPSG = int(self.cuenca.epsg),
            Dict = DicVar)
        
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

    def Basin_Update(self, PathNC):
        '''Actualiza el archivo nc de la cuenca con las variables agregadas o borradas de la misma'''
        #Lectura del archivo 
        g = netCDF4.Dataset(PathNC,'a')
        #Inclusion de nuevas variables
        for l in self.Nc2Save:
            #Selecciona el grupo del nc en donde va a meter la variable
            Grupo = self.DicBasinNc[l]['categoria']
            Group = g.groups[Grupo]
            #mira si el grupo tiene la dimension ncells en caso de que no, la crea
            try:
                pos = Group.dimensions.keys().index('ncell')
            except:
                DimNcell = Group.createDimension('ncell',self.cuenca.ncells)
            #Obtiene el nombre, tipo y variable a actualizar
            nombre = l
            tipo = self.DicBasinNc[l]['tipo']
            if tipo == 'float32':
                tipo = 'f4'
            elif tipo == 'int64':
                tipo = 'i4'
            Var = np.copy(self.DicBasinNc[l]['var'])
            #Trata de meter la variable como algo no existente 
            try:
                print 'variable nueva'
                VarName = Group.createVariable(nombre,tipo,('ncell',),zlib=True)
            #si ya existe la variable la sobre escribe 
            except:
                print 'variable vieja'
                VarName = Group.variables[nombre]
            #guarda la variable
            VarName[:] = Var
        self.Nc2Save = []
        #Cerrado del archivo nc
        g.close()

    def Basin_LoadBasin(self, PathNC):
        # Numero Total de Variables
        self.NumDicBasinNcVariables = 0
        # Numero Total de Variables Basicas
        self.NumDicBasinNcVariablesBasicas = 0
        #Cargar la cuenca y sus variables base a WMF 
        self.cuenca = wmf.SimuBasin(rute = PathNC)
        #Cargar las variables de la cuenca a un diccionario.
        g = netCDF4.Dataset(PathNC)
        for grupoKey in ['base','Geomorfo','Hidro']:        
            #Carga los grupos de variables en donde si se tengan variables
            if len(g.groups[grupoKey].variables.keys())>0:
                #itera
                for k in g.groups[grupoKey].variables.keys():
                    #Evalua si tiene la misma cantidad de celdas y puede ser un mapa
                    shape = g.groups[grupoKey].variables[k].shape
                    MapaRaster = False
                    for s in shape:
                        if s == self.cuenca.ncells:
                            MapaRaster = True
                        #Actualiza el diccionario
                        self.DicBasinNc.update({k:
                            {'nombre':k,
                            'tipo':g.groups[grupoKey].variables[k].dtype.name,
                            'shape':g.groups[grupoKey].variables[k].shape,
                            'raster':MapaRaster,
                            'basica': True,
                            'categoria': grupoKey,
                            'var': g.groups[grupoKey].variables[k][:]}})
                        self.NumDicBasinNcVariables = self.NumDicBasinNcVariables + 1
                        self.NumDicBasinNcVariablesBasicas = self.NumDicBasinNcVariablesBasicas + 1
        g.close()
        #Cargar la cuenca y sus variables base a WMF 
        self.cuenca = wmf.SimuBasin(rute = PathNC)
        self.NumDicBasinWMFVariables = 50
        self.cuenca.GetGeo_Cell_Basics()
        #Area de la cuenca y codigo EPSG  
        return self.cuenca.ncells*wmf.cu.dxp**2./1e6, self.cuenca.epsg, wmf.models.dxp, wmf.cu.nodata, self.cuenca.umbral
    
    def Basin_LoadBasinDivisory(self, PathDivisory):
        # Guarda los shapes de divisoria y de red hidrica.
        self.cuenca.Save_Basin2Map(PathDivisory, wmf.cu.dxp)
    
    def Basin_LoadBasinNetwork(self, PathNetwork):
        # Guarda los shapes de divisoria y de red hidrica.
        self.cuenca.Save_Net2Map(PathNetwork, wmf.cu.dxp, self.cuenca.umbral)
    
    def Basin_LoadVariableFromDicNC(self, VarName):
        '''Toma una variable del diccionario de variables basicas y la carga a Qgis'''
        #Transforma a un raster 
        rutaSalida = '/tmp/HydroSED/'+VarName+'.tiff'
        self.cuenca.Transform_Basin2Map(self.DicBasinNc[VarName]['var'],
            ruta = rutaSalida,
            EPSG = self.cuenca.epsg)
        return rutaSalida
    
    def Basin_LoadVariableFromDicWMF(self,VarName):
        '''Toma una variable del diccionario de variables basicas y la carga a Qgis'''
        #Transforma a un raster 
        rutaSalida = '/tmp/HydroSED/'+VarName+'.tiff'
        self.cuenca.Transform_Basin2Map(self.DicBasinWMF[VarName]['var'],
            ruta = rutaSalida,
            EPSG = self.cuenca.epsg)
        return rutaSalida
    
    def Basin_Raster2WMF(self, VarName, VarPath, VarGroup):
        '''toma una ruta y convierte un mapa raster a cuenca para luego ponerlo 
        en el diccionario de WMF'''
        #Lee la variable 
        Var, prop, epsg = wmf.read_map_raster(VarPath)
        if epsg == self.cuenca.epsg:
            #Convierte a cuenca
            Var = self.cuenca.Transform_Map2Basin(Var, prop)
            #Actualiza el diccionario de WMF
            self.DicBasinWMF.update({VarName:
                {'nombre':VarName,
                'tipo':Var.dtype.name,
                'shape':Var.shape,
                'raster':True,
                'basica': False,
                'categoria': VarGroup,
                'var': Var}})
            return 0
        else:
            return 1
    
    def Basin_GeoGetAcumSlope(self):
        #self.cuenca.GetGeo_Cell_Basics()
        self.DicBasinWMF.update({'Area':
            {'nombre':'Area',
            'tipo':'float32',
            'shape':self.cuenca.CellAcum.shape,
            'raster':True,
            'basica': False,
            'categoria': 'Geomorfo',
            'var': self.cuenca.CellAcum}})
        self.DicBasinWMF.update({'Pendiente':
            {'nombre':'Pendiente',
            'tipo':'float32',
            'shape':self.cuenca.CellSlope.shape,
            'raster':True,
            'basica': False,
            'categoria': 'Geomorfo',
            'var': self.cuenca.CellSlope}})
            
    def Basin_GeoGetOrder(self):
        self.cuenca.GetGeo_StreamOrder(umbral = self.cuenca.umbral)
        self.DicBasinWMF.update({'Order_hills':
            {'nombre':'Order_hills',
            'tipo':self.cuenca.CellHorton_Hill.dtype.name,
            'shape':self.cuenca.CellHorton_Hill.shape,
            'raster':True,
            'basica': False,
            'categoria': 'Geomorfo',
            'var': self.cuenca.CellHorton_Hill}})
        self.DicBasinWMF.update({'Order_channels':
            {'nombre':'Order_channels',
            'tipo':self.cuenca.CellHorton_Stream.dtype.name,
            'shape':self.cuenca.CellHorton_Stream.shape,
            'raster':True,
            'basica': False,
            'categoria': 'Geomorfo',
            'var': self.cuenca.CellHorton_Stream}})
    
    def Basin_GeoGetIT(self):
        IT = self.cuenca.GetGeo_IT()
        self.DicBasinWMF.update({'Topo_index':
            {'nombre':'Topo_index',
            'tipo':IT.dtype.name,
            'shape':IT.shape,
            'raster':True,
            'basica': False,
            'categoria': 'Geomorfo',
            'var': IT}})
    
    def Basin_GeoGetChannels(self):
        #self.cuenca.GetGeo_Cell_Basics()
        self.DicBasinWMF.update({'Channels':
            {'nombre':'Channels',
            'tipo':self.cuenca.CellCauce.dtype.name,
            'shape':self.cuenca.CellCauce.shape,
            'raster':True,
            'basica': False,
            'categoria': 'Geomorfo',
            'var': self.cuenca.CellCauce}})

    def Basin_GeoGetDist2Out(self):
        self.cuenca.GetGeo_WidthFunction(show = False)
        self.DicBasinWMF.update({'Dist2Out':
            {'nombre':'Dist2Out',
            'tipo':self.cuenca.CellDist2Out.dtype.name,
            'shape':self.cuenca.CellDist2Out.shape,
            'raster':True,
            'basica': False,
            'categoria': 'Geomorfo',
            'var': self.cuenca.CellDist2Out}})
    
    def Basin_GeoGetHAND(self):
        self.cuenca.GetGeo_HAND(umbral = self.cuenca.umbral)
        self.DicBasinWMF.update({'HAND':
            {'nombre':'HAND',
            'tipo':self.cuenca.CellHAND.dtype.name,
            'shape':self.cuenca.CellHAND.shape,
            'raster':True,
            'basica': False,
            'categoria': 'Geomorfo',
            'var': self.cuenca.CellHAND}})
        self.DicBasinWMF.update({'HDND':
            {'nombre':'HDND',
            'tipo':self.cuenca.CellHDND.dtype.name,
            'shape':self.cuenca.CellHDND.shape,
            'raster':True,
            'basica': False,
            'categoria': 'Geomorfo',
            'var': self.cuenca.CellHDND}})
        self.DicBasinWMF.update({'HAND_class':
            {'nombre':'HAND_class',
            'tipo':self.cuenca.CellHAND_class.dtype.name,
            'shape':self.cuenca.CellHAND_class.shape,
            'raster':True,
            'basica': False,
            'categoria': 'Geomorfo',
            'var': self.cuenca.CellHAND_class}})
            
    def Basin_GeoGetParameters(self):
        self.cuenca.GetGeo_Parameters()
        return self.cuenca.GeoParameters, self.cuenca.Tc

    def Interpol_GetFields(self, Path2Points):
        '''Entrega una lista con los nombres de los atributos de un shp de puntos'''
        #Lectura del archivo 
        Driver = ogr.Open(Path2Points)
        layer = Driver.GetLayer(0)
        #Obtiene los nombres de las columnas
        self.Interpol_Columnas = []
        ldefn = layer.GetLayerDefn()
        for n in range(ldefn.GetFieldCount()):
            fdefn = ldefn.GetFieldDefn(n)
            self.Interpol_Columnas.append(fdefn.name)
        #Cierra el archivos
        Driver.Destroy()
        
    def Interpol_GetDateTimeParams(self, Path2Excel):
        '''Encuentra las fechas y el paso de tiempo del archivo de excel que contiene los registros de lluvia'''
        #Abre el archivo y encuentra fechas
        Data = pd.read_excel(Path2Excel)
        self.InterpolData = Data.copy()
        self.Interpol_fi = Data.index[0].to_pydatetime()
        self.Interpol_ff = Data.index[-1].to_pydatetime()
        self.Interpol_fd = Data.index[1] - Data.index[0]
    
    def Interpol_GetInterpolation(self, Path2Shp,Campo2Read,fi,ff,fd,expo, PathOutput):
        '''Interpola los campos de precipitacion para el periodo seleccionado con las estaciones disponibles'''
        #Fechas en letras
        fi = fi.strftime('%Y-%m-%d-%H:%M')
        ff = ff.strftime('%Y-%m-%d-%H:%M')
        tr = '%dS'%fd
        #Obtiene los ids de excel 
        idExcel = self.InterpolData.columns.values.tolist()
        idExcel = [i.encode() for i in idExcel]
        #Lee el shp con los puntos y los ids 
        xy,idShape = wmf.read_map_points(Path2Shp,[Campo2Read.encode()])
        idShape = idShape[Campo2Read].astype(int).astype(str).tolist()
        #Organiza los datos para interpolar
        xyNew = []
        ColGood = []
        for i in idShape:
            try:
                pos = idExcel.index(i)
                xyNew.append(xy.T[pos].tolist())
                ColGood.append(i)
            except:
                pass
        xyNew = np.array(xyNew)
        Data = self.InterpolData[ColGood][fi:ff].resample(tr).sum()
        #Interpola
        self.cuenca.rain_interpolate_idw(xyNew.T, Data, PathOutput,p = expo)
        
        
        
        
