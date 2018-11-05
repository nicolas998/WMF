from __future__ import print_function
from __future__ import absolute_import
from builtins import zip
from builtins import str
from builtins import range
from builtins import object
import os.path

from qgis.core import QgsRasterLayer,  QgsProject, QgsVectorLayer, QgsFillSymbol,QgsRenderContext
from qgis.PyQt import QtGui, uic
import netCDF4
from wmf import wmf
import numpy as np 
import datetime as dt
import scipy.stats as stat
import pandas as pd
import osgeo.ogr as ogr
from . import HydroSEDPlots as HSplots

class controlHS(object):
    
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
        self.DicParameters = {}
        self.Nc2Save = []
        self.Interpol_Columnas = []

    def cargar_mapa_raster (self,pathMapaRaster):
    
        retornoCargaLayerMapaRaster = False
    
        pathMapaRaster = pathMapaRaster.strip ()
    
        if (os.path.exists (pathMapaRaster)):
    
            baseNameMapaRaster   = os.path.basename (pathMapaRaster)
            baseNameMapaRaster = os.path.splitext(baseNameMapaRaster)[0]
            layerMapaRaster = QgsRasterLayer (pathMapaRaster, baseNameMapaRaster)
            QgsProject.instance ().addMapLayer (layerMapaRaster)
    
            retornoCargaLayerMapaRaster = layerMapaRaster.isValid()
    
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
            QgsProject.instance().addMapLayer(layerMapaVector)

            if tipo_style == self.TIPO_STYLE_POLILINEA:

                symbols = layerMapaVector.renderer().symbols(QgsRenderContext())
                symbol = symbols[0]
                symbol.setColor(QtGui.QColor.fromRgb(color[0],color[1],color[2]))
                symbol.setWidth(width)

            #try:
            #    symbol.setWidth(width)
            #except:
            #    symbol.setBorderWidth(width)
            #if layerMapVector.geometryType() == QGis.Polygon:

            elif tipo_style == self.TIPO_STYLE_POLIGONO:

                Render = layerMapaVector.renderer()
                mySymbol1 = QgsFillSymbol.createSimple({'color':'blue', 
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
    
    def __pathAsynch2paths__(self, path):
        filename, ext = os.path.splitext(path)
        path_rvr = filename + '.rvr'
        path_prm = filename + '.prm'
        path_lookup = filename + '.lookup'
        return path_rvr, path_prm, path_lookup
        
    def trazador_cuenca(self,x,y, dxp,umbral,PathDiv, PathRed, PathNC,
        PathAsynch,PathDEM, PathDIR,TopoNodes = False, LastStream = True):
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
        if len(PathAsynch)>2:
            path_rvr, path_prm, path_lookup = self.__pathAsynch2paths__(PathAsynch)
            self.cuenca.Transform_Basin2Asynch(path_rvr, 
                path_prm,
                path_lookup)
    
    def BasinNc2Network(self,pathNetwork, names):
        '''Guarda una red hidrica con las variables seleccionadas del diccionario NC'''
        # genera el diccionario de las variables a guardar
        DicVar = {}
        for n in names:
            DicVar.update({n:self.DicBasinNc[n]['var'].data})
        # fix_print_with_import
        print(DicVar)
        #Guarda la red hidrica 
        self.cuenca.Save_Net2Map(pathNetwork, wmf.cu.dxp, 
            umbral = self.cuenca.umbral,
            EPSG = int(self.cuenca.epsg),
            Dict = DicVar)
    
    def BasinWMF2Network(self,pathNetwork, names):
        '''Guarda una red hidrica con las variables seleccionadas del diccionario NC'''
        # genera el diccionario de las variables a guardar
        DicVar = {}
        for n in names:
            DicVar.update({n.encode():self.DicBasinWMF[n]['var']})
        # fix_print_with_import
        print(DicVar)
        #Guarda la red hidrica 
        self.cuenca.Save_Net2Map(pathNetwork, wmf.cu.dxp, 
            umbral = self.cuenca.umbral,
            EPSG = int(self.cuenca.epsg),
            Dict = DicVar)
        
    def hidologia_balance(self, dxp, umbral, PathRain, PathETR):
        #Se fija si la lluvia es un path o un valor 
        try:
            Rain = float(PathRain)
        except:
            #Trata de sacarlo del WMF
            try:
                Rain = np.copy(self.DicBasinWMF[PathRain]['var'])
            except:
                try:
                    Rain = np.copy(self.DicBasinNc[PathRain]['var'])
                except:
                    return 1
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
            'var': self.cuenca.CellQmed,
            'saved':False}})
        self.DicBasinWMF.update({'ETR':
            {'nombre':'ETR',
            'tipo':'float32',
            'shape':self.cuenca.CellETR.shape,
            'raster':True,
            'basica': False,
            'categoria': 'Hidro',
            'var': self.cuenca.CellETR,
            'saved':False}})
        Runoff = Rain - self.cuenca.CellETR
        self.DicBasinWMF.update({'Runoff':
            {'nombre':'Runoff',
            'tipo':'float32',
            'shape':self.cuenca.CellETR.shape,
            'raster':True,
            'basica': False,
            'categoria': 'Hidro',
            'var': Runoff,
            'saved':False}})
        #Retorna el resultado a la salida 
        return 0,self.cuenca.CellQmed[-1]

    def hidrologia_extremos_regional(self, Qmed, CoefExpoList, Pdf2Use, MaxorMin):
        '''Obtiene caudales maximos o minimos para los periodos de retorno 
        de 2.33, 5, 10, 25, 50, 100, 500'''
        #Organiza cioeficientes y expo 
        Coef = [CoefExpoList[i] for i in [0, 2]]
        Expo = [CoefExpoList[i] for i in [1, 3]]
        #Trata de hacerlo
        #try:
        # fix_print_with_import
        print(Coef)
        # fix_print_with_import
        print(Expo)
        # fix_print_with_import
        print(Pdf2Use)
        #Hace de acuerdo a una cosa o la otra 
        if MaxorMin == 'QMax':
            Qext = self.cuenca.GetQ_Max(Qmed, Coef, Expo, Dist = Pdf2Use)
        elif MaxorMin == 'QMin':
            Qext = self.cuenca.GetQ_Min(Qmed, Coef, Expo, Dist = Pdf2Use)
        #Actualiza el diccionario 
        Nombre = MaxorMin
        #Actualiza diccionarios 
        for Tr,Q in zip([2.33, 5, 10, 25, 50, 100], Qext):
            nombre2 = Nombre + '_' + str(Tr)
            self.DicBasinWMF.update({nombre2:{'nombre':nombre2,
            'tipo':'float32',
            'shape':Q.shape,
            'raster':True,
            'basica': False,
            'categoria': 'Hidro',
            'var': np.copy(Q),
            'saved':False}})
          #  return 0
        #except:
         #   return 1

    def Basin_Update(self, PathNC):
        '''Actualiza el archivo nc de la cuenca con las variables agregadas o borradas de la misma'''
        #Lectura del archivo 
        g = netCDF4.Dataset(PathNC,'a')
        #Inclusion de nuevas variables
        for l in list(self.DicBasinNc.keys()):#self.Nc2Save:
            #Si la variable no esta guradada la actualiza en el nc
            if self.DicBasinNc[l]['saved'] is False:
                #Selecciona el grupo del nc en donde va a meter la variable
                Grupo = self.DicBasinNc[l]['categoria']
                Group = g.groups[Grupo]
                #mira si el grupo tiene la dimension ncells en caso de que no, la crea
                try:
                    pos = list(Group.dimensions.keys()).index('ncell')
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
                    #print 'variable nueva'
                    VarName = Group.createVariable(nombre,tipo,('ncell',),zlib=True)
                #si ya existe la variable la sobre escribe 
                except:
                    #print 'variable vieja'
                    VarName = Group.variables[nombre]
                #guarda la variable y la actualiza en estado a guardada
                VarName[:] = Var
                self.DicBasinNc[l]['saved'] = True
        #self.Nc2Save = []
        #Cerrado del archivo nc
        g.close()

    def Basin_LoadBasin(self, PathNC, LoadSim = False, LoadSed = False):
        '''Carga un proyecto cuenca de nc en la memoria de Qgis:
        LoadSim: Carga las variables de simulacion
        LoadSed: Carga variables de simulacion de sedimentos'''
        #Inicia de cero el diccionario de la cuenca 
        self.DicBasinNc = {}
        self.DicBasinWMF = {}
        # Numero Total de Variables
        self.NumDicBasinNcVariables = 0
        # Numero Total de Variables Basicas
        self.NumDicBasinNcVariablesBasicas = 0
        #Cargar la cuenca y sus variables base a WMF 
        self.cuenca = wmf.SimuBasin(rute = PathNC)
        #Cargar las variables de la cuenca a un diccionario.
        g = netCDF4.Dataset(PathNC)
        Basicas = [True, True, False]
        ListGrupos = ['base','Geomorfo','Hidro']
        if LoadSim:
            ListGrupos.append('SimHidro')
            Basicas.append(True)
        if LoadSed:
            ListGrupos.append('SimSediments')
            Basicas.append(True)
        for grupoKey,basicas in zip(ListGrupos, Basicas):  
            #Carga los grupos de variables en donde si se tengan variables
            if len(list(g.groups[grupoKey].variables.keys()))>0:
                # fix_print_with_import
                print(list(g.groups[grupoKey].variables.keys()))
                #itera
                for k in list(g.groups[grupoKey].variables.keys()):
                    #Evalua si tiene la misma cantidad de celdas y puede ser un mapa
                    shape = g.groups[grupoKey].variables[k].shape
                    MapaRaster = False
                    for s in shape:
                        if s == self.cuenca.ncells: MapaRaster = True
                    #Actualiza el diccionario
                    self.DicBasinNc.update({k:
                        {'nombre':k,
                        'tipo':g.groups[grupoKey].variables[k].dtype.name,
                        'shape':g.groups[grupoKey].variables[k].shape,
                        'raster':MapaRaster,
                        'basica': basicas,
                        'categoria': grupoKey,
                        'var': g.groups[grupoKey].variables[k][:],
                        'saved':True}})
                    self.NumDicBasinNcVariables = self.NumDicBasinNcVariables + 1
                    self.NumDicBasinNcVariablesBasicas = self.NumDicBasinNcVariablesBasicas + 1
        #Carga parametros de calib del modelo 
        if LoadSim:
            #Si tiene param
            if len(list(g.groups['Parametros'].variables.keys()))>0:
                for k in list(g.groups['Parametros'].variables.keys()):
                    self.DicParameters.update({k:{
                        'nombre': k,
                        'var': g.groups['Parametros'].variables[k][:]}})
        # fix_print_with_import
        print(self.DicParameters)
        #Cierra el archivo netCDF de la cuenca        
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
    
    def Basin_LoadVariableFromDicNC(self, VarName, capa = None):
        '''Toma una variable del diccionario de variables basicas y la carga a Qgis'''
        #Variable de acuerdo a si hay capa o no
        if capa is None:
            Variable = np.copy(self.DicBasinNc[VarName]['var'])
            rutaSalida = '/tmp/HydroSED/raster/'+VarName+'.tiff'
        else:
            Variable = np.copy(self.DicBasinNc[VarName]['var'][capa])
            rutaSalida = '/tmp/HydroSED/raster/'+VarName+'_'+str(capa+1)+'.tiff'
        # fix_print_with_import
        print(Variable)
        #Conversion
        self.cuenca.Transform_Basin2Map(Variable,
            ruta = rutaSalida,
            EPSG = self.cuenca.epsg)
        #Devuelve la ruta temporal donde esta el mapa raster con al variable
        return rutaSalida
    
    def Basin_LoadVariableFromDicWMF(self,VarName):
        '''Toma una variable del diccionario de variables basicas y la carga a Qgis'''
        #Transforma a un raster 
        rutaSalida = '/tmp/HydroSED/raster/'+VarName+'.tiff'
        self.cuenca.Transform_Basin2Map(self.DicBasinWMF[VarName]['var'],
            ruta = rutaSalida,
            EPSG = self.cuenca.epsg)
        return rutaSalida
    
    def Basin_Raster2WMF(self, VarName, VarPath, VarGroup, PathOrValue = 'Path', Value = None):
        '''toma una ruta y convierte un mapa raster a cuenca para luego ponerlo 
        en el diccionario de WMF'''
        Convierte = False
        #si es un mapa raster
        if PathOrValue == 'Path':
            Var, prop, epsg = wmf.read_map_raster(VarPath)
            if epsg == self.cuenca.epsg:
                #Convierte a cuenca
                Var = self.cuenca.Transform_Map2Basin(Var, prop)
                Convierte = True
            else:
                return 1
        #Si es una constante
        else:
            Var = np.ones(self.cuenca.ncells)*Value
            Convierte = True
        #Actualiza el diccionario de WMF
        if Convierte:
            self.DicBasinWMF.update({VarName:
                {'nombre':VarName,
                'tipo':Var.dtype.name,
                'shape':Var.shape,
                'raster':True,
                'basica': False,
                'categoria': VarGroup,
                'var': Var,
                'saved':False}})
            return 0
        else:
            return 1
    
    def Basin_Geo2Excel(self,ExcelPath):
        '''Guarda los parametros geomorfologicos de la cuenca en un archivo de excel'''
        #Obtiene param de curvas
        self.cuenca.GetGeo_Ppal_Hipsometric()
        #PArametros
        GeoData = pd.DataFrame.from_dict(self.cuenca.GeoParameters, 'index')
        GeoData = pd.DataFrame(GeoData.values, index=GeoData.index,columns=['Parametros'])
        TcData = pd.DataFrame.from_dict(self.cuenca.Tc, 'index')
        TcData = pd.DataFrame(TcData.values, TcData.index, columns=['T viaje [hrs]'])
        #perfil
        Perfil = pd.DataFrame(np.vstack([self.cuenca.ppal_stream[1]/1000., self.cuenca.ppal_stream[0][::-1]]).T,
            index=list(range(self.cuenca.ppal_stream.shape[1])), 
            columns=['Dist2Out[km]','Elevacion[m]'])
        #Curva hipsometrica
        CurvaHipso = pd.DataFrame(np.vstack([self.cuenca.hipso_basin[0]*wmf.cu.dxp**2/1e6, self.cuenca.hipso_basin[1]]).T,
            index = list(range(self.cuenca.hipso_basin.shape[1])),
            columns = ['Area[km2]','Elevacion[m]'])
        #Escritor
        W = pd.ExcelWriter(ExcelPath)
        #Escribe 
        GeoData.to_excel(W)
        TcData.to_excel(W, startcol=3, startrow=0)
        Perfil.to_excel(W, startcol=6, startrow=0)
        CurvaHipso.to_excel(W, startcol=10, startrow=0)
        #Cierra el archivo 
        W.close()
        return 0
        
    
    def Basin_GeoGetAcumSlope(self):
        #self.cuenca.GetGeo_Cell_Basics()
        self.DicBasinWMF.update({'Area':
            {'nombre':'Area',
            'tipo':'float32',
            'shape':self.cuenca.CellAcum.shape,
            'raster':True,
            'basica': False,
            'categoria': 'Geomorfo',
            'var': self.cuenca.CellAcum,
            'saved':False}})
        self.DicBasinWMF.update({'Pendiente':
            {'nombre':'Pendiente',
            'tipo':'float32',
            'shape':self.cuenca.CellSlope.shape,
            'raster':True,
            'basica': False,
            'categoria': 'Geomorfo',
            'var': self.cuenca.CellSlope,
            'saved':False}})
            
    def Basin_GeoGetOrder(self):
        self.cuenca.GetGeo_StreamOrder(umbral = self.cuenca.umbral)
        self.DicBasinWMF.update({'Order_hills':
            {'nombre':'Order_hills',
            'tipo':self.cuenca.CellHorton_Hill.dtype.name,
            'shape':self.cuenca.CellHorton_Hill.shape,
            'raster':True,
            'basica': False,
            'categoria': 'Geomorfo',
            'var': self.cuenca.CellHorton_Hill,
            'saved':False}})
        self.DicBasinWMF.update({'Order_channels':
            {'nombre':'Order_channels',
            'tipo':self.cuenca.CellHorton_Stream.dtype.name,
            'shape':self.cuenca.CellHorton_Stream.shape,
            'raster':True,
            'basica': False,
            'categoria': 'Geomorfo',
            'var': self.cuenca.CellHorton_Stream,
            'saved':False}})
    
    def Basin_GeoGetIT(self):
        IT = self.cuenca.GetGeo_IT()
        self.DicBasinWMF.update({'Topo_index':
            {'nombre':'Topo_index',
            'tipo':IT.dtype.name,
            'shape':IT.shape,
            'raster':True,
            'basica': False,
            'categoria': 'Geomorfo',
            'var': IT,
            'saved':False}})
    
    def Basin_GeoGetChannels(self):
        #self.cuenca.GetGeo_Cell_Basics()
        self.DicBasinWMF.update({'Channels':
            {'nombre':'Channels',
            'tipo':self.cuenca.CellCauce.dtype.name,
            'shape':self.cuenca.CellCauce.shape,
            'raster':True,
            'basica': False,
            'categoria': 'Geomorfo',
            'var': self.cuenca.CellCauce,
            'saved':False}})

    def Basin_GeoGetDist2Out(self):
        self.cuenca.GetGeo_WidthFunction(show = False)
        self.DicBasinWMF.update({'Dist2Out':
            {'nombre':'Dist2Out',
            'tipo':self.cuenca.CellDist2Out.dtype.name,
            'shape':self.cuenca.CellDist2Out.shape,
            'raster':True,
            'basica': False,
            'categoria': 'Geomorfo',
            'var': self.cuenca.CellDist2Out,
            'saved':False}})
    
    def Basin_GeoGetHAND(self):
        self.cuenca.GetGeo_HAND(umbral = self.cuenca.umbral)
        self.DicBasinWMF.update({'HAND':
            {'nombre':'HAND',
            'tipo':self.cuenca.CellHAND.dtype.name,
            'shape':self.cuenca.CellHAND.shape,
            'raster':True,
            'basica': False,
            'categoria': 'Geomorfo',
            'var': self.cuenca.CellHAND,
            'saved':False}})
        self.DicBasinWMF.update({'HDND':
            {'nombre':'HDND',
            'tipo':self.cuenca.CellHDND.dtype.name,
            'shape':self.cuenca.CellHDND.shape,
            'raster':True,
            'basica': False,
            'categoria': 'Geomorfo',
            'var': self.cuenca.CellHDND,
            'saved':False}})
        self.DicBasinWMF.update({'HAND_class':
            {'nombre':'HAND_class',
            'tipo':self.cuenca.CellHAND_class.dtype.name,
            'shape':self.cuenca.CellHAND_class.shape,
            'raster':True,
            'basica': False,
            'categoria': 'Geomorfo',
            'var': self.cuenca.CellHAND_class,
            'saved':False}})

    def Basin_GeoGetOCG(self):
        '''Obtiene el coeficiente y el exponente de OCG (velez, 2001) para la cuenca'''
        #Obtiene param basicos de la cuenca 
        self.cuenca.GetGeo_Cell_Basics()
        Area = self.cuenca.CellAcum*wmf.cu.dxp**2/1e6
        #Obtiene el coeficiente 
        Coef, Expo = wmf.OCG_param(pend=self.cuenca.CellSlope, area=Area)
        Expo = np.ones(self.cuenca.ncells)*Expo
        #Agrea los resultados al diccionario
        self.DicBasinWMF.update({'OCG_coef':
            {'nombre':'OCG_coef',
            'tipo':Coef.dtype.name,
            'shape':Coef.shape,
            'raster':True,
            'basica': False,
            'categoria': 'Geomorfo',
            'var': np.copy(Coef),
            'saved':False}})
        self.DicBasinWMF.update({'OCG_exp':
            {'nombre':'OCG_exp',
            'tipo':Expo.dtype.name,
            'shape':Expo.shape,
            'raster':True,
            'basica': False,
            'categoria': 'Geomorfo',
            'var': np.copy(Expo),
            'saved':False}})
    
    def Basin_GeoGetKubota(self):
        '''Obtiene el coeficiente de kubota y sivapalan para el flujo de agua en el suelo
        para el correcto calculo el susuario debe tener definidas las variables:
            - Hu: almacenamiento capilar maximo [mm] DicBasinNc['h1_max']
            - ks: conductividad saturada del suelo [mm/s] DicBasinNc['v_coef'][1]'''
        #Busca si los parametros de la cuenca estan definidos, si no los reemplaza con constantes para Hu y ks 
        if u'h3_max' in list(self.DicBasinNc.keys()): 
            Hg = np.copy(self.DicBasinNc['h3_max']['var'])
            Hg[Hg<=0] = Hg.mean()
        else: 
            Hg = 250    
        if u'v_coef' in list(self.DicBasinNc.keys()):  
            ks = np.copy(self.DicBasinNc['v_coef']['var'][1])
        else: 
            ks = 0.00128 
        self.cuenca.GetGeo_Cell_Basics()
        So = np.copy(self.cuenca.CellSlope)
        Factor = (wmf.cu.dxp**2.)/1000./1000. #[m3/mm]
        #Calculo del coeficiente de kubota
        Coef = (ks*So*(wmf.cu.dxp**2.))/(3*(Hg*Factor)**2.)
        Coef[Coef==0] = np.percentile(Coef, 50)
        #ksh=(Ks*cuSalgar.CellSlope*(12.7**2.0))/(3*(Hg*0.9/1000.0)**2)
        Coef[np.where(np.isinf(Coef))]=np.mean(Coef[np.where(np.isfinite(Coef))])
        # fix_print_with_import
        print(Coef)
        #Actualiza el diccionario 
        self.DicBasinWMF.update({'Kubota_coef':
            {'nombre':'Kubota_coef',
            'tipo':Coef.dtype.name,
            'shape':Coef.shape,
            'raster':True,
            'basica': False,
            'categoria': 'Geomorfo',
            'var': np.copy(Coef),
            'saved': False}})
        KubotaExpo = np.ones(self.cuenca.ncells)*2
        self.DicBasinWMF.update({'Kubota_exp':
            {'nombre':'Kubota_exp',
            'tipo':KubotaExpo.dtype.name,
            'shape':KubotaExpo.shape,
            'raster':True,
            'basica': False,
            'categoria': 'Geomorfo',
            'var': np.copy(KubotaExpo),
            'saved': False}})
        
    def Basin_GeoGetRunoff(self, e1, Epsilon):
        '''Obtiene el coeficiente de escorrentia para carcavas'''
        #Obtiene pendiente y param requeridos
        self.cuenca.GetGeo_Cell_Basics()
        So = np.copy(self.cuenca.CellSlope)
        #Busca si la variable esta guardada en los diccionarios, si no, la reemplaza con una constante. 
        if 'Manning' in list(self.DicBasinWMF.keys()) or 'Manning' in list(self.DicBasinNc.keys()):
            try: 
                Man = np.copy(self.DicBasinWMF['Manning']['var'])
            except: 
                Man = np.copy(self.DicBasinNc['Manning']['var'])
            Tipo_var = Man.dtype.name
        else: 
            Man = np.ones([self.cuenca.ncells])*0.05 
            Tipo_var = 'float'              
        #Calcula 
        Coef = (float(Epsilon)/Man)*(So**2.)
        Coef[np.where(np.isinf(Coef))]=np.mean(Coef[np.where(np.isfinite(Coef))])
        Coef[Coef<=0] = np.percentile(Coef, 50)
        # fix_print_with_import
        print(Coef)
        Expo = (2./3.)*e1*np.ones([self.cuenca.ncells])
        #Pone en los diccionarios 
        self.DicBasinWMF.update({'Runoff_coef':
            {'nombre':'Runoff_coef',
            'tipo':Coef.dtype.name,
            'shape':Coef.shape,
            'raster':True,
            'basica': False,
            'categoria': 'Geomorfo',
            'var': np.copy(Coef),
            'saved':False}})
        self.DicBasinWMF.update({'Runoff_exp':
            {'nombre':'Runoff_exp',
            'tipo':Coef.dtype.name,
            'shape':Coef.shape,
            'raster':True,
            'basica': False,
            'categoria': 'Geomorfo',
            'var': np.copy(Expo),
            'saved':False}})
          
    def Basin_GeoGetParameters(self):
        self.cuenca.GetGeo_Parameters()
        return self.cuenca.GeoParameters, self.cuenca.Tc

    def Basin_GeoAcumVar(self, VarName, VarAcumName, WhatAcumDic,VarMaskName,WhatMaskDic,umbral):
        '''Acumula una variable teniendo en cuenta una mascara y un umbral sobre la misma'''
        #Selecciona variable de acumulacion 
        if WhatAcumDic == 'NC':
            Var = np.copy(self.DicBasinNc[VarAcumName]['var'])
        elif WhatAcumDic == 'WMF':
            Var = np.copy(self.DicBasinWMF[VarAcumName]['var'])
        #Selecciona la mascara
        Mask = np.ones(self.cuenca.ncells)
        if VarMaskName is not None:
            #Si hay mascara busca en que diccionario esta y la toma 
            if WhatMaskDic == 'NC':
                MaskTemp = np.copy(self.DicBasinNc[VarMaskName]['var'])
            elif WhatMaskDic == 'WMF':
                MaskTemp = np.copy(self.DicBasinWMF[VarMaskName]['var'])
            #Con el umbral convierte a la mascara en una variable de ceros y unos 
            Mask[MaskTemp<=umbral] = 0
        #Variable a acumular y acumulacion
        Var = Var*Mask
        AcumVar = wmf.cu.basin_acum_var(self.cuenca.structure[0], 
            np.ones((1,self.cuenca.ncells))*Var,
            self.cuenca.ncells)
        AcumVar = AcumVar[0]
        #Mete la variable en el diccionario de WMF 
        self.DicBasinWMF.update({VarName:
            {'nombre':VarName,
            'tipo':AcumVar.dtype.name,
            'shape':AcumVar.shape,
            'raster':True,
            'basica': False,
            'categoria': 'Geomorfo',
            'var': np.copy(AcumVar),
            'saved':False}})
        return 0

        
            

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
        #Lee el shp con los puntos y los ids 
        xy,idShape = wmf.read_map_points(Path2Shp,[Campo2Read])
        idShape = idShape[Campo2Read].astype(int).tolist()
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
        
    def Interpol_GetRainfallAcum(self, path2bin, inicio, fin):
        '''Obtiene un campo acumulado de precipitacion para el periodo especifico'''
        #Lee el binario y los datos en el intervalo
        Vsum = np.zeros(self.cuenca.ncells)
        # fix_print_with_import
        print(path2bin)
        for i in range(inicio, fin+1):
            vect,res = wmf.models.read_int_basin(path2bin,i,self.cuenca.ncells)
            if res == 0:
                vect = vect.astype(float)/1000.
                Vsum+=vect
        # fix_print_with_import
        print('salida')
        #Pasa el acumulado al diccionario de WMF 
        self.DicBasinWMF.update({'Lluvia':
            {'nombre':'Lluvia',
            'tipo':Vsum.dtype.name,
            'shape':Vsum.shape,
            'raster':True,
            'basica': False,
            'categoria': 'Hidro',
            'var': Vsum,
            'saved':False}})
            
    def Sim_GetQobsInfo(self,PathExcelQobs):
         DataQ = pd.read_excel(PathExcelQobs)
         self.QobsData = DataQ.copy()
         self.DataQindex = DataQ.index 
         self.QobsData_fi = DataQ.index[0].to_pydatetime()
         self.QobsData_ff = DataQ.index[-1].to_pydatetime()
         self.QobsData_fd = DataQ.index[1] - DataQ.index[0]
         idExcelQ = self.QobsData.columns.values.tolist()
         return idExcelQ,DataQ,self.DataQindex
                
           
    def Sim_SaveParameters(self, PathNC, ParamName, scalarParam):
        '''Actualiza el nc con un conjunto de parametros escalares nuevo'''
        #Abre el archivo nc y apunta al grupo de parametros
        g = netCDF4.Dataset(PathNC,'a')
        GrupoParam = g.groups['Parametros']
        Ejecuto = 1
        if len(ParamName)>0:
            #mira si el grupo ya tiene parametros adentro, si no crea la dimension
            try:
                pos = list(GrupoParam.dimensions.keys()).index('nparam')
            except:
                DimNparam = GrupoParam.createDimension('nparam',15)
                #Crea el sub-grupo con el nombre de los parametros dentro del nc 
                SubGroup = g.groups['Parametros'].createGroup('NombreParam')
                nnames = SubGroup.createDimension('nnames',15)
                nombres = SubGroup.createVariable('nombres', str, ('nnames',))
                nombres[:] = np.array(['Hu','Hg','Evp','Inf','Per','Loss','vRun','vSub','vSup','vChannel',
                    'SedParam','Exp1','Exp2','Exp3','Exp4'], dtype = object)
            #Crea la variable en el grupo o la sobre-escribe
            try:
                Var = GrupoParam.createVariable(ParamName,'f4',('nparam',),zlib=True)
            except:
                Var = GrupoParam.variables[ParamName]
            Var[:] = scalarParam
            self.DicParameters.update({ParamName:
                {'nombre': ParamName,
                'var':scalarParam}})
            Ejecuto = 0
            
        #cierra el archivo 
        g.close()
        return Ejecuto
    
    def ExpressionParser(self, expresion):
        '''Toma una expresion enviada en string, extrae las variables del diccionario
        y evalua la expresion'''
        #Busca en el diccionario de NC las variables de la expresion
        Lista = list(self.DicBasinWMF.keys())
        Lista.extend(list(self.DicBasinNc.keys()))
        for k in Lista:
            #Tamano y pos de la letra
            largo = len(k)
            pos = -1
            try:
                pos = expresion.index(k)
                #Extrae el valor de la var del diccionario y define variables
                try:
                    var = self.DicBasinNc[k]['var']
                except:
                    var = self.DicBasinWMF[k]['var']
                ejec = expresion[pos:pos+largo]+'=var'
                exec(ejec)
            except:
                pass
        #Formula la expresion 
        try:
            AA = eval(expresion)
        except:
            exec(expresion)
            pos = expresion.index('=')
            VarName = expresion[:pos].strip()
            AA = eval(VarName)
        return AA
   
    def BasinConvert2HillsOrChannels(self, Var, Metodo, Agregado):
        '''Agrega una variable por canales o laderas de acuerdo a una metodologia'''
        #Funcion de metodologias
        def __fx__(x, metodo = 'media'):
            if metodo == 'media':
                return np.nanmean(x)
            elif metodo[0] == 'P':
                percentil = int(metodo[1:])
                return np.nanpercentile(x, percentil)
            elif metodo == 'min':
                return np.nanmin(x)
            elif metodo == 'max':
                return np.nanmax(x)
            elif metodo == 'moda':
                res = stat.mode(x)
                return res.mode[0]
        #Variable vacia nula
        VarAgregada = np.ones(self.cuenca.ncells)*wmf.cu.nodata
        #Itera por las laderas del elemento cuenca
        for i in range(1,self.cuenca.nhills+1):
            #Define posiciones de acuerdo a la metodologia
            hacer = True
            pos = np.where(self.cuenca.hills_own == i)[0]
            if Agregado == 'Canales':
                pos2 = np.where((self.cuenca.hills_own == i) & (self.cuenca.CellCauce == 1))[0]
                if len(pos2)==0:
                    hacer = False
                else:
                    pos = pos2
            #Agrega la variable
            if hacer:
                VarTemporal = __fx__(Var[pos], Metodo)
                VarAgregada[pos] = VarTemporal
            #correccion no data de canales 
        #if Agregado == 'Canales':
            #posMalos = np.where((self.cuenca.CellCauce == 1) & (VarAgregada == wmf.cu.nodata))[0]
            #print posMalos
                ##posBuenos = np.where((self.cuenca.CellCauce == 1) & (VarAgregada <> wmf.cu.nodata))[0]
                #VarAgregada[posMalos] = VarAgregada[posBuenos].mean()
                #print VarAgregada[posMalos].min()
        return VarAgregada
    
    def Radar_FechasProcess(self, Fi, Ff, TimeStep, FechasRadar, ListaRadar):
        '''Obtiene las fechas adecuadas y las posiciones para ls ttfmacion de radar a cuenca'''
        #Rango de fechas a convertir
        self.ConvertDates = pd.date_range(Fi,Ff,freq=TimeStep)
        #TextoFechas = [i.strftime('%Y-%m-%d %H:%M') for i in self.ConvertDates.to_pydatetime()]
        #Obtiene fechas equivalentes para las fechas a convertir
        DatesRadar = pd.to_datetime(FechasRadar)
        DatesRadar = DatesRadar.ceil(TimeStep)
        # fix_print_with_import
        print(DatesRadar[:5])
        # fix_print_with_import
        print(DatesRadar[-5:])
        #Slice de analisis
        Flag = True
        pos = 0
        while Flag:
            try:
                Ini = np.where(DatesRadar == self.ConvertDates[pos])[0].tolist()[0]
                Flag = False
                # fix_print_with_import
                print(Ini, pos)
            except:
                pos += 1
            if pos>DatesRadar.size: Flag = False
        Flag = True
        pos = -1
        while Flag:
            try:
                Fin = np.where(DatesRadar == self.ConvertDates[pos])[0].tolist()[0]
                Flag = False
                # fix_print_with_import
                print(Fin, pos)
            except:
                pos -= 1
            if pos>DatesRadar.size: Flag = False
        #Corte en vectores
        self.DatesRadar = DatesRadar[Ini:Fin]
        FechasRadar = FechasRadar[Ini:Fin]
        self.ListaRadar = ListaRadar[Ini:Fin]
        #self.TextoFechas = TextoFechas
    
    def Radar_Conver2Basin(self, PathRadarBasin, TimeStep, umbral = 0.01,
        verbose = True, old = False):
        '''Convierte barridos de radar a la cuenca'''
        #Se fija si ya existia un binario con campos 
        if old:
            self.cuenca.rain_radar2basin_from_array(status='old',ruta_out= PathRadarBasin)
        #Convierte para las fechas
        Rain = np.zeros(self.cuenca.ncells)
        for ft,date in zip(self.ConvertDates, self.ConvertDates.to_pydatetime()):
            #Posiciones de datos de radar
            ListaPos = np.where(self.DatesRadar == ft)[0].tolist()
            Entra = True
            if Rain.mean()>umbral and len(ListaPos) == 0:
                RainFin = np.copy(Rain)
                Entra = False
            elif len(ListaPos)>0:
                Entra = True
            #Acumula para ese periodo
            if Entra:        
                Rain = np.zeros(self.cuenca.ncells)
                for l in ListaPos:
                    #Lee el netCDF y lo transforma 
                    g = netCDF4.Dataset(self.ListaRadar[l])
                    RadProp = [g.ncols, g.nrows, g.xll, g.yll, g.dx, g.dx]
                    Rain += self.cuenca.Transform_Map2Basin(g.variables['Rain'][:].T/ (12*1000.0),RadProp)
                    g.close()
                RainFin = np.copy(Rain)
            #Actualiza el binario con datos de radar
            dentro = self.cuenca.rain_radar2basin_from_array(vec = RainFin,
                ruta_out = PathRadarBasin,
                fecha = date-dt.timedelta(hours = 5),
                dt = TimeStep,
                umbral = umbral)
            if verbose:
                # fix_print_with_import
                print(date, RainFin.mean(), ListaPos)
        #Cierra el binario y cea el encabezado
        self.cuenca.rain_radar2basin_from_array(status = 'close',ruta_out = PathRadarBasin)
        self.cuenca.rain_radar2basin_from_array(status = 'reset')


    def Sim_Basin(self, Start, Nsteps, Calibracion, PathRain, DeltaT, exponentes, PathStore):
        '''Hace la simulacion hidrologica de la cuenca'''        
        #Set de calibracion
        wmf.models.sed_factor = Calibracion[-1]
        Calibracion = Calibracion[:-1]
        Calibracion = Calibracion[2:] + Calibracion[:2]
        #Set del intervalo de tiempo de simulacion
        wmf.models.dt = DeltaT
        #SEtea el tipo de velocidad de acuerdo a los exponentes
        for c,i in enumerate(exponentes):
            if i!=1:
                wmf.models.speed_type[c] = 2
            else:
                wmf.models.speed_type[c] = 1
        #Simulacion de la cuenca
        if wmf.models.sim_sediments == 0:
            Results, Qsim = self.cuenca.run_shia(Calibracion, 
               PathRain, 
               Nsteps,
               Start,
               ruta_storage = PathStore)
        elif wmf.models.sim_sediments == 1:
            Results, Qsim, Qsed = self.cuenca.run_shia(Calibracion, 
               PathRain, 
               Nsteps,
               Start,
               ruta_storage = PathStore)
        #Obtiene resultados como cosas genericas
        self.Sim_index = Qsim.index
        self.Sim_Streamflow = Qsim.copy()
        self.Sim_Rainfall = pd.Series(Results['Rain_hietogram'][0], index = self.Sim_index)
        self.Sim_Balance = pd.Series(Results['Balance'][0], index = self.Sim_index)
        self.Sim_Storage = Results['Storage']
        self.Sim_RainfallField = np.copy(Results['Rain_Acum'])
        #Resultados opcionales
        if wmf.models.show_storage == 1:
            self.Sim_StorageSerie = pd.DataFrame(Results['Mean_Storage'].T, index = self.Sim_index)
        if wmf.models.show_mean_speed == 1:
            self.Sim_SpeedSerie = pd.DataFrame(Results['Mean_Speed'].T, index = self.Sim_index)
        if wmf.models.retorno == 1:
            self.Sim_RetornoSerie = pd.DataFrame(wmf.models.mean_retorno, index = self.Sim_index)
            self.Sim_RetornoMap = np.copy(wmf.models.retorned)
        #Resultados de sedimentos 
        if wmf.models.sim_sediments == 1:
            #Mapas de erosion y depositacion
            self.Sim_ErosionMap = np.copy(wmf.models.volero)
            self.Sim_DepositionMap = np.copy(wmf.models.voldepo)
            #SErie de sedimentos simulada 
            self.Sim_Sediments = Qsed.copy()
        
    def Sim_setStates_ConstValue(self, Tanque, Valor):
        '''Establece condiciones de almacenamiento constantes basado en un valor constante para toda la cuenca'''
        #Si es un almacenamiento del modelo 
        if Tanque < 6:
            self.cuenca.set_Storage(Valor, Tanque - 1)
        #Si es sedimentos de la cuenca.
        else:
            wmf.models.vd = np.ones(self.cuenca.ncells)*Valor
    
    def Sim_setStates_FileValue(self, Tanque, Path2Bin, Fecha):
        '''Establece las condiciones de almacenamiento basado en un binario con datos de almacenamiento anteriores'''
        #Si es un tanque del modelo 
        if Tanque < 6:
            #Obtiene el record.
            Fecha = Fecha.strftime('%Y-%m-%d %H:%M')
            Data = wmf.read_storage_struct(Path2Bin)
            record = Data.index.get_loc(Fecha)
            #lee el archivo 
            Valor, res = wmf.models.read_float_basin_ncol(Path2Bin, 
                record, self.cuenca.ncells, 5) 
            if res == 0:
                self.cuenca.set_Storage(Valor[Tanque-1], Tanque - 1)
                return Valor[Tanque-1].mean(), Valor[Tanque-1]
            else:
                return 0.0
    
    def Sim_setStates_NcValue(self, Tanque):
        '''Establece condiciones de un tanque en funcion del estado que se tiene en la tabla NC'''
        #Saca el valor del diccionario
        if Tanque < 6:
            Valor = self.DicBasinNc['storage']['var'][Tanque-1]
            self.cuenca.set_Storage(Valor, Tanque -1)
            return Valor.mean(), Valor
        
    def Sim_Series2Excel(self,ExcelPath,dataframe):
        '''Guarda los dataframes a un excel.'''
        W = pd.ExcelWriter(ExcelPath)
        #Escribe 
        dataframe.to_excel(W)
        #Cierra el archivo 
        W.close()
        return 0     
        
    def Convert2CDC (self,serie):
        '''Toma una serie y arroja como salidas la serie organizada y el porcentaje de 
           excedencia de cada valor'''
        serie = np.sort(serie)
        porcen_s=[]
        for i in range(len(serie)):
            porcen_s.append((len(serie[serie>serie[i]]))/float(len(serie))*100)
        return np.array(serie),np.array(porcen_s)
        
    def Convert_Df2CDC(self,df):
        '''Toma un dataframe que tenga como index fechas, y como columnas id de estaciones
            arroja como salida un dataframe con la series organizadas y los respectivos
            porcentajes de excelencia de cada valor'''
        ids = np.array(df.keys())
        tupla = []
        Data = []
        for ID in ids:
            tupla.append((str(ID),'Q_sort'))
            tupla.append((str(ID),'P_exc'))
            Data.extend([self.Convert2CDC(df[ID].values)[0],self.Convert2CDC(df[ID].values)[1]])
        Data = np.array(Data)
        index = pd.MultiIndex.from_tuples(tupla, names=['Tramo','CDC'])
        dfCDC = pd.DataFrame(Data.T,index = np.arange(len(df.index)),columns = index)
        return dfCDC        
        
    def Convert_DfAnualCaudales(self,df):
        '''Toma un dataframe que tenga como index fechas, y como columnas id de estaciones
            arroja como salida un dataframe con las medias mensuales multianuales por cada tramo'''
        ids = np.array(df.keys())
        tupla = []
        Data = []
        media= df.groupby(df.index.month).mean()
        desv = df.groupby(df.index.month).std()
        
        for ID in ids:
            tupla.append((str(ID),'Media'))
            tupla.append((str(ID),'Desv_Estandar'))
            Data.extend([media[ID],desv[ID]])
        Data = np.array(Data)
        index = pd.MultiIndex.from_tuples(tupla, names=['Tramo','Media_Mensual'])
        dfAnual = pd.DataFrame(Data.T,index = np.arange(1,13),columns = index)
        return dfAnual
            
    def CDC_Series2Excel(self,ExcelPath,df_sim=None,df_obs=None):
        '''Guarda los dataframes a un excel.'''
        W = pd.ExcelWriter(ExcelPath,engine='xlsxwriter')
        #Escribe 
        if df_sim is not None:
            df_sim.to_excel(W,sheet_name='CDC_simulados')
        if df_obs is not None:
            df_obs.to_excel(W,sheet_name='CDC_observados')
        #Cierra el archivo 
        W.close()
        return 0 
        
    def MediaMensual_Q2Excel(self,ExcelPath,df_sim=None,df_obs=None):
        '''Guarda los dataframes a un excel.'''
        W = pd.ExcelWriter(ExcelPath,engine='xlsxwriter')
        #Escribe 
        if df_sim is not None:
            df_sim.to_excel(W,sheet_name='Media_Mensual_Qsim')
        if df_obs is not None:
            df_obs.to_excel(W,sheet_name='Media_Mensual_Qobs')
        #Cierra el archivo 
        W.close()
        return 0        
            


    
            
            
            
            
            
            
            
            
            
            
