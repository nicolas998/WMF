import os.path

from qgis.core import QgsRasterLayer, QgsMapLayerRegistry, QgsVectorLayer, QgsFillSymbolV2
from PyQt4 import QtGui, uic
import netCDF4
from wmf import wmf
import numpy as np 
import scipy.stats as stat
import pandas as pd
import osgeo.ogr as ogr
import HydroSEDPlots as HSplots

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
            Rain, prop, epsg = wmf.read_map_raster(PathRain) 
            if epsg == self.cuenca.epsg:
                Rain = self.cuenca.Transform_Map2Basin(Rain, prop)
            else:
                return 1, 1
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
        # Guarda el resultado 
        if len(PathQmed)>2:
            self.cuenca.Save_Net2Map(PathQmed, dxp, umbral, qmed = self.cuenca.CellQmed)
        #Retorna el resultado a la salida 
        return 0,self.cuenca.CellQmed[-1]

    def Basin_Update(self, PathNC):
        '''Actualiza el archivo nc de la cuenca con las variables agregadas o borradas de la misma'''
        #Lectura del archivo 
        g = netCDF4.Dataset(PathNC,'a')
        #Inclusion de nuevas variables
        for l in self.DicBasinNc.keys():#self.Nc2Save:
            #Si la variable no esta guradada la actualiza en el nc
            if self.DicBasinNc[l]['saved'] is False:
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
            if len(g.groups[grupoKey].variables.keys())>0:
                print g.groups[grupoKey].variables.keys()
                #itera
                for k in g.groups[grupoKey].variables.keys():
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
            if len(g.groups['Parametros'].variables.keys())>0:
                for k in g.groups['Parametros'].variables.keys():
                    self.DicParameters.update({k:{
                        'nombre': k,
                        'var': g.groups['Parametros'].variables[k][:]}})
        print self.DicParameters
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
        print Variable
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
            index=range(self.cuenca.ppal_stream.shape[1]), 
            columns=['Dist2Out[km]','Elevacion[m]'])
        #Curva hipsometrica
        CurvaHipso = pd.DataFrame(np.vstack([self.cuenca.hipso_basin[0]*wmf.cu.dxp**2/1e6, self.cuenca.hipso_basin[1]]).T,
            index = range(self.cuenca.hipso_basin.shape[1]),
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
            'tipo':Coef.dtype.name,
            'shape':Coef.shape,
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
        if u'h1_max' in self.DicBasinNc.keys(): 
            Hu = np.copy(self.DicBasinNc['h1_max']['var'])
        else: 
            Hu = 100    
        if u'v_coef' in self.DicBasinNc.keys():  
            ks = np.copy(self.DicBasinNc['v_coef']['var'][1])
        else: 
            ks = 0.003 
        self.cuenca.GetGeo_Cell_Basics()
        So = np.copy(self.cuenca.CellSlope)
        Factor = (wmf.cu.dxp**2.)/1000. #[m3/mm]
        #Calculo del coeficiente de kubota
        Coef = (ks*So*(wmf.cu.dxp**2.))/(3*(Hu*Factor)**2.)
        Coef[np.where(np.isinf(Coef))]=np.mean(Coef[np.where(np.isfinite(Coef))])
        print Coef
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
        
    def Basin_GeoGetRunoff(self, e1, Epsilon):
        '''Obtiene el coeficiente de escorrentia para carcavas'''
        #Obtiene pendiente y param requeridos
        self.cuenca.GetGeo_Cell_Basics()
        So = np.copy(self.cuenca.CellSlope)
        #Busca si la variable esta guardada en los diccionarios, si no, la reemplaza con una constante. 
        if 'Manning' in self.DicBasinWMF.keys() or 'Manning' in self.DicBasinNc.keys():
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
        print Coef
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
        xy,idShape = wmf.read_map_points(Path2Shp,[Campo2Read.encode()])
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
        for i in range(inicio, fin+1):
            vect,res = wmf.models.read_int_basin(path2bin,i,self.cuenca.ncells)
            if res == 0:
                vect = vect.astype(float)/1000.
                Vsum+=vect
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
           
    def Sim_SaveParameters(self, PathNC, ParamName, scalarParam):
        '''Actualiza el nc con un conjunto de parametros escalares nuevo'''
        #Abre el archivo nc y apunta al grupo de parametros
        g = netCDF4.Dataset(PathNC,'a')
        GrupoParam = g.groups['Parametros']
        Ejecuto = 1
        if len(ParamName)>0:
            #mira si el grupo ya tiene parametros adentro, si no crea la dimension
            try:
                pos = GrupoParam.dimensions.keys().index('nparam')
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
        Lista = self.DicBasinWMF.keys()
        Lista.extend(self.DicBasinNc.keys())
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
            
                
        
