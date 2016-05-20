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

from cu import *
from models import *
import numpy as np
from mpl_toolkits.basemap import Basemap, addcyclic, shiftgrid, cm
import pylab as pl
import osgeo.ogr, osgeo.osr
import gdal
from scipy.spatial import Delaunay
from scipy.stats import norm
import os
import pandas as pd
import datetime as datetime
try:
	import netcdf as netcdf
except:
	try:
		import netCDF4 as netcdf
	except:
		print 'No netcdf en esta maquina, se desabilita la funcion SimuBasin.save_SimuBasin'
		pass
#-----------------------------------------------------------------------
#Ploteo de variables
#-----------------------------------------------------------------------
def plot_sim_single(Qs,Qo=None,mrain=None,Dates=None,ruta=None,
	figsize=(8,4.5),ids=None,legend=True,
	ylabel='Stream $[m^3/seg]$'):
	fig=pl.figure(facecolor='w',edgecolor='w',figsize=figsize)
	ax1=fig.add_subplot(111)
	#Fechas
	if Dates==None:
		if len(Qs.shape)>1:
			ejeX=range(len(Qs[0]))
		else:
			ejeX = range(len(Qs))
	else:
		ejeX=Dates
	#Grafica la lluvia
	if mrain<>None:
		ax1AX=pl.gca()
		ax2=ax1.twinx()
		ax2AX=pl.gca()
		ax2.fill_between(ejeX,0,mrain,alpha=0.4,color='blue',lw=0)
		ax2.set_ylabel('Precipitation [$mm$]',size=14)
		ax2AX.set_ylim(ax2AX.get_ylim() [::-1])    
	#grafica las hidrografas
	ColorSim=['r','g','k','c','y']
	if len(Qs.shape)>1:
		if ids is None:
			ids = np.arange(1,Qs.shape[0]+1)
		for i,c,d in zip(Qs,ColorSim,ids):
			ax1.plot(ejeX,i,c,lw=1.5,label='Simulated '+str(d))    
	else:
		ax1.plot(ejeX,Qs,'r',lw=1.5,label='Simulated ')    
	if Qo<>None: ax1.plot(ejeX,Qo,'b',lw=2,label='Observed')
	#Pone elementos en la figura
	ax1.set_xlabel('Time [$min$]',size=14)
	ax1.set_ylabel(ylabel,size=14)
	ax1.grid(True)
	if legend:
		lgn1=ax1.legend(loc='upper center', bbox_to_anchor=(0.5, -0.12),
			fancybox=True, shadow=True, ncol=4)
	if ruta<>None:
		pl.savefig(ruta, bbox_inches='tight')
	pl.show()

#-----------------------------------------------------------------------
#Lectura de informacion y mapas 
#-----------------------------------------------------------------------
def read_map_raster(ruta_map,isDEMorDIR=False,dxp=None):
	'Funcion: read_map\n'\
	'Descripcion: Lee un mapa raster soportado por GDAL.\n'\
	'Parametros Obligatorios:.\n'\
	'	-ruta_map: Ruta donde se encuentra el mapa.\n'\
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
		cu.ncols=ncols
		cu.nrows=nrows
		cu.nodata=noData
		cu.dx=dx
		cu.xll=xll
		cu.yll=yll
		if dxp==None:
			cu.dxp=30.0
		else:
			cu.dxp=dxp
		return Mapa.T
	else:
		return Mapa.T,[ncols,nrows,xll,yll,dx,noData]

def read_map_points(ruta_map, ListAtr = None):
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
	'	Escribe el mapa en la ruta especificada.\n'\
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
			if i[0]<>-9999.0 and i[1]<>-9999.0:
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
	'	-FechaI, FechaF: Fechas inicial y final en formato datetime.\n'\
	'	-fmt : Formato en el que se pasa FechaI a texto, debe couincidir con.\n'\
	'		el del texto que se tenga en los archivos.\n'\
	'	-exten : Extension de los archivos .asc, .nc, .bin ...\n'\
	'	-string : texto antes de la fecha.\n'\
	'	-dt : Intervalos de tiempo entre los eventos.\n'\
	'Retorno:.\n'\
	'	Lista : la lista de python con los nombres de los binarios.\n'\
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

def __Add_hdr_bin_2route__(rute):
	if rute.endswith('.bin') is False and rute.endswith('.hdr') is True:
		ruteBin = rute[:-3] + 'bin'
		ruteHdr = rute
	elif rute.endswith('.bin') is True and rute.endswith('.hdr') is False:
		ruteBin = rute
		ruteHdr = rute[:-3] + 'hdr'
	elif  rute.endswith('.bin') is False and rute.endswith('.hdr') is False:
		ruteBin = rute + '.bin'
		ruteHdr = rute + '.hdr'
	return ruteBin,ruteHdr

def read_mean_rain(ruta,Nintervals=None,FirstInt=None):
	#Abrey cierra el archivo plano
	Data = np.loadtxt(ruta,skiprows=6,usecols=(2,3),delimiter=',',dtype='str')
	Rain = np.array([float(i[0]) for i in Data])
	Dates = [datetime.datetime.strptime(i[1],' %Y-%m-%d-%H:%M  ') for i in Data]
	#Obtiene el pedazo
	if Nintervals is not None or FirstInt is not None:
		Dates = Dates[FirstInt:FirstInt+Nintervals]
		Rain = Rain[FirstInt:FirstInt+Nintervals]
	#Regresa el resultado de la funcion
	return pd.Series(Rain,index = pd.to_datetime(Dates))
	
	
#-----------------------------------------------------------------------
#Ecuaciones Que son de utilidad
#-----------------------------------------------------------------------
def OCG_param(alfa=[0.75,0.2],sigma=[0.0,0.225,0.225],
	c1=5.54,k=0.5,fhi=0.95,Omega=0.13,pend=None,area=None,
	):
	'Funcion: OCG_param\n'\
	'Descripcion: Calcula los parametros de la onda cinematica.\n'\
	'geomorfologica (Velez, 2001).\n'\
	'Parametros Opcionales:.\n'\
	'	-isDEMorDIR: Pasa las propiedades de los mapas al modulo cuencas \n'\
	'		escrito en fortran \n'\
	'Retorno:.\n'\
	'	Parametros: B, w1, w2 y w3.\n'\
	'		si se entregan los mapas de pend, aacum entrega h_coef(4,:) .\n' \
	'		se asume que w1 corresponde a h_exp(4,:) .\n' \
	#Calcula los parametros de la ecuacion de onda cinematica 
	B = Omega*(c1*k**(alfa[0]-alfa[1]))**((2.0/3.0)-alfa[1])
	eB=1.0/(1+alfa[1]*((2/3.0)-sigma[1]))
	w1=((2/3.0)-sigma[1])*(1.0-alfa[1])*eB
	w2=(1+alfa[1]*((2/3.0)-sigma[1]))/(fhi*((2/3.0)-sigma[1])*(alfa[0]-alfa[1])+sigma[0])
	w2=(fhi*(0.667-sigma[1])*(alfa[1]-alfa[0])+sigma[0])*eB
	w3=(0.5-sigma[2])*eB
	B=B**(-eB)
	if pend<>None and area<>None:
		var = B*(pend**w2)*(area**w3)
		return var,w1
	else:
		return B,w1,w2,w3		

def PotCritica(S,D,te = 0.056):
    ti = te * (D*1600*9.8)
    return ti *(8.2* (((ti/(1000*9.8*S))/D)**(1.0/6.0)) * np.sqrt(ti/1000.0))/9800.0

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
	Y = pd.rolling_mean(Y,window)
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
	def __init__(self,lat,lon,DEM,DIR,name='NaN',stream=None,
		umbral=1000, ruta = None):
		'Descripcion: Inicia la variable de la cuenca, y la traza \n'\
		'	obtiene las propiedades basicas de la cuenca. \n'\
		'\n'\
		'Parametros\n'\
		'----------\n'\
		'self : Inicia las variables vacias.\n'\
		'lat : Coordenada en X de la salida de la cuenca.\n'\
		'lon : Coordenada en Y de la salida de la cuenca.\n'\
		'name : Nombre con el que se va a conocer la cuenca.\n'\
		'stream : Opcional, si se coloca, las coordenadas no tienen.\n'\
		'	que ser exactas, estas se van a corregir para ubicarse.\n'\
		'	en el punto mas cercano dentro de la corriente, este.\n'\
		'	debe ser un objeto del tipo stream.\n'\
		'umbral : umbral minimo para la creacion de cauce (defecto =1000).\n'\
		'ruta : Ruta donde se encuentra un archivo binario de la cuenca.\n'\
		'\n'\
		'Retornos\n'\
		'----------\n'\
		'self : Con las variables iniciadas.\n'\
		#Relaciona con el DEM y el DIR
		self.DEM=DEM
		self.DIR=DIR
		#Si se entrega cauce corrige coordenadas
		if ruta == None:
			if stream<>None:
				error=[]
				for i in stream.structure.T:
					error.append( np.sqrt((lat-i[0])**2+(lon-i[1])**2) )
				loc=np.argmin(error)
				lat=stream.structure[0,loc]
				lon=stream.structure[1,loc]
			#copia la direccion de los mapas de DEM y DIR, para no llamarlos mas
			self.name=name
			#Traza la cuenca 
			self.ncells = cu.basin_find(lat,lon,DIR,
				cu.ncols,cu.nrows)
			self.structure = cu.basin_cut(self.ncells)
			self.umbral = umbral
		else:
			self.__Load_BasinNc(ruta)
	#Cargador de cuenca 
	def __Load_BasinNc(self,ruta,Var2Search=None):
		'Descripcion: Lee una cuenca posteriormente guardada\n'\
		'	La cuenca debio ser guardada con Basin.Save_Basin2nc\n'\
		'\n'\
		'Parametros\n'\
		'----------\n'\
		'self : Inicia las variables vacias.\n'\
		'ruta : ruta donde se encuentra ubicada la cuenca guardada\n'\
		'Retornos\n'\
		'----------\n'\
		'self : La cuenca con sus parametros ya cargada.\n'\
		#Abre el archivo binario de la cuenca 
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
		'	se guardan como flotantes.\n'\
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
		rutaTcPlot = None, figsize=(8,5)):
		'Descripcion: Obtiene los parametros geomorfologicos de la cuenca \n'\
		'	y los tiempos de concentracion calculados por diferentes metodologias. \n'\
		'\n'\
		'Parametros\n'\
		'	rutaParamASC: ruta del ascii donde se escriben los param.\n'\
		'	plotTc: Plotea o no los tiempos de concentracion.\n'\
		'	rutaTcPlot: Si se da se guarda la figura de tiempos de concentracion.\n'\
		'----------\n'\
		'\n'\
		'Retornos\n'\
		'----------\n'\
		'GeoParameters : Parametros de la cuenca calculados.\n'\
		'Tc :  Tiempo de concentracion calculado para la cuenca.\n'\
		#Calcula lo que se necesita para sacar los parametros
		acum,longCeld,S0,Elev=cu.basin_basics(self.structure,
			self.DEM,self.DIR,cu.ncols,cu.nrows,self.ncells)
		slope=cu.basin_arc_slope(self.structure,self.DEM,self.ncells,
			cu.ncols,cu.nrows)
		Lpma,puntto=cu.basin_findlong(self.structure,self.ncells)
		cauce,nodos,trazado,n_nodos,n_cauce = cu.basin_stream_nod(self.structure,
			acum,self.umbral,self.ncells)
		ppal_nceldas,punto = cu.basin_ppalstream_find(self.structure,
			nodos,longCeld,Elev,self.ncells)
		ppal = cu.basin_ppalstream_cut(ppal_nceldas,self.ncells)
		self.hipso_main,self.hipso_basin=cu.basin_ppal_hipsometric(
			self.structure,Elev,punto,30,ppal_nceldas,self.ncells)
		self.main_stream=ppal
		nperim = cu.basin_perim_find(self.structure,self.ncells)
		#Obtiene los parametros 
		Area=(self.ncells*cu.dxp**2)/1e6
		Perim=nperim*cu.dxp/1000.0
		Lcau=ppal[1,-1]/1000.0
		Scau=np.polyfit(ppal[1,::-1],ppal[0],1)[0]*100
		Scue=slope.mean()*100
		Hmin=Elev[-1]; Hmax=Elev[puntto]; Hmean=Elev.mean()
		HCmax=Elev[punto]
		x,y = cu.basin_coordxy(self.structure,self.ncells)
		CentXY = [np.median(x),np.median(y)]
		#Genera un diccionario con las propiedades de la cuenca 
		self.GeoParameters={'Area[km2]': Area,
			'Perimetro[km]':Perim,
			'Pend Cauce [%]':Scau,
			'Long Cau [km]': Lcau,
			'Pend Cuenca [%]': Scue,
			'Long Cuenca [km]': Lpma,
			'Hmax [m]':Hmax,
			'Hmin [m]':Hmin,
			'Hmean [m]':Hmean,
			'H Cauce Max [m]':HCmax,
			'Centro [X]': CentXY[0],
			'Centro [Y]': CentXY[1]}
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
	
	#Parametros por mapas (distribuidos)
	def GetGeo_Cell_Basics(self):
		'Descripcion: Obtiene: area acumulada, long de celdas, Pendiente \n'\
		'	y Elevacion en cada elemento de la cuenca. \n'\
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
			self.DEM,self.DIR,cu.ncols,cu.nrows,self.ncells)
		self.CellAcum=acum; self.CellLong=longCeld
		self.CellSlope=S0; self.CellHeight=Elev
		#Obtiene el canal en la cuenca 
		self.CellCauce = np.zeros(self.ncells)
		self.CellCauce[self.CellAcum>self.umbral]=1
	def GetGeo_IsoChrones(self,Tc,Niter=4):
		'Descripcion: Obtiene el tiempo de viaje aproximado de cada  \n'\
		'	celda a la salida de la cuenca, para eso usa el tiempo de . \n'\
		'	concentracion obtenido por la funcion GetGeo_Parameters . \n'\
		'\n'\
		'Parametros\n'\
		'----------\n'\
		'	self : no necesita nada es autocontenido.\n'\
		'	Tc : Valor escalar de tiempo de concentracion.\n'\
		'	Niter: Cantidad de iteraciones para aproximar vel, defecto 4.\n'\
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
	def GetGeo_Ppal_Hipsometric(self,umbral=1000,
		intervals = 30):
		'Descripcion: Calcula y grafica la curva hipsometrica de\n'\
		'	la cuenca.\n'\
		'\n'\
		'Parametros\n'\
		'----------\n'\
		'umbral : cantidad minima de celdas para el trazado.\n'\
		'intervals: Cantidad de intervalos en los cuales se haran .\n'\
		'	los muestreos de la curva hipsometrica.\n'\
		'\n'\
		'Retornos\n'\
		'----------\n'\
		'Hipso : Curva hipsometrica en formato vectorial, contiene.\n'\
		'	- Curva sobre cauce ppal.\n'\
		'	- Curva como histograma de la cuenca.\n'\
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
	
	def GetGeo_HAND(self,umbral=1000):
		'Descripcion: Calcula Height Above the Nearest Drainage (HAND) \n'\
		'	y Horizontal Distance to the Nearest Drainage (HDND) (Renno, 2008). \n'\
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
			self.DEM,self.DIR,cu.ncols,cu.nrows,self.ncells)		
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
		if CN<>None:
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
			#if Tr==None or len(Tr)<>IntTr.size():
			#	Lista=[2.33,5,10,25,50,100,500,1000]
			#	Tr=[l for l in Lista[:IntTr.size]]			
			fig=pl.figure(edgecolor='w',facecolor='w')
			ax=fig.add_subplot(111)
			X=np.linspace(0,Dur,lluviaTrEfect.shape[1])
			Grosor=np.arange(0.5,4,0.2)
			for l,le,t,g in zip(lluviaTr,lluviaTrEfect,Tr,Grosor):	
				ax.plot(X,l,c='b',lw=g,label=str(t))
			if CN<>None:
				for l,le,t,g in zip(lluviaTr,lluviaTrEfect,Tr,Grosor):	
					ax.plot(X,le,c='r',lw=g,label=str(t))
			ax.set_xlabel('Tiempo $[h]$',size=14)
			ax.set_ylabel('Precipitacion $[mm]$',size=14)	
			ax.grid(True)
			pl.legend(loc=0,ncol=2)
			if ruta<>None:
				pl.savefig(ruta,bbox_inches='tight')
			pl.show()
		if CN<>None:
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
		if ruta<>None:
			pl.savefig(ruta,bbox_inches='tight')
		pl.show()
	
	#------------------------------------------------------
	# Trabajo con mapas externos y variables fisicas
	#------------------------------------------------------
	def Transform_Map2Basin(self,Map,MapProp):
		'Descripcion: A partir de un mapa leido obtiene un vector \n'\
		'	con la forma de la cuenca, el cual luego puede ser agregado a esta. \n'\
		'\n'\
		'Parametros\n'\
		'----------\n'\
		'self : Inicia las variables vacias.\n'\
		'Map : Matriz con la informacion del mapa.\n'\
		'MapProp : Propiedades del mapa.\n'\
		'	1. Ncols Mapa.\n'\
		'	2. Nrows Mapa.\n'\
		'	3. Xll Mapa.\n'\
		'	4. Yll Mapa.\n'\
		'	5. dx Mapa.\n'\
		'\n'\
		'Retornos\n'\
		'----------\n'\
		'vecMap : Vector conla informacion del mapa al interio de la cuenca.\n'\
		#Comienza le codifgo 
		vec = cu.basin_map2basin(self.structure,
			Map,MapProp[2],MapProp[3],MapProp[4],
			cu.nodata,
			'fill_mean',
			self.ncells,
			MapProp[0],MapProp[1])
		return vec
	def Transform_Hills2Basin(self,HillsMap):
		'Descripcion: A partir de un vector con propiedades de las laderas\n'\
		'	obtiene un vector con las propiedades por celda, ojo estas \n'\
		'	quedan con las formas de las laderas y la variable queda agregada. \n'\
		'\n'\
		'Parametros\n'\
		'----------\n'\
		'self : la cuenca misma.\n'\
		'MapHills : Vector con las variables por laderas [nhills].\n'\
		'\n'\
		'Retornos\n'\
		'----------\n'\
		'CellMap : Vector con la variable agregada por laderas, pero .\n'\
		'	pasada a celdas.\n'\
		#Genera el mapa de basin vacio
		CellMap = np.ones(self.ncells)
		#itera por la cantidad de elementos y les va asignando
		for i,k in enumerate(HillsMap[::-1]):
			CellMap[self.hills_own==i+1] = k
		return CellMap
	def Transform_Basin2Hills(self,CellMap,mask=None,sumORmean=0):
		'Descripcion: A partir de un vector tipo Basin obtiene un\n'\
		'	vector del tipo laderas, en donde las propiedades se \n'\
		'	agregan para cada ladera. \n'\
		'\n'\
		'Parametros\n'\
		'----------\n'\
		'self : la cuenca misma.\n'\
		'CellMap : Vector con las propiedades por celdas [ncells].\n'\
		'mask : Celdas sobre las cuales se agrega la variable (1), y\n'\
		'	sobre las que no (0).\n'\
		'sumORmean : si la variable sera agregada como un promedio (0)\n'\
		'	o como una suma (1).\n'\
		'\n'\
		'Retornos\n'\
		'----------\n'\
		'HillsMap : Vector con las prop agregadas a laderas .\n'\
		#Si hay mascara la tiene en cuenta
		if mask<>None:
			Ma = np.zeros(self.ncells)
			if type(mask) is float or type(mask) is int:
				Ma[CellMap==mask] = 1
			elif type(mask) is np.ndarray:
				Ma = np.copy(mask)
		else:
			Ma = np.ones(self.ncells)
		#Pasa el mapa de celdas a mapa de laderas		
		HillsMap = cu.basin_subbasin_map2subbasin(self.hills_own,
			CellMap, self.nhills, sumORmean, self.ncells, Ma)
		return HillsMap
	#------------------------------------------------------
	# Trabajo con datos puntuales puntos 
	#------------------------------------------------------
	def Points_Points2Stream(self,coordXY,ids):
		'Descripcion: toma las coordenadas de puntos XY y las mueve\n'\
		'	hacia los cauces de la cuenca\n'\
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
		'	donde quedaron localizadas las coordenadas.\n'\
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
		return xy_new, basin_pts, basin_pts[basin_pts<>0]
	def Points_Points2Basin(self,coordXY,ids):
		'Descripcion: toma las coordenadas de puntos XY y las pone\n'\
		'	en la cuenca, no las mueve hacia los cuaces\n'\
		'\n'\
		'Parametros\n'\
		'----------\n'\
		'coordXY : coordenadas de los puntos [2,Ncoord].\n'\
		'ids : Identificacion de los puntos que se van a mover.\n'\
		'\n'\
		'Retornos\n'\
		'----------\n'\
		'basin_pts : Vector con la forma de la cuenca con los puntos.\n'\
		'	donde quedaron localizadas las coordenadas.\n'\
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
		return basin_pts,basin_pts[basin_pts<>0]
	#------------------------------------------------------
	# Caudales de largo plazo y regionalizacion
	#------------------------------------------------------
	#Caudal de largo plazo 
	def GetQ_Balance(self,Precipitation, Tipo_ETR = 1, mu_choud = 1.37):
		'Descripcion: Calcula el caudal medio por balance de largo plazo\n'\
		'	para ello requiere conocer la precipitacion y el metodo de\n'\
		'	estimacion de la evaporacion.\n'\
		'\n'\
		'Parametros\n'\
		'----------\n'\
		'self : no necesita nada es autocontenido.\n'\
		'Precipitation : Cantidad anual de lluvia, escalar o vector.\n'\
		'Elevacion : Elevacion en cada punto de la cuenca.\n'\
		'Tipo_ETR : Tipo de ecuacion para calcular la evaporacion.\n'\
		'	-1. Turc.\n'\
		'	-2. Cenicafe Budyko.\n'\
		'	-3. Choundry.\n'\
		'	Defecto: 1.\n'\
		'\n'\
		'Retornos\n'\
		'----------\n'\
		'CellQmed : Caudal medio calculado para toda la cuenca.\n'\
		'CellETR : ETR calculada para toda la cuenca.\n'\
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
		Tr = [2.33, 5, 10, 25, 50, 100], Dist = 'gumbel'):
		'Descripcion: Calcula el caudal maximo para diferentes\n'\
		'	periodos de retorno. Recibe como entrada el caudal medio \n'\
		'	calculado con GetQ_Balance.\n'\
		'	el calculo se hace a partir de la ecuacion:.\n'\
		'	MedMax = Coef[0] * Qmed ** Exp[0] .\n'\
		'	DesMax = Coef[1] * Qmed ** Exp[1] .\n'\
		'	Qmax = MedMax + K(Tr) * DesMax.\n'\
		'	Donde K es un coeficiente que depende de Tr y f(x).\n'\
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
		MedMax = Coef[0] * Qmed ** Expo[0]
		DesMax = Coef[1] * Qmed ** Expo[1]
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
		'	periodos de retorno. Recibe como entrada el caudal medio \n'\
		'	calculado con GetQ_Balance.\n'\
		'	el calculo se hace a partir de la ecuacion:.\n'\
		'	MedMin = Coef[0] * Qmed ** Exp[0] .\n'\
		'	DesMin = Coef[1] * Qmed ** Exp[1] .\n'\
		'	Qmin = MedMin + K(Tr) * DesMin.\n'\
		'	Donde K es un coeficiente que depende de Tr y f(x).\n'\
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
	def Save_Net2Map(self,ruta,dx=30.0,umbral=1000,
		qmed=None,Dict=None,DriverFormat='ESRI Shapefile',
		EPSG=4326):
		'Descripcion: Guarda la red hidrica simulada de la cuenca en .shp \n'\
		'	Puede contener un diccionario con propiedades de la red hidrica. \n'\
		'\n'\
		'Parametros\n'\
		'----------\n'\
		'self : no necesita nada es autocontenido.\n'\
		'ruta : Lugar y nombre donde se va a guardar la red hidrica.\n'\
		'dx : Longitud de las celdas planas.\n'\
		'umbral : cantidad de celdas necesarias para corriente.\n'\
		'qmed : caudal medio calculado por alguna metodologia.\n'\
		'Dict : Diccionario con parametros de la red hidrica que se quieren imprimir.\n'\
		'DriverFormat : nombre del tipo de archivo vectorial de salida (ver OsGeo).\n'\
		'EPSG : Codigo de proyeccion utilizada para los datos, defecto WGS84.\n'\
		'\n'\
		'Retornos\n'\
		'----------\n'\
		'Escribe un archivo vectorial con la estructura de la red hidrica y sus propiedades.\n'\
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
		if qmed<>None:
			netsize = cu.basin_netxy_find(self.structure,nodos,cauce*qmed,self.ncells)
			netQmed=cu.basin_netxy_cut(netsize,self.ncells)
		cortes=np.where(net[0,:]==-999)
		cortes=cortes[0].tolist()
		cortes.insert(0,0)
		#Escribe el shp de la red hidrica
		if ruta.endswith('.shp')==False:
			ruta=ruta+'.shp'
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
		if qmed<>None:
			new_field=osgeo.ogr.FieldDefn('Qmed[m3s]',osgeo.ogr.OFTReal)
			layer.CreateField(new_field)
		if Dict<>None:
			if type(Dict==dict):
				netDict=[]
				for k in Dict.keys():
					new_field=osgeo.ogr.FieldDefn(k[:10],osgeo.ogr.OFTReal)
					layer.CreateField(new_field)
					netsizeT = cu.basin_netxy_find(self.structure,nodos,cauce*Dict[k],self.ncells)
					netDict.append(cu.basin_netxy_cut(netsize,self.ncells))
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
			if qmed<>None:	
				feature.SetField('Qmed[m3s]',float(netQmed[0,i+1]))
			if Dict<>None:
				if type(Dict==dict):
					for n,k in zip(netDict,Dict.keys()):					
						feature.SetField(k[:10],float(n[0,i+1]))
			#featureFID+=1
			layer.CreateFeature(feature)
			line.Destroy()
			feature.Destroy()
		shapeData.Destroy()
	def Save_Basin2Map(self,ruta,dx=30.0,Param={},
		DriverFormat='ESRI Shapefile',EPSG=4326):
		'Descripcion: Guarda un archivo vectorial de la cuenca en .shp \n'\
		'	Puede contener un diccionario con propiedades. \n'\
		'\n'\
		'Parametros\n'\
		'----------\n'\
		'self : no necesita nada es autocontenido.\n'\
		'ruta : Lugar y nombre donde se va a guardar la cuenca.\n'\
		'dx : Longitud de las celdas planas.\n'\
		'DriverFormat : nombre del tipo de archivo vectorial de salida (ver OsGeo).\n'\
		'EPSG : Codigo de proyeccion utilizada para los datos, defecto WGS84.\n'\
		'\n'\
		'Retornos\n'\
		'----------\n'\
		'Escribe un archivo vectorial de la cuenca.\n'\
		#Obtiene el perimetro de la cuenca 
		nperim = cu.basin_perim_find(self.structure,self.ncells)
		basinPerim=cu.basin_perim_cut(nperim)
		#Param,Tc=basin_Tc(basin.structure,DEM,DIR,cu.dxp,basin.ncells,cu.ncols,cu.nrows)
		#Construye el shp 
		if ruta.endswith('.shp')==False:
			ruta=ruta+'.shp'
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
		#Calcula el tamano de la muestra
		ring = osgeo.ogr.Geometry(osgeo.ogr.wkbLinearRing)
		for i in basinPerim.T:
			ring.AddPoint(x=float(i[0]),y=float(i[1]))
		poly=osgeo.ogr.Geometry(osgeo.ogr.wkbPolygon)
		poly.AddGeometry(ring)
		feature = osgeo.ogr.Feature(layerDefinition)
		feature.SetGeometry(poly)
		feature.SetFID(0)
		for p in Param.keys():		
			#feature.SetField(p[:p.index('[')].strip()[:10],float("%.2f" % Param[p]))
			feature.SetField(p,float("%.2f" % Param[p]))
		layer.CreateFeature(feature)
		poly.Destroy()
		ring.Destroy()
		feature.Destroy()
		shapeData.Destroy()

	#------------------------------------------------------
	# Graficas de la cuenca
	#------------------------------------------------------
	def Plot_basin(self,vec=None,Min=None,
		Max=None,ruta=None,mostrar='si',barra='si',figsize=(10,8),
		ZeroAsNaN = 'no',extra_lat=0.0,extra_long=0.0,lines_spaces=0.02,
		xy=None,xycolor='b',colorTable=None,alpha=1.0,vmin=None,vmax=None,
		colorbar=True):
		#Plotea en la terminal como mapa un vector de la cuenca
		'Funcion: write_proyect_int_ext\n'\
		'Descripcion: Genera un plot del mapa entrgeado.\n'\
		'del mismo en forma de mapa \n'\
		'Parametros Obligatorios:.\n'\
		'	-basin: Vector con la forma de la cuenca.\n'\
		'	-vec: Vector con los valores a plotear.\n'\
		'Parametros Opcionales:.\n'\
		'	-Min: Valor minimo del plot, determina el rango de colores.\n'\
		'	-Max: Valor maximo del plot, determina el rango de colores.\n'\
		'	-ruta: Ruta en la cual se guarda la grafica.\n'\
		'	-mostrar: Muestra o no la grafica, defecto: si.\n'\
		'	-barra: Muestra o no la barra de colores, defecto: si.\n'\
		'	-figsize: tamano de la ventana donde se muestra la cuenca.\n'\
		'	-ZeroAsNaN: Convierte los valores de cero en NaN.\n'\
		'Retorno:.\n'\
		'	Actualizacion del binario .int\n'\
		#Plot en basemap
		Mcols,Mrows=cu.basin_2map_find(self.structure,self.ncells)
		Map,mxll,myll=cu.basin_2map(self.structure,self.structure[0]
			,Mcols,Mrows,self.ncells)
		longs=np.array([mxll+0.5*cu.dx+i*cu.dx for i in range(Mcols)])
		lats=np.array([myll+0.5*cu.dx+i*cu.dx for i in range(Mrows)])
		X,Y=np.meshgrid(longs,lats)
		Y=Y[::-1]
		fig=pl.figure(figsize=figsize)		
		m = Basemap(projection='merc',
			llcrnrlat=lats.min()-extra_lat,
			urcrnrlat=lats.max()+extra_lat,
			llcrnrlon=longs.min()-extra_long,
			urcrnrlon=longs.max()+extra_long,
			resolution='c')
		m.drawparallels(np.arange(lats.min(),
			lats.max(),lines_spaces),
			labels=[1,0,0,0],
			fmt="%.2f",
			rotation='vertical',
			xoffset=0.1)
		m.drawmeridians(np.arange(longs.min(),
			longs.max(),lines_spaces),
			labels=[0,0,1,0], 
			fmt="%.2f",
			yoffset=0.1)
		Xm,Ym=m(X,Y)
	    #Plotea el contorno de la cuenca y la red 
		nperim = cu.basin_perim_find(self.structure,self.ncells)
		perim = cu.basin_perim_cut(nperim)	
		xp,yp=m(perim[0],perim[1])  
		m.plot(xp, yp, color='r',lw=2)
	    #Si hay una variable la monta
		if vec<>None:
			if vmin is None:
				vmin = vec.min()
			if vmax is None:
				vmax = vec.max()
			MapVec,mxll,myll=cu.basin_2map(self.structure,vec,Mcols,Mrows,
				self.ncells)
			MapVec[MapVec==cu.nodata]=np.nan
			if ZeroAsNaN is 'si':
				MapVec[MapVec == 0] = np.nan
			if colorTable<>None:
				cs = m.contourf(Xm, Ym, MapVec.T, 25, alpha=alpha,cmap=colorTable,
					vmin=vmin,vmax=vmax)
			else:
				cs = m.contourf(Xm, Ym, MapVec.T, 25, alpha=alpha,
					vmin=vmin,vmax=vmax)
			if colorbar:
				cbar = m.colorbar(cs,location='bottom',pad="5%")	
		#Si hay coordenadas de algo las plotea
		if xy<>None:
			xc,yc=m(xy[0],xy[1])
			m.scatter(xc,yc,color=xycolor,
				s=30,
				linewidth=0.5,
				edgecolor='black')
		#Guarda
		if ruta<>None:
			pl.savefig(ruta, bbox_inches='tight',pad_inches = 0.25)
		pl.show()
	# Grafica barras de tiempos de concentracion
	def Plot_Tc(self,ruta=None,figsize=(8,6)):
		keys=self.Tc.keys()
		keys[2]=u'Carr Espana'
		Media=np.array(self.Tc.values()).mean()
		Desv=np.array(self.Tc.values()).std()
		Mediana=np.percentile(self.Tc.values(),50)
		rango=[Media-Desv,Media+Desv]
		colores=[]
		for t in self.Tc.values():
			if t>rango[0] and t<rango[1]:
				colores.append('b')
			else:
				colores.append('r')
		fig=pl.figure(edgecolor='w',facecolor='w',figsize=figsize)
		ax=fig.add_subplot(111)
		box = ax.get_position()
		ax.set_position([box.x0, box.y0 + box.height * 0.18,
			box.width, box.height * 0.9])
		ax.set_xlim(-0.4,len(keys)+1-0.8)
		ax.bar(range(len(keys)),self.Tc.values(),color=colores)
		ax.hlines(Media,-0.4,len(keys)+1-0.8,'k',lw=2)
		ax.hlines(Mediana,-0.4,len(keys)+1-0.8,'r',lw=2)
		ax.hlines(Media+Desv,-0.4,len(keys)+1-0.8,'b',lw=2)
		Texto='%.2f' % Media
		ax.text(len(keys)/3.0,Media+0.03,'$\\mu='+Texto+'$')
		Texto='%.2f' % Desv
		ax.text(len(keys)/2.0,Media+Desv+0.03,u'$\mu+\sigma='+Texto+'$')
		Texto='%.2f' % Mediana
		ax.text(len(keys)/2.0,Mediana+0.03,'$P_{50}='+Texto+'$')	
		ax.set_xticks(list(np.arange(1,len(keys)+1)-0.8))
		ax.set_xticklabels(keys,rotation=60)
		ax.set_ylabel(u'Tiempo de concentracion $T_c[hrs]$',size=14)
		ax.grid(True)
		if ruta<>None:
			pl.savefig(ruta,bbox_inches='tight')
		pl.show()
	#Plot de cuace ppal
	def PlotPpalStream(self,ruta = None, figsize = (8,6)):
		fig = pl.figure(figsize = figsize, edgecolor = 'w',
			facecolor = 'w')
		ax = fig.add_subplot(111)
		pl.scatter(self.ppal_stream[1]/1000.0,
			self.ppal_stream[0],
			c = self.ppal_slope,
			linewidth = 0)
		ax.grid(True)
		ax.set_xlabel('Distancia $[km]$',size = 16)
		ax.set_ylabel(u'Elevacion $[m.s.n.m]$',size = 16)
		pl.colorbar()
		if ruta<>None:
			pl.savefig(ruta,bbox_inches='tight')
		pl.show()
	#Plot de hidtograma de pendientes 
	def PlotSlopeHist(self,ruta=None,bins=[0,2,0.2],
		Nsize=1000, figsize = (8,6)):
		if Nsize>self.ncells:
			Nsize = self.ncells
		pos = np.random.choice(self.ncells,Nsize)
		h,b=np.histogram(self.CellSlope[pos],bins=np.arange(bins[0],bins[1],bins[2]))
		b=(b[:-1]+b[1:])/2.0
		h=h.astype(float)/h.astype(float).sum()
		fig=pl.figure(figsize = figsize, edgecolor='w',facecolor='w')
		ax=fig.add_subplot(111)
		ax.plot(b,h,lw=2)
		ax.grid(True)
		ax.set_xlabel('Pendiente',size=14)
		ax.set_ylabel('$pdf [\%]$',size=14)
		if ruta<>None:
			pl.savefig(ruta,bbox_inches='tight')
		pl.show()
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
		if ruta<>None:
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
		box = ax.get_position()
		ax.set_position([box.x0, box.y0 + box.height * 0.1,
			box.width, box.height * 0.9])
		ax.plot(ppal_acum,elevPpal,c='b',lw=2,label='Cacuce Principal')
		ax.plot(basin_acum,elevBasin,c='r',lw=2,label='Cuenca')
		ax.grid()
		ax.set_xlabel('Porcentaje Area Acumulada $[\%]$',size=14)
		if normed==False:
			ax.set_ylabel('Elevacion $[m.s.n.m]$',size=14)
		elif normed==True:
			ax.set_ylabel('Elevacion $[\%]$',size=14)
		lgn1=ax.legend(loc='upper center', bbox_to_anchor=(0.5, -0.12),
			fancybox=True, shadow=True, ncol=2)
		if ruta<>None:
			pl.savefig(ruta, bbox_inches='tight')
			pl.close('all')
		else:
			pl.show()
class SimuBasin(Basin):
	
	def __init__(self,lat,lon,DEM,DIR,name='NaN',stream=None,umbral=500,
		noData=-999,modelType='cells',SimSed='no',SimSlides='no',dt=60,
		SaveStorage='no',SaveSpeed='no',rute = None, retorno = 0,
		SeparateFluxes = 'no'):
		'Descripcion: Inicia un objeto para simulacion \n'\
		'	el objeto tiene las propieades de una cuenca con. \n'\
		'	la diferencia de que inicia las variables requeridas. \n'\
		'	para simular. \n'\
		'\n'\
		'Parametros\n'\
		'----------\n'\
		'self : Inicia las variables vacias.\n'\
		'lat : Coordenada en X de la salida de la cuenca.\n'\
		'lon : Coordenada en Y de la salida de la cuenca.\n'\
		'name : Nombre con el que se va a conocer la cuenca.\n'\
		'	(defecto = NaN).\n'\
		'stream : Opcional, si se coloca, las coordenadas no tienen.\n'\
		'	que ser exactas, estas se van a corregir para ubicarse.\n'\
		'	en el punto mas cercano dentro de la corriente, este.\n'\
		'	debe ser un objeto del tipo stream.\n'\
		'umbral : Cantidad minima de celdas para la creacion de cauces.\n'\
		'	(defecto = 500 ).\n'\
		'noData : Valor correspondiente a valores sin dato (defecto = -999).\n'\
		'modelType : Tipo de modelo, por celdas o por laderas (defecto = cells).\n'\
		'	opciones: .\n'\
		'		cells => modela por celdas.\n'\
		'		hills => modela por laderas.\n'\
		'SimSed : Simula si, o no simula sedimentos no.\n'\
		'SimSlides : Simula si, o no simula deslizamientos no.\n'\
		'dt : Tamano del intervlao de tiempo en que trabaj el modelo (defecto=60seg) [secs].\n'\
		'SaveStorage : Guarda o no el almacenamiento.\n'\
		'SaveSpeed : Guarda o no mapas de velocidad.\n'\
		'rute : por defecto es None, si se coloca la ruta, el programa no.\n'\
		'	hace una cuenca, si no que trata de leer una en el lugar especificado.\n'\
		'	Ojo: la cuenca debio ser guardada con la funcion: Save_SimuBasin().\n'\
		'retorno : (defecto = 0), si es cero no se considera alm maximo en .\n'\
		'	el tanque 3, si es 1, si se considera.\n'\
		'SeparateFluxes : Separa el flujo en base, sub-superficial y escorrentia.\n'\
		'SeparateRain : Separa el flujo proveniente de convectivas y de estratiformes.\n'\
		'ShowStorage : Muestra en la salida del modelo el alm promedio en cada uno de los tanques.\n'\
		'\n'\
		'Retornos\n'\
		'----------\n'\
		'self : Con las variables iniciadas.\n'\
		#Variables de radar
		self.radarDates = []
		self.radarPos = []
		self.radarMeanRain = []
		self.radarCont = 1
		#Si no hay ruta traza la cuenca
		if rute is None:
			#Si se entrega cauce corrige coordenadas
			if stream<>None:
				error=[]
				for i in stream.structure.T:
					error.append( np.sqrt((lat-i[0])**2+(lon-i[1])**2) )
				loc=np.argmin(error)
				lat=stream.structure[0,loc]
				lon=stream.structure[1,loc]
			#copia la direccion de los mapas de DEM y DIR, para no llamarlos mas
			self.name=name
			self.DEM=DEM
			self.DIR=DIR
			self.modelType=modelType
			self.nodata=noData
			self.umbral = umbral
			#Traza la cuenca 
			self.ncells = cu.basin_find(lat,lon,DIR,
				cu.ncols,cu.nrows)
			self.structure = cu.basin_cut(self.ncells)
			#traza las sub-cuencas
			acum=cu.basin_acum(self.structure,self.ncells)
			cauce,nodos,self.nhills = cu.basin_subbasin_nod(self.structure
				,acum,umbral,self.ncells)
			self.hills_own,sub_basin = cu.basin_subbasin_find(self.structure,
				nodos,self.nhills,self.ncells)
			self.hills = cu.basin_subbasin_cut(self.nhills)
			models.drena=self.structure		
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
			models.retorno = 0
			models.verbose = 0
			#Define los puntos de control		
			models.control = np.zeros((1,N))
			models.control_h = np.zeros((1,N))
			#Define las simulaciones que se van a hacer 
			models.sim_sediments=0
			if SimSed is 'si':
				models.sim_sediments=1
			models.sim_slides=0
			if SimSlides is 'si':
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
		# si hay tura lee todo lo de la cuenca
		elif rute is not None:
			self.__Load_SimuBasin(rute)
	
	def __Load_SimuBasin(self,ruta):
		'Descripcion: Lee una cuenca posteriormente guardada\n'\
		'	La cuenca debio ser guardada con SimuBasin.save_SimuBasin\n'\
		'\n'\
		'Parametros\n'\
		'----------\n'\
		'self : Inicia las variables vacias.\n'\
		'ruta : ruta donde se encuentra ubicada la cuenca guardada\n'\
		'Retornos\n'\
		'----------\n'\
		'self : La cuenca con sus parametros ya cargada.\n'\
		#Abre el archivo binario de la cuenca 
		gr = netcdf.Dataset(ruta,'a')
		#obtiene las prop de la cuenca
		self.name = gr.nombre
		self.modelType = gr.modelType.encode()
		self.nodata = gr.noData
		self.umbral = gr.umbral
		self.ncells = gr.ncells
		self.nhills = gr.nhills
		models.dt = gr.dt
		models.retorno = gr.retorno
		#Asigna dem y DIr a partir de la ruta 
		try:
			DEM = read_map_raster(gr.DEM,True,gr.dxp)
			DIR = read_map_raster(gr.DIR,True,gr.dxp)
			cu.nodata = -9999.0
			DIR[DIR<=0]=cu.nodata.astype(int)
			DIR=cu.dir_reclass(DIR,cu.ncols,cu.nrows)			
			self.DEM = DEM
			self.DIR = DIR
		except:
			print 'me jodi'
			pass
		#de acuerdo al tipo de modeloe stablece numero de elem
		if self.modelType[0] is 'c':
			N = self.ncells
		elif self.modelType[0] is 'h':
			N = self.nhills
		#Obtiene las variables vectoriales 
		self.structure = gr.variables['structure'][:]
		self.hills = gr.variables['hills'][:]
		self.hills_own = gr.variables['hills_own'][:]
		#obtiene las propieades del modelo 
		models.h_coef = np.ones((4,N)) * gr.variables['h_coef'][:]
		models.v_coef = np.ones((4,N)) * gr.variables['v_coef'][:]
		models.h_exp = np.ones((4,N)) * gr.variables['h_exp'][:]
		models.v_exp = np.ones((4,N)) * gr.variables['v_exp'][:]
		models.max_capilar = np.ones((1,N)) * gr.variables['h1_max'][:]
		models.max_gravita = np.ones((1,N)) * gr.variables['h3_max'][:]
		
		if self.modelType[0] is 'c':			
			models.drena = np.ones((3,N)) *gr.variables['drena'][:]
		elif self.modelType[0] is 'h':
			models.drena = np.ones((1,N)) * gr.variables['drena'][:]			
		models.unit_type = np.ones((1,N)) * gr.variables['unit_type'][:]
		models.hill_long = np.ones((1,N)) * gr.variables['hill_long'][:]
		models.hill_slope = np.ones((1,N)) * gr.variables['hill_slope'][:]
		models.stream_long = np.ones((1,N)) * gr.variables['stream_long'][:]
		models.stream_slope = np.ones((1,N)) * gr.variables['stream_slope'][:]
		models.stream_width = np.ones((1,N)) * gr.variables['stream_width'][:]
		models.elem_area = np.ones((1,N)) * gr.variables['elem_area'][:]
		models.speed_type = np.ones((3)) * gr.variables['speed_type'][:]
		models.storage = np.ones((5,N)) * gr.variables['storage'][:]
		
		#propiedades de puntos de control
		models.control = np.ones((1,N)) * gr.variables['control'][:]
		models.control_h = np.ones((1,N)) * gr.variables['control_h'][:]
		#Cierra el archivo 
		gr.close()
		
	#------------------------------------------------------
	# Subrutinas de lluvia, interpolacion, lectura, escritura
	#------------------------------------------------------	
	def rain_interpolate_mit(self,coord,registers,ruta):
		'Descripcion: Interpola la lluvia mediante una malla\n'\
		'	irregular de triangulos, genera campos que son. \n'\
		'	guardados en un binario para luego ser leido por el. \n'\
		'	modelo en el momento de simular. \n'\
		'\n'\
		'Parametros\n'\
		'----------\n'\
		'self : .\n'\
		'coord : Array (2,Ncoord) con las coordenadas de estaciones.\n'\
		'registers : Array (Nest,Nregisters) con los registros de lluvia.\n'\
		'ruta : Ruta con nombre en donde se guardara el binario con.\n'\
		'	la informacion de lluvia.\n'\
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
		TIN_perte = models.rain_pre_mit(xy_basin,TIN_mesh,coord,self.ncells,
			TIN_mesh.shape[1],coord.shape[1]) 	 			
		#Genera las interpolaciones para el rango de datos 		
		meanRain = models.rain_mit(xy_basin,coord,reg,TIN_mesh,
			TIN_perte,ruta,self.ncells,coord.shape[1],
			TIN_mesh.shape[1],reg.shape[1])
		#Guarda un archivo con informacion de la lluvia 
		f=open(ruta[:-3]+'hdr','w')
		f.write('Numero de celdas: %d \n' % self.ncells)
		f.write('Numero de laderas: %d \n' % self.nhills)
		f.write('Numero de registros: %d \n' % reg.shape[1])
		f.write('Tipo de interpolacion: TIN\n')
		f.write('Record, Fecha \n')
		if isPandas:
			dates=registers.index.to_pydatetime()
			for c,d in enumerate(dates):
				f.write('%d, %s \n' % (c,d.strftime('%Y-%m-%d-%H:%M')))
		f.close()
		return meanRain
			
	def rain_interpolate_idw(self,coord,registers,ruta,p=1,umbral=0.0):
		'Descripcion: Interpola la lluvia mediante la metodologia\n'\
		'	del inverso de la distancia ponderado. \n'\
		'\n'\
		'Parametros\n'\
		'----------\n'\
		'self : .\n'\
		'coord : Array (2,Ncoord) con las coordenadas de estaciones.\n'\
		'registers : Array (Nest,Nregisters) con los registros de lluvia.\n'\
		'p :  exponente para la interpolacion de lluvia.\n'\
		'ruta : Ruta con nombre en donde se guardara el binario con.\n'\
		'	la informacion de lluvia.\n'\
		'umbral : Umbral de suma total de lluvia bajo el cual se considera\n'\
		'	que un intervalo tiene suficiente agua como para generar reaccion\n'\
		'	(umbral = 0.0) a medida que incremente se generaran archivos mas\n'\
		'	livianos, igualmente existe la posibilidad de borrar informacion.\n'\
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
		status='update',umbral = 0.01):
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
		'	update: (Defecto) Con este estado se genera un binario y se agregan nuevos campos.\n'\
		'	old: Estado para abrir y tomar las propiedades de self.radar.. para la generacion de un binario.\n'\
		'	close: Cierra un binario que se ha generado mediante update.\n'\
		'	reset: Reinicia las condiciones de self.radar... para la creacion de un campo nuevo.\n'\
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
		if ruta_out <> None:
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
			if vec.mean() > umbral:
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
			f=open(ruta_hdr[:-3]+'.hdr','r')
			Lista = f.readlines()
			self.radarCont = int(Lista[3].split()[-1])
			f.close()
			#Abre con numpy para simplificar las cosas 
			a = np.loadtxt(ruta_hdr,skiprows=6,dtype='str').T
			self.radarPos = [int(i.split(',')[0]) for i in a[1]]
			self.radarMeanRain = [float(i.split(',')[0]) for i in a[2]]
			for i in a[3]:
				d = datetime.datetime.strptime(i,'%Y-%m-%d-%H:%M')
				self.radarDates.append(d)
		return actualizo 
		
	#------------------------------------------------------
	# Subrutinas para preparar modelo 
	#------------------------------------------------------	
	def set_Geomorphology(self,umbrales=[30,500],stream_width=None):
		'Descripcion: calcula las propiedades geomorfologicas necesarias \n'\
		'	para la simulacion. \n'\
		'\n'\
		'Parametros\n'\
		'----------\n'\
		'self : Inicia las variables vacias.\n'\
		'umbrales : Lista con la cantidad de celdas necesarias .\n'\
		'	para que una celda sea: ladera, carcava o cauce .\n'\
		'stream_width = Ancho del canal en cada tramo (opcional).\n'\
		'\n'\
		'Retornos\n'\
		'----------\n'\
		'self : Con las variables geomorfologicas de simulacion iniciadas.\n'\
		'	models.drena : Numero de celda o ladera destino. \n'\
		'	models.nceldas : Numero de celdas o laderas. \n'\
		'	models.unit_type : tipo de celda, en el caso de ladera no aplica.\n'\
		'		1: Celda tipo ladera.\n'\
		'		2: Celda tipo carcava.\n'\
		'		3: Celda tipo cauce.\n'\
		'	models.hill_long : Longitud de la ladera (o celda). \n'\
		'	models.hill_slope : Pendiente de cada ladera (o celda).\n'\
		'	models.stream_long : Longitud de cada tramo de cuace. \n'\
		'	models.stream_slope : Pendiente de cada tramo de cauce. \n'\
		'	models.stream_width : Ancho de cada tramo de cauce. \n'\
		'	models.elem_area : Area de cada celda o ladera. \n'\
		#Obtiene lo basico para luego pasar argumentos
		acum,hill_long,pend,elev = cu.basin_basics(self.structure,
			self.DEM,self.DIR,cu.ncols,cu.nrows,self.ncells)
		#Obtiene parametros para la cuenca como si fuera celdas
		#Obtiene el tipo de celdas
		unit_type = cu.basin_stream_type(self.structure,
			acum,umbrales,len(umbrales),self.ncells)
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
		if stream_width==None:
			stream_width=np.ones(self.ncells)
		#De acuerdo a si el modelo es por laderas o por celdas agrega lass varaibeles 
		if self.modelType[0]=='c':
			models.drena = np.ones((1,self.ncells))*self.structure
			models.nceldas = self.ncells
			models.unit_type = np.ones((1,self.ncells))*unit_type
			models.hill_long = np.ones((1,self.ncells))*hill_long
			models.hill_slope = np.ones((1,self.ncells))*pend
			models.stream_long = np.ones((1,self.ncells))*hill_long
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
			models.stream_width = np.ones((1,N))*cu.basin_subbasin_map2subbasin(
				self.hills_own,stream_width,self.nhills,0,self.ncells,self.CellCauce)
			models.elem_area = np.ones((1,N))*np.array([self.hills_own[self.hills_own==i].shape[0] for i in range(1,self.hills.shape[1]+1)])*cu.dxp**2.0			
	def set_Speed_type(self,types=np.ones(3)):
		'Descripcion: Especifica el tipo de velocidad a usar en cada \n'\
		'	nivel del modelo. \n'\
		'\n'\
		'Parametros\n'\
		'----------\n'\
		'self : Inicia las variables vacias.\n'\
		'types : tipos de velocidad .\n'\
		'	1. Velocidad tipo embalse lineal, no se especifica h_exp.\n'\
		'	2. Velocidad onda cinematica, se debe especificar h_exp.\n'\
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
	def set_PhysicVariables(self,modelVarName,var,pos,mask=None):
		'Descripcion: Coloca las variables fisicas en el modelo \n'\
		'	Se debe assignarel nombre del tipo de variable, la variable\n'\
		'	y la posicion en que esta va a ser insertada\n'\
		'\n'\
		'Parametros\n'\
		'----------\n'\
		'self : Objeto de la cuenca 2que luego se va a simular.\n'\
		'modelVarName : Nombre de la variable del modelo que se va a incertar.\n'\
		'	- h_coef.\n'\
		'	- h_exp.\n'\
		'	- v_coef.\n'\
		'	- v_exp.\n'\
		'	- capilar.\n'\
		'	- gravit.\n'\
		'var : variable que ingresa en el modelo, esta puede ser:.\n'\
		'	- Ruta: una ruta del tipo string.\n'\
		'	- Escalar : Un valor escalar que se asignara a toda la cuenca.\n'\
		'	- Vector : Un vector con la informacion leida (1,ncells).\n'\
		'pos : Posicion de insercion, aplica para : h_coef, v_coef,.\n'\
		'	h_exp, v_exp.\n'\
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
	def set_Storage(self,var,pos):
		'Descripcion: \n'\
		'	Establece el almacenamiento inicial del modelo\n'\
		'	la variable puede ser un valor, una ruta o un vector.\n'\
		'\n'\
		'Parametros\n'\
		'----------\n'\
		'var : Variable conla cual se va a iniciar el almancenamiento.\n'\
		'	- ruta : es una ruta a un archivo binario de almacenamiento.\n'\
		'	- escalar : valor de almacenamiento constate para toda la cuenca.\n'\
		'	- vector : Vector con valores para cadda unidad de la cuenca.\n'\
		'pos : Posicion de insercion,.\n'\
		'	- 0 :  alm cpailar.\n'\
		'	- 1 :  alm superficial.\n'\
		'	- 2 :  alm sub-superficial.\n'\
		'	- 3 :  alm subterraneo.\n'\
		'	- 4 :  alm cauce.\n'\
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
			#Si es un string lee el binario de almacenamiento alojado en esa ruta 
			Vec,res = models.read_float_basin(var,pos+1,N)
			isVec=True
		elif type(var) is int or float:
			Vec = np.ones((1,N))*var
			isVec=True
		elif type(var) is np.ndarray and var.shape[0] == N:			
			Vec = var
			isVec=True
		#Aloja ese almacenamiento en la cuenca 
		models.storage[pos] = Vec
	def set_Control(self,coordXY,ids,tipo = 'Q'):
		'Descripcion: \n'\
		'	Establece los puntos deonde se va a realizar control del caudal\n'\
		'	y de la humedad simulada.\n'\
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
		'	- models.control : control de caudal y sedimentos.\n'\
		'	- models.control_h : control de la humedad del suelo.\n'\
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
				IdsConvert = basinPts[basinPts<>0]
			elif self.modelType[0] is 'h':
				unitario = basinPts / basinPts
				pos = self.hills_own * self.CellCauce * unitario
				posGrande = self.hills_own * self.CellCauce * basinPts
				IdsConvert = posGrande[posGrande<>0] / pos[pos<>0]
				models.control[0][pos[pos<>0].astype(int).tolist()] = IdsConvert 			
		elif tipo is 'H':
			xyNew = coordXY
			basinPts, order = self.Points_Points2Basin(coordXY,ids)
			if self.modelType is 'cells':
				models.control_h[0] = basinPts
				IdsConvert = basinPts[basinPts<>0]
			elif self.modelType is 'hills':
				unitario = basinPts / basinPts
				pos = self.hills_own * unitario
				posGrande = self.hills_own * basinPts
				IdsConvert = posGrande[posGrande<>0] / pos[pos<>0]
				models.control_h[0][pos[pos<>0].astype(int).tolist()] = IdsConvert 
		return IdsConvert,xyNew
	
	
	#def set_sediments(self,var,varName):
		
		
	def set_slides(self,var,VarName):
		'Descripcion: Alojas las variables requeridas para la ejecucion\n'\
		'	del modelo de sedimentos.\n'\
		'\n'\
		'Parametros\n'\
		'----------\n'\
		'var : Variable que describe la propiedad (constante, ruta, mapa o vector).\n'\
		'varName: Nombre de la variable a ingresar en el modelo.\n'\
		'	Zs : Profundidad del suelo.\n'\
		'	GammaSoil: Densidad del suelo.\n'\
		'	Cohesion: Cohesion del suelo.\n'\
		'	FrictionAngle : Angulo de friccion del suelo.\n'\
		'	FS: Factor de Seguridad, en este caso se envia una constante.\n'\
		'	RadSlope : .\n'\
		'\n'\
		'Retornos\n'\
		'----------\n'\
		'self : variables iniciadas en el modelo bajo los nombres de:.\n'\
		'	wmf.models.gammas.\n'\
		'	wmf.models.cohesion.\n'\
		'	wmf.models.frictionangle.\n'\
		'	wmf.models.zs.\n'\
		#Obtiene el vector que va a alojar en el modelo
		if VarName <> 'FS':
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
			models.gammas = np.ones((1,N))*Vec
		elif VarName is 'Cohesion':
			models.cohesion = np.ones((1,N))*Vec
		elif VarName is 'FrictionAngle':
			models.frictionangle = np.ones((1,N))*Vec
		elif VarName is 'Zs':
			models.zs = np.ones((1,N))*Vec
		elif VarName is 'FS':
			models.fs = var
		elif Varname is 'Slope':
			models.radslope = np.ones((1,N))*Vec
	#------------------------------------------------------
	# Guardado y Cargado de modelos de cuencas preparados 
	#------------------------------------------------------	
	def Save_SimuBasin(self,ruta,ruta_dem = None,ruta_dir = None):
		'Descripcion: guarda una cuenca previamente ejecutada\n'\
		'\n'\
		'Parametros\n'\
		'----------\n'\
		'ruta : Ruta donde la cuenca sera guardada.\n'\
		'ruta_dem : direccion donde se aloja el DEM (se recomienda absoluta).\n'\
		'ruta_dir : direccion donde se aloja el DIR (se recomienda absoluta).\n'\
		'\n'\
		'Retornos\n'\
		'----------\n'\
		'self : Con las variables iniciadas.\n'\
		#Guarda la cuenca
		if self.modelType[0] is 'c':
			N = self.ncells
		elif self.modelType[0] is 'h':
			N = self.nhills
		#Determina las rutas
		if ruta_dem is None:
			ruta_dem = 'not rute'
		if ruta_dir is None:
			ruta_dir = 'not rute'
		Dict = {'nombre':self.name,'DEM':ruta_dem,
			'DIR':ruta_dir,
		    'modelType':self.modelType,'noData':self.nodata,'umbral':self.umbral,
		    'ncells':self.ncells,'nhills':self.nhills,
		    'dt':models.dt,'Nelem':N,'dxp':cu.dxp,'retorno':models.retorno}
		#abre el archivo 
		gr = netcdf.Dataset(ruta,'w',format='NETCDF4')
		#Establece tamano de las variables 
		DimNcell = gr.createDimension('ncell',self.ncells)
		DimNhill = gr.createDimension('nhills',self.nhills)
		DimNelem = gr.createDimension('Nelem',N)
		DimCol3 = gr.createDimension('col3',3)
		DimCol2 = gr.createDimension('col2',2)
		DimCol4 = gr.createDimension('col4',4)
		DimCol5 = gr.createDimension('col5',5)		
		#Crea variables
		VarStruc = gr.createVariable('structure','i4',('col3','ncell'),zlib=True)
		VarHills = gr.createVariable('hills','i4',('col2','nhills'),zlib=True)
		VarHills_own = gr.createVariable('hills_own','i4',('ncell',),zlib=True)
		VarH_coef = gr.createVariable('h_coef','f4',('col4','Nelem'),zlib=True)
		VarV_coef = gr.createVariable('v_coef','f4',('col4','Nelem'),zlib=True)
		VarH_exp = gr.createVariable('h_exp','f4',('col4','Nelem'),zlib=True)
		VarV_exp = gr.createVariable('v_exp','f4',('col4','Nelem'),zlib=True)
		Var_H1max = gr.createVariable('h1_max','f4',('Nelem',),zlib = True)
		Var_H3max = gr.createVariable('h3_max','f4',('Nelem',),zlib = True)
		Control = gr.createVariable('control','i4',('Nelem',),zlib = True)
		ControlH = gr.createVariable('control_h','i4',('Nelem',),zlib = True)
		if self.modelType[0] is 'c':
			drena = gr.createVariable('drena','i4',('col3','Nelem'),zlib = True)
		elif self.modelType[0] is 'h':
			drena = gr.createVariable('drena','i4',('Nelem'),zlib = True)
		unitType = gr.createVariable('unit_type','i4',('Nelem',),zlib = True)
		hill_long = gr.createVariable('hill_long','f4',('Nelem',),zlib = True)
		hill_slope = gr.createVariable('hill_slope','f4',('Nelem',),zlib = True)
		stream_long = gr.createVariable('stream_long','f4',('Nelem',),zlib = True)
		stream_slope = gr.createVariable('stream_slope','f4',('Nelem',),zlib = True)
		stream_width = gr.createVariable('stream_width','f4',('Nelem',),zlib = True)
		elem_area = gr.createVariable('elem_area','f4',('Nelem',),zlib = True)
		speed_type = gr.createVariable('speed_type','i4',('col3',),zlib = True)
		storage = gr.createVariable('storage','i4',('col5','Nelem'),zlib = True)
		
		#Asigna valores a las variables
		VarStruc[:] = self.structure
		VarHills[:] = self.hills
		VarHills_own[:] = self.hills_own
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
		
		#asigna las prop a la cuenca 
		gr.setncatts(Dict)
		#Cierra el archivo 
		gr.close()
		#Sale del programa
		return 

	#------------------------------------------------------
	# Ejecucion del modelo
	#------------------------------------------------------	
	def run_shia(self,Calibracion,
		rain_rute, N_intervals, start_point = 1, ruta_storage = None,
		ruta_conv = None, ruta_stra = None):
		'Descripcion: Ejecuta el modelo una ves este es preparado\n'\
		'	Antes de su ejecucion se deben tener listas todas las . \n'\
		'	variables requeridas . \n'\
		'\n'\
		'Parametros\n'\
		'----------\n'\
		'self : Cuenca a ejecutar con todo listo para ser ejecutada.\n'\
		'Calibracion : Parametros de calibracion del modelo, orden:.\n'\
		'	- Evaporacion.\n'\
		'	- Infiltracion.\n'\
		'	- Percolacion.\n'\
		'	- Perdidas.\n'\
		'	- Vel Superficial .\n'\
		'	- Vel Sub-superficial.\n'\
		'	- Vel Subterranea.\n'\
		'	- Vel Cauce.\n'\
		'	- Max Capilar.\n'\
		'	- Max Gravitacional.\n'\
		'rain_rute : Ruta donde se encuentra el archivo binario de lluvia:.\n'\
		'	generado por rain_interpolate_* o por rain_radar2basin.\n'\
		'N_intervals : Numero de intervalos de tiempo.\n'\
		'start_point : Punto donde comienza a usar registros de lluvia.\n'\
		'	los binarios generados por rain_* generan un archivo de texto.\n'\
		'	que contiene fechas par aayudar a ubicar el punto de inicio deseado.\n'\
		'ruta_storage : Ruta donde se guardan los estados del modelo en cada intervalo.\n'\
		'	de tiempo, esta es opcional, solo se guardan si esta variable es asignada.\n'\
		'ruta_conv : Ruta al binario y hdr indicando las nubes que son convectivas.\n'\
		'ruta_stra : Ruta al binario y hdr indicando las nubes que son estratiformes.\n'\
		'\n'\
		'Retornos\n'\
		'----------\n'\
		'Qsim : Caudal simulado en los puntos de control.\n'\
		'Hsim : Humedad simulada en los puntos de control.\n'\
		#genera las rutas 
		rain_ruteBin,rain_ruteHdr = __Add_hdr_bin_2route__(rain_rute)
		# De acuerdo al tipo de modelo determina la cantidad de elementos
		if self.modelType[0] is 'c':
			N = self.ncells
		elif self.modelType[0] is 'h':
			N = self.nhills
		#prepara variables globales
		models.rain_first_point = start_point
		#Prepara terminos para control
		if np.count_nonzero(models.control) == 0 :
			NcontrolQ = 1
		else:
			NcontrolQ = np.count_nonzero(models.control)+1
		if np.count_nonzero(models.control_h) == 0 :
			NcontrolH = 1
		else:
			NcontrolH = np.count_nonzero(models.control_h)
		#Prepara variables de guardado de variables 
		if ruta_storage is not None:
			models.save_storage = 1
		else:
			ruta_storage = 'no_guardo_nada.bin'
		#Variables de separacion de flujo por tipo de lluvia 
		if models.separate_rain == 1 and ruta_conv <> None and ruta_stra <> None:
			ruta_binConv,ruta_hdrConv = __Add_hdr_bin_2route__(ruta_conv)
			ruta_binStra,ruta_hdrStra = __Add_hdr_bin_2route__(ruta_stra)
		else:
			models.separate_rain = 0
			ruta_binConv = 'none'
			ruta_hdrConv = 'none'
			ruta_binStra = 'none'
			ruta_hdrStra = 'none'
			
		# Ejecuta el modelo 
		Qsim,Qsed,Qseparated,Humedad,Balance,Alm,Qsep_byrain = models.shia_v1(
			rain_ruteBin,
			rain_ruteHdr,
			Calibracion,
			N,
			NcontrolQ,
			NcontrolH,
			N_intervals,
			ruta_storage,
			ruta_binConv,
			ruta_binStra,
			ruta_hdrConv,
			ruta_hdrStra)
		#Retorno de variables de acuerdo a lo simulado 
		Retornos={'Qsim' : Qsim}
		Retornos.update({'Balance' : Balance})
		Retornos.update({'Storage' : Alm})
		if np.count_nonzero(models.control_h)>0:
			Retornos.update({'Humedad' : Humedad})
		if models.sim_sediments == 1:
			Retornos.update({'Sediments' : Qsed})
		if models.separate_fluxes == 1:
			Retornos.update({'Fluxes' : Qseparated})
		if models.separate_rain == 1:
			Retornos.update({'Rain_sep' : Qsep_byrain})
		if models.show_storage == 1:
			Retornos.update({'Mean_Storage' : np.copy(models.mean_storage)})
		return Retornos
		
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
		if ruta<>None:
			pl.savefig(ruta,bbox_inches='tight')
		pl.show()
	#def Plot_Map(self,ruta=None):
