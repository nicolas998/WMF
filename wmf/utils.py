#!/usr/bin/env python
#---------------------------------------------------------------------------------------------------
#Descripcion
#---------------------------------------------------------------------------------------------------
# utils.py: desarrollado como un modulo orientado a guardar datos en el formato de la cuenca 
# trazada, este es usado desde python
#utils.py: Modulo para el manejo de datos propios de cuencas.f90 y modelos.f90
#Copyright (C) <2014>  <Nicolas Velasquez Giron>

#This program is free software: you can redistribute it and/or modify
#it under the terms of the GNU General Public License as published by
#the Free Software Foundation, either version 3 of the License, or
#(at your option) any later version.

#This program is distributed in the hope that it will be useful,
#but WITHOUT ANY WARRANTY; without even the implied warranty of
#MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#GNU General Public License for more details.

#You should have received a copy of the GNU General Public License
#along with this program.  If not, see <http://www.gnu.org/licenses/>.

#---------------------------------------------------------------------------------------------------
#Modulos a usar
#---------------------------------------------------------------------------------------------------
import numpy as np
import sys 
import osgeo.ogr, osgeo.osr
import os
from modelacion import cuencas as cu
from modelacion import lluvia as ll
from scipy.spatial import Delaunay
import pylab as pl
import gdal
from osgeo import ogr
from matplotlib.colors import LogNorm, LightSource
from mpl_toolkits.basemap import Basemap, addcyclic, shiftgrid, cm
from osgeo import ogr,osr
import pandas 
import scipy.stats as stat
import osgeo 

#---------------------------------------------------------------------------------------------------
#Guardado y lectura de cuenca
#---------------------------------------------------------------------------------------------------
#Funciones de guardado
def write_proyect_int(basin,nceldas,ruta,proyecto=None,ncols=None,nrows=None,xll=None,yll=None,noData=None,dx=None,dxp=None,nombre=None):
    #Guarda el archivo base de la cuenca 
    'Funcion: write_proyect_int\n'\
    'Descripcion: Guarda la cuenca trazada como un proyecto compuesto por un archivo header de texto plano\n'\
    'y un archivo base binario que contiene lo basico de la cuenca trazada. En el header se indican las \n'\
    'propiedades del mapa base de la cuenca, y la cantidad de celdas que la componen, en el binario se \n'\
    'almacenan en la primera columna el punto al cual drena cada celda y en la 2da y 3ra columnas la fila y columna\n'\
    'Parametros Obligatorios:.\n'\
    '	-basin: vector obtenido con cu.basin_cut, contiene 1: drena id, 2: col, 3: fil.\n'\
    '	-nceldas: Cantidad de filas que componen el vector basin.\n'\
    '	-ruta: Ruta en la cual se escribe el proyecto de la cuenca.\n'\
    'Parametros Opcionales:.\n'\
    '	-proyecto: Nombre de la clase cuencas usada para el trabajo (ej: from fmod import cuencas as cu1).\n'\
    '		Si se indica esta clase no es necesario indicar los demas argumentos.\n'\
    '	-ncols: Numero de columnas del mapa base.\n'\
    '	-nrows: Numero de filas del mapa base.\n'\
    '	-xll: Coordenada en X de la esquina inferior izquierda de la celda inferior izquierda.\n'\
    '	-yll: Coordenada en Y de la esquina inferior izquierda de la celda inferior izquierda.\n'\
    '	-noData: Numero indicador de no data.\n'\
    '	-dx: Longitud de un lado de la celda en las coordenadas originales del mapa.\n'\
    '	-dxp: .\n'\
    'Retorno:.\n'\
    '	-ruta: ruta en la cual se escribiran el archivo binario .int y su respectivo encabezado .int_hdr \n'\
     'Lista Nombres reconocidos por el modelo:.\n'\
    '	-Tipo Cauce: Tipo de celda, para el modelo.\n'\
    '	-Pts Control: Vector con las celdas diferentes a cero donde hay puntos de control.\n'\
    #Evalua si el archivo ya existe o no
    if os.path.exists(ruta) and basin.ndim==1 and nombre<>None:
	#Si ya existe y basin es el vector de la cuenca, y se asigno nombre, escribe la nueva entrada
	write_proyect_int_ext(basin,nombre,nceldas,ruta)
    else:
	#Si ya existe pero basin es de 3 columnas, o no se da el nombre, sobre escribe el archivo de la cuenca
	#Define el nombre de la ruta
	if ruta.endswith('.int'):
	    ruta_hdr=ruta+'_hdr'
	else:
	    ruta=ruta+'.int'
	    ruta_hdr=ruta+'_hdr'
	#define las propieades
	if proyecto<>None:
	    ncols=proyecto.ncols; nrows=proyecto.nrows; xll=proyecto.xll; yll=proyecto.yll
	    dx=proyecto.dx; dxp=proyecto.dxp; noData=proyecto.nodata
	#Guarda el header
	f=open(ruta_hdr,'w')
	f.write("ncols:           %i \n" % ncols)
	f.write("nrows:           %i \n" % nrows)
	f.write("xllcorner:       %f \n" % xll)
	f.write("yllcorner:       %f \n" % yll)
	f.write("cellsize:        %f \n" % dx)
	f.write("cellsize_plana:  %f \n" % dxp)
	f.write("noData_Value:    %i \n" % noData)
	f.write("Numero Columnas: %i \n" % 3)
	f.write("Numero de filas: %i \n" % nceldas)
	f.write("1: DrenaID \n")
	f.write("2: Columnas \n")
	f.write("3: Filas \n")
	f.close()
	#Guarda el binario
	cu.write_int_basin(ruta,basin[0,:],1,'replace',nceldas)
	for i in range(1,3):
	    cu.write_int_basin(ruta,basin[i,:],i+1,'old',nceldas)
def write_proyect_int_ext(vect,nombre,nceldas,ruta):
    #Guarda vectores extras en el archivo de la cuenca
    'Funcion: write_proyect_int_ext\n'\
    'Descripcion: Ingresa nuevas columnas en el proyecto .int a partir de vectores de enteros \n'\
    'Parametros Obligatorios:.\n'\
    '	-vect: Vector con la variable entera a ingresar, debe tener la estructura de la cuenca.\n'\
    '	-nombre: nombre con el que sera reconocida la variable .\n'\
    '	-nceldas: Cantidad de filas que componen el vector basin.\n'\
    '	-ruta: Ruta en la cual se escribe el proyecto de la cuenca.\n'\
    'Retorno:.\n'\
    '	Actualizacion del binario .int\n'\
    #Define nombre de la ruta
    if ruta.endswith('.int'):
	ruta_hdr=ruta+'_hdr'
    else:
	sys.exit('Error: Archivo de cuenca no reconocido')
    #Abre header para conocer la cantidad de columnas
    f=open(ruta_hdr,'r')
    Linea=f.readlines()
    N_entradas=int(Linea[7].split()[-1])
    N_filas=int(Linea[8].split()[-1])
    f.close()
    #Verifica que el vector nuevo si tenga el tamano de la cuenca
    if N_filas<>nceldas:
	sys.exit('Error: Vector tiene %i filas, cuenca tiene %i filas' % (nceldas,N_filas))
    #Verifica que la variable a agregar si corresponda al tipo de archivo .int o .float
    if ruta.endswith('.int'):
	terminacion=ruta[-3:]
    if str(vect.dtype)[:-2]<>terminacion:
	sys.exit('Error: Vector a agregar tipo %s, archivo contenedor de variables tipo %s' % (str(vect.dtype)[:-2],terminacion))			
    #si pasa las dos verificaciones abre header para actualizar
    f=open(ruta_hdr,'r')
    Lineas=f.readlines()
    f.close()
    #Busca linea por linea si esa variable ya ha sido ingresada
    N_col=-9999
    N_col=read_proyect_identify(ruta_hdr,nombre)
    #Depende de si esta o no actualiza la misma linea o crea una nueva
    if N_col>0:
	#Ya existia solo sobre-escribe
	cu.write_int_basin(ruta,vect,N_col,'old',nceldas)
    else:    	
	#No existia, le toca escribir nueva	
	Lineas[7]='Numero Columnas: '+str(N_entradas+1)+' \n'	
	f=open(ruta_hdr,'w')
	f.writelines(Lineas)
	f.write("%i: %s \n" % (N_entradas+1,str(nombre)))
	f.close()
	#Escribe sobre el binario
	cu.write_int_basin(ruta,vect,N_entradas+1,'old',nceldas)	
def write_proyect_float(vect,nombre,nceldas,ruta):
    'Funcion: write_proyect_float\n'\
    'Descripcion: Toma como base el proyecto .int escrito por la funcion write_proyect_in, y escribe \n'\
    'sobre este valores de mapas reales. \n'\
    'Parametros Obligatorios:.\n'\
    '	-vect: Vector con la variable real a ingresar, debe tener la estructura de la cuenca.\n'\
    '	-nombre: nombre con el que sera reconocida la variable .\n'\
    '	-nceldas: Cantidad de filas que componen el vector basin.\n'\
    '	-ruta: Ruta en la cual se escribe el proyecto de la cuenca, debe tener el nombre del proyecto\n'\
    '		asociado entero "nombre_proy.int", pero sin la extension o bien con la extension ".float".\n'\
    'Retorno:.\n'\
    '	Escribe el archibo binario flotante con la variable ingresada\n'\
    'Lista Nombres reconocidos por el modelo:.\n'\
    '	-Longitud: Longitud de cada una de las celdas de la cuenca.\n'\
    '	-Pendiente Laderas: Pendiente de cada una de las celdas sin ser consideradas cauce.\n'\
    '	-Pendiente Cauce: Pendiente de las celdas pertenecientes a cada uno de los tramos de la red.\n'\
    '	-hu: Almacenamiento maximo capilar.\n'\
    '	-ks: Conductividad saturada sub-superficial.\n'\
    '	-kp: Conductividad saturada subterranea.\n'\
    '	-N: Valor de rugosidad de Manning de la superficie.\n'\
    '	-evp: Mapa de evaporacion potencial.\n'\
    '	-Coord X: Coordenada X de cada una de las celdas.\n'\
    '	-Coord Y: Coordenada Y de cada una de las celdas.\n'\
    #chequea si termina en float si no lo pone
    if ruta.endswith('.float')==False:
	ruta=ruta+'.float'
    else:
	sys.exit('Error: El archivo debe terminar en .float')
    ruta_hdr=ruta+'_hdr'
    #chequea si el archivo float ya existe
    if os.path.exists(ruta):
	write_proyect_float_ext(vect,nombre,nceldas,ruta)
    #si no existe crea uno nuevo
    else:
	#Verifica que la variable a agregar si corresponda al tipo de archivo .int o .float
	terminacion='float'
	if str(vect.dtype)[:-2]<>terminacion:
	    sys.exit('Error: Vector a agregar tipo %s, archivo contenedor de variables tipo %s' % (str(vect.dtype)[:-2],terminacion))			
	#si pasa las dos verificaciones abre header por primera ves
	f=open(ruta_hdr,'w')
	f.write("Numero Columnas: %i \n" % 1)
	f.write("Numero de filas: %i \n" % nceldas)
	f.write("1: %s \n" % (str(nombre)))
	f.close()
	#Graba en el binario
	cu.write_float_basin(ruta,vect,1,'replace',nceldas)		
def write_proyect_float_ext(vect,nombre,nceldas,ruta):
    #Guarda vectores extras en el archivo de la cuenca
    'Funcion: write_proyect_float_ext\n'\
    'Descripcion: Ingresa nuevas columnas en el proyecto .float a partir de vectores de reales \n'\
    'Parametros Obligatorios:.\n'\
    '	-vect: Vector con la variable entera a ingresar, debe tener la estructura de la cuenca.\n'\
    '	-nombre: nombre con el que sera reconocida la variable .\n'\
    '	-nceldas: Cantidad de filas que componen el vector basin.\n'\
    '	-ruta: Ruta en la cual se escribe el proyecto de la cuenca.\n'\
    'Retorno:.\n'\
    '	Actualizacion del binario .float\n'\
    #Define nombre de la ruta
    if ruta.endswith('.float'):
	ruta_hdr=ruta+'_hdr'
    else:
	sys.exit('Error: Archivo de cuenca no reconocido')
    #Abre header para conocer la cantidad de columnas
    f=open(ruta_hdr,'r')
    Linea=f.readlines()
    N_entradas=int(Linea[0].split()[-1])
    N_filas=int(Linea[1].split()[-1])
    f.close()
    #Verifica que el vector nuevo si tenga el tamano de la cuenca
    if N_filas<>nceldas:
	sys.exit('Error: Vector tiene %i filas, cuenca tiene %i filas' % (nceldas,N_filas))
    #Verifica que la variable a agregar si corresponda al tipo de archivo .int o .float
    if ruta.endswith('.int'):
	terminacion=ruta[-3:]
    elif ruta.endswith('.float'):
	terminacion=ruta[-5:]
    if str(vect.dtype)[:-2]<>terminacion:
	sys.exit('Error: Vector a agregar tipo %s, archivo contenedor de variables tipo %s' % (str(vect.dtype)[:-2],terminacion))			
    #si pasa las dos verificaciones abre header para actualizar
    f=open(ruta_hdr,'r')
    Lineas=f.readlines()
    f.close()
    #Busca linea por linea si esa variable ya ha sido ingresada
    N_col=-9999
    N_col=read_proyect_identify(ruta_hdr,nombre)
    #Depende de si esta o no actualiza la misma linea o crea una nueva
    if N_col>0:
	#Ya existia solo sobre-escribe
	cu.write_float_basin(ruta,vect,N_col,'old',nceldas)
    else:    	
	#No existia, le toca escribir nueva	
	Lineas[0]='Numero Columnas: '+str(N_entradas+1)+' \n'	
	f=open(ruta_hdr,'w')
	f.writelines(Lineas)
	f.write("%i: %s \n" % (N_entradas+1,str(nombre)))
	f.close()
	#Escribe sobre el binario
	cu.write_float_basin(ruta,vect,N_entradas+1,'old',nceldas)		
def write_proyect_storage(vect,nceldas,ruta,associated=None,calibracion=None,model='SHIA',nombre=None,columna=None):
    #Guarda las condiciones iniciales de un modelo
    'Funcion: write_proyect_storage\n'\
    'Descripcion: Crea un archivo binario y un header con el almacenamiento del modelo seleccionado\n'\
    'Parametros Obligatorios:.\n'\
    '	-vect: Vector con las entradas de cada uno de los almacenamientos del modelo, debe tener.\n'\
    '	tantas entradas como almacenamientos tenga el modelo seleccionado.\n'\
    '	-nceldas: Cantidad de filas que componen el vector basin.\n'\
    '	-ruta: Ruta en la cual se escribe el proyecto de la cuenca.\n'\
    'Parametros Opcionales:.\n'\
    '	-associated: Ruta del proyecto .int asociado al archivo de almacenamiento a crear.\n'\
    '	-nombre: Lista con los nombres de cada uno de los almacenamientos, si el modelo esta en la.\n'\
    '	lista se puede omitir esta accion y el selecciona los nombres de acuerdo al modelo.\n'\
    '	-calibracion: Parametros de calibracion utilizados para la corrida del modelo.\n'\
    '	-model: Modelo para el cual se estan guardando los estados de almacenamiento (Defecto: SHIA).\n'\
    'Retorno:.\n'\
    '	Escribe un archivo binario .store con los almacenamientos, y un header .store_hdr.\n'\
    #Chequea si ya existe o no el archivo de almacenamiento
    if os.path.exists(ruta):
	#si ya existe actualiza con el nombre o con la columna
	if nombre<>None:
	    write_proyect_storage_ext(vect,nceldas,ruta,nombre=nombre)
	elif columna<>None:
	    write_proyect_storage_ext(vect,nceldas,ruta,columna=columna)
    #Si no existe crea un archivo nuevo
    else:
	#Chequea si la ruta termina con la extension
	if ruta.endswith('.store')==False:
	    ruta=ruta+'.store'
	#Chequea que el archivo asociado tenga la extension correcta
	if associated.endswith('.int_hdr')==True or associated.endswith('.int')==True:
	    #si el associado no termina en .int_hdr le agrega eso
	    if associated.endswith('.int_hdr')==False and associated.endswith('.int')==True:
		associated=associated+'_hdr'
	    #Chequea que la cantidad de celdas sea igual a la que esta en la tabla
	    nceldas_file=read_proyect_hdr(associated)
	    nceldas_file=int(nceldas_file['Numero de filas'])
	    if nceldas<>nceldas_file:
		sys.exit('Error: Cantidad de filas diferente a la cantidad de filas del proyecto seleccionado')
	else:
	    sys.exit('Error: Archivo asociado no termina con .int_hdr o .int')	
	#Depende del tipo de modelo crea los nombres de las entradas 
	if model=='SHIA':
	    names=['S1','S2','S3','S4','S5']
	    #Si es el modelo shia el vector debe tener 5 entradas
	    if vect.shape[0]<>5:
		sys.exit('Error: Cantidad de almacenamientos diferente de 5')
	else:
	    sys.exit('Error: Nombre de modelo no implementado')
	#Crea el archivo header y escribe en el 
	f=open(ruta+'_hdr','w')
	f.write("Archivo asociado: %s \n" % os.path.abspath(associated))
	f.write("Modelo: %s \n" % model)
	f.write("Parametros Calibracion: ")
	if calibracion<>None:
	    for i in calibracion:
		f.write("%.3f, " % i)
	f.write("\n")
	f.write("Numero de columnas: %d \n" % vect.shape[0])
	f.write("Numero de filas: %d \n" % nceldas)
	cont=1
	for i in names:
	    f.write("%d: %s \n" % (cont,i))
	    cont+=1
	f.close()
	#Crea el binario correspondiente
	cu.write_float_basin(ruta,vect[0],1,'replace',nceldas)
	cont=2
	for i in vect[1:]:
	    cu.write_float_basin(ruta,i,cont,'old',nceldas)
	    cont+=1
def write_proyect_storage_ext(vect,nceldas,ruta,nombre=None,columna=None):
    #Guarda o modifica una de las entradas de almacenamiento de los archivos .storage
    'Funcion: write_proyect_storage_ext\n'\
    'Descripcion: Actualiza el binario y el header del almacenamiento del modelo\n'\
    'Parametros Obligatorios:.\n'\
    '	-vect: Vector con las entradas del almacenamiento a modificar.\n'\
    '	-nceldas: Cantidad de filas que componen el vector basin.\n'\
    '	-ruta: Ruta en la cual se escribe el proyecto de la cuenca.\n'\
    'Parametros Opcionales:.\n'\
    '	-nombre: Nombre de la variable a actualizar.\n'\
    '	-columna: numero de la columna a actualizar. \n'\
    'Retorno:.\n'\
    '	Actualiza el archivo binario .store con los almacenamientos, y el header .store_hdr.\n'\
    #Verifica que el nombre de la ruta si sea coherente
    if ruta.endswith('.store_hdr')==False and ruta.endswith('.store')==False:
	sys.exit('Error: Nombre de archivo no corresponde con un archivo de estados iniciales')
    else:
	#si la ruta no termina con el hdr, se lo agrega
	if ruta.endswith('.store'):
	    ruta=ruta+'_hdr'
    #Lee la cantidad de columnas que hay en el archivo 
    f=open(ruta,r)
    Lineas=f.readlines()
    f.close()
    ncols=int(Lineas[2].split()[-1])
    asociado=Lineas[0].split(':')[:-1]
    #Verifica que se indique un nombre o la columna para especificar cual entrada va a ser editada
    if nombre<>None and columna==None:
	col=read_proyect_identify(ruta,nombre)
    elif columna<>None and nombre==None:
	if columna>ncols or columna<0:
	    sys.exit('Error: Columna especificada mayor a la cantidad de columnas o menor a cero')
	else:
	    col=columna
    else:
	sys.exit('Error: Se debe especificar o nombre o columna')
    #Verifica que el tamano del vector si cumpla 
    nceldas_asociado=read_proyect_hdr(asociado)
    nceldas_asociado=nceldas_asociado['Numero de filas']
    if nceldas<>nceldas_asociado:
	sys.exit('Error: Numero de filas del vector diferente al numero de filas del proyecto')
    #Si todo lo anterior se cumplio ingresa el vector en la nueva entrada
    cu.write_float_basin(ruta[:-4],vect,col,'old',nceldas)
def write_proyect_rain(x,y,proy_rain,nceldas,ruta,pp=1,date='0000-00-00-00-00',dt='x [min]',tipo='idw',associated=None):
    #Guarda las condiciones iniciales de un modelo
    'Funcion: write_proyect_rain\n'\
    'Descripcion: Crea un archivo binario y header de los campos de lluvia en el tiempo de acuerdo a la\n'\
    'topologia de la cuenca. \n'\
    'Parametros Obligatorios:.\n'\
    '	-x,y: vectores con las coordenadas X y Y de las celdas (x,y=cu.basin_coordxy(basin,nceldas)).\n'\
    '	-proy_rain: Nombre del modulo de lluvia importado (from fwm.modelacion import lluvia).\n'\
    '		Debe tener cargadas las variables: punt_lluvia, punt_coord.\n'\
    '	-nceldas: Cantidad de filas que componen el vector basin.\n'\
    '	-ruta: Ruta en la cual se escribe el campo de lluvia de la cuenca.\n'\
    'Parametros Opcionales:.\n'\
    '	-pp: Exponente del calculo de pesos (1/d(x0,x1)**pp), defecto PP=1.\n'\
    '	-date: Fecha y hora del comienzo del evento, YYYY-MM-DD-HH-mm.\n'\
    '	-dt: Paso del tiempo entre intervalos Val [Unidades].\n'\
    '	-tipo: tipo de campo: idw: Inverso de la distancia, tin: triangulos irregulares, rad: radar.\n'\
    '	-associated: Archivo de cuenca asociado al campo de lluvia interpolado.\n'\
    'Retorno:.\n'\
    '	Escribe un archivo binario .rain con los almacenamientos, y un header .rain_hdr.\n'\
    #Evalua si la ruta si tiene la extension que es
    if ruta.endswith('rain')==False and ruta.endswith('rain_hdr')==False:
	ruta=ruta+'.rain'
	ruta_hdr=ruta+'_hdr'
    elif ruta.endswith('rain'):
	ruta_hdr=ruta+'_hdr'
    elif ruta.endswith('rain_hdr'):
	ruta_hdr=ruta
	ruta=ruta[:-4]
    N_reg=proy_rain.punt_lluvia.shape[1]
    N_est=proy_rain.punt_lluvia.shape[0]
    #Escribe el archivo header
    f=open(ruta_hdr,'w')
    f.write('nceldas:		%i \n' % nceldas)
    f.write('Numero Registros:	%i \n' % N_reg)
    if associated<>None:
	f.write('Asociado:		%s \n' % os.path.abspath(associated))
    else: 
	f.write('Asociado:		%s \n' % 'Null')
    f.write('Tipo de campo:		%s \n' % tipo)
    f.write('Fecha inicial:		%s \n' % date)
    f.write('Delta t:		%s \n' % dt)
    f.write('Exponente pp:		%s \n' % pp)
    f.close()
    #Escribe el archivo binario
    if tipo=='idw':
	r=proy_rain.interpolation_idw(x,y,N_reg,N_est,pp,nceldas,ruta)
	return r
    if tipo=='tin':
	if proy_rain.tin==None:
	    find_tin(proy_rain.punt_coord,proy_rain)
	    r,tin_perte=proy_rain.pre_tin_points_select(x,y,nceldas)
	r=proy_rain.interpolation_tin(x,y,tin_perte,proy_rain.num_reg,nceldas,ruta)
	return r
#funciones de lectura
def read_proyect_hdr(ruta,proyecto=None):
    'Funcion: read_proyect_hdr\n'\
    'Descripcion: Lee el encabezado de un proyecto de cuenca .int o .int_hdr.\n'\
    'Parametros Obligatorios:.\n'\
    '	-ruta: Ruta del proyecto .int o .int_hdr a leer.\n'\
    'Parametros Opcionales:.\n'\
    '	-proyecto: Si se especifica un proyecto importado como: from fwm.modealcion import cuencas as cu.\n'\
    '	como proyecto=cu, entonces cu copiara todas las propieades del header del proyecto.\n'\
    'Retorno:.\n'\
    '	Orden: Diccionario con la informacion del proyecto.\n'\
    '	proyecto: Actualizacion del modulo ingresado.\n'\
    #Determina si se trata del header si no lo hace y lo busca
    if ruta.endswith('.int'):
	ruta=ruta+'_hdr'
	f=open(ruta,'r')
    elif ruta.endswith('.int_hdr'):
	f=open(ruta,'r')
    else:		
	sys.exit("Error: No se encuentra el header del archivo base de la cuenca")
    #Determina cantidad de filas en la cuenca y cantidad de columnas
    Lista=f.readlines()
    Orden={}
    for i in Lista:
	try:    
	    nombre=int(i.split(':')[0])
	except:
	    nombre=i.split(':')[0]
	try:
	    val=float(i.split(':')[1])
	except:
	    val=i.split(':')[1]
	Orden.update({nombre:val})
    nceldas=int(Orden['Numero de filas'])
    #Si se asigna un proyecto de cuenca a trabajar le asigna las propiedades
    if proyecto<>None:
	proyecto.xll=Orden['xllcorner']
	proyecto.yll=Orden['yllcorner']
	proyecto.ncols=int(Orden['ncols'])
	proyecto.nrows=int(Orden['nrows'])
	proyecto.nodata=int(Orden['noData_Value'])
	proyecto.dx=Orden['cellsize']
	proyecto.dxp=Orden['cellsize_plana']
	return nceldas,Orden,proyecto
    return nceldas,Orden
def read_proyect_identify(ruta,name2find):
    #Funcion para identificar en donde esta algo dentro de un proyecto
    'Funcion: read_proyect_identify\n'\
    'Descripcion: Lee el encabezado de algun proyecto (.int, .float o .store) y de acuerdo al nombre.\n'\
    'Entrega la columna en la cual se localiza eso buscado.\n'\
    'Parametros Obligatorios:.\n'\
    '	-ruta: Ruta del proyecto .int o .int_hdr a leer.\n'\
    'Parametros Opcionales:.\n'\
    '	-name2find: Nombre de la variable a buscar.\n'\
    'Retorno:.\n'\
    '	Num_col: Numero de la columna donde se encuentra lo buscado, si Num_col=-9999 no se ha encontrado.\n'\
    #Carga el header de la cuenca e identifica en que columna esta lo buscado
    if ruta.endswith('.int') or ruta.endswith('.float') or ruta.endswith('.store'):
	ruta=ruta+'_hdr'
	f=open(ruta,'r')
    elif ruta.endswith('.int_hdr') or ruta.endswith('.float_hdr') or ruta.endswith('.store_hdr'):
	f=open(ruta,'r')
    else:		
	sys.exit("Error: No se encuentra el header del archivo base de la cuenca")
    Lineas=f.readlines()
    #Obtiene la primera entrada de cada linea
    Num_col=-9999
    for i in Lineas:
	if i.find(name2find)>0:
	    try:
		Num_col=int(i.split(':')[0])
	    except:
		pass
    f.close()
    return Num_col
def read_for_model(ruta_int,ruta_float,proyecto,ruta_store=None,model='SHIA'):
    #Funcion para identificar en donde esta algo dentro de un proyecto
    'Funcion: read_for_model\n'\
    'Descripcion: Lee los datos necesarios para la ejecucion del modelo a partir de la informacion\n'\
    'depositada en los binarios .int, .float y .store del proyecto.\n'\
    'Parametros Obligatorios:.\n'\
    '	-ruta_int: Ruta del proyecto .int.\n'\
    '	-ruta_float: Ruta del proyecto .float.\n'\
    '	-proyecto: Proyecto al cual se le asignan las propiedades, importado como:\n'\
    '	 from fwm.modelacion import modelos as mod, en el caso donde proyecto=mod.\n'\
    'Parametros Opcionales:.\n'\
    '	-ruta_store: Ruta de las condiciones iniciales de almacenamiento.\n'\
    'Retorno:.\n'\
    '	proyecto: El proyecto de modelacion con la informacion cargada.\n'\
    #Revisa cada una de las rutas
    if ruta_int.endswith('int_hdr') or ruta_int.endswith('int'):
	if ruta_int.endswith('int_hdr'): ruta_int=ruta_int[:-4]
    else:
	sys.exit('Error: La ruta del proyecto de enteros debe terminar en .int o .int_hdr')
    if ruta_float.endswith('float_hdr') or ruta_float.endswith('float'):
	if ruta_int.endswith('float_hdr'): ruta_float=ruta_float[:-4]
    else:
	sys.exit('Error: La ruta del proyecto de enteros debe terminar en .float o .float_hdr')
    if ruta_store<>None:
	if ruta_store.endswith('store_hdr') or ruta_store.endswith('store'):
	    if ruta_store.endswith('store_hdr'): ruta_store=ruta_store[:-4]
	else:
	    sys.exit('Error: La ruta del proyecto de enteros debe terminar en .store o .store_hdr')
    #Determina el numero de celdas yel nodata para el modelo
    Orden=read_proyect_hdr(ruta_int)[1]
    proyecto.nceldas=int(Orden['Numero de filas'])
    proyecto.nodata=int(Orden['noData_Value'])
    proyecto.dxp=float(Orden['cellsize_plana'])    
    nceldas=proyecto.nceldas
    proyecto.drena=cu.read_int_basin(ruta_int,1,nceldas,1)
    #De acuerdo al tipo de modelo determina que se va a buscar 
    if model=='SHIA':
	Cols=[read_proyect_identify(ruta_int,i) for i in ['Tipo Cauce','Pts Control']]
	Cols.extend([read_proyect_identify(ruta_float,i) for i in ['Longitud','Pendiente Laderas','Pendiente Cauce','hu','ks','kp','N','evp','Coord X','Coord Y']])
	flag=True
	for i in Cols:
	    if i<0: flag=False
	if flag:    	    
	    proyecto.tipo_celda=cu.read_int_basin(ruta_int,Cols[0],nceldas,1)
	    proyecto.control=cu.read_int_basin(ruta_int,Cols[1],nceldas,1)
	    proyecto.l_celda=cu.read_float_basin(ruta_float,Cols[2],nceldas,1)
	    proyecto.pend_celda=cu.read_float_basin(ruta_float,Cols[3],nceldas,1)
	    proyecto.pend_cauce=cu.read_float_basin(ruta_float,Cols[4],nceldas,1)
	    proyecto.hu=cu.read_float_basin(ruta_float,Cols[5],nceldas,1)
	    proyecto.ks=cu.read_float_basin(ruta_float,Cols[6],nceldas,1)
	    proyecto.kp=cu.read_float_basin(ruta_float,Cols[7],nceldas,1)
	    proyecto.man=cu.read_float_basin(ruta_float,Cols[8],nceldas,1)
	    proyecto.evp=cu.read_float_basin(ruta_float,Cols[9],nceldas,1)
	    proyecto.xy_celdas=cu.read_float_basin(ruta_float,[Cols[10],Cols[11]],nceldas,2)
	#Si se incluyo ruta de los almacenamientos genera los almacenamientos del modelo
	if ruta_store<>None:
	    Cols=[read_proyect_identify(ruta_store,i) for i in ['S1','S2','S3','S4','S5']]	    
	    proyecto.s=cu.read_float_basin(ruta_store,Cols,nceldas,5)
	else:
	    proyecto.s=np.zeros((5,nceldas))
	#Vuelve el proyecto con los datos actualizados
	return proyecto
def read_map_raster(ruta_map,proyecto=None):
    'Funcion: read_map\n'\
    'Descripcion: Lee un mapa raster soportado por GDAL.\n'\
    'Parametros Obligatorios:.\n'\
    '	-ruta_map: Ruta donde se encuentra el mpaa.\n'\
    'Parametros Opcionales:.\n'\
    '	-proyecto: Si se especifica un proyecto importado como: from fwm.modealcion import cuencas as cu.\n'\
    '	como proyecto=cu, entonces cu copiara todas las propieades del header del proyecto, exepcto dxp.\n'\
    'Retorno:.\n'\
    '	Si se especifica el proyecto: actualiza el proyecto y devuelve la matriz del mapa.\n'\
    '	Si no se especifica: devuelve la matriz del mapa y una lista con las caracteristicas:.\n'\
    '		En el siguiente orden: ncols,nrows,xll,yll,dx,nodata.\n'\
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
    if proyecto<>None:
	proyecto.ncols=ncols
	proyecto.nrows=nrows
	proyecto.nodata=noData
	proyecto.dx=dx
	proyecto.xll=xll
	proyecto.yll=yll
	return Mapa.T
    else:
	return Mapa.T,[ncols,nrows,xll,yll,dx,noData]
def read_coord_from_vector(ruta_map,col_id=None):
    'Funcion: read_coord_from_vector\n'\
    'Descripcion: Lee un mapa de puntos soportado por osgeo para obtener a partir de este las coordenadas.\n'\
    'Parametros Obligatorios:.\n'\
    '	-ruta_map: Ruta del archivo vectorial GIS donde se encuentran las estaciones.\n'\
    'Parametros Opcionales:.\n'\
    '	-col_id: Valor entero indicando la columna en la cual esta el id de cada estacion.\n'\
    '		comienza en 0.\n'\
    'Retorno:.\n'\
    '	-coordenadas x y y de cada punto.\n'\
    '	-ids: ID de cada una de las estaciones.\n'\
    #Abre el archivo
    d=ogr.Open(ruta_map)
    lyr=d.GetLayer()    
    #Evalua si va a sacar los ids
    f=lyr.GetFeature(0); flag=0
    if col_id<>None:
	if type(col_id)==int:
	    if col_id<f.GetFieldCount()-1:
		flag=1
	elif type(col_id)==str:
	    col_id=f.GetFieldIndex(col_id)
	    flag=1
    #Obtiene las coordenadas
    coord=[]; ids=[]
    for i in range(lyr.GetFeatureCount()):
	f=lyr.GetFeature(i)
	geo=f.GetGeometryRef()
	coord.append([geo.GetX(),geo.GetY()])
	if flag==1:
	    ids.append(f.GetField(col_id))
    if flag==1:
	return np.array(coord).T,ids
    else:
	return np.array(coord).T
	
#---------------------------------------------------------------------------------------------------
#Guardado y lectura de lluvia
#---------------------------------------------------------------------------------------------------	
#funciones de guardado
def write_rain_rain(ruta,rain,dt,ddmmaa,hhmmss):
    #Funcion desarrollada para guardar informacion de lluvia en formato del modelo
    f=open(ruta,'w')
    f.write("Cantidad Registros:   %d \n" % np.size(rain,axis=1))
    f.write("Cantidad estaciones:    %d \n" % np.size(rain,axis=0))
    f.write("dt:   %1.2f \n" % dt)
    f.write("No Data:   %d \n" % rain[rain<0][0])
    f.write("Fecha Inicio: %d:%d:%d \n" % (ddmmaa[0],ddmmaa[1],ddmmaa[2]))
    f.write("Hora Inicio: %d:%d:%d \n" % (hhmmss[0],hhmmss[1],hhmmss[2]))
    for i in range(np.size(rain,axis=0)):
	f.write("      %3d" % (i+1))
    f.write("\n")
    for i in rain.T:
	for j in i:
	    if j<0:
		f.write("     %3d" % j)
	    else:
		f.write("     %1.2f" % j)
	f.write("\n")
    f.close()
def write_rain_coord(ruta,ids,coord):
    #Escribe un archivo con estaciones
    f=open(ruta,'w')
    f.write('Numero de estaciones \n')
    f.write('%i \n' % np.size(coord,axis=1))
    f.write('Coordenadas \n')
    f.write('ID  Lat     Long \n')
    cont=0
    for i in coord.T:
	f.write('%i  %3.4f  %3.4f \n' % (ids[cont],float(i[0]),float(i[1])))
	cont+=1
    f.close()
#funciones de lectura
def read_rain_coord(ruta):
    #lee las coordenadas de lass estaciones de lluvia
    ids,y,x=np.loadtxt(ruta,skiprows=4,unpack=True,usecols=(0,1,2))
    N_est=np.size(ids)
    ids=ids.astype(int)
    est_coord=np.array([x,y])
    return ids,est_coord,N_est
def read_rain_rain(ruta):
    #lee los registros de lluvia
    f=open(ruta,'r')
    f.readline()
    f.readline()
    #Lee delta de tiempo y no data
    dt=int(f.readline().split()[-1])
    noData=int(f.readline().split()[-1])
    #Lee la fecha de inicio de la simulacion
    ddmmaa=f.readline().split()[-1]
    ddmmaa=ddmmaa.split(':')
    ddmmaa=np.array([int(i) for i in ddmmaa])
    hhmmss=f.readline().split()[-1]
    hhmmss=hhmmss.split(':')
    hhmmss=np.array([int(i) for i in hhmmss])
    f.close()
    #lee los registros
    lluvia=np.loadtxt(ruta,skiprows=9,unpack=True)
    #devuelve los datos leidos
    return dt,noData,ddmmaa,hhmmss,lluvia

    
#funciones propias de tin
def find_tin(coordXY,proy_rain=None):
    #Guarda las condiciones iniciales de un modelo
    'Funcion: find_tin\n'\
    'Descripcion: A partir de pares coordenados encuentra los triangulos conformados.\n'\
    'Parametros Obligatorios:.\n'\
    '	-coordXY: Vector con las coordenadas X y Y de las estaciones de precipitacion.\n'\
    'Parametros Opcionales:.\n'\
    '	-proy_rain: Modulo de lluvia (from fwm.modelacion import lluvia), para asignarle \n'\
    '		a este el valor de lluvia.tin.'
    'Retorno:.\n'\
    '	Si no se asigna proy_rain, retorna la matriz conformate de los triangulos.\n'\
    #calcula la malla de triangulos
    if coordXY.shape[1]<>2: 
	coordXY=coordXY.T
	if coordXY.shape[1]<>2:
	    sys.exit('Error: La matriz debe contener dos columnas.')
    T=Delaunay(coordXY)
    T=T.vertices.T
    #si se asigno un proyecto de lluvia le asigna el valor de lluvia.tin
    if proy_rain<>None:
	proy_rain.tin=T+1
    else:
	return T

#---------------------------------------------------------------------------------------------------
#Guardado y lectura de Calibracion
#---------------------------------------------------------------------------------------------------
#funciones de lectura
def read_shia_calib(ruta):
    #Lee los parametros de calibracion del modelo hidrologico
    #Abre el archivo y lo lee
    f=open(ruta,'r')
    Lista=f.readlines()
    f.close()
    #toma los parametros de claibracion
    Calib=[]
    for i in Lista[:22]:
	try:
	    Calib.append(float(i))
	except:
	    pass
    #toma los paramtros de la OCG
    OCG=[]
    for i in Lista[23:]:
	try:
	    OCG.append(float(i))
	except:
	    pass
    return np.array(Calib[:-1]),np.array(OCG),Calib[-1]
	
#---------------------------------------------------------------------------------------------------
#Guardado y grabado shapes
#---------------------------------------------------------------------------------------------------
def save_shp_cuenca(basin,nceldas,ruta,dx=30.0):
	#Obtiene el perimetro de la cuenca 
	nperim = cu.basin_perim_find(basin,nceldas)
	basinPerim=cu.basin_perim_cut(nperim)
	#Construye el shp 
	if ruta.endswith('.shp')==False:
		ruta=ruta+'.shp'
	#Genera el shapefile
	spatialReference = osgeo.osr.SpatialReference()
	spatialReference.ImportFromProj4('+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs')
	driver = osgeo.ogr.GetDriverByName('ESRI Shapefile')
	if os.path.exists(ruta):
	     driver.DeleteDataSource(ruta)
	shapeData = driver.CreateDataSource(ruta)
	layer = shapeData.CreateLayer('layer1', spatialReference, osgeo.ogr.wkbPolygon)
	layerDefinition = layer.GetLayerDefn()
	new_field=osgeo.ogr.FieldDefn('Area[km2]',osgeo.ogr.OFTReal)
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
	feature.SetField('Area[km2]',nceldas*dx**2/1e6)
	layer.CreateFeature(feature)
	poly.Destroy()
	ring.Destroy()
	feature.Destroy()
	shapeData.Destroy()
def save_shp_network(basin,nceldas,ruta,dx=30.0,umbral=1000):
	#division de la cuenca 
	acum=cu.basin_acum(basin,nceldas)
	cauce,nod_f,n_nodos=cu.basin_subbasin_nod(basin,acum,umbral,nceldas)
	sub_pert,sub_basin=cu.basin_subbasin_find(basin,nod_f,n_nodos,nceldas)
	sub_basins=cu.basin_subbasin_cut(n_nodos)
	sub_horton,nod_hort=cu.basin_subbasin_horton(sub_basins,nceldas,n_nodos)
	sub_hort=cu.basin_subbasin_find(basin,nod_hort,n_nodos,nceldas)[0]
	cauceHorton=sub_hort*cauce
	#Obtiene la red en manera vectorial 
	nodos = cu.basin_stream_nod(basin,acum,umbral,nceldas)[1]
	netsize = cu.basin_netxy_find(basin,nodos,cauceHorton,nceldas)
	net=cu.basin_netxy_cut(netsize,nceldas)
	cortes=np.where(net[0,:]==-999)
	cortes=cortes[0].tolist()
	cortes.insert(0,0)
	#Escribe el shp de la red hidrica
	if ruta.endswith('.shp')==False:
		ruta=ruta+'.shp'
	spatialReference = osgeo.osr.SpatialReference()
	spatialReference.ImportFromProj4('+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs')
	driver = osgeo.ogr.GetDriverByName('ESRI Shapefile')
	if os.path.exists(ruta):
	     driver.DeleteDataSource(ruta)
	shapeData = driver.CreateDataSource(ruta)
	layer = shapeData.CreateLayer('layer1', spatialReference, osgeo.ogr.wkbLineString)
	layerDefinition = layer.GetLayerDefn()
	new_field=osgeo.ogr.FieldDefn('Long[km]',osgeo.ogr.OFTReal)
	layer.CreateField(new_field)
	new_field=osgeo.ogr.FieldDefn('Horton',osgeo.ogr.OFTInteger)
	layer.CreateField(new_field)
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
		#featureFID+=1
		layer.CreateFeature(feature)
		line.Destroy()
		feature.Destroy()
	shapeData.Destroy()
def save_shp_estaciones(ids,est_coord,nombre):
    #guarda un shapefile de las estaciones de trabajo, de momento solo usa wgs84
    #Se fija si el nombre termina con .shp
    if nombre.endswith('.shp')==False:
	nombre=nombre+'.shp'
    #Genera el shapefile
    spatialReference = osgeo.osr.SpatialReference()
    spatialReference.ImportFromProj4('+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs')
    driver = osgeo.ogr.GetDriverByName('ESRI Shapefile')
    shapeData = driver.CreateDataSource(nombre)
    layer = shapeData.CreateLayer('layer1', spatialReference, osgeo.ogr.wkbPoint)
    layerDefinition = layer.GetLayerDefn()
    new_field=osgeo.ogr.FieldDefn('Estacion',osgeo.ogr.OFTInteger)
    layer.CreateField(new_field)
    #Mete todos los puntos
    featureIndex=0
    contador=0
    #Calcula el tamano de la muestra
    N=np.size(est_coord,axis=1)
    if N>1:
	for i in est_coord.T:
	    point = osgeo.ogr.Geometry(osgeo.ogr.wkbPoint)
	    point.SetPoint(0, float(i[0]), float(i[1]))
	    feature = osgeo.ogr.Feature(layerDefinition)
	    feature.SetGeometry(point)
	    feature.SetFID(featureIndex)
	    feature.SetField('Estacion',ids[contador])
	    contador+=1
	    layer.CreateFeature(feature)
	    point.Destroy()
	    feature.Destroy()
    else:
	point = osgeo.ogr.Geometry(osgeo.ogr.wkbPoint)
	point.SetPoint(0, float(est_coord[0,0]), float(est_coord[1,0]))
	feature = osgeo.ogr.Feature(layerDefinition)
	feature.SetGeometry(point)
	feature.SetFID(featureIndex)
	feature.SetField('Estacion',ids)
	contador+=1
	layer.CreateFeature(feature)
	point.Destroy()
	feature.Destroy()
    shapeData.Destroy()	
def save_shp_triangulos(TIN,est_coord,nombre):
    #guarda un shapefile de las estaciones de trabajo, de momento solo usa wgs84
    #Se fija si el nombre termina con .shp
    if nombre.endswith('.shp')==False:
	    nombre=nombre+'.shp'
    #Genera el shapefile
    spatialReference = osgeo.osr.SpatialReference()
    spatialReference.ImportFromProj4('+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs')
    driver = osgeo.ogr.GetDriverByName('ESRI Shapefile')
    shapeData = driver.CreateDataSource(nombre)
    layer = shapeData.CreateLayer('layer1', spatialReference, osgeo.ogr.wkbPolygon)
    layerDefinition = layer.GetLayerDefn()
    new_field=osgeo.ogr.FieldDefn('TIN',osgeo.ogr.OFTInteger)
    layer.CreateField(new_field)
    #Mete todos los puntos
    featureIndex=0
    contador=0
    TIN=TIN.T-1
    for i in TIN:    
	ring = osgeo.ogr.Geometry(osgeo.ogr.wkbLinearRing)
	for j in i:
	    ring.AddPoint(est_coord[0,j], est_coord[1,j])
	poly=osgeo.ogr.Geometry(osgeo.ogr.wkbPolygon)
	poly.AddGeometry(ring)
	feature = osgeo.ogr.Feature(layerDefinition)
	feature.SetGeometry(poly)
	feature.SetFID(featureIndex)
	feature.SetField('TIN',contador+1)
	contador+=1
	layer.CreateFeature(feature)
	poly.Destroy()
	ring.Destroy()
	feature.Destroy()
    shapeData.Destroy()

#---------------------------------------------------------------------------------------------------
#Graficas
#---------------------------------------------------------------------------------------------------
def plot_basin_map(basin,vec,nceldas,Min=None,Max=None,ruta=None,mostrar='si',barra='si'):
    #Plotea en la terminal como mapa un vector de la cuenca
    'Funcion: write_proyect_int_ext\n'\
    'Descripcion: Genera un plot del mapa entrgeado.\n'\
    'del mismo en forma de mapa \n'\
    'Parametros Obligatorios:.\n'\
    '	-basin: Vector con la forma de la cuenca.\n'\
    '	-vec: Vector con los valores a plotear.\n'\
    '	-nceldas: Cantidad de filas que componen el vector basin.\n'\
    'Parametros Opcionales:.\n'\
    '	-Min: Valor minimo del plot, determina el rango de colores.\n'\
    '	-Max: Valor maximo del plot, determina el rango de colores.\n'\
    '	-ruta: Ruta en la cual se guarda la grafica.\n'\
    '	-mostrar: Muestra o no la grafica, defecto: si.\n'\
    '	-barra: Muestra o no la barra de colores, defecto: si.\n'\
    'Retorno:.\n'\
    '	Actualizacion del binario .int\n'\
    #Delimita el mapa
    Mcols,Mrows=cu.basin_2map_find(basin,nceldas)
    m,mxll,myll=cu.basin_2map(basin,vec,Mcols,Mrows,nceldas)
    m[m==cu.nodata]=np.nan
    if Min<>None and Max<>None:
	pl.imshow(m.T,norm=LogNorm(vmin=Min,vmax=Max))
    else:
	pl.imshow(m.T)
    if barra=='si':
	pl.colorbar()
    if ruta<>None and type(ruta)==str:
	pl.savefig(ruta,bbox_inches='tight')
    if mostrar=='si':
	pl.show()
def plot_basin_elegant(basin,nceldas,dx,dem_rute=None,vec=None,ruta=None,lines_spaces=0.02,
	extra_lat=0.0,extra_long=0.0,alpha=0.5,umbral_red=1000,horton='si',red='si',
	figsize=None,ppal_stream='no',longCeldas=None,elevation=None,xy=None,
	xycolor='red',colorTable=None):
	#Plotea en la terminal como mapa un vector de la cuenca, en este caso lo hace bonito, recibe
    #varios mapas
	'Funcion: plot_basin_elegant\n'\
	'Descripcion: Genera un plot del mapa elegido, de una ves lo hace elegante.\n'\
	'Tiene que estar en wgs84 para que funcione de manera correcta.\n'\
	'Parametros Obligatorios:.\n'\
	'	-basin: Vector con la forma de la cuenca.\n'\
	'	-nceldas: Cantidad de filas que componen el vector basin.\n'\
	'	-dx: Distancia entre celdas.\n'\
	'Parametros Opcionales:.\n'\
	'	-vec: Vector con la variable que se quiere plotear.\n'\
	'	-colorTable: Paleta de colores para plotear, defecto=None.\n'\
	'	-dem_rute: direccion al modelo de elevacion de la zona, si lo colocan lo pone de fondo.\n'\
	'	-ruta: Ruta en la cual se guarda la grafica.\n'\
	'	-lines_spaces: Espacios entre las lineas de la cuadricula.\n'\
	'	-extra_lat: Valor fuera de los margenes de la latitud para plotear.\n'\
	'	-extra_long: Valor fuera d elos margenes de la longitud para plotear.\n'\
	'	-alpha: Nivel de transparecnia de la variable a plotear (Defecto: 0.5).\n'\
	'	-umbral_red: Umbral para plotar la red (cantidad de celdas), Defecto=1000.\n'\
	'	-horton: Pintar red por orden de horton, Defecto: si.\n'\
	'	-figsize: Determina el tamano del mapa a plotear, Defecto: None.\n'\
	'	-ppal_stream: Pinta encima el cauce principal, Defecto: no.\n'\
	'	-longCeldas: Distancia plana entre celdas, solo sera util si ppal_stream=si.\n'\
	'	-elevation: Elevacion de cada celda, solo sera util si ppal_stream=si.\n'\
	'	-xy: Cordenadas de puntos que se graficaran sobre el mapa, Defecto: None.\n'\
	'		-xycolor: Color de los puntos, Defecto: Rojo.\n'\
    #Obtiene el borde de la cuenca
	nperim = cu.basin_perim_find(basin,nceldas)
	perim=cu.basin_perim_cut(nperim)	
	#Obtiene la red hidrica con orden de horton
	acum=cu.basin_acum(basin,nceldas)	
	cauce,nod_f,n_nodos=cu.basin_subbasin_nod(basin,acum,umbral_red,nceldas)
	if horton=='si':
		sub_pert,sub_basin=cu.basin_subbasin_find(basin,nod_f,n_nodos,nceldas)
		sub_basins=cu.basin_subbasin_cut(n_nodos)
		sub_horton,nod_hort=cu.basin_subbasin_horton(sub_basins,nceldas,n_nodos)
		sub_hort=cu.basin_subbasin_find(basin,nod_hort,n_nodos,nceldas)[0]
	else:
		sub_hort=1
	cauce,nodos,trazado,n_nodos,n_cauce = cu.basin_stream_nod(basin,acum,umbral_red,nceldas)
	cauce=sub_hort*cauce
	netsize = cu.basin_netxy_find(basin,nodos,cauce,nceldas)	
	net=cu.basin_netxy_cut(netsize,nceldas)
	cortes=np.where(net[0,:]==-999)
	cortes=cortes[0].tolist()
	cortes.insert(0,0)
	#Cauce principal si si 
	if ppal_stream=='si':
		ppal_nceldas,punto = cu.basin_ppalstream_find(basin,nodos,longCeldas,elevation,nceldas)
		ppal = cu.basin_ppalstream_cut(ppal_nceldas,nceldas)
	#Obtiene los parametros del mapa a partir de la cuenca
	Mcols,Mrows=cu.basin_2map_find(basin,nceldas)
	Map,mxll,myll=cu.basin_2map(basin,basin[0],Mcols,Mrows,nceldas)
	#Monta la base del mapa a partir de las prop de la cuenca
	longs=np.array([mxll+0.5*dx+i*cu.dx for i in range(Mcols)])
	lats=np.array([myll+0.5*dx+i*cu.dx for i in range(Mrows)])
	X,Y=np.meshgrid(longs,lats)
	Y=Y[::-1]
	if figsize<>None:
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
	#Si hay mapa de elevaciones lo plotea
	if dem_rute<>None:
		DataSet=gdal.Open(dem_rute)
		#Cantidad de columnas y de filas
		Xsize=DataSet.RasterXSize
		Ysize=DataSet.RasterYSize
		#Ubicacion del raster
		GeoTransform=DataSet.GetGeoTransform()
		Xcorner=GeoTransform[0]
		Ycorner=GeoTransform[3]
		dem_dx=GeoTransform[1]
		dem_dy=GeoTransform[-1]
		#Lee la informacion y enmascara los datos nulos
		Data=DataSet.ReadAsArray()
		Data=np.ma.array(Data,mask=Data==-9999.0)
		Data=np.ma.array(Data,mask=Data>4000.0)
		#Generacion de una matriz de coordenadas
		longsD=np.array([Xcorner+0.5*dem_dx+i*dem_dx for i in range(Xsize)])
		latsD=np.array([Ycorner+0.5*dem_dy+i*dem_dy for i in range(Ysize)])
		dem_X,dem_Y=np.meshgrid(longsD,latsD)
		demXm,demYm=m(dem_X,dem_Y)
		#Imprime el dem de fondo
		m.contourf(demXm, demYm, Data, 100, cmap=pl.cm.terrain,alpha=0.95)
	#Si hay una variable la monta
	if vec<>None:
		MapVec,mxll,myll=cu.basin_2map(basin,vec,Mcols,Mrows,nceldas)
		MapVec[MapVec==cu.nodata]=np.nan
		if colorTable<>None:
			cs = m.contourf(Xm, Ym, MapVec.T, 25, alpha=alpha,cmap=colorTable)
		else:
			cs = m.contourf(Xm, Ym, MapVec.T, 25, alpha=alpha)
		cbar = m.colorbar(cs,location='bottom',pad="5%")	
	#Plotea el contorno de la cuenca y la red 
	xp,yp=m(perim[0],perim[1])  
	m.plot(xp, yp, color='r',lw=2)
	#red
	if red=='si':
		if horton=='si':
			for i,j in zip(cortes[:-1],cortes[1:]):
			    if net[0,i+1]==1: co=(0.4,0.4,1.0); w=1
			    if net[0,i+1]==2: co=(0.2,0.2,1.0); w=1.5
			    if net[0,i+1]==3: co=(0,0,1.0); w=2
			    if net[0,i+1]==4: co=(0,0,0.8); w=2.5
			    if net[0,i+1]==5: co=(0,0,0.6); w=3.0
			    if net[0,i+1]==6: co=(0,0,0.4); w=3.5
			    xn,yn=m(net[1,i+1:j],net[2,i+1:j])
			    m.plot(xn,yn,c=co,lw=w)
		else:
			for i,j in zip(cortes[:-1],cortes[1:]):
				xn,yn=m(net[1,i+1:j],net[2,i+1:j])
				m.plot(xn,yn,c='k',lw=2)
	#Cauce principal
	if ppal_stream=='si':
		xp,yp=m(ppal[2],ppal[3])
		m.plot(xp,yp,c=(0,0,0.3),lw=3.5)
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

def plot_basin_map_noBorder(basin,vec,nceldas,ruta,Min=None,Max=None):
    #Plotea en la terminal como mapa un vector de la cuenca
	'Funcion: write_proyect_int_ext\n'\
	'Descripcion: Genera un plot del mapa entrgeado.\n'\
	'del mismo en forma de mapa \n'\
	'Parametros Obligatorios:.\n'\
	'	-basin: Vector con la forma de la cuenca.\n'\
	'	-vec: Vector con los valores a plotear.\n'\
	'	-nceldas: Cantidad de filas que componen el vector basin.\n'\
	'	-ruta: Ruta en la cual se guarda la grafica.\n'\
	'Parametros Opcionales:.\n'\
	'	-Min: Valor minimo del plot, determina el rango de colores.\n'\
	'	-Max: Valor maximo del plot, determina el rango de colores.\n'\
	'Retorno:.\n'\
	'	Actualizacion del binario .int\n'\
    #Delimita el mapa y guarda el mapa
	Mcols,Mrows=cu.basin_2map_find(basin,nceldas)
	m,mxll,myll=cu.basin_2map(basin,vec,Mcols,Mrows,nceldas)
	m[m==cu.nodata]=np.nan
    #Plot de la barra
	fig=pl.figure(frameon=False)
	ax=fig.add_subplot(111)
	ax.axis('off')    
	try:
		pl.imshow(m.T,norm=LogNorm(vmin=Min,vmax=Max))
		pl.colorbar(orientation="horizontal")
	except:
		pl.imshow(m.T)   
	pl.colorbar(orientation="horizontal")
	pl.gca().set_visible(False)    	
	pl.savefig(ruta+"_bar.png",bbox_inches='tight',transparent=True)	
    #Plot del mapa
	fig=pl.figure(frameon=False)
	ax=fig.add_subplot(111)
	ax.axis('off')
	ax.get_xaxis().set_visible(False)
	ax.get_yaxis().set_visible(False)    	
	m[m==0]=np.nan
	try:
		pl.imshow(m.T,norm=LogNorm(vmin=Min,vmax=Max))
	except:
		pl.imshow(m.T)   
	pl.tight_layout()
	pl.savefig(ruta+"_map.png",bbox_inches='tight',transparent=True,pad_inches = 0)
	#Obtiene los cuatro bordes
	f=open(ruta+"_extent.txt","w")
	f.write("%.8f \n" % mxll)
	mxll=mxll+cu.dx*m.shape[0]
	f.write("%.8f \n" % mxll)
	f.write("%.8f \n" % myll)
	myll=myll+cu.dx*m.shape[1]
	f.write("%.8f \n" % myll)
	f.close()
def plot_sim_single(Qs,Qo=None,mrain=None,Dates=None,ruta=None,Mostrar='si',figTam=(12,10)):
    fig=pl.figure(facecolor='w',edgecolor='w',figsize=figTam)
    ax1=fig.add_subplot(111)
    #Fechas
    if Dates==None:
	ejeX=range(len(Qs[0]))
    else:
	ejeX=Dates
    #Grafica la lluvia
    if mrain<>None:
	ax1AX=pl.gca()
	ax2=ax1.twinx()
	ax2AX=pl.gca()
	ax2.fill_between(ejeX,0,mrain,alpha=0.4,color='blue',lw=0)
	ax2.set_ylabel('Precipitacion [mm]',size=14)
	ax2AX.set_ylim(ax2AX.get_ylim() [::-1])    
    #grafica las hidrografas
    ColorSim=['r','g','k','c','y']
    #for pos,i in enumerate(Qs):
    ax1.plot(ejeX,Qs,'r',lw=1.5,label='Simulado ')    
    if Qo<>None: ax1.plot(ejeX,Qo,'b',lw=2,label='Obsevado')
    #Pone elementos en la figura
    ax1.set_xlabel('Tiempo [$min$]',size=14)
    ax1.set_ylabel('Caudal [$m^3/seg$]',size=14)
    lgn1=ax1.legend(bbox_to_anchor=(0.01, 0.99), loc=2, borderaxespad=0., prop={'size':12})
    if Mostrar=='si':
	pl.show()
    if ruta<>None:
	pl.savefig(ruta, bbox_inches='tight')
#---------------------------------------------------------------------------------------------------
#Funciones de plot de geomorfologia
#---------------------------------------------------------------------------------------------------
def plot_ppal_stream(dist,elev,ruta=None,ventana=10,normed=False):
	if normed==True:
		elev=elev-elev.min(); elev=(elev/elev.max())*100.0
		dist=dist-dist.min(); dist=(dist/dist.max())*100.0
	#Calcula cositas antes
	dist=dist[::-1]
	elevSuave=pandas.Series(elev)
	elevSuave=pandas.rolling_mean(elevSuave,ventana)
	r=np.polyfit(dist,elev,1)
	#Figura
	fig=pl.figure(edgecolor='w',facecolor='w')
	ax=fig.add_subplot(111)
	box = ax.get_position()
	ax.set_position([box.x0, box.y0 + box.height * 0.1,
		box.width, box.height * 0.9])	
	ax.plot(dist,elev,c='b',lw=2,label='Perfil')
	ax.plot(dist,elevSuave,c='g',lw=1.5,label='Perfil Suavizado')
	ax.plot(dist,dist*r[0]+r[1],c='r',lw=1.5,label='Pendiente')
	ax.grid()
	if normed==False:
		ax.set_xlabel('Distancia a la salida $[m]$',size=14)
		ax.set_ylabel('Elevacion $[m.s.n.m]$',size=14)
	if normed==True:
		ax.set_xlabel('Distancia a la salida $[\%]$',size=14)
		ax.set_ylabel('Elevacion $[\%]$',size=14)
	r=r[0]*100; s="%4.2f" % r
	l=dist[0]/1000.0; l="%4.2f" % l  
	ax.text(0.7,0.2,'$S_{0}='+s+'[\%]$ \n $L='+l+'[Km]$',size=16,transform=ax.transAxes)
	lgn1=ax.legend(loc='upper center', bbox_to_anchor=(0.5, -0.12),
		fancybox=True, shadow=True, ncol=3)
	if ruta<>None:
		pl.savefig(ruta, bbox_inches='tight')
		pl.close('all')
	else:
		pl.show()
	return r/100.0,dist[0]/1000.0
def plot_hipsometric(ppal_hipso,basin_hipso,ruta=None,ventana=10,normed=False):
	#Suaviza la elevacion en el cuace ppal
	elevPpal=pandas.Series(ppal_hipso[1])
	elevBasin=basin_hipso[1]
	if normed==True:
		elevPpal=elevPpal-elevPpal.min()
		elevPpal=(elevPpal/elevPpal.max())*100.0
		elevBasin=elevBasin-elevBasin.min()
		elevBasin=(elevBasin/elevBasin.max())*100.0
	elevPpal=pandas.rolling_mean(elevPpal,ventana)
	ppal_acum=(ppal_hipso[0]/ppal_hipso[0,-1])*100
	basin_acum=(basin_hipso[0]/basin_hipso[0,0])*100
	#Genera el plot
	fig=pl.figure(edgecolor='w',facecolor='w')
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
		
#---------------------------------------------------------------------------------------------------
#Funciones varias
#---------------------------------------------------------------------------------------------------
def control_minimize(control):
    #Funcion para tomar el vector de los puntos de control y entregar un vector con la posicion y el id de los puntos de control
    #encuentra la posicion
    posicion=np.where(control>0)[1]
    ids=control[control>0]
    #arma con eso un vector
    control_min=np.zeros((2,len(posicion)),dtype=int)
    control_min[0]=posicion+1
    control_min[1]=ids
    return control_min,len(posicion)
def sort_rain(id_dst,id_src,rain_src,coords=None,pp=1,nodata=-999.0):
    #toma los ids del orden de las coordenadas y los ids en el roden del registro los reordena 
    #para entregar en el orden de las coordenadas y si hay estaciones ficticias las llena
    'Funcion: sort_rain\n'\
    'Descripcion: Reorganiza la lluvia leida de acuerdo a los ids de las coordenadas, si .\n'\
    'debe interpolar lluvia sobre estaciones ficticias lo hace.\n'\
    'Parametros Obligatorios:.\n'\
    '	-id_dst: Lista con los ids tomados de las coordenadas.\n'\
    '	-id_src: Lista con los ids tomados de los registros de lluvia (orden de la matriz de lluvia).\n'\
    '	-rain_src: Matriz con los registros de lluvia leidos.\n'\
    'Parametros Opcionales:.\n'\
    '	-coord: Vector con las coordenadas leidas, si se asigna es utilizado para interpolar.\n'\
    '	-pp: potencia utilizada para sesgar las distancias en la interpolacion idw, defecto pp=1.\n'\
    '	-nodata: Valor que representa los datos faltantes en las series de lluvia.\n'\
    'Retorno:.\n'\
    '	-rain_dst: Matriz de lluvia organizada de acuerdo a las coordenadas y con valores interpolados.\n'\
    #Genera una matriz de lluvia vacia
    rain_dst=np.zeros((len(id_dst),rain_src.shape[1]))
    encontrado=np.zeros(len(id_dst),dtype=int)
    #Itera sobre los ids de las coordenadas para buscarlos sobre los ids de la info de lluvia
    for pos,i_dts in enumerate(id_dst):
	#Determina si esa estacion esta o no esta en ambas listas
	try:
	    pos_id=id_src.index(i_dts)
	    encontrado[pos]=1	    
	    #Si esta en la lista la copia
	    rain_dst[pos]=rain_src[pos_id]
	except:
	    encontrado[pos]=0
    #Si se dieron las coordenadas interpola sobre las estaciones vacias
    m1=np.isnan(rain_dst)
    m2=rain_dst<0
    r_aux=np.ma.array(rain_dst,mask=m1+m2)
    print '----------------'
    if coords<>None:
	for pos,i in enumerate(id_dst):
	    if encontrado[pos]==0:
		#Si la posicion no se ha llenado interpola
		W=np.zeros(len(id_dst))
		#Encuentra los pesos
		for pos2,coord in enumerate(coords):		    
		    if pos2<>pos:
			#Solo calcula el peso para todos menos el mismo
			W[pos2]=1/(np.sqrt((coords[0][pos]-coord[0])**2+(coords[1][pos]-coord[1])**2)**pp)
		Wm=np.repeat(W,rain_dst.shape[1]).reshape(len(W),rain_dst.shape[1])
		r=r_aux*Wm
		rain_dst[pos]=r.sum(axis=0)/W.sum()    
    rain_dst[np.isnan(rain_dst)]=nodata
    return rain_dst

#---------------------------------------------------------------------------------------------------
#Modelos y geomorfologia agregados
#---------------------------------------------------------------------------------------------------	
def shia_lag_day(rain,Param,Watershed,Storage=None,Ext_param=[1.5,0.5]):
    #corre el modelo shia a escala diaria para unos parametros y una lluvia dada
    'Funcion: shia_lag_day\n'\
    'Descripcion: Ejecuta el modelo agregado shia a escala diaria.\n'\
    'Parametros Obligatorios:.\n'\
    '	-rain: Vector numpy con la lluvia en los diferentes intervalos de tiempo.\n'\
    '	-Param: Lista con los parametros de calibracion del modelo [Hu, Ks, Kp, Kpp, Ts, Tss, Tb].\n'\
    '	-Watershed: Lista con los parametros de la cuenca [Area, Elevacion].\n'\
    'Parametros Opcionales:.\n'\
    '	-Storage: Vector numpy con el estado de almacenamiento inicial de los 4 tanques [S1,S2,S3,S4].\n'\
    '	 si no se indica, comienza en cero.\n'\
    '	-Ext_param: Lista con parametros extras [Exp Inf=1.5, Exp EVP=0.5].\n'\
    'Retorno:.\n'\
    '	-Caudal: Caudal simulado a escala diaria.\n'\
    #Realiza verificaciones
    if len(Param)<7 or len(Param)>7: sys.exit('Error: Se deben indicar los 7 parametros de calibracion')
    if len(Ext_param)<>2: sys.exit('Error: Se deben indicar los 2 parametros extra de calibracion')
    if Storage==None:
	Storage=np.zeros(4)
    #calcula la evp
    Evp=4.658*np.exp(-0.0002*Watershed[1])
    #Itera para cada intervalo de tiempo
    Q=np.zeros(rain.shape[0])
    Param=np.array(Param,dtype=float)
    for pos,p in enumerate(rain):
	if np.isnan(p): p=0.0
	#Tanque 1: Almacenamiento capilar
	Entra_1=np.min([p*(1-(Storage[0]/Param[0])**Ext_param[0]),Param[0]-Storage[0]])
	Storage[0]+=Entra_1
	Salida_1=np.min([Evp*(Storage[0]/Param[0])**Ext_param[1],Storage[0]])
	Storage[0]-=Salida_1
	#tanque 2: Almacenamiento superficial
	Entra_2=np.max([p-Entra_1-Param[1],0])
	Storage[1]+=Entra_2
	Salida_2=(1/Param[4])*Storage[1]
	Storage[1]-=Salida_2
	#tanque 3: Almacenamiento sub-superficial
	Entra_3=np.max([p-Entra_1-Entra_2-Param[2],0])
	Storage[2]+=Entra_3
	Salida_3=(1/Param[5])*Storage[2]
	Storage[2]-=Salida_3
	#tanque 4: Almacenamiento base
	Entra_4=np.max([p-Entra_1-Entra_2-Entra_3-Param[3],0])
	Storage[3]+=Entra_4
	Salida_4=(1/Param[6])*Storage[3]
	Storage[3]-=Salida_4
	#Calcula el caudal del intervalo de tiempo
	Q[pos]=((Salida_2+Salida_3+Salida_4)*1000*Watershed[0])/86400
    return Q,Storage   
def basin_Tc(basin,DEM,DIR,dx,nceldas,ncols,nrows):
	#Calcula lo que se necesita para sacar los parametros
	acum,longCeld,S0,Elev=cu.basin_basics(basin,DEM,DIR,cu.ncols,cu.nrows,nceldas)
	slope=cu.basin_arc_slope(basin,DEM,nceldas,ncols,nrows)
	Lpma,puntto=cu.basin_findlong(basin,nceldas)
	cauce,nodos,trazado,n_nodos,n_cauce = cu.basin_stream_nod(basin,acum,1000,nceldas)
	ppal_nceldas,punto = cu.basin_ppalstream_find(basin,nodos,longCeld,Elev,nceldas)
	ppal = cu.basin_ppalstream_cut(ppal_nceldas,nceldas)
	nperim = cu.basin_perim_find(basin,nceldas)
	#Obtiene los parametros 
	Area=(nceldas*dx**2)/1e6
	Perim=nperim*dx/1000.0
	Lcau=ppal[1,-1]/1000.0
	Scau=np.polyfit(ppal[1,::-1],ppal[0],1)[0]*100
	Scue=slope.mean()*100
	Hmin=Elev[-1]; Hmax=Elev[puntto]; Hmean=Elev.mean()
	HCmax=Elev[punto]
	#Genera un diccionario con las propiedades de la cuenca 
	CuencaParam={'Area[km2]': Area,
		'Perimetro[km]':Perim,
		'Pend Cauce [%]':Scau,
		'Long Cau [km]': Lcau,
		'Pend Cuenca [%]': Scue,
		'Long Cuenca [km]': Lpma,
		'Hmax [m]':Hmax,
		'Hmin [m]':Hmin,
		'Hmean [m]':Hmean,
		'H Cauce Max [m]':HCmax}
	#Calcula los tiempos de concentracion
	Tiempos={}
	Tc=0.3*(Lcau/(Scue**0.25))**0.75
	Tiempos.update({'US Army': Tc})
	Tc=0.3*(Lcau/((Hmax-Hmin)/Lcau)**0.25)**0.75
	Tiempos.update({'Direccion Carreteras Espana': Tc})
	Tc=(0.02*Lpma**0.77)/((Scau)**0.385)
	Tiempos.update({'Kiprich': Tc})
	Tc=8.157*((Area**0.316)/(((Scau*100)**0.17)*Scue**0.565))
	Tiempos.update({'Campo y Munera': Tc})
	Tc=(4*np.sqrt(Area)+1.5*Lcau)/(0.8*np.sqrt(Hmean))
	Tiempos.update({'Giandotti':Tc})
	Tc=5*((Lcau*0.62)/np.sqrt(Scau*10))**0.5
	Tiempos.update({'John Stone': Tc})
	Tc=(Lcau/Area)*(np.sqrt(Area)/Scau)
	Tiempos.update({'Ventura': Tc})
	Tc=0.3*(Lcau/(((HCmax-Hmin)/Lcau)*100)**0.25)**0.75
	Tiempos.update({'Temez': Tc})
	return CuencaParam,Tiempos
		    
#---------------------------------------------------------------------------------------------------
#Evaluacion de eficiencia
#---------------------------------------------------------------------------------------------------
def eval_nash(s_o,s_s):
    med_s_o=np.nansum(s_o)/np.sum(np.isfinite(s_o))
    Dif_sim=(s_o-s_s[:len(s_o)])**2
    Dif_med=(s_o-med_s_o)**2
    E=1-(np.nansum(Dif_sim)/np.nansum(Dif_med))
    return E
def eval_rmse(So,Ss):
    Dif=(So-Ss)**2
    Dif=np.ma.array(Dif,mask=np.isnan(Dif))
    return np.sqrt(Dif.sum()/Dif.shape[0])
def eval_rmse_log(So,Ss):
    So=np.ma.array(So,mask=np.isnan(So))
    Ss=np.ma.array(Ss,mask=np.isnan(Ss))
    So=np.log(So); Ss=np.log(Ss)
    Dif=(So-Ss)**2
    Dif=np.ma.array(Dif,mask=np.isnan(Dif))
    return np.sqrt(Dif.sum()/Dif.shape[0])
def eval_t_pico(s_o,s_s,dt):    
    max_o=np.argmax(np.ma.array(s_o,mask=np.isnan(s_o)))
    max_s=np.argmax(np.ma.array(s_s,mask=np.isnan(s_s)))
    dif_tpico=(max_o-max_s)*dt
    return dif_tpico
def eval_q_pico(s_o,s_s):
    max_o=np.argmax(np.ma.array(s_o,mask=np.isnan(s_o)))
    max_s=np.argmax(np.ma.array(s_s,mask=np.isnan(s_s)))
    Qo_max=s_o[max_o]
    Qs_max=s_s[max_s]
    dif_qpico=((Qo_max-Qs_max)/Qo_max)*100
    return dif_qpico
def eval_pareto2d(vec1,vec2,Nsuperficies=1,mejor1='inf',mejor2='inf'): #Determina las superficies de valores optimos de pareto
    # Descripcion
    'Funcion: pareto2d\n'\
    'Descripcion: Dadas dos variables determina los grupos pertenecientes a las diferentes superficies optimas\n'\
    'Parametros Obligatorios:.\n'\
    '	-vec1: vector con el desempeno 1.\n'\
    '	-vec2: vector con el desempeno 2.\n'\
    'Parametros Opcionales:.\n'\
    '	-Nsuperficies: cantidad de superficies optimas que se van a encontrar (menor a la cantidad de datos), defecto Nsuperficies=1.\n'\
    '	-mejor1: si mejor1=inf (defecto) valores bajos buen desempeno, si mejor1=sup valores altos buen desempeno.\n'\
    '	-mejor2: si mejor1=inf (defecto) valores bajos buen desempeno, si mejor2=sup valores altos buen desempeno.\n'\
    #Si la medida de desempeno mejora a medida que es mayor la cambia
    if mejor1=='sup':
	vec1=vec1[:]*-1.0
    if mejor2=='sup':
	vec2=vec2[:]*-1.0
    #Crea la lista vacia de superficies y por primera vez el vector de donde se tomaran los que son optimos
    superficies=list()
    evaluar=np.zeros((len(vec1),1),dtype=int)
    #crea la lista de numpy donde va a meter los indices de quienes son optimos y atriz vacia de NxN
    superficie=np.array(1,dtype=int)
    #Itera desde la posicion 1 hasta la ultima
    for j in range(len(vec1)):
	#Evalua si esa columna ya fue eliminada o no
	k=j+1
	while evaluar[j]==0 and k<len(vec1):
	    #Compara y asigna un valor a la matriz de acuerdo a una de las 3 reglas
	    if vec1[j]<vec1[k] and vec2[j]<vec2[k]: #Regla 2: j supera a k (j es posible pareto)	
		evaluar[k]=1 #no se evalua esa posicion k (ha sido eliminada como un posible pareto)
	    elif vec1[k]<vec1[j] and vec2[k]<vec2[j]: #Regla 3: k supera a j(j no es pareto)	
		evaluar[j]=1
	    k=k+1
    #Determina finalmente quienes pertenecen al conjunto N de optimos
    for j in range(len(evaluar)):
	if evaluar[j]==0:
	    superficie=np.append(superficie,j) 
    #Agrega ese conjunto a la lista de superficies y hace nan esas entradas en los vectores de evalucaion
    superficies.append(superficie[1:])
    evaluar=np.zeros((len(vec1),1),dtype=int)
    evaluar[superficie[1:]]=1 #Ya no tomara en cuenta a quienes fueron elegidos para un grupo
    #Una vez termina de evaluar la cantidad de grupos entrge la lista con los indices de los que son optimos
    return superficies

#---------------------------------------------------------------------------------------------------
#Caudales maximos 
#---------------------------------------------------------------------------------------------------
def Rain2IDF(data,dates,duraciones=[1,2,4,6,8,12,24],Tr=[2.33,5,10,25,50,75,100],
	figura='no',ruta='no'):
	Fechas=pandas.to_datetime(dates)
	Rain=pandas.Series(data,index=Fechas)
	Anos=Rain.resample('A',how='sum').index.year
	#Obtiene los maximos para cada duracion
	Duracion=[]
	for d in duraciones:	
		Rain1=Rain.resample(str(d)+'H',how='sum')
		intensidad=[]
		for y in Anos:
			intensidad.append(Rain1[str(y)].max())
		Duracion.append(np.array(intensidad)/float(d))	
	Intensidad=np.array(Duracion)
	#construye las intensidades
	param=[stat.gumbel_r.fit(d) for d in Intensidad]
	P=np.array([1-1/float(i) for i in Tr])
	IntAdjusted=[]
	for t in P:
		I=[]
		for p in param:
			gum=stat.gumbel_r.freeze(p)
			I.append(gum.ppf(t)[0])
		IntAdjusted.append(I)
	IntAdjusted=np.array(IntAdjusted)
	#Ajuste numerico
	Text=np.log(np.array([Tr for i in range(len(duraciones))]).T.reshape(len(duraciones)*len(Tr)))
	Dext=np.log(np.array([duraciones for i in range(len(Tr))]).reshape(len(duraciones)*len(Tr)))
	Ones=np.ones(Dext.shape)
	X=np.vstack((Ones,Dext,Text)).T
	Y=np.log(np.reshape(IntAdjusted,(IntAdjusted.size)))
	param,residual=np.linalg.lstsq(X,Y)[:2]
	k=np.exp(param[0]); n=-1*param[1]; m=param[2]
	R2 = 1 - residual / (Y.size * Y.var())
	#construccion de las intensidades simuladas
	IntS=[]
	for t in Tr:
		Is=[]
		for d in duraciones:
			Is.append((k*t**m)/(d**n))
		IntS.append(Is)
	IntSim=np.array(IntS)
	if figura=='si':	
		fig=pl.figure(edgecolor='w',facecolor='w')
		ax=fig.add_subplot(111)
		colores=['b','r','k','c','y','m','g']
		for i,col,t in zip(IntSim,colores,Tr):
			ax.plot(duraciones,i,c=col,label='Tr: '+str(t))
		ktext='%.2f' % k
		mtext='%.2f' % m
		ntext='%.2f' % n
		ax.text(0.2,0.8,
			r'$I_{T_{r},D}=\frac{'+ktext+'T_{r}^{'+mtext+'}}{D^{'+ntext+'}}$',
			transform=ax.transAxes,
			size=24)
		ax.grid(True)
		ax.set_xlabel('Duracion $[h]$',size=14)
		ax.set_ylabel('Intensidad $[mm/h]$',size=14)
		ax.legend(loc=0)
		if ruta<>None:
			pl.savefig(ruta,bbox_inches='tight')
		pl.show()
	return IntSim,k,m,n,R2
def DesingStorm(IntTr,Dur,CN=None,plot='no',ruta=None,Tr=None):	
	CurvaHuff=np.array([0.18,0.47,0.65,0.74,0.80,0.85,0.89,0.92,0.94,1.0])
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
def QmaxGardex(Tr,Area,Tc,Qp,Tp=2.33,alfaP=10,ecuacion=1):
	if ecuacion==1:
		alfaQ=Area*alfaP/86.4
	elif ecuacion==2:
		alfaQ=Area*alfaP/(3.6*Tc)
	Qtr=[]
	for t in Tr:
		Qtr.append(Qp+alfaQ*np.log(t/Tp))
	return np.array(Qtr)
def QmaxRacional(IntTr,A,Coef):
	Qtr=IntTr*A*Coef/3.6
	return Qtr
def HidrogramaSnyder(Area,Tc,Cp,Fc=2.9,pequena='si'):
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
def HidrogramaWilliams(Area,LongCuenca,PendCauce,Tc):
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
def HidrogramaSCS(Area,Tc,N=25):
	#Parametros fijos
	ttp=np.array([0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.1,
	1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2.0,2.2,2.4,2.6,2.8,3.0,3.2,
	3.4,3.6,3.8,4.0,4.5,5.0])
	UUp=np.array([0.0,0.03,0.1,0.19,0.31,0.47,0.66,0.82,0.93,0.99,1.0,0.99,0.93,0.86,
	0.78,0.68,0.56,0.46,0.39,0.33,0.28,0.21,0.15,0.11,0.08,0.06,0.04,0.03,0.02,0.02,0.01,
	0.01,0.0])
	#Parametros de la hidrografa
	Trezago=(3.0/5.0)*Tc
	Tlluvia=0.133*Tc
	Tpico=Tlluvia/2+Trezago
	Up=484*Area*0.386/Tpico
	Upm=Up*0.3048**3/25.4
	#Obtiene el HU
	HU=[]
	for t,u in zip(ttp,UUp):
		uTemp=u*Up
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
def Convolution(Tiempo,Qhu,lluvEfect):
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
	
