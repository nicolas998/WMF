___
# wmf
___

Este es un modulo orientado al trabajo con cuencas hidrograficas 
y a la ejecucion de modelos hidrologicos, ademas contiene herramientas
de vision de datos y analisis geomorfologico de cuencas.

La manera de importar el modulo es:

    #!/usr/bin/env python
	>- from wmf import wmf 
	
   
## Modulos
___

Contiene modulos provenientes de fortran, compilados mediante 
gfortran e integrados a python mediante f2py.

**cuencas**

Contiene las herramientas para determinar la cuenca y varios 
parametros asociados a la misma, cuenta igualmente con herramientas
para determinar el cauce y algunas propiedades geomorfologicas.


**modelos**

Modulo orientado a la preparación y ejecución del modelo SHIA o bien
modelos que posteriormente se agreguen al sistema.

## wmf.py
___

Reune todo el código y las funciones definicas en **cuencas**
y en **modelos**

## Instalación:
___

Para instalar el paquete se deben seguir los siguientes pasos:

>1. Anclarse al repositorio: https://github.com/nicolas998/WMF.git
o bien descargar el archivo **WMF-master.zip**.
>2. En el caso de descargar el .zip este debe ser descomprimido, 
luego se debe mover la terminal hasta el directorio 
**cd ~/wmf**
>3. se instala el paquete: **sudo python setup.py install**

### Requerimientos:

>- Compilador de fortran (gfortran probado).
>- Python 2.7
>- Sistema Unix (puede compilar en Windows con Cwing)
>- Debe tener las siguiente dependencias de python.
	- numpy 
	- glob
	- mpl_toolkits.basemap  	
	- netCDF4
	- osgeo
	- gdal
	- scipy
	- os
	- pandas
	- datetime
	- matplotlib
