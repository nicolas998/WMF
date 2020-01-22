___
# wmf
___

WMF (Watershed Modeling Framework) is a module to design to work with hydrographic watersheds and the execution of hydrological models. Additionally, it contains tools for data visualization, analysis of variables, and geomorphological analysis.

In a Linux machine you just need you can install the module by typing:

>- sudo python3 setup.py install 

or 

>- python3 setup.py install --user

To import the module:

>- from wmf import wmf


## Modules
___

WMF is the result of two Fortran modules and one python script.

**cuencas.f90**

This module has the tools to extract streams and watersheds from a given DIR and DEM maps.  Also, it has multiple functions to perform different tasks based on the topology of the watershed. Some of the tasks are:

- Obtain the area, slope, TI, HAND, Width function, Horton Order.
- Extract the network, identify hillslopes and links. 
- Convert raster to the topology of the watershed.  It can also convert variables to hills or vice-versa.

**modelos.f90**

This module has the TETIS model (Velez, 2001) written from scratch. It also contains several sub-models and has functions to write and read binary data related to the input and output variables of the hydrological model.  The sub-models are:

- HydroFlash: A flash flood 1D hydraulic model to obtain flood plains during execution based on the DEM and the results of the hydrological model.

- SHIA-landslide: An adaptation of the landslide model developed by Aristizabal (2014).

- SED: An adaptation of the sediment production model developed for the CASC2D-SED model.

**wmf.py**

This is the base script that merges cuencas.f90 and modelos.f90.  It has defined several classes such as the **SimuBasin** class that could be considered the heart of WMF.  In this module, we have many functionalities done as an interface to the Fortran modules.

## Requirements:

>- A Fortran compiler such as **gfortran**, also the user must have the python3-dev tools.
>- Python 3.6 
>- A Unix machine
>- You must have these packages.
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

I also have tried WMF in **Google Colab** just type:

!pip install git+https://github.com/nicolas998/WMF.git

