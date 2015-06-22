===========
wmf
===========

Este es un modulo orientado al trabajo con cuencas hidrograficas 
y a la ejecucion de modelos hidrologicos, ademas contiene herramientas
de vision de datos y analisis geomorfologico de cuencas.

La manera de importar el modulo es:

    #!/usr/bin/env python

    from wmf.modelacion import cuencas
    from wmf.modelacion import modelos 
    from wmf import utils
   
modulos
=========

Contiene modulos provenientes de fortran, compilados mediante 
gfortran e integrados a python mediante f2py.

cuencas
-------

Contiene las herramientas para determinar la cuenca y varios 
parametros asociados a la misma, cuenta igualmente con herramientas
para determinar el cauce y algunas propiedades geomorfologicas.


modelos
-------

Modulo orientado a la preparación y ejecución del modelo SHIA o bien
modelos que posteriormente se agreguen al sistema.

utils
=====

Se compone de diferentes funciones orientadas a la lectura y escritura 
de datos propios del modelo.
