py.install_sources(
  ['__init__.py', 'wmf.py'],
  subdir : 'wmf')

# f2py generated wrappers; the output lists match exactly what f2py emits
# for each source (cuencas.f90 also produces F77-style wrappers)
cu_gen = custom_target('cu-f2py',
  input : ['cuencas.f90'],
  output : ['cumodule.c', 'cu-f2pywrappers.f', 'cu-f2pywrappers2.f90'],
  command : [py, '-m', 'numpy.f2py', '@INPUT@', '-m', 'cu', '--lower',
             '--build-dir', '@OUTDIR@'])

models_gen = custom_target('models-f2py',
  input : ['modelosv2.f90'],
  output : ['modelsmodule.c', 'models-f2pywrappers2.f90'],
  command : [py, '-m', 'numpy.f2py', '@INPUT@', '-m', 'models', '--lower',
             '--build-dir', '@OUTDIR@'])

py.extension_module('cu',
  ['cuencas.f90', cu_gen, fortranobject_c],
  c_args : np_inc_args,
  dependencies : py_dep,
  fortran_args : fortran_args,
  link_args : fortran_link_args,
  link_language : 'fortran',
  subdir : 'wmf',
  install : true)

py.extension_module('models',
  ['modelosv2.f90', models_gen, fortranobject_c],
  c_args : np_inc_args,
  dependencies : py_dep,
  fortran_args : fortran_args,
  link_args : fortran_link_args,
  link_language : 'fortran',
  subdir : 'wmf',
  install : true)
