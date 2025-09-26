# meson.build (root)
project('wmf', 'fortran', version: '1.3.0', default_options: ['warning_level=1'])

py = import('python').find_installation(pure: false)

# Build the Fortran extensions *inside* the package namespace:
py.extension_module('wmf.cu',     'wmf/cuencas.f90')
py.extension_module('wmf.models', 'wmf/modelosv2.f90')

# Install your pure-Python sources (exclude the .f90 files)
py.install_sources(
  ['wmf/__init__.py', 'wmf/wmf.py'],
  subdir: 'wmf',
  pure: true
)