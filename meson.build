project('wmf', 'c', 'fortran',
  version : '1.1.0',
  license : 'GPL-3.0-or-later',
  meson_version : '>= 1.1.0',
  default_options : ['buildtype=release', 'fortran_std=legacy'])

py = import('python').find_installation(pure: false)
py_dep = py.dependency()

ff = meson.get_compiler('fortran')

# The Fortran sources predate the strict checks of modern gfortran
fortran_args = []
foreach arg : ['-fallow-argument-mismatch', '-fallow-invalid-boz']
  if ff.has_argument(arg)
    fortran_args += arg
  endif
endforeach

# On Windows (MinGW gfortran) link the GNU runtimes statically so the
# extensions do not depend on MinGW DLLs being on the PATH at import time
fortran_link_args = []
if host_machine.system() == 'windows' and ff.get_id() == 'gcc'
  foreach arg : ['-static-libgfortran', '-static-libgcc', '-static-libquadmath']
    if ff.has_link_argument(arg)
      fortran_link_args += arg
    endif
  endforeach
  # Also satisfy the pthread dependency from the static libwinpthread.a
  fortran_link_args += ['-Wl,-Bstatic', '-lwinpthread', '-Wl,-Bdynamic']
endif

incdir_numpy = run_command(py,
  ['-c', 'import numpy; print(numpy.get_include())'],
  check : true).stdout().strip()
incdir_f2py = run_command(py,
  ['-c', 'import numpy.f2py; print(numpy.f2py.get_include())'],
  check : true).stdout().strip()
# Plain -I flags instead of include_directories() so the build also works
# when numpy lives in a venv inside the source tree
np_inc_args = ['-I' + incdir_numpy, '-I' + incdir_f2py]
fortranobject_c = incdir_f2py / 'fortranobject.c'

subdir('wmf')
