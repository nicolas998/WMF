#!/usr/bin/env python
from setuptools import setup, Extension
import numpy
import os

# Define Fortran extensions using numpy's f2py
from numpy.distutils.core import Extension as NumpyExtension
from numpy.distutils.core import setup as numpy_setup



ext1 = Extension(name = 'cu',
                 sources = ['wmf/cuencas.f90'])
ext2 = Extension(name = 'models',
                 sources = ['wmf/modelosv2.f90'])

setup(
    name='wmf',
    version='1.0',
    author='Nicolas Velasquez G',
    author_email='nicolas.velasquezgiron@gmail.com',    
    packages=['wmf'],
    package_data={'wmf':['cu.so','models.so']},
    url='https://github.com/nicolas998/WMF.git',
    license='LICENSE.txt',
    description='Watershed Modelling Framework',
    long_description=open('README.txt').read(),
    install_requires=[ ],
    ext_modules=[ext1, ext2],
	)
