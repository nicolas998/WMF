#!/usr/bin/env python
import os
from numpy.distutils.core import setup, Extension

ext1 = Extension(name = 'cu',
                 sources = ['wmf/cuencas.f90'])
ext2 = Extension(name = 'models',
                 sources = ['wmf/modelosv2.f90'])

setup(
    name='wmf',
    version='0.2.11',
    author='Nicolas Velasquez G',
    author_email='nicolas.velasquezgiron@gmail.com',    
    packages=['wmf'],
    package_data={'wmf':['cu.so','models.so']},
    url='https://github.com/nicolas998/WMF.git',
    license='LICENSE.txt',
    description='.',
    long_description=open('README.txt').read(),
    install_requires=[ ],
    ext_modules=[ext1, ext2],
	)
