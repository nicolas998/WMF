#!/usr/bin/env python
import os
from numpy.distutils.core import setup, Extension

setup(
    name='wmf',
    version='0.1.0',
    author='Nicolas Velasquez G',
    author_email='nicolas.velasquezgiron@gmail.com',
    packages=['wmf'],
    package_data={'wmf':['modelacion.so']},
    url='http://pypi.python.org/pypi/wmf/',
    license='LICENSE.txt',
    description='.',
    long_description=open('README.txt').read(),
    install_requires=[
        "numpy >= 1.6.1",
        "osgeo >= 2.7.3",
	"matplotlib >= 0.15",
    ],
)
