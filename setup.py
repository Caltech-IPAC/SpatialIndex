import os

from setuptools import setup
from distutils.extension import Extension
from Cython.Build import cythonize

sptind_extension = Extension(
    name="spatial_index",
    sources = ["spatial_index.pyx"],
    extra_objects = ["lib/sptQueryLib.o"],
    libraries = ["htm"],
    library_dirs = ["lib"],
    include_dirs = ["lib/include"]
)
setup(
    name = "spatial_index",
    version = '1.1.0',
    author = 'John Good',
    author_email = 'jcg@ipac.caltech.edu',
    description = 'Library for generating DBMS spatial index constraints.',
    long_description = open('README.md').read(),
    license = 'LICENSE.txt',
    keywords = 'astronomy database spatial-index',
    url = 'https://github.com/Caltech-IPAC/SpatialIndex',
    packages = ['SpatialIndex'],
    ext_modules = cythonize([sptind_extension])
)
