import os

from setuptools import setup
from distutils.extension import Extension
from Cython.Build import cythonize

objs = []
 
objs.append('SpatialIndex/lib/src/sptQueryLib.o')
objs.append('SpatialIndex/lib/src/tinyhtm/src/htm.o')
objs.append('SpatialIndex/lib/src/tinyhtm/src/tree.o')
objs.append('SpatialIndex/lib/src/tinyhtm/src/common.o')
objs.append('SpatialIndex/lib/src/tinyhtm/src/id_list.o')
objs.append('SpatialIndex/lib/src/tinyhtm/src/tree_count.o')
objs.append('SpatialIndex/lib/src/tinyhtm/src/tree_gen.o')
objs.append('SpatialIndex/lib/src/tinyhtm/src/geometry.o')
objs.append('SpatialIndex/lib/src/tinyhtm/src/select.o')


extensions = Extension(
    sources = ["SpatialIndex/spatial_index.pyx"],
    include_dirs = ["SpatialIndex/lib/include"],
    extra_objects = objs,
)

setup(
    packages = ['SpatialIndex'],

    ext_modules = cythonize(extensions,
                            compiler_directives={'language_level' : '3str'})
)
