from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext
import numpy as np

ext_modules = [
    Extension('VocalTractLabApi',
              ['VocalTractLabApi.pyx'],
              # Note here that the C++ language was specified
              # The default language is C
              language="c",  
              libraries=['VocalTractLabApi'],
              library_dirs=['.'])
    ]

setup(
    name = 'VocalTractLabApi',
    cmdclass = {'build_ext': build_ext},
    ext_modules = ext_modules,
    include_dirs=[np.get_include()]  # This gets all the required Numpy core files
)