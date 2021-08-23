#!/usr/bin/env python

import logging
import sys
import pprint
from setuptools import setup, find_packages
from setuptools.extension import Extension
import numpy as np

# Set up the logging environment
logging.basicConfig()
log = logging.getLogger()

# Handle the -W all flag
if 'all' in sys.warnoptions:
    log.level = logging.DEBUG

# Parse the verison from the PyVTL module
with open('PyVTL/__init__.py') as f:
    for line in f:
        if line.find('__version__') >= 0:
            version = line.split('=')[1].strip()
            version = version.strip('"')
            version = version.strip("'")
            continue

#with open('VERSION.txt', 'w') as f:
#    f.write(version)

# Use Cython if available
try:
    from Cython.Build import cythonize
except:
    log.critical(
        'Cython.Build.cythonize not found. '
        'Cython is required to build from a repo.')
    sys.exit(1)


# Extension options
#include_dirs = []
#try:
#    import numpy as np
#    include_dirs.append(numpy.get_include())
#except ImportError:
#    log.critical('Numpy and its headers are required to run setup(). Exiting')
#    sys.exit(1)



# Build extension modules 
EXT_MODULES = cythonize(
    [
        Extension( 'VocalTractLabApi',
              ['PyVTL/VocalTractLabApi.pyx'],
              language="c",  
              libraries=['PyVTL/VocalTractLabApi'],
              library_dirs=['.'],
              include_dirs=[np.get_include()]
              )
    ]
)

# Dependencies
DEPENDENCIES = [
    'librosa>=0.8.1',
    'matplotlib>=3.4.3',
    'numpy>=1.21.2',
    'pandas>=1.3.2',
    'tqdm>=4.62.1',
]


CLASSIFIERS = """\
Development Status :: 3 - Alpha
Intended Audience :: Science/Research
Intended Audience :: Developers
License :: OSI Approved :: GPL-3.0 License
Programming Language :: C++
Programming Language :: Python
Programming Language :: Python :: 3
Programming Language :: Python :: 3.8
Programming Language :: Python :: 3.9
Programming Language :: Python :: 3.10
Programming Language :: Python :: Implementation :: CPython
Topic :: Software Development
Topic :: Scientific/Engineering
Typing :: Typed
Operating System :: Microsoft :: Windows
Operating System :: POSIX
Operating System :: Unix
"""


setup_args = dict(
    name='PyVTL',
    version=version,
    description='Articulatory (text-to-) speech synthesis for Python',
    long_description=open('README.md').read(),
    url='http://pypi.python.org/pypi/PackageName/',
    author='Paul Konstantin Krug',
    author_email='paul_konstantin.krug@tu-dresden.de',
    license='GPL-3.0',
    classifiers = [_f for _f in CLASSIFIERS.split('\n') if _f],
    keywords=[ 'text-to-speech', 'speech synthesis', 'articulatory synthesis', 'vocal tract' ],
    ext_modules=EXT_MODULES,
    packages=find_packages(),
    package_dir={'PyVTL': 'PyVTL'},
    #package_data= {'PyVTL': ['API/*', 'Models/*', 'Speaker/*', 'Data/*']},
    include_package_data = True,
    install_requires=DEPENDENCIES,
    #zip_safe= False,
)






setup(**setup_args)



















from setuptools import setup

setup(
name='PyVTL',
version='0.1.1',
author='Paul Konstantin Krug',
author_email='paul_konstantin.krug@tu-dresden.de',
packages=['PyVTL'],
#scripts=['PyVTL/Models/PCA_nc_10.joblib'],
#data_files = ['./PyVTL/Models/PCA_nc_10.joblib'],
package_dir={'PyVTL': 'PyVTL'},
package_data= {'PyVTL': ['API/*', 'Models/*', 'Speaker/*', 'Data/*']},
include_package_data = True,
url='http://pypi.python.org/pypi/PackageName/',
license='LICENSE.txt',
description='Articulatory (text-to-) speech synthesis for Python',
long_description=open('README.md').read(),
install_requires=[],
#zip_safe= False,
)