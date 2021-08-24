#!/usr/bin/env python

import os, sys
import logging
import subprocess
import shutil

from setuptools import setup, find_packages
from setuptools.command.build_py import build_py
from setuptools.extension import Extension

import numpy as np
from Cython.Build import cythonize
import cmake


WORKING_PATH = os.getcwd()

class Build_VTL(build_py):
    """Build VocalTractLab API"""
    def run( self ):
        print( 'Building VocalTractLab API using cmake:' )
        os.chdir( 'PyVTL/src' )
        #with TemporaryDirectory() as tmpdir:
        #    os.chdir(tmpdir)
        subprocess.check_call( [ 'cmake', '.' ] )
        subprocess.check_call( [ 'cmake', '--build', '.' ] )
        vtl_api_name = 'VocalTractLabApi'
        if sys.platform == 'win32':
            file_extension = '.dll'
        else:
            file_extension = '.so'
        shutil.move( os.path.join( 'Debug', vtl_api_name + file_extension ), os.path.join( WORKING_PATH, 'PyVTL' ) )
        shutil.move( os.path.join( 'Debug', vtl_api_name + '.lib' ), os.path.join( WORKING_PATH, 'PyVTL' ) )
        #shutil.move( os.path.join( '', vtl_api_name + '.h' ), os.path.join( WORKING_PATH, 'PyVTL' ) )
        print( ' chir dir: ' )
        print( os.listdir( os.getcwd() ) )
        os.chdir( WORKING_PATH )
        print( 'working dir:' )
        print( os.listdir( os.getcwd() ) )
        print( 'PyVTL dir:' )
        print( os.listdir( os.getcwd()+'/PyVTL' ) )
        build_py.run( self )







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
#try:
#    from Cython.Build import cythonize
#except:
#    log.critical(
#        'Cython.Build.cythonize not found. '
#        'Cython is required to build from a repo.')
#    sys.exit(1)


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
        Extension( 'PyVTL.VocalTractLabApi',
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
    cmdclass = {"build_py": Build_VTL},
    include_dirs=np.get_include(),
    packages=find_packages(),
    package_dir={'PyVTL': 'PyVTL'},
    #package_data= {'PyVTL': ['API/*', 'Models/*', 'Speaker/*', 'Data/*']},
    package_data= {'PyVTL': [ os.path.join( WORKING_PATH, 'PyVTL/Speaker/*'), os.path.join( WORKING_PATH,'./PyVTL/*' ) ]},
    include_package_data = True,
    install_requires=DEPENDENCIES,
    #zip_safe= False,
)






setup(**setup_args)