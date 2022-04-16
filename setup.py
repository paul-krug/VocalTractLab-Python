#!/usr/bin/env python
DOCLINES = 'This Python library provides an implementation of the articulatory synthesizer VocalTractLab.\n\
\n\
Core features are:\
- a one dimensional aero-acoustic simulation of vocal tract dynamics using a\
  vocal tract model based on MRI scans of a human vocal tract\
- synthesis of artificial speech through high-level gestural control or low-\
  level motor control of individual articulators through dynamic targets\
- extensive visualization options\
- estimation of articulatory targets from audio files directly (pitch target-\
  estimation) or from arbitrary data (e.g. articulatory measurements)\
- phoneme-to-speech functionality and a basic text-to-speech pipeline\
- and much more\
\
Besides the scientific purpose, this module is well suited for any production\
environments that need to access the VocalTractLab-backend.\
VocalTractLab and its Python module are licensed under GPL-3.0.'


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

class Build_Target_Optimizer( build_py ):
    """Build TargetOptimizer-Backend API"""
    def run( self ):
        print( 'Building Target_Optimizer-Backend using cmake:' )
        os.chdir( 'VocalTractLab/src/targetoptimizer-backend' )
        #with TemporaryDirectory() as tmpdir:
        #    os.chdir(tmpdir)
        subprocess.check_call( [ 'cmake', '.' ] )
        subprocess.check_call( [ 'cmake', '--build', '.', '--config', 'Release' ] )
        api_name = 'TargetOptimizerApi'
        if sys.platform == 'win32':
            file_extension = '.dll'
            shutil.move( os.path.join( 'Release', api_name + file_extension ), os.path.join( WORKING_PATH, 'VocalTractLab' ) )
            shutil.move( os.path.join( 'release', api_name + '.lib' ), os.path.join( WORKING_PATH, 'VocalTractLab' ) )
        else:
            file_extension = '.so'
            shutil.move( 'lib' + api_name + file_extension, os.path.join( WORKING_PATH, 'VocalTractLab' ) )
        #print( os.listdir( os.getcwd() ) )
        #print( os.listdir( os.path.join( os.getcwd(), 'CMakeFiles' ) ) )
        shutil.move( os.path.join( '', api_name + '.h' ), os.path.join( WORKING_PATH, 'VocalTractLab' ) )
        shutil.move( os.path.join( '', 'Data.h' ), os.path.join( WORKING_PATH, 'VocalTractLab' ) )
        os.chdir( WORKING_PATH )
        #build_py.run( self )

class Build_VTL( build_py ):
    """Build VocalTractLab-Backend API"""
    def run( self ):
        print( 'Building VocalTractLab-Backend using cmake:' )
        os.chdir( 'VocalTractLab/src/vocaltractlab-backend' )
        #with TemporaryDirectory() as tmpdir:
        #    os.chdir(tmpdir)
        subprocess.check_call( [ 'cmake', '.' ] )
        subprocess.check_call( [ 'cmake', '--build', '.', '--config', 'Release' ] )
        api_name = 'VocalTractLabApi'
        if sys.platform == 'win32':
            file_extension = '.dll'
            shutil.move( os.path.join( 'Release', api_name + file_extension ), os.path.join( WORKING_PATH, 'VocalTractLab' ) )
            shutil.move( os.path.join( 'Release', api_name + '.lib' ), os.path.join( WORKING_PATH, 'VocalTractLab' ) )
        else:
            file_extension = '.so'
            shutil.move( 'lib' + api_name + file_extension, os.path.join( WORKING_PATH, 'VocalTractLab' ) )
        shutil.move( os.path.join( '', api_name + '.h' ), os.path.join( WORKING_PATH, 'VocalTractLab' ) )
        #print( ' chir dir: ' )
        #print( os.listdir( os.getcwd() ) )
        os.chdir( WORKING_PATH )
        #print( 'working dir:' )
        #print( os.listdir( os.getcwd() ) )
        #print( 'VocalTractLab dir:' )
        #print( os.listdir( os.getcwd()+'/VocalTractLab' ) )
        #build_py.run( self )

class Build_Backends( build_py ):
    def run(self):
        #self.run_command( 'build_target_optimizer' )
        self.run_command( 'build_vtl' )
        build_py.run(self)




# Set up the logging environment
logging.basicConfig()
log = logging.getLogger()

# Handle the -W all flag
if 'all' in sys.warnoptions:
    log.level = logging.DEBUG

# Get version from the VocalTractLab module
with open('VocalTractLab/__init__.py') as f:
    for line in f:
        if line.find('__version__') >= 0:
            version = line.split('=')[1].strip()
            version = version.strip('"')
            version = version.strip("'")
            continue


# Build extension modules 
EXT_MODULES = cythonize(
    [
        #Extension( 'VocalTractLab.target_estimation',
        #      ['VocalTractLab/target_estimation.pyx'],
        #      language="c++",
        #      libraries=['VocalTractLab/TargetOptimizerApi'],
        #      library_dirs=['.'],#, './src/', './src/targetoptimizer-backend/'],
        #      include_dirs=[np.get_include()],#, './src/', './src/targetoptimizer-backend/']
        #      ),
        Extension( 'VocalTractLab.VocalTractLabApi',
              ['VocalTractLab/VocalTractLabApi.pyx'],
              language="c",
              libraries=['VocalTractLab/VocalTractLabApi'],
              library_dirs=['.'],
              include_dirs=[np.get_include()]
              ),
    ]
)

# Dependencies
DEPENDENCIES = [
    'librosa>=0.8.1',
    'matplotlib>=3.4.3',
    'numpy>=1.22.0',
    'pandas>=1.3.2',
    'praat-parselmouth>=0.4.0',
    'tqdm>=4.62.1',
]


CLASSIFIERS = """\
Development Status :: 3 - Alpha
Intended Audience :: Science/Research
Intended Audience :: Developers
License :: OSI Approved :: GNU General Public License v3 (GPLv3)
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


cmdclass = dict( #build_target_optimizer = Build_Target_Optimizer,
                 build_vtl = Build_VTL,
                 build_py = Build_Backends,
                 )
#cmdclass['build_target_optimizer'] = Build_Target_Optimizer
#cmdclass['build_vtl'] = Build_VTL
#cmdclass['build_py'] = my_build

setup_args = dict(
    name='VocalTractLab',
    version=version,
    description='Articulatory (text-to-) speech synthesis for Python',
    long_description= DOCLINES,
    url='https://github.com/paul-krug/VocalTractLab',
    #download_url='https://github.com/paul-krug/VocalTractLab/archive/v_0.3.tar.gz',
    author='Paul Krug',
    author_email='paul_konstantin.krug@tu-dresden.de',
    license='GPL-3.0',
    classifiers = [_f for _f in CLASSIFIERS.split('\n') if _f],
    keywords=[ 'text-to-speech', 'speech synthesis', 'articulatory synthesis', 'vocal tract' ],
    ext_modules=EXT_MODULES,
    cmdclass = cmdclass,
    include_dirs=np.get_include(),
    packages=find_packages(),
    package_dir={'VocalTractLab': 'VocalTractLab'},
    #package_data= {'VocalTractLab': ['API/*', 'Models/*', 'Speaker/*', 'Data/*']},
    package_data= {'VocalTractLab': [ os.path.join( WORKING_PATH, 'VocalTractLab/speaker/*'),
                              os.path.join( WORKING_PATH, 'VocalTractLab/data/*'), 
                              os.path.join( WORKING_PATH, 'VocalTractLab/data/dictionaries/phonecodes/*'),
                              os.path.join( WORKING_PATH, 'VocalTractLab/src/vocaltractlab-backend/*' ),
                              os.path.join( WORKING_PATH, 'VocalTractLab/src/targetoptimizer-backend/*' ),
                              os.path.join( WORKING_PATH, 'VocalTractLab/src/targetoptimizer-backend/dlib/*' ),
                              os.path.join( WORKING_PATH,'./VocalTractLab/*' ) ]},
    include_package_data = True,
    install_requires=DEPENDENCIES,
    use_scm_version=True,
    setup_requires=['setuptools_scm'],
    #zip_safe= False,
)






setup(**setup_args)