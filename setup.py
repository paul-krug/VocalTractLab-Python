#!/usr/bin/env python
'''
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
'''

from setuptools import setup

install_requires = [
    #'beautifulsoup4>=4.11.1',
    'matplotlib',
    'numpy',
    'pandas',
    'tools-mp',
    'vocaltractlab-cython',
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

setup_args = dict(
    name='VocalTractLab',
    version='0.9.0',
    description='High-performance articulatory speech synthesis in Python',
    long_description_content_type='text/markdown',
    long_description=open('README.md').read(),
    #long_description= DOCLINES,
    url='https://github.com/paul-krug/VocalTractLab-Python',
    author='Paul Krug',
    author_email='paul_konstantin.krug@tu-dresden.de',
    license='GPL-3.0',
    #classifiers = [_f for _f in CLASSIFIERS.split('\n') if _f],
    keywords=[
        'text-to-speech',
        'speech synthesis',
        'articulatory synthesis',
        'vocal tract',
        'speech production',
        'vocoder',
        ],
    packages=['vocaltractlab'],
    install_requires=install_requires

    #use_scm_version=True,
    #setup_requires=['setuptools_scm'],
    #zip_safe= False,
)

setup(**setup_args)