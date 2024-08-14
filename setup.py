


from setuptools import setup



install_requires = [
    'numpy',
    'target-approximation>=0.0.4',
    'tools-mp',
    'torch',
    'torchaudio',
    'vocaltractlab-cython==0.0.13',
    ]



#CLASSIFIERS = """\
#Development Status :: 3 - Alpha
#Intended Audience :: Science/Research
#Intended Audience :: Developers
#License :: OSI Approved :: GNU General Public License v3 (GPLv3)
#Programming Language :: C++
#Programming Language :: Python
#Programming Language :: Python :: 3
#Programming Language :: Python :: 3.8
#Programming Language :: Python :: 3.9
#Programming Language :: Python :: Implementation :: CPython
#Topic :: Software Development
#Topic :: Scientific/Engineering
#Typing :: Typed
#Operating System :: Microsoft :: Windows
#Operating System :: POSIX
#Operating System :: Unix
#"""

setup_args = dict(
    name='VocalTractLab',
    version='0.5.2',
    description='High-performance articulatory speech synthesis in Python',
    long_description_content_type='text/markdown',
    long_description=open('README.md').read(),
    url='https://github.com/paul-krug/VocalTractLab-Python',
    author='Paul Krug',
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
)

setup(**setup_args)