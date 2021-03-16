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
package_data= {'PyVTL': ['API/*', 'Models/*', 'Speaker/*']},
include_package_data = True,
url='http://pypi.python.org/pypi/PackageName/',
license='LICENSE.txt',
description='An awesome package that does something3',
long_description=open('README.md').read(),
install_requires=[],
#zip_safe= False,
)