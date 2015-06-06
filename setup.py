#from setuptools import setup
import setuptools # I think this is needed for the following
from numpy.distutils.core import Extension,setup
import os
file_list = filter(lambda x: x[x.find('.'):] in ['.f90','.h','.c'],os.listdir('.'))
file_list.remove('pltx.c')# for future interfacing with qt -- using dummy now
bad_startstring = 'fortrancore' # generated from the .pyf file, so don't include them a second time!
file_list = filter(lambda x: not (len(x) > len(bad_startstring) and x[:len(bad_startstring)] == bad_startstring),file_list)

ext_nlsl = Extension(name = 'nlsl.fortrancore',
        sources = ['pynlsl.pyf',]+file_list,
        define_macros = [('ADD_UNDERSCORE',None)],
        )

setup(
    name='NLSL',
    author='NLSL development team',
    version='0.1.0',
    packages=['nlsl'],
    license='LICENSE.md',
    description='the ACERT NLSL package, now wrapped up in python!',
    long_description=open('README.rst').read(),
    install_requires=[
        "sympy",
        "numpy",
        "scipy",
        "matplotlib",
        "tables",
        ],
    ext_modules = [ext_nlsl],
)
