import os
import subprocess
from glob import glob
from setuptools import setup, Extension

file_list = [x for x in os.listdir('.') if os.path.splitext(x)[1] in ['.f90', '.h', '.c']]
file_list.remove('pltx.c')  # for future interfacing with qt -- using dummy now
bad_startstring = 'fortrancore'  # generated from the .pyf file, so don't include them a second time!
file_list = [x for x in file_list if not x.startswith(bad_startstring)]

print(
    "first, I'm going to run make, because setup.py doesn't do build order correctly -- if this fails, you need to be sure to call from within the git bash shell"
)
subprocess.call('make')


def build_extension():
    """Build the Fortran extension with f2py if it isn't already present."""
    existing = glob(os.path.join('nlsl', 'fortrancore*.so'))
    if existing:
        print(f"Found prebuilt extension: {existing[0]}")
        return
    cmd = ['f2py', '-c', '-m', 'nlsl.fortrancore', 'pynlsl.pyf'] + file_list
    print('Building Fortran extension with', ' '.join(cmd))
    try:
        subprocess.check_call(cmd)
    except subprocess.CalledProcessError as err:
        print('Warning: f2py failed to build the extension:', err)
        print('Continuing without a compiled fortrancore module.')


build_extension()

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
    ext_modules=[],
)
