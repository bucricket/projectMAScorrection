#!/usr/bin/env python
import subprocess
import os
import shutil

prefix  = os.environ.get('PREFIX')
processDi = os.path.abspath(os.path.join(prefix,os.pardir))
processDir = os.path.join(processDi,'work')
binDir = os.path.join(prefix,'bin')
libDir = os.path.join(processDir,'source','lib')
srcDir = os.path.join(processDir,'source','src')

#=====write Makefile.local===========
makeFilename = os.path.join(processDir,'source','build',"Makefile.local")

file = open(makeFilename,"w") 
 
file.write("HDF5_PREFIX  = %s\n" % prefix) 
file.write("FFLAGS_HDF5  = -D_RTTOV_HDF $(FFLAG_MOD)$(HDF5_PREFIX)/include\n") 
file.write("LDFLAGS_HDF5 = -L$(HDF5_PREFIX)/lib -lhdf5hl_fortran -lhdf5_hl -lhdf5_fortran -lhdf5\n") 
file.write("FFLAGS_EXTERN  = $(FFLAGS_NETCDF)  $(FFLAGS_HDF5)  $(FFLAGS_DRHOOK)\n") 
file.write("LDFLAGS_EXTERN = $(LDFLAGS_NETCDF) $(LDFLAGS_HDF5) $(LDFLAGS_DRHOOK)")
file.close() 


#=====compile rttov=================
os.chdir(srcDir)
subprocess.call(["../build/rttov_compile.sh"])

#====moving shared library to bin ===========

# get Anaconda root location
p = subprocess.Popen(["conda", "info", "--root"],stdout=subprocess.PIPE)
out = p.communicate()
condaPath = out[0][:-1]

shutil.copyfile(os.path.join(libDir,'rttov_wrapper_f2py.so'),os.path.join(processDir,'pyrttov','rttov_wrapper_f2py.so'))

try:
    from setuptools import setup
#    setup_kwargs = {'entry_points': {'console_scripts':['pyrttov=pyrttov.__init__:Rttov']}}
except ImportError:
    from distutils.core import setup
#    setup_kwargs = {'scripts': ['bin/pyrttov']}
#    
#from pyrttov import __version__

setup(
    name="pyrttov",
    version="0.1.0",
    description="pythonic wrapper for rttov",
    author="Mitchell Schull",
    author_email="mitch.schull@noaa.gov",
    packages= ['pyrttov'],
    package_data = {'pyrttov':['rttov_wrapper_f2py.so']},
    platforms='Posix; MacOS X; Windows',
    license='BSD 3-Clause',
    classifiers=[
        'Development Status :: 2 - Pre-Alpha',
        'Intended Audience :: Developers',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: MIT License',
        'Programming Language :: Python :: 2',
        # Uses dictionary comprehensions ==> 2.7 only
        'Programming Language :: Python :: 2.7',
        'Topic :: Scientific/Engineering :: GIS',
    ],  
#    **setup_kwargs
)

