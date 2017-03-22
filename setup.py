#!/usr/bin/env python
import subprocess
import os
import shutil
import tarfile
from setuptools import setup



def untar(fname, fpath):
    if fname.endswith('tar.gz') or fname.endswith('tar.bz') or fname.endswith('tar'):
        tar = tarfile.open(fname)
        tar.extractall(path=fpath)
        tar.close()
        os.remove(fname)
base = os.getcwd()
prefix = os.environ.get('PREFIX')
processDi = os.path.abspath(os.path.join(prefix,os.pardir))
processDir = os.path.join(processDi,'work')
binDir = os.path.join(prefix,'bin')
libDir = os.path.join(processDir,'source','lib')
srcDir = os.path.join(processDir,'source','src')

#=====write Makefile.local===========
makeFilename = os.path.join(processDir,'source','build',"Makefile.local")

fn = open(makeFilename,"w")

fn.write("HDF5_PREFIX  = %s\n" % prefix)
fn.write("FFLAGS_HDF5  = -D_RTTOV_HDF $(FFLAG_MOD)$(HDF5_PREFIX)/include\n")
fn.write("LDFLAGS_HDF5 = -L$(HDF5_PREFIX)/lib -lhdf5hl_fortran -lhdf5_hl -lhdf5_fortran -lhdf5\n")
fn.write("FFLAGS_EXTERN  = $(FFLAGS_NETCDF)  $(FFLAGS_HDF5)  $(FFLAGS_DRHOOK)\n")
fn.write("LDFLAGS_EXTERN = $(LDFLAGS_NETCDF) $(LDFLAGS_HDF5) $(LDFLAGS_DRHOOK)")
fn.close()


#=====compile rttov=================
rttovPath = os.path.join(prefix,'share','rttov')
os.makedirs(rttovPath)
rttovEmisPath = os.path.join(rttovPath,'emis_data')
os.makedirs(rttovEmisPath)
rttovBRDFPath = os.path.join(rttovPath,'brdf_data')
os.makedirs(rttovBRDFPath)

os.chdir(srcDir)
subprocess.call(["../build/rttov_compile.sh"])

#====moving shared library to bin ===========
p = subprocess.Popen(["conda", "info", "--root"],stdout=subprocess.PIPE)
out = p.communicate()
condaPath = out[0][:-1]

os.chdir(base)

shutil.copyfile(os.path.join(libDir,'rttov_wrapper_f2py.so'),
                os.path.join(prefix,'lib','python2.7',
                             'site-packages','rttov_wrapper_f2py.so'))

shutil.copyfile(os.path.join(processDir,'source','rtcoef_rttov11',
                             'rttov7pred54L','rtcoef_landsat_8_tirs.dat'),
                os.path.join(rttovPath,'rtcoef_landsat_8_tirs.dat'))

setup(
    name="pyrttov",
    version="0.1.0",
    description="pythonic wrapper for rttov",
    author="Mitchell Schull",
    author_email="mitch.schull@noaa.gov",
    packages=['pyrttov'],
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
)

