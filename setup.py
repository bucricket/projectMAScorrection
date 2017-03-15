#!/usr/bin/env python
import subprocess
import os
import shutil
import wget
import tarfile



def untar(fname, fpath):
    if (fname.endswith('tar.gz') or fname.endswith('tar.bz') or fname.endswith('tar')):
        tar = tarfile.open(fname)
        tar.extractall(path = fpath)
        tar.close()
        os.remove(fname)
        
base = os.getcwd()
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

#======download, untar and move atlases and coefficients=======================
attempts =0
while attempts < 10:
    try:
        wget.download('https://nwpsaf.eu/downloads/emis_data/uw_ir_emis_atlas_hdf5.tar')
        break
    except :
        attempts += 1
untar('uw_ir_emis_atlas_hdf5.tar',base)
source = os.listdir(base)
for files in source:
    if files.endswith('.H5'):
        shutil.move(os.path.join(base,files), os.path.join(rttovEmisPath,files))
        
attempts =0
while attempts < 5:
    try:        
        wget.download('https://nwpsaf.eu/downloads/emis_data/uw_ir_emis_atlas_covariances_hdf5.tar')
        break
    except :
        attempts += 1
untar('uw_ir_emis_atlas_covariances_hdf5.tar',base)
sourcePath = os.path.join(base,'uw_ir_emis_atlas_covariances_hdf5.tar')
source = os.listdir(base)
for files in source:
    if files.endswith('.H5'):
        shutil.move(os.path.join(base,files), os.path.join(rttovEmisPath,files))
        
wget.download('https://nwpsaf.eu/downloads/emis_data/uw_ir_emis_atlas_angcorr_hdf5.tar')
untar('uw_ir_emis_atlas_angcorr_hdf5.tar',base)
sourcePath = os.path.join(base,'uw_ir_emis_atlas_angcorr_hdf5.tar')
source = os.listdir(base)
for files in source:
    if files.endswith('.H5'):
        shutil.move(os.path.join(base,files), os.path.join(rttovEmisPath,files))
#=========BRDF=================================================================
attempts =0
while attempts < 5:
    try:     
        wget.download('https://nwpsaf.eu/site/download/rttov_downloads/brdf_data/cms_brdf_atlas_hdf5.tar')
        break
    except :
        attempts += 1
untar('cms_brdf_atlas_hdf5.tar',base)
source = os.listdir(base)
for files in source:
    if files.endswith('.H5'):
        shutil.move(os.path.join(base,files), os.path.join(rttovBRDFPath,files))

shutil.copyfile(os.path.join(libDir,'rttov_wrapper_f2py.so'),os.path.join(prefix,'lib','python2.7','site-packages','rttov_wrapper_f2py.so'))

shutil.copyfile(os.path.join(processDir,'source','rtcoef_rttov11','rttov7pred54L','rtcoef_landsat_8_tirs.dat'),os.path.join(rttovPath,'rtcoef_landsat_8_tirs.dat'))
#try:
from setuptools import setup
#    setup_kwargs = {'entry_points': {'console_scripts':['pyrttov=pyrttov.__init__:Rttov']}}
#except ImportError:
#    from distutils.core import setup
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

