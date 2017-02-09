#! /usr/bin/env python
# -*- coding: utf-8 -*-

"""
This module tests the rttov Fortran
Notice that you should have the PYTHONPATH environment variable
containig the directorys that holds the rttov Python package
This test does:
- load US76 profile
- update latitude/logitude/surface_type in order to be over land
- load MSG3 coefficients for RTTOV11. Solar and Themal channels.
- retrieve surafce emissivity/reflectance for all channels
- run forward model
- display surface emissivity/reflectance
- display radiances
- display brightness temperatures
- display relectances
"""

import numpy
import h5py
import os
from rttov import emissivity, reflectance, profile, \
   option, radiance, getcoefval
from rttov import rttov_gui_f2py

file = "../rttov_tests/us_standard_76.H5"

# Read US76 profile
nprof = profile.getNumberOfProfiles(file)
print "number of profiles found :", nprof

p1 = profile.Profile()
o1 = option.Option()

iprof = 6
p1.loadProfileNumber(file, iprof)

f = h5py.File(file, 'r')
h5 = f['/OPTIONS/']
o1.loadh5(h5)
f.close()

# ADD solar calculations
o1['ADDSOLAR'] = True
o1['ADDINTERP'] = True

# change latitude longitude date and surface type in order to
# put the profile over land in France
# Black body cloud 800hPa cloud fraction 0.40
p1['LATITUDE'] = 45.0
p1['LONGITUDE'] = 5.0
p1['SKIN']['SURFTYPE'] = 0
p1['AZANGLE'] = 45.0
p1['SUNAZANGLE'] = 100.0
p1['SUNZENANGLE'] = 30.0
p1['ZENANGLE'] = 30.0
p1['CFRACTION'] = 0.4
p1['CTP'] = 800.0
o1["CO_DATA"] = False
o1["CO2_DATA"] = False
o1["N2O_DATA"] = False
o1["CH4_DATA"] = False
o1["ADDSOLAR"] = False


# Ecriture du fichier resultat "profile.H5" avec 1 profil et 1 option
# Output p1 to HDF5 file
ofile = "../rttov_tests/profile.H5"
of = h5py.File(ofile, 'w')
p1.saveh5(of, "/PROFILES/0001/")
o1.saveh5(of, '/OPTIONS/')
of.close()


# Read RTTOV v11 coefficinets
nchannels, err = rttov_gui_f2py.rttov_gui_load(
    [0],
    "../rttov_tests/rtcoef_msg_3_seviri.dat", "", "", "")
print "nchannels", nchannels
listchan = getcoefval.rttov_get_coef_val_i1('FF_ORI_CHN')

print listchan

e1 = emissivity.Emissivity()
r1 = reflectance.Reflectance()
e1.setEmissivity(nchannels)
r1.setReflectance(nchannels)
ofile = "../rttov_tests/surface.H5"
of = h5py.File(ofile, 'w')
e1.saveh5(of, '/EMISSIVITY/')
r1.saveh5(of, '/REFLECTANCE/')
of.close()

# Get surface emissivity/reflectance from Atlas
err = rttov_gui_f2py.rttov_gui_get_emisbrdf(os.environ["RTTOV_GUI_EMISS_DIR"],
                                            "../rttov_tests/profile.H5",
                                            "../rttov_tests/surface.H5",
                                            100, 100, 100)


# Read emissivity and Reflectance; print and plot
fileName = "../rttov_tests/surface.H5"

f = h5py.File(fileName, 'r')
h5 = f['/EMISSIVITY/']
e1 = emissivity.Emissivity()
e1.loadh5(h5)
h5 = f['/REFLECTANCE/']
r1 = reflectance.Reflectance()
r1.loadh5(h5)
f.close()

print "\n \n e1 read from surface.H5 "
e1.display()
print "\n"
print "\n \n r1 read from surface.H5 "
r1.display()
print "\n"

# Now RUN RTTOV Direct model
nthreads = 2
err = rttov_gui_f2py.rttov_gui_run("../rttov_tests/profile.H5",
                                   "../rttov_tests/surface.H5",
                                   "../rttov_tests/radr.H5",
                                   "../rttov_tests/kmat.H5",
                                   "../rttov_tests/trns.H5",
                                   "DIRECT", nthreads)

# At present time the matplotlib.pyplot import has conflicts
# with the F90 interface if placed before
# calls to RTTOV routines.
# So do not move this import at the top of the module
import matplotlib.pyplot as plt

# Read emissivity and Reflectance; print and plot
fileName = "../rttov_tests/surface.H5"

f = h5py.File(fileName, 'r')
h5 = f['/EMISSIVITY/']
e1 = emissivity.Emissivity()
e1.loadh5(h5)
h5 = f['/REFLECTANCE/']
r1 = reflectance.Reflectance()
r1.loadh5(h5)
f.close()
print "\n \n e1 read from surface.H5 after RTTOV run "
e1.display()
print "\n"
print "\n \n r1 read from surface.H5 after RTTOV run  "
r1.display()
print "\n"

# Plots emissivity
plt.plot(listchan, e1['EMIS_OUT'], 'b-', marker='^', markersize=4)
plt.plot(listchan, r1['REFL_OUT'], 'r-', marker='v', markersize=4)
plt.xlabel('channel number')
plt.ylabel(e1['EMIS_OUT_ATTRIBUTE']['COMMENT'])
plt.title(file + '  Emissivity/Reflectance _OUT')
# plt.show()


# Read Radiance results
file = '../rttov_tests/radr.H5'
f = h5py.File(file, 'r')
h5 = f['/RADIANCE/']

# rad is a rttov_hdf_mod radiance instance
rad = radiance.Radiance()
rad.loadh5(h5)
f.close()
rad.display()


# Plots brightness temperature, radiance, reflectance

colors = ['b', 'g', 'r', 'c', 'm', 'y', 'k']
i = 0
for v in ['BT', 'BT_CLEAR']:
    color = colors[i]
    i += 1
    plt.plot(listchan, rad[v], color + '-', label=v, marker='o', markersize=8)
plt.xlabel('channel')
plt.ylabel('Brightness Temperatures' +
           '   (' + rad[v + '_ATTRIBUTE']['UNITS'] + ')')
plt.ylim(200, 320)
plt.legend()
plt.title(file)
plt.show()


for v in ['TOTAL', 'CLEAR']:
    color = colors[i]
    i += 1
    plt.plot(listchan, rad[v], color + '-', label=v, marker='o', markersize=6)
plt.xlabel('channel')
plt.ylabel('Radiances' + '   (' + rad[v + '_ATTRIBUTE']['UNITS'] + ')')
plt.yscale('log')
plt.ylim(1, 150)
plt.legend()
plt.title(file)
plt.show()

for v in ['REFL', 'REFL_CLEAR']:
    color = colors[i]
    i += 1
    plt.plot(listchan, rad[v], color + '-', label=v, marker='o', markersize=6)
plt.xlabel('channel')
plt.ylabel('Reflectances' + '   (' + rad[v + '_ATTRIBUTE']['UNITS'] + ')')
plt.legend()
plt.title(file)
plt.show()
print ">>>>>>>>>>TEST OK"
