#!/usr/bin/env python

# ============================================================================
# README
# ============================================================================

# This script can generate all the test profile datasets for the RTTOV test
# suite using the profiles defined in the associated data Python script.
# (This doesn't include the diverse profile datasets).

# The option --gui allows to generate the Python files for the RTTOV_GUI.

# The purpose of scripting this is to make it much easier to modify the test
# profiles, for example to account for missed cases or new code features.

# The code at the bottom of this script demonstrates how each type of profile
# dataset is created. The code is quite general so should accommodate changes
# in array dimensions such as the addition of new profiles, changes in
# profile levels, etc, as long as they are implemented consistently across
# the input data.

# Gas units in the profile data are assumed to be ppmv over dry air.
# Test suite profiles are generally converted to ppmv over moist air.

# The associated data script should define the following variables:
# NPROF - the number of profiles to create per dataset
# NLEV  - the number of levels in the underlying profile set (e.g. US76)
#
# p, t, q, o3, co2, co,  - arrays of profiles in underlying set (e.g. US76)
# n2o, ch4
#
# clw                    - cloud liquid water array on same levels as
#                          basis profile set
#
# cld101, cfrac101,      - additional 101L arrays for cloud/aerosol profiles
# icede101, aer101         (can also define profiles on different numbers of
#                          levels if required - these profiles are not tied
#                          to particular pressure levels and are not
#                          currently interpolated by the code below)
#
# s2m, skin, datetime,   - remaining profile variables
# satzen, satazi,
# sunzen, sunazi,
# lat, elev, ctp,
# cfraction,
# be, cosbk, ish, idg
#
# p101, p54, p51         - arrays containing RTTOV standard pressures
#
# ch4ref101, co2ref101,  - trace gas reference profiles on 101L
# coref101, n2oref101,
# o3ref101, qref101

# co2ref54, o3ref54,     - trace gas reference profiles on 54L
# qref54
#
# co2ref51, qref51       - trace gas reference profiles on 51L
#
# Trace gases and standard pressures may be supplied for other levels.
#
# ============================================================================

from prof_gen_data import *
import numpy as np
import matplotlib.pyplot as plt
import os
import sys
import pprint
import argparse

UNIT_PPMVDRY = 0
UNIT_KGKGWET = 1
UNIT_PPMVWET = 2

# These masses should match those in rttov_const.F90
DRYAIRMASS = np.float64(28.9644)
GASMASSLIST = {'q'   : np.float64(18.01528),
               'o3'  : np.float64(47.9982),
               'co2' : np.float64(44.0095),
               'co'  : np.float64(28.0101),
               'ch4' : np.float64(16.04246),
               'n2o' : np.float64(44.0128)}

# ============================================================================
# Helper functions to create profiles and write profile datasets
# ============================================================================

# ----------------------------------------------------------------------------
# Conversion functions for gas units
# ----------------------------------------------------------------------------

def ppmvdry2kgkgwet(d, gasname):
    """Convert d ppmv dry to kgkg wet for gas gasname"""
    r = d * np.float64(1.E-06) * GASMASSLIST[gasname] / DRYAIRMASS
    return r / (1. + r)

def ppmvdry2ppmvwet(d, gasname, h2o_dry):
    """Convert d ppmv dry to ppmv wet for gas gasname given water vapour h2o_dry in ppmv dry"""
    if gasname == 'q':
        return d / (1. + np.float64(1.E-06) * d)
    else:
        return d * (1. - np.float64(1.E-06) * h2o_dry / (1. + np.float64(1.E-06) * h2o_dry))

def kgkgwet2ppmvdry(d, gasname):
    """Convert d kgkg wet to ppmv dry for gas gasname"""
    return np.float64(1.E06) * (d / (1. - d)) * DRYAIRMASS / GASMASSLIST[gasname]

def kgkgwet2ppmvwet(d, gasname, h2o_wet):
    """Convert d kg/kg wet to ppmv wet for gas gasname given water vapour h2o_wet in kg/kg wet"""
    ppmvdry = kgkgwet2ppmvdry(d, gasname)
    h2oppmvdry = kgkgwet2ppmvdry(h2o_wet, 'q')
    return ppmvdry2ppmvwet(ppmvdry, gasname, h2oppmvdry)

#def kgkgdry2ppmvwet(d, gasname, h2o_dry):
    #"""Convert d kg/kg dry to ppmv wet for gas gasname given water vapour h2o_dry in kg/kg dry"""
    #ppmvdry = np.float64(1.E06) * d * DRYAIRMASS / GASMASSLIST[gasname]
    #h2oppmvdry = np.float64(1.E06) * h2o_dry * DRYAIRMASS / GASMASSLIST['q']
    #return ppmvdry2ppmvwet(ppmvdry, gasname, h2oppmvdry)

# Simple tests to ensure consistency between conversion functions
assert(kgkgwet2ppmvwet(.2, 'o3', .3) == ppmvdry2ppmvwet(kgkgwet2ppmvdry(.2, 'o3'), 'o3', kgkgwet2ppmvdry(.3, 'q')))
assert(kgkgwet2ppmvdry(ppmvdry2kgkgwet(1000., 'o3'), 'o3') == 1000.)


# ----------------------------------------------------------------------------
# These two functions generate 101L profiles using variations on the AIRS formula
# ----------------------------------------------------------------------------

def make_p(nlev, i_fix, p_fix, e):
    """Make new pressure profiles using the AIRS 101L formula"""
    coef = np.zeros([3,3])
    v = np.array(map(lambda x:x**(1./e), p_fix))

    for r in range(3):
        for c in range(3):
          coef[r,c] = i_fix[r]**(2-c)

    abc = np.linalg.solve(coef, v)
    f = lambda i: (abc[0] * i**2 + abc[1] * i + abc[2])**e
    newp = map(f, range(nlev))
    newp.reverse()
    return np.array(newp)

def make_101L():
    """Create NPROF 101L pressure profiles using slightly variations on the AIRS formula
       This is fairly arbitrary, but it makes a selection of nice smooth profiles which
       differ significantly from one another"""
    plist = []
    for (i, f) in enumerate(np.linspace(0.93, 1.03, NPROF)):
        newp = make_p(101, [0, 35 + i, 100], [f * 1060., f * 300., f * 0.005], (7. + i)/2.)
        if np.max(newp) > 1100.:
            print 'ERROR: maximum generated pressure is > 1100 hPa'
            sys.exit(1)
        plist.append(newp)
        #plt.plot(newp)
    #plt.show()
    return plist

# ----------------------------------------------------------------------------
# Function to create the directories for a new profile dataset
# ----------------------------------------------------------------------------

def make_dirs(d, np, gui=False):
    """Create the necessary directories under the given dirname"""
    if not os.path.exists(d):
        os.makedirs(d)
    for i in range(1, np + 1):
        if gui:
            fn = d + '/{:03d}.py'.format(i)
            with open(fn, 'w') as f:
                f.write('"""\n    Profile {}\n'.format(fn))
                f.write('        file automaticaly created by prof_gen.py script\n')
                f.write('""" \n\n')
                f.write('self["ID"] = "{}"'.format(fn))

        else:
            for sd in ['atm', 'ground']:
                profdir = (d + '/{:03d}/' + sd).format(i)
                if not os.path.exists(profdir):
                    os.makedirs(profdir)

# ----------------------------------------------------------------------------
# The following functions write out different kinds of profile data
# ----------------------------------------------------------------------------

def write_list(fn, data, gui=False, name='', gas_unit=UNIT_PPMVWET, h2o_dry=None):
    """Write a 1-D profile to the given filename"""
    dat = data
    if name.lower() in GASMASSLIST:
        if gas_unit == UNIT_KGKGWET:
            dat = ppmvdry2kgkgwet(data, name.lower())
        elif gas_unit == UNIT_PPMVWET:
            dat = ppmvdry2ppmvwet(data, name.lower(), h2o_dry)
    if gui:
        with open(fn, 'a') as f:
            f.write('\nself["{}"] = numpy.'.format(name))
            pprint.pprint(dat, stream=f)
    else:
        with open(fn, 'w') as f:
            for v in dat:
                f.write('{:12.6E}'.format(v) + '\n')


def write_refgas(dirname, dlist, flist, gui=False, gas_unit=UNIT_PPMVWET, h2o_dry=None):
    """Write out the ref gas profiles from the dlist to the files in flist"""
    if h2o_dry is None: h2o_dry = [None] * NPROF
    if gui:
        for i in range(NPROF):
            for (d, fn) in zip(dlist, flist):
                gn = fn.replace('.txt','').upper()
                profdir = dirname + '/{:03d}.py'.format(i + 1)
                write_list(profdir, d, gui, gn, gas_unit, h2o_dry[i])

    else:
        for i in range(NPROF):
            for (d, fn) in zip(dlist, flist):
                gn = fn.replace('.txt','').upper()
                profdir = dirname + '/{:03d}/atm/'.format(i + 1)
                write_list(profdir + fn, d, gui, gn, gas_unit, h2o_dry[i])


def write_common(dirname, gui=False, gas_unit=UNIT_PPMVWET):
    """Writes the common files for all profiles in directory dirname:
       atm/cloud0.txt, atm/aerosli.txt, ground/*.txt, angles.txt, be.txt, datetime.txt"""

    if gui:
        for i in range(NPROF):
            d = dirname + '/{:03d}.py'.format(i + 1)
            with open(d , 'a') as f:
                f.write('\nself["GAS_UNITS"] = {}'.format(gas_unit))

                f.write('\nself["CTP"] = {}'.format(ctp[i]))
                f.write('\nself["CFRACTION"] = {}'.format(cfraction[i]))

                f.write('\nself["IDG"] = {}'.format(idg[i]))
                f.write('\nself["ISH"] = {}'.format(ish[i]))

                f.write('\nself["S2M"]["T"] = {}'.format(s2m['t'][i]))
                if gas_unit == UNIT_KGKGWET:
                    f.write('\nself["S2M"]["Q"] = {}'.format(ppmvdry2kgkgwet(s2m['q'][i], 'q')))
                    f.write('\nself["S2M"]["O"] = {}'.format(ppmvdry2kgkgwet(s2m['o'][i], 'o3')))
                elif gas_unit == UNIT_PPMVWET:
                    f.write('\nself["S2M"]["Q"] = {}'.format(ppmvdry2ppmvwet(s2m['q'][i], 'q', s2m['q'][i])))
                    f.write('\nself["S2M"]["O"] = {}'.format(ppmvdry2ppmvwet(s2m['o'][i], 'o3', s2m['q'][i])))
                else:
                    f.write('\nself["S2M"]["Q"] = {}'.format(s2m['q'][i]))
                    f.write('\nself["S2M"]["O"] = {}'.format(s2m['o'][i]))
                f.write('\nself["S2M"]["P"] = {}'.format(s2m['p'][i]))
                f.write('\nself["S2M"]["U"] = {}'.format(s2m['u'][i]))
                f.write('\nself["S2M"]["V"] = {}'.format(s2m['v'][i]))
                f.write('\nself["S2M"]["WFETC"] = {}'.format(s2m['wfetc'][i]))

                f.write('\nself["SKIN"]["SURFTYPE"] = {}'.format(skin['surftype'][i]))
                f.write('\nself["SKIN"]["WATERTYPE"] = {}'.format(skin['watertype'][i]))
                f.write('\nself["SKIN"]["T"] = {}'.format(skin['t'][i]))
                f.write('\nself["SKIN"]["SALINITY"] = {}'.format(skin['salinity'][i]))
                f.write('\nself["SKIN"]["FOAM_FRACTION"] = {}'.format(skin['foam_fraction'][i]))

                f.write('\nself["ELEVATION"] = {}'.format(elev[i]))
                f.write('\nself["ZENANGLE"] = {}'.format(satzen[i]))
                f.write('\nself["AZANGLE"] = {}'.format(satazi[i]))
                f.write('\nself["SUNZENANGLE"] = {}'.format(sunzen[i]))
                f.write('\nself["SUNAZANGLE"] = {}'.format(sunazi[i]))
                f.write('\nself["LATITUDE"] = {}'.format(lat[i]))

                f.write('\nself["BE"] = {}'.format(be[i]))
                f.write('\nself["COSBK"] = {}'.format(cosbk[i]))

            write_list(d, np.array(skin['fastem'][i][:]), gui, 'SKIN"]["FASTEM')
            write_list(d, np.array(datetime[i][0:3]), gui, 'DATE')
            write_list(d, np.array(datetime[i][3:6]), gui, 'TIME')

    else:
        for i in range(NPROF):
            d = dirname + '/{:03d}/atm/'.format(i + 1)
            with open(d + 'cloud0.txt', 'w') as f:
                f.write('&cloud\n')
                f.write('  ctp       = ' + str(ctp[i]) + '\n')
                f.write('  cfraction = ' + str(cfraction[i]) + '\n')
                f.write('/\n')

            with open(d + 'aerosli.txt', 'w') as f:
                f.write('&aerosli\n')
                f.write('  idg = ' + str(idg[i]) + '\n')
                f.write('  ish = ' + str(ish[i]) + '\n')
                f.write('/\n')

            d = dirname + '/{:03d}/ground/'.format(i + 1)
            with open(d + 's2m.txt', 'w') as f:
                f.write('&s2m\n')
                f.write('  s0%t     = ' + str(s2m['t'][i]) + '\n')
                if gas_unit == UNIT_KGKGWET:
                    f.write('  s0%q     = ' + str(ppmvdry2kgkgwet(s2m['q'][i], 'q')) + '\n')
                    f.write('  s0%o     = ' + str(ppmvdry2kgkgwet(s2m['o'][i], 'o3')) + '\n')
                elif gas_unit == UNIT_PPMVWET:
                    f.write('  s0%q     = ' + str(ppmvdry2ppmvwet(s2m['q'][i], 'q', s2m['q'][i])) + '\n')
                    f.write('  s0%o     = ' + str(ppmvdry2ppmvwet(s2m['o'][i], 'o3', s2m['q'][i])) + '\n')
                else:
                    f.write('  s0%q     = ' + str(s2m['q'][i]) + '\n')
                    f.write('  s0%o     = ' + str(s2m['o'][i]) + '\n')
                f.write('  s0%p     = ' + str(s2m['p'][i]) + '\n')
                f.write('  s0%u     = ' + str(s2m['u'][i]) + '\n')
                f.write('  s0%v     = ' + str(s2m['v'][i]) + '\n')
                f.write('  s0%wfetc = ' + str(s2m['wfetc'][i]) + '\n')
                f.write('/\n')

            with open(d + 'skin.txt', 'w') as f:
                f.write('&skin\n')
                f.write('  k0%surftype        = ' + str(skin['surftype'][i]) + '\n')
                f.write('  k0%watertype       = ' + str(skin['watertype'][i]) + '\n')
                f.write('  k0%t               = ' + str(skin['t'][i]) + '\n')
                f.write('  k0%salinity        = ' + str(skin['salinity'][i]) + '\n')
                f.write('  k0%foam_fraction   = ' + str(skin['foam_fraction'][i]) + '\n')
                f.write('  k0%fastem(1)       = ' + str(skin['fastem'][i][0]) + '\n')
                f.write('  k0%fastem(2)       = ' + str(skin['fastem'][i][1]) + '\n')
                f.write('  k0%fastem(3)       = ' + str(skin['fastem'][i][2]) + '\n')
                f.write('  k0%fastem(4)       = ' + str(skin['fastem'][i][3]) + '\n')
                f.write('  k0%fastem(5)       = ' + str(skin['fastem'][i][4]) + '\n')
                f.write('/\n')

            with open(d + 'elevation.txt', 'w') as f:
                f.write('&elev\n')
                f.write('  elevation = ' + str(elev[i]) + '\n')
                f.write('/\n')

            d = dirname + '/{:03d}/'.format(i + 1)
            with open(d + 'angles.txt', 'w') as f:
                f.write('&angles\n')
                f.write('  zenangle     = ' + str(satzen[i]) + '\n')
                f.write('  azangle      = ' + str(satazi[i]) + '\n')
                f.write('  sunzenangle  = ' + str(sunzen[i]) + '\n')
                f.write('  sunazangle   = ' + str(sunazi[i]) + '\n')
                f.write('  latitude     = ' + str(lat[i]) + '\n')
                f.write('/\n')

            with open(d + 'be.txt', 'w') as f:
                f.write(str(be[i]) + '   ' + str(cosbk[i]) + '\n')

            with open(d + 'datetime.txt', 'w') as f:
                f.write('   '.join(map(str, datetime[i])) + '\n')

            with open(d + 'gas_units.txt', 'w') as f:
                f.write('&units\n')
                f.write('  gas_units = ' + str(gas_unit) + '\n')
                f.write('/\n')


def write_seaonly_skin(dirname, gui=False):
    """Overwrite skin.txt for all-sea profile sets"""
    if gui:
        for i in range(NPROF):
            d = dirname + '/{:03d}.py'.format(i + 1)
            with open(d , 'a') as f:
                f.write('\nself["SKIN"]["SURFTYPE"] = {}'.format(skin['surftype_allsea'][i]))
                f.write('\nself["SKIN"]["WATERTYPE"] = {}'.format(skin['watertype'][i]))
                f.write('\nself["SKIN"]["T"] = {}'.format(skin['t'][i]))
                f.write('\nself["SKIN"]["SALINITY"] = {}'.format(skin['salinity'][i]))
                f.write('\nself["SKIN"]["FOAM_FRACTION"] = {}'.format(skin['foam_fraction'][i]))
            write_list(d, np.array(skin['fastem'][i][:]), gui, 'SKIN"]["FASTEM')

    else:
        for i in range(NPROF):
            d = dirname + '/{:03d}/ground/'.format(i + 1)
            with open(d + 'skin.txt', 'w') as f:
                f.write('&skin\n')
                f.write('  k0%surftype        = ' + str(skin['surftype_allsea'][i]) + '\n')
                f.write('  k0%watertype       = ' + str(skin['watertype'][i]) + '\n')
                f.write('  k0%t               = ' + str(skin['t'][i]) + '\n')
                f.write('  k0%salinity        = ' + str(skin['salinity'][i]) + '\n')
                f.write('  k0%foam_fraction   = ' + str(skin['foam_fraction'][i]) + '\n')
                f.write('  k0%fastem(1)       = ' + str(skin['fastem'][i][0]) + '\n')
                f.write('  k0%fastem(2)       = ' + str(skin['fastem'][i][1]) + '\n')
                f.write('  k0%fastem(3)       = ' + str(skin['fastem'][i][2]) + '\n')
                f.write('  k0%fastem(4)       = ' + str(skin['fastem'][i][3]) + '\n')
                f.write('  k0%fastem(5)       = ' + str(skin['fastem'][i][4]) + '\n')
                f.write('/\n')


def write_cloud(dirname, cld, cfrac, icede, nlev, gui=False):
    """Write out cloud.txt, cfrac.txt and icede.txt"""
    ntyp = cld.shape[1] / (nlev - 1)
    if cld.shape[1] % ntyp > 0:
        print 'ERROR: number of cloud levels does not match supplied nlev, ', nlev
        sys.exit(1)

    cloud_list = ['STCO', 'STMA', 'CUCC', 'CUCP', 'CUMA', 'CIRR']

    form = lambda x:'{:6.3f}'.format(x)

    if gui:
        for i in range(NPROF):
            d = dirname + '/{:03d}.py'.format(i + 1)
            acld = cld.reshape((NPROF, nlev-1, ntyp))
            for ityp in range(ntyp):
                if( any(acld[i,:,ityp] > 0) ):
                    write_list(d, acld[i,:,ityp], gui, cloud_list[ityp])
            write_list(d, cfrac[i,:], gui, "CFRAC")

            if np.any(icede[i,:] > 0):
                write_list(d, icede[i,:], gui, "ICEDE")

    else:
        for i in range(NPROF):
            d = dirname + '/{:03d}/atm/'.format(i + 1)
            with open(d + 'cloud.txt', 'w') as f:
                for l in range(nlev - 1):
                    f.write(' '.join(map(form, cld[i,ntyp*l:ntyp*(l+1)])) + '\n')

            with open(d + 'cfrac.txt', 'w') as f:
                for l in range(nlev - 1):
                    f.write(form(cfrac[i,l]) + '\n')

            if np.any(icede[i,:] < 0): continue
            with open(d + 'icede.txt', 'w') as f:
                for l in range(nlev - 1):
                    f.write(form(icede[i,l]) + '\n')


def write_aerosol(dirname, aer, nlev, gui=False):
    """Write out aerosl.txt"""
    ntyp = aer.shape[1] / (nlev - 1)
    if aer.shape[1] % ntyp > 0:
        print 'ERROR: number of aerosol levels does not match supplied nlev, ', nlev
        sys.exit(1)

    aerosol_list = ['INSO','WASO','SOOT','SSAM','SSCM','MINM','MIAM','MICM','MITR',  \
                    'SUSO','VOLA','VAPO','ASDU']

    form = lambda x:'{:9.4f}'.format(x)

    if gui:
        for i in range(NPROF):
            d = dirname + '/{:03d}.py'.format(i + 1)
            aaer = aer.reshape((NPROF, nlev-1, ntyp))
            for ityp in range(ntyp):
                if( any(aaer[i,:,ityp] > 0) ):
                    write_list(d, aaer[i,:,ityp], gui, aerosol_list[ityp] )

    else:
        for i in range(NPROF):
            d = dirname + '/{:03d}/atm/'.format(i + 1)
            with open(d + 'aerosl.txt', 'w')  as f:
                for l in range(nlev - 1):
                    f.write(' '.join(map(form, aer[i,ntyp*l:ntyp*(l+1)])) + '\n')


# ----------------------------------------------------------------------------
# The following two functions are the main ones to call for generating profile
# datasets from the underlying set (e.g. US76)
# ----------------------------------------------------------------------------

def make_interp_profs(dirname, gui, pout, g, gas_unit=UNIT_PPMVWET):
    """Write out the NPROF underlying profiles interpolated onto the p levels in each member of pout"""

    make_dirs(dirname, NPROF, gui)

    h2o_dry = [None] * NPROF
    for i in range(NPROF):
        if gui:
            proffn = dirname + '/{:03d}.py'.format(i + 1)
            write_list(proffn, pout[i], True, 'P')
        else:
            profdir = dirname + '/{:03d}/atm/'.format(i + 1)
            write_list(profdir + 'p.txt', pout[i])

        logpin = np.log(p[i,:])
        logpout = np.log(pout[i])

        for (din, fn) in zip(DATALIST[g], FILELIST[g]):
            dout = np.interp(x=logpout, xp=logpin[::-1], fp=din[i,::-1], left=-999., right=-999.)
            for j in range(len(dout)):
                if dout[j] >= 0.:
                    continue
                dout[j] = dout[j-1] + (logpout[j] - logpout[j-1]) * (dout[j-1] - dout[j-2]) / (logpout[j-1] - logpout[j-2])
            if fn == 'q.txt': h2o_dry[i] = dout
            if gui:
                gn = fn.replace('.txt','').upper()
                write_list(proffn, dout, gui, gn, gas_unit, h2o_dry[i])
            else:
                gn = fn.replace('.txt','').upper()
                write_list(profdir + fn, dout, gui, gn, gas_unit, h2o_dry[i])

    write_common(dirname, gui, gas_unit)

    return h2o_dry # For ref_gas generation

def make_nointerp_profs(dirname, gui, nlevels, g, gas_unit=UNIT_PPMVWET):
    """Write out the bottom nlevels of the NPROF underlying profiles without any interpolation (1 <= nlevels <= NLEV)"""

    make_dirs(dirname, NPROF, gui)

    nlevels = max(1, min(nlevels, NLEV))
    for i in range(NPROF):
        if gui:
            proffn = dirname + '/{:03d}.py'.format(i + 1)
            write_list(proffn, p[i,nlevels-1::-1], True, 'P')
        else:
            profdir = dirname + '/{:03d}/atm/'.format(i + 1)
            write_list(profdir + 'p.txt', p[i,nlevels-1::-1])

        for (din, fn) in zip(DATALIST[g], FILELIST[g]):
            if gui:
                gn = fn.replace('.txt','').upper()
                write_list(proffn, din[i,nlevels-1::-1], gui, gn, gas_unit, q[i,nlevels-1::-1])
            else:
                gn = fn.replace('.txt','').upper()
                write_list(profdir + fn, din[i,nlevels-1::-1], gui, gn, gas_unit, q[i,nlevels-1::-1])

    write_common(dirname, gui, gas_unit)


# ============================================================================
# Define gas combinations for output profile datasets
# ============================================================================

# Each "DATALIST" contains lists of profile data in the data module for each gas combination
# Each "FILELIST" contains corresponding names of output files

DATALIST = {'allgas' : [t, q, o3, co2, co, n2o, ch4],
            'nogas'  : [t, q],
            'co2'    : [t, q, co2],
            'o3'     : [t, q, o3],
            'co2o3'  : [t, q, o3, co2],
            'clw'    : [t, q, clw]}
FILELIST = {'allgas' : ['t.txt', 'q.txt', 'o3.txt', 'co2.txt', 'co.txt', 'n2o.txt', 'ch4.txt'],
            'nogas'  : ['t.txt', 'q.txt'],
            'co2'    : ['t.txt', 'q.txt', 'co2.txt'],
            'o3'     : ['t.txt', 'q.txt', 'o3.txt'],
            'co2o3'  : ['t.txt', 'q.txt', 'o3.txt', 'co2.txt'],
            'clw'    : ['t.txt', 'q.txt', 'clw.txt']}
GASLIST = DATALIST.keys()

REFGAS101DATALIST = {'allgasref' : [o3ref101, co2ref101, coref101, n2oref101, ch4ref101]}
REFGAS101FILELIST = {'allgasref' : ['o3.txt', 'co2.txt', 'co.txt', 'n2o.txt', 'ch4.txt']}
REFGAS101LIST = REFGAS101DATALIST.keys()

REFGAS54DATALIST = {'co2ref'    : [co2ref54],
                    'o3ref'     : [o3ref54],
                    'co2o3ref'  : [o3ref54, co2ref54]}
REFGAS54FILELIST = {'co2ref'    : ['co2.txt'],
                    'o3ref'     : ['o3.txt'],
                    'co2o3ref'  : ['o3.txt', 'co2.txt']}
REFGAS54LIST = REFGAS54DATALIST.keys()

REFGAS51DATALIST = {'co2ref'    : [co2ref51]}
REFGAS51FILELIST = {'co2ref'    : ['co2.txt']}
REFGAS51LIST = REFGAS51DATALIST.keys()


# ============================================================================
# Make the profile sets
# ============================================================================

# ----------------------------------------------------------------------------
# Standard and varying
# ----------------------------------------------------------------------------

def make_std_varying(gui, gas_unit=UNIT_PPMVWET):
    for g in GASLIST:
        # standard 54L
        dirname = 'standard54lev_' + g
        make_interp_profs(dirname, gui, [p54] * NPROF, g, gas_unit)

        # standard 101L
        dirname = 'standard101lev_' + g
        make_interp_profs(dirname, gui, [p101] * NPROF, g, gas_unit)

        # varying 101L
        dirname = 'varying101lev_' + g
        plist = make_101L()
        make_interp_profs(dirname, gui, plist, g, gas_unit)

    for g in ['co2']:
        # standard 51L
        dirname = 'standard51lev_' + g
        make_interp_profs(dirname, gui, [p51] * NPROF, g, gas_unit)

    # We need some additional 'nogas' profiles for the refgas_vs_nogas test
    g = 'nogas'
    dirname = 'standard54lev_' + g + '_ppmvdry'
    make_interp_profs(dirname, gui, [p54] * NPROF, g, UNIT_PPMVDRY)

    dirname = 'standard101lev_' + g + '_ppmvdry'
    make_interp_profs(dirname, gui, [p101] * NPROF, g, UNIT_PPMVDRY)


# ----------------------------------------------------------------------------
# Sea-surface only
# ----------------------------------------------------------------------------

# Make the profiles as above and then overwrite skin.txt

def make_sea_only(gui, gas_unit=UNIT_PPMVWET):
    for g in ['allgas']:
        # standard 54L
        dirname = 'standard54lev_' + g + '_seaonly'
        make_interp_profs(dirname, gui, [p54] * NPROF, g, gas_unit)
        write_seaonly_skin(dirname, gui)

        # standard 101L
        dirname = 'standard101lev_' + g + '_seaonly'
        make_interp_profs(dirname, gui, [p101] * NPROF, g, gas_unit)
        write_seaonly_skin(dirname, gui)

        # varying 101L
        dirname = 'varying101lev_' + g + '_seaonly'
        plist = make_101L()
        make_interp_profs(dirname, gui, plist, g, gas_unit)
        write_seaonly_skin(dirname, gui)

    for g in ['o3']:
        # varying 101L
        dirname = 'varying101lev_' + g + '_seaonly'
        plist = make_101L()
        make_interp_profs(dirname, gui, plist, g, gas_unit)
        write_seaonly_skin(dirname, gui)

# ----------------------------------------------------------------------------
# Cloud and aerosol
# ----------------------------------------------------------------------------

# Note that the cloud/aerosol profiles are not interpolated
# and are not associated with particular pressure levels:
# the script creates varying101L profiles as above and then
# "drops" the cld/aer profiles in.

def make_cld_aer(gui, gas_unit=UNIT_PPMVWET):
    for g in ['allgas', 'o3']:
        # cloud 101L
        dirname = 'cld101lev_' + g
        plist = make_101L()
        make_interp_profs(dirname, gui, plist, g, gas_unit)
        write_cloud(dirname, cld101, cfrac101, icede101, 101, gui)

        # aerosol 101L
        dirname = 'aer101lev_' + g
        plist = make_101L()
        make_interp_profs(dirname, gui, plist, g, gas_unit)
        write_aerosol(dirname, aer101, 101, gui)

        # cloud+aerosol 101L
        dirname = 'cldaer101lev_' + g
        plist = make_101L()
        make_interp_profs(dirname, gui, plist, g, gas_unit)
        write_cloud(dirname, cld101, cfrac101, icede101, 101, gui)
        write_aerosol(dirname, aer101, 101, gui)

    # Also make a sea-only allgas profile for cloudy PC-RTTOV
    # NB profiles 5 and 6 have ish=3 which shouldn't be used with PC-RTTOV
    for g in ['allgas']:
        # cloud 101L
        dirname = 'cld101lev_' + g + '_seaonly'
        plist = make_101L()
        make_interp_profs(dirname, gui, plist, g, gas_unit)
        write_cloud(dirname, cld101, cfrac101, icede101, 101, gui)
        write_seaonly_skin(dirname, gui)

# ----------------------------------------------------------------------------
# Truncated (low model top)
# ----------------------------------------------------------------------------

# Exactly the same as above, but the pressure profile is truncated

def make_trunc(gui, gas_unit=UNIT_PPMVWET):
    for g in ['allgas', 'o3']:
        # standard 54L truncated (low model top)
        dirname = 'standard54lev_' + g + '_trunc'
        make_interp_profs(dirname, gui, [p54[19:]] * NPROF, g, gas_unit)

# ----------------------------------------------------------------------------
# With reference gas profiles
# ----------------------------------------------------------------------------

# Create the profile on standard levels as above and then
# overwrite the gas profiles with the references

def make_refgas(gui, gas_unit=UNIT_PPMVDRY):
    for g in REFGAS101LIST:
        dirname = 'standard101lev_' + g
        h2o_dry = make_interp_profs(dirname, gui, [p101] * NPROF, g[:-3], gas_unit)
        write_refgas(dirname, REFGAS101DATALIST[g], REFGAS101FILELIST[g], gui, gas_unit, h2o_dry)

    for g in REFGAS54LIST:
        dirname = 'standard54lev_' + g
        h2o_dry = make_interp_profs(dirname, gui, [p54] * NPROF, g[:-3], gas_unit)
        write_refgas(dirname, REFGAS54DATALIST[g], REFGAS54FILELIST[g], gui, gas_unit, h2o_dry)

    #for g in REFGAS51LIST:
        #dirname = 'standard51lev_' + g
        #h2o_dry = make_interp_profs(dirname, gui, [p51] * NPROF, g[:-3], gas_unit)
        #write_refgas(dirname, REFGAS51DATALIST[g], REFGAS51FILELIST[g], gui, gas_unit, h2o_dry)

# ----------------------------------------------------------------------------
# Uninterpolated basis profiles
# ----------------------------------------------------------------------------

def make_noninterp(gui, gas_unit=UNIT_PPMVWET):
    for g in ['allgas']:
        # US76 on 43L (up to ~0.005hPa)
        dirname = 'us76_43lev_' + g
        make_nointerp_profs(dirname, gui, 43, g, gas_unit)

        # US76 on 50L
        dirname = 'us76_50lev_' + g
        make_nointerp_profs(dirname, gui, NLEV, g, gas_unit)

# ----------------------------------------------------------------------------
# Call required functions
# ----------------------------------------------------------------------------

def parse_args():
    parser = argparse.ArgumentParser(description='Generate RTTOV test profiles', conflict_handler='resolve')
    parser.add_argument('-g', '--gui', dest='gui', action='store_true', help='generates Python file for RTTOV_GUI')
    return parser.parse_args()

args = parse_args()
gui = args.gui

make_std_varying(gui)
make_sea_only(gui)
make_cld_aer(gui)
make_trunc(gui)
make_refgas(gui)
make_noninterp(gui)

