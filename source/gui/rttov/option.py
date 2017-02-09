#! /usr/bin/env python
# -*- coding: utf-8 -*-

import h5py
import re

from core import _V, loadDset_arr, loadDset_log, getAttributes, saveAttributes, checkAttributes


class Option(dict, _V):
    """
    The Option class
    A Python dictionary which contains a copy of the RTTOV v11
    rttov_options structure
    The creation method fills the class with the same defaults
    as the Fortran structure
    Inherits the _V class for some default methods
    """

    options_list_logical = ['ADDAEROSL', 'ADDCLOUDS', 'CLDSTR_SIMPLE',
                            'ADDINTERP',
                            'ADDPC', 'ADDRADREC', 'ADDREFRAC',
                            'ADDSOLAR', 'APPLY_REG_LIMITS', 'CH4_DATA',
                            'CLW_DATA',
                            'CO2_DATA', 'CO_DATA', 'DO_CHECKINPUT',
                            'DO_LAMBERTIAN',
                            'LGRADP', 'N2O_DATA', 'OZONE_DATA',
                            'SPACETOP', 'SWITCHRAD', 'USE_Q2M', 'VERBOSE',
                            'USER_AER_OPT_PARAM', 'USER_CLD_OPT_PARAM',
                            'USE_Q2M',
                            'DO_NLTE_CORRECTION', 'REG_LIMIT_EXTRAP',
                            'SUPPLY_FOAM_FRACTION']
    options_list = options_list_logical + ['CLDSTR_THRESHOLD',
                                           'FASTEM_VERSION',
                                           'IPCREG', 'IPCBND',
                                           'INTERP_MODE']

    def __init__(self):
        dict.__init__(self)
        self.default()

    def default(self):
        for key in self.options_list_logical:
            self[key] = False
        self['SPACETOP'] = True
        self['VERBOSE'] = True
        self['USE_Q2M'] = True
        self['DO_CHECKINPUT'] = True
        self['CLDSTR_THRESHOLD'] = -1.0
        self['IPCREG'] = -1
        self['IPCBND'] = -1
        self['FASTEM_VERSION'] = 5
        self['OZONE_DATA'] = True
        self['ADDINTERP'] = True
        self['INTERP_MODE'] = 5
        self['REG_LIMIT_EXTRAP'] = True
        self['APPLY_REG_LIMITS'] = True
        self['CLDSTR_SIMPLE'] = False  # hidden option so not put in the list


# ADDAEROSL_ATTRIBUTE      =>    {u'COMMENT': 'Account for scattering due toaerosols', u'UNITS': 'n/a'} <type 'dict'>
# ADDCLOUDS_ATTRIBUTE      =>    {u'COMMENT': 'account for scattering due to clouds', u'UNITS': 'n/a'} <type 'dict'>
# ADDINTERP_ATTRIBUTE      =>    {u'COMMENT': 'Input profiles may be supplied on user-defined levels', u'UNITS': 'n/a'} <type 'dict'>
# ADDPC_ATTRIBUTE      =>    {u'COMMENT': 'Carry out Principal Components calculations', u'UNITS': 'n/a'} <type 'dict'>
# ADDRADREC_ATTRIBUTE      =>    {u'COMMENT': 'If addpc is true, the PC calculations will return reconstructed radiances', u'UNITS': 'n/a'} <type 'dict'>
# ADDREFRAC_ATTRIBUTE      =>    {u'COMMENT': 'RTTOV calculations account for atmospheric refraction', u'UNITS': 'n/a'} <type 'dict'>
# ADDSOLAR_ATTRIBUTE      =>    {u'COMMENT': 'Enable solar calculations', u'UNITS': 'n/a'} <type 'dict'>
# APPLY_REG_LIMITS_ATTRIBUTE      =>    {u'COMMENT': 'Input profiles outside the limits specified in the coefficient files are reset to the min/max', u'UNITS': 'n/a'} <type 'dict'>
# CH4_DATA_ATTRIBUTE      =>    {u'COMMENT': 'Profile contains CH4 data', u'UNITS': 'n/a'} <type 'dict'>
# CLDSTR_THRESHOLD_ATTRIBUTE      =>    {u'COMMENT': 'Cloud threshold, recommended to set this negative for TL/AD/K', u'UNITS': 'n/a'} <type 'dict'>
# CLW_DATA_ATTRIBUTE      =>    {u'COMMENT': 'Profile contains Cloud Liquid Water (micro wave) data', u'UNITS': 'n/a'} <type 'dict'>
# CO2_DATA_ATTRIBUTE      =>    {u'COMMENT': 'Profile contains CO2 data', u'UNITS': 'n/a'} <type 'dict'>
# CO_DATA_ATTRIBUTE      =>    {u'COMMENT': 'Profile contains CO data', u'UNITS': 'n/a'} <type 'dict'>
# DO_CHECKINPUT_ATTRIBUTE      =>    {u'COMMENT': 'Checks whether input profiles are within both absolute and regression limits', u'UNITS': 'n/a'} <type 'dict'>
# DO_LAMBERTIAN_ATTRIBUTE      =>    {u'COMMENT': 'Activate treatment of surface as Lambertian instead of specular reflector for downwelling emitted radiation', u'UNITS': 'n/a'} <type 'dict'>
# DO_NLTE_CORRECTION_ATTRIBUTE      =>    {u'COMMENT': 'Enable Non Local Thermal Equilibrium (NLTE) calculations', u'UNITS': 'n/a'} <type 'dict'>
# FASTEM_VERSION_ATTRIBUTE      =>    {u'COMMENT': 'Fastem version; Valid range: 1-5. Otherwise version taken from coef file', u'UNITS': 'n/a'} <type 'dict'>
# INTERP_MODE_ATTRIBUTE      =>    {u'COMMENT': 'Interpolation mode (see user guide for options)', u'UNITS': 'n/a'} <type 'dict'>
# IPCBND_ATTRIBUTE      =>    {u'COMMENT': 'The index of the required spectral band for PC-RTTOV 1-3 for IASI clear-sky, 1 otherwise', u'UNITS': 'n/a'} <type 'dict'>
# IPCREG_ATTRIBUTE      =>    {u'COMMENT': 'The index of the required set of PC predictors 1-4 for IASI, 1-3 for AIRS', u'UNITS': 'n/a'} <type 'dict'>
# LGRADP_ATTRIBUTE      =>    {u'COMMENT': 'Allow TL/AD of user pressure levels if addinterp is true', u'UNITS': 'n/a'} <type 'dict'>
# N2O_DATA_ATTRIBUTE      =>    {u'COMMENT': 'Profile contains N2O data', u'UNITS': 'n/a'} <type 'dict'>
# OZONE_DATA_ATTRIBUTE      =>    {u'COMMENT': 'Profile contains Ozone data', u'UNITS': 'n/a'} <type 'dict'>
# REG_LIMIT_EXTRAP_ATTRIBUTE      =>    {u'COMMENT': 'extrapolate top profile levels using regression limits where input top level is below coef top level', u'UNITS': 'n/a'} <type 'dict'>
# SPACETOP_ATTRIBUTE      =>    {u'COMMENT': 'Treat user s model-top as space boundary', u'UNITS': 'n/a'} <type 'dict'>
# SWITCHRAD_ATTRIBUTE      =>    {u'COMMENT': 'Determines input perturbation for AD/K: if true, radiance_k%bt is used, otherwise radiance_k%total', u'UNITS': 'n/a'} <type 'dict'>
# USER_AER_OPT_PARAM_ATTRIBUTE      =>    {u'COMMENT': 'User input for aerosols optical parameters', u'UNITS': 'n/a'} <type 'dict'>
# USER_CLD_OPT_PARAM_ATTRIBUTE      =>    {u'COMMENT': 'User input for cloud optical parameters', u'UNITS': 'n/a'} <type 'dict'>
# USE_Q2M_ATTRIBUTE      =>    {u'COMMENT': 'Activate use of surface humidity', u'UNITS': 'n/a'} <type 'dict'>
# VERBOSE_ATTRIBUTE      =>    {u'COMMENT': 'verbose', u'UNITS': 'n/a'}
# <type 'dict'>
        for key in self.options_list:
            self[key + '_ATTRIBUTE'] = {'COMMENT': 'none', 'UNITS': 'n/a'}

    def loadh5(self, h5):
        # Load is driven by the content of the file

        for cle in list(h5):
            if(cle in self.options_list_logical):
                # Logicals are specific
                self[cle] = loadDset_log(h5[cle])
            else:
                # Reals, Integers
                self[cle] = loadDset_arr(h5[cle])
            self[cle + '_ATTRIBUTE'] = getAttributes(h5[cle])
            checkAttributes(self[cle + '_ATTRIBUTE'])

    def saveh5(self, h5, path):
        """Save option class to an HDF5 file h5 under path directory
        Attributes are saved at the same time datasets are saved
        """
        for cle in self.keys():

            if(re.search('_ATTRIBUTE', cle) or self[cle] is None):
                continue

            if(cle in self.options_list_logical):
                # Logicals are specific
                if(self[cle]):
                    val = 1
                else:
                    val = 0
                dset = h5.create_dataset(path + cle, data=val)
            else:
                dset = h5.create_dataset(path + cle, data=self[cle])

            saveAttributes(h5[path + cle], self[cle + '_ATTRIBUTE'])


if __name__ == '__main__':

    # Open file ReadOnly
    fileName = '../rttov_tests/options.H5'
    f = h5py.File(fileName, 'r')

    # get the Dataset
    h5 = f['/OPTIONS/']

    # p1 is a rttov_hdf_mod Option instance
    o1 = Option()

    # defaut
    o1.default()
    import sys
    sys.exit(0)

    # Load option
    o1.loadh5(h5)

    # Close HDF file
    f.close()

    # Display option
    o1.display()

    # Save option o1 to HDF5
    ofile = '../rttov_tests/options_res.H5'
    of = h5py.File(ofile, 'w')

    o1.saveh5(of, '/OPTIONS/')

    of.close()
