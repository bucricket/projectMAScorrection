#! /usr/bin/env python
# -*- coding: utf-8 -*-
try:
    import h5py
except ImportError:
    import sys
    sys.stderr.write('ERROR: h5py is not installed\n')
    sys.exit(1)
try:
    import numpy
except ImportError:
    import sys
    sys.stderr.write('ERROR: numpy is not installed\n')
    sys.exit(1)
import re

from core import loadDset_arr, loadDset_str, \
    getAttributes, saveAttributes, checkAttributes, cleanpath
from rttov_gui_f2py import rttov_gui_aer_clim_prof

# conversion constantes for gas content :
q_mixration_to_ppmv = 1.60771704e+6
o3_mixratio_to_ppmv = 6.03504e+5
co2_mixratio_to_ppmv = 6.58114e+5
co_mixratio_to_ppmv = 1.0340699e+6
n2o_mixratio_to_ppmv = 6.58090e+5
ch4_mixratio_to_ppmv = 1.80548e+6
q_mixratio_to_ppmv = 1.60771704e+6
mixratio_to_ppmv = {"Q": q_mixration_to_ppmv,
                    "O3": o3_mixratio_to_ppmv,
                    "CO2": co2_mixratio_to_ppmv,
                    "CO": co_mixratio_to_ppmv,
                    "NO2": n2o_mixratio_to_ppmv,
                    "CH4": ch4_mixratio_to_ppmv,
                    "Q": q_mixratio_to_ppmv,
                    }


def fct_mixratio_to_ppmv(value, gas):
    return value * mixratio_to_ppmv[gas]


def fct_ppmv_to_mixratio(value, gas):
    return value / mixratio_to_ppmv[gas]


def getNumberOfProfiles(fileName):
    """
    Get the number of profiles in and HDF5 file
    If the /PROFILE is found number of profiles is 1
    If the /PROFILES is found then the number of profiles
      is the number of subdirectories which names are formed of 4 digits
    """
    f = h5py.File(fileName, 'r')
    nprof = 0
    keys = f.keys()
    if("PROFILE" in keys):
        nprof = 1
    elif("PROFILES" in keys):
        for name, value in f['PROFILES'].iteritems():
            if (re.search('[0-9]{4}', name)):
                nprof += 1
    f.close()
    return nprof


class Profile(dict):
    """
    The Profile class
    A Python dictionary which contains a copy of the RTTOV v11
    profile_type structure
    All members are present, for example if the profile does not
    contain any CO2
    the CO2 key is present and is set to the Python None
    Methods are:
    - loadh5 which returns a full profile from a given dataset
    - loadprofilenumber which returns a profile defined by its number
      from a HDF5 file containing several profiles
    - display display profile to terminal
    - saveh5 saves a profile class to a HDF5 file
    - addgas will add a new gas to the profile filled witha  default
      concentration or a user value
    - removegas remove a gas from the profile
    - removeallbutgas remove all gases except the given one
    """
    # Default gas concentrations values for 2012
    defaultppmv = {'Q':    2e4, 'O3': 0.5,   'CO2': 392.0,
                   'CH4': 1809.0, 'CO': 0.174, 'N2O': 325.0}
    defaultaerdensity = 1.e0
    defaultclddensity = 1.e-01

    # All known profile_type members
    gas_list = ['Q', 'O3', 'CO2', 'CH4', 'CO', 'N2O']
    profile_list = ['AEROSOLS', 'AZANGLE', 'BE', 'CFRAC', 'CFRACTION', 'CH4',
                    'CLOUD',
                    'CLW', 'CO', 'CO2', 'COSBK', 'CTP',
                    'DATE', 'ELEVATION', 'GAS_UNITS', 'ICEDE', 'ID', 'IDG',
                    'ISH', 'LATITUDE', 'LONGITUDE', 'N2O',
                    'NLAYERS', 'NLEVELS', 'O3', 'P', 'Q',
                    'SNOW_FRAC', 'SOIL_MOISTURE',
                    'SUNAZANGLE', 'SUNZENANGLE', 'T', 'TIME', 'ZENANGLE']
    sskin_list = ['SURFTYPE', 'WATERTYPE', 'T',
                  'SALINITY', 'FOAM_FRACTION', 'FASTEM']
    s2m_list = ['T', 'Q', 'O', 'P', 'U', 'V', 'WFETC']
    aerosol_list = ['INSO', 'WASO', 'SOOT', 'SSAM', 'SSCM', 'MINM', 'MIAM',
                    'MICM', 'MITR',
                    'SUSO', 'VOLA', 'VAPO', 'ASDU']

    aerosol_long_list = ['Insoluble', 'Water soluble', 'Soot',
                         'Sea salt (acc mode)', 'Sea salt (coa mode)',
                         'Mineral (nuc mode)', 'Mineral (acc mode)',
                         'Mineral (coa mode)',
                         'Mineral transported', 'Sulphated droplets',
                         'OPAC Volcanic ash',
                         'New Volcanic ash', 'Asian dust']

    cloud_list = ['STCO', 'STMA', 'CUCC', 'CUCP', 'CUMA', 'CIRR']
    cloud_long_list = ['Stratus Continental', 'Stratus Maritime',
                       'Cumulus Continental Clean',
                       'Cumulus Continental Polluted', 'Cumulus Maritime',
                       'Cirrus']
    gas_units = {-1: "ppmv over dry air", 0: "ppmv",
                 1: "kg/kg over moist air", 2: "ppmv over moist air"}

    # add by Pascale

    profileParamList = ['DATE', 'TIME', 'ID', 'LATITUDE', 'LONGITUDE',
                        'ELEVATION', 'AZANGLE', 'ZENANGLE', 'SUNAZANGLE',
                        'SUNZENANGLE', 'CFRACTION',
                        'BE', 'COSBK', 'CTP', 'IDG', 'ISH', 'SNOW_FRAC',
                        'SOIL_MOISTURE']

    maxValue = {'LATITUDE': 90., 'LONGITUDE': 180., 'ELEVATION': 9000.,
                'AZANGLE': 180., 'ZENANGLE': 90.0, 'SUNAZANGLE': 180.,
                'SUNZENANGLE': 90., 'CFRACTION': 1.,
                'BE': 0.7, 'COSBK': 1., 'CTP': 1100., 'IDG': 4, 'ISH': 3,
                'SNOW_FRAC': 1., 'SOIL_MOISTURE': 1.,
                'S2M_T': 400., 'S2M_Q': 0.60E+06, 'S2M_O': 1000.,
                'S2M_P': 1100.0, 'S2M_U': 100., 'S2M_V': 100.,
                "S2M_WFETC": 100000.,
                "SKIN_SALINITY": 100., "SKIN_FOAM_FRACTION": 1.,
                "SKIN_SURFTYPE": 2, 'SKIN_T': 400., 'SKIN_WATERTYPE': 1
                }
    minValue = {'LATITUDE': -90., 'LONGITUDE': -180.,
                'ELEVATION': -500., 'AZANGLE': -180., 'ZENANGLE': 0.0,
                'SUNAZANGLE': -180, 'SUNZENANGLE': 0.0, 'CFRACTION': 0.0,
                'BE': 0.201, 'COSBK': 0., 'CTP': 50., 'IDG': 1, 'ISH': 1,
                'SNOW_FRAC': 0., 'SOIL_MOISTURE': 0.,
                'S2M_T': 90., 'S2M_Q': 0.1E-10, 'S2M_O': 0.1E-10,
                'S2M_P': 400.0, 'S2M_U': -100., 'S2M_V': -100.,
                "S2M_WFETC": 0.,
                "SKIN_SALINITY": 0., "SKIN_FOAM_FRACTION": 0.,
                "SKIN_SURFTYPE": 0, 'SKIN_T': 90., 'SKIN_WATERTYPE': 0
                }

    typeListe = {"IDG": {1: "Ou and Liou", 2: "Wyser", 3: "Boudala",
                         4: "McFarquhar"},
                 "ISH": {1: "Hexagonal crystals", 2: "Ice aggregates",
                         3: "Baran 2013", 4: "Baran 2014"},
                 "SURFTYPE": {0: "land", 1: "sea", 2: "sea-ice"},
                 "WATERTYPE": {0: "fresh water", 1: "ocean water"}}

    def __init__(self):
        dict.__init__(self)
        self['S2M'] = {}
        self['SKIN'] = {}

        for key in self.profile_list:
            self[key] = None
            self[key + '_ATTRIBUTE'] = {'COMMENT': 'none', 'UNITS': 'n/a'}
        for key in self.s2m_list:
            self['S2M'][key] = None
            self['S2M'][
                key + '_ATTRIBUTE'] = {'COMMENT': 'none', 'UNITS': 'n/a'}
        for key in self.sskin_list:
            self['SKIN'][key] = None
            self['SKIN'][
                key + '_ATTRIBUTE'] = {'COMMENT': 'none', 'UNITS': 'n/a'}
        for key in self.aerosol_list:
            self[key] = None
            self[key + '_ATTRIBUTE'] = {'COMMENT': 'none', 'UNITS': 'n/a'}
        for key in self.cloud_list:
            self[key] = None
            self[key + '_ATTRIBUTE'] = {'COMMENT': 'none', 'UNITS': 'n/a'}

    def setDefaultAttributes(self):
        """
         Set the default values for ALL the profile variables even if
         variables are set to None
        """
        self["ASDU_ATTRIBUTE"] = {
            'COMMENT': 'Asian dust', 'UNITS': 'number density (cm-3)'}
        self["AZANGLE_ATTRIBUTE"] = {
            'COMMENT': 'Satellite azimuth angle',
            'UNITS': 'degree (0-360 deg; east=90)'}
        self["BE_ATTRIBUTE"] = {
            'COMMENT': 'Earth magnetic field strength', 'UNITS': 'Gauss'}
        self["CFRACTION_ATTRIBUTE"] = {
         'COMMENT': 'Black body cloud fraction (0 - 1) 1 for 100% cloud cover',
         'UNITS': '0:1'}
        self["CFRAC_ATTRIBUTE"] = {
            'COMMENT': 'Cloud fraction (0 - 1) 1 for 100% cloud cover',
            'UNITS': '0:1'}
        self["CH4_ATTRIBUTE"] = {'COMMENT': 'CH4 (ppmv)', 'UNITS': 'ppmv'}
        self["CIRR_ATTRIBUTE"] = {'COMMENT': 'Cirrus', 'UNITS': 'g/m-3'}
        self["CLW_ATTRIBUTE"] = {'COMMENT': 'none', 'UNITS': 'n/a'}
        self["CO2_ATTRIBUTE"] = {
            'COMMENT': 'Carbon dioxide (ppmv)', 'UNITS': 'ppmv'}
        self["COSBK_ATTRIBUTE"] = {
            'COMMENT': ('Cosine of the angle between the Earth '
                        'magnetic field and wave propagation direction'),
            'UNITS': 'n/a'}
        self["CO_ATTRIBUTE"] = {'COMMENT': 'CO (ppmv)', 'UNITS': 'ppmv'}
        self["CTP_ATTRIBUTE"] = {
            'COMMENT': 'Black body cloud top pressure  (hPa)', 'UNITS': 'hPa'}
        self["CUCC_ATTRIBUTE"] = {
            'COMMENT': 'Cumulus Continental Clean', 'UNITS': 'g/m-3'}
        self["CUCP_ATTRIBUTE"] = {
            'COMMENT': 'Cumulus Continental Polluted', 'UNITS': 'g/m-3'}
        self["CUMA_ATTRIBUTE"] = {
            'COMMENT': 'Cumulus Maritime', 'UNITS': 'g/m-3'}
        self["DATE_ATTRIBUTE"] = {'COMMENT': 'Year Month Day'}
        self["ELEVATION_ATTRIBUTE"] = {
            'COMMENT': 'Surface elevaton', 'UNITS': 'km'}
        self["GAS_UNITS_ATTRIBUTE"] = {
            'COMMENT': 'Units for gas profiles', 'UNITS': 'n/a'}
        self["ICEDE_ATTRIBUTE"] = {
            'COMMENT': 'Ice crystals diameter', 'UNITS': 'microns'}
        self["IDG_ATTRIBUTE"] = {
            'COMMENT': ('Scheme for IWC to eff diameter,Dg 1=Ou and Liou; 2='
                        'Wyser et al.; 3=Boudala et al; 4=McFarquhar et al.'),
            'UNITS': 'n/a'}
        self["ID_ATTRIBUTE"] = {
            'COMMENT': 'User may give text ID to each profile', 'UNITS': 'n/a'}
        self["INSO_ATTRIBUTE"] = {
            'COMMENT': 'Insoluble', 'UNITS': 'number density (cm-3)'}
        self["ISH_ATTRIBUTE"] = {
            'COMMENT': ('Choose the shape of ice crystals, 1=Hexagonal ice '
                        'crystals; 2=Ice aggregates; 3=Baran 2013;'
                        ' 4=Baran 2014'),
            'UNITS': 'n/a'}
        self["LATITUDE_ATTRIBUTE"] = {
            'COMMENT': 'Latitude (deg)', 'UNITS': 'degree'}
        self["LONGITUDE_ATTRIBUTE"] = {
            'COMMENT': 'Longitude (deg)', 'UNITS': 'degree (0-360)'}
        self["MIAM_ATTRIBUTE"] = {
            'COMMENT': 'Mineral (acc mode)', 'UNITS': 'number density (cm-3)'}
        self["MICM_ATTRIBUTE"] = {
            'COMMENT': 'Mineral (coa mode)', 'UNITS': 'number density (cm-3)'}
        self["MINM_ATTRIBUTE"] = {
            'COMMENT': 'Mineral (nuc mode)', 'UNITS': 'number density (cm-3)'}
        self["MITR_ATTRIBUTE"] = {
            'COMMENT': 'Mineral transported', 'UNITS': 'number density (cm-3)'}
        self["N2O_ATTRIBUTE"] = {'COMMENT': 'N2O (ppmv)', 'UNITS': 'ppmv'}
        self["NLAYERS_ATTRIBUTE"] = {
            'COMMENT': 'Number of atmospheric layers', 'UNITS': 'n/a'}
        self["NLEVELS_ATTRIBUTE"] = {
            'COMMENT': 'Number of atmospheric levels', 'UNITS': 'n/a'}
        self["O3_ATTRIBUTE"] = {'COMMENT': 'Ozone (ppmv)', 'UNITS': 'ppmv'}
        self["P_ATTRIBUTE"] = {'COMMENT': 'Pressure (hPa)', 'UNITS': 'hPa'}
        self["Q_ATTRIBUTE"] = {
            'COMMENT': 'Water vapour (ppmv)', 'UNITS': 'ppmv'}
        self["S2M"]["WFETC_ATTRIBUTE"] = {
            'COMMENT': 'Wind fetch (metres)', 'UNITS': 'm'}
        self["S2M"]["P_ATTRIBUTE"] = {
            'COMMENT': 'Surface pressure (hPa)', 'UNITS': 'hPa'}
        self["S2M"]["Q_ATTRIBUTE"] = {
            'COMMENT': 'Water vapour (ppmv)', 'UNITS': 'ppmv'}
        self["S2M"]["U_ATTRIBUTE"] = {
            'COMMENT': 'U 10m wind component (m/s)', 'UNITS': 'm/s'}
        self["S2M"]["V_ATTRIBUTE"] = {
            'COMMENT': 'V 10m wind component (m/s)', 'UNITS': 'm/s'}
        self["S2M"]["T_ATTRIBUTE"] = {
            'COMMENT': 'Temperature (K)', 'UNITS': 'K'}
        self["S2M"]["O_ATTRIBUTE"] = {
            'COMMENT': 'Ozone (ppmv)', 'UNITS': 'ppmv'}
        self["SKIN"]["FASTEM_ATTRIBUTE"] = {
            u'COMMENT': 'Land/sea-ice surface parameters for fastem',
            'UNITS': 'n/a'}
        self["SKIN"]["WATERTYPE_ATTRIBUTE"] = {
            u'COMMENT': '0=fresh water, 1=ocean water', 'UNITS': 'n/a'}
        self["SKIN"]["T_ATTRIBUTE"] = {
            u'COMMENT': 'Radiative skin temperature (K)', u'UNITS': 'K'}
        self["SKIN"]["SURFTYPE_ATTRIBUTE"] = {
            u'COMMENT': '0=land, 1=sea, 2=sea-ice', 'UNITS': 'n/a'}
        self["SKIN"]["SALINITY_ATTRIBUTE"] = {
            u'COMMENT': 'Practical salinity unit %o - FASTEM-4/5 only',
            'UNITS': 'n/a'}
        self["SKIN"]["FOAM_FRACTION_ATTRIBUTE"] = {
            u'COMMENT': 'Ocean foam coverage fraction passed to FASTEM',
            'UNITS': '0:1'}
        self["SNOW_FRAC_ATTRIBUTE"] = {
            'COMMENT': ('Surface snow coverage fraction (0-1).'
                        ' Used only by IR emissivity atlas'), 'UNITS': '0:1'}
        self["SOIL_MOISTURE_ATTRIBUTE"] = {
            'COMMENT': 'soil moisture', 'UNITS': 'm^3/m^3'}
        self["SOOT_ATTRIBUTE"] = {'COMMENT': 'Soot',
                                  'UNITS': 'number density (cm-3)'}
        self["SSAM_ATTRIBUTE"] = {
            'COMMENT': 'Sea salt (acc mode)', 'UNITS': 'number density (cm-3)'}
        self["SSCM_ATTRIBUTE"] = {
            'COMMENT': 'Sea salt (coa mode)', 'UNITS': 'number density (cm-3)'}
        self["STCO_ATTRIBUTE"] = {
            'COMMENT': 'Stratus Continental', 'UNITS': 'g/m-3'}
        self["STMA_ATTRIBUTE"] = {
            'COMMENT': 'Stratus Maritime', 'UNITS': 'g/m-3'}
        self["SUNAZANGLE_ATTRIBUTE"] = {
            'COMMENT': 'Sun azimuth angle',
            'UNITS': 'degree (0-360 deg; east=90)'}
        self["SUNZENANGLE_ATTRIBUTE"] = {
            'COMMENT': 'Sun zenith angle', 'UNITS': 'degree'}
        self["SUSO_ATTRIBUTE"] = {
            'COMMENT': 'Sulphated droplets', 'UNITS': 'number density (cm-3)'}
        self["TIME_ATTRIBUTE"] = {'COMMENT': 'Hour Minute Second'}
        self["T_ATTRIBUTE"] = {'COMMENT': 'Temperature', 'UNITS': 'K'}
        self["VAPO_ATTRIBUTE"] = {
            'COMMENT': 'New Volcanic ash', 'UNITS': 'number density (cm-3)'}
        self["VOLA_ATTRIBUTE"] = {
            'COMMENT': 'OPAC Volcanic ash', 'UNITS': 'number density (cm-3)'}
        self["WASO_ATTRIBUTE"] = {
            'COMMENT': 'Water soluble', 'UNITS': 'number density (cm-3)'}
        self["ZENANGLE_ATTRIBUTE"] = {
            'COMMENT': 'Satellite zenith angle', 'UNITS': 'degree'}

    def setDefaultProfileAsciiInput(self):
        """Put here default values for all mandatory profile variables
           The profile shall be already filled with minimum values like
           array of pressure, temperature
           water vapour
           Defaults values such as:
           - Pressure, Temperatur, WaterVapour units
           - surface types and values (taken from bottom level"""
        if self["IDG"] is None:
            self["IDG"] = 1   # Scheme for Ice water content
        if self["ISH"] is None:
            self["ISH"] = 1   # Ice cristal shape

        if self["GAS_UNITS"] == None:
            self["GAS_UNITS"] = 2  # ppmv over moist air

        " Skin variables "
        if self["SKIN"]["T"] is None:
            self["SKIN"]["T"] = self["T"][-1]  # (K)
        if self["SKIN"]["SURFTYPE"] is None:
            self["SKIN"]["SURFTYPE"] = 1  # (0=Land, 1=Sea, 2=sea-ice)
        if self["SKIN"]["WATERTYPE"] is None:
            self["SKIN"]["WATERTYPE"] = 1  # (0=fresh water, 1=ocean water)
        if self["SKIN"]["SALINITY"] is None:
            self["SKIN"]["SALINITY"] = 37  # (%o)
        if self["SKIN"]["FASTEM"] is None:
            self["SKIN"]["FASTEM"] = numpy.array(
                [3., 5., 15., 0.1, 0.30000001])  # (5 parameters Land/sea-ice)
        if self["SKIN"]["FOAM_FRACTION"] is None:
            self["SKIN"]["FOAM_FRACTION"] = 0.0

        " 2m and 10m air variables "
        if self["S2M"]["T"] is None:
            self["S2M"]["T"] = self["T"][-1]  # (K)
        if self["S2M"]["Q"] is None:
            self["S2M"]["Q"] = self["Q"][-1]  # (ppmv)
        if self["S2M"]["P"] is None:
            self["S2M"]["P"] = self["P"][-1]  # (hPa)
        if self["S2M"]["U"] is None:
            self["S2M"]["U"] = 0  # (m/s)
        if self["S2M"]["V"] is None:
            self["S2M"]["V"] = 0  # (m/s)
        if self["S2M"]["WFETC"] is None:
            self["S2M"]["WFETC"] = 100000  # (m)

        " Simple cloud "
        if self["CTP"] is None:
            self["CTP"] = 500.0  # (hPa)
        if self["CFRACTION"] is None:
            self["CFRACTION"] = 0.0   # [0,1]    Clear sky is the default

        " Viewing geometry "
        if self["AZANGLE"] is None:
            self["AZANGLE"] = 0.  # (deg)
        if self["ELEVATION"] is None:
            self["ELEVATION"] = 0.  # (km)
        if self["SUNAZANGLE"] is None:
            self["SUNAZANGLE"] = 0.  # (deg)
        if self["SUNZENANGLE"] is None:
            self["SUNZENANGLE"] = 0.  # (deg)
        if self["ZENANGLE"] is None:
            self["ZENANGLE"] = 0.  # (deg)
        if self["SNOW_FRAC"] is None:
            self["SNOW_FRAC"] = 0.  # [0,1]
        if self["SOIL_MOISTURE"] == None:
            self["SOIL_MOISTURE"] = 0.  # [0,1]

        # Set latitude/longitude in between Lannion and Exeter,
        # in the "Manche"/"Channel"
        # Lannion is 48.750,  -3.470
        # Exeter is  50.726,  -3.476
        if self["LATITUDE"] is None:
            self["LATITUDE"] = 49.738  # (deg)
        if self["LONGITUDE"] is None:
            # (deg)    Lat/lon compatible with suftype==ocean
            self["LONGITUDE"] = -3.473

        " Magnetic field "
        if self["BE"] is None:
            self["BE"] = 0.3  # (Gauss)
        if self["COSBK"] is None:
            self["COSBK"] = 1.

        " Mislaneous "
        if self["ID"] is None:
            self["ID"] = "This is my profile"
        if self["DATE"] is None:
            self["DATE"] = numpy.array([2014, 04, 30])  # Year, Month, Day
        if self["TIME"] is None:
            self["TIME"] = numpy.array([12, 0, 0])     # Hour, Minute, Second
        self.setDefaultAttributes()
        if self["GAS_UNITS"] is None:
            self["GAS_UNITS"] = 2  # defaut case

    def checkProfileAsciiInput(self):
        if (self["P"] is None or self["T"] is None or self["Q"] is None):
            print "At least one of P, T, Q arrays is missing"
            raise IOError
        self["NLEVELS"] = numpy.shape(self["P"])[0]
        self["NLAYERS"] = self["NLEVELS"] - 1
        for cle in self.gas_list:
            if(self[cle] is not None):
                if (numpy.shape(self[cle])[0] != self["NLEVELS"]):
                    print "Bad array size for", cle
                    raise IOError
        for cle in self.aerosol_list:
            if(self[cle] is not None):
                if (numpy.shape(self[cle])[0] != self["NLAYERS"]):
                    print "Bad array size for", cle
                    raise IOError
        for cle in self.cloud_list:
            if(self[cle] is not None):
                if (numpy.shape(self[cle])[0] != self["NLAYERS"]):
                    print "Bad array size for", cle
                    raise IOError
        if self.anyCloud():
            if(self["CFRAC"] is not None):
                if (numpy.shape(self["CFRAC"])[0] != self["NLAYERS"]):
                    print "Bad array size for", cle
                    raise IOError
            else:
                print "Cloud fraction missing "
                raise IOError
            if(self["ICEDE"] is None):
                # ICEDE array is mandatory for cloudy profiles
                self["ICEDE"] = numpy.zeros(
                    self["NLAYERS"], dtype=self["P"].dtype)
        if "GAS_UNITS" not in self:
            self["GAS_UNITS"] = 2  # defaut case

    def loadh5(self, h5):
        # Load driven by the file content
        for cle in list(h5):
            if(cle in ['S2M', 'SKIN']):

                for scle in list(h5[cle]):
                    self[cle][scle + '_ATTRIBUTE'] = {}
                    self[cle][scle] = loadDset_arr(h5[cle + '/' + scle])
                    self[cle][
                        scle + '_ATTRIBUTE'] = getAttributes(h5[
                                                            cle + '/' + scle])
                    checkAttributes(self[cle][scle + '_ATTRIBUTE'])
            elif(cle == 'ID'):

                self[cle] = loadDset_str(h5[cle])
                self[cle + '_ATTRIBUTE'] = getAttributes(h5[cle])
                checkAttributes(self[cle + '_ATTRIBUTE'])
            else:

                self[cle] = loadDset_arr(h5[cle])
                self[cle + '_ATTRIBUTE'] = getAttributes(h5[cle])
                checkAttributes(self[cle + '_ATTRIBUTE'])

        if(self['AEROSOLS'] is not None):
            na = self['AEROSOLS'].shape[1]
            if(self['AEROSOLS'].ndim == 2):
                for i in range(na):
                    if (any(self['AEROSOLS'][:, i] != 0.)):
                        self[self.aerosol_list[i]] = self['AEROSOLS'][:, i]
                        self[self.aerosol_list[i] + '_ATTRIBUTE'] = {
                            'COMMENT': self.aerosol_long_list[i],
                            'UNITS': 'number density (cm-3)'}
                        if ('UNITS' in self['AEROSOLS_ATTRIBUTE'].keys()):
                            self[self.aerosol_list[
                                i] + '_ATTRIBUTE'] = {
                                        'COMMENT': self.aerosol_long_list[i]}
                            self[self.aerosol_list[
                                i] + '_ATTRIBUTE']['UNITS'] = self[
                                                'AEROSOLS_ATTRIBUTE']['UNITS']
                        else:
                            self[self.aerosol_list[i] + '_ATTRIBUTE'] = {
                                'COMMENT': self.aerosol_long_list[i],
                                'UNITS': 'number density (cm-3)'}
            else:
                for i in range(na):
                    if (numpy.any(self['AEROSOLS'][:, i, :] != 0.)):
                        self[self.aerosol_list[i]] = self['AEROSOLS'][:, i, :]
                        if ('UNITS' in self['AEROSOLS_ATTRIBUTE'].keys()):
                            self[self.aerosol_list[
                                i] + '_ATTRIBUTE'] = {
                                        'COMMENT': self.aerosol_long_list[i]}
                            self[self.aerosol_list[
                                i] + '_ATTRIBUTE']['UNITS'] = self[
                                                'AEROSOLS_ATTRIBUTE']['UNITS']
                        else:
                            self[self.aerosol_list[i] + '_ATTRIBUTE'] = {
                                'COMMENT': self.aerosol_long_list[i],
                                'UNITS': 'number density (cm-3)'}
            self['AEROSOLS'] = None
            self['AEROSOLS_ATTRIBUTE'] = {'COMMENT': 'none', 'UNITS': 'n/a'}

        if(self['CLOUD'] is not None):
            na = self['CLOUD'].shape[1]
            if(self['CLOUD'].ndim == 2):
                for i in range(na):
                    if (any(self['CLOUD'][:, i] != 0.)):
                        self[self.cloud_list[i]] = self['CLOUD'][:, i]
                        if ('UNITS' in self['CLOUD_ATTRIBUTE'].keys()):
                            self[self.cloud_list[
                                i] + '_ATTRIBUTE'] = {
                                        'COMMENT': self.cloud_long_list[i]}
                            self[self.cloud_list[
                                i] + '_ATTRIBUTE']['UNITS'] = self[
                                                    'CLOUD_ATTRIBUTE']['UNITS']
                        else:
                            self[self.cloud_list[i] + '_ATTRIBUTE'] = {
                                'COMMENT': self.cloud_long_list[i],
                                'UNITS': 'layer mean content(g/m3)'}
            else:
                for i in range(na):
                    if (numpy.any(self['CLOUD'][:, i, :] != 0.)):
                        self[self.cloud_list[i]] = self['CLOUD'][:, i, :]
                        if ('UNITS' in self['CLOUD_ATTRIBUTE'].keys()):
                            self[self.cloud_list[
                                i] + '_ATTRIBUTE'] = {
                                            'COMMENT': self.cloud_long_list[i]}
                            self[self.cloud_list[
                                i] + '_ATTRIBUTE']['UNITS'] = self[
                                                    'CLOUD_ATTRIBUTE']['UNITS']
                        else:
                            self[self.cloud_list[i] + '_ATTRIBUTE'] = {
                                'COMMENT': self.cloud_long_list[i],
                                'UNITS': 'layer mean content(g/m3)'}
            self['CLOUD'] = None
            self['CLOUD_ATTRIBUTE'] = {'COMMENT': 'none', 'UNITS': 'n/a'}
            if(self['ICEDE'] is None):
                # ICEDE array is mandatory for cloudy profiles
                self['ICEDE'] = numpy.zeros(
                    self['NLAYERS'], dtype=self['P'].dtype)
                self['ICEDE_ATTRIBUTE'] = {'COMMENT': 'Ice crystals diameter',
                                           'UNITS': 'microns', 'LBOUND': 1,
                                           'UBOUND': self['NLAYERS']}

        if self["GAS_UNITS"] is None:
            self["GAS_UNITS"] = 2  # defaut case
            self['GAS_UNITS_ATTRIBUTE'] = {
                'COMMENT': 'Units for gas profiles', 'UNITS': 'n/a'}

    def controlProfile(self):
        """ control some values of the profile regarding to min et max
            values """
        for label in self.minValue.keys():
            if label[:3] == "S2M":
                if self["S2M"][label[4:]] < self.minValue[label]:
                    self["S2M"][label[4:]] = self.minValue[label]
                if self["S2M"][label[4:]] > self.maxValue[label]:
                    self["S2M"][label[4:]] = self.maxValue[label]
            else:
                if label[:4] == "SKIN":
                    if self["SKIN"][label[5:]] < self.minValue[label]:
                        self["SKIN"][label[5:]] = self.minValue[label]
                    if self["SKIN"][label[5:]] > self.maxValue[label]:
                        self["SKIN"][label[5:]] = self.maxValue[label]
                else:
                    if self[label] < self.minValue[label]:
                        self[label] = self.minValue[label]
                    if self[label] > self.maxValue[label]:
                        self[label] = self.maxValue[label]

    def loadProfileNumber(self, file, iprof):
        """
        Load profile number iprof from file
        Opens and closes the file
        """
        # opens the HDF5 file read only
        f = h5py.File(file, 'r')
        dir = '/PROFILES/{:0>4}/'.format(iprof)

        # get the Dataset
        h5 = f[dir]

        # Load profile
        self.loadh5(h5)

        # Close HDF file
        f.close()

    def loadProfileAscii(self, fname):
        """
        Load an ASCII profile from file fname
        """
        try:
            execfile(fname)
        except SyntaxError:
            print "Error Reading ASCII profile file"
            return 1
        try:
            self.checkProfileAsciiInput()
        except IOError:
            print "Error Checking ASCII profile file"
            return 1
        self.setDefaultProfileAsciiInput()

        return 0

    def display(self):
        # Display profile content to terminal
        for cle in sorted(self.keys()):
            print cle, '=>', self[cle], type(self[cle])

    def saveh5(self, h5, path):
        """Save profile class to an HDF5 file h5 under path directory
        Do not save empty arrays
        Attributes are saved at the same time datasets are saved
        """
        list_of_names = []
        h5.visit(list_of_names.append)

        if self.anyAerosol():
            na = len(self.aerosol_list)
            self['AEROSOLS'] = numpy.zeros(
                (self['NLAYERS'], len(self.aerosol_list)),
                dtype=self['P'].dtype)
            self['AEROSOLS_ATTRIBUTE'] = {
                'COMMENT': 'Aerosols', 'UNITS': 'number density (cm-3)'}

            for cle in self.aerosol_list:
                if(self[cle] is None):
                    continue
                else:
                    i = self.aerosol_list.index(cle)
                    self['AEROSOLS'][:, i] = self[cle]

        if self.anyCloud():
            na = len(self.cloud_list)
            self['CLOUD'] = numpy.zeros(
                (self['NLAYERS'], len(self.cloud_list)), dtype=self['P'].dtype)
            self['CLOUD_ATTRIBUTE'] = {
                'COMMENT': 'Cloud water/ice - IR only', 'UNITS': 'g/m-3'}

            for cle in self.cloud_list:
                if(self[cle] is None):
                    continue
                else:
                    i = self.cloud_list.index(cle)
                    self['CLOUD'][:, i] = self[cle]
            if(self['ICEDE'] is None):
                # ICEDE array is mandatory for cloudy profiles
                self['ICEDE'] = numpy.zeros(
                    self['NLAYERS'], dtype=self['P'].dtype)
                self['ICEDE_ATTRIBUTE'] = {'COMMENT': 'Ice crystals diameter',
                                           'UNITS': 'microns', 'LBOUND': 1,
                                           'UBOUND': self['NLAYERS']}

        for cle in self.keys():

            if(re.search('_ATTRIBUTE', cle) or self[cle] is None):
                # do not save empty members
                continue

            if(cle in ['S2M', 'SKIN']):
                # S2M and SKIN are HDF5 sub-directories
                for scle in self[cle].keys():
                    if(
                     re.search('_ATTRIBUTE', scle) or self[cle][scle] is None):
                        continue
                    # create the dataset and save attributes
                    cpath = cleanpath(path + '/' + cle + '/' + scle)
                    if cpath in list_of_names:
                        del h5[cpath]
                    dset = h5.create_dataset(cpath, data=self[cle][scle])
                    saveAttributes(h5[cpath], self[cle][scle + '_ATTRIBUTE'])

            elif(cle in self.aerosol_list or cle in self.cloud_list):
                # do not save each individual aerosol array but "aerosol" 2
                # dimensions array
                # HDF5 file only contains the RTTOV arrays
                continue

            else:
                # create the dataset and save attributes
                cpath = cleanpath(path + '/' + cle)
                if cpath in list_of_names:
                    del h5[cpath]
                dset = h5.create_dataset(cpath, data=self[cle])
                saveAttributes(h5[cpath], self[cle + '_ATTRIBUTE'])

        # ici on peut supprimer le tableau a deux dimensions des aerosols et
        # nuages
        self['AEROSOLS'] = None
        self['AEROSOLS_ATTRIBUTE'] = {'COMMENT': 'none', 'UNITS': 'n/a'}
        self['CLOUD'] = None
        self['CLOUD_ATTRIBUTE'] = {'COMMENT': 'none', 'UNITS': 'n/a'}

    def addGas(self, gas, *conc):
        """
        on ajoute un gas de nom gas a la bonne dimension nlevels
        et avec le meme type de données que le tableau de pression
        """
        # create numpy array filled with zeros, size is based on number of
        # pressure levels
        self[gas] = numpy.zeros(self['NLEVELS'], dtype=self['P'].dtype)
        self[gas + '_ATTRIBUTE'] = {'COMMENT': gas, 'UNITS': 'ppmv',
                                    'LBOUND': 1, 'UBOUND': self['NLEVELS']}
        # Fills with given concentration or default value
        if(conc):
            self[gas] += conc
        else:
            if self['GAS_UNITS'] in [-1, 0, 2]:
                self[gas] += self.defaultppmv[gas]
            else:
                self[gas] += fct_ppmv_to_mixratio(self.defaultppmv[gas], gas)

    def removeGas(self, gas):
        """
        on enleve un gas de nom gas
        """
        self[gas] = None
        self[gas + '_ATTRIBUTE'] = {'COMMENT': 'none', 'UNITS': 'n/a'}

    def removeAllButGas(self, gas):
        """
        on enleve tous les gas sauf celui de nom gas
        """
        for g in self.gas_list:
            if(not g == gas):
                self[g] = None
                self[g + '_ATTRIBUTE'] = {'COMMENT': 'none', 'UNITS': 'n/a'}

    def removeAllAerosol(self):
        """
        on enleve tous les aerosols
        """
        for i in self.aerosol_list:
            self[i] = None
            self[i + '_ATTRIBUTE'] = {'COMMENT': 'none', 'UNITS': 'n/a'}

        self['AEROSOLS'] = None
        self['AEROSOLS_ATTRIBUTE'] = {'COMMENT': 'none', 'UNITS': 'n/a'}

    def removeAerosol(self, aerosol):
        """
        on enleve un aerosol
        """
        self[aerosol] = None
        self[aerosol + '_ATTRIBUTE'] = {'COMMENT': 'none', 'UNITS': 'n/a'}

    def addAerosol(self, aerosol, *conc):
        """
        on ajoute un aerosol de nom aerosol a la bonne dimension nlayers
        et avec le meme type de données que le tableau de pression
        """
        # create numpy array filled with zeros, size is based on number of
        # pressure levels
        i = self.aerosol_list.index(aerosol)
        self[aerosol] = numpy.zeros(self['NLAYERS'], dtype=self['P'].dtype)
        self[aerosol + '_ATTRIBUTE'] = {'COMMENT': self.aerosol_long_list[i],
                                        'UNITS': 'number density (cm-3)',
                                        'LBOUND': 1, 'UBOUND': self['NLAYERS']}
        if(aerosol == 'VOLA' or aerosol == 'VAPO'):
            # volcanic ashes
            pmin = 100.
            pmax = 300.

        else:
            pmin = 600.
            pmax = 1000.

        for i in range(self['NLAYERS']):
            # here we use presure levels indices to fill cloud LAYERS.
            if(self['P'][i] > pmin and self['P'][i] < pmax):
                if(conc):
                    self[aerosol][i] += conc
                else:
                    self[aerosol][i] += self.defaultaerdensity

    def replaceByAerosolClim(self, clim):
        """
        on remplace les aerosols par une valeur climatologique
        avec clim de 1 a 10:
        1  -->Continental clean
        2  -->Continental average
        3  -->Continental polluted
        4  -->Urban
        5  -->Desert
        6  -->Maritime clean
        7  -->Maritime polluted
        8  -->Maritime tropical
        9  -->Arctic
        10 -->Antarctic
        """
        # create numpy array filled with zeros, size is based on number of
        # pressure levels
        tabaer = None
        tabaer, err = rttov_gui_aer_clim_prof(self['P'], self['T'], self['Q'],
                                              self['NLEVELS'],
                                              self['LATITUDE'],
                                              self['ELEVATION'],
                                              1.0,
                                              self['NLEVELS'])

        self.removeAllAerosol()
        self['AEROSOLS'] = tabaer[:, clim - 1, :]
        na = self['AEROSOLS'].shape[1]
        for i in range(na):
            if (any(self['AEROSOLS'][:, i] != 0.)):
                self[self.aerosol_list[i]] = self['AEROSOLS'][:, i]
                self[self.aerosol_list[i] + '_ATTRIBUTE'] = {
                    'COMMENT': self.aerosol_long_list[i],
                    'UNITS': 'number density (cm-3)'}
        # ici on peut supprimer le tableau a deux dimensions des aerosols
        self['AEROSOLS'] = None
        self['AEROSOLS_ATTRIBUTE'] = {'COMMENT': 'none', 'UNITS': 'n/a'}

    def anyAerosol(self):
        for cle in self.aerosol_list:
            if(self[cle] is None):
                continue
            else:
                return True
        return False

    def removeCloud(self, cloud):
        """
        remove a cloud
        """
        self[cloud] = None
        self[cloud + '_ATTRIBUTE'] = {'COMMENT': 'none', 'UNITS': 'n/a'}
        if not self.anyCloud():
            self['CFRAC'] = None
            self['CFRAC_ATTRIBUTE'] = {'COMMENT': 'none', 'UNITS': 'n/a',
                                       'LBOUND': 0, 'UBOUND': 0}
        if not self.anyCloud() and self['CLW'] is not None:
            self['CLW'] = None
            self['CLW_ATTRIBUTE'] = {'COMMENT': 'none', 'UNITS': 'n/a',
                                     'LBOUND': 0, 'UBOUND': 0}
        if not self.anyCloud() and self['ICEDE'] is not None:
            self['ICEDE'] = None
            self['ICEDE_ATTRIBUTE'] = {'COMMENT': 'none', 'UNITS': 'n/a',
                                       'LBOUND': 0, 'UBOUND': 0}

    def addCloud(self, cloud, *conc):
        """
        on ajoute un cloud de nom cloud a la bonne dimension nlayers
        et avec le meme type de données que le tableau de pression
        """
        # create numpy array filled with zeros, size is based on number of
        # pressure levels
        i = self.cloud_list.index(cloud)
        self[cloud] = numpy.zeros(self['NLAYERS'], dtype=self['P'].dtype)
        self[cloud + '_ATTRIBUTE'] = {'COMMENT': self.cloud_long_list[i],
                                      'UNITS': 'layer maen content(g/m3)',
                                      'LBOUND': 1, 'UBOUND': self['NLAYERS']}

        if(self['CFRAC'] is None):
            self['CFRAC'] = numpy.zeros(self['NLAYERS'], dtype=self['P'].dtype)
            self['CFRAC_ATTRIBUTE'] = {
                                'COMMENT': 'Cloud fractional cover (0:1)',
                                'UNITS': '0:1',
                                'LBOUND': 1, 'UBOUND': self['NLAYERS']}

        if(self['ICEDE'] is None):
                # ICEDE array is mandatory for cloudy profiles
            self['ICEDE'] = numpy.zeros(self['NLAYERS'], dtype=self['P'].dtype)
            self['ICEDE_ATTRIBUTE'] = {'COMMENT': 'Ice crystals diameter',
                                       'UNITS': 'microns',
                                       'LBOUND': 1, 'UBOUND': self['NLAYERS']}

        if(cloud == 'STCO' or cloud == 'STMA'):
            # Stratus
            pmin = 850.
            pmax = 1000.
        elif(cloud == 'CUCC' or cloud == 'CUCP' or cloud == 'CUMA'):
            # Cumulus
            pmin = 600.
            pmax = 800.
        else:
            # Cirrus
            pmin = 450.
            pmax = 550.

        for i in range(self['NLAYERS']):
            # here we use presure levels indices to fill cloud LAYERS.
            if(self['P'][i] > pmin and self['P'][i] < pmax):
                self[cloud][i] = self.defaultclddensity
                if (self['CFRAC'][i] == 0.):
                    self['CFRAC'][i] = 0.20

    def anyCloud(self):
        for cle in self.cloud_list:
            if(self[cle] is None):
                continue
            else:
                return True
        return False

    def to1DvarVect(self):
        """ return an 1Dvar vector """
        """ 1Dvar vector are only on 43, 51 or 54 levels"""
        """ a vector contains only T values, lnq bottom lnq vales ,
            Tsurf, lnq surf and Tskin"""
        """ temperature, and ln q are stored from top of atmosphere
            to ground """

        nlevels = len(self['T'])
        if (nlevels) not in (43, 51, 54):
            print ("WARNING :", "nlevels " +
                   str(nlevels) + "not in (43,51,54) ")
        veclen = nlevels + 29 + 3
        # 1-54=temp (K), 55-83=lnq (g/kg) (bottom 29 levels only), 84=Tsurf,
        # 85=lnq surf, 86=Tskin
        vector = numpy.zeros((veclen), dtype=float)
        for i in range(nlevels):
            vector[i] = self['T'][i]
        for i in range(29):
            vector[nlevels + i] = self['Q'][nlevels - 29 + i]
        vector[nlevels + 29] = self['S2M']['T']
        vector[nlevels + 29 + 1] = self['S2M']['Q']
        vector[nlevels + 29 + 2] = self['SKIN']['T']

        return vector

    def to1DvarBackgroundFile(self, filename):
        """ write an 1Dvar Background File """
        """ 1Dvar work only on 43, 51 or 54 levels"""
        """ must respect File description """
        """ https://nwpsaf.eu/deliverables/nwpsaf_1dvar/
            nwpsaf-mo-ud-032_NWPSAF_1DVar_Manual.html#aux """
        f = open(filename, "w")
        nlevels = len(self['T'])
        if (nlevels) not in (43, 51, 54):
            print ("WARNING :", "nlevels " +
                   str(nlevels) + "not in (43,51,54) \n")
        f.write("Backgrounnd file generated from\n")
        f.write(self["ID"] + "\n")
        f.write(str(nlevels) + " fixed layers")
        for k in range(0, 6):
            f.write("\n")
        f.write(
            "x------------------- End of Header---------------------------x\n")
        f.write("\n")
        f.write("No. Background Profiles:            1\n")
        f.write("No. of Levels/Profile:             " + str(nlevels) + "\n")
        f.write("Unit for q:                     1  (1=ppmv, 2=kg/kg, 3=RH)\n")
        f.write(
            "x------------------------------------------------------------x\n")
        f.write("Profile #    1 Follows\n")
        f.write(
            "x------------------------------------------------------------x\n")
        for i in range(nlevels - 1, -1, -1):
            sp = ("%10.5f" % self["P"][i]).rjust(13)
            st = ("%9.5f" % self["T"][i]).rjust(13)
            sq = ("%f" % self["Q"][i]).rjust(13)
            so3 = ("%f" % self["O3"][i]).rjust(13)
            f.write(sp + st + sq + so3)
            f.write("\n")
        f.write("Surface Temperature (K):" +
                ("%9.5f" % self["S2M"]["T"]).rjust(18))
        f.write("\n")
        f.write("Surface Humidity (ppmv)  :" +
                ("%f" % self["S2M"]["Q"]).rjust(18))
        f.write("\n")
        f.write("Skin Temperature (K):" + ("%9.5f" %
                                           self["SKIN"]["T"]).rjust(18))
        f.write("\n")
        f.write("Surface Pressure (hPa):" + ("%10.5f" %
                                             self["S2M"]["P"]).rjust(18))
        f.write("\n")
        f.write("10m U-Wind (m/s):" + ("%f" % self["S2M"]["U"]).rjust(18))
        f.write("\n")
        f.write("10m U-Wind (m/s):" + ("%f" % self["S2M"]["V"]).rjust(18))
        f.write("\n")
        f.close()
if __name__ == '__main__':

    file = "../rttov_tests/cldaer101lev_allgas.H5"

    # Open file ReadOnly
    f = h5py.File(file, 'r')

    # get the Dataset
    h5 = f['/PROFILES/0001/']

    # p1 is a rttov_hdf_mod Profile instance
    p1 = Profile()

    # Load profile
    p1.loadh5(h5)

    # Close HDF file
    f.close()

    # Display profile
    p1.display()

    ofile = "../rttov_tests/output_profile.H5"
    of = h5py.File(ofile, 'w')
    p1.saveh5(of, "/PROFILES/0001/")
    of.close()

    # p1.replaceByAerosolClim(7)
    # p1.addAerosol('SSAM')
    # p1.addCloud('CUMA')

    # Plots temperature and gas concentrations
    import matplotlib.pyplot as plt
    plt.plot(p1['T'], p1['P'], 'ro-')
    plt.xlabel(p1['T_ATTRIBUTE']['COMMENT'] +
               '   (' + p1['T_ATTRIBUTE']['UNITS'] + ')')
    plt.ylim(1100, 0.01)
    plt.yscale('log')
    plt.ylabel(p1['P_ATTRIBUTE']['COMMENT'] +
               '  ' + p1['P_ATTRIBUTE']['UNITS'])
    plt.title(p1['ID'] + ' TEMPERATURE')
    plt.show()

    if(p1.anyCloud()):
        nlevels = p1['NLEVELS']
        player = (p1['P'][0:nlevels - 1] + p1['P'][1:nlevels]) / 2.0
        colors = ['b', 'g', 'r', 'c', 'm', 'y', 'k']
        i = 0
        for cle in p1.cloud_list:
            if(p1[cle] is None):
                continue
            else:
                color = colors[i]
                x = p1[cle]
                plt.plot(x, player, color + 'o-',
                         label=p1[cle + '_ATTRIBUTE']['COMMENT'])
                i += 1

        x = p1['CFRAC'][:]
        plt.plot(x, player, 'ko-', label='Cloud Fraction')
        plt.xlabel('Cloud Liquid water/ Cloud Ice Water')
        plt.xscale('log')
        plt.ylim(1100, 10)
        plt.yscale('log')
        plt.ylabel(p1['P_ATTRIBUTE']['COMMENT'] +
                   '   (' + p1['P_ATTRIBUTE']['UNITS'] + ')')
        plt.legend()
        plt.title(p1['ID'] + ' CLOUDS')
        plt.show()

    if (p1.anyAerosol()):

        nlevels = p1['NLEVELS']
        player = (p1['P'][0:nlevels - 1] + p1['P'][1:nlevels]) / 2.0
        colors = ['b', 'g', 'r', 'c', 'm', 'y']
        i = 0
        for cle in p1.aerosol_list:
            if(p1[cle] is None):
                continue
            else:
                color = colors[i % len(colors)]
                if (i >= len(colors)):
                    pattern = 'o'
                else:
                    pattern = 'v'
                x = p1[cle]
                plt.plot(x, player, color + pattern + '-',
                         label=p1[cle + '_ATTRIBUTE']['COMMENT'])
                i += 1
        plt.xlabel('Aerosols (number density cm^-3)')
        plt.xscale('log')
        plt.ylim(1100, 5)
        plt.yscale('log')
        plt.ylabel(p1['P_ATTRIBUTE']['COMMENT'] +
                   '   (' + p1['P_ATTRIBUTE']['UNITS'] + ')')
        plt.legend()
        plt.title(p1['ID'] + ' AEROSOLS')
        plt.show()

    colors = ['b', 'g', 'r', 'c', 'm', 'y', 'k']
    i = 0
    for gas in Profile.gas_list:
        if(p1[gas] is None):
            continue
        color = colors[i]
        i += 1
        plt.plot(p1[gas], p1['P'], color + 'o-', label=gas)
    plt.xlabel('gas concentrations (ppmv)')
    plt.xscale('log')
    plt.ylim(1100, 0.01)
    plt.yscale('log')
    plt.ylabel(p1['P_ATTRIBUTE']['COMMENT'] +
               '   (' + p1['P_ATTRIBUTE']['UNITS'] + ')')
    plt.legend()
    plt.title(p1['ID'] + ' GASES')
    plt.show()
