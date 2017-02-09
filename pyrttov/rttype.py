'''
:Copyright: 2015, EUMETSAT, All Rights Reserved.

This software was developed within the context of
the EUMETSAT Satellite Application Facility on
Numerical Weather Prediction (NWP SAF), under the
Cooperation Agreement dated 25 November 1998, between
EUMETSAT and the Met Office, UK, by one or more partners
within the NWP SAF. The partners in the NWP SAF are
the Met Office, ECMWF, KNMI and MeteoFrance.

This module defines abstract classes to generate data types.
'''

import numpy as np


# Numpy types expected by the RTTOV f2py wrapper
wrapfloat = np.float64
wrapint = np.int32


class _genericEnum(int):
    '''Type that do something close to a C/C++ enum.'''

    _enum_dict = {}

    def __new__(cls, value):
        if (isinstance(value, (int, np.int16, np.int32, np.int64)) and
                int(value) in cls._enum_dict.itervalues()):
            return int.__new__(cls, value)
        if isinstance(value, basestring) and value in cls._enum_dict:
            return int.__new__(cls, cls._enum_dict[value])
        raise ValueError("Incorect value for a {}: {!s}".format(cls.__name__,
                                                                value))

    @property
    def name(self):
        for key, item in self._enum_dict.iteritems():
            if item == self:
                return key


class gasUnitType(_genericEnum):
    '''Type that represents a gas unit.'''

    _enum_dict = {'unknown': -2, 'ppmv_dry': -1,
                  'compatibility_mode': 0,
                  'kg_per_kg': 1, 'ppmv_wet': 2}


class itemIdType(_genericEnum):
    '''
    Type that represents a mapping between gas/item name and its internal id.
    '''

    _enum_dict = {'Q': 1, 'O3': 2, 'CO2': 3, 'N2O': 4, 'CO': 5, 'CH4': 6,
                  'CLW': 15, 'CFRAC': 20, 'STCO': 21, 'STMA': 22, 'CUCC': 23,
                  'CUCP': 24, 'CUMA': 25, 'CIRR': 30, 'ICEDE': 31, 'INSO': 41,
                  'WASO': 42, 'SOOT': 43, 'SSAM': 44, 'SSCM': 45, 'MINM': 46,
                  'MIAM': 47, 'MICM': 48, 'MITR': 49, 'SUSO': 50, 'VOLA': 51,
                  'VAPO': 52, 'ASDU': 53}
