'''
:Copyright: 2015, EUMETSAT, All Rights Reserved.

This software was developed within the context of
the EUMETSAT Satellite Application Facility on
Numerical Weather Prediction (NWP SAF), under the
Cooperation Agreement dated 25 November 1998, between
EUMETSAT and the Met Office, UK, by one or more partners
within the NWP SAF. The partners in the NWP SAF are
the Met Office, ECMWF, KNMI and MeteoFrance.
'''

from __future__ import absolute_import, print_function
import collections
from .decorator import add_descriptors_opts, lock_attributes

(BAND1, BAND2, BAND3) = range(1, 4)
(INDEX1, INDEX2, INDEX3, INDEX4) = range(1, 5)
(FASTEM1, FASTEM2, FASTEM3, FASTEM4, FASTEM5, FASTEM6) = range(1, 7)


@lock_attributes
@add_descriptors_opts(['wrapper', 'config', 'interpolation', 'rt_all', 'rt_mw',
                       'rt_ir', 'rt_ir_pc'])
class Options(object):
    '''The Options class holds and manages the RTTOV options.'''

    # Defaults for wrapper options
    _defaults_wrapper = {
        'nthreads': 1,
        'nprofs_per_call': 1,
        'verbose_wrapper': False,
        'check_opts': True,
        'store_rad': False,
        'store_rad2': False,
        'store_trans': False,
    }
    # Defaults for RTTOV configuration: config
    _defaults_config = {
        'apply_reg_limits': False,
        'verbose': True,
        'do_checkinput': True,
    }
    # Defaults for RTTOV configuration: interpolation
    _defaults_interpolation = {
        'addinterp': True,
        'spacetop': True,
        'lgradp': False,
        'reg_limit_extrap': False,
        'interp_mode': 1,
    }
    # Defaults for RTTOV configuration: rt_all
    _defaults_rt_all = {
        'addrefrac': False,
        'switchrad': False,
        'use_q2m': True,
        'do_lambertian': False,
    }
    # Defaults for RTTOV configuration: rt_mw
    _defaults_rt_mw = {
        'clw_data': False,
        'supply_foam_fraction': False,
        'fastem_version': FASTEM5,
    }
    # Defaults for RTTOV configuration: rt_ir
    _defaults_rt_ir = {
        'addsolar': False,
        'do_nlte_correction': False,
        'addaerosl': False,
        'addclouds': False,
        'cldstr_threshold': -1.,
        'cldstr_simple': False,
        'user_aer_opt_param': False,
        'user_cld_opt_param': False,
        'ozone_data': False,
        'co2_data': False,
        'n2o_data': False,
        'co_data': False,
        'ch4_data': False,
    }
    # Defaults for RTtOV configuration: rt_pc
    _defaults_rt_ir_pc = {
        'addpc': False,
        'addradrec': False,
        'ipcreg': BAND1,
        'ipcbnd': INDEX1,
    }
    # Override the naming convention for some of the descriptors
    _naming_override = {
        'addpc': 'AddPC',
        'clw_data': 'CLWData',
        'co2_data': 'CO2Data',
        'ch4_data': 'CH4Data',
        'co_data': 'COData',
        'n2o_data': 'N2OData',
    }

    def __init__(self):
        # Call the function define by the decorator in order to initialise
        # the defaults
        self._initDefaults()

    def _initDefaults(self):
        '''May be overwritten by the class decorator.'''
        pass

    @staticmethod
    def _prepareValue(val):
        '''Transform option's values if needed.'''
        if isinstance(val, bool):
            return str(int(val))
        else:
            return str(val)

    def _prepareOpts(self):
        '''Prepare the option string.'''
        opts_processed = collections.OrderedDict()
        for confkey, thedict in {'opts%config': self._config,
                                 'opts%interpolation': self._interpolation,
                                 'opts%rt_all': self._rt_all,
                                 'opts%rt_mw': self._rt_mw,
                                 'opts%rt_ir': self._rt_ir,
                                 'opts%rt_ir%pc': self._rt_ir_pc}.iteritems():
            for varkey in sorted(thedict.keys()):
                opts_processed[confkey + '%' + varkey] = \
                    self._prepareValue(thedict[varkey])
        for varkey in sorted(self._wrapper.keys()):
            opts_processed[varkey] = self._prepareValue(self._wrapper[varkey])
        if not self.AddClouds:
            del opts_processed['opts%rt_ir%cldstr_threshold']
        if not self.AddPC:
            del opts_processed['opts%rt_ir%pc%ipcbnd']
            del opts_processed['opts%rt_ir%pc%ipcreg']
        return opts_processed

    def defineStrOptions(self):
        '''Returns an option string understandable by RTTOV's wrapper.'''
        return ' '.join(['{} {!s}'.format(key, val)
                         for key, val in self._prepareOpts().iteritems()])

    def __str__(self):
        return '\n'.join(['{:40s}  {!s}'.format(key, val)
                          for key, val in self._prepareOpts().iteritems()])
