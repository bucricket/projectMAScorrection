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
import sys
import copy
import collections
import numpy as np

import rttov_wrapper_f2py as rtwrap

from pyrttov.option import Options
from pyrttov.profile import Profiles
from pyrttov.profile import ItemDescriptorNaming as _ItemDescriptorNaming
from pyrttov.decorator import property_ro, lock_attributes, add_descriptors_gasesK
from pyrttov.descriptor import FilepathDescriptorRWD as _FPathDescRWD
from pyrttov.descriptor import DirpathDescriptorRWD as _DPathDescRWD
from pyrttov.descriptor import GenericDescriptorRO as _DescRO
from pyrttov.rttype import wrapfloat, wrapint, itemIdType


# No automatic exports
__all__ = []
__version__ = 0.1.0


class RttovError(Exception):
    pass


class RttovInternalError(RttovError):
    pass


class RttovWrapperError(RttovError):
    pass


@lock_attributes
@add_descriptors_gasesK(_ItemDescriptorNaming)
class Rttov(object):

    '''Main wrapper class for RTTOV.'''

    # Because dict is a mutable type, it will be shared between
    # child and ancestor classes
    _shared_classvar = {
        'instance_maxcount': 100,
        'instance_count': 0}

    @classmethod
    def _instance_add(cls):
        if (cls._shared_classvar['instance_count'] ==
                cls._shared_classvar['instance_maxcount']):
            errmsg = "Maximum count of RTTOV instances was reached."
            raise RttovInternalError(errmsg)
        cls._shared_classvar['instance_count'] += 1

    @classmethod
    def _instance_del(cls):
        cls._shared_classvar['instance_count'] -= 1

    def __new__(self, *args, **kwargs):
        self._instance_add()
        return object.__new__(self, *args, **kwargs)

    def __init__(self):
        self._coefs_conf = {}
        self._nchannels = 0
        self._instid = -1
        self._debug = False
        self._options = Options()
        self._brdfatlasloaded = False
        self._iratlasloaded = False
        self._mwatlasloaded = False
        del self.SurfEmisRefl
        self._profiles = None
        self._bBasics = dict()
        self._bTrans = dict()
        self._bRad = dict()
        self._bRad2 = dict()

    def __del__(self):
        self.deallocBrdfAtlas()
        self.deallocIrAtlas()
        self.deallocMwAtlas()
        if self._instid > 0:
            self._printVerbWrap("Deallocating this inst_id.")
            err = rtwrap.rttov_drop_inst(self._instid)
            self._errCheck(err, 'Error in rttov_drop_inst')
        self._instance_del()

    def _printDbg(self, msg):
        if self._debug:
            sys.stderr.write("DEBUG: " + msg + "\n")

    def _printVerbWrap(self, msg):
        if self.Options.VerboseWrapper or self._debug:
            sys.stderr.write(msg + "\n")

    @staticmethod
    def _errCheck(err, msg, exc=RttovInternalError):
        if err != 0:
            raise exc(msg)

    def _get_options(self):
        return self._options

    def _set_options(self, obj):
        if not isinstance(obj, Options):
            raise ValueError("Wrong type for options.")
        self._options = obj

    dstr = (":class:`.Options` object where RTTOV options " +
            "can be tuned." +
            "\n\nThis is a Read/Write attribute of type :class:`.Options`.")
    Options = property(_get_options, _set_options, doc=dstr)

    @property_ro
    def CoeffsLoaded(self):
        """Whether or not :meth:`loadInst` has been called successfully."""
        return self._instid >= 0

    @property_ro
    def BrdfAtlasLoaded(self):
        """Whether or not a BRDF atlas is loaded."""
        return self._brdfatlasloaded

    @property_ro
    def IrAtlasLoaded(self):
        """Whether or not an IR atlas is loaded."""
        return self._iratlasloaded

    @property_ro
    def MwAtlasLoaded(self):
        """Whether or not a MW atlas is loaded."""
        return self._mwatlasloaded

    @property_ro
    def Nchannels(self):
        """The number of channels available to this instance."""
        return self._nchannels

    FileCoef = _FPathDescRWD('coefs_conf', 'file_coef',
                             doc="Path to the RTTOV coefficients.")

    FileSccld = _FPathDescRWD('coefs_conf', 'file_sccld',
                              doc="Path to the RTTOV cloud coefficients.")

    FileScaer = _FPathDescRWD('coefs_conf', 'file_scaer',
                              doc="Path to the RTTOV aerosol coefficients.")

    FilePcCoef = _FPathDescRWD('coefs_conf', 'file_pccoef',
                               doc="Path to the RTTOV PC coefficients.")

    BrdfAtlasPath = _DPathDescRWD('coefs_conf', 'brdf_atlas_path',
                                  doc="Path to the BRDF atlas files.")

    EmisAtlasPath = _DPathDescRWD('coefs_conf', 'emis_atlas_path',
                                  doc="Path to the emissivity atlas files.")

    def loadInst(self, channels=None):
        '''Load instrument for a list of channels.

        At least :data:`FileCoef` must have been set previously.

        :param channels: list of channels to be loaded (all if omitted).
        :type channels: List of int or numpy.ndarray
        '''
        optionsloc = copy.deepcopy(self.Options)

        if self.FileCoef is None:
            raise RttovInternalError('The FileCoef attribute have to be set')
        OptStack = ['file_coef', self.FileCoef]

        if self.FileSccld is not None:
            OptStack.extend(['file_sccld', self.FileSccld])
            optionsloc.AddClouds = True

        if self.FileScaer is not None:
            OptStack.extend(['file_scaer', self.FileScaer])
            optionsloc.AddAerosl = True

        if self.FilePcCoef is not None:
            OptStack.extend(['file_pccoef', self.FilePcCoef])
            optionsloc.AddPC = True

        if self.CoeffsLoaded:
            errmsg = "The coefficients where already loaded once..."
            raise RttovInternalError(errmsg)

        OptStack.append(optionsloc.defineStrOptions())
        self._printDbg(' '.join(OptStack))

        if channels is None:
            npchannels = np.array([0, ], dtype=wrapint)
        elif isinstance(channels, collections.Iterable):
            npchannels = np.array(channels, dtype=wrapint)
        else:
            errmsg = "Incorrect type {!s} for channels".format(type(channels))
            raise TypeError(errmsg)

        self._instid = rtwrap.rttov_load_inst(' '.join(OptStack), npchannels)
        if self._instid > 0:
            err, self._nchannels = rtwrap.rttov_get_coef_val_i0(self._instid,
                                                                "NCHANNELS")
            if err != 0:
                raise RttovInternalError('Error getting number of channels')
        else:
            raise RttovInternalError('Error in rttov_load_inst')

        tmpmsg = "Load successful >>>>> inst_id : {0:d}, nchannels : {1:d}."
        self._printVerbWrap(tmpmsg.format(self._instid, self._nchannels))

    @property_ro
    def CoeffsNlevels(self):
        """Number of levels of the coefficient file."""
        if self.CoeffsLoaded:
            err, l_nlevels = rtwrap.rttov_get_coef_val_i0(self._instid,
                                                          "SIZE_REF_PRFL_P")
            self._errCheck(err, 'Error getting the number of levels')
            return l_nlevels
        else:
            return 0

    @property_ro
    def RefPressures(self):
        """Pressure levels of the coefficient file."""
        if self.CoeffsLoaded:
            err, l_p = rtwrap.rttov_get_coef_val_r1(self._instid,
                                                    "REF_PRESSURE",
                                                    self.CoeffsNlevels)
            self._errCheck(err, 'Error getting pressures from coefficient')
            return l_p
        else:
            return None

    @property_ro
    def WaveNumbers(self):
        """Channel's central wavenumbers of the coefficient file."""
        if self.CoeffsLoaded:
            err, l_p = rtwrap.rttov_get_coef_val_r1(self._instid,
                                                    "WAVENUMBERS",
                                                    self.Nchannels)
            self._errCheck(err, 'Error getting channel central wavenumbers')
            return l_p
        else:
            return None

    def updateOptions(self):
        """Update RTTOV options for the currently loaded instrument."""
        if not self.CoeffsLoaded:
            raise RttovWrapperError('Coefficients not loaded')

        # This may prevent crash
        if not self.FileSccld:
            self.Options.AddClouds = False
        if not self.FileScaer:
            self.Options.AddAerosl = False

        str_options = self.Options.defineStrOptions()
        self._printDbg('Set RTTOV options: {}'.format(str_options))
        err = rtwrap.rttov_set_options(self._instid, str_options)
        self._errCheck(err, 'Error in rttov_set_options')

    def printOptions(self):
        """Print RTTOV options for the currently loaded instrument."""
        if not self.CoeffsLoaded:
            raise RttovWrapperError('Coefficients not loaded')
        err = rtwrap.rttov_print_options(self._instid)
        self._errCheck(err, 'Error in rttov_print_options')

    def _GenericAtlasSetupCheck(self, atlas, atlas_p, month):
        atlas_uc = atlas.upper()
        atlas_ucf = atlas.capitalize()

        if month < 1 or month > 12:
            self._printVerbWrap('Warning: Incorrect month number.')
            return False

        if getattr(self, '{}AtlasLoaded'.format(atlas_ucf)):
            errmsg = 'Warning: {} atlas already initialised.'.format(atlas_uc)
            self._printVerbWrap(errmsg)
            return False

        if not self.CoeffsLoaded:
            self._printVerbWrap('Warning: instrument must be loaded first, ' +
                                '{} atlas not initialised'.format(atlas_uc))
            return False

        atlas_path = getattr(self, '{}AtlasPath'.format(atlas_p))
        self._printDbg('{}AtlasPath={}'.format(atlas_p, atlas_path))
        if not atlas_path:
            self._printVerbWrap(('Warning: {}AtlasPath not set. ' +
                                 '{} atlas not initialised').format(atlas_p,
                                                                    atlas_uc))
            return False

        return True

    def _GenericAtlasSetupMonthFix(self, atlas, month):
        atlas = atlas.upper()
        if month is not None:
            return month
        if not self._profiles:
            errmsg = ("Warning: profiles must be set with valid month, " +
                      "{} atlas will not initialise").format(atlas)
            self._printVerbWrap(errmsg)
            return 0
        return self.Profiles.DateTimes[0, 1]

    def _GenericAtlasSetupFinalise(self, atlas, err):
        atlas_uc = atlas.upper()
        if err:
            errmsg = 'Warning: {} atlas not (re-)initialised'.format(atlas_uc)
            self._printVerbWrap(errmsg)
            return False
        else:
            errmsg = '{} atlas loaded successfully'.format(atlas_uc)
            self._printVerbWrap(errmsg)
            setattr(self, '_{}atlasloaded'.format(atlas), True)
            return True

    def brdfAtlasSetup(self, month=None, single_inst=False, version=-1):
        """Initialise the BRDF atlas (all options available).

        :param int month: month (1-12 => Jan-Dec) for which atlas should be
            initialised, optional
        :param bool single_inst: initialise atlas for a single instrument,
            optional
        :param int version: version number for the BRDF atlas, optional
        """
        month = self._GenericAtlasSetupMonthFix('brdf', month)
        if self._GenericAtlasSetupCheck('brdf', 'Brdf', month):
            tmp_instid = self._instid if single_inst else 0
            err = rtwrap.rttov_brdf_atlas_setup(self.BrdfAtlasPath,
                                                month, version, tmp_instid)
            return self._GenericAtlasSetupFinalise('brdf', err)
        else:
            return False

    def irEmisAtlasSetup(self, month=None, ang_corr=False, single_inst=False,
                         version=-1):
        """Initialise the IR atlas (all options available).

        :param int month: month (1-12 => Jan-Dec) for which atlas should be
            initialised, optional
        :param bool ang_corr: apply zenith angle correction to IR emissivities,
            optional
        :param bool single_inst: initialise atlas for a single instrument,
            optional
        :param int version: version number for IR atlas, optional
        """
        month = self._GenericAtlasSetupMonthFix('ir', month)
        if self._GenericAtlasSetupCheck('ir', 'Emis', month):
            tmp_instid = self._instid if single_inst else 0
            err = rtwrap.rttov_ir_emis_atlas_setup(self.EmisAtlasPath,
                                                   month, version, tmp_instid,
                                                   int(ang_corr))
            return self._GenericAtlasSetupFinalise('ir', err)
        else:
            return False

    def mwEmisAtlasSetup(self, month=None, version=-1):
        """Initialise the MW atlas.

        :param int month: month (1-12 => Jan-Dec) for which atlas should be
            initialised, optional
        :param int version: version number for MW atlas, optional
        """
        month = self._GenericAtlasSetupMonthFix('mw', month)
        if self._GenericAtlasSetupCheck('mw', 'Emis', month):
            err = rtwrap.rttov_mw_emis_atlas_setup(self.EmisAtlasPath,
                                                   month, version, self._instid)
            return self._GenericAtlasSetupFinalise('mw', err)
        else:
            return False

    def deallocBrdfAtlas(self):
        """Deallocate memory for the BRDF atlas."""
        if self.BrdfAtlasLoaded:
            self._printVerbWrap("Deallocating the BRDF Atlas.")
            rtwrap.rttov_brdf_atlas_dealloc()

    def deallocIrAtlas(self):
        """Deallocate memory for the IR atlas."""
        if self.IrAtlasLoaded:
            self._printVerbWrap("Deallocating the IR Atlas.")
            rtwrap.rttov_ir_emis_atlas_dealloc()

    def deallocMwAtlas(self):
        """Deallocate memory for the MW atlas."""
        if self.MwAtlasLoaded:
            self._printVerbWrap("Deallocating the MW Atlas.")
            rtwrap.rttov_mw_emis_atlas_dealloc()

    def _getSurfEmisRefl(self):
        if self._surfemisrefl is not None:
            return self._surfemisrefl
        else:
            raise RttovInternalError("SurfEmisRefl not yet initialised.")

    def _setSurfEmisRefl(self, a):
        if (isinstance(a, np.ndarray) and a.shape[0] == 2 and
                len(a.shape) == 3):
            if a.dtype is not np.dtype(wrapfloat) or np.isfortran(a):
                a = np.array(a, order='C', dtype=wrapfloat)
            self._surfemisrefl_manual = True
            self._surfemisrefl = a
        else:
            raise TypeError("Incorrect value type/shape for SurfEmisRefl")

    def _delSurfEmisRefl(self):
        self._surfemisrefl_manual = False
        self._surfemisrefl = None

    dstring = """Array containing input/output surface emissivity and reflectance values.

    This must be a ``numpy.ndarray`` of shape [2, nprofiles, nchannels], the
    size of the last two dimensions cannot be checked here. However, it will
    eventually be checked during the :meth:`runDirect` or :meth:`runK` calls.

    The `dtype` of the ``numpy.ndarray`` must be {!r}. If the input array has
    another `dtype`, a conversion will be attempted.

    Use `del` to empty the SurfEmisRefl array.

    If this array is not set by the user (or has been emptied unsig `del`), it
    will be initialised to -1 in the :meth:`runDirect` or :meth:`runK` calls.
    """.format(wrapfloat)
    SurfEmisRefl = property(_getSurfEmisRefl, _setSurfEmisRefl,
                            _delSurfEmisRefl, doc=dstring)

    def _getProfiles(self):
        if self._profiles is not None:
            return self._profiles
        else:
            raise RttovInternalError("Profiles not yet initialised.")

    def _setProlifes(self, p):
        if not isinstance(p, Profiles):
            raise TypeError("Incorrect value for Profiles")
        if not p.check():
            errmsg = "Error: some mandatory profile fields are missing"
            raise RttovWrapperError(errmsg)
        self._profiles = p

    dstr = (":class:`.Profiles` object currently associated " +
            "with this Rttov instance." +
            "\n\nThis is a Read/Write attribute of type :class:`.Profiles`.")
    Profiles = property(_getProfiles, _setProlifes, doc=dstr)

    def _pressure_fixer(self):
        '''Returns the pressure array. It may come from the coefficient file'''
        if self.Profiles.DefaultPressureLevels:
            if not self.CoeffsLoaded:
                errmsg = ("Error: instrument not loaded, " +
                          "cannot use coefficient file pressure levels")
                raise RttovWrapperError(errmsg)
            if self.CoeffsNlevels != self.Profiles.Nlevels:
                errmsg = ("Error: number of levels differs to coefficient " +
                          "file, cannot use coefficient file pressure levels")
                raise RttovWrapperError(errmsg)
            self._printVerbWrap("Using pressure levels from coefficient file.")
            p = np.empty((self.Profiles.Nprofiles, self.CoeffsNlevels),
                         dtype=wrapfloat)
            p[...] = self.RefPressures
            return p
        else:
            return self.Profiles.P

    def printSurfEmisRefl(self, pfunction=print):
        """Print the :data:`SurfEmisRefl` content on the standard output."""
        fmt_e = 'Channel #{:5d} e={:6.3f} r={:6.3f}'
        try:
            emisrefl = self.SurfEmisRefl
        except RttovInternalError:
            emisrefl = None
        if emisrefl is not None:
            for prof in range(self.Profiles.Nprofiles):
                pfunction(">> Surfemiss/Refl for profile #{:d}".format(prof))
                for ich in range(emisrefl.shape[-1]):
                    pfunction(fmt_e.format(ich + 1,
                                           emisrefl[0, prof, ich],
                                           emisrefl[1, prof, ich]))
        else:
            pfunction(">> Surfemiss/Refl is not (yet?) initialised.")

    def _prerun_commons(self, channels):
        if channels is None:
            channels = np.arange(1, self.Nchannels + 1, dtype=wrapint)
        if (not isinstance(channels, np.ndarray) or
                channels.dtype is not wrapint):
            channels = np.array(channels, dtype=wrapint)
        if not self.CoeffsLoaded:
            raise RttovWrapperError("Error: coefficients not loaded")
        if self._profiles is None:
            raise RttovWrapperError("Error: profiles are not set yet")

        # update the options which are locally stored in the class instance
        # for RTTOV
        self.updateOptions()

        # deal with SurfEmisRefl
        surfemisrefl_shape = (2, self.Profiles.Nprofiles, len(channels))
        if self._surfemisrefl_manual:
            # Check that the dimensions are fine
            if self.SurfEmisRefl.shape != surfemisrefl_shape:
                raise RttovWrapperError("Incorrect shape for SurfEmisRefl")
        else:
            self._printVerbWrap("No surface emissivity/reflectance supplied:" +
                                " setting calcemis and calcrefl to true")
            self._surfemisrefl = - np.ones(surfemisrefl_shape, dtype=wrapfloat)

        # rads and btrefl are always required
        self._printDbg(("Allocating btrefl and rads for {:d} profiles and " +
                        "{:d} channels").format(self.Profiles.Nprofiles,
                                                len(channels)))
        # First clear _bBasics of the previous results
        self._bBasics = dict()
        self._bBasics['btrefl'] = np.empty((self.Profiles.Nprofiles,
                                            len(channels)), dtype=wrapfloat)
        self._bBasics['rads'] = np.empty((self.Profiles.Nprofiles,
                                          len(channels)), dtype=wrapfloat)

        # A lot of prints...
        if self._debug:
            self.Profiles.printProfiles(pfunction=self._printDbg)
            self.Profiles.printGases(pfunction=self._printDbg)
            self.printSurfEmisRefl(pfunction=self._printDbg)

        return (channels, self._pressure_fixer())

    def runDirect(self, channels=None):
        """Run the RTTOV direct model.

        :param channels: list of channels to simulate (all the channels if
                         omitted)
        :type channels: list of int or numpy.ndarray
        """
        channels, p = self._prerun_commons(channels)

        err = rtwrap.rttov_call_direct(
            self._instid,
            channels,
            self.Profiles.DateTimes.transpose(),    # profile dates/times                                     [nprofiles][6]
            self.Profiles.Angles.transpose(),       # satzen, satazi, sunzen, sunazi angles                   [nprofiles][4]
            self.Profiles.SurfGeom.transpose(),     # lat, lon, elevation                                     [nprofiles][3]
            self.Profiles.SurfType.transpose(),     # surftype, watertype                                     [nprofiles][2]
            self.Profiles.Skin.transpose(),         # skin T, salinity, snow_frac, foam_frac, fastem_coefsx5  [nprofiles][9]
            self.Profiles.S2m.transpose(),          # 2m p, 2m t, 2m q, 10m wind u, v, wind-fetch             [nprofiles][6]
            self.Profiles.SimpleCloud.transpose(),  # ctp, cfraction                                          [nprofiles][2]
            self.Profiles.IceCloud.transpose(),     # ish, idg                                                [nprofiles][2]
            self.Profiles.Zeeman.transpose(),       # Be, cosbk                                               [nprofiles][2]
            p.transpose(),                          # pressure                                                [nprofiles][nlevels]
            self.Profiles.T.transpose(),            # temperature                                             [nprofiles][nlevels]
            self.Profiles.GasUnits,                 # units for gas profiles
            self.Profiles.GasId.transpose(),        # gas ID list                                             [ngases]
            self.Profiles.Gases.transpose(),        # gas profiles                                            [ngases][nprofiles][nlevels]
            self._surfemisrefl.transpose(),         # input/output surface emissivities/BRDFs                 [2][nprofiles][nchannels]
            self._bBasics['btrefl'].transpose(),    # output BTs/refls (for thermal/solar chans)              [nprofiles][nchannels]
            self._bBasics['rads'].transpose(),      # output radiances                                        [nprofiles][nchannels]
            )
        self._errCheck(err, 'Error in rttov_call_direct')

        self._doStoreTrans(len(channels))
        self._doStoreRad(len(channels))
        self._doStoreRad2(len(channels))

    def runK(self, channels=None):
        """Run the RTTOV K model.

        :param channels: list of channels to simulate (all the channels if
                         omitted)
        :type channels: list of int or numpy.ndarray
        """
        channels, p = self._prerun_commons(channels)

        self._printDbg("Allocating the many _k arrays")

        def ini_k_like(tab):
            theshape = list(tab.shape)
            # An extra dimension is needed for K computations
            theshape.insert(-1, len(channels))
            return np.empty(theshape, dtype=tab.dtype)

        self._bBasics['skin_k'] = ini_k_like(self.Profiles.Skin)
        self._bBasics['s2m_k'] = ini_k_like(self.Profiles.S2m)
        self._bBasics['simplecloud_k'] = ini_k_like(self.Profiles.SimpleCloud)
        self._bBasics['p_k'] = ini_k_like(p)
        self._bBasics['t_k'] = ini_k_like(self.Profiles.T)
        self._bBasics['gases_k'] = ini_k_like(self.Profiles.Gases)

        self._bBasics['surfemisrefl_k'] = np.zeros_like(self._surfemisrefl)
        self._bBasics['bt_k'] = np.ones_like(self._bBasics['btrefl'])
        self._bBasics['rad_k'] = np.ones_like(self._bBasics['rads'])

        err = rtwrap.rttov_call_k(
            self._instid,
            channels,
            self.Profiles.DateTimes.transpose(),        # profile dates/times                                     [nprofiles][6]
            self.Profiles.Angles.transpose(),           # satzen, satazi, sunzen, sunazi angles                   [nprofiles][4]
            self.Profiles.SurfGeom.transpose(),         # lat, lon, elevation                                     [nprofiles][3]
            self.Profiles.SurfType.transpose(),         # surftype, watertype                                     [nprofiles][2]
            self.Profiles.Skin.transpose(),             # skin T, salinity, snow_frac, foam_frac, fastem_coefsx5  [nprofiles][9]
            self._bBasics['skin_k'].transpose(),        # skin T, salinity, snow_frac, foam_frac, fastem_coefsx5  [nprofiles][nchannels][9]
            self.Profiles.S2m.transpose(),              # 2m p, 2m t, 2m q, 10m wind u, v, wind-fetch             [nprofiles][6]
            self._bBasics['s2m_k'].transpose(),         # 2m p, 2m t, 2m q, 10m wind u, v, wind-fetch             [nprofiles][nchannels][6]
            self.Profiles.SimpleCloud.transpose(),      # ctp, cfraction                                          [nprofiles][2]
            self._bBasics['simplecloud_k'].transpose(), # ctp, cfraction                                          [nprofiles][nchannels][2]
            self.Profiles.IceCloud.transpose(),         # ish, idg                                                [nprofiles][2]
            self.Profiles.Zeeman.transpose(),           # Be, cosbk                                               [nprofiles][2]
            p.transpose(),                              # pressure                                                [nprofiles][nlevels]
            self._bBasics['p_k'].transpose(),           # pressure                                                [nprofiles][nchannels][nlevels]
            self.Profiles.T.transpose(),                # temperature                                             [nprofiles][nlevels]
            self._bBasics['t_k'].transpose(),           # temperature                                             [nprofiles][nchannels][nlevels]
            self.Profiles.GasUnits,                     # units for gas profiles
            self.Profiles.GasId.transpose(),            # gas ID list                                             [ngases]
            self.Profiles.Gases.transpose(),            # gas profiles                                            [ngases][nprofiles][nlevels]
            self._bBasics['gases_k'].transpose(),       # gas profiles                                            [ngases][nprofiles][nchannels][nlevels]
            self._surfemisrefl.transpose(),             # input/output surface emissivities/BRDFs                 [2][nprofiles][nchannels]
            self._bBasics['surfemisrefl_k'].transpose(),# input/output surface emissivities/BRDFs                 [2][nprofiles][nchannels]
            self._bBasics['btrefl'].transpose(),        # output BTs/refls (for thermal/solar chans)              [nprofiles][nchannels]
            self._bBasics['rads'].transpose(),          # output radiances                                        [nprofiles][nchannels]
            self._bBasics['bt_k'].transpose(),          # output BTs/refls (for thermal/solar chans)              [nprofiles][nchannels]
            self._bBasics['rad_k'].transpose(),         # output radiances                                        [nprofiles][nchannels]
            )
        self._errCheck(err, 'Error in rttov_call_k')

        # Save the gas_id list for a later use
        self._bBasics['gases_k_id'] = self.Profiles.GasId

        self._doStoreTrans(len(channels))
        self._doStoreRad(len(channels))
        self._doStoreRad2(len(channels))

    def _1DstoreGeneric(self, outd, storeid, callback, nchannels):
        outd[storeid] = np.empty((self.Profiles.Nprofiles * nchannels),
                                 dtype=wrapfloat)
        err = callback(self._instid, outd[storeid])
        self._errCheck(err, 'Error when getting {}'.format(storeid))
        outd[storeid] = outd[storeid].reshape((self.Profiles.Nprofiles,
                                               nchannels))

    def _2DstoreGeneric(self, outd, storeid, callback, nchannels,
                        onlayers=False):
        trailind_dim = self.Profiles.Nlevels - int(onlayers)
        outd[storeid] = np.empty((self.Profiles.Nprofiles * nchannels,
                                  trailind_dim), dtype=wrapfloat)
        err = callback(self._instid, outd[storeid].transpose())
        self._errCheck(err, 'Error when getting {}'.format(storeid))
        outd[storeid] = outd[storeid].reshape((self.Profiles.Nprofiles,
                                               nchannels,
                                               trailind_dim))

    def _doStoreTrans(self, nchannels):
        if not self.Options.StoreTrans:
            self._bTrans = dict()
            return
        self._1DstoreGeneric(self._bTrans, 'tau_total',
                             rtwrap.rttov_get_tau_total, nchannels)
        self._2DstoreGeneric(self._bTrans, 'tau_levels',
                             rtwrap.rttov_get_tau_levels, nchannels)
        self._1DstoreGeneric(self._bTrans, 'tausun_total_path1',
                             rtwrap.rttov_get_tausun_total_path1, nchannels)
        self._2DstoreGeneric(self._bTrans, 'tausun_levels_path1',
                             rtwrap.rttov_get_tausun_levels_path1, nchannels)
        self._1DstoreGeneric(self._bTrans, 'tausun_total_path2',
                             rtwrap.rttov_get_tausun_total_path2, nchannels)
        self._2DstoreGeneric(self._bTrans, 'tausun_levels_path2',
                             rtwrap.rttov_get_tausun_levels_path2, nchannels)

    def _doStoreRad(self, nchannels):
        if not self.Options.StoreRad:
            self._bRad = dict()
            return
        self._1DstoreGeneric(self._bRad, 'rad_clear',
                             rtwrap.rttov_get_rad_clear, nchannels)
        self._1DstoreGeneric(self._bRad, 'rad_total',
                             rtwrap.rttov_get_rad_total, nchannels)
        self._1DstoreGeneric(self._bRad, 'bt_clear',
                             rtwrap.rttov_get_bt_clear, nchannels)
        self._1DstoreGeneric(self._bRad, 'bt',
                             rtwrap.rttov_get_bt, nchannels)
        self._1DstoreGeneric(self._bRad, 'refl_clear',
                             rtwrap.rttov_get_refl_clear, nchannels)
        self._1DstoreGeneric(self._bRad, 'refl',
                             rtwrap.rttov_get_refl, nchannels)
        self._1DstoreGeneric(self._bRad, 'rad_cloudy',
                             rtwrap.rttov_get_rad_cloudy, nchannels)
        self._2DstoreGeneric(self._bRad, 'overcast',
                             rtwrap.rttov_get_overcast, nchannels,
                             onlayers=True)

    def _doStoreRad2(self, nchannels):
        if not self.Options.StoreRad2:
            self._bRad2 = dict()
            return
        self._1DstoreGeneric(self._bRad2, 'rad2_up_clear',
                             rtwrap.rttov_get_rad2_upclear, nchannels)
        self._1DstoreGeneric(self._bRad2, 'rad2_dn_clear',
                             rtwrap.rttov_get_rad2_dnclear, nchannels)
        self._1DstoreGeneric(self._bRad2, 'rad2_refldn_clear',
                             rtwrap.rttov_get_rad2_refldnclear, nchannels)
        self._2DstoreGeneric(self._bRad2, 'rad2_up',
                             rtwrap.rttov_get_rad2_up, nchannels,
                             onlayers=True)
        self._2DstoreGeneric(self._bRad2, 'rad2_down',
                             rtwrap.rttov_get_rad2_down, nchannels,
                             onlayers=True)
        self._2DstoreGeneric(self._bRad2, 'rad2_surf',
                             rtwrap.rttov_get_rad2_surf, nchannels,
                             onlayers=True)

    dstring = ("Computed brightness temperatures or reflectances from the " +
               "pevious run. ``numpy.ndarray`` of shape " +
               "[nprofiles, nchannels].")
    BtRefl = _DescRO('bBasics', 'btrefl', doc=dstring)

    dstring = ("Computed radiances from the pevious run. " +
               "``numpy.ndarray`` of shape [nprofiles, nchannels].")
    Rads = _DescRO('bBasics', 'rads', doc=dstring)

    dstring = ("Computed pressure jacobians from the pevious run. " +
               "``numpy.ndarray`` of shape [nprofiles, nchannels, nlevels].")
    PK = _DescRO('bBasics', 'p_k', doc=dstring)

    dstring = ("Computed temperature jacobians from the pevious run. " +
               "``numpy.ndarray`` of shape [nprofiles, nchannels, nlevels].")
    TK = _DescRO('bBasics', 't_k', doc=dstring)

    dstring = ("Computed skin variables jacobians from the pevious run. " +
               "``numpy.ndarray`` of shape [nprofiles, nchannels, 9].")
    SkinK = _DescRO('bBasics', 'skin_k', doc=dstring)

    dstring = ("Computed 2m variables jacobians from the pevious run. " +
               "``numpy.ndarray`` of shape [nprofiles, nchannels, 6].")
    S2mK = _DescRO('bBasics', 's2m_k', doc=dstring)

    dstring = ("Computed simple cloud jacobians from the pevious run. " +
               "``numpy.ndarray`` of shape [nprofiles, nchannels, 2].")
    SimpleCloudK = _DescRO('bBasics', 'simplecloud_k', doc=dstring)

    dstring = ("Computed gas, cloud and aerosol jacobians from the pevious " +
               "run. ``numpy.ndarray`` of shape " +
               "[ngases, nprofiles, nchannels, nlevels].")
    GasesK = _DescRO('bBasics', 'gases_k', doc=dstring)

    def getItemK(self, gas_id):
        """Computed jacobian for a given gas, cloud or aerosol.

        If the jacobian was not computed, ``None`` will be returned.

        :param gas_id: the gas_id to look for
        :type gas_id: int, str or itemIdType
        :return: computed jacobian
        :rtype: numpy.ndarray of shape [nprofiles, nchannels, nlevels].
        """
        g_ids = self._bBasics.get('gases_k_id', None)
        g = self._bBasics.get('gases_k', None)
        if g_ids is not None and g is not None:
            g_idx = np.where(g_ids == itemIdType(gas_id))
            if len(g_idx[0]):
                return g[g_idx[0][0], ...]
        # The gas was not found... (or no computation was done
        return None

    @property_ro
    def SurfEmisK(self):
        """Computed surface emissivity jacobians from the pevious run. ``numpy.ndarray`` of shape [nprofiles, nchannels]."""
        res = self._bBasics.get('surfemisrefl_k', None)
        if res is not None:
            res = res[0, ...]
        return res

    @property_ro
    def SurfReflK(self):
        """Computed surface reflectivity jacobians from the pevious run. ``numpy.ndarray`` of shape [nprofiles, nchannels]."""
        res = self._bBasics.get('surfemisrefl_k', None)
        if res is not None:
            res = res[1, ...]
        return res

    # Outputs under the store_trans key
    fdstring = ('RTTOV transmission {} output ``numpy.ndarray`` of shape' +
                '[nprofiles, nchannels{}]. Requires store_trans=True.')

    dstring = fdstring.format('tau_total', '')
    TauTotal = _DescRO('bTrans', 'tau_total', doc=dstring)

    dstring = fdstring.format('tau_levels', ', nlevels')
    TauLevels = _DescRO('bTrans', 'tau_levels', doc=dstring)

    dstring = fdstring.format('tausun_total_path1', '')
    TauSunTotalPath1 = _DescRO('bTrans', 'tausun_total_path1', doc=dstring)

    dstring = fdstring.format('tausun_levels_path1', ', levels')
    TauSunLevelsPath1 = _DescRO('bTrans', 'tausun_levels_path1', doc=dstring)

    dstring = fdstring.format('tausun_total_path1', '')
    TauSunTotalPath2 = _DescRO('bTrans', 'tausun_total_path2', doc=dstring)

    dstring = fdstring.format('tausun_levels_path2', ', levels')
    TauSunLevelsPath2 = _DescRO('bTrans', 'tausun_levels_path2', doc=dstring)

    # Outputs under the store_rad key
    fdstring = ('RTTOV radiance {} output ``numpy.ndarray`` of shape' +
                '[nprofiles, nchannels{}]. Requires store_rad=True.')

    dstring = fdstring.format('clear', '')
    RadClear = _DescRO('bRad', 'rad_clear', doc=dstring)

    dstring = fdstring.format('total', '')
    RadTotal = _DescRO('bRad', 'rad_total', doc=dstring)

    dstring = fdstring.format('bt_clear', '')
    BtClear = _DescRO('bRad', 'bt_clear', doc=dstring)

    dstring = fdstring.format('bt', '')
    Bt = _DescRO('bRad', 'bt', doc=dstring)

    dstring = fdstring.format('refl_clear', '')
    ReflClear = _DescRO('bRad', 'refl_clear', doc=dstring)

    dstring = fdstring.format('refl', '')
    Refl = _DescRO('bRad', 'refl', doc=dstring)

    dstring = fdstring.format('cloudy', '')
    RadCloudy = _DescRO('bRad', 'rad_cloudy', doc=dstring)

    dstring = fdstring.format('overcast', '')
    Overcast = _DescRO('bRad', 'overcast', doc=dstring)

    # Outputs under the store_rad2 key
    fdstring = ('RTTOV radiance2 {} output ``numpy.ndarray`` of shape' +
                '[nprofiles, nchannels{}]. Requires store_rad2=True.')

    dstring = fdstring.format('upclear', '')
    Rad2UpClear = _DescRO('bRad2', 'rad2_up_clear', doc=dstring)

    dstring = fdstring.format('dnclear', '')
    Rad2DnClear = _DescRO('bRad2', 'rad2_dn_clear', doc=dstring)

    dstring = fdstring.format('refldnclear', '')
    Rad2ReflDnClear = _DescRO('bRad2', 'rad2_refldn_clear', doc=dstring)

    dstring = fdstring.format('up', ', nlayers')
    Rad2Up = _DescRO('bRad2', 'rad2_up', doc=dstring)

    dstring = fdstring.format('down', ', nlayers')
    Rad2Down = _DescRO('bRad2', 'rad2_down', doc=dstring)

    dstring = fdstring.format('surf', ', nlayers')
    Rad2Surf = _DescRO('bRad2', 'rad2_surf', doc=dstring)
