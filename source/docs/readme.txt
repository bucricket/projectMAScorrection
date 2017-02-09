##############################################################################
#                                                                            #
#                               RTTOV v11.3                                  #
#         Release Notes for version 11.3 FCM version number 1835             #
#                          28 September 2015                                 #
#                                                                            #
##############################################################################

1. Licence
To use this software, users need to have registered for RTTOV v11 with the
NWP SAF (http://nwpsaf.eu/), and to have agreed to the terms of the RTTOV v11
license agreement.


2. Updates from RTTOV v11.2


New features/improvements:

- new Python and C++ interface to RTTOV: allows much RTTOV functionality
  to be called directly from Python or C++.

- updates to the RTTOV GUI including new 1DVar functionality.

- new option to specify input units for gas profiles (ppmv over moist air,
  ppmv over dry air, kg/kg over moist air).

- PC-RTTOV can now be run for all surface types (requires new PC
  coefficient files available on the website).

- option to treat surfaces as Lambertian reflectors for the downwelling
  surface-reflected emission can now be applied to IR as well as MW sensors.

- option to supply ocean surface foam coverage fraction to FASTEM instead
  of using the internal foam parameterisation.

- improved treatment of snow in land surface BRDF atlas (requires a new
  eigenvector data file to be downloaded from the website).

- IR emissivity atlas has new option to include a correction for zenith
  angle (requires new additional atlas data files available on the
  website).

- IR emissivity and BRDF atlases now also available in HDF5 format with
  much reduced file size due to compression, and avoids the need for
  compiling against the netCDF library.

- example programs simplified and new examples added for IR cloudy
  simulations and RTTOV-SCATT.

- four new subroutines which can be used to allocate any or all
  input/output arrays and structures for the RTTOV direct, TL, AD
  and K models simplifying code which calls RTTOV.

- new interactive script to make compiling RTTOV easier.

- new interactive script to download coefficients from website.

- RTTOV libraries no longer have minor version numbers in their names.



Bug fixes (see the RTTOV v11 web page for more details):

- the treatment of gas units was not consistent between IR and MW sensors
  (see http://nwpsaf.eu/deliverables/rtm/docs_rttov11/rttov_v11.3_gas_units
  _and_new_coefs.pdf).

- fix to prevent fatal errors occurring under certain circumstances when
  reg_limit_extrap and do_checkinput options are both true.

- fix run-time error in TL which occurs when cloudy IR simulations are
  run with all cfrac values set to zero.

- fix to make treatment of optical depths at the top of the atmosphere for
  interpolation modes 4 and 5 consistent with modes 1-3 when spacetop option
  is true (only has a significant impact if input pressure profile does not
  span weighting functions of channels being simulated).

- when apply_reg_limits option is true, ensure this is applied to the level
  just below the surface as this can contribute to the optical depth
  calculation.

- BRDF atlas was artificially capping BRDF values at 1/pi.

- errors in emissivity_ad/k for channels with calcemis set to false when
  there are also some channels in the same call with calcemis set to true.



Updates from RTTOV v11.1

New features/improvements:

- RTTOV graphical user interface (GUI) is included in the package: allows you
  to run direct and K models and visualise results immediately.

- new interpolation options are available which eliminate the jaggedness
  observed in temperature Jacobians under certain circumstances.

- new option for profile extrapolation at top of atmosphere based on
  regression limits.

- new FASTEM-6 option which has improved treatment of azimuthal dependence.

- updated Baran ice optical property parameterisation for cloudy IR
  simulations.

- apply_reg_limits is now independent of do_checkinput so the former can be
  applied when do_checkinput is false.

- HDF5 coef I/O now allows for channel selection when reading/writing files.

- new rttov_coef_info.exe executable which can be used to print out
  information about a coefficient file (of any format) on commandline.

- new example_k.F90 program which gives an example of calling the K model.

- example program src/other/rttov_obs_to_pc.F90 which demonstrates conversion
  of observations to PC space for PC assimilation applications.

- optimisation of PC code.

- optimisation of sea-surface refl code for solar calculations (especially for
  K model).

- optimisation of IR cloud and aerosol scattering code.

- test suite profile sets are now based on the six US76 standard atmosphere
  profiles.


Bug fixes:

- when extracting a subset of channels from a coef file using
  rttov_conv_coef.exe, NLTE coefficients are now extracted only for the
  requested channels. NB If you have extracted a subset of channels from any
  hi-res sounder coef file containing NLTE coefs to a new coef file, this must
  be recreated using rttov_conv_coef.exe.

- running an NLTE simulation for IASI channel 7021, AIRS channel 2122 and CrIS
  channel 1245 (i.e. the channels immediately following the block of
  NLTE-affected channels for each instrument) causes an array bounds error.

- for NLTE simulations fix minor bug whereby rad%cloudy, rad_tl%cloudy,
  rad_tl%total, rad_tl%bt and profiles_ad/k%cfraction are in error with the
  simple cloud scheme (cfraction > 0).

- add missing code to calculate profiles_k_rec%icede for PC cloudy simulations.

- array bounds errors can occur under certain conditions when reading PC
  coefficients.

- the TL radiances and surface emissivities and the AD/K 10m wind u and v
  components and wind fetch may be in error for PC calculations.

- the TL radiances and surface reflectance and the AD/K 10m wind u and v
  components and wind fetch may be in error for solar simulations.

- division by zero is possible with solar simulations when 10m u wind speed is
  exactly zero and 10m v wind speed is negative.

- solar simulations may be in error when the profile month is set to January.

- rarely-occurring bug which can result in small errors in simulated
  radiances, reflectances and BTs (and also in TL/AD/K output) for
  solar-affected channels (requires multiple profiles with different dates to
  be passed to parallel interface).

- BRDF atlas should output negative values for channels with wavelengths
  beyond the range of the dataset instead of zeros.

- bug in FASTEM-3 K code related to the azimuthal correction (errors are
  small in magnitude).

- ice effective diameter limits are not applied in very rare cases when
  calculated effective diameter is negative.

- bug in rttov_alloc_aux_prof whereby code would fail due to an array
  bounds error if run with 12 or 13 predefined aerosol types present.

- bug rttov_checkinput.F90 subroutine which means that RTTOV fails to
  throw a fatal error in the case where an element of an input temperature,
  water vapour or trace gas profiles exceeds the hard profile limits.

- checkinput and apply_reg_limits subroutines now check all levels down to
  the one below the surface (previously they stopped at the level above the
  surface).

- when calling rttov_read_coefs to read HDF5-formatted coefficients the code
  may fail because the HDF5 library has not been initialised.



Updates from RTTOV v10.2

New features:

- capability to simulate visible/near-IR channels.

- land surface bi-directional reflectance factor (BRDF) atlas.

- two new aerosol particle types for volcanic ash and dust.

- new ice cloud parameterisation scheme for cloudy IR simulations.

- capability to specify aerosol and cloud optical parameters directly instead
  of supplying number density profiles for pre-defined particle types for
  visible/IR scattering simulations.

- capability to include a correction for non-local thermodynamic equilibrium
  effects for hi-res sounders.

- option to treat surfaces as Lambertian reflectors when calculating surface-
  reflected down-welling component of atmospheric emission.

- Principle Components calculations may be carried out for clear and cloudy
  profiles over sea surfaces, and for IASI it is possible to calculate PCs
  for limited spectral bands.

- additional options for internal interpolation.

- updates to the IR land surface emissivity atlas to improve speed and memory
  usage.

- optional capability to read/write coefficients in HDF5 format.

- capability to convert MW scattering coefficient files (mietables) to binary
  format for improved performance.

- updates to SSU coefficients to account for time-variation of cell pressure.


3. Installation instructions

This readme file refers to the version v11.3 of RTTOV v11. This is the third
release of the RTTOV v11 code. The entire v11.3 package is available via the
NWP SAF helpdesk. This file is called rttov113.tar.gz. New RTTOV v11 users MUST
obtain this file by requesting it via the helpdesk.

Alternatively, if you already have RTTOV v11.2, the updates for RTTOV v11.3
are contained in a single file rttov113_update.tar.gz. This file is available
from the RTTOV v11 web page:
http://nwpsaf.eu/deliverables/rtm/rtm_rttov11.html

If you have RTTOV v11.1 you can upgrade to v11.3 by applying the v11.2 and
v11.3 upgrades: see the documentation associated with the
rttov112_update.tar.gz file on the RTTOV web page.

Only updated source code, documentatation, compilation flags, test scripts,
test input datasets and test reference outputs are provided in the update, no
coefficient files. Updates to coefficient files are available on the RTTOV v11
web page. Instructions on how to apply the update may be found below.

NB For the release of v11.3 most coefficient files were regenerated so users
   are advised to read the associated documentation and download the latest
   coefficients from the web site.

Coefficient files are available on the RTTOV v11 web page at:
http://nwpsaf.eu/deliverables/rtm/rtm_rttov11.html

The docs/ directory contains the following files:
- a copy of this file (readme.txt)
- user guide (users_guide_11_v1.4.pdf - NWPSAF-MO-UD-028)
- test suite guide (rttov-test.pdf - NWPSAF-MO-TV-031)
- user guide for the GUI (rttov_gui_for_rttov_v11.3.pdf - NWPSAF-MF-UD-010)
- user guide for the wrapper (rttov-wrapper.pdf - NWPSAF-MO-UD-035)

New users are advised to read the user guide which gives all the details
necessary to compile and run the code. Updated versions of this may also be
obtained from the above link.

RTTOV v11 is designed for UNIX/Linux systems. The software is now successfully
tested on the following architectures and compilers: IBM; Cray; Intel systems
with gfortran, ifort, NAG, and pgf90; MacOS-X with gfortran. The compilers tested
are represented by the compiler flag files available in the build/arch/
directory.

The following system components are needed before running RTTOV v11:
    * UNIX or Linux operating system
    * Fortran 90 compiler
    * make utilities
    * gzip and gunzip
    * about 100 MBytes of free disk space (several GBytes if you require
      hyperspectral IR sounder coefficient files)

RTTOV v11 will not work with older versions of some compilers. The following list
gives the versions of several common compilers known to be compatible:
  gfortran - v4.4.7 and later
  ifort - v12.0.4 and later
  NAG - v5.2, v5.3
  pgf90 - v11.7 and later
  IBM - xlf95 v12.1 and v13.1
  Cray Fortran - v8.3.4

Known compiler issues:
  - NAG v6.0: cannot read nested derived types in namelists and so cannot
    run the test suite. The test and example programs outside of the test
    harness (rttov_test.pl) run fine which suggests the core of RTTOV is
    compiled successfully and can be used with this compiler.

Some basic information on installing the RTTOV v11 Fortran 90 code in a UNIX or
LINUX environment follows.

------------------------------------------------------------------------------

For new users:

Instructions for installing the FULL RTTOV v11.3 package obtained via the
NWP SAF helpdesk.

The file name should be rttov113.tar.gz and this file should be copied to
the location in which you wish to install RTTOV (e.g. ~user/rttov113/).

Extract the tarball in this top-level RTTOV v11.3 directory:
$ tar -zxvf rttov113.tar.gz

This creates a number of folders in the top-level directory.
Larger coefficient files such as those for hi-res IR sounders are not
included in the package so if you wish to simulate these instruments you
must download the relevant coefficient files from the web site.

The rtcoef_rttov11/ directory contains various sub-directories for
different types of coefficient file. Coefficient files should be
downloaded from the RTTOV v11 web page and placed in the appropriate
locations in order to successfully run all tests.

The script rtcoef_rttov11/rttov_coef_download.sh can be run to download the
coefficients to their default locations within rtcoef_rttov11/.

The web page for RTTOV v11 coefficients is:
http://nwpsaf.eu/deliverables/rtm/rttov11_coefficients.html

------------------------------------------------------------------------------

For existing RTTOV v11.1 and v11.2 users:

Instructions for installing the UPDATE package obtained from the RTTOV v11 web
page to upgrade from RTTOV v11.2 to v11.3.

The file name should be rttov113_update.tar.gz

Navigate to your top-level RTTOV v11.2 directory. For example:
$ cd ~user/rttov11

Remove some files which are no longer needed:
$ rm -f src/other/us76_to_hdf5.F90
$ rm -f src/test/rttov_zutility_test.F90
$ rm -f src/test/test_weighting_fn_dev.F90
$ rm -fr rttov_test/*

Copy the file rttov113_update.tar.gz into the top-level RTTOV directory and
extract the updates:
$ tar zxf rttov113_update.tar.gz

Note: to update from v11.1 first you must download the rttov112_update.tar.gz
package from the website and follow the instructions contained in the
readme.txt file. You CANNOT update from v11.1 to v11.3 directly using the
rttov113_update.tar.gz package.

As noted above no coefficient files are included in the update tarball. You
should download the latest coefficient files from the RTTOV v11 web page and
place them in the rtcoef_rttov11/ directory. The updated test reference data
sets were generated with the latest coefficients.

The script rtcoef_rttov11/rttov_coef_download.sh can be run to download the
coefficients to their default locations within rtcoef_rttov11/.

------------------------------------------------------------------------------

The fortran code is organised within a number of subdirectories within
the src/ directory. The code consists of subroutine files and top
level programs for testing the RTTOV v11 code:
- test/rttov_test.F90              Test program for RTTOV (to be called via
                                   rttov_test/rttov_test.pl)
- test/example_*.F90               Example source code provided to aid user
- brdf_atlas/example_atlas_fwd.F90
- mw_scatt/rttovscatt_test.F90     Test program for scattering code
- coef_io/rttov_conv_coef.F90      Program to convert ascii coeff files to
                                   binary or HDF5 format.

The easiest way to compile the code is to run the following interactive
script:

$ cd src/
$ ../build/rttov_compile.sh

This script prompts you for information and describes which features of
RTTOV are available given the input you provide.

The directory build/arch/ contains a list of files containing compilation
flags for a variety of architectures and compilers. You should give the
name of one of these files when prompted. You can create your own using
the existing files as a template. If you want to take advantage of
multi-threaded execution via the RTTOV parallel interface you must choose
one of the "openmp" compiler flag files.

It is recommended to compile RTTOV with the HDF5 library: to do this you
must edit the file ../build/Makefile.local so that the HDF5 lines are
uncommented and the HDF5_PREFIX variable points to your installation of
the HDF5 library (v1.8.8 or later). You should do this *before* running
the rttov_compile.sh script.

If you compile RTTOV with the HDF5 library and you wish to use the land
surface IR emissivity and BRDF atlases you MUST download the new HDF5
format files from the website.

The IR emissivity and BRDF atlas files are still available in netCDF
format: if you do not compile RTTOV with HDF5 you can edit the file
../build/Makefile.local with details of your netCDF installation and
you can then use the netCDF format atlas files.

If you want to run the RTTOV GUI or you want to call RTTOV from Python
you must have f2py installed on your system.

It is possible to compile RTTOV manually as with previous versions: see
the user guide for more details.


Once the code is compiled you will find bin/ and lib/ directories in your
top-level RTTOV directory containing the RTTOV binaries and libraries. One
library is created for each subfolder within src/ and you should link all
required libraries in your application (at the very least librttov11_main.a
and librttov11_coef_io - see the user guide for more information).

Example stand-alone Makefiles are included for the example_*fwd.F90
demonstration programs. These are src/test/Makefile_examples and
src/brdf_atlas/Makefile_examples. Each has a section at the top which must be
edited with paths appropriate the the user's system. These demonstrate how to
link the RTTOV libraries with separate application code.


The following executables are created in the bin directory (found in the
top-level directory by default):
bin/rttov_test.exe            (basic RTTOV v11 test program)
bin/example_*.exe             (compilation of example source code provided)
bin/rttovscatt_test.exe       (test program for RTTOV-SCATT)
bin/rttov_*atlas_test.exe     (test programs for emissivity atlas code)
bin/rttov_brdf_atlas_test.exe (test program for BRDF atlas code)
and some other utilities for reading coefficient files etc

The rttov_test/ directory contains input data and scripts for running the
various test executables listed above. The script rttov_test.pl is used to run
rttov_test.exe and this is fully documented in the docs/ directory. The
test_rttov11.sh script can be used to run a series of tests to verify your
installation. Separate tests exist for the MW scattering and emissivity and BRDF
atlas code, and these are documented in the user guide.

Test reference data was generated using the latest coefficient files from
the website.

If the test scripts complete successfully or if any observed differences 
with the reference data are as expected (see user guide), the code is ready 
to be incorporated in your code.
Report all bugs to the NWP SAF Helpdesk: http://nwpsaf.eu/feedback.html


The contents of rttov113.tar.gz for RTTOV v11.3 are:

readme.txt
docs/
docs/NWPSAFLogo_gradient_S.png
docs/doxygen_config_dev
docs/doxygen_config_user
docs/readme.txt
docs/rttov-test.pdf
docs/rttov-wrapper.pdf
docs/rttov11_svr.pdf
docs/rttov_doxygen_readme.dox
docs/rttov_gui_for_rttov_v11.3.pdf
docs/users_guide_11_v1.4.pdf
src/
src/Makefile
src/brdf_atlas/
src/brdf_atlas/Makefile
src/brdf_atlas/Makefile_examples
src/brdf_atlas/example_atlas_fwd.F90
src/brdf_atlas/mod_brdf_atlas.F90
src/brdf_atlas/mod_rttov_brdf_atlas.F90
src/brdf_atlas/rttov_brdf_atlas_test.F90
src/brdf_atlas/rttov_deallocate_brdf_atlas.F90
src/brdf_atlas/rttov_get_brdf.F90
src/brdf_atlas/rttov_setup_brdf_atlas.F90
src/coef_io/
src/coef_io/Makefile
src/coef_io/rttov_channel_extract_coef.F90
src/coef_io/rttov_channel_extract_pccoef.F90
src/coef_io/rttov_channel_extract_scaercoef.F90
src/coef_io/rttov_channel_extract_sccldcoef.F90
src/coef_io/rttov_check_channels_pc.F90
src/coef_io/rttov_cmpuc.F90
src/coef_io/rttov_coef_io_mod.F90
src/coef_io/rttov_coeffname.F90
src/coef_io/rttov_conv_coef.F90
src/coef_io/rttov_dealloc_coef.F90
src/coef_io/rttov_dealloc_coef_pccomp.F90
src/coef_io/rttov_dealloc_coef_scatt_ir.F90
src/coef_io/rttov_dealloc_coefs.F90
src/coef_io/rttov_dealloc_optpar_ir.F90
src/coef_io/rttov_deletecomment.F90
src/coef_io/rttov_findnextsection.F90
src/coef_io/rttov_get_pc_predictindex.F90
src/coef_io/rttov_init_coef.F90
src/coef_io/rttov_init_coef_optpar_ir.F90
src/coef_io/rttov_init_coef_pccomp.F90
src/coef_io/rttov_init_coefs.F90
src/coef_io/rttov_nullify_coef.F90
src/coef_io/rttov_nullify_coef_pccomp.F90
src/coef_io/rttov_nullify_coef_scatt_ir.F90
src/coef_io/rttov_nullify_coefs.F90
src/coef_io/rttov_nullify_optpar_ir.F90
src/coef_io/rttov_opencoeff.F90
src/coef_io/rttov_q2v.F90
src/coef_io/rttov_read_ascii_coef.F90
src/coef_io/rttov_read_ascii_pccoef.F90
src/coef_io/rttov_read_ascii_scaercoef.F90
src/coef_io/rttov_read_ascii_sccldcoef.F90
src/coef_io/rttov_read_binary_coef.F90
src/coef_io/rttov_read_binary_pccoef.F90
src/coef_io/rttov_read_binary_scaercoef.F90
src/coef_io/rttov_read_binary_sccldcoef.F90
src/coef_io/rttov_read_coefs.F90
src/coef_io/rttov_skipcommentline.F90
src/coef_io/rttov_test_get_pc_predictindex.F90
src/coef_io/rttov_v2q.F90
src/coef_io/rttov_write_ascii_coef.F90
src/coef_io/rttov_write_ascii_pccoef.F90
src/coef_io/rttov_write_ascii_scaercoef.F90
src/coef_io/rttov_write_ascii_sccldcoef.F90
src/coef_io/rttov_write_binary_coef.F90
src/coef_io/rttov_write_binary_pccoef.F90
src/coef_io/rttov_write_binary_scaercoef.F90
src/coef_io/rttov_write_binary_sccldcoef.F90
src/coef_io/rttov_write_coefs.F90
src/coef_io_789/
src/coef_io_789/Makefile
src/coef_io_789/rttov789_coeffname.F90
src/coef_io_789/rttov789_conv_coef.F90
src/coef_io_789/rttov789_conv_coef_11to9.F90
src/coef_io_789/rttov789_opencoeff.F90
src/coef_io_789/rttov789_readcoeffs.F90
src/coef_io_789/rttov789_readcoeffs_ascii.F90
src/coef_io_789/rttov789_write_ascii_coef.F90
src/emis_atlas/
src/emis_atlas/Makefile
src/emis_atlas/mod_cnrm_mw_atlas.F90
src/emis_atlas/mod_iratlas.F90
src/emis_atlas/mod_mwatlas.F90
src/emis_atlas/mod_rttov_emis_atlas.F90
src/emis_atlas/rttov_cnrmmwatlas_test.F90
src/emis_atlas/rttov_deallocate_emis_atlas.F90
src/emis_atlas/rttov_get_emis.F90
src/emis_atlas/rttov_iratlas_test.F90
src/emis_atlas/rttov_mwatlas_test.F90
src/emis_atlas/rttov_setup_emis_atlas.F90
src/gui/
src/gui/Makefile
src/gui/f2py_f2cmap
src/gui/rttov_gui_context.F90
src/gui/rttov_gui_f2py.F90
src/gui/rttov_gui_handle.F90
src/gui/rttov_gui_test_run.F90
src/hdf/
src/hdf/Makefile
src/hdf/rttov_hdf_chanprof_io.F90
src/hdf/rttov_hdf_coefs.F90
src/hdf/rttov_hdf_emissivity_io.F90
src/hdf/rttov_hdf_load.F90
src/hdf/rttov_hdf_mod.F90
src/hdf/rttov_hdf_opt_param_io.F90
src/hdf/rttov_hdf_options_config_io.F90
src/hdf/rttov_hdf_options_interp_io.F90
src/hdf/rttov_hdf_options_io.F90
src/hdf/rttov_hdf_options_pc_io.F90
src/hdf/rttov_hdf_options_rt_all_io.F90
src/hdf/rttov_hdf_options_rt_ir_io.F90
src/hdf/rttov_hdf_options_rt_mw_io.F90
src/hdf/rttov_hdf_pccomp_io.F90
src/hdf/rttov_hdf_profile_io.F90
src/hdf/rttov_hdf_profiles.F90
src/hdf/rttov_hdf_radiance2_io.F90
src/hdf/rttov_hdf_radiance_io.F90
src/hdf/rttov_hdf_reflectance_io.F90
src/hdf/rttov_hdf_rttov_coef_io.F90
src/hdf/rttov_hdf_rttov_coef_pcc1_io.F90
src/hdf/rttov_hdf_rttov_coef_pcc2_io.F90
src/hdf/rttov_hdf_rttov_coef_pcc_io.F90
src/hdf/rttov_hdf_rttov_fast_coef_io.F90
src/hdf/rttov_hdf_rttov_nlte_coef_io.F90
src/hdf/rttov_hdf_s2m_io.F90
src/hdf/rttov_hdf_save.F90
src/hdf/rttov_hdf_sskin_io.F90
src/hdf/rttov_hdf_transmission_io.F90
src/main/
src/main/Makefile
src/main/lapack.f
src/main/mod_rttov_baran2013_icldata.F90
src/main/mod_rttov_baran2014_icldata.F90
src/main/mod_rttov_fastem3_coef.F90
src/main/mod_rttov_fastem5_coef.F90
src/main/parkind1.F90
src/main/rttov_ad.F90
src/main/rttov_add_aux_prof.F90
src/main/rttov_add_opdp_path.F90
src/main/rttov_add_prof.F90
src/main/rttov_add_raytracing.F90
src/main/rttov_alloc_ad.F90
src/main/rttov_alloc_aux_prof.F90
src/main/rttov_alloc_auxrad.F90
src/main/rttov_alloc_auxrad_stream.F90
src/main/rttov_alloc_direct.F90
src/main/rttov_alloc_ircld.F90
src/main/rttov_alloc_k.F90
src/main/rttov_alloc_opdp_path.F90
src/main/rttov_alloc_opt_param.F90
src/main/rttov_alloc_pc_dimensions.F90
src/main/rttov_alloc_pccomp.F90
src/main/rttov_alloc_predictor.F90
src/main/rttov_alloc_prof.F90
src/main/rttov_alloc_rad.F90
src/main/rttov_alloc_raytracing.F90
src/main/rttov_alloc_sunglint.F90
src/main/rttov_alloc_tl.F90
src/main/rttov_alloc_traj.F90
src/main/rttov_alloc_traj_dyn.F90
src/main/rttov_alloc_traj_sta.F90
src/main/rttov_alloc_trans_scatt_ir.F90
src/main/rttov_alloc_transmission.F90
src/main/rttov_alloc_transmission_aux.F90
src/main/rttov_apply_reg_limits.F90
src/main/rttov_apply_reg_limits_ad.F90
src/main/rttov_apply_reg_limits_k.F90
src/main/rttov_apply_reg_limits_tl.F90
src/main/rttov_baran2013_calc_optpar.F90
src/main/rttov_baran2013_calc_optpar_ad.F90
src/main/rttov_baran2013_calc_optpar_tl.F90
src/main/rttov_baran2014_calc_optpar.F90
src/main/rttov_baran2014_calc_optpar_ad.F90
src/main/rttov_baran2014_calc_optpar_tl.F90
src/main/rttov_baran_calc_phase.F90
src/main/rttov_baran_calc_phase_ad.F90
src/main/rttov_baran_calc_phase_tl.F90
src/main/rttov_calc_solar_spec_esd.F90
src/main/rttov_calcbt.F90
src/main/rttov_calcbt_ad.F90
src/main/rttov_calcbt_pc.F90
src/main/rttov_calcbt_pc_ad.F90
src/main/rttov_calcbt_pc_tl.F90
src/main/rttov_calcbt_tl.F90
src/main/rttov_calcemis_ir.F90
src/main/rttov_calcemis_ir_ad.F90
src/main/rttov_calcemis_ir_k.F90
src/main/rttov_calcemis_ir_tl.F90
src/main/rttov_calcemis_mw.F90
src/main/rttov_calcemis_mw_ad.F90
src/main/rttov_calcemis_mw_k.F90
src/main/rttov_calcemis_mw_tl.F90
src/main/rttov_calcrad.F90
src/main/rttov_calcrad_ad.F90
src/main/rttov_calcrad_k.F90
src/main/rttov_calcrad_tl.F90
src/main/rttov_calcsatrefl.F90
src/main/rttov_calcsatrefl_ad.F90
src/main/rttov_calcsatrefl_tl.F90
src/main/rttov_calcsurfrefl.F90
src/main/rttov_calcsurfrefl_ad.F90
src/main/rttov_calcsurfrefl_k.F90
src/main/rttov_calcsurfrefl_tl.F90
src/main/rttov_check_traj.F90
src/main/rttov_checkinput.F90
src/main/rttov_checkpcchan.F90
src/main/rttov_cldstr.F90
src/main/rttov_cldstr_ad.F90
src/main/rttov_cldstr_k.F90
src/main/rttov_cldstr_tl.F90
src/main/rttov_const.F90
src/main/rttov_convert_gas_dry.F90
src/main/rttov_convert_gas_dry_ad.F90
src/main/rttov_convert_gas_dry_k.F90
src/main/rttov_convert_gas_dry_tl.F90
src/main/rttov_copy_aux_prof.F90
src/main/rttov_copy_opdp_path.F90
src/main/rttov_copy_pccomp.F90
src/main/rttov_copy_prof.F90
src/main/rttov_copy_rad.F90
src/main/rttov_copy_raytracing.F90
src/main/rttov_direct.F90
src/main/rttov_erfcx.F90
src/main/rttov_errorhandling.F90
src/main/rttov_errorreport.F90
src/main/rttov_fastem5.F90
src/main/rttov_fastem5_ad.F90
src/main/rttov_fastem5_k.F90
src/main/rttov_fastem5_tl.F90
src/main/rttov_fresnel.F90
src/main/rttov_fresnel_ad.F90
src/main/rttov_fresnel_k.F90
src/main/rttov_fresnel_tl.F90
src/main/rttov_getoptions.F90
src/main/rttov_global.F90
src/main/rttov_init_aux_prof.F90
src/main/rttov_init_auxrad_stream.F90
src/main/rttov_init_ircld.F90
src/main/rttov_init_opdp_path.F90
src/main/rttov_init_opt_param.F90
src/main/rttov_init_pccomp.F90
src/main/rttov_init_predictor.F90
src/main/rttov_init_prof.F90
src/main/rttov_init_rad.F90
src/main/rttov_init_raytracing.F90
src/main/rttov_init_sunglint.F90
src/main/rttov_init_trans_scatt_ir.F90
src/main/rttov_init_transmission.F90
src/main/rttov_init_transmission_aux.F90
src/main/rttov_intavg_chan.F90
src/main/rttov_intavg_chan_ad.F90
src/main/rttov_intavg_chan_k.F90
src/main/rttov_intavg_chan_tl.F90
src/main/rttov_intavg_prof.F90
src/main/rttov_intavg_prof_ad.F90
src/main/rttov_intavg_prof_k.F90
src/main/rttov_intavg_prof_tl.F90
src/main/rttov_integrate.F90
src/main/rttov_integrate_ad.F90
src/main/rttov_integrate_k.F90
src/main/rttov_integrate_tl.F90
src/main/rttov_k.F90
src/main/rttov_layeravg.F90
src/main/rttov_layeravg_ad.F90
src/main/rttov_layeravg_k.F90
src/main/rttov_layeravg_tl.F90
src/main/rttov_locpat.F90
src/main/rttov_locpat_ad.F90
src/main/rttov_locpat_k.F90
src/main/rttov_locpat_tl.F90
src/main/rttov_math_mod.F90
src/main/rttov_mult_profiles_k.F90
src/main/rttov_nlte_bias_correction.F90
src/main/rttov_nlte_bias_correction_ad.F90
src/main/rttov_nlte_bias_correction_k.F90
src/main/rttov_nlte_bias_correction_tl.F90
src/main/rttov_opdep.F90
src/main/rttov_opdep_9.F90
src/main/rttov_opdep_9_ad.F90
src/main/rttov_opdep_9_k.F90
src/main/rttov_opdep_9_tl.F90
src/main/rttov_opdep_ad.F90
src/main/rttov_opdep_k.F90
src/main/rttov_opdep_tl.F90
src/main/rttov_opdpscattir.F90
src/main/rttov_opdpscattir_ad.F90
src/main/rttov_opdpscattir_k.F90
src/main/rttov_opdpscattir_tl.F90
src/main/rttov_opts_eq.F90
src/main/rttov_pcscores.F90
src/main/rttov_pcscores_ad.F90
src/main/rttov_pcscores_k.F90
src/main/rttov_pcscores_rec_k.F90
src/main/rttov_pcscores_tl.F90
src/main/rttov_profaux.F90
src/main/rttov_profaux_ad.F90
src/main/rttov_profaux_k.F90
src/main/rttov_profaux_tl.F90
src/main/rttov_reconstruct.F90
src/main/rttov_reconstruct_ad.F90
src/main/rttov_reconstruct_k.F90
src/main/rttov_reconstruct_tl.F90
src/main/rttov_refsun.F90
src/main/rttov_refsun_ad.F90
src/main/rttov_refsun_k.F90
src/main/rttov_refsun_tl.F90
src/main/rttov_setgeometry.F90
src/main/rttov_setgeometry_ad.F90
src/main/rttov_setgeometry_k.F90
src/main/rttov_setgeometry_tl.F90
src/main/rttov_setpredictors_7.F90
src/main/rttov_setpredictors_7_ad.F90
src/main/rttov_setpredictors_7_k.F90
src/main/rttov_setpredictors_7_tl.F90
src/main/rttov_setpredictors_8.F90
src/main/rttov_setpredictors_8_ad.F90
src/main/rttov_setpredictors_8_k.F90
src/main/rttov_setpredictors_8_tl.F90
src/main/rttov_setpredictors_9.F90
src/main/rttov_setpredictors_9_ad.F90
src/main/rttov_setpredictors_9_k.F90
src/main/rttov_setpredictors_9_tl.F90
src/main/rttov_sublayer.F90
src/main/rttov_sublayer_ad.F90
src/main/rttov_sublayer_k.F90
src/main/rttov_sublayer_tl.F90
src/main/rttov_tl.F90
src/main/rttov_transmit.F90
src/main/rttov_transmit_9_solar.F90
src/main/rttov_transmit_9_solar_ad.F90
src/main/rttov_transmit_9_solar_k.F90
src/main/rttov_transmit_9_solar_tl.F90
src/main/rttov_transmit_ad.F90
src/main/rttov_transmit_k.F90
src/main/rttov_transmit_tl.F90
src/main/rttov_types.F90
src/main/rttov_unix_env.F90
src/main/rttov_user_options_checkinput.F90
src/main/rttov_user_profile_checkinput.F90
src/main/throw.h
src/main/yomhook.F90
src/mw_scatt/
src/mw_scatt/Makefile
src/mw_scatt/example_rttovscatt.F90
src/mw_scatt/mod_rttovscatt_test.F90
src/mw_scatt/rttov_alloc_scatt_prof.F90
src/mw_scatt/rttov_boundaryconditions.F90
src/mw_scatt/rttov_boundaryconditions_ad.F90
src/mw_scatt/rttov_boundaryconditions_tl.F90
src/mw_scatt/rttov_dealloc_scattcoeffs.F90
src/mw_scatt/rttov_eddington.F90
src/mw_scatt/rttov_eddington_ad.F90
src/mw_scatt/rttov_eddington_tl.F90
src/mw_scatt/rttov_hydro.F90
src/mw_scatt/rttov_hydro_ad.F90
src/mw_scatt/rttov_hydro_tl.F90
src/mw_scatt/rttov_iniedd.F90
src/mw_scatt/rttov_iniedd_ad.F90
src/mw_scatt/rttov_iniedd_tl.F90
src/mw_scatt/rttov_iniscatt.F90
src/mw_scatt/rttov_iniscatt_ad.F90
src/mw_scatt/rttov_iniscatt_tl.F90
src/mw_scatt/rttov_init_scatt_prof.F90
src/mw_scatt/rttov_integratesource.F90
src/mw_scatt/rttov_integratesource_ad.F90
src/mw_scatt/rttov_integratesource_tl.F90
src/mw_scatt/rttov_mieproc.F90
src/mw_scatt/rttov_mieproc_ad.F90
src/mw_scatt/rttov_mieproc_tl.F90
src/mw_scatt/rttov_read_scattcoeffs.F90
src/mw_scatt/rttov_scatt.F90
src/mw_scatt/rttov_scatt_ad.F90
src/mw_scatt/rttov_scatt_setupindex.F90
src/mw_scatt/rttov_scatt_tl.F90
src/mw_scatt/rttov_setpressure.F90
src/mw_scatt/rttovscatt_test.F90
src/mw_scatt/rttovscatt_test_one.F90
src/mw_scatt_coef/
src/mw_scatt_coef/Makefile
src/mw_scatt_coef/channels.dat
src/mw_scatt_coef/channels.dat_all
src/mw_scatt_coef/channels.dat_amsua
src/mw_scatt_coef/channels.dat_atovs
src/mw_scatt_coef/channels.dat_debug
src/mw_scatt_coef/channels.dat_mwhs2
src/mw_scatt_coef/channels.dat_mwts2
src/mw_scatt_coef/channels.dat_smos
src/mw_scatt_coef/channels.dat_ssmis
src/mw_scatt_coef/convert_mietable.F90
src/mw_scatt_coef/density_all.F90
src/mw_scatt_coef/gamma_dsd.F90
src/mw_scatt_coef/get_dia.F90
src/mw_scatt_coef/ice_density.F90
src/mw_scatt_coef/liu_dda.F90
src/mw_scatt_coef/liu_density.F90
src/mw_scatt_coef/melting_layer.F90
src/mw_scatt_coef/mg_ellips.F90
src/mw_scatt_coef/mie_coated_sphere.F90
src/mw_scatt_coef/mie_one_temp.F90
src/mw_scatt_coef/mie_one_wc.F90
src/mw_scatt_coef/mie_sphere.F90
src/mw_scatt_coef/mie_table_generation.ksh
src/mw_scatt_coef/mod_gamma_dsd.F90
src/mw_scatt_coef/mod_mie.F90
src/mw_scatt_coef/n0_t.F90
src/mw_scatt_coef/perm_ice.F90
src/mw_scatt_coef/perm_melt.F90
src/mw_scatt_coef/perm_water.F90
src/mw_scatt_coef/permittivity.F90
src/mw_scatt_coef/predict_mom07.F90
src/mw_scatt_coef/predict_psd.F90
src/mw_scatt_coef/predict_psd_F07.F90
src/mw_scatt_coef/readme.txt
src/mw_scatt_coef/rttov_ascii2bin_scattcoef.F90
src/mw_scatt_coef/rttov_scatt_make_coef.F90
src/mw_scatt_coef/scat_db2.dda
src/mw_scatt_coef/scatdb.c
src/mw_scatt_coef/scattering.F90
src/mw_scatt_coef/set_spectra.F90
src/mw_scatt_coef/vol_fracs.F90
src/other/
src/other/Makefile
src/other/create_aer_clim_prof.F90
src/other/rttov_aer_clim_prof.F90
src/other/rttov_bpr_calc.F90
src/other/rttov_bpr_dealloc.F90
src/other/rttov_bpr_init.F90
src/other/rttov_bpr_mod.F90
src/other/rttov_calc_weighting_fn.F90
src/other/rttov_coef_info.F90
src/other/rttov_mie_params_aer.F90
src/other/rttov_mie_params_cld.F90
src/other/rttov_obs_to_pc.F90
src/other/rttov_print_info.F90
src/other/rttov_print_opts.F90
src/other/rttov_print_profile.F90
src/other/rttov_scattering_mod.F90
src/other/rttov_us76_prof.F90
src/other/rttov_zutility.F90
src/parallel/
src/parallel/Makefile
src/parallel/rttov_parallel_ad.F90
src/parallel/rttov_parallel_direct.F90
src/parallel/rttov_parallel_k.F90
src/parallel/rttov_parallel_tl.F90
src/test/
src/test/Makefile
src/test/Makefile_examples
src/test/example_aer_file_fwd.F90
src/test/example_aer_param_fwd.F90
src/test/example_cld_file_fwd.F90
src/test/example_cld_param_fwd.F90
src/test/example_fwd.F90
src/test/example_k.F90
src/test/example_pc_fwd.F90
src/test/example_rttovscatt_fwd.F90
src/test/rttov_chain.F90
src/test/rttov_k_ad.F90
src/test/rttov_k_bf.F90
src/test/rttov_k_tl.F90
src/test/rttov_lun.F90
src/test/rttov_make_pccomp_inc.F90
src/test/rttov_make_profile_inc.F90
src/test/rttov_make_radiance_inc.F90
src/test/rttov_scale_pccomp_inc.F90
src/test/rttov_scale_profile_inc.F90
src/test/rttov_scale_radiance_inc.F90
src/test/rttov_test.F90
src/test/rttov_test_k_mod.F90
src/wrapper/
src/wrapper/Makefile
src/wrapper/f2py_f2cmap
src/wrapper/rttov_c_interface.h
src/wrapper/rttov_wrapper_f2py.F90
src/wrapper/rttov_wrapper_handle.F90
src/wrapper/rttov_wrapper_transfer.F90
brdf_data/
build/
build/Makefile.PL
build/Makefile.inc
build/Makefile.local
build/arch/
build/arch/aix
build/arch/aix-debug
build/arch/aix-ops
build/arch/cray-ecmwf
build/arch/cray-mo-ifort
build/arch/gfortran
build/arch/gfortran-debug
build/arch/gfortran-openmp
build/arch/ifort
build/arch/ifort-debug
build/arch/ifort-mf
build/arch/ifort-openmp
build/arch/ifort-ops
build/arch/mpaix
build/arch/nag5.1
build/arch/nag5.1-debug
build/arch/nagfor
build/arch/nagfor-debug
build/arch/nagfor-openmp
build/arch/nec-meteofrance
build/arch/pgf90
build/arch/pgf90-debug
build/arch/pgf90-openmp
build/arch/solaris
build/arch/solaris-debug
build/arch/sunstudio-debug
build/cpinch.pl
build/mkintf.pl
build/mvdmod.pl
build/mypcpp.pl
build/rttov_compile.sh
data/
data/Be_LUT.2007.txt
data/asdu00
data/iasi_pc_band_2_chans.txt
data/iasi_pc_band_3_chans.txt
data/plevs.dat
data/prof.dat
data/prof_aerosl_cl.dat
data/vapo00
emis_data/
gui/
gui/README
gui/doc/
gui/doc/helpDiffRad.html
gui/doc/helpKMatrixFrame.html
gui/doc/helpKPC.html
gui/doc/helpKpcMatrixFrame.html
gui/doc/helpOptions.html
gui/doc/helpPC.html
gui/doc/helpR1DVAR.html
gui/doc/helpRadianceFrame.html
gui/icons/
gui/icons/CH4.png
gui/icons/CO.png
gui/icons/CO2.png
gui/icons/N2O.png
gui/icons/NO2.png
gui/icons/O3.png
gui/icons/Q.png
gui/icons/exit.png
gui/icons/fileclose.png
gui/icons/hand.png
gui/icons/help.png
gui/icons/k10.png
gui/icons/kP.png
gui/icons/kmat_O3.png
gui/icons/kmat_Q.png
gui/icons/kmat_T.png
gui/icons/kmatrix_toolbar.png
gui/icons/matplotlib_toolbar.png
gui/icons/reset.png
gui/icons/right.png
gui/icons/stock_down.png
gui/icons/stock_up.png
gui/r1Dvar/
gui/r1Dvar/__init__.py
gui/r1Dvar/data/
gui/r1Dvar/data/AIRS_COEFFS_DIR/
gui/r1Dvar/data/AIRS_COEFFS_DIR/ChannelChoice_orig.dat
gui/r1Dvar/data/AIRS_COEFFS_DIR/Rmatrix_orig
gui/r1Dvar/data/ATMS_COEFFS_DIR/
gui/r1Dvar/data/ATMS_COEFFS_DIR/ChannelChoice_orig.dat
gui/r1Dvar/data/ATMS_COEFFS_DIR/Rmatrix_orig
gui/r1Dvar/data/ATOVS_CLOUDY_COEFFS_DIR/
gui/r1Dvar/data/ATOVS_CLOUDY_COEFFS_DIR/ChannelChoice_orig.dat
gui/r1Dvar/data/ATOVS_CLOUDY_COEFFS_DIR/Rmatrix_orig
gui/r1Dvar/data/ATOVS_COEFFS_DIR/
gui/r1Dvar/data/ATOVS_COEFFS_DIR/ChannelChoice_orig.dat
gui/r1Dvar/data/ATOVS_COEFFS_DIR/Rmatrix_orig
gui/r1Dvar/data/CrIS_COEFFS_DIR/
gui/r1Dvar/data/CrIS_COEFFS_DIR/ChannelChoice_orig.dat
gui/r1Dvar/data/CrIS_COEFFS_DIR/Rmatrix_orig
gui/r1Dvar/data/IASI_COEFFS_DIR/
gui/r1Dvar/data/IASI_COEFFS_DIR/Bmatrix
gui/r1Dvar/data/IASI_COEFFS_DIR/ChannelChoice.dat
gui/r1Dvar/data/IASI_COEFFS_DIR/ChannelChoice_orig.dat
gui/r1Dvar/data/IASI_COEFFS_DIR/Rmatrix
gui/r1Dvar/data/IASI_COEFFS_DIR/Rmatrix_orig
gui/r1Dvar/data/SSMIS_COEFFS_DIR/
gui/r1Dvar/data/SSMIS_COEFFS_DIR/ChannelChoice_orig.dat
gui/r1Dvar/data/SSMIS_COEFFS_DIR/Rmatrix_orig
gui/r1Dvar/data/Sample_Background/
gui/r1Dvar/data/Sample_Background/BACKGROUND_43L.dat
gui/r1Dvar/data/Sample_Background/BACKGROUND_51L.dat
gui/r1Dvar/data/Sample_Background/BACKGROUND_54L.dat
gui/r1Dvar/data/Sample_Background/BACKGROUND_with_CLW.dat
gui/r1Dvar/data/Sample_Background/truth_43L.dat
gui/r1Dvar/data/Sample_Background/truth_51L.dat
gui/r1Dvar/data/Sample_Background/truth_54L.dat
gui/r1Dvar/data/Sample_Bmatrices/
gui/r1Dvar/data/Sample_Bmatrices/Bmatrix_43L
gui/r1Dvar/data/Sample_Bmatrices/Bmatrix_51L
gui/r1Dvar/data/Sample_Bmatrices/Bmatrix_54L
gui/r1Dvar/r1dvar.py
gui/r1Dvar/r1dvarObjects.py
gui/rcontroller/
gui/rcontroller/__init__.py
gui/rcontroller/controller.py
gui/rcontroller/optionctrl.py
gui/rcontroller/profilectrl.py
gui/rcontroller/r1dvarController.py
gui/rcontroller/surfacectrl.py
gui/rcontroller/util.py
gui/rmodel/
gui/rmodel/__init__.py
gui/rmodel/config.py
gui/rmodel/project.py
gui/rttov/
gui/rttov/__init__.py
gui/rttov/chanprof.py
gui/rttov/core.py
gui/rttov/core.pyc
gui/rttov/emissivity.py
gui/rttov/getcoefval.py
gui/rttov/kmatrix.py
gui/rttov/kpcmatrix.py
gui/rttov/list_of_profile_vars.py
gui/rttov/misc.py
gui/rttov/option.py
gui/rttov/pccomp.py
gui/rttov/profile.py
gui/rttov/profile.pyc
gui/rttov/radiance.py
gui/rttov/reflectance.py
gui/rttov/transmission.py
gui/rttov_gui.env
gui/rttov_gui_f2py.so
gui/rttovgui
gui/run
gui/rview/
gui/rview/__init__.py
gui/rview/coeff.py
gui/rview/colors.py
gui/rview/console.py
gui/rview/helpframe.py
gui/rview/kmatrixframe.py
gui/rview/kpcView.py
gui/rview/kpcmatrixframe.py
gui/rview/kprofileframe.py
gui/rview/layeritem.py
gui/rview/option.py
gui/rview/option_help.html
gui/rview/pcView.py
gui/rview/profileframe.py
gui/rview/profileframeutils.py
gui/rview/r1dvarView.py
gui/rview/r1dvarprofileframe.py
gui/rview/rBtView.py
gui/rview/radianceframe.py
gui/rview/surface.py
gui/rview/surfedit.py
gui/rview/util.py
gui/rview/wxmpl.py
gui/test/
gui/test/__init__.py
gui/test/rttovgui_unittest_class.py
gui/test/test_1dvar.py
gui/test/test_aerosols.py
gui/test/test_atlas.py
gui/test/test_clouds.py
gui/test/test_full.py
gui/test/test_nlte.py
gui/test/test_options.py
gui/test/test_ozone.py
gui/test/test_retrieve.py
gui/test/test_retrieve_iasi.py
gui/test/test_run.py
gui/test/test_runpc.py
gui/test/unitest_project.py
rtcoef_rttov11/
rtcoef_rttov11/cldaer/
rtcoef_rttov11/make_metopb_iasi.pl
rtcoef_rttov11/mietable/
rtcoef_rttov11/pc/
rtcoef_rttov11/rttov7pred101L/
rtcoef_rttov11/rttov7pred54L/
rtcoef_rttov11/rttov7pred54L/rtcoef_calipso_1_iir.dat
rtcoef_rttov11/rttov7pred54L/rtcoef_coms_1_mi.dat
rtcoef_rttov11/rttov7pred54L/rtcoef_coriolis_1_windsat.dat
rtcoef_rttov11/rttov7pred54L/rtcoef_dmsp_10_ssmi.dat
rtcoef_rttov11/rttov7pred54L/rtcoef_dmsp_11_ssmi.dat
rtcoef_rttov11/rttov7pred54L/rtcoef_dmsp_11_ssmt2.dat
rtcoef_rttov11/rttov7pred54L/rtcoef_dmsp_12_ssmi.dat
rtcoef_rttov11/rttov7pred54L/rtcoef_dmsp_12_ssmt2.dat
rtcoef_rttov11/rttov7pred54L/rtcoef_dmsp_13_ssmi.dat
rtcoef_rttov11/rttov7pred54L/rtcoef_dmsp_14_ssmi.dat
rtcoef_rttov11/rttov7pred54L/rtcoef_dmsp_14_ssmt2.dat
rtcoef_rttov11/rttov7pred54L/rtcoef_dmsp_15_ssmi.dat
rtcoef_rttov11/rttov7pred54L/rtcoef_dmsp_15_ssmt2.dat
rtcoef_rttov11/rttov7pred54L/rtcoef_dmsp_16_ssmis.dat
rtcoef_rttov11/rttov7pred54L/rtcoef_dmsp_17_ssmis.dat
rtcoef_rttov11/rttov7pred54L/rtcoef_dmsp_18_ssmis.dat
rtcoef_rttov11/rttov7pred54L/rtcoef_dmsp_19_ssmis.dat
rtcoef_rttov11/rttov7pred54L/rtcoef_dmsp_8_ssmi.dat
rtcoef_rttov11/rttov7pred54L/rtcoef_dmsp_9_ssmi.dat
rtcoef_rttov11/rttov7pred54L/rtcoef_envisat_1_atsr-shifted.dat
rtcoef_rttov11/rttov7pred54L/rtcoef_envisat_1_atsr.dat
rtcoef_rttov11/rttov7pred54L/rtcoef_eos_1_aster.dat
rtcoef_rttov11/rttov7pred54L/rtcoef_eos_1_modis-shifted.dat
rtcoef_rttov11/rttov7pred54L/rtcoef_eos_1_modis.dat
rtcoef_rttov11/rttov7pred54L/rtcoef_eos_2_amsre.dat
rtcoef_rttov11/rttov7pred54L/rtcoef_eos_2_amsua.dat
rtcoef_rttov11/rttov7pred54L/rtcoef_eos_2_hsb.dat
rtcoef_rttov11/rttov7pred54L/rtcoef_eos_2_modis-shifted.dat
rtcoef_rttov11/rttov7pred54L/rtcoef_eos_2_modis.dat
rtcoef_rttov11/rttov7pred54L/rtcoef_ers_1_atsr.dat
rtcoef_rttov11/rttov7pred54L/rtcoef_ers_2_atsr.dat
rtcoef_rttov11/rttov7pred54L/rtcoef_fy1_3_mvisr.dat
rtcoef_rttov11/rttov7pred54L/rtcoef_fy1_4_mvisr.dat
rtcoef_rttov11/rttov7pred54L/rtcoef_fy2_2_vissr.dat
rtcoef_rttov11/rttov7pred54L/rtcoef_fy2_3_vissr.dat
rtcoef_rttov11/rttov7pred54L/rtcoef_fy2_4_vissr.dat
rtcoef_rttov11/rttov7pred54L/rtcoef_fy3_1_iras.dat
rtcoef_rttov11/rttov7pred54L/rtcoef_fy3_1_mwhs.dat
rtcoef_rttov11/rttov7pred54L/rtcoef_fy3_1_mwri.dat
rtcoef_rttov11/rttov7pred54L/rtcoef_fy3_1_mwts.dat
rtcoef_rttov11/rttov7pred54L/rtcoef_fy3_2_mwhs.dat
rtcoef_rttov11/rttov7pred54L/rtcoef_fy3_2_mwri.dat
rtcoef_rttov11/rttov7pred54L/rtcoef_fy3_2_mwts.dat
rtcoef_rttov11/rttov7pred54L/rtcoef_fy3_3_mwhs2.dat
rtcoef_rttov11/rttov7pred54L/rtcoef_fy3_3_mwri.dat
rtcoef_rttov11/rttov7pred54L/rtcoef_fy3_3_mwts2.dat
rtcoef_rttov11/rttov7pred54L/rtcoef_gcom-w_1_amsr2.dat
rtcoef_rttov11/rttov7pred54L/rtcoef_gms_5_imager.dat
rtcoef_rttov11/rttov7pred54L/rtcoef_goes_10_imager.dat
rtcoef_rttov11/rttov7pred54L/rtcoef_goes_10_sounder.dat
rtcoef_rttov11/rttov7pred54L/rtcoef_goes_11_imager.dat
rtcoef_rttov11/rttov7pred54L/rtcoef_goes_11_sounder.dat
rtcoef_rttov11/rttov7pred54L/rtcoef_goes_12_imager.dat
rtcoef_rttov11/rttov7pred54L/rtcoef_goes_12_sounder.dat
rtcoef_rttov11/rttov7pred54L/rtcoef_goes_13_imager.dat
rtcoef_rttov11/rttov7pred54L/rtcoef_goes_13_sounder.dat
rtcoef_rttov11/rttov7pred54L/rtcoef_goes_14_imager.dat
rtcoef_rttov11/rttov7pred54L/rtcoef_goes_14_sounder.dat
rtcoef_rttov11/rttov7pred54L/rtcoef_goes_15_imager.dat
rtcoef_rttov11/rttov7pred54L/rtcoef_goes_15_sounder.dat
rtcoef_rttov11/rttov7pred54L/rtcoef_goes_16_abi.dat
rtcoef_rttov11/rttov7pred54L/rtcoef_goes_4_sounder.dat
rtcoef_rttov11/rttov7pred54L/rtcoef_goes_5_sounder.dat
rtcoef_rttov11/rttov7pred54L/rtcoef_goes_6_sounder.dat
rtcoef_rttov11/rttov7pred54L/rtcoef_goes_7_sounder.dat
rtcoef_rttov11/rttov7pred54L/rtcoef_goes_8_imager.dat
rtcoef_rttov11/rttov7pred54L/rtcoef_goes_8_sounder.dat
rtcoef_rttov11/rttov7pred54L/rtcoef_goes_9_imager.dat
rtcoef_rttov11/rttov7pred54L/rtcoef_goes_9_sounder.dat
rtcoef_rttov11/rttov7pred54L/rtcoef_gpm_1_gmi.dat
rtcoef_rttov11/rttov7pred54L/rtcoef_himawari_8_ahi.dat
rtcoef_rttov11/rttov7pred54L/rtcoef_jason_2_amr.dat
rtcoef_rttov11/rttov7pred54L/rtcoef_jpss_0_atms.dat
rtcoef_rttov11/rttov7pred54L/rtcoef_jpss_0_viirs.dat
rtcoef_rttov11/rttov7pred54L/rtcoef_landsat_4_tm.dat
rtcoef_rttov11/rttov7pred54L/rtcoef_landsat_5_tm.dat
rtcoef_rttov11/rttov7pred54L/rtcoef_landsat_7_tm.dat
rtcoef_rttov11/rttov7pred54L/rtcoef_landsat_8_tirs.dat
rtcoef_rttov11/rttov7pred54L/rtcoef_meghatr_1_madras.dat
rtcoef_rttov11/rttov7pred54L/rtcoef_meghatr_1_saphir.dat
rtcoef_rttov11/rttov7pred54L/rtcoef_meteor-m_1_msumr.dat
rtcoef_rttov11/rttov7pred54L/rtcoef_meteor-m_2_mtvzagy.dat
rtcoef_rttov11/rttov7pred54L/rtcoef_meteosat_1_mviri.dat
rtcoef_rttov11/rttov7pred54L/rtcoef_meteosat_2_mviri.dat
rtcoef_rttov11/rttov7pred54L/rtcoef_meteosat_3_mviri.dat
rtcoef_rttov11/rttov7pred54L/rtcoef_meteosat_4_mviri.dat
rtcoef_rttov11/rttov7pred54L/rtcoef_meteosat_5_mviri.dat
rtcoef_rttov11/rttov7pred54L/rtcoef_meteosat_6_mviri.dat
rtcoef_rttov11/rttov7pred54L/rtcoef_meteosat_7_mviri.dat
rtcoef_rttov11/rttov7pred54L/rtcoef_metop_1_amsua.dat
rtcoef_rttov11/rttov7pred54L/rtcoef_metop_1_avhrr.dat
rtcoef_rttov11/rttov7pred54L/rtcoef_metop_1_hirs.dat
rtcoef_rttov11/rttov7pred54L/rtcoef_metop_1_mhs.dat
rtcoef_rttov11/rttov7pred54L/rtcoef_metop_2_amsua.dat
rtcoef_rttov11/rttov7pred54L/rtcoef_metop_2_avhrr.dat
rtcoef_rttov11/rttov7pred54L/rtcoef_metop_2_hirs-shifted.dat
rtcoef_rttov11/rttov7pred54L/rtcoef_metop_2_hirs.dat
rtcoef_rttov11/rttov7pred54L/rtcoef_metop_2_mhs.dat
rtcoef_rttov11/rttov7pred54L/rtcoef_metopsg_1_ici.dat
rtcoef_rttov11/rttov7pred54L/rtcoef_metopsg_1_metimage.dat
rtcoef_rttov11/rttov7pred54L/rtcoef_metopsg_1_mwi.dat
rtcoef_rttov11/rttov7pred54L/rtcoef_metopsg_1_mws.dat
rtcoef_rttov11/rttov7pred54L/rtcoef_msg_1_seviri.dat
rtcoef_rttov11/rttov7pred54L/rtcoef_msg_2_seviri.dat
rtcoef_rttov11/rttov7pred54L/rtcoef_msg_3_seviri.dat
rtcoef_rttov11/rttov7pred54L/rtcoef_msg_4_seviri.dat
rtcoef_rttov11/rttov7pred54L/rtcoef_mtg_1_fci.dat
rtcoef_rttov11/rttov7pred54L/rtcoef_mtsat_1_imager.dat
rtcoef_rttov11/rttov7pred54L/rtcoef_mtsat_2_imager.dat
rtcoef_rttov11/rttov7pred54L/rtcoef_nimbus_6_hirs.dat
rtcoef_rttov11/rttov7pred54L/rtcoef_nimbus_7_smmr.dat
rtcoef_rttov11/rttov7pred54L/rtcoef_noaa_10_avhrr.dat
rtcoef_rttov11/rttov7pred54L/rtcoef_noaa_10_hirs-shifted.dat
rtcoef_rttov11/rttov7pred54L/rtcoef_noaa_10_hirs.dat
rtcoef_rttov11/rttov7pred54L/rtcoef_noaa_10_msu.dat
rtcoef_rttov11/rttov7pred54L/rtcoef_noaa_11_avhrr.dat
rtcoef_rttov11/rttov7pred54L/rtcoef_noaa_11_hirs-shifted.dat
rtcoef_rttov11/rttov7pred54L/rtcoef_noaa_11_hirs.dat
rtcoef_rttov11/rttov7pred54L/rtcoef_noaa_11_msu.dat
rtcoef_rttov11/rttov7pred54L/rtcoef_noaa_12_avhrr.dat
rtcoef_rttov11/rttov7pred54L/rtcoef_noaa_12_hirs-shifted.dat
rtcoef_rttov11/rttov7pred54L/rtcoef_noaa_12_hirs.dat
rtcoef_rttov11/rttov7pred54L/rtcoef_noaa_12_msu.dat
rtcoef_rttov11/rttov7pred54L/rtcoef_noaa_13_avhrr.dat
rtcoef_rttov11/rttov7pred54L/rtcoef_noaa_14_avhrr.dat
rtcoef_rttov11/rttov7pred54L/rtcoef_noaa_14_hirs-shifted.dat
rtcoef_rttov11/rttov7pred54L/rtcoef_noaa_14_hirs.dat
rtcoef_rttov11/rttov7pred54L/rtcoef_noaa_14_msu.dat
rtcoef_rttov11/rttov7pred54L/rtcoef_noaa_15_amsua.dat
rtcoef_rttov11/rttov7pred54L/rtcoef_noaa_15_amsub.dat
rtcoef_rttov11/rttov7pred54L/rtcoef_noaa_15_avhrr.dat
rtcoef_rttov11/rttov7pred54L/rtcoef_noaa_15_hirs-shifted.dat
rtcoef_rttov11/rttov7pred54L/rtcoef_noaa_15_hirs.dat
rtcoef_rttov11/rttov7pred54L/rtcoef_noaa_16_amsua.dat
rtcoef_rttov11/rttov7pred54L/rtcoef_noaa_16_amsub.dat
rtcoef_rttov11/rttov7pred54L/rtcoef_noaa_16_avhrr.dat
rtcoef_rttov11/rttov7pred54L/rtcoef_noaa_16_hirs-shifted.dat
rtcoef_rttov11/rttov7pred54L/rtcoef_noaa_16_hirs.dat
rtcoef_rttov11/rttov7pred54L/rtcoef_noaa_17_amsua.dat
rtcoef_rttov11/rttov7pred54L/rtcoef_noaa_17_amsub.dat
rtcoef_rttov11/rttov7pred54L/rtcoef_noaa_17_avhrr.dat
rtcoef_rttov11/rttov7pred54L/rtcoef_noaa_17_hirs-shifted.dat
rtcoef_rttov11/rttov7pred54L/rtcoef_noaa_17_hirs.dat
rtcoef_rttov11/rttov7pred54L/rtcoef_noaa_18_amsua.dat
rtcoef_rttov11/rttov7pred54L/rtcoef_noaa_18_avhrr.dat
rtcoef_rttov11/rttov7pred54L/rtcoef_noaa_18_hirs-shifted.dat
rtcoef_rttov11/rttov7pred54L/rtcoef_noaa_18_hirs.dat
rtcoef_rttov11/rttov7pred54L/rtcoef_noaa_18_mhs.dat
rtcoef_rttov11/rttov7pred54L/rtcoef_noaa_19_amsua.dat
rtcoef_rttov11/rttov7pred54L/rtcoef_noaa_19_avhrr.dat
rtcoef_rttov11/rttov7pred54L/rtcoef_noaa_19_hirs-shifted.dat
rtcoef_rttov11/rttov7pred54L/rtcoef_noaa_19_hirs.dat
rtcoef_rttov11/rttov7pred54L/rtcoef_noaa_19_mhs.dat
rtcoef_rttov11/rttov7pred54L/rtcoef_noaa_1_vtpr1.dat
rtcoef_rttov11/rttov7pred54L/rtcoef_noaa_2_vtpr1.dat
rtcoef_rttov11/rttov7pred54L/rtcoef_noaa_3_vtpr1.dat
rtcoef_rttov11/rttov7pred54L/rtcoef_noaa_4_vtpr1.dat
rtcoef_rttov11/rttov7pred54L/rtcoef_noaa_5_avhrr.dat
rtcoef_rttov11/rttov7pred54L/rtcoef_noaa_5_hirs.dat
rtcoef_rttov11/rttov7pred54L/rtcoef_noaa_5_msu.dat
rtcoef_rttov11/rttov7pred54L/rtcoef_noaa_6_avhrr.dat
rtcoef_rttov11/rttov7pred54L/rtcoef_noaa_6_hirs.dat
rtcoef_rttov11/rttov7pred54L/rtcoef_noaa_6_msu.dat
rtcoef_rttov11/rttov7pred54L/rtcoef_noaa_7_avhrr.dat
rtcoef_rttov11/rttov7pred54L/rtcoef_noaa_7_hirs.dat
rtcoef_rttov11/rttov7pred54L/rtcoef_noaa_7_msu.dat
rtcoef_rttov11/rttov7pred54L/rtcoef_noaa_8_avhrr.dat
rtcoef_rttov11/rttov7pred54L/rtcoef_noaa_8_hirs.dat
rtcoef_rttov11/rttov7pred54L/rtcoef_noaa_8_msu.dat
rtcoef_rttov11/rttov7pred54L/rtcoef_noaa_9_avhrr.dat
rtcoef_rttov11/rttov7pred54L/rtcoef_noaa_9_hirs-shifted.dat
rtcoef_rttov11/rttov7pred54L/rtcoef_noaa_9_hirs.dat
rtcoef_rttov11/rttov7pred54L/rtcoef_noaa_9_msu.dat
rtcoef_rttov11/rttov7pred54L/rtcoef_saral_1_altika.dat
rtcoef_rttov11/rttov7pred54L/rtcoef_sentinel3_1_slstr.dat
rtcoef_rttov11/rttov7pred54L/rtcoef_trmm_1_tmi.dat
rtcoef_rttov11/rttov8pred101L/
rtcoef_rttov11/rttov8pred51L/
rtcoef_rttov11/rttov8pred51L/rtcoef_noaa_11_ssu_1st.dat
rtcoef_rttov11/rttov8pred51L/rtcoef_noaa_11_ssu_1st_pmcshift.dat
rtcoef_rttov11/rttov8pred51L/rtcoef_noaa_11_ssu_2nd.dat
rtcoef_rttov11/rttov8pred51L/rtcoef_noaa_11_ssu_2nd_pmcshift.dat
rtcoef_rttov11/rttov8pred51L/rtcoef_noaa_14_ssu.dat
rtcoef_rttov11/rttov8pred51L/rtcoef_noaa_14_ssu_pmcshift.dat
rtcoef_rttov11/rttov8pred51L/rtcoef_noaa_5_ssu.dat
rtcoef_rttov11/rttov8pred51L/rtcoef_noaa_5_ssu_pmcshift.dat
rtcoef_rttov11/rttov8pred51L/rtcoef_noaa_6_ssu.dat
rtcoef_rttov11/rttov8pred51L/rtcoef_noaa_6_ssu_pmcshift.dat
rtcoef_rttov11/rttov8pred51L/rtcoef_noaa_7_ssu.dat
rtcoef_rttov11/rttov8pred51L/rtcoef_noaa_7_ssu_pmcshift.dat
rtcoef_rttov11/rttov8pred51L/rtcoef_noaa_8_ssu.dat
rtcoef_rttov11/rttov8pred51L/rtcoef_noaa_8_ssu_pmcshift.dat
rtcoef_rttov11/rttov8pred51L/rtcoef_noaa_9_ssu.dat
rtcoef_rttov11/rttov8pred51L/rtcoef_noaa_9_ssu_pmcshift.dat
rtcoef_rttov11/rttov8pred54L/
rtcoef_rttov11/rttov8pred54L/rtcoef_metop_1_hirs.dat
rtcoef_rttov11/rttov8pred54L/rtcoef_metop_2_hirs-shifted.dat
rtcoef_rttov11/rttov8pred54L/rtcoef_metop_2_hirs.dat
rtcoef_rttov11/rttov8pred54L/rtcoef_nimbus_3_mrir.dat
rtcoef_rttov11/rttov8pred54L/rtcoef_nimbus_6_hirs.dat
rtcoef_rttov11/rttov8pred54L/rtcoef_noaa_10_hirs-shifted.dat
rtcoef_rttov11/rttov8pred54L/rtcoef_noaa_10_hirs.dat
rtcoef_rttov11/rttov8pred54L/rtcoef_noaa_11_hirs-shifted.dat
rtcoef_rttov11/rttov8pred54L/rtcoef_noaa_11_hirs.dat
rtcoef_rttov11/rttov8pred54L/rtcoef_noaa_12_hirs-shifted.dat
rtcoef_rttov11/rttov8pred54L/rtcoef_noaa_12_hirs.dat
rtcoef_rttov11/rttov8pred54L/rtcoef_noaa_14_hirs-shifted.dat
rtcoef_rttov11/rttov8pred54L/rtcoef_noaa_14_hirs.dat
rtcoef_rttov11/rttov8pred54L/rtcoef_noaa_15_hirs-shifted.dat
rtcoef_rttov11/rttov8pred54L/rtcoef_noaa_15_hirs.dat
rtcoef_rttov11/rttov8pred54L/rtcoef_noaa_16_hirs-shifted.dat
rtcoef_rttov11/rttov8pred54L/rtcoef_noaa_16_hirs.dat
rtcoef_rttov11/rttov8pred54L/rtcoef_noaa_17_hirs-shifted.dat
rtcoef_rttov11/rttov8pred54L/rtcoef_noaa_17_hirs.dat
rtcoef_rttov11/rttov8pred54L/rtcoef_noaa_18_hirs-shifted.dat
rtcoef_rttov11/rttov8pred54L/rtcoef_noaa_18_hirs.dat
rtcoef_rttov11/rttov8pred54L/rtcoef_noaa_19_hirs-shifted.dat
rtcoef_rttov11/rttov8pred54L/rtcoef_noaa_19_hirs.dat
rtcoef_rttov11/rttov8pred54L/rtcoef_noaa_5_hirs.dat
rtcoef_rttov11/rttov8pred54L/rtcoef_noaa_6_hirs.dat
rtcoef_rttov11/rttov8pred54L/rtcoef_noaa_7_hirs.dat
rtcoef_rttov11/rttov8pred54L/rtcoef_noaa_8_hirs.dat
rtcoef_rttov11/rttov8pred54L/rtcoef_noaa_9_hirs-shifted.dat
rtcoef_rttov11/rttov8pred54L/rtcoef_noaa_9_hirs.dat
rtcoef_rttov11/rttov9pred101L/
rtcoef_rttov11/rttov9pred54L/
rtcoef_rttov11/rttov9pred54L/rtcoef_coms_1_mi.dat
rtcoef_rttov11/rttov9pred54L/rtcoef_envisat_1_atsr-shifted.dat
rtcoef_rttov11/rttov9pred54L/rtcoef_envisat_1_atsr.dat
rtcoef_rttov11/rttov9pred54L/rtcoef_eos_1_aster.dat
rtcoef_rttov11/rttov9pred54L/rtcoef_eos_1_modis-shifted.dat
rtcoef_rttov11/rttov9pred54L/rtcoef_eos_1_modis.dat
rtcoef_rttov11/rttov9pred54L/rtcoef_eos_2_modis-shifted.dat
rtcoef_rttov11/rttov9pred54L/rtcoef_eos_2_modis.dat
rtcoef_rttov11/rttov9pred54L/rtcoef_ers_1_atsr.dat
rtcoef_rttov11/rttov9pred54L/rtcoef_ers_2_atsr.dat
rtcoef_rttov11/rttov9pred54L/rtcoef_fy2_3_vissr.dat
rtcoef_rttov11/rttov9pred54L/rtcoef_fy2_4_vissr.dat
rtcoef_rttov11/rttov9pred54L/rtcoef_fy3_1_iras.dat
rtcoef_rttov11/rttov9pred54L/rtcoef_goes_13_imager.dat
rtcoef_rttov11/rttov9pred54L/rtcoef_goes_14_imager.dat
rtcoef_rttov11/rttov9pred54L/rtcoef_goes_15_imager.dat
rtcoef_rttov11/rttov9pred54L/rtcoef_goes_16_abi.dat
rtcoef_rttov11/rttov9pred54L/rtcoef_himawari_8_ahi.dat
rtcoef_rttov11/rttov9pred54L/rtcoef_jpss_0_viirs.dat
rtcoef_rttov11/rttov9pred54L/rtcoef_landsat_8_oli.dat
rtcoef_rttov11/rttov9pred54L/rtcoef_metop_1_avhrr.dat
rtcoef_rttov11/rttov9pred54L/rtcoef_metop_2_avhrr.dat
rtcoef_rttov11/rttov9pred54L/rtcoef_metopsg_1_metimage.dat
rtcoef_rttov11/rttov9pred54L/rtcoef_msg_1_seviri.dat
rtcoef_rttov11/rttov9pred54L/rtcoef_msg_2_seviri.dat
rtcoef_rttov11/rttov9pred54L/rtcoef_msg_3_seviri.dat
rtcoef_rttov11/rttov9pred54L/rtcoef_msg_4_seviri.dat
rtcoef_rttov11/rttov9pred54L/rtcoef_mtg_1_fci.dat
rtcoef_rttov11/rttov9pred54L/rtcoef_mtsat_1_imager.dat
rtcoef_rttov11/rttov9pred54L/rtcoef_mtsat_2_imager.dat
rtcoef_rttov11/rttov9pred54L/rtcoef_noaa_10_avhrr.dat
rtcoef_rttov11/rttov9pred54L/rtcoef_noaa_11_avhrr.dat
rtcoef_rttov11/rttov9pred54L/rtcoef_noaa_12_avhrr.dat
rtcoef_rttov11/rttov9pred54L/rtcoef_noaa_13_avhrr.dat
rtcoef_rttov11/rttov9pred54L/rtcoef_noaa_14_avhrr.dat
rtcoef_rttov11/rttov9pred54L/rtcoef_noaa_15_avhrr.dat
rtcoef_rttov11/rttov9pred54L/rtcoef_noaa_16_avhrr.dat
rtcoef_rttov11/rttov9pred54L/rtcoef_noaa_17_avhrr.dat
rtcoef_rttov11/rttov9pred54L/rtcoef_noaa_18_avhrr.dat
rtcoef_rttov11/rttov9pred54L/rtcoef_noaa_19_avhrr.dat
rtcoef_rttov11/rttov9pred54L/rtcoef_noaa_5_avhrr.dat
rtcoef_rttov11/rttov9pred54L/rtcoef_noaa_6_avhrr.dat
rtcoef_rttov11/rttov9pred54L/rtcoef_noaa_7_avhrr.dat
rtcoef_rttov11/rttov9pred54L/rtcoef_noaa_8_avhrr.dat
rtcoef_rttov11/rttov9pred54L/rtcoef_noaa_9_avhrr.dat
rtcoef_rttov11/rttov9pred54L/rtcoef_sentinel3_1_slstr.dat
rtcoef_rttov11/rttov_coef_download.sh
rtcoef_rttov11/vtpr.pl
rttov_test/
rttov_test/arch/
rttov_test/profile-datasets-hdf/*
rttov_test/profile-datasets-py/*
rttov_test/profile-datasets/*
rttov_test/rttov_test.pl
rttov_test/rttov_test_plot.py
rttov_test/rttov_test_plot_mod.py
rttov_test/run_example_aer_file_fwd.sh
rttov_test/run_example_aer_param_fwd.sh
rttov_test/run_example_atlas_fwd.sh
rttov_test/run_example_cld_file_fwd.sh
rttov_test/run_example_cld_param_fwd.sh
rttov_test/run_example_fwd.sh
rttov_test/run_example_k.sh
rttov_test/run_example_pc_fwd.sh
rttov_test/run_example_rttovscatt_fwd.sh
rttov_test/test_brdf_atlas.1/
rttov_test/test_brdf_atlas.1/profiles_visnir
rttov_test/test_brdf_atlas.2/*
rttov_test/test_brdf_atlas.sh
rttov_test/test_cnrmmwatlas.sh
rttov_test/test_coef_io.sh
rttov_test/test_coef_io_hdf.sh
rttov_test/test_core.sh
rttov_test/test_cpu.sh
rttov_test/test_emisatlas.1/
rttov_test/test_emisatlas.1/profiles_ir
rttov_test/test_emisatlas.1/profiles_mw
rttov_test/test_emisatlas.2/*
rttov_test/test_example_aer_file_fwd.1/
rttov_test/test_example_aer_file_fwd.1/aer_prof.dat
rttov_test/test_example_aer_file_fwd.1/prof.dat
rttov_test/test_example_aer_file_fwd.2/*
rttov_test/test_example_aer_param_fwd.1/
rttov_test/test_example_aer_param_fwd.1/aer_opt_param_avhrr.dat
rttov_test/test_example_aer_param_fwd.1/prof.dat
rttov_test/test_example_aer_param_fwd.2/*
rttov_test/test_example_atlas_fwd.1/
rttov_test/test_example_atlas_fwd.1/prof.dat
rttov_test/test_example_atlas_fwd.2/*
rttov_test/test_example_cld_file_fwd.1/
rttov_test/test_example_cld_file_fwd.1/cld_prof.dat
rttov_test/test_example_cld_file_fwd.1/prof.dat
rttov_test/test_example_cld_file_fwd.2/*
rttov_test/test_example_cld_param_fwd.1/
rttov_test/test_example_cld_param_fwd.1/cld_opt_param_avhrr.dat
rttov_test/test_example_cld_param_fwd.1/prof.dat
rttov_test/test_example_cld_param_fwd.2/*
rttov_test/test_example_fwd.1/
rttov_test/test_example_fwd.1/prof.dat
rttov_test/test_example_fwd.2/*
rttov_test/test_example_k.1/
rttov_test/test_example_k.1/prof.dat
rttov_test/test_example_k.2/*
rttov_test/test_example_pc_fwd.1/
rttov_test/test_example_pc_fwd.1/prof.dat
rttov_test/test_example_pc_fwd.1/radrec.dat
rttov_test/test_example_pc_fwd.2/*
rttov_test/test_example_rttovscatt_fwd.1/
rttov_test/test_example_rttovscatt_fwd.1/prof.dat
rttov_test/test_example_rttovscatt_fwd.2/*
rttov_test/test_fwd.2/*
rttov_test/test_fwd.sh
rttov_test/test_iratlas.sh
rttov_test/test_multi_instrument.2/*
rttov_test/test_multi_instrument.sh
rttov_test/test_mwatlas.sh
rttov_test/test_pc.2/*
rttov_test/test_pc.sh
rttov_test/test_rttov11.2/*
rttov_test/test_rttov11.sh
rttov_test/test_rttov11_hires.2/*
rttov_test/test_rttov11_hires.sh
rttov_test/test_rttovscatt.1/
rttov_test/test_rttovscatt.1/example_rttovscatt.asc
rttov_test/test_rttovscatt.1/profiles2_fmt
rttov_test/test_rttovscatt.2/*
rttov_test/test_rttovscatt.sh
rttov_test/test_solar.2/*
rttov_test/test_solar.sh
rttov_test/test_zeeman.2/*
rttov_test/test_zeeman.sh
rttov_test/tests.0/*
wrapper/
wrapper/Makefile
wrapper/Options.cpp
wrapper/Options.h
wrapper/Profile.cpp
wrapper/Profile.h
wrapper/Profiles.cpp
wrapper/Profiles.h
wrapper/Rttov.cpp
wrapper/Rttov.h
wrapper/RttovSafe.cpp
wrapper/RttovSafe.h
wrapper/RttovSafe_example.cpp
wrapper/Rttov_common.h
wrapper/Rttov_example.cpp
wrapper/doxygen_config_wrapper
wrapper/example_data.h
wrapper/example_data.py
wrapper/example_data_cpp.h
wrapper/interface_example_c.c
wrapper/interface_example_cpp.cpp
wrapper/interface_example_python.py
wrapper/rttov_c_interface.h
wrapper/rttov_cc_interface.h


James Hocking, NWP SAF
28 September 2015
