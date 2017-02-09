!
! Description
!   User-level subroutines for RTTOV wrapper:
!     rttov_load_inst     - call to initialise an instrument (set options, read coefs)
!     rttov_set_options   - call to change options for an instrument (after reading coefs)
!     rttov_print_options - call to print current options for an instrument
!     rttov_call_direct   - call RTTOV direct
!     rttov_call_k        - call RTTOV K
!     rttov_drop_inst     - deallocation of a single instrument
!     rttov_drop_all      - deallocation of all instruments and atlases
!
!     rttov_ir_emis_atlas_setup   - initialise IR emissivity atlas
!     rttov_mw_emis_atlas_setup   - initialise MW emissivity atlas
!     rttov_brdf_atlas_setup      - initialise BRDF atlas
!     rttov_ir_emis_atlas_dealloc - deallocate IR emissivity atlas
!     rttov_mw_emis_atlas_dealloc - deallocate MW emissivity atlas
!     rttov_brdf_atlas_dealloc    - deallocate BRDF atlas
!
!     rttov_get_*         - call to return members of radiance,
!                           radiance2 and transmission structures
!
!     rttov_get_coef_val_* - call to return variables from RTTOV coefficients
!                            structure; not intended for general use
!
! Copyright:
!    This software was developed within the context of
!    the EUMETSAT Satellite Application Facility on
!    Numerical Weather Prediction (NWP SAF), under the
!    Cooperation Agreement dated 25 November 1998, between
!    EUMETSAT and the Met Office, UK, by one or more partners
!    within the NWP SAF. The partners in the NWP SAF are
!    the Met Office, ECMWF, KNMI and MeteoFrance.
!
!    Copyright 2015, EUMETSAT, All Rights Reserved.

SUBROUTINE rttov_wrapper_f2py
! This subroutine is only here to ensure the Makefile.PL will create dependencies
!INTF_END

  USE parkind1
  USE rttov_types
  USE rttov_const
  USE rttov_wrapper_handle
  USE rttov_wrapper_transfer

END SUBROUTINE rttov_wrapper_f2py


SUBROUTINE rttov_load_inst( &
    inst_id,         &
    opts_str,        &
    nchannels,       &
    channels)
!f2py intent(hide):: nchannels=len(channels)
!
! Set up wrapper:
!   - initialise options
!   - load coefficients
!   - initialise atlases (if requested)
!
#include "throw.h"

  USE parkind1
  USE rttov_wrapper_handle
  USE rttov_getoptions

  IMPLICIT NONE

  INTEGER(jpim),    INTENT(OUT) :: inst_id
  CHARACTER(LEN=*), INTENT(IN)  :: opts_str
  INTEGER(jpim),    INTENT(IN)  :: nchannels
  INTEGER(jpim),    INTENT(IN)  :: channels(nchannels)

  INTEGER(jpim) :: err

#include "rttov_user_options_checkinput.interface"
#include "rttov_errorreport.interface"

!------------------------------------------------------------------------------

  TRY

  DO inst_id = 1, max_ninst
    IF (.NOT. the_rth(inst_id)%init) EXIT
  ENDDO

  IF (inst_id > max_ninst) THEN
    err = errorstatus_fatal
    THROWM(err.NE.0, 'Reached maximum number of instruments: deallocate one and try again')
  ENDIF

  ! Reset values to defaults (in case this inst_id was deallocated and reallocated)
  CALL rttov_wrapper_handle_defaults(inst_id)

  CALL rttov_set_options(err, inst_id, opts_str(1:LEN(opts_str)))
  THROWM(err.NE.0, 'Error setting options')

  CALL rttov_wrapper_handle_load(err, the_rth(inst_id), channels)
  THROW(err.NE.0)

  CALL rttov_wrapper_nullify_structs(the_rth(inst_id))

  the_rth(inst_id)%init = .TRUE.

  IF (the_rth(inst_id)%check_opts) THEN
    CALL rttov_user_options_checkinput(err, the_rth(inst_id)%opts, the_rth(inst_id)%coefs)
    THROW(err.NE.0)
  ENDIF

  IF (the_rth(inst_id)%verbose) THEN
    INFO("Load coefficients: ")
    INFO(the_rth(inst_id)%file_coef)
    IF (TRIM(the_rth(inst_id)%file_scaer) .NE. '') &
      INFO(the_rth(inst_id)%file_scaer)
    IF (TRIM(the_rth(inst_id)%file_sccld) .NE. '') &
      INFO(the_rth(inst_id)%file_sccld)
    IF (TRIM(the_rth(inst_id)%file_pccoef) .NE. '') &
      INFO(the_rth(inst_id)%file_pccoef)
  ENDIF

  CATCH
  the_rth(inst_id)%init = .FALSE.
  inst_id = -1
END SUBROUTINE rttov_load_inst


SUBROUTINE rttov_drop_inst(err, inst_id)
!
! Deallocation of one instrument
!
#include "throw.h"

  USE parkind1
  USE rttov_wrapper_handle

  IMPLICIT NONE
  INTEGER(jpim), INTENT(OUT) :: err
  INTEGER(jpim), INTENT(IN)  :: inst_id

#include "rttov_errorreport.interface"

!------------------------------------------------------------------------------

  TRY

  IF (inst_id < 1 .OR. inst_id > max_ninst) THEN
    err = errorstatus_fatal
    THROWM(err.NE.0, 'inst_id outside valid range')
  ENDIF

  IF (.NOT. the_rth(inst_id)%init) THEN
    IF (the_rth(inst_id)%verbose) INFO('No coefficients loaded')
    RETURN
  ENDIF

  ! Note that handle_drop drops coefficients and atlases
  CALL rttov_wrapper_handle_drop(err, the_rth(inst_id))
  THROWM(err.NE.0, 'An error occured while dropping coefficients')

  IF (the_rth(inst_id)%verbose) THEN
    INFO('Drop coefficients: ')
    INFO(the_rth(inst_id)%file_coef)
    IF (TRIM(the_rth(inst_id)%file_scaer) .NE. '') &
      INFO(the_rth(inst_id)%file_scaer)
    IF (TRIM(the_rth(inst_id)%file_sccld) .NE. '') &
      INFO(the_rth(inst_id)%file_sccld)
    IF (TRIM(the_rth(inst_id)%file_pccoef) .NE. '') &
      INFO(the_rth(inst_id)%file_pccoef)
  ENDIF

  CATCH
END SUBROUTINE rttov_drop_inst


SUBROUTINE rttov_drop_all(err)
!
! Deallocation of all instruments and atlases
!
#include "throw.h"

  USE parkind1
  USE rttov_wrapper_handle

  IMPLICIT NONE
  INTEGER(jpim), INTENT(OUT) :: err

#include "rttov_errorreport.interface"

!------------------------------------------------------------------------------

  TRY

  CALL rttov_ir_emis_atlas_dealloc()
  CALL rttov_mw_emis_atlas_dealloc()
  CALL rttov_brdf_atlas_dealloc()

  CALL rttov_wrapper_handle_drop_all(err)
  THROW(err.NE.0)

  CATCH
END SUBROUTINE rttov_drop_all


SUBROUTINE rttov_ir_emis_atlas_setup(err, path, month, version, inst_id, ang_corr)
!
! Initialise the IR emissivity atlas
!
#include "throw.h"

  USE parkind1
  USE rttov_types, ONLY : rttov_options, rttov_coefs
  USE rttov_wrapper_handle
  USE rttov_const, ONLY : sensor_id_mw, sensor_id_po, sensor_id_ir
  USE mod_iratlas, ONLY : ir_atlas_version
  USE mod_rttov_emis_atlas, ONLY : ir_atlas_init

  IMPLICIT NONE
  INTEGER(jpim),    INTENT(OUT) :: err
  CHARACTER(LEN=*), INTENT(IN)  :: path
  INTEGER(jpim),    INTENT(IN)  :: month
  INTEGER(jpim),    INTENT(IN)  :: version
  INTEGER(jpim),    INTENT(IN)  :: inst_id
  INTEGER(jpim),    INTENT(IN)  :: ang_corr

  LOGICAL(jplm)       :: do_ang_corr
  INTEGER(jpim)       :: use_version
  TYPE(rttov_options) :: opts
  TYPE(rttov_coefs)   :: coefs

#include "rttov_setup_emis_atlas.interface"
#include "rttov_errorreport.interface"

!------------------------------------------------------------------------------

  TRY

  IF (ir_atlas_init) THEN
    err = errorstatus_fatal
    THROWM(err.NE.0, 'IR emissivity atlas already initialised')
  ENDIF

  IF (month < 1 .OR. month > 12) THEN
    err = errorstatus_fatal
    THROWM(err.NE.0, 'Invalid month')
  ENDIF

  IF (inst_id > max_ninst) THEN
    err = errorstatus_fatal
    THROWM(err.NE.0, 'inst_id outside valid range')
  ENDIF

  IF (version > 0) THEN
    use_version = version
  ELSE
    use_version = ir_atlas_version
  ENDIF

  do_ang_corr = (ang_corr .NE. 0)

  IF (inst_id > 0) THEN
    IF (.NOT. the_rth(inst_id)%init) THEN
      err = errorstatus_fatal
      THROWM(err.NE.0, 'inst_id is not valid: must be zero or the ID for an initialised instrument')
    ENDIF
    IF (the_rth(inst_id)%coefs%coef%id_sensor == sensor_id_mw .OR. &
        the_rth(inst_id)%coefs%coef%id_sensor == sensor_id_po) THEN
      err = errorstatus_fatal
      THROWM(err.NE.0, 'Specified instrument is not an IR sensor')
    ENDIF

    ! Initialise atlas for the specific instrument inst_id
    CALL rttov_setup_emis_atlas(                     &
           err,                                      &
           the_rth(inst_id)%opts,                    &
           month,                                    &
           the_rth(inst_id)%coefs,                   &
           path(1:LEN(path)),                        &
           ir_atlas_read_std          = .FALSE.,     &
           ir_atlas_single_instrument = .TRUE.,      &
           ir_atlas_ang_corr          = do_ang_corr, &
           ir_atlas_ver               = use_version)
  ELSE
    ! Initialise atlas for multiple instruments using dummy opts and coefs
    coefs%coef%id_sensor = sensor_id_ir
    opts%config%verbose = .FALSE.
    CALL rttov_setup_emis_atlas(                     &
           err,                                      &
           opts,                                     &
           month,                                    &
           coefs,                                    &
           path(1:LEN(path)),                        &
           ir_atlas_read_std          = .FALSE.,     &
           ir_atlas_single_instrument = .FALSE.,     &
           ir_atlas_ang_corr          = do_ang_corr, &
           ir_atlas_ver               = use_version)
  ENDIF
  THROW(err.NE.0)

  CATCH
END SUBROUTINE rttov_ir_emis_atlas_setup


SUBROUTINE rttov_mw_emis_atlas_setup(err, path, month, version, inst_id)
!
! Initialise the MW emissivity atlas
!
#include "throw.h"

  USE parkind1
  USE rttov_types, ONLY : rttov_options, rttov_coefs
  USE rttov_wrapper_handle
  USE rttov_const, ONLY : sensor_id_ir, sensor_id_hi, sensor_id_mw
  USE mod_mwatlas, ONLY : mw_atlas_version
  USE mod_cnrm_mw_atlas, ONLY : cnrm_mw_atlas_version => mw_atlas_version
  USE mod_rttov_emis_atlas, ONLY : mw_atlas_init

  IMPLICIT NONE
  INTEGER(jpim),    INTENT(OUT) :: err
  CHARACTER(LEN=*), INTENT(IN)  :: path
  INTEGER(jpim),    INTENT(IN)  :: month
  INTEGER(jpim),    INTENT(IN)  :: version
  INTEGER(jpim),    INTENT(IN)  :: inst_id

  INTEGER(jpim)       :: use_version
  TYPE(rttov_options) :: opts
  TYPE(rttov_coefs)   :: coefs

#include "rttov_setup_emis_atlas.interface"
#include "rttov_errorreport.interface"

!------------------------------------------------------------------------------

  TRY

  IF (mw_atlas_init) THEN
    err = errorstatus_fatal
    THROWM(err.NE.0, 'MW emissivity atlas already initialised')
  ENDIF

  IF (month < 1 .OR. month > 12) THEN
    err = errorstatus_fatal
    THROWM(err.NE.0, 'Invalid month')
  ENDIF

  IF (inst_id > max_ninst) THEN
    err = errorstatus_fatal
    THROWM(err.NE.0, 'inst_id outside valid range')
  ENDIF

  IF (version > 0) THEN
    use_version = version
  ELSE
    use_version = mw_atlas_version
  ENDIF

  IF (inst_id > 0) THEN
    IF (.NOT. the_rth(inst_id)%init) THEN
      err = errorstatus_fatal
      THROWM(err.NE.0, 'inst_id is not valid: must be zero or the ID for an initialised instrument')
    ENDIF
    IF (the_rth(inst_id)%coefs%coef%id_sensor == sensor_id_ir .OR. &
        the_rth(inst_id)%coefs%coef%id_sensor == sensor_id_hi) THEN
      err = errorstatus_fatal
      THROWM(err.NE.0, 'Specified instrument is not a MW sensor')
    ENDIF

    ! Initialise atlas using the specific instrument inst_id
    ! This is required for CNRM atlas; makes no difference for TELSEM
    CALL rttov_setup_emis_atlas(   &
           err,                    &
           the_rth(inst_id)%opts,  &
           month,                  &
           the_rth(inst_id)%coefs, &
           path(1:LEN(path)),      &
           mw_atlas_ver = use_version)
  ELSE

    IF (use_version == cnrm_mw_atlas_version) THEN
      err = errorstatus_fatal
      THROWM(err.NE.0, 'inst_id for loaded instrument is required to initialise CNRM MW atlas')
    ENDIF

    ! Initialise TELSEM atlas using dummy coefs and opts
    coefs%coef%id_sensor = sensor_id_mw
    opts%config%verbose = .FALSE.
    CALL rttov_setup_emis_atlas(   &
           err,                    &
           opts,                   &
           month,                  &
           coefs,                  &
           path(1:LEN(path)),      &
           mw_atlas_ver = use_version)
  ENDIF
  THROW(err.NE.0)

  CATCH
END SUBROUTINE rttov_mw_emis_atlas_setup


SUBROUTINE rttov_brdf_atlas_setup(err, path, month, version, inst_id)
!
! Initialise the BRDF atlas
!
#include "throw.h"

  USE parkind1
  USE rttov_types, ONLY : rttov_options, rttov_coefs
  USE rttov_wrapper_handle
  USE rttov_const, ONLY : sensor_id_mw, sensor_id_po, sensor_id_ir
  USE mod_brdf_atlas, ONLY : vn_atlas_version
  USE mod_rttov_brdf_atlas, ONLY : vn_atlas_init

  IMPLICIT NONE
  INTEGER(jpim),    INTENT(OUT) :: err
  CHARACTER(LEN=*), INTENT(IN)  :: path
  INTEGER(jpim),    INTENT(IN)  :: month
  INTEGER(jpim),    INTENT(IN)  :: version
  INTEGER(jpim),    INTENT(IN)  :: inst_id

  INTEGER(jpim)       :: use_version
  TYPE(rttov_options) :: opts
  TYPE(rttov_coefs)   :: coefs

#include "rttov_setup_brdf_atlas.interface"
#include "rttov_errorreport.interface"

!------------------------------------------------------------------------------

  TRY

  IF (vn_atlas_init) THEN
    err = errorstatus_fatal
    THROWM(err.NE.0, 'BRDF atlas already initialised')
  ENDIF

  IF (month < 1 .OR. month > 12) THEN
    err = errorstatus_fatal
    THROWM(err.NE.0, 'Invalid month')
  ENDIF

  IF (inst_id > max_ninst) THEN
    err = errorstatus_fatal
    THROWM(err.NE.0, 'inst_id outside valid range')
  ENDIF

  IF (version > 0) THEN
    use_version = version
  ELSE
    use_version = vn_atlas_version
  ENDIF

  IF (inst_id > 0) THEN
    IF (.NOT. the_rth(inst_id)%init) THEN
      err = errorstatus_fatal
      THROWM(err.NE.0, 'inst_id is not valid: must be zero or the ID for an initialised instrument')
    ENDIF
    IF (the_rth(inst_id)%coefs%coef%id_sensor == sensor_id_mw .OR. &
        the_rth(inst_id)%coefs%coef%id_sensor == sensor_id_po) THEN
      err = errorstatus_fatal
      THROWM(err.NE.0, 'Specified instrument is not a VIS/IR sensor')
    ENDIF

    ! Initialise atlas for the specific instrument inst_id
    CALL rttov_setup_brdf_atlas(                     &
           err,                                      &
           the_rth(inst_id)%opts,                    &
           month,                                    &
           the_rth(inst_id)%coefs,                   &
           path(1:LEN(path)),                        &
           brdf_atlas_single_instrument = .TRUE.,    &
           vn_atlas_ver                 = use_version)
  ELSE
    ! Initialise atlas for multiple instruments using dummy opts and coefs
    coefs%coef%id_sensor = sensor_id_ir
    opts%config%verbose = .FALSE.
    CALL rttov_setup_brdf_atlas(                     &
           err,                                      &
           opts,                                     &
           month,                                    &
           coefs,                                    &
           path(1:LEN(path)),                        &
           brdf_atlas_single_instrument = .FALSE.,   &
           vn_atlas_ver                 = use_version)
  ENDIF
  THROW(err.NE.0)

  CATCH
END SUBROUTINE rttov_brdf_atlas_setup



SUBROUTINE rttov_ir_emis_atlas_dealloc()
!
! Deallocate the IR emissivity atlas
!
  USE rttov_types, ONLY : rttov_coefs
  USE rttov_const, ONLY : sensor_id_ir

  IMPLICIT NONE

  TYPE(rttov_coefs) :: coefs

#include "rttov_deallocate_emis_atlas.interface"

!------------------------------------------------------------------------------

  coefs%coef%id_sensor = sensor_id_ir
  CALL rttov_deallocate_emis_atlas(coefs)

END SUBROUTINE rttov_ir_emis_atlas_dealloc

SUBROUTINE rttov_mw_emis_atlas_dealloc()
!
! Deallocate the MW emissivity atlas
!
  USE rttov_types, ONLY : rttov_coefs
  USE rttov_const, ONLY : sensor_id_mw

  IMPLICIT NONE

  TYPE(rttov_coefs) :: coefs

#include "rttov_deallocate_emis_atlas.interface"

!------------------------------------------------------------------------------

  coefs%coef%id_sensor = sensor_id_mw
  CALL rttov_deallocate_emis_atlas(coefs)

END SUBROUTINE rttov_mw_emis_atlas_dealloc

SUBROUTINE rttov_brdf_atlas_dealloc()
!
! Deallocate the BRDF atlas
!
  USE rttov_types, ONLY : rttov_coefs
  USE rttov_const, ONLY : sensor_id_ir

  IMPLICIT NONE

  TYPE(rttov_coefs) :: coefs

#include "rttov_deallocate_brdf_atlas.interface"

!------------------------------------------------------------------------------

  coefs%coef%id_sensor = sensor_id_ir
  CALL rttov_deallocate_brdf_atlas(coefs)

END SUBROUTINE rttov_brdf_atlas_dealloc


SUBROUTINE rttov_call_direct( &
    err,               &
    inst_id,           &
    channel_list,      &
    datetimes,         &    ! profile dates/times                                     (6,nprofiles)
    angles,            &    ! satzen, satazi, sunzen, sunazi angles                   (4,nprofiles)
    surfgeom,          &    ! lat, lon, elevation                                     (3,nprofiles)
    surftype,          &    ! surftype, watertype                                     (2,nprofiles)
    skin,              &    ! skin T, salinity, snow_frac, foam_frac, fastem_coefsx5  (9,nprofiles)
    s2m,               &    ! 2m p, 2m t, 2m q, 10m wind u, v, wind-fetch             (6,nprofiles)
    simplecloud,       &    ! ctp, cfraction                                          (2,nprofiles)
    icecloud,          &    ! ish, idg                                                (2,nprofiles)
    zeeman,            &    ! Be, cosbk                                               (2,nprofiles)
    p,                 &    ! pressure                                                (nlevels,nprofiles)
    t,                 &    ! temperature                                             (nlevels,nprofiles)
    gas_units,         &    ! units for gas profiles
    gas_id,            &    ! gas ID list                                             (ngases)
    gases,             &    ! gas profiles                                            (nlevels,nprofiles,ngases)
    surfemisrefl,      &    ! surface emissivities and BRDFs                          (nchannels,nprofiles,2)
    btrefl,            &    ! output BTs/refls (for thermal/solar chans)              (nchannels,nprofiles)
    rads,              &    ! output radiances                                        (nchannels,nprofiles)
    nchannels, ngases, nlevels, nprofiles)
!f2py threadsafe
!f2py intent(hide):: nchannels=len(channel_list)
!f2py intent(hide):: ngases=len(gas_id)
!f2py intent(hide):: nlevels=shape(p, 1)
!f2py intent(hide):: nprofiles=shape(p, 2)

!
! Prepares input profiles, calls RTTOV and returns radiances
!

#include "throw.h"

  USE rttov_wrapper_handle

  USE rttov_wrapper_transfer

  USE parkind1, ONLY : jpim, jprb, jplm

  USE rttov_types, ONLY : &
      rttov_chanprof,     &
      rttov_emissivity,   &
      rttov_reflectance,  &
      profile_type,       &
      radiance_type,      &
      transmission_type

  USE rttov_const, ONLY : &
      surftype_sea,       &
      sensor_id_ir,       &
      sensor_id_hi,       &
      sensor_id_mw,       &
      sensor_id_po

  USE mod_rttov_emis_atlas, ONLY : ir_atlas_init, mw_atlas_init
  USE mod_rttov_brdf_atlas, ONLY : vn_atlas_init

  IMPLICIT NONE

  INTEGER(jpim),    INTENT(OUT)   :: err
  INTEGER(jpim),    INTENT(IN)    :: nchannels, ngases, nlevels, nprofiles

  INTEGER(jpim),    INTENT(IN)    :: inst_id
  INTEGER(jpim),    INTENT(IN)    :: channel_list(nchannels)
  INTEGER(jpim),    INTENT(IN)    :: datetimes(6,nprofiles)
  REAL(jprb),       INTENT(IN)    :: angles(4,nprofiles)
  REAL(jprb),       INTENT(IN)    :: surfgeom(3,nprofiles)
  INTEGER(jpim),    INTENT(IN)    :: surftype(2,nprofiles)
  REAL(jprb),       INTENT(IN)    :: skin(9,nprofiles)
  REAL(jprb),       INTENT(IN)    :: s2m(6,nprofiles)
  REAL(jprb),       INTENT(IN)    :: simplecloud(2,nprofiles)
  INTEGER(jpim),    INTENT(IN)    :: icecloud(2,nprofiles)
  REAL(jprb),       INTENT(IN)    :: zeeman(2,nprofiles)
  REAL(jprb),       INTENT(IN)    :: p(nlevels,nprofiles)
  REAL(jprb),       INTENT(IN)    :: t(nlevels,nprofiles)
  INTEGER(jpim),    INTENT(IN)    :: gas_units
  INTEGER(jpim),    INTENT(IN)    :: gas_id(ngases)
  REAL(jprb),       INTENT(IN)    :: gases(nlevels,nprofiles,ngases)

  REAL(jprb),       INTENT(INOUT) :: surfemisrefl(nchannels,nprofiles,2)

  REAL(jprb),       INTENT(INOUT) :: btrefl(nchannels,nprofiles)
  REAL(jprb),       INTENT(INOUT) :: rads(nchannels,nprofiles)


  INTEGER(jpim)      :: asw
  INTEGER(jpim)      :: i, k, iprof, iprof1, lochan, hichan
  INTEGER(jpim)      :: nchanprof, thisnchanprof, nprof, thisnprof, ncalls
  LOGICAL(jplm)      :: use_emis_atlas, use_brdf_atlas
  CHARACTER(LEN=256) :: msg

  TYPE(rttov_chanprof),    POINTER :: chanprof(:)    => NULL()
  LOGICAL(jplm),           POINTER :: calcemis(:)    => NULL()
  TYPE(rttov_emissivity),  POINTER :: emissivity(:)  => NULL()
  LOGICAL(jplm),           POINTER :: calcrefl(:)    => NULL()
  TYPE(rttov_reflectance), POINTER :: reflectance(:) => NULL()
  TYPE(profile_type),      POINTER :: profiles(:)    => NULL()
  TYPE(transmission_type)          :: transmission
  TYPE(radiance_type)              :: radiance
  TYPE(radiance2_type)             :: radiance2

#include "rttov_alloc_direct.interface"
#include "rttov_get_emis.interface"
#include "rttov_get_brdf.interface"
#include "rttov_direct.interface"
#include "rttov_parallel_direct.interface"
#include "rttov_errorreport.interface"

!------------------------------------------------------------------------------

  TRY

  CALL rttov_wrapper_handle_check(err, inst_id)
  THROW(err.NE.0)

  ! Deallocate (if necessary) and re-allocate the output structures
  CALL rttov_wrapper_handle_alloc(err, the_rth(inst_id), nprofiles, nchannels, nlevels, 0_jpim)
  THROW(err.NE.0)
  CALL rttov_wrapper_handle_alloc(err, the_rth(inst_id), nprofiles, nchannels, nlevels, 1_jpim)
  THROW(err.NE.0)

  ! nchannels     - input argument, number of channels simulated for every profile
  ! nprofiles     - input argument, total number of profiles passed to this subroutine
  ! ncalls        - number of calls to RTTOV
  ! nprof         - number of profiles per call to RTTOV
  ! thisnprof     - number of profiles in current call (may differ to nprof on last call)
  ! nchanprof     - size of chanprof per call to RTTOV
  ! thisnchanprof - size of chanprof for current call (may differ to nchanprof on last call)

  nprof = MIN(nprofiles, the_rth(inst_id)%nprofs_per_call)
  nchanprof = nprof * nchannels

  ! Allocate RTTOV data
  asw = 1_jpim
  CALL rttov_alloc_direct( &
              err,                         &
              asw,                         &
              nprof,                       &
              nchanprof,                   &
              nlevels,                     &
              chanprof,                    &
              the_rth(inst_id)%opts,       &
              profiles,                    &
              the_rth(inst_id)%coefs,      &
              transmission,                &
              radiance,                    &
              radiancedata2 = radiance2,   &
              calcemis      = calcemis,    &
              emissivity    = emissivity,  &
              calcrefl      = calcrefl,    &
              reflectance   = reflectance, &
              init=.TRUE._jplm)
  THROWM(err.NE.0, 'Error allocating RTTOV structures')

  ! Populate input arrays which don't change
  DO i = 1, nprof
    lochan = (i - 1) * nchannels + 1
    hichan = i * nchannels
    chanprof(lochan:hichan)%prof = i
    chanprof(lochan:hichan)%chan = channel_list(:)
  ENDDO

  ! Simple cloud scheme cloud top albedo cannot currently be set by wrapper
  reflectance%refl_cloud_top = 0._jprb

  ! Flags to indicate whether atlases should be used
  use_emis_atlas = (ir_atlas_init .AND. &
                     (the_rth(inst_id)%coefs%coef%id_sensor == sensor_id_ir .OR. &
                      the_rth(inst_id)%coefs%coef%id_sensor == sensor_id_hi)) &
                   .OR. &
                   (mw_atlas_init .AND. &
                     (the_rth(inst_id)%coefs%coef%id_sensor == sensor_id_mw .OR. &
                      the_rth(inst_id)%coefs%coef%id_sensor == sensor_id_po))

  use_brdf_atlas = the_rth(inst_id)%opts%rt_ir%addsolar .AND. vn_atlas_init .AND. &
                     (the_rth(inst_id)%coefs%coef%id_sensor == sensor_id_ir .OR. &
                      the_rth(inst_id)%coefs%coef%id_sensor == sensor_id_hi)

  IF (the_rth(inst_id)%verbose) THEN
    WRITE(msg, '(a,i4,a,i8)') &
      'Running RTTOV using nthreads = ', the_rth(inst_id)%nthreads, &
      ' and nprofs_per_call = ', the_rth(inst_id)%nprofs_per_call
    INFO(msg)
  ENDIF

  ! Main loop over batches of profiles

  ncalls = nprofiles / the_rth(inst_id)%nprofs_per_call
  IF (MOD(nprofiles, the_rth(inst_id)%nprofs_per_call) /= 0) ncalls = ncalls + 1
  thisnchanprof = nchanprof ! This may be overwritten in the final loop below

  DO i = 1, ncalls

    ! Calculate low/high profile and chanprof indices
    iprof1 = (i - 1) * the_rth(inst_id)%nprofs_per_call + 1
    IF (iprof1 + the_rth(inst_id)%nprofs_per_call > nprofiles) THEN
      thisnprof = nprofiles - iprof1 + 1
      thisnchanprof = thisnprof * nchannels
    ELSE
      thisnprof = the_rth(inst_id)%nprofs_per_call
    ENDIF

    ! Populate the profile structure
    CALL rttov_copy_to_profiles(  &
        the_rth(inst_id),         &
        profiles(1:thisnprof),    &
        iprof1,                   &
        datetimes,                &
        angles, surfgeom,         &
        surftype, skin, s2m,      &
        simplecloud, icecloud,    &
        zeeman,                   &
        p, t, gas_units, gas_id, gases)


    ! Set surface emissivities and BRDFs: default to calcemis/calcrefl TRUE
    calcemis(:) = .TRUE.
    calcrefl(:) = .TRUE.

    ! If atlas was initialised and there are any non-sea profiles then call the atlas

    ! Emissivities
    IF (use_emis_atlas .AND. ANY(profiles(:)%skin%surftype /= surftype_sea)) THEN
      CALL rttov_get_emis(              &
                err,                    &
                the_rth(inst_id)%opts,  &
                chanprof,               &
                profiles,               &
                the_rth(inst_id)%coefs, &
                emissivity = emissivity(:)%emis_in)
      THROWM(err.NE.0, 'Error reading emissivity atlas')

      ! Set calcemis TRUE where atlas has no data
      calcemis(:) = (emissivity(:)%emis_in <= 0._jprb)

      ! Ensure calcemis is TRUE over sea
      DO iprof = 1, thisnprof
        IF (profiles(iprof)%skin%surftype == surftype_sea) &
          calcemis((iprof-1)*nchannels+1:iprof*nchannels) = .TRUE.
      ENDDO
    ENDIF

    ! BRDFs
    IF (use_brdf_atlas .AND. ANY(profiles(:)%skin%surftype /= surftype_sea)) THEN
      CALL rttov_get_brdf(              &
                err,                    &
                chanprof,               &
                profiles,               &
                the_rth(inst_id)%coefs, &
                brdf = reflectance(:)%refl_in)
      THROWM(err.NE.0, 'Error reading BRDF atlas')

      ! Set calcrefl TRUE where atlas has no data
      calcrefl(:) = (reflectance(:)%refl_in < 0._jprb)

      ! Ensure calcrefl is TRUE over sea
      DO iprof = 1, thisnprof
        IF (profiles(iprof)%skin%surftype == surftype_sea) &
          calcrefl((iprof-1)*nchannels+1:iprof*nchannels) = .TRUE.
      ENDDO
    ENDIF

    ! If user passed any non-negative emissivities/BRDFs then use these
    ! and set calcemis/calcrefl FALSE
    DO iprof = 1, thisnprof
      k = iprof1 + iprof - 1
      lochan = (iprof - 1) * nchannels + 1
      hichan = iprof * nchannels
      WHERE (surfemisrefl(:,k,1) >= 0._jprb)
        emissivity(lochan:hichan)%emis_in = surfemisrefl(:,k,1)
        calcemis(lochan:hichan) = .FALSE.
      ENDWHERE
      WHERE (surfemisrefl(:,k,2) >= 0._jprb)
        reflectance(lochan:hichan)%refl_in = surfemisrefl(:,k,2)
        calcrefl(lochan:hichan) = .FALSE.
      ENDWHERE
    ENDDO

    ! Run RTTOV
    IF (the_rth(inst_id)%nthreads <= 1) THEN
      CALL rttov_direct( &
          err,                                        &
          chanprof(1:thisnchanprof),                  &
          the_rth(inst_id)%opts,                      &
          profiles(1:thisnprof),                      &
          the_rth(inst_id)%coefs,                     &
          transmission,                               &
          radiance,                                   &
          radiance2,                                  &
          calcemis    = calcemis(1:thisnchanprof),    &
          emissivity  = emissivity(1:thisnchanprof),  &
          calcrefl    = calcrefl(1:thisnchanprof),    &
          reflectance = reflectance(1:thisnchanprof))
    ELSE
      CALL rttov_parallel_direct( &
          err,                                        &
          chanprof(1:thisnchanprof),                  &
          the_rth(inst_id)%opts,                      &
          profiles(1:thisnprof),                      &
          the_rth(inst_id)%coefs,                     &
          transmission,                               &
          radiance,                                   &
          radiance2,                                  &
          calcemis    = calcemis(1:thisnchanprof),    &
          emissivity  = emissivity(1:thisnchanprof),  &
          calcrefl    = calcrefl(1:thisnchanprof),    &
          reflectance = reflectance(1:thisnchanprof), &
          nthreads    = the_rth(inst_id)%nthreads)
    ENDIF
    THROWM(err.NE.0, 'Error running RTTOV')


    ! Store results
    DO iprof = 1, thisnprof
      k = iprof1 + iprof - 1
      lochan = (iprof - 1) * nchannels + 1
      hichan = iprof * nchannels
      WHERE (the_rth(inst_id)%coefs%coef%ss_val_chn(channel_list(:)) < 2)
        btrefl(:,k) = radiance%bt(lochan:hichan)
      ELSEWHERE
        btrefl(:,k) = radiance%refl(lochan:hichan)
      ENDWHERE
      rads(:,k) = radiance%total(lochan:hichan)
      surfemisrefl(:,k,1) = emissivity(lochan:hichan)%emis_out
      surfemisrefl(:,k,2) = reflectance(lochan:hichan)%refl_out
    ENDDO

    lochan = (iprof1 - 1) * nchannels + 1
    hichan = (iprof1 + thisnprof - 1) * nchannels

    IF (the_rth(inst_id)%store_trans) THEN
      the_rth(inst_id)%transmission%tau_total(lochan:hichan)             = transmission%tau_total(1:thisnchanprof)
      the_rth(inst_id)%transmission%tau_levels(:,lochan:hichan)          = transmission%tau_levels(:,1:thisnchanprof)
      the_rth(inst_id)%transmission%tausun_total_path2(lochan:hichan)    = transmission%tausun_total_path2(1:thisnchanprof)
      the_rth(inst_id)%transmission%tausun_levels_path2(:,lochan:hichan) = transmission%tausun_levels_path2(:,1:thisnchanprof)
      the_rth(inst_id)%transmission%tausun_total_path1(lochan:hichan)    = transmission%tausun_total_path1(1:thisnchanprof)
      the_rth(inst_id)%transmission%tausun_levels_path1(:,lochan:hichan) = transmission%tausun_levels_path1(:,1:thisnchanprof)
    ENDIF

    IF (the_rth(inst_id)%store_rad) THEN
      the_rth(inst_id)%radiance%clear(lochan:hichan)      = radiance%clear(1:thisnchanprof)
      the_rth(inst_id)%radiance%total(lochan:hichan)      = radiance%total(1:thisnchanprof)
      the_rth(inst_id)%radiance%bt_clear(lochan:hichan)   = radiance%bt_clear(1:thisnchanprof)
      the_rth(inst_id)%radiance%bt(lochan:hichan)         = radiance%bt(1:thisnchanprof)
      the_rth(inst_id)%radiance%refl_clear(lochan:hichan) = radiance%refl_clear(1:thisnchanprof)
      the_rth(inst_id)%radiance%refl(lochan:hichan)       = radiance%refl(1:thisnchanprof)
      the_rth(inst_id)%radiance%cloudy(lochan:hichan)     = radiance%cloudy(1:thisnchanprof)
      the_rth(inst_id)%radiance%overcast(:,lochan:hichan) = radiance%overcast(:,1:thisnchanprof)
    ENDIF

    IF (the_rth(inst_id)%store_rad2) THEN
      the_rth(inst_id)%radiance2%upclear(lochan:hichan)     = radiance2%upclear(1:thisnchanprof)
      the_rth(inst_id)%radiance2%dnclear(lochan:hichan)     = radiance2%dnclear(1:thisnchanprof)
      the_rth(inst_id)%radiance2%refldnclear(lochan:hichan) = radiance2%refldnclear(1:thisnchanprof)
      the_rth(inst_id)%radiance2%up(:,lochan:hichan)        = radiance2%up(:,1:thisnchanprof)
      the_rth(inst_id)%radiance2%down(:,lochan:hichan)      = radiance2%down(:,1:thisnchanprof)
      the_rth(inst_id)%radiance2%surf(:,lochan:hichan)      = radiance2%surf(:,1:thisnchanprof)
    ENDIF

  ENDDO ! loop over calls to RTTOV

  ! Deallocate RTTOV data
  asw = 0_jpim
  CALL rttov_alloc_direct( &
              err,                         &
              asw,                         &
              nprof,                       &
              nchanprof,                   &
              nlevels,                     &
              chanprof,                    &
              the_rth(inst_id)%opts,       &
              profiles,                    &
              the_rth(inst_id)%coefs,      &
              transmission,                &
              radiance,                    &
              radiancedata2 = radiance2,   &
              calcemis      = calcemis,    &
              emissivity    = emissivity,  &
              calcrefl      = calcrefl,    &
              reflectance   = reflectance)
  THROWM(err.NE.0, 'Error deallocating RTTOV structures')

  CATCH
END SUBROUTINE rttov_call_direct


SUBROUTINE rttov_call_k( &
    err,               &
    inst_id,           &
    channel_list,      &
    datetimes,         &    ! profile dates/times                                     (6,nprofiles)
    angles,            &    ! satzen, satazi, sunzen, sunazi angles                   (4,nprofiles)
    surfgeom,          &    ! lat, lon, elevation                                     (3,nprofiles)
    surftype,          &    ! surftype, watertype                                     (2,nprofiles)
    skin,              &    ! skin T, salinity, snow_frac, foam_frac, fastem_coefsx5  (9,nprofiles)
    skin_k,            &    ! skin K                                                  (9,nchannels,nprofiles)
    s2m,               &    ! 2m p, 2m t, 2m q, 10m wind u, v, wind-fetch             (6,nprofiles)
    s2m_k,             &    ! 2m K                                                    (6,nchannels,nprofiles)
    simplecloud,       &    ! ctp, cfraction                                          (2,nprofiles)
    simplecloud_k,     &    ! ctp, cfraction K                                        (2,nchannels,nprofiles)
    icecloud,          &    ! ish, idg                                                (2,nprofiles)
    zeeman,            &    ! Be, cosbk                                               (2,nprofiles)
    p,                 &    ! pressure                                                (nlevels,nprofiles)
    p_k,               &    ! output pressure K                                       (nlevels,nchannels,nprofiles)
    t,                 &    ! temperature                                             (nlevels,nprofiles)
    t_k,               &    ! output temperature K                                    (nlevels,nchannels,nprofiles)
    gas_units,         &    ! units for gas profiles
    gas_id,            &    ! gas ID list                                             (ngases)
    gases,             &    ! gas profiles                                            (nlevels,nprofiles,ngases)
    gases_k,           &    ! output gas profiles K                                   (nlevels,nchannels,nprofiles,ngases)
    surfemisrefl,      &    ! surface emissivities and BRDFs                          (nchannels,nprofiles,2)
    surfemisrefl_k,    &    ! output surface emissivities and BRDFs K                 (nchannels,nprofiles,2)
    btrefl,            &    ! output BTs/refls (for thermal/solar chans)              (nchannels,nprofiles)
    rads,              &    ! output radiances                                        (nchannels,nprofiles)
    bt_k,              &    ! input BT perturbations (thermal chans only)             (nchannels,nprofiles)
    rads_k,            &    ! input radiance perturbations                            (nchannels,nprofiles)
    nchannels, ngases, nlevels, nprofiles)
!f2py threadsafe
!f2py intent(hide):: nchannels=len(channel_list)
!f2py intent(hide):: ngases=len(gas_id)
!f2py intent(hide):: nlevels=shape(p, 1)
!f2py intent(hide):: nprofiles=shape(p, 2)

!
! Prepares input profiles, calls RTTOV K and returns radiances and profiles_k
!

#include "throw.h"

  USE rttov_wrapper_handle

  USE rttov_wrapper_transfer

  USE parkind1, ONLY : jpim, jprb, jplm

  USE rttov_types, ONLY : &
      rttov_chanprof,     &
      rttov_emissivity,   &
      rttov_reflectance,  &
      profile_type,       &
      radiance_type,      &
      transmission_type

  USE rttov_const, ONLY : &
      surftype_sea,       &
      sensor_id_ir,       &
      sensor_id_hi,       &
      sensor_id_mw,       &
      sensor_id_po

  USE mod_rttov_emis_atlas, ONLY : ir_atlas_init, mw_atlas_init
  USE mod_rttov_brdf_atlas, ONLY : vn_atlas_init

  IMPLICIT NONE

  INTEGER(jpim),    INTENT(OUT)   :: err

  INTEGER(jpim),    INTENT(IN)    :: nchannels, ngases, nlevels, nprofiles

  INTEGER(jpim),    INTENT(IN)    :: inst_id
  INTEGER(jpim),    INTENT(IN)    :: channel_list(nchannels)
  INTEGER(jpim),    INTENT(IN)    :: datetimes(6,nprofiles)
  REAL(jprb),       INTENT(IN)    :: angles(4,nprofiles)
  REAL(jprb),       INTENT(IN)    :: surfgeom(3,nprofiles)
  INTEGER(jpim),    INTENT(IN)    :: surftype(2,nprofiles)
  REAL(jprb),       INTENT(IN)    :: skin(9,nprofiles)
  REAL(jprb),       INTENT(INOUT) :: skin_k(9,nchannels,nprofiles)
  REAL(jprb),       INTENT(IN)    :: s2m(6,nprofiles)
  REAL(jprb),       INTENT(INOUT) :: s2m_k(6,nchannels,nprofiles)
  REAL(jprb),       INTENT(IN)    :: simplecloud(2,nprofiles)
  REAL(jprb),       INTENT(INOUT) :: simplecloud_k(2,nchannels,nprofiles)
  INTEGER(jpim),    INTENT(IN)    :: icecloud(2,nprofiles)
  REAL(jprb),       INTENT(IN)    :: zeeman(2,nprofiles)
  REAL(jprb),       INTENT(IN)    :: p(nlevels,nprofiles)
  REAL(jprb),       INTENT(INOUT) :: p_k(nlevels,nchannels,nprofiles)
  REAL(jprb),       INTENT(IN)    :: t(nlevels,nprofiles)
  REAL(jprb),       INTENT(INOUT) :: t_k(nlevels,nchannels,nprofiles)
  INTEGER(jpim),    INTENT(IN)    :: gas_units
  INTEGER(jpim),    INTENT(IN)    :: gas_id(ngases)
  REAL(jprb),       INTENT(IN)    :: gases(nlevels,nprofiles,ngases)
  REAL(jprb),       INTENT(INOUT) :: gases_k(nlevels,nchannels,nprofiles,ngases)

  REAL(jprb),       INTENT(INOUT) :: surfemisrefl(nchannels,nprofiles,2)
  REAL(jprb),       INTENT(INOUT) :: surfemisrefl_k(nchannels,nprofiles,2)

  REAL(jprb),       INTENT(INOUT) :: btrefl(nchannels,nprofiles)
  REAL(jprb),       INTENT(INOUT) :: rads(nchannels,nprofiles)
  REAL(jprb),       INTENT(IN)    :: bt_k(nchannels,nprofiles)
  REAL(jprb),       INTENT(IN)    :: rads_k(nchannels,nprofiles)


  INTEGER(jpim)      :: asw
  INTEGER(jpim)      :: i, j, k, iprof, iprof1, lochan, hichan
  INTEGER(jpim)      :: nchanprof, thisnchanprof, nprof, thisnprof, ncalls
  LOGICAL(jplm)      :: use_emis_atlas, use_brdf_atlas
  CHARACTER(LEN=128) :: msg

  TYPE(rttov_chanprof),    POINTER :: chanprof(:)      => NULL()
  LOGICAL(jplm),           POINTER :: calcemis(:)      => NULL()
  TYPE(rttov_emissivity),  POINTER :: emissivity(:)    => NULL()
  TYPE(rttov_emissivity),  POINTER :: emissivity_k(:)  => NULL()
  LOGICAL(jplm),           POINTER :: calcrefl(:)      => NULL()
  TYPE(rttov_reflectance), POINTER :: reflectance(:)   => NULL()
  TYPE(rttov_reflectance), POINTER :: reflectance_k(:) => NULL()
  TYPE(profile_type),      POINTER :: profiles(:)      => NULL()
  TYPE(profile_type),      POINTER :: profiles_k(:)    => NULL()
  TYPE(transmission_type)          :: transmission
  TYPE(transmission_type)          :: transmission_k
  TYPE(radiance_type)              :: radiance
  TYPE(radiance_type)              :: radiance_k

#include "rttov_alloc_k.interface"
#include "rttov_get_emis.interface"
#include "rttov_get_brdf.interface"
#include "rttov_init_rad.interface"
#include "rttov_init_transmission.interface"
#include "rttov_init_prof.interface"
#include "rttov_k.interface"
#include "rttov_parallel_k.interface"
#include "rttov_errorreport.interface"

!------------------------------------------------------------------------------

  TRY

  CALL rttov_wrapper_handle_check(err, inst_id)
  THROW(err.NE.0)

  ! Deallocate (if necessary) and re-allocate the output structures
  CALL rttov_wrapper_handle_alloc(err, the_rth(inst_id), nprofiles, nchannels, nlevels, 0_jpim)
  THROW(err.NE.0)
  CALL rttov_wrapper_handle_alloc(err, the_rth(inst_id), nprofiles, nchannels, nlevels, 1_jpim)
  THROW(err.NE.0)

  ! nchannels     - input argument, number of channels simulated for every profile
  ! nprofiles     - input argument, total number of profiles passed to this subroutine
  ! ncalls        - number of calls to RTTOV
  ! nprof         - number of profiles per call to RTTOV
  ! thisnprof     - number of profiles in current call (may differ to nprof on last call)
  ! iprof1        - index of first profile in current batch
  ! nchanprof     - size of chanprof per call to RTTOV
  ! thisnchanprof - size of chanprof for current call (may differ to nchanprof on last call)

  nprof = MIN(nprofiles, the_rth(inst_id)%nprofs_per_call)
  nchanprof = nprof * nchannels

  ! Allocate RTTOV data
  asw = 1_jpim
  CALL rttov_alloc_k( &
              err,                           &
              asw,                           &
              nprof,                         &
              nchanprof,                     &
              nlevels,                       &
              chanprof,                      &
              the_rth(inst_id)%opts,         &
              profiles,                      &
              profiles_k,                    &
              the_rth(inst_id)%coefs,        &
              transmission,                  &
              transmission_k,                &
              radiance,                      &
              radiance_k,                    &
              calcemis      = calcemis,      &
              emissivity    = emissivity,    &
              emissivity_k  = emissivity_k,  &
              calcrefl      = calcrefl,      &
              reflectance   = reflectance,   &
              reflectance_k = reflectance_k, &
              init=.TRUE._jplm)
  THROWM(err.NE.0, 'Error allocating RTTOV structures')

  ! Populate input arrays which don't change
  DO j = 1, nprof
    lochan = (j - 1) * nchannels + 1
    hichan = j * nchannels
    chanprof(lochan:hichan)%prof = j
    chanprof(lochan:hichan)%chan = channel_list(:)
  ENDDO

  ! Simple cloud scheme cloud top albedo cannot currently be set by wrapper
  reflectance%refl_cloud_top = 0._jprb

  ! Flags to indicate whether atlases should be used
  use_emis_atlas = (ir_atlas_init .AND. &
                     (the_rth(inst_id)%coefs%coef%id_sensor == sensor_id_ir .OR. &
                      the_rth(inst_id)%coefs%coef%id_sensor == sensor_id_hi)) &
                   .OR. &
                   (mw_atlas_init .AND. &
                     (the_rth(inst_id)%coefs%coef%id_sensor == sensor_id_mw .OR. &
                      the_rth(inst_id)%coefs%coef%id_sensor == sensor_id_po))

  use_brdf_atlas = the_rth(inst_id)%opts%rt_ir%addsolar .AND. vn_atlas_init .AND. &
                     (the_rth(inst_id)%coefs%coef%id_sensor == sensor_id_ir .OR. &
                      the_rth(inst_id)%coefs%coef%id_sensor == sensor_id_hi)

  IF (the_rth(inst_id)%verbose) THEN
    WRITE(msg, '(a,i4,a,i8)') &
      'Running RTTOV K using nthreads = ', the_rth(inst_id)%nthreads, &
      ' and nprofs_per_call = ', the_rth(inst_id)%nprofs_per_call
    INFO(msg)
  ENDIF

  ! Main loop over batches of profiles

  ncalls = nprofiles / the_rth(inst_id)%nprofs_per_call
  IF (MOD(nprofiles, the_rth(inst_id)%nprofs_per_call) /= 0) ncalls = ncalls + 1
  thisnchanprof = nchanprof
  thisnprof = the_rth(inst_id)%nprofs_per_call

  DO i = 1, ncalls

    ! Calculate first profile index of this batch
    iprof1 = (i - 1) * the_rth(inst_id)%nprofs_per_call + 1

    ! Last loop may have fewer profiles
    IF (iprof1 + the_rth(inst_id)%nprofs_per_call > nprofiles) THEN
      thisnprof = nprofiles - iprof1 + 1
      thisnchanprof = thisnprof * nchannels
    ENDIF

    ! Populate the profile structure
    CALL rttov_copy_to_profiles(  &
        the_rth(inst_id),         &
        profiles(1:thisnprof),    &
        iprof1,                   &
        datetimes,                &
        angles, surfgeom,         &
        surftype, skin, s2m,      &
        simplecloud, icecloud,    &
        zeeman,                   &
        p, t, gas_units, gas_id, gases)


    ! Set surface emissivities and BRDFs: default to calcemis/calcrefl TRUE
    calcemis(:) = .TRUE.
    calcrefl(:) = .TRUE.

    ! If atlas was initialised and there are any non-sea profiles then call the atlas

    ! Emissivities
    IF (use_emis_atlas .AND. ANY(profiles(:)%skin%surftype /= surftype_sea)) THEN
      CALL rttov_get_emis(                 &
                err,                       &
                the_rth(inst_id)%opts,     &
                chanprof(1:thisnchanprof), &
                profiles(1:thisnprof),     &
                the_rth(inst_id)%coefs,    &
                emissivity = emissivity(:)%emis_in)
      THROWM(err.NE.0, 'Error reading emissivity atlas')

      ! Set calcemis TRUE where atlas has no data
      calcemis(:) = (emissivity(:)%emis_in <= 0._jprb)

      ! Ensure calcemis is TRUE over sea
      DO iprof = 1, thisnprof
        IF (profiles(iprof)%skin%surftype == surftype_sea) &
          calcemis((iprof-1)*nchannels+1:iprof*nchannels) = .TRUE.
      ENDDO
    ENDIF

    ! BRDFs
    IF (use_brdf_atlas .AND. ANY(profiles(:)%skin%surftype /= surftype_sea)) THEN
      CALL rttov_get_brdf(                 &
                err,                       &
                chanprof(1:thisnchanprof), &
                profiles(1:thisnprof),     &
                the_rth(inst_id)%coefs,    &
                brdf = reflectance(:)%refl_in)
      THROWM(err.NE.0, 'Error reading BRDF atlas')

      ! Set calcrefl TRUE where atlas has no data
      calcrefl(:) = (reflectance(:)%refl_in < 0._jprb)

      ! Ensure calcrefl is TRUE over sea
      DO iprof = 1, thisnprof
        IF (profiles(iprof)%skin%surftype == surftype_sea) &
          calcrefl((iprof-1)*nchannels+1:iprof*nchannels) = .TRUE.
      ENDDO
    ENDIF

    ! If user passed any non-negative emissivities/BRDFs then use these
    ! and set calcemis/calcrefl FALSE
    DO iprof = 1, thisnprof
      k = iprof1 + iprof - 1
      lochan = (iprof - 1) * nchannels + 1
      hichan = iprof * nchannels
      WHERE (surfemisrefl(:,k,1) >= 0._jprb)
        emissivity(lochan:hichan)%emis_in = surfemisrefl(:,k,1)
        calcemis(lochan:hichan) = .FALSE.
      ENDWHERE
      WHERE (surfemisrefl(:,k,2) >= 0._jprb)
        reflectance(lochan:hichan)%refl_in = surfemisrefl(:,k,2)
        calcrefl(lochan:hichan) = .FALSE.
      ENDWHERE
    ENDDO

    ! Populate the radiance_k structure
    CALL rttov_init_rad(radiance_k)
    CALL rttov_copy_to_radiance_k( &
        radiance_k,                &
        iprof1,                    &
        thisnprof,                 &
        bt_k, rads_k)

    ! Initialise other K arrays
    CALL rttov_init_prof(profiles_k)
    CALL rttov_init_transmission(transmission_k)
    emissivity_k(:)%emis_in = 0._jprb
    reflectance_k(:)%refl_in = 0._jprb
    emissivity_k(:)%emis_out = 0._jprb
    reflectance_k(:)%refl_out = 0._jprb

    ! Run RTTOV K
    IF (the_rth(inst_id)%nthreads <= 1) THEN
      CALL rttov_k( &
          err,                                            &
          chanprof(1:thisnchanprof),                      &
          the_rth(inst_id)%opts,                          &
          profiles(1:thisnprof),                          &
          profiles_k(1:thisnchanprof),                    &
          the_rth(inst_id)%coefs,                         &
          transmission,                                   &
          transmission_k,                                 &
          radiance,                                       &
          radiance_k,                                     &
          calcemis      = calcemis(1:thisnchanprof),      &
          emissivity    = emissivity(1:thisnchanprof),    &
          emissivity_k  = emissivity_k(1:thisnchanprof),  &
          calcrefl      = calcrefl(1:thisnchanprof),      &
          reflectance   = reflectance(1:thisnchanprof),   &
          reflectance_k = reflectance_k(1:thisnchanprof))
    ELSE
      CALL rttov_parallel_k( &
          err,                                            &
          chanprof(1:thisnchanprof),                      &
          the_rth(inst_id)%opts,                          &
          profiles(1:thisnprof),                          &
          profiles_k(1:thisnchanprof),                    &
          the_rth(inst_id)%coefs,                         &
          transmission,                                   &
          transmission_k,                                 &
          radiance,                                       &
          radiance_k,                                     &
          calcemis      = calcemis(1:thisnchanprof),      &
          emissivity    = emissivity(1:thisnchanprof),    &
          emissivity_k  = emissivity_k(1:thisnchanprof),  &
          calcrefl      = calcrefl(1:thisnchanprof),      &
          reflectance   = reflectance(1:thisnchanprof),   &
          reflectance_k = reflectance_k(1:thisnchanprof), &
          nthreads      = the_rth(inst_id)%nthreads)
    ENDIF
    THROWM(err.NE.0, 'Error running RTTOV')

    ! Store results
    DO iprof = 1, thisnprof
      k = iprof1 + iprof - 1
      lochan = (iprof - 1) * nchannels + 1
      hichan = iprof * nchannels
      WHERE (the_rth(inst_id)%coefs%coef%ss_val_chn(channel_list(:)) < 2)
        btrefl(:,k) = radiance%bt(lochan:hichan)
      ELSEWHERE
        btrefl(:,k) = radiance%refl(lochan:hichan)
      ENDWHERE
      rads(:,k) = radiance%total(lochan:hichan)
      surfemisrefl(:,k,1) = emissivity(lochan:hichan)%emis_out
      surfemisrefl(:,k,2) = reflectance(lochan:hichan)%refl_out
      surfemisrefl_k(:,k,1) = emissivity_k(lochan:hichan)%emis_in
      surfemisrefl_k(:,k,2) = reflectance_k(lochan:hichan)%refl_in
    ENDDO
    CALL rttov_copy_from_profiles_k( &
        the_rth(inst_id),            &
        profiles_k(1:thisnchanprof), &
        iprof1,                      &
        skin_k, s2m_k,               &
        simplecloud_k,               &
        p_k, t_k, gas_id, gases_k)

    lochan = (iprof1 - 1) * nchannels + 1
    hichan = (iprof1 + thisnprof - 1) * nchannels

    IF (the_rth(inst_id)%store_trans) THEN
      the_rth(inst_id)%transmission%tau_total(lochan:hichan)             = transmission%tau_total(1:thisnchanprof)
      the_rth(inst_id)%transmission%tau_levels(:,lochan:hichan)          = transmission%tau_levels(:,1:thisnchanprof)
      the_rth(inst_id)%transmission%tausun_total_path2(lochan:hichan)    = transmission%tausun_total_path2(1:thisnchanprof)
      the_rth(inst_id)%transmission%tausun_levels_path2(:,lochan:hichan) = transmission%tausun_levels_path2(:,1:thisnchanprof)
      the_rth(inst_id)%transmission%tausun_total_path1(lochan:hichan)    = transmission%tausun_total_path1(1:thisnchanprof)
      the_rth(inst_id)%transmission%tausun_levels_path1(:,lochan:hichan) = transmission%tausun_levels_path1(:,1:thisnchanprof)
    ENDIF

    IF (the_rth(inst_id)%store_rad) THEN
      the_rth(inst_id)%radiance%clear(lochan:hichan)      = radiance%clear(1:thisnchanprof)
      the_rth(inst_id)%radiance%total(lochan:hichan)      = radiance%total(1:thisnchanprof)
      the_rth(inst_id)%radiance%bt_clear(lochan:hichan)   = radiance%bt_clear(1:thisnchanprof)
      the_rth(inst_id)%radiance%bt(lochan:hichan)         = radiance%bt(1:thisnchanprof)
      the_rth(inst_id)%radiance%refl_clear(lochan:hichan) = radiance%refl_clear(1:thisnchanprof)
      the_rth(inst_id)%radiance%refl(lochan:hichan)       = radiance%refl(1:thisnchanprof)
      the_rth(inst_id)%radiance%cloudy(lochan:hichan)     = radiance%cloudy(1:thisnchanprof)
      the_rth(inst_id)%radiance%overcast(:,lochan:hichan) = radiance%overcast(:,1:thisnchanprof)
    ENDIF
  ENDDO ! loop over calls to RTTOV

  ! Deallocate RTTOV data
  asw = 0_jpim
  CALL rttov_alloc_k( &
              err,                           &
              asw,                           &
              nprof,                         &
              nchanprof,                     &
              nlevels,                       &
              chanprof,                      &
              the_rth(inst_id)%opts,         &
              profiles,                      &
              profiles_k,                    &
              the_rth(inst_id)%coefs,        &
              transmission,                  &
              transmission_k,                &
              radiance,                      &
              radiance_k,                    &
              calcemis      = calcemis,      &
              emissivity    = emissivity,    &
              emissivity_k  = emissivity_k,  &
              calcrefl      = calcrefl,      &
              reflectance   = reflectance,   &
              reflectance_k = reflectance_k)
  THROWM(err.NE.0, 'Error deallocating RTTOV structures')

  CATCH
END SUBROUTINE rttov_call_k


SUBROUTINE rttov_set_options(err, inst_id, opts_str)
!
! Set RTTOV options (once instrument has been initialised
!   with a call to rttov_load_inst)
!
#include "throw.h"

  USE rttov_wrapper_handle
  USE rttov_getoptions
  USE parkind1, ONLY : jprb, jpim

  IMPLICIT NONE

  INTEGER(jpim),    INTENT(OUT) :: err
  INTEGER(jpim),    INTENT(IN)  :: inst_id
  CHARACTER(LEN=*), INTENT(IN)  :: opts_str

  INTEGER(jpim)       :: ival
  REAL(jprb)          :: rval
  CHARACTER(LEN=256)  :: strval
  LOGICAL(jplm)       :: exists

#include "rttov_user_options_checkinput.interface"
#include "rttov_errorreport.interface"

!------------------------------------------------------------------------------

  TRY

  IF (inst_id < 1 .OR. inst_id > max_ninst) THEN
    err = errorstatus_fatal
    THROWM(err.NE.0, 'inst_id outside valid range')
  ENDIF

  CALL initoptionsstr(opts_str(1:LEN(opts_str)))


  ! General config options

  CALL getoption('opts%config%apply_reg_limits', ival, exists = exists)
  IF (exists) the_rth(inst_id)%opts%config%apply_reg_limits = (ival .NE. 0)

  CALL getoption('opts%config%verbose', ival, exists = exists)
  IF (exists) the_rth(inst_id)%opts%config%verbose = (ival .NE. 0)

  CALL getoption('opts%config%do_checkinput', ival, exists = exists)
  IF (exists) the_rth(inst_id)%opts%config%do_checkinput = (ival .NE. 0)


  ! Interpolation options

  CALL getoption('opts%interpolation%addinterp', ival, exists = exists)
  IF (exists) the_rth(inst_id)%opts%interpolation%addinterp = (ival .NE. 0)

  CALL getoption('opts%interpolation%interp_mode', ival, exists = exists)
  IF (exists) the_rth(inst_id)%opts%interpolation%interp_mode = ival

  CALL getoption('opts%interpolation%lgradp', ival, exists = exists)
  IF (exists) the_rth(inst_id)%opts%interpolation%lgradp = (ival .NE. 0)

  CALL getoption('opts%interpolation%spacetop', ival, exists = exists)
  IF (exists) the_rth(inst_id)%opts%interpolation%spacetop = (ival .NE. 0)

  CALL getoption('opts%interpolation%reg_limit_extrap', ival, exists = exists)
  IF (exists) the_rth(inst_id)%opts%interpolation%reg_limit_extrap = (ival .NE. 0)


  ! General RT options

  CALL getoption('opts%rt_all%addrefrac', ival, exists = exists)
  IF (exists) the_rth(inst_id)%opts%rt_all%addrefrac = (ival .NE. 0)

  CALL getoption('opts%rt_all%switchrad', ival, exists = exists)
  IF (exists) the_rth(inst_id)%opts%rt_all%switchrad = (ival .NE. 0)

  CALL getoption('opts%rt_all%use_q2m', ival, exists = exists)
  IF (exists) the_rth(inst_id)%opts%rt_all%use_q2m = (ival .NE. 0)

  CALL getoption('opts%rt_all%do_lambertian', ival, exists = exists)
  IF (exists) the_rth(inst_id)%opts%rt_all%do_lambertian = (ival .NE. 0)


  ! MW RT options

  CALL getoption('opts%rt_mw%fastem_version', ival, exists = exists)
  IF (exists) the_rth(inst_id)%opts%rt_mw%fastem_version = ival

  CALL getoption('opts%rt_mw%supply_foam_fraction', ival, exists = exists)
  IF (exists) the_rth(inst_id)%opts%rt_mw%supply_foam_fraction = (ival .NE. 0)

  CALL getoption('opts%rt_mw%clw_data', ival, exists = exists)
  IF (exists) the_rth(inst_id)%opts%rt_mw%clw_data = (ival .NE. 0)

  CALL getoption('opts%rt_mw%do_lambertian', ival, exists = exists)
  IF (exists) the_rth(inst_id)%opts%rt_mw%do_lambertian = (ival .NE. 0)


  ! VIS/IR RT options

  CALL getoption('opts%rt_ir%ozone_data', ival, exists = exists)
  IF (exists) the_rth(inst_id)%opts%rt_ir%ozone_data = (ival .NE. 0)

  CALL getoption('opts%rt_ir%co2_data', ival, exists = exists)
  IF (exists) the_rth(inst_id)%opts%rt_ir%co2_data = (ival .NE. 0)

  CALL getoption('opts%rt_ir%n2o_data', ival, exists = exists)
  IF (exists) the_rth(inst_id)%opts%rt_ir%n2o_data = (ival .NE. 0)

  CALL getoption('opts%rt_ir%co_data', ival, exists = exists)
  IF (exists) the_rth(inst_id)%opts%rt_ir%co_data = (ival .NE. 0)

  CALL getoption('opts%rt_ir%ch4_data', ival, exists = exists)
  IF (exists) the_rth(inst_id)%opts%rt_ir%ch4_data = (ival .NE. 0)


  CALL getoption('opts%rt_ir%addsolar', ival, exists = exists)
  IF (exists) the_rth(inst_id)%opts%rt_ir%addsolar = (ival .NE. 0)

  CALL getoption('opts%rt_ir%do_nlte_correction', ival, exists = exists)
  IF (exists) the_rth(inst_id)%opts%rt_ir%do_nlte_correction = (ival .NE. 0)

  CALL getoption('opts%rt_ir%addclouds', ival, exists = exists)
  IF (exists) the_rth(inst_id)%opts%rt_ir%addclouds = (ival .NE. 0)

  CALL getoption('opts%rt_ir%addaerosl', ival, exists = exists)
  IF (exists) the_rth(inst_id)%opts%rt_ir%addaerosl = (ival .NE. 0)

  CALL getoption('opts%rt_ir%user_aer_opt_param', ival, exists = exists)
  IF (exists) the_rth(inst_id)%opts%rt_ir%user_aer_opt_param = (ival .NE. 0)

  CALL getoption('opts%rt_ir%user_cld_opt_param', ival, exists = exists)
  IF (exists) the_rth(inst_id)%opts%rt_ir%user_cld_opt_param = (ival .NE. 0)

  CALL getoption('opts%rt_ir%cldstr_threshold', rval, exists = exists)
  IF (exists) the_rth(inst_id)%opts%rt_ir%cldstr_threshold = rval

  CALL getoption('opts%rt_ir%cldstr_simple', ival, exists = exists)
  IF (exists) the_rth(inst_id)%opts%rt_ir%cldstr_simple = (ival .NE. 0)

  CALL getoption('opts%rt_ir%do_lambertian', ival, exists = exists)
  IF (exists) the_rth(inst_id)%opts%rt_ir%do_lambertian = (ival .NE. 0)


  ! PC options

  CALL getoption('opts%rt_ir%pc%addpc', ival, exists = exists)
  IF (exists) the_rth(inst_id)%opts%rt_ir%pc%addpc = (ival .NE. 0)

  CALL getoption('opts%rt_ir%pc%addradrec', ival, exists = exists)
  IF (exists) the_rth(inst_id)%opts%rt_ir%pc%addradrec = (ival .NE. 0)

  CALL getoption('opts%rt_ir%pc%ipcbnd', ival, exists = exists)
  IF (exists) the_rth(inst_id)%opts%rt_ir%pc%ipcbnd = ival

  CALL getoption('opts%rt_ir%pc%ipcreg', ival, exists = exists)
  IF (exists) the_rth(inst_id)%opts%rt_ir%pc%ipcreg = ival


  ! Wrapper options

  CALL getoption('nthreads', ival, exists = exists)
  IF (exists) the_rth(inst_id)%nthreads = ival

  CALL getoption('nprofs_per_call', ival, exists = exists)
  IF (exists) the_rth(inst_id)%nprofs_per_call = ival
  IF (the_rth(inst_id)%nprofs_per_call < 1) the_rth(inst_id)%nprofs_per_call = 1

  CALL getoption('verbose_wrapper', ival, exists = exists)
  IF (exists) the_rth(inst_id)%verbose = (ival .NE. 0)

  CALL getoption('check_opts', ival, exists = exists)
  IF (exists) the_rth(inst_id)%check_opts = (ival .NE. 0)

  CALL getoption('store_trans', ival, exists = exists)
  IF (exists) the_rth(inst_id)%store_trans = (ival .NE. 0)

  CALL getoption('store_rad', ival, exists = exists)
  IF (exists) the_rth(inst_id)%store_rad = (ival .NE. 0)

  CALL getoption('store_rad2', ival, exists = exists)
  IF (exists) the_rth(inst_id)%store_rad2 = (ival .NE. 0)

  ! Radiance structure must be allocated if radiance2 is allocated
  IF (the_rth(inst_id)%store_rad2) the_rth(inst_id)%store_rad = .TRUE.


  CALL getoption('file_coef', strval, exists = exists)
  IF (exists) the_rth(inst_id)%file_coef = TRIM(strval)
  CALL getoption('file_scaer', strval, exists = exists)
  IF (exists) the_rth(inst_id)%file_scaer = TRIM(strval)
  CALL getoption('file_sccld', strval, exists = exists)
  IF (exists) the_rth(inst_id)%file_sccld = TRIM(strval)
!   CALL getoption('file_pccoef', strval, exists = exists)
!   IF (exists) the_rth(inst_id)%file_pccoef = TRIM(strval)
!   CALL getoption('file_mietable', strval, exists = exists)
!   IF (exists) the_rth(inst_id)%file_mietable = TRIM(strval)


  CALL checkoptions()

  ! When this is called from rttov_load_inst, the coefs haven't been loaded (init = FALSE)
  ! so we shouldn't call rttov_user_options_checkinput at this point: only call it when
  ! user is calling rttov_set_options after initialising an intrument.
  IF (the_rth(inst_id)%init .AND. the_rth(inst_id)%check_opts) THEN
    CALL rttov_user_options_checkinput(err, the_rth(inst_id)%opts, the_rth(inst_id)%coefs)
    THROW(err.NE.0)
  ENDIF

  CATCH
END SUBROUTINE rttov_set_options

SUBROUTINE rttov_print_options(err, inst_id)
!
! Print current RTTOV options for an instrument
!
#include "throw.h"

  USE rttov_wrapper_handle
  USE parkind1, ONLY : jpim
  USE rttov_global, ONLY : error_unit

  IMPLICIT NONE

  INTEGER(jpim),    INTENT(OUT) :: err
  INTEGER(jpim),    INTENT(IN)  :: inst_id

  CHARACTER(LEN=20)   :: tmp_text ! temporary string for formatting output

#include "rttov_print_opts.interface"
#include "rttov_errorreport.interface"

!------------------------------------------------------------------------------

  TRY

  CALL rttov_wrapper_handle_check(err, inst_id)
  THROW(err.NE.0)

  CALL rttov_print_opts(the_rth(inst_id)%opts)

  WRITE(error_unit,'(a)')           "Wrapper options"
  WRITE(error_unit,'(2x,a,a)')      "Coef file            ", TRIM(the_rth(inst_id)%file_coef)
  IF (the_rth(inst_id)%file_scaer .NE. '') &
    WRITE(error_unit,'(2x,a,a)')      "Aerosol coef file    ", TRIM(the_rth(inst_id)%file_scaer)
  IF (the_rth(inst_id)%file_sccld .NE. '') &
    WRITE(error_unit,'(2x,a,a)')      "Cloud coef file      ", TRIM(the_rth(inst_id)%file_sccld)
!   IF (the_rth(inst_id)%file_pccoef .NE. '') &
!     WRITE(error_unit,'(2x,a,a)')      "PC-RTTOV coef file   ", the_rth(inst_id)%file_pccoef
!   IF (the_rth(inst_id)%file_mietable .NE. '') &
!     WRITE(error_unit,'(2x,a,a)')      "RTTOV-SCATT mietable ", the_rth(inst_id)%file_mietable
  WRITE(tmp_text,'(i6)') the_rth(inst_id)%nprofs_per_call
  WRITE(error_unit,'(2x,a,a)')      "nprofs_per_call      ", TRIM(ADJUSTL(tmp_text))
  WRITE(tmp_text,'(i5)') the_rth(inst_id)%nthreads
  WRITE(error_unit,'(2x,a,a)')      "nthreads             ", TRIM(ADJUSTL(tmp_text))
  WRITE(error_unit,'(2x,a,l1)')     "verbose_wrapper      ", the_rth(inst_id)%verbose
  WRITE(error_unit,'(2x,a,l1)')     "check_opts           ", the_rth(inst_id)%check_opts
  WRITE(error_unit,'(2x,a,l1)')     "store_trans          ", the_rth(inst_id)%store_trans
  WRITE(error_unit,'(2x,a,l1)')     "store_rad            ", the_rth(inst_id)%store_rad
  WRITE(error_unit,'(2x,a,l1)')     "store_rad2           ", the_rth(inst_id)%store_rad2
  WRITE(error_unit,'(a)',advance='yes')

  CATCH
END SUBROUTINE rttov_print_options


! The following subroutines return RTTOV output arrays: they can be called
! after RTTOV has been called.

SUBROUTINE rttov_get_rad_clear(err, inst_id, rad_clear, nchanprof)
!f2py intent(hide):: nchanprof=len(rad_clear)
#include "throw.h"
  USE rttov_wrapper_handle
  USE parkind1, ONLY : jpim, jprb

  IMPLICIT NONE

  INTEGER(jpim), INTENT(OUT)   :: err
  INTEGER(jpim), INTENT(IN)    :: inst_id
  INTEGER(jpim), INTENT(IN)    :: nchanprof
  REAL(jprb),    INTENT(INOUT) :: rad_clear(nchanprof)

#include "rttov_errorreport.interface"

  TRY
  CALL rttov_wrapper_handle_check(err, inst_id)
  THROW(err.NE.0)

  IF (.NOT. ASSOCIATED(the_rth(inst_id)%radiance%clear)) THEN
    err = errorstatus_fatal
    THROWM(err.NE.0, 'Wrapper radiance not allocated: has RTTOV been called?')
  ENDIF

  rad_clear(:) = the_rth(inst_id)%radiance%clear(:)
  CATCH
END SUBROUTINE rttov_get_rad_clear

SUBROUTINE rttov_get_rad_total(err, inst_id, rad_total, nchanprof)
!f2py intent(hide):: nchanprof=len(rad_total)
#include "throw.h"
  USE rttov_wrapper_handle
  USE parkind1, ONLY : jpim, jprb

  IMPLICIT NONE

  INTEGER(jpim), INTENT(OUT)   :: err
  INTEGER(jpim), INTENT(IN)    :: inst_id
  INTEGER(jpim), INTENT(IN)    :: nchanprof
  REAL(jprb),    INTENT(INOUT) :: rad_total(nchanprof)

#include "rttov_errorreport.interface"

  TRY
  CALL rttov_wrapper_handle_check(err, inst_id)
  THROW(err.NE.0)

  IF (.NOT. ASSOCIATED(the_rth(inst_id)%radiance%total)) THEN
    err = errorstatus_fatal
    THROWM(err.NE.0, 'Wrapper radiance not allocated: has RTTOV been called?')
  ENDIF

  rad_total(:) = the_rth(inst_id)%radiance%total(:)
  CATCH
END SUBROUTINE rttov_get_rad_total

SUBROUTINE rttov_get_rad_cloudy(err, inst_id, rad_cloudy, nchanprof)
!f2py intent(hide):: nchanprof=len(rad_cloudy)
#include "throw.h"
  USE rttov_wrapper_handle
  USE parkind1, ONLY : jpim, jprb

  IMPLICIT NONE

  INTEGER(jpim), INTENT(OUT)   :: err
  INTEGER(jpim), INTENT(IN)    :: inst_id
  INTEGER(jpim), INTENT(IN)    :: nchanprof
  REAL(jprb),    INTENT(INOUT) :: rad_cloudy(nchanprof)

#include "rttov_errorreport.interface"

  TRY
  CALL rttov_wrapper_handle_check(err, inst_id)
  THROW(err.NE.0)

  IF (.NOT. ASSOCIATED(the_rth(inst_id)%radiance%cloudy)) THEN
    err = errorstatus_fatal
    THROWM(err.NE.0, 'Wrapper radiance not allocated: has RTTOV been called?')
  ENDIF

  rad_cloudy(:) = the_rth(inst_id)%radiance%cloudy(:)
  CATCH
END SUBROUTINE rttov_get_rad_cloudy

SUBROUTINE rttov_get_bt_clear(err, inst_id, bt_clear, nchanprof)
!f2py intent(hide):: nchanprof=len(bt_clear)
#include "throw.h"
  USE rttov_wrapper_handle
  USE parkind1, ONLY : jpim, jprb

  IMPLICIT NONE

  INTEGER(jpim), INTENT(OUT)   :: err
  INTEGER(jpim), INTENT(IN)    :: inst_id
  INTEGER(jpim), INTENT(IN)    :: nchanprof
  REAL(jprb),    INTENT(INOUT) :: bt_clear(nchanprof)

#include "rttov_errorreport.interface"

  TRY
  CALL rttov_wrapper_handle_check(err, inst_id)
  THROW(err.NE.0)

  IF (.NOT. ASSOCIATED(the_rth(inst_id)%radiance%bt_clear)) THEN
    err = errorstatus_fatal
    THROWM(err.NE.0, 'Wrapper radiance not allocated: has RTTOV been called?')
  ENDIF

  bt_clear(:) = the_rth(inst_id)%radiance%bt_clear(:)
  CATCH
END SUBROUTINE rttov_get_bt_clear

SUBROUTINE rttov_get_bt(err, inst_id, bt, nchanprof)
!f2py intent(hide):: nchanprof=len(bt)
#include "throw.h"
  USE rttov_wrapper_handle
  USE parkind1, ONLY : jpim, jprb

  IMPLICIT NONE

  INTEGER(jpim), INTENT(OUT)   :: err
  INTEGER(jpim), INTENT(IN)    :: inst_id
  INTEGER(jpim), INTENT(IN)    :: nchanprof
  REAL(jprb),    INTENT(INOUT) :: bt(nchanprof)

#include "rttov_errorreport.interface"

  TRY
  CALL rttov_wrapper_handle_check(err, inst_id)
  THROW(err.NE.0)

  IF (.NOT. ASSOCIATED(the_rth(inst_id)%radiance%bt)) THEN
    err = errorstatus_fatal
    THROWM(err.NE.0, 'Wrapper radiance not allocated: has RTTOV been called?')
  ENDIF

  bt(:) = the_rth(inst_id)%radiance%bt(:)
  CATCH
END SUBROUTINE rttov_get_bt

SUBROUTINE rttov_get_refl_clear(err, inst_id, refl_clear, nchanprof)
!f2py intent(hide):: nchanprof=len(refl_clear)
#include "throw.h"
  USE rttov_wrapper_handle
  USE parkind1, ONLY : jpim, jprb

  IMPLICIT NONE

  INTEGER(jpim), INTENT(OUT)   :: err
  INTEGER(jpim), INTENT(IN)    :: inst_id
  INTEGER(jpim), INTENT(IN)    :: nchanprof
  REAL(jprb),    INTENT(INOUT) :: refl_clear(nchanprof)

#include "rttov_errorreport.interface"

  TRY
  CALL rttov_wrapper_handle_check(err, inst_id)
  THROW(err.NE.0)

  IF (.NOT. ASSOCIATED(the_rth(inst_id)%radiance%refl_clear)) THEN
    err = errorstatus_fatal
    THROWM(err.NE.0, 'Wrapper radiance not allocated: has RTTOV been called?')
  ENDIF

  refl_clear(:) = the_rth(inst_id)%radiance%refl_clear(:)
  CATCH
END SUBROUTINE rttov_get_refl_clear

SUBROUTINE rttov_get_refl(err, inst_id, refl, nchanprof)
!f2py intent(hide):: nchanprof=len(refl)
#include "throw.h"
  USE rttov_wrapper_handle
  USE parkind1, ONLY : jpim, jprb

  IMPLICIT NONE

  INTEGER(jpim), INTENT(OUT)   :: err
  INTEGER(jpim), INTENT(IN)    :: inst_id
  INTEGER(jpim), INTENT(IN)    :: nchanprof
  REAL(jprb),    INTENT(INOUT) :: refl(nchanprof)

#include "rttov_errorreport.interface"

  TRY
  CALL rttov_wrapper_handle_check(err, inst_id)
  THROW(err.NE.0)

  IF (.NOT. ASSOCIATED(the_rth(inst_id)%radiance%refl)) THEN
    err = errorstatus_fatal
    THROWM(err.NE.0, 'Wrapper radiance not allocated: has RTTOV been called?')
  ENDIF

  refl(:) = the_rth(inst_id)%radiance%refl(:)
  CATCH
END SUBROUTINE rttov_get_refl

SUBROUTINE rttov_get_overcast(err, inst_id, overcast, nchanprof, nlayers)
!f2py intent(hide):: nchanprof=len(overcast)
!f2py intent(hide):: nlayers=len(overcast[0])
#include "throw.h"
  USE rttov_wrapper_handle
  USE parkind1, ONLY : jpim, jprb

  IMPLICIT NONE

  INTEGER(jpim), INTENT(OUT)   :: err
  INTEGER(jpim), INTENT(IN)    :: inst_id
  INTEGER(jpim), INTENT(IN)    :: nchanprof
  INTEGER(jpim), INTENT(IN)    :: nlayers
  REAL(jprb),    INTENT(INOUT) :: overcast(nlayers,nchanprof)

#include "rttov_errorreport.interface"

  TRY
  CALL rttov_wrapper_handle_check(err, inst_id)
  THROW(err.NE.0)

  IF (.NOT. ASSOCIATED(the_rth(inst_id)%radiance%overcast)) THEN
    err = errorstatus_fatal
    THROWM(err.NE.0, 'Wrapper radiance not allocated: has RTTOV been called?')
  ENDIF

  overcast(:,:) = the_rth(inst_id)%radiance%overcast(:,:)
  CATCH
END SUBROUTINE rttov_get_overcast


SUBROUTINE rttov_get_rad2_upclear(err, inst_id, rad2_upclear, nchanprof)
!f2py intent(hide):: nchanprof=len(rad2_upclear)
#include "throw.h"
  USE rttov_wrapper_handle
  USE parkind1, ONLY : jpim, jprb

  IMPLICIT NONE

  INTEGER(jpim), INTENT(OUT)   :: err
  INTEGER(jpim), INTENT(IN)    :: inst_id
  INTEGER(jpim), INTENT(IN)    :: nchanprof
  REAL(jprb),    INTENT(INOUT) :: rad2_upclear(nchanprof)

#include "rttov_errorreport.interface"

  TRY
  CALL rttov_wrapper_handle_check(err, inst_id)
  THROW(err.NE.0)

  IF (.NOT. ASSOCIATED(the_rth(inst_id)%radiance2%upclear)) THEN
    err = errorstatus_fatal
    THROWM(err.NE.0, 'Wrapper radiance2 not allocated: has RTTOV been called?')
  ENDIF

  rad2_upclear(:) = the_rth(inst_id)%radiance2%upclear(:)
  CATCH
END SUBROUTINE rttov_get_rad2_upclear

SUBROUTINE rttov_get_rad2_dnclear(err, inst_id, rad2_dnclear, nchanprof)
!f2py intent(hide):: nchanprof=len(rad2_dnclear)
#include "throw.h"
  USE rttov_wrapper_handle
  USE parkind1, ONLY : jpim, jprb

  IMPLICIT NONE

  INTEGER(jpim), INTENT(OUT)   :: err
  INTEGER(jpim), INTENT(IN)    :: inst_id
  INTEGER(jpim), INTENT(IN)    :: nchanprof
  REAL(jprb),    INTENT(INOUT) :: rad2_dnclear(nchanprof)

#include "rttov_errorreport.interface"

  TRY
  CALL rttov_wrapper_handle_check(err, inst_id)
  THROW(err.NE.0)

  IF (.NOT. ASSOCIATED(the_rth(inst_id)%radiance2%dnclear)) THEN
    err = errorstatus_fatal
    THROWM(err.NE.0, 'Wrapper radiance2 not allocated: has RTTOV been called?')
  ENDIF

  rad2_dnclear(:) = the_rth(inst_id)%radiance2%dnclear(:)
  CATCH
END SUBROUTINE rttov_get_rad2_dnclear

SUBROUTINE rttov_get_rad2_refldnclear(err, inst_id, rad2_refldnclear, nchanprof)
!f2py intent(hide):: nchanprof=len(rad2_refldnclear)
#include "throw.h"
  USE rttov_wrapper_handle
  USE parkind1, ONLY : jpim, jprb

  IMPLICIT NONE

  INTEGER(jpim), INTENT(OUT)   :: err
  INTEGER(jpim), INTENT(IN)    :: inst_id
  INTEGER(jpim), INTENT(IN)    :: nchanprof
  REAL(jprb),    INTENT(INOUT) :: rad2_refldnclear(nchanprof)

#include "rttov_errorreport.interface"

  TRY
  CALL rttov_wrapper_handle_check(err, inst_id)
  THROW(err.NE.0)

  IF (.NOT. ASSOCIATED(the_rth(inst_id)%radiance2%refldnclear)) THEN
    err = errorstatus_fatal
    THROWM(err.NE.0, 'Wrapper radiance2 not allocated: has RTTOV been called?')
  ENDIF

  rad2_refldnclear(:) = the_rth(inst_id)%radiance2%refldnclear(:)
  CATCH
END SUBROUTINE rttov_get_rad2_refldnclear

SUBROUTINE rttov_get_rad2_up(err, inst_id, rad2_up, nchanprof, nlayers)
!f2py intent(hide):: nchanprof=len(rad2_up)
!f2py intent(hide):: nlayers=len(rad2_up[0])
#include "throw.h"
  USE rttov_wrapper_handle
  USE parkind1, ONLY : jpim, jprb

  IMPLICIT NONE

  INTEGER(jpim), INTENT(OUT)   :: err
  INTEGER(jpim), INTENT(IN)    :: inst_id
  INTEGER(jpim), INTENT(IN)    :: nchanprof
  INTEGER(jpim), INTENT(IN)    :: nlayers
  REAL(jprb),    INTENT(INOUT) :: rad2_up(nlayers,nchanprof)

#include "rttov_errorreport.interface"

  TRY
  CALL rttov_wrapper_handle_check(err, inst_id)
  THROW(err.NE.0)

  IF (.NOT. ASSOCIATED(the_rth(inst_id)%radiance2%up)) THEN
    err = errorstatus_fatal
    THROWM(err.NE.0, 'Wrapper radiance2 not allocated: has RTTOV been called?')
  ENDIF

  rad2_up(:,:) = the_rth(inst_id)%radiance2%up(:,:)
  CATCH
END SUBROUTINE rttov_get_rad2_up

SUBROUTINE rttov_get_rad2_down(err, inst_id, rad2_down, nchanprof, nlayers)
!f2py intent(hide):: nchanprof=len(rad2_down)
!f2py intent(hide):: nlayers=len(rad2_down[0])
#include "throw.h"
  USE rttov_wrapper_handle
  USE parkind1, ONLY : jpim, jprb

  IMPLICIT NONE

  INTEGER(jpim), INTENT(OUT)   :: err
  INTEGER(jpim), INTENT(IN)    :: inst_id
  INTEGER(jpim), INTENT(IN)    :: nchanprof
  INTEGER(jpim), INTENT(IN)    :: nlayers
  REAL(jprb),    INTENT(INOUT) :: rad2_down(nlayers,nchanprof)

#include "rttov_errorreport.interface"

  TRY
  CALL rttov_wrapper_handle_check(err, inst_id)
  THROW(err.NE.0)

  IF (.NOT. ASSOCIATED(the_rth(inst_id)%radiance2%down)) THEN
    err = errorstatus_fatal
    THROWM(err.NE.0, 'Wrapper radiance2 not allocated: has RTTOV been called?')
  ENDIF

  rad2_down(:,:) = the_rth(inst_id)%radiance2%down(:,:)
  CATCH
END SUBROUTINE rttov_get_rad2_down

SUBROUTINE rttov_get_rad2_surf(err, inst_id, rad2_surf, nchanprof, nlayers)
!f2py intent(hide):: nchanprof=len(rad2_surf)
!f2py intent(hide):: nlayers=len(rad2_surf[0])
#include "throw.h"
  USE rttov_wrapper_handle
  USE parkind1, ONLY : jpim, jprb

  IMPLICIT NONE

  INTEGER(jpim), INTENT(OUT)   :: err
  INTEGER(jpim), INTENT(IN)    :: inst_id
  INTEGER(jpim), INTENT(IN)    :: nchanprof
  INTEGER(jpim), INTENT(IN)    :: nlayers
  REAL(jprb),    INTENT(INOUT) :: rad2_surf(nlayers,nchanprof)

#include "rttov_errorreport.interface"

  TRY
  CALL rttov_wrapper_handle_check(err, inst_id)
  THROW(err.NE.0)

  IF (.NOT. ASSOCIATED(the_rth(inst_id)%radiance2%surf)) THEN
    err = errorstatus_fatal
    THROWM(err.NE.0, 'Wrapper radiance2 not allocated: has RTTOV been called?')
  ENDIF

  rad2_surf(:,:) = the_rth(inst_id)%radiance2%surf(:,:)
  CATCH
END SUBROUTINE rttov_get_rad2_surf


SUBROUTINE rttov_get_tau_total(err, inst_id, tau_total, nchanprof)
!f2py intent(hide):: nchanprof=len(tau_total)
#include "throw.h"
  USE rttov_wrapper_handle
  USE parkind1, ONLY : jpim, jprb

  IMPLICIT NONE

  INTEGER(jpim), INTENT(OUT)   :: err
  INTEGER(jpim), INTENT(IN)    :: inst_id
  INTEGER(jpim), INTENT(IN)    :: nchanprof
  REAL(jprb),    INTENT(INOUT) :: tau_total(nchanprof)

#include "rttov_errorreport.interface"

  TRY
  CALL rttov_wrapper_handle_check(err, inst_id)
  THROW(err.NE.0)

  IF (.NOT. ASSOCIATED(the_rth(inst_id)%transmission%tau_total)) THEN
    err = errorstatus_fatal
    THROWM(err.NE.0, 'Wrapper transmission not allocated: has RTTOV been called?')
  ENDIF

  tau_total(:) = the_rth(inst_id)%transmission%tau_total(:)
  CATCH
END SUBROUTINE rttov_get_tau_total

SUBROUTINE rttov_get_tau_levels(err, inst_id, tau_levels, nchanprof, nlevels)
!f2py intent(hide):: nchanprof=len(tau_levels)
!f2py intent(hide):: nlevels=len(tau_levels[0])
#include "throw.h"
  USE rttov_wrapper_handle
  USE parkind1, ONLY : jpim, jprb

  IMPLICIT NONE

  INTEGER(jpim), INTENT(OUT)   :: err
  INTEGER(jpim), INTENT(IN)    :: inst_id
  INTEGER(jpim), INTENT(IN)    :: nchanprof
  INTEGER(jpim), INTENT(IN)    :: nlevels
  REAL(jprb),    INTENT(INOUT) :: tau_levels(nlevels,nchanprof)

#include "rttov_errorreport.interface"

  TRY
  CALL rttov_wrapper_handle_check(err, inst_id)
  THROW(err.NE.0)

  IF (.NOT. ASSOCIATED(the_rth(inst_id)%transmission%tau_levels)) THEN
    err = errorstatus_fatal
    THROWM(err.NE.0, 'Wrapper transmission not allocated: has RTTOV been called?')
  ENDIF

  tau_levels(:,:) = the_rth(inst_id)%transmission%tau_levels(:,:)
  CATCH
END SUBROUTINE rttov_get_tau_levels

SUBROUTINE rttov_get_tausun_total_path2(err, inst_id, tausun_total_path2, nchanprof)
!f2py intent(hide):: nchanprof=len(tausun_total_path2)
#include "throw.h"
  USE rttov_wrapper_handle
  USE parkind1, ONLY : jpim, jprb

  IMPLICIT NONE

  INTEGER(jpim), INTENT(OUT)   :: err
  INTEGER(jpim), INTENT(IN)    :: inst_id
  INTEGER(jpim), INTENT(IN)    :: nchanprof
  REAL(jprb),    INTENT(INOUT) :: tausun_total_path2(nchanprof)

#include "rttov_errorreport.interface"

  TRY
  CALL rttov_wrapper_handle_check(err, inst_id)
  THROW(err.NE.0)

  IF (.NOT. ASSOCIATED(the_rth(inst_id)%transmission%tausun_total_path2)) THEN
    err = errorstatus_fatal
    THROWM(err.NE.0, 'Wrapper transmission not allocated: has RTTOV been called?')
  ENDIF

  tausun_total_path2(:) = the_rth(inst_id)%transmission%tausun_total_path2(:)
  CATCH
END SUBROUTINE rttov_get_tausun_total_path2

SUBROUTINE rttov_get_tausun_levels_path2(err, inst_id, tausun_levels_path2, nchanprof, nlevels)
!f2py intent(hide):: nchanprof=len(tausun_levels_path2)
!f2py intent(hide):: nlevels=len(tausun_levels_path2[0])
#include "throw.h"
  USE rttov_wrapper_handle
  USE parkind1, ONLY : jpim, jprb

  IMPLICIT NONE

  INTEGER(jpim), INTENT(OUT)   :: err
  INTEGER(jpim), INTENT(IN)    :: inst_id
  INTEGER(jpim), INTENT(IN)    :: nchanprof
  INTEGER(jpim), INTENT(IN)    :: nlevels
  REAL(jprb),    INTENT(INOUT) :: tausun_levels_path2(nlevels,nchanprof)

#include "rttov_errorreport.interface"

  TRY
  CALL rttov_wrapper_handle_check(err, inst_id)
  THROW(err.NE.0)

  IF (.NOT. ASSOCIATED(the_rth(inst_id)%transmission%tausun_levels_path2)) THEN
    err = errorstatus_fatal
    THROWM(err.NE.0, 'Wrapper transmission not allocated: has RTTOV been called?')
  ENDIF

  tausun_levels_path2(:,:) = the_rth(inst_id)%transmission%tausun_levels_path2(:,:)
  CATCH
END SUBROUTINE rttov_get_tausun_levels_path2

SUBROUTINE rttov_get_tausun_total_path1(err, inst_id, tausun_total_path1, nchanprof)
!f2py intent(hide):: nchanprof=len(tausun_total_path1)
#include "throw.h"
  USE rttov_wrapper_handle
  USE parkind1, ONLY : jpim, jprb

  IMPLICIT NONE

  INTEGER(jpim), INTENT(OUT)   :: err
  INTEGER(jpim), INTENT(IN)    :: inst_id
  INTEGER(jpim), INTENT(IN)    :: nchanprof
  REAL(jprb),    INTENT(INOUT) :: tausun_total_path1(nchanprof)

#include "rttov_errorreport.interface"

  TRY
  CALL rttov_wrapper_handle_check(err, inst_id)
  THROW(err.NE.0)

  IF (.NOT. ASSOCIATED(the_rth(inst_id)%transmission%tausun_total_path1)) THEN
    err = errorstatus_fatal
    THROWM(err.NE.0, 'Wrapper transmission not allocated: has RTTOV been called?')
  ENDIF

  tausun_total_path1(:) = the_rth(inst_id)%transmission%tausun_total_path1(:)
  CATCH
END SUBROUTINE rttov_get_tausun_total_path1

SUBROUTINE rttov_get_tausun_levels_path1(err, inst_id, tausun_levels_path1, nchanprof, nlevels)
!f2py intent(hide):: nchanprof=len(tausun_levels_path1)
!f2py intent(hide):: nlevels=len(tausun_levels_path1[0])
#include "throw.h"
  USE rttov_wrapper_handle
  USE parkind1, ONLY : jpim, jprb

  IMPLICIT NONE

  INTEGER(jpim), INTENT(OUT)   :: err
  INTEGER(jpim), INTENT(IN)    :: inst_id
  INTEGER(jpim), INTENT(IN)    :: nchanprof
  INTEGER(jpim), INTENT(IN)    :: nlevels
  REAL(jprb),    INTENT(INOUT) :: tausun_levels_path1(nlevels,nchanprof)

#include "rttov_errorreport.interface"

  TRY
  CALL rttov_wrapper_handle_check(err, inst_id)
  THROW(err.NE.0)

  IF (.NOT. ASSOCIATED(the_rth(inst_id)%transmission%tausun_levels_path1)) THEN
    err = errorstatus_fatal
    THROWM(err.NE.0, 'Wrapper transmission not allocated: has RTTOV been called?')
  ENDIF

  tausun_levels_path1(:,:) = the_rth(inst_id)%transmission%tausun_levels_path1(:,:)
  CATCH
END SUBROUTINE rttov_get_tausun_levels_path1


! The following subroutines return information from an RTTOV coefficients
! structure once an instrument has been initialised. These are not intended
! for general use.

SUBROUTINE RTTOV_GET_COEF_VAL_I0(err, inst_id, varch, I0)
!
#include "throw.h"

USE parkind1, ONLY : jpim
!INTF_OFF
USE rttov_wrapper_handle
USE rttov_const, ONLY : ngases_max
!INTF_ON
  IMPLICIT NONE
  INTEGER(jpim),      INTENT(OUT) :: err
  INTEGER(jpim),      INTENT(IN)  :: inst_id
  CHARACTER(LEN=*),   INTENT(IN)  :: varch
  INTEGER(jpim),      INTENT(OUT) :: I0
!INTF_END

#include "rttov_errorreport.interface"
!
TRY

  CALL rttov_wrapper_handle_check(err, inst_id)
  THROW(err.NE.0)

  SELECT CASE (TRIM(VARCH))
    CASE ("ID_PLATFORM")
      I0 = the_rth(inst_id)%COEFS%COEF%ID_PLATFORM
    CASE ("ID_SAT")
      I0 = the_rth(inst_id)%COEFS%COEF%ID_SAT
    CASE ("ID_INST")
      I0 = the_rth(inst_id)%COEFS%COEF%ID_INST
    CASE ("ID_SENSOR")
      I0 = the_rth(inst_id)%COEFS%COEF%ID_SENSOR
    CASE ("ID_COMP_LVL")
      I0 = the_rth(inst_id)%COEFS%COEF%ID_COMP_LVL
    CASE ("ID_COMP_PC")
      I0 = the_rth(inst_id)%COEFS%COEF%ID_COMP_PC
    CASE ("FMV_MODEL_VER")
      I0 = the_rth(inst_id)%COEFS%COEF%FMV_MODEL_VER
    CASE ("FMV_CHN")
      I0 = the_rth(inst_id)%COEFS%COEF%FMV_CHN
    CASE ("FMV_GAS")
      I0 = the_rth(inst_id)%COEFS%COEF%FMV_GAS
    CASE ("NLEVELS")
      I0 = the_rth(inst_id)%COEFS%COEF%NLEVELS
    CASE ("NLAYERS")
      I0 = the_rth(inst_id)%COEFS%COEF%NLAYERS
    CASE ("NMIXED")
      I0 = the_rth(inst_id)%COEFS%COEF%NMIXED
    CASE ("NWATER")
      I0 = the_rth(inst_id)%COEFS%COEF%NWATER
    CASE ("NOZONE")
      I0 = the_rth(inst_id)%COEFS%COEF%NOZONE
    CASE ("NWVCONT")
      I0 = the_rth(inst_id)%COEFS%COEF%NWVCONT
    CASE ("NCO2")
      I0 = the_rth(inst_id)%COEFS%COEF%NCO2
    CASE ("NN2O")
      I0 = the_rth(inst_id)%COEFS%COEF%NN2O
    CASE ("NCO")
      I0 = the_rth(inst_id)%COEFS%COEF%NCO
    CASE ("NCH4")
      I0 = the_rth(inst_id)%COEFS%COEF%NCH4
    CASE ("INCZEEMAN")
      IF( the_rth(inst_id)%COEFS%COEF%INCZEEMAN ) THEN
        I0 = 1
      ELSE
        I0 = 0
      ENDIF
    CASE ("FMV_PC_BANDS")
      I0 = the_rth(inst_id)%COEFS%COEF_PCCOMP%FMV_PC_BANDS
    CASE ("FMV_PC_CLD")
      I0 = the_rth(inst_id)%COEFS%COEF_PCCOMP%FMV_PC_CLD
    CASE ("FMV_PC_MNUM")
      I0 = the_rth(inst_id)%COEFS%COEF_PCCOMP%FMV_PC_MNUM
    CASE ("SOLARCOEF")
      IF( the_rth(inst_id)%COEFS%COEF%SOLARCOEF ) THEN
        I0 = 1
      ELSE
        I0 = 0
      ENDIF
    CASE ("NLTECOEF")
      IF( the_rth(inst_id)%COEFS%COEF%NLTECOEF) THEN
        I0 = 1
      ELSE
        I0 = 0
      ENDIF
    CASE ("PMC_SHIFT")
      IF( the_rth(inst_id)%COEFS%COEF%PMC_SHIFT ) THEN
        I0 = 1
      ELSE
        I0 = 0
      ENDIF

! water optical constant
! wave spectrum
! fastem
! ssirem                 all omited

! aliases
    CASE ("NCHANNELS")
      I0 = the_rth(inst_id)%COEFS%COEF%FMV_CHN

! array sizes
    CASE ("SIZE_FMV_GAS_ID" ,&
        & "SIZE_FMV_VAR"    ,&
        & "SIZE_FMV_COE"    ,&
        & "SIZE_FMV_INT"    ,&
        & "SIZE_FMV_LVL"    ,&
        & "SIZE_GAZ_UNITS")
      I0 = the_rth(inst_id)%COEFS%COEF%FMV_GAS

    CASE ("SIZE_FMV_GAS_POS")
      I0 = ngases_max

    CASE ("SIZE_FF_ORI_CHN" ,&
        & "SIZE_FF_VAL_CHN" ,&
        & "SIZE_FF_CWN" ,&
        & "SIZE_FF_BCO" ,&
        & "SIZE_FF_BCS" ,&
        & "SIZE_FF_GAM" ,&
        & "SIZE_TT_CHN"     ,&
        & "SIZE_TT_VAL_CHN" ,&
        & "SIZE_TT_CWN" ,&
        & "SIZE_TT_A0"  ,&
        & "SIZE_TT_A1"  ,&
        & "SIZE_PW_CHN"     ,&
        & "SIZE_PW_VAL_CHN" ,&
        & "SIZE_SS_CHN"     ,&
        & "SIZE_SS_VAL_CHN" ,&
        & "SIZE_SS_CWN"     ,&
        & "SIZE_SS_SOLAR_SPECTRUM" )
      I0 = the_rth(inst_id)%COEFS%COEF%FMV_CHN

    CASE ("SIZE_REF_PRFL_P" ,&
        & "SIZE_LIM_PRFL_P" ,&
        & "SIZE_LIM_PRFL_TMAX" ,&
        & "SIZE_LIM_PRFL_TMIN" )
      I0 = the_rth(inst_id)%COEFS%COEF%NLEVELS

    CASE ("SIZE_NOISE_IN")
      I0=the_rth(inst_id)%COEFS%COEF_PCCOMP%FMV_PC_NCHN
    CASE ("SIZE_FMV_PC_SETS")
      I0 = the_rth(inst_id)%COEFS%COEF_PCCOMP%FMV_PC_BANDS

!    CASE ("SIZE_LINE_BY_LINE")
!      I0 = 20

! 2 dimensions arrays, should be specified in I0 and I1
    CASE ("SIZE_FMV_PC_NPRED")
      I0 = 2

! array size for aliases
    CASE ("SIZE_WAVENUMBERS")
      I0 = the_rth(inst_id)%COEFS%COEF%FMV_CHN

    CASE ("SIZE_REF_PRESSURE" ,&
        & "SIZE_REF_TEMPERATURE" ,&
        & "SIZE_REF_WATERVAPOR" ,&
        & "SIZE_REF_OZONE" ,&
        & "SIZE_REF_CO2" ,&
        & "SIZE_REF_N2O" ,&
        & "SIZE_REF_CO" ,&
        & "SIZE_REF_CH4" ,&
        & "SIZE_MIN_TEMPERATURE" ,&
        & "SIZE_MIN_WATERVAPOR" ,&
        & "SIZE_MIN_OZONE" ,&
        & "SIZE_MIN_CO2" ,&
        & "SIZE_MIN_N2O" ,&
        & "SIZE_MIN_CO" ,&
        & "SIZE_MIN_CH4" ,&
        & "SIZE_MAX_TEMPERATURE" ,&
        & "SIZE_MAX_WATERVAPOR" ,&
        & "SIZE_MAX_OZONE" ,&
        & "SIZE_MAX_CO2" ,&
        & "SIZE_MAX_N2O" ,&
        & "SIZE_MAX_CO" ,&
        & "SIZE_MAX_CH4")
      I0 = the_rth(inst_id)%COEFS%COEF%NLEVELS

    CASE DEFAULT
      err = errorstatus_fatal
      THROWM(ERR.NE.0,"INTERGER SCALAR VARIABLE "//TRIM(VARCH)//" NOT FOUND")

  END SELECT

CATCH

END SUBROUTINE

SUBROUTINE RTTOV_GET_COEF_VAL_R0(err, inst_id, varch, R0)
!
#include "throw.h"

USE parkind1, ONLY : jpim, jprb
!INTF_OFF
USE rttov_wrapper_handle
USE rttov_const, ONLY : speedl
!INTF_ON
  IMPLICIT NONE
  INTEGER(jpim),      INTENT(OUT) :: err
  INTEGER(jpim),      INTENT(IN)  :: inst_id
  CHARACTER(LEN=*),   INTENT(IN)  :: varch
  REAL(jprb),         INTENT(OUT) :: R0
!INTF_END

#include "rttov_errorreport.interface"
!
TRY

  CALL rttov_wrapper_handle_check(err, inst_id)
  THROW(err.NE.0)

  SELECT CASE (TRIM(VARCH))

    CASE ("FC_SPEEDL")
      R0 = speedl
    CASE ("FC_PLANCK_C1")
      R0 = the_rth(inst_id)%COEFS%COEF%FC_PLANCK_C1
    CASE ("FC_PLANCK_C2")
      R0 = the_rth(inst_id)%COEFS%COEF%FC_PLANCK_C2
    CASE ("FC_SAT_HEIGHT")
      R0 = the_rth(inst_id)%COEFS%COEF%FC_SAT_HEIGHT

    CASE DEFAULT
      err = errorstatus_fatal
      THROWM(ERR.NE.0,"REAL SCALAR VARIABLE "//TRIM(VARCH)//" NOT FOUND")

  END SELECT

CATCH

END SUBROUTINE

SUBROUTINE RTTOV_GET_COEF_VAL_C0(err, inst_id, varch, C0)
!
#include "throw.h"

USE parkind1, ONLY : jpim
!INTF_OFF
USE rttov_wrapper_handle
!INTF_ON
  IMPLICIT NONE
  INTEGER(jpim),      INTENT(OUT) :: err
  INTEGER(jpim),      INTENT(IN)  :: inst_id
  CHARACTER(LEN=*),   INTENT(IN)  :: varch
  CHARACTER(LEN=80),  INTENT(OUT) :: C0
!INTF_END

#include "rttov_errorreport.interface"
!
TRY

  CALL rttov_wrapper_handle_check(err, inst_id)
  THROW(err.NE.0)

  SELECT CASE (TRIM(VARCH))

    CASE ("ID_CREATION")
      C0 = the_rth(inst_id)%COEFS%COEF%ID_CREATION
    CASE ("ID_COMMON_NAME")
      C0 = the_rth(inst_id)%COEFS%COEF%ID_COMMON_NAME
    CASE ("FMV_MODEL_DEF")
      C0 = the_rth(inst_id)%COEFS%COEF%FMV_MODEL_DEF

    CASE DEFAULT
      err = errorstatus_fatal
      THROWM(ERR.NE.0,"CHARACTER STRING VARIABLE "//TRIM(VARCH)//" NOT FOUND")

  END SELECT

CATCH

END SUBROUTINE

SUBROUTINE RTTOV_GET_COEF_VAL_I1(err, inst_id, varch, M, I1)
!
#include "throw.h"

USE parkind1, ONLY : jpim
!INTF_OFF
USE rttov_wrapper_handle
!INTF_ON
  IMPLICIT NONE
  INTEGER(jpim),      INTENT(OUT) :: err
  INTEGER(jpim),      INTENT(IN)  :: inst_id
  CHARACTER(LEN=*),   INTENT(IN)  :: varch
  INTEGER(jpim),      INTENT(IN)  :: M
  INTEGER(jpim),      INTENT(OUT) :: I1(M)
!INTF_END

#include "rttov_errorreport.interface"
!
TRY

  CALL rttov_wrapper_handle_check(err, inst_id)
  THROW(err.NE.0)

  SELECT CASE (TRIM(VARCH))
    CASE ("FMV_GAS_ID")
      I1 = the_rth(inst_id)%COEFS%COEF%FMV_GAS_ID
    CASE ("FMV_GAS_POS")
      I1 = the_rth(inst_id)%COEFS%COEF%FMV_GAS_POS
    CASE ("FMV_VAR")
      I1 = the_rth(inst_id)%COEFS%COEF%FMV_VAR
    CASE ("FMV_COE")
      I1 = the_rth(inst_id)%COEFS%COEF%FMV_COE
    CASE ("FMV_INT")
      I1 = the_rth(inst_id)%COEFS%COEF%FMV_INT
    CASE ("FMV_LVL")
      I1 = the_rth(inst_id)%COEFS%COEF%FMV_LVL
    CASE ("GAZ_UNITS")
      I1 = the_rth(inst_id)%COEFS%COEF%GAZ_UNITS
    CASE ("FF_ORI_CHN")
      I1 = the_rth(inst_id)%COEFS%COEF%FF_ORI_CHN
    CASE ("FF_VAL_CHN")
      I1 = the_rth(inst_id)%COEFS%COEF%FF_VAL_CHN
    CASE ("TT_CHN")
      I1 = the_rth(inst_id)%COEFS%COEF%TT_CHN
    CASE ("TT_VAL_CHN")
      I1 = the_rth(inst_id)%COEFS%COEF%TT_VAL_CHN
    CASE ("PW_CHN")
      I1 = the_rth(inst_id)%COEFS%COEF%PW_CHN
    CASE ("PW_VAL_CHN")
      I1 = the_rth(inst_id)%COEFS%COEF%PW_VAL_CHN
    CASE ("SS_CHN")
      I1 = the_rth(inst_id)%COEFS%COEF%SS_CHN
    CASE ("SS_VAL_CHN")
      I1 = the_rth(inst_id)%COEFS%COEF%SS_VAL_CHN
    CASE ("FMV_PC_SETS")
      I1 = the_rth(inst_id)%COEFS%COEF_PCCOMP%FMV_PC_SETS

    CASE ("SIZE_FMV_PC_NPRED")
      I1 = (/the_rth(inst_id)%COEFS%COEF_PCCOMP%FMV_PC_BANDS ,the_rth(inst_id)%COEFS%COEF_PCCOMP%FMV_PC_SETS(1)/)

    CASE DEFAULT
      err = errorstatus_fatal
      THROWM(ERR.NE.0,"INTEGER ARRAY VARIABLE "//TRIM(VARCH)//" NOT FOUND")

  END SELECT

CATCH

END SUBROUTINE

SUBROUTINE RTTOV_GET_COEF_VAL_I2(err, inst_id, varch, M, N, I2)
!
#include "throw.h"

USE parkind1, ONLY : jpim
!INTF_OFF
USE rttov_wrapper_handle
!INTF_ON
  IMPLICIT NONE
  INTEGER(jpim),      INTENT(OUT) :: err
  INTEGER(jpim),      INTENT(IN)  :: inst_id
  CHARACTER(LEN=*),   INTENT(IN)  :: varch
  INTEGER(jpim),      INTENT(IN)  :: M
  INTEGER(jpim),      INTENT(IN)  :: N
  INTEGER(jpim),      INTENT(OUT) :: I2(M, N)
!INTF_END

#include "rttov_errorreport.interface"
!
TRY

  CALL rttov_wrapper_handle_check(err, inst_id)
  THROW(err.NE.0)

  SELECT CASE (TRIM(VARCH))

    CASE ("FMV_PC_NPRED")
      I2(:,:) = the_rth(inst_id)%COEFS%COEF_PCCOMP%PCREG(:,:)%FMV_PC_NPRED

    CASE DEFAULT
      err = errorstatus_fatal
      THROWM(ERR.NE.0,"INTEGER ARRAY VARIABLE "//TRIM(VARCH)//" NOT FOUND")

  END SELECT

CATCH

END SUBROUTINE
!
SUBROUTINE RTTOV_GET_COEF_VAL_R1(err, inst_id, varch, M, R1)

#include "throw.h"

USE parkind1, ONLY : jpim, jprb
!INTF_OFF
USE rttov_wrapper_handle
USE rttov_const, ONLY : GID_MIXED => GAS_ID_MIXED, &
                     &  GID_OZONE => GAS_ID_OZONE, &
                     &  GID_WATERVAPOUR => GAS_ID_WATERVAPOUR, &
                     &  GID_CO2 => GAS_ID_CO2,     &
                     &  GID_N2O => GAS_ID_N2O,     &
                     &  GID_CO  => GAS_ID_CO,      &
                     &  GID_CH4 => GAS_ID_CH4
!INTF_ON
  IMPLICIT NONE
  INTEGER(jpim),      INTENT(OUT) :: err
  INTEGER(jpim),      INTENT(IN)  :: inst_id
  CHARACTER(LEN=*),   INTENT(IN)  :: varch
  INTEGER(jpim),      INTENT(IN)  :: M
  REAL(jprb),         INTENT(OUT) :: R1(M)
!INTF_END

#include "rttov_errorreport.interface"
  INTEGER(jpim) :: N
!
TRY

  CALL rttov_wrapper_handle_check(err, inst_id)
  THROW(err.NE.0)

  SELECT CASE (TRIM(VARCH))

    CASE ("FF_CWN")
      R1 = the_rth(inst_id)%COEFS%COEF%FF_CWN
    CASE ("FF_BCO")
      R1 = the_rth(inst_id)%COEFS%COEF%FF_BCO
    CASE ("FF_BCS")
      R1 = the_rth(inst_id)%COEFS%COEF%FF_BCS
    CASE ("FF_GAM")
      R1 = the_rth(inst_id)%COEFS%COEF%FF_GAM
    CASE ("TT_CWN")
      R1 = the_rth(inst_id)%COEFS%COEF%TT_CWN
    CASE ("TT_A0")
      R1 = the_rth(inst_id)%COEFS%COEF%TT_A0
    CASE ("TT_A1")
      R1 = the_rth(inst_id)%COEFS%COEF%TT_A1
    CASE ("SS_CWN")
      R1 = the_rth(inst_id)%COEFS%COEF%SS_CWN
    CASE ("SS_SOLAR_SPECTRUM")
      R1 = the_rth(inst_id)%COEFS%COEF%SS_SOLAR_SPECTRUM
    CASE ("REF_PRFL_P")
      R1 = the_rth(inst_id)%COEFS%COEF%REF_PRFL_P
    CASE ("LIM_PRFL_P")
      R1 = the_rth(inst_id)%COEFS%COEF%LIM_PRFL_P
    CASE ("LIM_PRFL_TMAX")
      R1 = the_rth(inst_id)%COEFS%COEF%LIM_PRFL_TMAX
    CASE ("LIM_PRFL_TMIN")
      R1 = the_rth(inst_id)%COEFS%COEF%LIM_PRFL_TMIN
    CASE ("NOISE_IN")
      R1 = the_rth(inst_id)%COEFS%COEF_PCCOMP%NOISE_IN

! Aliases
    CASE ("WAVENUMBERS")
      R1 = the_rth(inst_id)%COEFS%COEF%FF_CWN

    CASE ("REF_PRESSURE")
      R1 = the_rth(inst_id)%COEFS%COEF%REF_PRFL_P

    CASE ("REF_TEMPERATURE")
      N  = the_rth(inst_id)%COEFS%COEF%FMV_GAS_POS(GID_MIXED)
      R1 = the_rth(inst_id)%COEFS%COEF%REF_PRFL_T(:,N)

    CASE ("REF_WATERVAPOR")
      IF( SIZE(R1) .NE. the_rth(inst_id)%COEFS%COEF%NLEVELS) THEN
        err = errorstatus_fatal
      ENDIF
      THROWM(ERR.NE.0,"REAL ARRAY ARGUMENT DOES NOT HAVE THE CORRECT DIMENSION")
      N  = the_rth(inst_id)%COEFS%COEF%FMV_GAS_POS(GID_WATERVAPOUR)
      R1 = the_rth(inst_id)%COEFS%COEF%REF_PRFL_MR(:,N)

     CASE ("REF_OZONE")
      IF( the_rth(inst_id)%COEFS%COEF%NOZONE .LE. 0 ) THEN
        err = errorstatus_fatal
      ENDIF
      THROWM(ERR.NE.0,"DOES NOT CONTAIN OZONE")
      N  = the_rth(inst_id)%COEFS%COEF%FMV_GAS_POS(GID_OZONE)
      R1 = the_rth(inst_id)%COEFS%COEF%REF_PRFL_MR(:,N)

     CASE ("REF_CO2")
      IF( the_rth(inst_id)%COEFS%COEF%NCO2 .LE. 0 ) THEN
        err = errorstatus_fatal
      ENDIF
      THROWM(ERR.NE.0,"DOES NOT CONTAIN CO2")
      N  = the_rth(inst_id)%COEFS%COEF%FMV_GAS_POS(GID_CO2)
      R1 = the_rth(inst_id)%COEFS%COEF%REF_PRFL_MR(:,N)

     CASE ("REF_N2O")
      IF( the_rth(inst_id)%COEFS%COEF%NN2O .LE. 0 ) THEN
        err = errorstatus_fatal
      ENDIF
      THROWM(ERR.NE.0,"DOES NOT CONTAIN N2O")
      N  = the_rth(inst_id)%COEFS%COEF%FMV_GAS_POS(GID_N2O)
      R1 = the_rth(inst_id)%COEFS%COEF%REF_PRFL_MR(:,N)

     CASE ("REF_CO")
      IF( the_rth(inst_id)%COEFS%COEF%NCO .LE. 0 ) THEN
        err = errorstatus_fatal
      ENDIF
      THROWM(ERR.NE.0,"DOES NOT CONTAIN CO")
      N  = the_rth(inst_id)%COEFS%COEF%FMV_GAS_POS(GID_CO)
      R1 = the_rth(inst_id)%COEFS%COEF%REF_PRFL_MR(:,N)

     CASE ("REF_CH4")
      IF( the_rth(inst_id)%COEFS%COEF%NCH4 .LE. 0 ) THEN
        err = errorstatus_fatal
      ENDIF
      THROWM(ERR.NE.0,"DOES NOT CONTAIN CH4")
      N  = the_rth(inst_id)%COEFS%COEF%FMV_GAS_POS(GID_CH4)
      R1 = the_rth(inst_id)%COEFS%COEF%REF_PRFL_MR(:,N)


    CASE ("MIN_TEMPERATURE")
      R1 = the_rth(inst_id)%COEFS%COEF%LIM_PRFL_TMIN(:)

    CASE ("MIN_WATERVAPOR")
      IF( SIZE(R1) .NE. the_rth(inst_id)%COEFS%COEF%NLEVELS) THEN
        err = errorstatus_fatal
      ENDIF
      THROWM(ERR.NE.0,"REAL ARRAY ARGUMENT DOES NOT HAVE THE CORRECT DIMENSION")
      N  = the_rth(inst_id)%COEFS%COEF%FMV_GAS_POS(GID_WATERVAPOUR)
      R1 = the_rth(inst_id)%COEFS%COEF%LIM_PRFL_GMIN(:,N)

     CASE ("MIN_OZONE")
      IF( the_rth(inst_id)%COEFS%COEF%NOZONE .LE. 0 ) THEN
        err = errorstatus_fatal
      ENDIF
      THROWM(ERR.NE.0,"DOES NOT CONTAIN OZONE")
      N  = the_rth(inst_id)%COEFS%COEF%FMV_GAS_POS(GID_OZONE)
      R1 = the_rth(inst_id)%COEFS%COEF%LIM_PRFL_GMIN(:,N)

     CASE ("MIN_CO2")
      IF( the_rth(inst_id)%COEFS%COEF%NCO2 .LE. 0 ) THEN
        err = errorstatus_fatal
      ENDIF
      THROWM(ERR.NE.0,"DOES NOT CONTAIN CO2")
      N  = the_rth(inst_id)%COEFS%COEF%FMV_GAS_POS(GID_CO2)
      R1 = the_rth(inst_id)%COEFS%COEF%LIM_PRFL_GMIN(:,N)

     CASE ("MIN_N2O")
      IF( the_rth(inst_id)%COEFS%COEF%NN2O .LE. 0 ) THEN
        err = errorstatus_fatal
      ENDIF
      THROWM(ERR.NE.0,"DOES NOT CONTAIN N2O")
      N  = the_rth(inst_id)%COEFS%COEF%FMV_GAS_POS(GID_N2O)
      R1 = the_rth(inst_id)%COEFS%COEF%LIM_PRFL_GMIN(:,N)

     CASE ("MIN_CO")
      IF( the_rth(inst_id)%COEFS%COEF%NCO .LE. 0 ) THEN
        err = errorstatus_fatal
      ENDIF
      THROWM(ERR.NE.0,"DOES NOT CONTAIN CO")
      N  = the_rth(inst_id)%COEFS%COEF%FMV_GAS_POS(GID_CO)
      R1 = the_rth(inst_id)%COEFS%COEF%LIM_PRFL_GMIN(:,N)

     CASE ("MIN_CH4")
      IF( the_rth(inst_id)%COEFS%COEF%NCH4 .LE. 0 ) THEN
        err = errorstatus_fatal
      ENDIF
      THROWM(ERR.NE.0,"DOES NOT CONTAIN CH4")
      N  = the_rth(inst_id)%COEFS%COEF%FMV_GAS_POS(GID_CH4)
      R1 = the_rth(inst_id)%COEFS%COEF%LIM_PRFL_GMIN(:,N)


    CASE ("MAX_TEMPERATURE")
      R1 = the_rth(inst_id)%COEFS%COEF%LIM_PRFL_TMAX(:)

    CASE ("MAX_WATERVAPOR")
      IF( SIZE(R1) .NE. the_rth(inst_id)%COEFS%COEF%NLEVELS) THEN
        err = errorstatus_fatal
      ENDIF
      THROWM(ERR.NE.0,"REAL ARRAY ARGUMENT DOES NOT HAVE THE CORRECT DIMENSION")
      N  = the_rth(inst_id)%COEFS%COEF%FMV_GAS_POS(GID_WATERVAPOUR)
      R1 = the_rth(inst_id)%COEFS%COEF%LIM_PRFL_GMAX(:,N)

     CASE ("MAX_OZONE")
      IF( the_rth(inst_id)%COEFS%COEF%NOZONE .LE. 0 ) THEN
        err = errorstatus_fatal
      ENDIF
      THROWM(ERR.NE.0,"DOES NOT CONTAIN OZONE")
      N  = the_rth(inst_id)%COEFS%COEF%FMV_GAS_POS(GID_OZONE)
      R1 = the_rth(inst_id)%COEFS%COEF%LIM_PRFL_GMAX(:,N)

     CASE ("MAX_CO2")
      IF( the_rth(inst_id)%COEFS%COEF%NCO2 .LE. 0 ) THEN
        err = errorstatus_fatal
      ENDIF
      THROWM(ERR.NE.0,"DOES NOT CONTAIN CO2")
      N  = the_rth(inst_id)%COEFS%COEF%FMV_GAS_POS(GID_CO2)
      R1 = the_rth(inst_id)%COEFS%COEF%LIM_PRFL_GMAX(:,N)

     CASE ("MAX_N2O")
      IF( the_rth(inst_id)%COEFS%COEF%NN2O .LE. 0 ) THEN
        err = errorstatus_fatal
      ENDIF
      THROWM(ERR.NE.0,"DOES NOT CONTAIN N2O")
      N  = the_rth(inst_id)%COEFS%COEF%FMV_GAS_POS(GID_N2O)
      R1 = the_rth(inst_id)%COEFS%COEF%LIM_PRFL_GMAX(:,N)

     CASE ("MAX_CO")
      IF( the_rth(inst_id)%COEFS%COEF%NCO .LE. 0 ) THEN
        err = errorstatus_fatal
      ENDIF
      THROWM(ERR.NE.0,"DOES NOT CONTAIN CO")
      N  = the_rth(inst_id)%COEFS%COEF%FMV_GAS_POS(GID_CO)
      R1 = the_rth(inst_id)%COEFS%COEF%LIM_PRFL_GMAX(:,N)

     CASE ("MAX_CH4")
      IF( the_rth(inst_id)%COEFS%COEF%NCH4 .LE. 0 ) THEN
        err = errorstatus_fatal
      ENDIF
      THROWM(ERR.NE.0,"DOES NOT CONTAIN CH4")
      N  = the_rth(inst_id)%COEFS%COEF%FMV_GAS_POS(GID_CH4)
      R1 = the_rth(inst_id)%COEFS%COEF%LIM_PRFL_GMAX(:,N)

    CASE DEFAULT
      err = errorstatus_fatal
      THROWM(ERR.NE.0,"REAL ARRAY VARIABLE "//TRIM(VARCH)//" NOT FOUND")

  END SELECT

CATCH

END SUBROUTINE
