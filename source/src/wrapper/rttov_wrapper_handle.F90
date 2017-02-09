MODULE rttov_wrapper_handle
! Description
!   Module for the RTTOV wrapper: internal subroutines
!     for allocation/deallocation and reading coefs.
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

#include "throw.h"

  USE parkind1
  USE rttov_types
  USE rttov_const, ONLY : naer_max, nwcl_max

  IMPLICIT NONE

  ! Maximum number of concurrent instruments
  INTEGER(jpim), PARAMETER :: max_ninst = 100


  INTEGER(jpim) :: iii

  ! Gas ID definitions
  INTEGER(jpim), PARAMETER :: gas_id_q   = 1
  INTEGER(jpim), PARAMETER :: gas_id_o3  = 2
  INTEGER(jpim), PARAMETER :: gas_id_co2 = 3
  INTEGER(jpim), PARAMETER :: gas_id_n2o = 4
  INTEGER(jpim), PARAMETER :: gas_id_co  = 5
  INTEGER(jpim), PARAMETER :: gas_id_ch4 = 6

  INTEGER(jpim), PARAMETER :: gas_id_clw = 15

  INTEGER(jpim), PARAMETER :: gas_id_cfrac = 20
  INTEGER(jpim), PARAMETER :: gas_id_lwc(1:nwcl_max) = (/ (iii, iii = 21, 25) /)
  INTEGER(jpim), PARAMETER :: gas_id_iwc   = 30
  INTEGER(jpim), PARAMETER :: gas_id_icede = 31

  INTEGER(jpim), PARAMETER :: gas_id_aer(1:naer_max) = (/ (iii, iii = 41, 53) /)

  TYPE rttovwrapperhandle_type
    LOGICAL(jplm)           :: init = .FALSE.

    CHARACTER(LEN=256)      :: file_coef
    CHARACTER(LEN=256)      :: file_scaer
    CHARACTER(LEN=256)      :: file_sccld
    CHARACTER(LEN=256)      :: file_pccoef
!     CHARACTER(LEN=256)      :: file_mietable

    LOGICAL(jplm)           :: verbose         = .FALSE.
    LOGICAL(jplm)           :: check_opts      = .TRUE.
    LOGICAL(jplm)           :: store_rad       = .FALSE.
    LOGICAL(jplm)           :: store_rad2      = .FALSE.
    LOGICAL(jplm)           :: store_trans     = .FALSE.
    INTEGER(jpim)           :: nthreads        = 1
    INTEGER(jpim)           :: nprofs_per_call = 1

    TYPE(rttov_options)     :: opts
    TYPE(rttov_coefs)       :: coefs

    TYPE(transmission_type) :: transmission
    TYPE(radiance_type)     :: radiance
    TYPE(radiance2_type)    :: radiance2
  ENDTYPE rttovwrapperhandle_type

  TYPE(rttovwrapperhandle_type), SAVE :: the_rth(max_ninst)

CONTAINS

SUBROUTINE rttov_wrapper_handle_defaults(inst_id)
  !
  ! Set handle variables to default values
  !

  INTEGER(jpim), INTENT(IN)  :: inst_id

  TYPE(rttov_options) :: opts

  the_rth(inst_id)%init            = .FALSE.
  the_rth(inst_id)%file_coef       = ''
  the_rth(inst_id)%file_scaer      = ''
  the_rth(inst_id)%file_sccld      = ''
  the_rth(inst_id)%file_pccoef     = ''
  the_rth(inst_id)%verbose         = .FALSE.
  the_rth(inst_id)%check_opts      = .TRUE.
  the_rth(inst_id)%store_rad       = .FALSE.
  the_rth(inst_id)%store_rad2      = .FALSE.
  the_rth(inst_id)%store_trans     = .FALSE.
  the_rth(inst_id)%nthreads        = 1
  the_rth(inst_id)%nprofs_per_call = 1
  the_rth(inst_id)%opts            = opts
END SUBROUTINE rttov_wrapper_handle_defaults

SUBROUTINE rttov_wrapper_handle_check(err, inst_id)
  !
  ! Check instrument has been initialised
  !

  INTEGER(jpim), INTENT(OUT) :: err
  INTEGER(jpim), INTENT(IN)  :: inst_id

  TRY

  IF (inst_id < 1 .OR. inst_id > max_ninst) THEN
    err = errorstatus_fatal
    THROWM(err.NE.0, 'inst_id outside valid range')
  ENDIF

  IF (.NOT. the_rth(inst_id)%init) THEN
    err = errorstatus_fatal
    THROWM(err.NE.0, 'No coefficients file loaded for this inst_id: call rttov_wrapper_load first')
  ENDIF

  CATCH
END SUBROUTINE rttov_wrapper_handle_check


SUBROUTINE rttov_wrapper_handle_load(err, rth, channels)
  !
  ! Read coefficients and check consistency with options
  !

  INTEGER(jpim),                 INTENT(OUT)   :: err
  TYPE(rttovwrapperhandle_type), INTENT(INOUT) :: rth
  INTEGER(jpim),                 INTENT(IN)    :: channels(:)

#include "rttov_read_coefs.interface"

  LOGICAL(jplm) :: allchannels, addclouds, addaerosl

  TRY

  allchannels = (SIZE(channels) == 1 .AND. channels(1) <= 0)

  ! Ensure cloud/aerosol coefs are read only if filenames were specified
  addclouds = rth%opts%rt_ir%addclouds
  addaerosl = rth%opts%rt_ir%addaerosl
  rth%opts%rt_ir%addclouds = (rth%file_sccld .NE. '')
  rth%opts%rt_ir%addaerosl = (rth%file_scaer .NE. '')

  IF (allchannels) THEN
    CALL rttov_read_coefs(err, rth%coefs, rth%opts, &
                          file_coef   = rth%file_coef,  &
                          file_scaer  = rth%file_scaer, &
                          file_sccld  = rth%file_sccld, &
                          file_pccoef = rth%file_pccoef)
  ELSE
    CALL rttov_read_coefs(err, rth%coefs, rth%opts, &
                          file_coef   = rth%file_coef,  &
                          file_scaer  = rth%file_scaer, &
                          file_sccld  = rth%file_sccld, &
                          file_pccoef = rth%file_pccoef,&
                          channels = channels)
  ENDIF

  ! Revert the options to original values
  rth%opts%rt_ir%addclouds = addclouds
  rth%opts%rt_ir%addaerosl = addaerosl

  THROWM(err.NE.0, 'Error reading coefficient files')

  CATCH
END SUBROUTINE rttov_wrapper_handle_load


SUBROUTINE rttov_wrapper_handle_alloc(err, rth, nprofiles, nchannels, nlevels, asw)
  !
  ! Allocate/deallocate wrapper handle output structures
  !

  INTEGER(jpim),                 INTENT(OUT)   :: err
  TYPE(rttovwrapperhandle_type), INTENT(INOUT) :: rth
  INTEGER(jpim),                 INTENT(IN)    :: nprofiles
  INTEGER(jpim),                 INTENT(IN)    :: nchannels
  INTEGER(jpim),                 INTENT(IN)    :: nlevels
  INTEGER(jpim),                 INTENT(IN)    :: asw

  INTEGER(jpim) :: nchanprof, nlayers

#include "rttov_alloc_rad.interface"
#include "rttov_alloc_transmission.interface"

  TRY

  nchanprof = nchannels * nprofiles
  nlayers = nlevels - 1_jpim

  IF (asw == 1_jpim) THEN

    IF (rth%store_rad) THEN
      IF (rth%store_rad2) THEN
        CALL rttov_alloc_rad(err, nchanprof, rth%radiance, nlayers, asw, rth%radiance2)
      ELSE
        CALL rttov_alloc_rad(err, nchanprof, rth%radiance, nlayers, asw)
      ENDIF
      THROWM(err.NE.0, 'Error allocating wrapper radiance')
    ENDIF

    IF (rth%store_trans) THEN
      CALL rttov_alloc_transmission(err, rth%transmission, nlayers, nchanprof, asw)
      THROWM(err.NE.0, 'Error allocating wrapper transmission')
    ENDIF

  ELSE IF (asw == 0_jpim) THEN

    IF (rth%store_rad) THEN
      IF (ASSOCIATED(rth%radiance%clear)) THEN
        IF (rth%store_rad2) THEN
          CALL rttov_alloc_rad(err, nchanprof, rth%radiance, nlayers, asw, rth%radiance2)
        ELSE
          CALL rttov_alloc_rad(err, nchanprof, rth%radiance, nlayers, asw)
        ENDIF
        THROWM(err.NE.0, 'Error deallocating wrapper radiance')
      ENDIF
    ENDIF

    IF (rth%store_trans) THEN
      IF (ASSOCIATED(rth%transmission%tau_total)) THEN
        CALL rttov_alloc_transmission(err, rth%transmission, nlayers, nchanprof, asw)
        THROWM(err.NE.0, 'Error deallocating wrapper transmission')
      ENDIF
    ENDIF

  ENDIF

  CATCH
END SUBROUTINE rttov_wrapper_handle_alloc


SUBROUTINE rttov_wrapper_nullify_structs(rth)
  !
  ! Nullify output structures
  !
  TYPE(rttovwrapperhandle_type), INTENT(INOUT) :: rth

  NULLIFY(rth%radiance%clear,      &
          rth%radiance%total,      &
          rth%radiance%bt_clear,   &
          rth%radiance%bt,         &
          rth%radiance%refl,       &
          rth%radiance%refl_clear, &
          rth%radiance%overcast)
  NULLIFY(rth%radiance2%up,      &
          rth%radiance2%down,    &
          rth%radiance2%surf,    &
          rth%radiance2%upclear, &
          rth%radiance2%dnclear, &
          rth%radiance2%refldnclear)
  NULLIFY(rth%transmission%tau_total,           &
          rth%transmission%tau_levels,          &
          rth%transmission%tausun_total_path2,  &
          rth%transmission%tausun_levels_path2, &
          rth%transmission%tausun_total_path1,  &
          rth%transmission%tausun_levels_path1)
END SUBROUTINE


SUBROUTINE rttov_wrapper_handle_drop(err, rth)
  !
  ! Deallocate wrapper handle structure
  !

  INTEGER(jpim),                 INTENT(OUT)   :: err
  TYPE(rttovwrapperhandle_type), INTENT(INOUT) :: rth

#include "rttov_dealloc_coefs.interface"

  TRY

  IF (rth%init) THEN
    rth%init = .FALSE.

    CALL rttov_dealloc_coefs(err, rth%coefs)
    THROWM(err.NE.0, 'Error deallocating coefficients')

    ! Deallocate RTTOV output structures - for deallocation the actual values
    !   of nprofiles, nchannels and nlevels arguments don't matter
    CALL rttov_wrapper_handle_alloc(err, rth, 1_jpim, 1_jpim, 1_jpim, 0_jpim)
    THROW(err.NE.0)
  ENDIF

  CATCH
END SUBROUTINE rttov_wrapper_handle_drop


SUBROUTINE rttov_wrapper_handle_drop_all(err)
  !
  ! Deallocate all allocated wrapper handles
  !

  INTEGER(jpim), INTENT(OUT)   :: err

  INTEGER(jpim)      :: inst_id
  CHARACTER(LEN=128) :: msg

#include "rttov_dealloc_coefs.interface"

  TRY

  DO inst_id = 1, max_ninst
    IF (the_rth(inst_id)%init) THEN
      CALL rttov_wrapper_handle_drop(err, the_rth(inst_id))
      IF (err .NE. 0) THEN
        WRITE(msg,'(a,i4)') 'Error deallocating wrapper for inst_id ', inst_id
        THROWM(err.NE.0, TRIM(msg))
      ENDIF
    ENDIF
  ENDDO

  CATCH
END SUBROUTINE rttov_wrapper_handle_drop_all


END MODULE rttov_wrapper_handle
