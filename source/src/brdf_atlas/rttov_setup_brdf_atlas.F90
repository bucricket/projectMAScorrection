! Description:
!> @file
!!   Sets up the BRDF atlas. Data is loaded for the
!!   specified month.
!
!> @brief
!!   Sets up the BRDF atlas. Data is loaded for the
!!   specified month.
!!
!! @details
!!   This subroutine should be called after the coefficient
!!   file has been read into the coefs structure. It should
!!   only be called for v9 predictor coefficient files which
!!   include visible/near-IR channels.
!!
!!   If calling the atlas for multiple instruments this
!!   subroutine only needs to be called once.
!!
!!   It is possible to choose between different atlas versions.
!!   See the user guide for the relevant values to pass. Note
!!   that only one version of the atlas may be initialised at
!!   any given time.
!!
!!   If you are only calling the BRDF atlas for one instrument
!!   you should set brdf_atlas_single_instrument=.TRUE. The
!!   calls to rttov_get_brdf are then much faster.
!!
!! @param[out]  err                           status on exit
!! @param[in]   opts                          options to configure the simulations
!! @param[in]   imonth                        month (1-12 => Jan-Dec) for which atlas should be initialised
!! @param[in]   coefs                         coefficients structure for instrument to simulate
!! @param[in]   path                          path to directory containing BRDF atlas data, optional
!! @param[in]   brdf_atlas_single_instrument  initialise BRDF atlas for a single instrument, optional
!! @param[in]   vn_atlas_ver                  version number for BRDF atlas, optional
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
!
SUBROUTINE rttov_setup_brdf_atlas( &
             &  err,                          &! out
             &  opts,                         &! in
             &  imonth,                       &! in
             &  coefs,                        &! in
             &  path,                         &! in, optional
             &  brdf_atlas_single_instrument, &! in, optional
             &  vn_atlas_ver)                  ! in, optional
!INTF_OFF
#include "throw.h"
!INTF_ON

  USE parkind1, ONLY : jpim, jplm

  USE rttov_types, ONLY : &
        & rttov_coefs,    &
        & rttov_options
!INTF_OFF

  USE rttov_const, ONLY : &
        & sensor_id_mw,   &
        & sensor_id_po

  USE mod_brdf_atlas, ONLY :      &
        & rttov_visnirbrdf_init, &
        & cms_vn_atlas_version => vn_atlas_version

  USE mod_rttov_brdf_atlas, ONLY : &
        & vn_atlas_version,        &
        & vn_atlas_init,           &
        & brdf_atlas_single_inst
!INTF_ON

  IMPLICIT NONE

  TYPE(rttov_options), INTENT(IN)          :: opts
  INTEGER(KIND=jpim), INTENT(IN)           :: imonth       ! month for which atlas data required
  TYPE(rttov_coefs),   INTENT(IN)          :: coefs        ! RTTOV instrument coefficients
  CHARACTER(LEN=*),   INTENT(IN), OPTIONAL :: path         ! path to atlas data (if not the default)
  LOGICAL(KIND=jplm), INTENT(IN), OPTIONAL :: brdf_atlas_single_instrument ! flag to indicate BRDF atlas to be init
                                                                           ! for a single instrument only (faster)
  INTEGER(KIND=jpim), INTENT(IN), OPTIONAL :: vn_atlas_ver ! version of vis/nir atlas to use

  INTEGER(KIND=jpim), INTENT(OUT)          :: err          ! output error status

!INTF_END

  CHARACTER(LEN=300) :: fpath
  CHARACTER(LEN=300) :: msg

!----------------------------------------------------------------------------
  TRY

  IF (coefs%coef%id_sensor == sensor_id_mw .OR. &
      coefs%coef%id_sensor == sensor_id_po .OR. &
      ALL(coefs%coef%ss_val_chn(:) == 0_jpim)) THEN
    err = errorstatus_fatal
    THROWM(err .NE. 0, 'BRDF atlas only valid for sensors with VIS/NIR channels')
  ENDIF

  IF (PRESENT(path)) THEN
    fpath = TRIM(path)//'/'
  ELSE
    fpath = './'
  END IF

  IF (PRESENT(brdf_atlas_single_instrument)) THEN
    brdf_atlas_single_inst = brdf_atlas_single_instrument
  ELSE
    brdf_atlas_single_inst = .FALSE.
  ENDIF

  ! default values for Atlas versions
  vn_atlas_version = cms_vn_atlas_version

  IF (PRESENT(vn_atlas_ver)) vn_atlas_version = vn_atlas_ver

  IF (vn_atlas_init) THEN
    IF (opts%config%verbose) WARN('vis/nir brdf atlas already initialised.')
  ELSE
    IF (vn_atlas_version == cms_vn_atlas_version) THEN
      IF (brdf_atlas_single_inst) THEN
        CALL rttov_visnirbrdf_init(TRIM(fpath),imonth,opts%config%verbose,err,coefs%coef%ff_cwn(:))
      ELSE
        CALL rttov_visnirbrdf_init(TRIM(fpath),imonth,opts%config%verbose,err)
      ENDIF
      THROWM(err .NE. 0, 'Error initialising vis/nir brdf atlas.')
    ELSE
      WRITE(msg,'(a,i5)') 'Unknown vis/nir brdf atlas version: ', vn_atlas_version
      err = errorstatus_fatal
      THROWM(err .NE. 0, msg)
    END IF
    IF (err == 0) vn_atlas_init = .TRUE.
  END IF

  CATCH

END SUBROUTINE rttov_setup_brdf_atlas
