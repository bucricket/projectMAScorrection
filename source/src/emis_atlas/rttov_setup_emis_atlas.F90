! Description:
!> @file
!!   Sets up the emissivity atlas appropriate to the sensor
!!   type. Data is loaded for the specified month.
!
!> @brief
!!   Sets up the emissivity atlas appropriate to the sensor
!!   type. Data is loaded for the specified month.
!!
!! @details
!!   This subroutine should be called after the coefficient
!!   file has been read into the coefs structure. The atlas
!!   corresponding to the instrument type (IR or MW) will
!!   be initialised. The IR and MW atlases can be initialised
!!   using multiple calls with different coefs arguments.
!!
!!   If calling the atlas for multiple IR instruments this
!!   subroutine only needs to be called once, and similarly
!!   for multiple MW instruments.
!!
!!   It is possible to choose between different atlas versions.
!!   See the user guide for the relevant values to pass. Note
!!   that only one IR and one MW atlas may be initialised at
!!   any given time.
!!
!!   There are additional options for the IR emissivity atlas:
!!   - the standard deviation information can optionally be
!!     returned by the atlas: do not select this if it is not
!!     required as the memory requirements are relatively large
!!   - if you are only calling the IR atlas for one instrument
!!     you should set ir_atlas_single_instrument=.TRUE. The
!!     calls to rttov_get_emis are then much faster.
!!   - the IR emissivites can optionally include a correction
!!     for the zenith angle. If this option is selected you
!!     must download the additional angular correction atlas
!!     files from the RTTOV website.
!!
!! @param[out]  err                         status on exit
!! @param[in]   opts                        options to configure the simulations
!! @param[in]   imonth                      month (1-12 => Jan-Dec) for which atlas should be initialised
!! @param[in]   coefs                       coefficients structure for instrument to simulate
!! @param[in]   path                        path to directory containing emissivity atlas data, optional
!! @param[in]   ir_atlas_read_std           return standard deviations for IR emissivity atlas, optional
!! @param[in]   ir_atlas_single_instrument  initialise IR emissivity atlas for a single instrument, optional
!! @param[in]   ir_atlas_ang_corr           apply zenith angle correction to IR emissivities, optional
!! @param[in]   ir_atlas_ver                version number for IR atlas, optional
!! @param[in]   mw_atlas_ver                version number for MW atlas, optional
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
SUBROUTINE rttov_setup_emis_atlas(          &
             &  err,                        &! out
             &  opts,                       &! in
             &  imonth,                     &! in
             &  coefs,                      &! in
             &  path,                       &! in, optional
             &  ir_atlas_read_std,          &! in, optional
             &  ir_atlas_single_instrument, &! in, optional
             &  ir_atlas_ang_corr,          &! in, optional
             &  ir_atlas_ver,               &! in, optional
             &  mw_atlas_ver)                ! in, optional
!INTF_OFF
#include "throw.h"
!INTF_ON
  USE parkind1, ONLY : jpim, jplm

  USE rttov_types, ONLY : &
        & rttov_coefs,    &
        & rttov_options
!INTF_OFF
  USE rttov_const, ONLY :      &
        & sensor_id_mw,        &
        & sensor_id_po,        &
        & errorstatus_success, &
        & errorstatus_fatal

  USE mod_iratlas, ONLY :      &
        & rttov_uwiremis_init, &
        & uw_ir_atlas_version => ir_atlas_version

  USE mod_mwatlas, ONLY :     &
        & rttov_readmw_atlas, &
        & telsem_mw_atlas_version => mw_atlas_version

  USE mod_cnrm_mw_atlas, ONLY :     &
        & rttov_cnrmmwemis_init,    &
        & cnrm_mw_atlas_version => mw_atlas_version

  USE mod_rttov_emis_atlas, ONLY : &
        & ir_atlas_version,        &
        & mw_atlas_version,        &
        & ir_atlas_init,           &
        & mw_atlas_init,           &
        & ir_atlas_std_init,       &
        & ir_atlas_single_inst,    &
        & ir_atlas_do_ang_corr
!INTF_ON
  IMPLICIT NONE

  TYPE(rttov_options), INTENT(IN)          :: opts
  INTEGER(KIND=jpim), INTENT(IN)           :: imonth                     ! month for which atlas data required
  TYPE(rttov_coefs),  INTENT(IN)           :: coefs                      ! RTTOV instrument coefficients
  CHARACTER(LEN=*),   INTENT(IN), OPTIONAL :: path                       ! path to atlas data (if not the default)
  LOGICAL(KIND=jplm), INTENT(IN), OPTIONAL :: ir_atlas_read_std          ! flag to indicate IR stdv are required
  LOGICAL(KIND=jplm), INTENT(IN), OPTIONAL :: ir_atlas_single_instrument ! flag to indicate IR atlas to be initialised
                                                                         !   for a single instrument only (faster)
  LOGICAL(KIND=jplm), INTENT(IN), OPTIONAL :: ir_atlas_ang_corr          ! flag to indicate IR angular correction
                                                                         !   should be included
  INTEGER(KIND=jpim), INTENT(IN), OPTIONAL :: ir_atlas_ver               ! version of IR atlas to use
  INTEGER(KIND=jpim), INTENT(IN), OPTIONAL :: mw_atlas_ver               ! version of MW atlas to use

  INTEGER(KIND=jpim), INTENT(OUT)          :: err                        ! output error status

!INTF_END

  CHARACTER(LEN=300) :: fpath
  CHARACTER(LEN=300) :: msg

!----------------------------------------------------------------------------
  TRY

  IF (PRESENT(path)) THEN
    fpath = TRIM(path)//'/'
  ELSE
    fpath = './'
  END IF

  IF (PRESENT(ir_atlas_read_std)) THEN
    ir_atlas_std_init = ir_atlas_read_std
  ELSE
    ir_atlas_std_init = .FALSE.
  ENDIF

  IF (PRESENT(ir_atlas_single_instrument)) THEN
    ir_atlas_single_inst = ir_atlas_single_instrument
  ELSE
    ir_atlas_single_inst = .FALSE.
  ENDIF

  IF (PRESENT(ir_atlas_ang_corr)) THEN
    ir_atlas_do_ang_corr = ir_atlas_ang_corr
  ELSE
    ir_atlas_do_ang_corr = .FALSE.
  ENDIF

  ! default values for atlas versions
  mw_atlas_version = telsem_mw_atlas_version
  ir_atlas_version = uw_ir_atlas_version

  IF (PRESENT(ir_atlas_ver)) ir_atlas_version = ir_atlas_ver
  IF (PRESENT(mw_atlas_ver)) mw_atlas_version = mw_atlas_ver

  IF (coefs%coef%id_sensor == sensor_id_mw .OR. coefs%coef%id_sensor == sensor_id_po) THEN

    ! MW atlas

    IF (mw_atlas_init) THEN
      IF (opts%config%verbose) WARN('MW emissivity atlas already initialised.')
    ELSE
      IF (mw_atlas_version == telsem_mw_atlas_version) THEN
        CALL rttov_readmw_atlas(TRIM(fpath),imonth,opts%config%verbose,err)
        THROWM(err .NE. 0, 'Error initialising TELSEM MW emissivity atlas.')
      ELSEIF (mw_atlas_version == cnrm_mw_atlas_version) THEN
        CALL rttov_cnrmmwemis_init(TRIM(fpath),imonth,coefs%coef,opts%config%verbose,err)
        THROWM(err .NE. 0, 'Error initialising CNRM MW emissivity atlas.')
      ELSE
        WRITE(msg,'(a,i5)') 'Unknown MW atlas version: ', mw_atlas_version
        err = errorstatus_fatal
        THROWM(err .NE. 0, msg)
      END IF
      IF (err == 0) mw_atlas_init = .TRUE.
    END IF
  ELSE

    ! IR atlas

    IF (ir_atlas_init) THEN
      IF (opts%config%verbose) WARN('IR emissivity atlas already initialised.')
    ELSE
      IF (ir_atlas_version == uw_ir_atlas_version) THEN
        IF (ir_atlas_single_inst) THEN
          CALL rttov_uwiremis_init(TRIM(fpath),imonth,opts%config%verbose,err,coefs%coef%ff_cwn(:))
        ELSE
          CALL rttov_uwiremis_init(TRIM(fpath),imonth,opts%config%verbose,err)
        ENDIF
        THROWM(err .NE. 0, 'Error initialising IR emissivity atlas.')
      ELSE
        WRITE(msg,'(a,i5)') 'Unknown IR atlas version: ', ir_atlas_version
        err = errorstatus_fatal
        THROWM(err .NE. 0, msg)
      END IF
      IF (err == 0) ir_atlas_init = .TRUE.
    END IF
  END IF

  CATCH

END SUBROUTINE rttov_setup_emis_atlas
