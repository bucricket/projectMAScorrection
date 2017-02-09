! Description:
!> @file
!!   Deallocate memory for emissivity atlas.
!
!> @brief
!!   Deallocate memory for emissivity atlas.
!!
!! @details
!!   If rttov_setup_emis_atlas was called for both an IR
!!   and a MW instrument you should call this once with each
!!   coefs structure to deallocate both the IR and MW arrays.
!!
!! @param[in]   coefs   coefficients structure for instrument to simulate
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
SUBROUTINE rttov_deallocate_emis_atlas(coefs)

  USE rttov_types, ONLY : &
        & rttov_coefs
!INTF_OFF
  USE rttov_const, ONLY :      &
        & sensor_id_mw,        &
        & sensor_id_po

  USE mod_iratlas, ONLY :             &
        & rttov_uwiremis_close_atlas, &
        & uw_ir_atlas_version => ir_atlas_version

  USE mod_mwatlas, ONLY :      &
        & rttov_closemw_atlas, &
        & telsem_mw_atlas_version => mw_atlas_version

  USE mod_cnrm_mw_atlas, ONLY :     &
        & cnrm_mw_atlas_version => mw_atlas_version

  USE mod_rttov_emis_atlas, ONLY : &
        & ir_atlas_version,       &
        & mw_atlas_version,       &
        & ir_atlas_init,          &
        & mw_atlas_init
!INTF_ON
  IMPLICIT NONE

  TYPE(rttov_coefs), INTENT(IN) :: coefs ! RTTOV instrument coefficients

!INTF_END

  IF (coefs%coef%id_sensor == sensor_id_mw .OR. coefs%coef%id_sensor == sensor_id_po) THEN

    ! MW atlas

    IF (mw_atlas_init) THEN
      IF (mw_atlas_version == telsem_mw_atlas_version) THEN
        CALL rttov_closemw_atlas
        mw_atlas_init = .FALSE.
      END IF

      IF (mw_atlas_version == cnrm_mw_atlas_version) THEN
        ! nothing to do for CNRM MW atlas
        mw_atlas_init = .FALSE.
      END IF
    END IF

  ELSE

    ! IR atlas

    IF (ir_atlas_init) THEN
      CALL rttov_uwiremis_close_atlas
      ir_atlas_init = .FALSE.
    END IF

  END IF

END SUBROUTINE rttov_deallocate_emis_atlas
