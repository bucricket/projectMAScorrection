! Description:
!> @file
!!   Deallocate memory for BRDF atlas.
!
!> @brief
!!   Deallocate memory for BRDF atlas.
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
SUBROUTINE rttov_deallocate_brdf_atlas(coefs)

  USE rttov_types, ONLY : &
        & rttov_coefs
!INTF_OFF
  USE rttov_const, ONLY : &
        & sensor_id_ir

  USE mod_brdf_atlas, ONLY :      &
        & rttov_visnirbrdf_close_atlas, &
        & cms_vn_atlas_version => vn_atlas_version

  USE mod_rttov_brdf_atlas, ONLY : &
        & vn_atlas_version,        &
        & vn_atlas_init
!INTF_ON
  IMPLICIT NONE

  TYPE(rttov_coefs), INTENT(IN) :: coefs ! RTTOV instrument coefficients

!INTF_END

 IF (coefs%coef%id_sensor == sensor_id_ir) THEN

    IF (vn_atlas_init) THEN
      CALL rttov_visnirbrdf_close_atlas
      vn_atlas_init = .FALSE.
    END IF

 END IF

END SUBROUTINE rttov_deallocate_brdf_atlas
