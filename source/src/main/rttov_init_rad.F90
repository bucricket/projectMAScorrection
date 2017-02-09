! Description:
!> @file
!!   Initialise a radiance structure.
!
!> @brief
!!   Initialise a radiance structure and optionally also a
!!   secondary radiance structure.
!!
!! @param[in,out]  rad   Radiance structure to initialise
!! @param[in,out]  rad2  Secondary radiance structure to initialise, optional
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
SUBROUTINE rttov_init_rad(rad, rad2)

  USE rttov_types, ONLY : radiance_type, radiance2_type
!INTF_OFF
  USE parkind1, ONLY : jprb
!INTF_ON

  IMPLICIT NONE

  TYPE(radiance_type),  INTENT(INOUT)           :: rad
  TYPE(radiance2_type), INTENT(INOUT), OPTIONAL :: rad2
!INTF_END

  rad%clear       = 0._jprb
  rad%total       = 0._jprb
  rad%bt_clear    = 0._jprb
  rad%bt          = 0._jprb
  rad%refl_clear  = 0._jprb
  rad%refl        = 0._jprb
  rad%cloudy      = 0._jprb
  rad%overcast    = 0._jprb

  IF (PRESENT(rad2)) THEN
    rad2%upclear     = 0._jprb
    rad2%dnclear     = 0._jprb
    rad2%refldnclear = 0._jprb
    rad2%up          = 0._jprb
    rad2%down        = 0._jprb
    rad2%surf        = 0._jprb
  ENDIF
END SUBROUTINE
