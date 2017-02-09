! Description:
!> @file
!!   Initialise a transmission structure.
!
!> @brief
!!   Initialise a transmission structure.
!!
!! @param[in,out]  transmission   Transmission structure to initialise
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
SUBROUTINE rttov_init_transmission(transmission)

  USE rttov_types, ONLY : transmission_type
!INTF_OFF
  USE parkind1, ONLY : jprb
!INTF_ON

  IMPLICIT NONE

  TYPE(transmission_type), INTENT(INOUT) :: transmission
!INTF_END

  transmission%tau_total  = 0._jprb
  transmission%tau_levels = 0._jprb
  transmission%tausun_total_path2  = 0._jprb
  transmission%tausun_levels_path2 = 0._jprb
  transmission%tausun_total_path1  = 0._jprb
  transmission%tausun_levels_path1 = 0._jprb
END SUBROUTINE 
