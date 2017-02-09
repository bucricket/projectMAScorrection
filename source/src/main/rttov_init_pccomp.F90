! Description:
!> @file
!!   Initialise a pccomp structure.
!
!> @brief
!!   Initialise a pccomp structure.
!!
!! @param[in,out]  pccomp   PC-RTTOV pccomp structure to initialise
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
SUBROUTINE rttov_init_pccomp(pccomp)

  USE rttov_types, ONLY : rttov_pccomp
!INTF_OFF
  USE parkind1, ONLY : jprb
!INTF_ON

  IMPLICIT NONE

  TYPE(rttov_pccomp), INTENT(INOUT) :: pccomp

!INTF_END

  IF (ASSOCIATED(pccomp%pcscores))     pccomp%pcscores     = 0._jprb
  IF (ASSOCIATED(pccomp%bt_pccomp))    pccomp%bt_pccomp    = 0._jprb
  IF (ASSOCIATED(pccomp%total_pccomp)) pccomp%total_pccomp = 0._jprb
END SUBROUTINE 
