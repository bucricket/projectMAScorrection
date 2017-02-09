SUBROUTINE rttov_init_sunglint(sunglint)
! Description:
!   Initialise sunglint structure
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
!    Copyright 2010, EUMETSAT, All Rights Reserved.
!
! Method:
!
! Current Code Owner: SAF NWP
!
! History:
! Version   Date     Comment
! -------   ----     -------
! 2010/03 Code cleaning Pascal Brunel, Philippe Marguinaud
!
!
! Code Description:
!   Language:           Fortran 90.
!   Software Standards: "European Standards for Writing and
!     Documenting Exchangeable Fortran 90 Code".
  USE rttov_types, ONLY : sunglint_type
!INTF_OFF
  USE parkind1, ONLY : jprb
!INTF_ON
  TYPE(sunglint_type), INTENT(INOUT) :: sunglint
!INTF_END
  sunglint%s%CSI         = 0._jprb
  sunglint%s%ALFA        = 0._jprb
  sunglint%s%C_SHAD      = 0._jprb
  sunglint%s%P_PRIME     = 0._jprb
  sunglint%s%PXY_GAMMAXY = 0._jprb
  sunglint%s%GAMMA_O     = 0._jprb
  sunglint%s%GAMMA_P     = 0._jprb
  sunglint%s%G_SHAD      = 0._jprb
  sunglint%s%GAMMAX      = 0._jprb
  sunglint%s%Q_SHAD      = 0._jprb
  sunglint%s%ZENSAT      = 0._jprb
  sunglint%s%ZENSUN      = 0._jprb
  sunglint%s%DAZNG       = 0._jprb
  sunglint%s%FAC1        = 0._jprb
  sunglint%s%A_SHAD      = 0._jprb
  sunglint%s%B_SHAD      = 0._jprb
  sunglint%s%LAMBDA_A    = 0._jprb
  sunglint%s%LAMBDA_B    = 0._jprb
  sunglint%s%X_U         = 0._jprb
  sunglint%s%ALFA1       = 0._jprb
  sunglint%s%OMEGA_M     = 0._jprb
  sunglint%s%WINDSP      = 0._jprb
  sunglint%s%WANGL       = 0._jprb
  sunglint%s%GAMMA_SQ    = 0._jprb
  sunglint%s%GLINT       = 0._jprb
  sunglint%s%OMEGA       = 0._jprb
  ! beta and psi are direct-model-only variables
!   sunglint%BETA          = 0._jprb
!   sunglint%PSI           = 0._jprb
END SUBROUTINE 
