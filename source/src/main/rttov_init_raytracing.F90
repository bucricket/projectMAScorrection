SUBROUTINE rttov_init_raytracing(addsolar, raytracings)
! Description:
!   Initialise raytracing structure
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
  USE rttov_types, ONLY : raytracing_Type
  USE parkind1, ONLY : jplm
!INTF_OFF
  USE parkind1, ONLY : jprb
!INTF_ON
  IMPLICIT NONE
  LOGICAL(jplm),         INTENT(IN)    :: addsolar
  TYPE(raytracing_Type), INTENT(INOUT) :: raytracings
!INTF_END
  raytracings%LTICK    = 0._jprb
  raytracings%HGPL     = 0._jprb
  raytracings%DMAIR_R  = 0._jprb
  raytracings%ZASAT    = 0._jprb
  raytracings%PATHSAT  = 0._jprb
  raytracings%CO2_CM   = 0._jprb
  IF (addsolar) THEN
    raytracings%ZASUN    = 0._jprb
    raytracings%PATHSUN  = 0._jprb
    raytracings%PATHEFF  = 0._jprb
  ENDIF
END SUBROUTINE 
