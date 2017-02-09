SUBROUTINE rttov_copy_raytracing(addsolar, raytracing1, raytracing2)
! Description:
!   Copy raytracing structure
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
  USE rttov_types, ONLY : raytracing_type
  USE parkind1, ONLY : jplm
  IMPLICIT NONE
  LOGICAL(jplm),         INTENT(IN)    :: addsolar
  TYPE(raytracing_type), INTENT(INOUT) :: raytracing1
  TYPE(raytracing_type), INTENT(IN)    :: raytracing2
!INTF_END

  ! Only copy those members which are used outside of locpat
  raytracing1%ZASAT    = raytracing2%ZASAT
  raytracing1%PATHSAT  = raytracing2%PATHSAT
  raytracing1%LTICK    = raytracing2%LTICK
  raytracing1%CO2_CM   = raytracing2%CO2_CM
  IF (addsolar) THEN
    raytracing1%ZASUN    = raytracing2%ZASUN
    raytracing1%PATHSUN  = raytracing2%PATHSUN
    raytracing1%PATHEFF  = raytracing2%PATHEFF
  ENDIF
END SUBROUTINE 
