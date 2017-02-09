SUBROUTINE rttov_scale_radiance_inc( &
            & radiance_inc, &
            & factor,    &
            & opts)
! Description:
!   Scales radiance or brightness temperature increment
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
!
  USE parkind1, ONLY : jprb
  USE rttov_types, ONLY : radiance_type, rttov_options
  IMPLICIT NONE
  TYPE(radiance_type), INTENT(INOUT) :: radiance_inc
  REAL(KIND=jprb)    , INTENT(IN)    :: factor
  TYPE(rttov_options), INTENT(IN)    :: opts
!INTF_END

  IF (opts%rt_all%switchrad) radiance_inc%bt = radiance_inc%bt * factor

  radiance_inc%total = radiance_inc%total * factor

END SUBROUTINE 
