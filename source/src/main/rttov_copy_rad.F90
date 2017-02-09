!
SUBROUTINE rttov_copy_rad( radiance1, radiance2 )
! Description:
! allocation/deallocation of a radiance structure
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
!    Copyright 2007, EUMETSAT, All Rights Reserved.
!
! Method:
!
! Current Code Owner: SAF NWP
!
! History:
! Version   Date     Comment
! -------   ----     -------
!  1.0       01/12/2002  New F90 code
!  2.0       02/12/2009  Marco matricardi:Added Principal component capability
!
! 2010/03 Code cleaning Pascal Brunel, Philippe Marguinaud
!
!
! Code Description:
!   Language:           Fortran 90.
!   Software Standards: "European Standards for Writing and
!     Documenting Exchangeable Fortran 90 Code".
!
! Declarations:
! Modules used:
! Imported Parameters:
! Imported Type Definitions:
  USE rttov_types, ONLY : radiance_type
  IMPLICIT NONE
! subroutine arguments
! scalar arguments with intent(out):
  TYPE(radiance_Type), INTENT(INOUT) :: radiance1 ! radiances
  TYPE(radiance_Type), INTENT(IN)    :: radiance2 ! radiances
!INTF_END
!
  radiance1%clear       = radiance2%clear
  radiance1%total       = radiance2%total
  radiance1%bt_clear    = radiance2%bt_clear
  radiance1%bt          = radiance2%bt
  radiance1%refl_clear  = radiance2%refl_clear
  radiance1%refl        = radiance2%refl
  radiance1%cloudy      = radiance2%cloudy
  radiance1%overcast    = radiance2%overcast
END SUBROUTINE
