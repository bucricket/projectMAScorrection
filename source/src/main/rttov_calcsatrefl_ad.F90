!
SUBROUTINE rttov_calcsatrefl_ad(chanprof, profiles, solar_spectrum, thermal, solar, rad_ad)
! Description:
! To convert an array of radiances in many channels
! to reflectance values.
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
!    Copyright 2011, EUMETSAT, All Rights Reserved.
!
! Method: Uses ToA solar irradiances for each channel
!         read in from RT coefficient file.
!
! Current Code Owner: SAF NWP
!
! History:
! Version   Date     Comment
! -------   ----     -------
! 1.0    07/10/2011  Created (J Hocking)
!
! Code Description:
!   Language:           Fortran 90.
!   Software Standards: "European Standards for Writing and
!     Documenting Exchangeable Fortran 90 Code".
!
! Declarations:
! Modules used:
! Imported Type Definitions:
  USE rttov_types, ONLY : rttov_chanprof, profile_type, radiance_type
  USE parkind1, ONLY : jprb, jplm
!INTF_OFF
  USE rttov_const, ONLY : pi_r, deg2rad
  USE parkind1, ONLY : jpim
!INTF_ON
  IMPLICIT NONE
!subroutine arguments:
  TYPE(rttov_chanprof), INTENT(IN)    :: chanprof(:)
  TYPE(profile_type),   INTENT(IN)    :: profiles(:)
  REAL(KIND=jprb),      INTENT(IN)    :: solar_spectrum(size(chanprof))
  LOGICAL(KIND=jplm),   INTENT(IN)    :: thermal(size(chanprof))
  LOGICAL(KIND=jplm),   INTENT(IN)    :: solar(size(chanprof))
  TYPE(radiance_type),  INTENT(INOUT) :: rad_ad
!INTF_END
! radiances are expressed in mw/cm-1/ster/sq.m
! reflectances have no unit
!local variables:
  REAL   (KIND=jprb) :: toarad_norm
  INTEGER(KIND=jpim) :: nchanprof
  INTEGER(KIND=jpim) :: i, prof
!- End of header ------------------------------------------------------
  nchanprof = size(chanprof)

  DO i = 1, nchanprof
    IF (solar(i) .AND. .NOT. thermal(i)) THEN
      prof = chanprof(i)%prof

      toarad_norm = solar_spectrum(i) * pi_r * cos(profiles(prof)%sunzenangle * deg2rad)

      rad_ad%clear(i) = rad_ad%refl_clear(i) / toarad_norm
      rad_ad%total(i) = rad_ad%refl(i) / toarad_norm
    ENDIF
  ENDDO

END SUBROUTINE rttov_calcsatrefl_ad
