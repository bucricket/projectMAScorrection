!
SUBROUTINE rttov_calc_solar_spec_esd(coef, chanprof, profiles, solar_spec_esd)
! Description:
! To calculate top-of-atmosphere solar spectral irradiance
!   adjusted for Earth-sun distance based on profile date.
!   If the date is not valid the distance is taken as 1 AU.
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
! 1.0    12/10/2011  Created (J Hocking)
!
! Code Description:
!   Language:           Fortran 90.
!   Software Standards: "European Standards for Writing and
!     Documenting Exchangeable Fortran 90 Code".
!
! Declarations:
! Modules used:
! Imported Type Definitions:
  USE rttov_types, ONLY : rttov_chanprof, rttov_coef, profile_type
  USE parkind1, ONLY : jprb
!INTF_OFF
  USE rttov_const, ONLY : pi
  USE parkind1, ONLY : jpim
!INTF_ON
  IMPLICIT NONE
!subroutine arguments:
  TYPE(rttov_coef),     INTENT(IN)    :: coef                           ! Coefficients
  TYPE(rttov_chanprof), INTENT(IN)    :: chanprof(:)                    ! Array of channel indices
  TYPE(profile_type),   INTENT(IN)    :: profiles(:)                    ! Profiles
  REAL(KIND=jprb),      INTENT(OUT)   :: solar_spec_esd(size(chanprof)) ! Output irradiances
!INTF_END
!local variables:
  INTEGER(KIND=jpim) :: nchanprof, nprofiles
  INTEGER(KIND=jpim) :: i, chan, prof
  INTEGER(KIND=jpim) :: year, month, day, day_of_year
  INTEGER(KIND=jpim) :: days(12)
  REAL(KIND=jprb)    :: d
  REAL(KIND=jprb)    :: esdsq_r(size(profiles))  ! Reciprocal of Earth-sun distance (in AU) squared
!- End of header ------------------------------------------------------
  nchanprof = size(chanprof)
  nprofiles = size(profiles)

  days = (/31,28,31,30,31,30,31,31,30,31,30,31/)

  DO prof = 1, nprofiles

    year  = profiles(prof) % date(1)
    month = profiles(prof) % date(2)
    day   = profiles(prof) % date(3)

    ! If the date is not valid, assume Earth-sun distance is 1 AU.
    esdsq_r(prof) = 1.0_jprb

    IF (year > 1950) THEN  ! Default profile year is 1950
      IF (month > 0 .AND. month < 13) THEN

        IF (MOD(year, 4_jpim) == 0 .AND. (.NOT. MOD(year, 100_jpim) == 0 .OR. MOD(year, 400_jpim) == 0)) THEN
          days(2) = 29
        ELSE
          days(2) = 28
        ENDIF

        IF (day > 0 .AND. day <= days(month)) THEN

          day_of_year = SUM(days(1:month-1)) + day - 1

          ! JAH - this approximation is common in textbooks/literature.
          d             = 2.0_jprb * pi * day_of_year/365.0_jprb
          esdsq_r(prof) = 1.00011_jprb + 0.034221_jprb * cos(d) + 0.00128_jprb * sin(d) + &
                          0.000719_jprb * cos(d*2.0_jprb) + 0.000077_jprb * sin(d*2.0_jprb)
        ENDIF
      ENDIF
    ENDIF
  ENDDO

  DO i = 1, nchanprof
    chan = chanprof(i)%chan
    prof = chanprof(i)%prof

    solar_spec_esd(i) = coef%ss_solar_spectrum(chan) * esdsq_r(prof)
  ENDDO

END SUBROUTINE rttov_calc_solar_spec_esd
