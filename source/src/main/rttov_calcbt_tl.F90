SUBROUTINE rttov_calcbt_tl( &
            & chanprof, &
            & coeffs,   &
            & thermal,  &
            & rad,      &
            & rad_tl)
!
! Description:
! To convert an array of radiances in many channels
! to equivalent black body brightness temperatures.
! Planck function is applied with a "central wave number"
! Temperature is corrected by "band corrections" coefficients
! derivative of inverse Planck function with respect to radiance is
!
!                             C1 * C2 * Nu**4
! B-1'(R,Nu) = --------------------------------------------- dR
!                     (    C1 * Nu**3) (     C1 * Nu**3 )**2
!               R**2 *(1 + ----------) ( Ln( ---------- )
!                     (         R    ) (        R       )
!
! which can be reduced to the following, with
!  C1 = C1 * Nu**3
!  C2 = C2 * Nu
!
!                  C1 * B-1(R,Nu)**2
! B-1'(R,Nu) = ----------------------- dR
!                  C2 * R * (R + C1)
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
!    Copyright 2002, EUMETSAT, All Rights Reserved.
!
! Method: Uses band correction coefficients for each channel
!         read in from RT coefficient file.
!
! Current Code Owner: SAF NWP
!
! History:
! Version   Date     Comment
! -------   ----     -------
! 1.0    01/12/2002  New F90 code with structures (P Brunel A Smith)
!                      based on BRIGV of previous RTTOV versions
! 1.1    02/01/2003  Comments added (R Saunders)
! 1.2    08/01/2004  Polarisation added (S English)
! 1.3    29/03/2005  Add end of header comment (J. Cameron)
! 1.4    05/02/2007  Remove polarisation (R Saunders)
! 1.5    09/11/2010  Bug fix for pol sensors (J Hocking)
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
! Imported Type Definitions:
  USE rttov_types, ONLY : rttov_chanprof, rttov_coef, radiance_type
  USE parkind1, ONLY : jplm
!INTF_OFF
  USE rttov_const, ONLY : sensor_id_po
  USE parkind1, ONLY : jprb, jpim
  USE rttov_math_mod, ONLY : inv_planck_tl
!INTF_ON
  IMPLICIT NONE
!subroutine arguments:
  TYPE(rttov_chanprof), INTENT(IN)    :: chanprof(:)
  TYPE(rttov_coef    ), INTENT(IN)    :: coeffs
  LOGICAL(KIND=jplm  ), INTENT(IN)    :: thermal(size(chanprof))
  TYPE(radiance_type ), INTENT(IN)    :: rad        ! rad%total rad%clear
  TYPE(radiance_type ), INTENT(INOUT) :: rad_tl
! input  rad_tl%total  and rad_tl%clear
! output rad_tl%bt     and rad_tl%bt_clear
!INTF_END
!local variables:
  REAL   (KIND=jprb) :: tstar    , tstar1   , tstar2
  REAL   (KIND=jprb) :: tstar_tl , tstar1_tl, tstar2_tl
  REAL   (KIND=jprb) :: radtotal , radclear
  INTEGER(KIND=jpim) :: chan     , i, pol_id
  INTEGER(KIND=jpim) :: nchannels
!- End of header --------------------------------------------------------
  nchannels = SIZE(chanprof)
  IF (coeffs%id_sensor /= sensor_id_po) THEN! All no polarimetric sensors
    DO i = 1, nchannels
      IF (thermal(i)) THEN
        chan = chanprof(i)%chan
! Clear+cloudy radiance
! T star for direct model
        tstar              = coeffs%ff_bco(chan) + coeffs%ff_bcs(chan) * rad%bt(i)
! TL
        CALL inv_planck_tl(coeffs%planck1(chan), coeffs%planck2(chan), rad%total(i), rad_tl%total(i), tstar, tstar_tl)
        rad_tl%bt(i)       = tstar_tl / coeffs%ff_bcs(chan)
!clear radiance
! T star for direct model
        tstar              = coeffs%ff_bco(chan) + coeffs%ff_bcs(chan) * rad%bt_clear(i)
! TL
        CALL inv_planck_tl(coeffs%planck1(chan), coeffs%planck2(chan), rad%clear(i), rad_tl%clear(i), tstar, tstar_tl)
        rad_tl%bt_clear(i) = tstar_tl / coeffs%ff_bcs(chan)
      ELSE
        rad_tl%bt(i)       = 0._jprb
        rad_tl%bt_clear(i) = 0._jprb
      ENDIF
    ENDDO
  ELSE! Special case for polarimetric radiometers
! Add average of 1st two elements of Stokes vector to differences in
! 3rd/4th before conversion to BT and subtract after conversion
    DO i = 1, nchannels
      IF (thermal(i)) THEN
        chan   = chanprof(i)%chan
        pol_id = coeffs%fastem_polar(chan) + 1_jpim
        IF (pol_id == 6_jpim) THEN! 3rd stoke element
          tstar1             = coeffs%ff_bco(chan) + coeffs%ff_bcs(chan) * rad%bt(i)
          tstar1             = tstar1 + coeffs%ff_bcs(chan) * 0.5_jprb * (rad%bt(i - 2) + rad%bt(i - 1))
          tstar2             = coeffs%ff_bco(chan) + coeffs%ff_bcs(chan) * rad%bt_clear(i)
          tstar2             = tstar2 + coeffs%ff_bcs(chan) * 0.5_jprb * (rad%bt_clear(i - 2) + rad%bt_clear(i - 1))
          radtotal           = rad%total(i) + 0.5 * (rad%total(i - 2) + rad%total(i - 1))
          radclear           = rad%clear(i) + 0.5 * (rad%clear(i - 2) + rad%clear(i - 1))
          CALL inv_planck_tl(coeffs%planck1(chan), coeffs%planck2(chan), radtotal, rad_tl%total(i), tstar1, tstar1_tl)
          rad_tl%bt(i)       = tstar1_tl / coeffs%ff_bcs(chan)
          CALL inv_planck_tl(coeffs%planck1(chan), coeffs%planck2(chan), radclear, rad_tl%clear(i), tstar2, tstar2_tl)
          rad_tl%bt_clear(i) = tstar2_tl / coeffs%ff_bcs(chan)
        ELSE IF (pol_id == 7_jpim) THEN! 4th stoke element
          tstar1             = coeffs%ff_bco(chan) + coeffs%ff_bcs(chan) * rad%bt(i)
          tstar1             = tstar1 + coeffs%ff_bcs(chan) * 0.5_jprb * (rad%bt(i - 3) + rad%bt(i - 2))
          tstar2             = coeffs%ff_bco(chan) + coeffs%ff_bcs(chan) * rad%bt_clear(i)
          tstar2             = tstar2 + coeffs%ff_bcs(chan) * 0.5_jprb * (rad%bt_clear(i - 3) + rad%bt_clear(i - 2))
          radtotal           = rad%total(i) + 0.5 * (rad%total(i - 3) + rad%total(i - 2))
          radclear           = rad%clear(i) + 0.5 * (rad%clear(i - 3) + rad%clear(i - 2))
          CALL inv_planck_tl(coeffs%planck1(chan), coeffs%planck2(chan), radtotal, rad_tl%total(i), tstar1, tstar1_tl)
          rad_tl%bt(i)       = tstar1_tl / coeffs%ff_bcs(chan)
          CALL inv_planck_tl(coeffs%planck1(chan), coeffs%planck2(chan), radclear, rad_tl%clear(i), tstar2, tstar2_tl)
          rad_tl%bt_clear(i) = tstar2_tl / coeffs%ff_bcs(chan)
        ELSE! 1st/2nd stoke elements
! Clear+cloudy radiance
! T star for direct model
          tstar              = coeffs%ff_bco(chan) + coeffs%ff_bcs(chan) * rad%bt(i)
! TL
          CALL inv_planck_tl(coeffs%planck1(chan), coeffs%planck2(chan), rad%total(i), rad_tl%total(i), tstar, tstar_tl)
          rad_tl%bt(i)       = tstar_tl / coeffs%ff_bcs(chan)
!clear radiance
! T star for direct model
          tstar              = coeffs%ff_bco(chan) + coeffs%ff_bcs(chan) * rad%bt_clear(i)
! TL
          CALL inv_planck_tl(coeffs%planck1(chan), coeffs%planck2(chan), rad%clear(i), rad_tl%clear(i), tstar, tstar_tl)
          rad_tl%bt_clear(i) = tstar_tl / coeffs%ff_bcs(chan)
        ENDIF
      ELSE
        rad_tl%bt(i)       = 0._jprb
        rad_tl%bt_clear(i) = 0._jprb
      ENDIF
    ENDDO
  ENDIF
END SUBROUTINE rttov_calcbt_tl
