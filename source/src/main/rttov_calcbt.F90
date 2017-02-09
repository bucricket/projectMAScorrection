!
SUBROUTINE rttov_calcbt(chanprof, coeffs, thermal, rad)
! Description:
! To convert an array of radiances in many channels
! to equivalent black body brightness temperatures.
! Planck function is applied with a "central wave number"
! Temperature is corrected by "band corrections" coefficients
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
! 1.3    28/02/2004  Improved vectorisation (D Dent)
! 1.4    26/01/2007  Removed polarisation (R Saunders)
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
  USE rttov_math_mod, ONLY : inv_planck
!INTF_ON
  IMPLICIT NONE
!subroutine arguments:
  TYPE(rttov_chanprof), INTENT(IN)    :: chanprof(:)! Array of channel indices.
  TYPE(rttov_coef    ), INTENT(IN)    :: coeffs     ! Coefficients
  LOGICAL(KIND=jplm  ), INTENT(IN)    :: thermal(size(chanprof))
  TYPE(radiance_type ), INTENT(INOUT) :: rad        ! input radiances and output BT
!INTF_END
! radiances are expressed in mw/cm-1/ster/sq.m
! and temperatures in Kelvin
!local variables:
  REAL   (KIND=jprb) :: tstore1, tstore2, radtotal, radclear
  INTEGER(KIND=jpim) :: nchannels
  INTEGER(KIND=jpim) :: chan, i, pol_id
!- End of header ------------------------------------------------------
  nchannels = SIZE(chanprof)
  IF (coeffs%id_sensor /= sensor_id_po) THEN! All no polarimetric sensors
    DO i = 1, nchannels
      IF (thermal(i)) THEN
        chan            = chanprof(i)%chan
        CALL inv_planck(coeffs%planck1(chan), coeffs%planck2(chan), rad%total(i), tstore1)
        CALL inv_planck(coeffs%planck1(chan), coeffs%planck2(chan), rad%clear(i), tstore2)
        rad%bt(i)       = (tstore1 - coeffs%ff_bco(chan)) / coeffs%ff_bcs(chan)
        rad%bt_clear(i) = (tstore2 - coeffs%ff_bco(chan)) / coeffs%ff_bcs(chan)
      ELSE
        rad%bt(i)       = 0._jprb
        rad%bt_clear(i) = 0._jprb
      ENDIF
    ENDDO
  ELSE! Special case for polarimetric radiometers
! Add average of 1st two elements of Stokes vector to differences in
! 3rd/4th before conversion to BT and subtract after conversion
    DO i = 1, nchannels
      IF (thermal(i)) THEN
        chan   = chanprof(i)%chan
        pol_id = coeffs%fastem_polar(chan) + 1_jpim
        IF (pol_id == 6_jpim) THEN! 3rd element
          radtotal        = rad%total(i) + 0.5_jprb * (rad%total(i - 2) + rad%total(i - 1))
          radclear        = rad%clear(i) + 0.5_jprb * (rad%clear(i - 2) + rad%clear(i - 1))
          CALL inv_planck(coeffs%planck1(chan), coeffs%planck2(chan), radtotal, tstore1)
          CALL inv_planck(coeffs%planck1(chan), coeffs%planck2(chan), radclear, tstore2)
          rad%bt(i)       = (tstore1 - coeffs%ff_bco(chan)) / coeffs%ff_bcs(chan)
          rad%bt_clear(i) = (tstore2 - coeffs%ff_bco(chan)) / coeffs%ff_bcs(chan)
          rad%bt(i)       = rad%bt(i) - 0.5_jprb * (rad%bt(i - 2) + rad%bt(i - 1))
          rad%bt_clear(i) = rad%bt_clear(i) - 0.5_jprb * (rad%bt_clear(i - 2) + rad%bt_clear(i - 1))
        ELSE IF (pol_id == 7_jpim) THEN! 4th element
          radtotal        = rad%total(i) + 0.5_jprb * (rad%total(i - 3) + rad%total(i - 2))
          radclear        = rad%clear(i) + 0.5_jprb * (rad%clear(i - 3) + rad%clear(i - 2))
          CALL inv_planck(coeffs%planck1(chan), coeffs%planck2(chan), radtotal, tstore1)
          CALL inv_planck(coeffs%planck1(chan), coeffs%planck2(chan), radclear, tstore2)
          rad%bt(i)       = (tstore1 - coeffs%ff_bco(chan)) / coeffs%ff_bcs(chan)
          rad%bt_clear(i) = (tstore2 - coeffs%ff_bco(chan)) / coeffs%ff_bcs(chan)
          rad%bt(i)       = rad%bt(i) - 0.5_jprb * (rad%bt(i - 3) + rad%bt(i - 2))
          rad%bt_clear(i) = rad%bt_clear(i) - 0.5_jprb * (rad%bt_clear(i - 3) + rad%bt_clear(i - 2))
        ELSE! 1st or 2nd element
          CALL inv_planck(coeffs%planck1(chan), coeffs%planck2(chan), rad%total(i), tstore1)
          CALL inv_planck(coeffs%planck1(chan), coeffs%planck2(chan), rad%clear(i), tstore2)
          rad%bt(i)       = (tstore1 - coeffs%ff_bco(chan)) / coeffs%ff_bcs(chan)
          rad%bt_clear(i) = (tstore2 - coeffs%ff_bco(chan)) / coeffs%ff_bcs(chan)
        ENDIF
      ELSE
        rad%bt(i)       = 0._jprb
        rad%bt_clear(i) = 0._jprb
      ENDIF
    ENDDO
  ENDIF
END SUBROUTINE rttov_calcbt
