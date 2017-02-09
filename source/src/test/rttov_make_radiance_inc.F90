SUBROUTINE rttov_make_radiance_inc( &
            & coef,         &
            & radiance_inc, &
            & channels,     &
            & nchannels,    &
            & opts)
! Description:
!   Computes a sensible radiance variation, either
!   in brightness temperature or in radiance
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
  USE parkind1, ONLY : jpim
  USE rttov_types, ONLY : rttov_coef, radiance_type, rttov_options
!INTF_OFF
  USE parkind1, ONLY : jprb
  USE rttov_const, ONLY : c1, c2
  USE rttov_math_mod, ONLY : planck, planck_tl
!INTF_ON
  IMPLICIT NONE
  TYPE(rttov_coef   ), INTENT(IN)    :: coef
  TYPE(radiance_type), INTENT(INOUT) :: radiance_inc
  INTEGER(KIND=jpim) , INTENT(IN)    :: nchannels
  INTEGER(KIND=jpim) , INTENT(IN)    :: channels(nchannels)
  TYPE(rttov_options), INTENT(IN)    :: opts
!INTF_END
  INTEGER(KIND=jpim) :: ichan
  REAL(KIND=jprb)    :: b

  IF (opts%rt_all%switchrad) radiance_inc%bt = .1_jprb

  DO ichan = 1, nchannels
    IF (coef%ss_val_chn(channels(ichan)) < 2) THEN
      CALL planck(coef%ff_cwn(channels(ichan)), c1, c2, 280._jprb, b)
      CALL planck_tl(coef%ff_cwn(channels(ichan)), c1, c2, 280._jprb, 0.1_jprb, b, radiance_inc%total(ichan))
    ELSE
      radiance_inc%total(ichan) = .1_jprb
    ENDIF
  ENDDO

END SUBROUTINE
