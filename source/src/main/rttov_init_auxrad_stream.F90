SUBROUTINE rttov_init_auxrad_stream(auxrad_stream)
! Description:
!   Initialise auxiliary radiance structure for streams
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
  USE rttov_types, ONLY : radiance_aux
!INTF_OFF
  USE parkind1, ONLY : jprb
!INTF_ON
  IMPLICIT NONE
  TYPE(radiance_aux ), INTENT(INOUT) :: auxrad_stream
!INTF_END

  IF (Associated (auxrad_stream%up))                 auxrad_stream%up                 = 0._jprb
  IF (Associated (auxrad_stream%down))               auxrad_stream%down               = 0._jprb
  IF (Associated (auxrad_stream%up_solar))           auxrad_stream%up_solar           = 0._jprb
  IF (Associated (auxrad_stream%down_solar))         auxrad_stream%down_solar         = 0._jprb
  IF (Associated (auxrad_stream%down_ref))           auxrad_stream%down_ref           = 0._jprb
  IF (Associated (auxrad_stream%down_ref_solar))     auxrad_stream%down_ref_solar     = 0._jprb
  IF (Associated (auxrad_stream%meanrad_up))         auxrad_stream%meanrad_up         = 0._jprb
  IF (Associated (auxrad_stream%meanrad_down))       auxrad_stream%meanrad_down       = 0._jprb
  IF (Associated (auxrad_stream%meanrad_up_solar))   auxrad_stream%meanrad_up_solar   = 0._jprb
  IF (Associated (auxrad_stream%meanrad_down_solar)) auxrad_stream%meanrad_down_solar = 0._jprb
  IF (Associated (auxrad_stream%cloudy))             auxrad_stream%cloudy             = 0._jprb
  IF (Associated (auxrad_stream%FAC1_2))             auxrad_stream%FAC1_2             = 0._jprb
  IF (Associated (auxrad_stream%FAC2_2))             auxrad_stream%FAC2_2             = 0._jprb
  IF (Associated (auxrad_stream%FAC3_2))             auxrad_stream%FAC3_2             = 0._jprb
  IF (Associated (auxrad_stream%FAC4_2))             auxrad_stream%FAC4_2             = 0._jprb
  IF (Associated (auxrad_stream%FAC5_2))             auxrad_stream%FAC5_2             = 0._jprb
  IF (Associated (auxrad_stream%FAC6_2))             auxrad_stream%FAC6_2             = 0._jprb
  IF (Associated (auxrad_stream%FAC7_2))             auxrad_stream%FAC7_2             = 0._jprb
  IF (Associated (auxrad_stream%FAC1_3))             auxrad_stream%FAC1_3             = 0._jprb
  IF (Associated (auxrad_stream%FAC2_3))             auxrad_stream%FAC2_3             = 0._jprb
  IF (Associated (auxrad_stream%FAC3_3))             auxrad_stream%FAC3_3             = 0._jprb
  IF (Associated (auxrad_stream%FAC4_3))             auxrad_stream%FAC4_3             = 0._jprb
  IF (Associated (auxrad_stream%FAC5_3))             auxrad_stream%FAC5_3             = 0._jprb
  IF (Associated (auxrad_stream%FAC6_3))             auxrad_stream%FAC6_3             = 0._jprb
  IF (Associated (auxrad_stream%FAC7_3))             auxrad_stream%FAC7_3             = 0._jprb

END SUBROUTINE 
