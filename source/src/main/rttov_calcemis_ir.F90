!
SUBROUTINE rttov_calcemis_ir( &
            & profiles,    &
            & geometry,    &
            & coef,        &
            & addpc,       &
            & coef_pccomp, &
            & thermal,     &
            & chanprof,    &
            & calcemis,    &
            & emissivity,  &
            & err)
! Description:
! To compute IR surface emissivities for all channels and all
! profiles if desired
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
! Method:
! RTTOV-6 IR surface emissivity report, V. Sherlock at:
! http://www.metoffice.com/research/interproj/nwpsaf/rtm/papers/isem6.pdf
!
! Current Code Owner: SAF NWP
!
! History:
! Version   Date     Comment
! -------   ----     -------
!  1.0       01/12/2002  New F90 code with structures (P Brunel A Smith)
!  1.1       02/01/2003  Comments added (R Saunders)
!  1.2       29/03/2005  Add end of header comment (J. Cameron)
!  1.3       02/12/2009  Introduced RTIASI emissivity model to be used with
!                        principal components. Marco Matricardi. ECMWF
!  1.4       18/10/2010  Always compute emissivity for PC-RTTOV regardless
!                        of calcemis (J Hocking)
!
! 2010/03 Code cleaning Pascal Brunel, Philippe Marguinaud
!
! Code Description:
!   Language:           Fortran 90.
!   Software Standards: "European Standards for Writing and
!     Documenting Exchangeable Fortran 90 Code".
!
!INTF_OFF
#include "throw.h"
!INTF_ON
  USE rttov_types, ONLY :  &
       & rttov_chanprof,    &
       & rttov_coef,        &
       & rttov_coef_pccomp, &
       & profile_type,      &
       & geometry_type
  USE parkind1, ONLY : jprb, jpim, jplm
!INTF_OFF
  USE rttov_const, ONLY :     &
       & surftype_sea,        &
       & surftype_land,       &
       & surftype_seaice,     &
       & errorstatus_success, &
       & errorstatus_fatal
  USE yomhook, ONLY : LHOOK, DR_HOOK
!INTF_ON
  IMPLICIT NONE
  LOGICAL(KIND=jplm)     , INTENT(IN)    :: addpc
  TYPE(rttov_chanprof   ), INTENT(IN)    :: chanprof(:)
  TYPE(profile_type     ), INTENT(IN)    :: profiles(:)
  TYPE(geometry_type    ), INTENT(IN)    :: geometry(:)
  TYPE(rttov_coef       ), INTENT(IN)    :: coef
  TYPE(rttov_coef_pccomp), INTENT(IN)    :: coef_pccomp
  LOGICAL(KIND=jplm)     , INTENT(IN)    :: thermal(SIZE(chanprof))
  LOGICAL(KIND=jplm)     , INTENT(IN)    :: calcemis(SIZE(chanprof))
  REAL(KIND=jprb)        , INTENT(INOUT) :: emissivity(SIZE(chanprof))
  INTEGER(KIND=jpim)     , INTENT(OUT)   :: err
!INTF_END
  INTEGER(KIND=jpim) :: j, chan, iprof
  REAL   (KIND=jprb) :: windsp
  REAL   (KIND=jprb) :: aems
  REAL   (KIND=jprb) :: bems
  REAL   (KIND=jprb) :: cems
  INTEGER(KIND=jpim) :: nchannels
  REAL   (KIND=jprb) :: ZHOOK_HANDLE
!- End of header --------------------------------------------------------
  TRY
  IF (LHOOK) CALL DR_HOOK('RTTOV_CALCEMIS_IR', 0_jpim, ZHOOK_HANDLE)
  nchannels = SIZE(chanprof)
  DO j = 1, nchannels

    IF (.NOT. thermal(j)) CYCLE
    IF (.NOT. calcemis(j)) CYCLE

    chan  = chanprof(j)%chan
    iprof = chanprof(j)%prof

    IF (addpc .AND. profiles(iprof)%skin%surftype == surftype_sea) THEN

      windsp  = SQRT(profiles(iprof)%s2m%u ** 2 + profiles(iprof)%s2m%v ** 2)
      aems    = coef_pccomp%emiss_c1(chan) + &
                coef_pccomp%emiss_c2(chan) * windsp + &
                coef_pccomp%emiss_c3(chan) * windsp ** 2_jpim
      bems    = coef_pccomp%emiss_c4(chan) + &
                coef_pccomp%emiss_c5(chan) * windsp + &
                coef_pccomp%emiss_c6(chan) * windsp ** 2_jpim
      cems    = coef_pccomp%emiss_c7(chan) + coef_pccomp%emiss_c8(chan) * windsp
      emissivity(j) = aems + (bems - aems) * &
                EXP(((coef_pccomp%emiss_c9(chan) - 60._jprb) ** 2_jpim - &
                     (profiles(iprof)%zenangle - coef_pccomp%emiss_c9(chan)) ** 2_jpim) / cems)

    ELSE ! Not using PC-RTTOV sea surface emissivity model

      !-----------------------------------------
      !1. Use a fixed value over land and seaice
      !-----------------------------------------
      IF (profiles(iprof)%skin%surftype == surftype_land) THEN
        emissivity(j) = 0.98_jprb
      ELSE IF (profiles(iprof)%skin%surftype == surftype_seaice) THEN
        emissivity(j) = 0.99_jprb
      ELSE
      !------------------------------------------------------------------
      !2. Over sea, emissivity is a polynomial in normalized zenith angle
      ! ISEM6 model
      !------------------------------------------------------------------
        emissivity(j) = coef%ssirem_a0(chan) - &
                        coef%ssirem_a1(chan) * geometry(iprof)%normzen ** coef%ssirem_xzn1(chan) -  &
                        coef%ssirem_a2(chan) * geometry(iprof)%normzen ** coef%ssirem_xzn2(chan)
      ENDIF
    ENDIF
  ENDDO
  IF (LHOOK) CALL DR_HOOK('RTTOV_CALCEMIS_IR', 1_jpim, ZHOOK_HANDLE)
  CATCH
  IF (LHOOK) CALL DR_HOOK('RTTOV_CALCEMIS_IR', 1_jpim, ZHOOK_HANDLE)
END SUBROUTINE rttov_calcemis_ir
