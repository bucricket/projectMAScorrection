!
SUBROUTINE rttov_calcemis_ir_tl( &
            & profiles,      &
            & profiles_tl,   &
            & coef,          &
            & addpc,         &
            & coef_pccomp,   &
            & thermal,       &
            & chanprof,      &
            & calcemis,      &
            & emissivity_tl)
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
  USE rttov_types, ONLY :  &
       & rttov_chanprof,    &
       & rttov_coef,        &
       & rttov_coef_pccomp, &
       & profile_type
  USE parkind1, ONLY : jprb, jplm
!INTF_OFF
  USE rttov_const, ONLY : surftype_sea
  USE yomhook, ONLY : LHOOK, DR_HOOK
  USE parkind1, ONLY : jpim
!INTF_ON
  IMPLICIT NONE
  LOGICAL(KIND=jplm)     , INTENT(IN)    :: addpc
  TYPE(rttov_chanprof   ), INTENT(IN)    :: chanprof(:)
  TYPE(profile_type     ), INTENT(IN)    :: profiles_tl(:)
  TYPE(profile_type     ), INTENT(IN)    :: profiles(:)
  TYPE(rttov_coef       ), INTENT(IN)    :: coef
  TYPE(rttov_coef_pccomp), INTENT(IN)    :: coef_pccomp
  LOGICAL(KIND=jplm)     , INTENT(IN)    :: thermal(SIZE(chanprof))
  LOGICAL(KIND=jplm)     , INTENT(IN)    :: calcemis(SIZE(chanprof))
  REAL(KIND=jprb)        , INTENT(INOUT) :: emissivity_tl(SIZE(chanprof))
!INTF_END
  INTEGER(KIND=jpim) :: j, chan, iprof
  REAL   (KIND=jprb) :: windsp
  REAL   (KIND=jprb) :: windsp_tl
  REAL   (KIND=jprb) :: aems
  REAL   (KIND=jprb) :: bems
  REAL   (KIND=jprb) :: cems
  REAL   (KIND=jprb) :: expf
  REAL   (KIND=jprb) :: fac
  REAL   (KIND=jprb) :: aems_tl
  REAL   (KIND=jprb) :: bems_tl
  REAL   (KIND=jprb) :: cems_tl
  INTEGER(KIND=jpim) :: nchannels
  REAL   (KIND=jprb) :: ZHOOK_HANDLE
!- End of header --------------------------------------------------------
! Loop on all channels
  IF (LHOOK) CALL DR_HOOK('RTTOV_CALCEMIS_IR_TL', 0_jpim, ZHOOK_HANDLE)
  nchannels = SIZE(chanprof)
  DO j = 1, nchannels

    IF (.NOT. thermal(j)) CYCLE
    IF (.NOT. calcemis(j)) CYCLE

    chan  = chanprof(j)%chan
    iprof = chanprof(j)%prof

    IF (addpc .AND. profiles(iprof)%skin%surftype == surftype_sea) THEN
      !------------------------------------------------------------------
      !1. Over sea, emissivity is computed using RTIASI model
      !------------------------------------------------------------------
      windsp = SQRT(profiles(iprof)%s2m%u ** 2_jpim + profiles(iprof)%s2m%v ** 2_jpim)
      IF (profiles(iprof)%s2m%u /= 0._jprb .OR. profiles(iprof)%s2m%v /= 0._jprb) THEN
        windsp_tl = (profiles_tl(iprof)%s2m%u * profiles(iprof)%s2m%u + &
                     profiles_tl(iprof)%s2m%v * profiles(iprof)%s2m%v) / windsp
      ELSE
        windsp_tl = (profiles_tl(iprof)%s2m%u + profiles_tl(iprof)%s2m%v) / SQRT(2._jprb)
      ENDIF
      aems    = coef_pccomp%emiss_c1(chan) + coef_pccomp%emiss_c2(chan) * windsp + &
                coef_pccomp%emiss_c3(chan) * windsp ** 2_jpim
      aems_tl = coef_pccomp%emiss_c2(chan) * windsp_tl + &
                2._jprb * coef_pccomp%emiss_c3(chan) * windsp * windsp_tl
      bems    = coef_pccomp%emiss_c4(chan) + coef_pccomp%emiss_c5(chan) * windsp + &
                coef_pccomp%emiss_c6(chan) * windsp ** 2_jpim
      bems_tl = coef_pccomp%emiss_c5(chan) * windsp_tl + &
                2._jprb * coef_pccomp%emiss_c6(chan) * windsp * windsp_tl
      cems    = coef_pccomp%emiss_c7(chan) + coef_pccomp%emiss_c8(chan) * windsp
      cems_tl = coef_pccomp%emiss_c8(chan) * windsp_tl
      expf    = EXP( ((coef_pccomp%emiss_c9(chan) - 60._jprb) ** 2_jpim - &
                           (profiles(iprof)%zenangle - coef_pccomp%emiss_c9(chan)) ** 2_jpim) / cems )
      fac     = - ((coef_pccomp%emiss_c9(chan) - 60._jprb) ** 2_jpim - &
                   (profiles(iprof)%zenangle - coef_pccomp%emiss_c9(chan)) ** 2_jpim) / cems ** 2_jpim
      emissivity_tl(j) = aems_tl + bems_tl * expf - aems_tl * expf + &
                         cems_tl * bems * expf * fac - cems_tl * aems * expf * fac

    ELSE ! Not using PC-RTTOV sea surface emissivity model

! In the manual, it is written that the user give RTTOV an emissivity variation "around"
! the calculated value
!         emissivity_tl(j) = 0._jprb
    ENDIF
  ENDDO
  IF (LHOOK) CALL DR_HOOK('RTTOV_CALCEMIS_IR_TL', 1_jpim, ZHOOK_HANDLE)
END SUBROUTINE rttov_calcemis_ir_tl
