!
SUBROUTINE rttov_calcemis_ir_k( &
            & profiles,     &
            & profiles_k,   &
            & coef,         &
            & addpc,        &
            & coef_pccomp,  &
            & thermal,      &
            & chanprof,     &
            & calcemis,     &
            & emissivity_k)
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
  TYPE(profile_type     ), INTENT(INOUT) :: profiles_k(SIZE(chanprof))
  TYPE(profile_type     ), INTENT(IN)    :: profiles(:)
  TYPE(rttov_coef       ), INTENT(IN)    :: coef
  TYPE(rttov_coef_pccomp), INTENT(IN)    :: coef_pccomp
  LOGICAL(KIND=jplm)     , INTENT(IN)    :: thermal(SIZE(chanprof))
  LOGICAL(KIND=jplm)     , INTENT(IN)    :: calcemis(SIZE(chanprof))
  REAL(KIND=jprb)        , INTENT(INOUT) :: emissivity_k(SIZE(chanprof))
!INTF_END
  INTEGER(KIND=jpim) :: j, chan, iprof
  REAL   (KIND=jprb) :: windsp
  REAL   (KIND=jprb) :: windsp_k
  REAL   (KIND=jprb) :: aems
  REAL   (KIND=jprb) :: bems
  REAL   (KIND=jprb) :: cems
  REAL   (KIND=jprb) :: expf
  REAL   (KIND=jprb) :: fac
  REAL   (KIND=jprb) :: aems_k
  REAL   (KIND=jprb) :: bems_k
  REAL   (KIND=jprb) :: cems_k
  INTEGER(KIND=jpim) :: nchannels
  REAL   (KIND=jprb) :: ZHOOK_HANDLE
!- End of header --------------------------------------------------------
! Loop on all channels
  IF (LHOOK) CALL DR_HOOK('RTTOV_CALCEMIS_IR_K', 0_jpim, ZHOOK_HANDLE)
  nchannels = SIZE(chanprof)
  aems_k = 0._jprb
  bems_k = 0._jprb
  cems_k = 0._jprb
  windsp_k = 0._jprb
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
      aems   = coef_pccomp%emiss_c1(chan) + coef_pccomp%emiss_c2(chan) * windsp + &
               coef_pccomp%emiss_c3(chan) * windsp ** 2_jpim
      bems   = coef_pccomp%emiss_c4(chan) + coef_pccomp%emiss_c5(chan) * windsp + &
               coef_pccomp%emiss_c6(chan) * windsp ** 2_jpim
      cems   = coef_pccomp%emiss_c7(chan) + coef_pccomp%emiss_c8(chan) * windsp
      expf   = EXP( ((coef_pccomp%emiss_c9(chan) - 60._jprb) ** 2_jpim - &
                    (profiles(iprof)%zenangle - coef_pccomp%emiss_c9(chan)) ** 2_jpim) / cems )
      fac    = - ((coef_pccomp%emiss_c9(chan) - 60._jprb) ** 2_jpim -    &
                     (profiles(iprof)%zenangle - coef_pccomp%emiss_c9(chan)) ** 2_jpim) / cems ** 2_jpim
      aems_k = aems_k + emissivity_k(j)
      bems_k = bems_k + emissivity_k(j) * expf
      cems_k = cems_k + bems * expf * emissivity_k(j) * fac
      aems_k = aems_k - emissivity_k(j) * expf
      cems_k = cems_k - aems * expf * emissivity_k(j) * fac
      windsp_k  = windsp_k + coef_pccomp%emiss_c8(chan) * cems_k
      windsp_k  = windsp_k + (coef_pccomp%emiss_c5(chan) + &
                              2._jprb * coef_pccomp%emiss_c6(chan) * windsp) * bems_k
      windsp_k  = windsp_k + (coef_pccomp%emiss_c2(chan) + &
                              2._jprb * coef_pccomp%emiss_c3(chan) * windsp) * aems_k
      aems_k = 0._jprb
      bems_k = 0._jprb
      cems_k = 0._jprb
      IF (profiles(iprof)%s2m%u /= 0._jprb .OR. profiles(iprof)%s2m%v /= 0._jprb) THEN
        profiles_k(j)%s2m%u = profiles_k(j)%s2m%u + profiles(iprof)%s2m%u * windsp_k / windsp
        profiles_k(j)%s2m%v = profiles_k(j)%s2m%v + profiles(iprof)%s2m%v * windsp_k / windsp
      ELSE
        profiles_k(j)%s2m%u = profiles_k(j)%s2m%u + windsp_k / SQRT(2._jprb)
        profiles_k(j)%s2m%v = profiles_k(j)%s2m%v + windsp_k / SQRT(2._jprb)
      ENDIF
      windsp_k            = 0._jprb
    END IF
  ENDDO
  IF (LHOOK) CALL DR_HOOK('RTTOV_CALCEMIS_IR_K', 1_jpim, ZHOOK_HANDLE)
END SUBROUTINE rttov_calcemis_ir_k
