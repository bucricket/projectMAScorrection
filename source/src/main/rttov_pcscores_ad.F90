!    Compute principal component scores
SUBROUTINE rttov_pcscores_ad( &
            & opts,            &
            & chanprof,        &
            & chanprof_pc,     &
            & pccomp,          &
            & pccomp_ad,       &
            & coef_pccomp,     &
            & radiancedata_ad)
!     Description:
!     To compute principal component scores
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
!     Method:
!     See: User's manual and scientific report for PC_RTTOV
!          (Available from EUMETSAT)
!
!     Owner:
!     EUMETSAT
!
!     History:
!     Version      Date        Comment
!     1           02/12/2009   Original version: Marco Matricardi. ECMWF.
!
! 2010/03 Code cleaning Pascal Brunel, Philippe Marguinaud
! 2011/11 Introduced PC multiple band and cloudy computations, Marco Matricardi, ECMWF
!
!     Code description:
!       Language:              Fortran 90.
!       Software Standards:    "European Standards for Writing and Documenting
!                              Exchangeable Fortran 90 code".
!
!     Module used:
  USE rttov_types, ONLY :  &
       & rttov_options,     &
       & rttov_chanprof,    &
       & rttov_pccomp,      &
       & rttov_coef_pccomp, &
       & radiance_Type
!INTF_OFF
  USE yomhook, ONLY : LHOOK, DR_HOOK
  USE parkind1, ONLY : jpim, jprb
!INTF_ON
  IMPLICIT NONE
!subroutine arguments:
  TYPE(rttov_options    ), INTENT(IN)    :: opts
  TYPE(rttov_chanprof   ), INTENT(IN)    :: chanprof   (:) ! Channel indices
  TYPE(rttov_chanprof   ), INTENT(IN)    :: chanprof_pc(:) ! Channel indices
  TYPE(rttov_pccomp     ), INTENT(INOUT) :: pccomp_ad
  TYPE(rttov_pccomp     ), INTENT(INOUT) :: pccomp
  TYPE(rttov_coef_pccomp), INTENT(IN)    :: coef_pccomp
  TYPE(radiance_Type    ), INTENT(INOUT) :: radiancedata_ad
!INTF_END
  INTEGER(KIND=jpim) :: i, k, chan, prof
  INTEGER(KIND=jpim) :: nchannels, nchannels_p    ! Number of radiances computed
  INTEGER(KIND=jpim) :: npcscores, npcscores_p
  INTEGER(KIND=jpim) :: nprofiles

  REAL   (KIND=JPRB) :: noise_r(size(coef_pccomp%noise_r))

  REAL   (KIND=JPRB) :: ZHOOK_HANDLE
!-----End of header-------------------------------------------------------------
!-------------------------------------------------------------------------------
!         1.   CALCULATE PRINCIPAL COMPONENT SCORES
!-------------------------------------------------------------------------------
  IF (LHOOK) CALL DR_HOOK('RTTOV_PCSCORES_AD', 0_jpim, ZHOOK_HANDLE)
  nchannels = size(chanprof)
  npcscores = size(chanprof_pc) 
  nprofiles = chanprof(nchannels)%prof

  nchannels_p = nchannels / nprofiles 
  npcscores_p = npcscores / nprofiles

! DAR: pre-calculate required channel noise
  DO i = 1, nchannels_p
    chan = chanprof(i)%chan
    noise_r(i) = coef_pccomp%noise_r(chan)
  ENDDO

! DAR: RTTOV 11.2 - Improve pcscores performance by reducing memory access and
! by improving use of vectorisation
  DO prof = 1, nprofiles
    k = 1
    DO i = (prof-1) * nchannels_p + 1, prof * nchannels_p
      radiancedata_ad%total(i) = radiancedata_ad%total(i) + &
         DOT_PRODUCT(&
         coef_pccomp%pcreg(opts%rt_ir%pc%ipcbnd, opts%rt_ir%pc%ipcreg)% &
         coefficients_t(1:npcscores_p, k), &
         pccomp_ad%pcscores((prof-1)*npcscores_p + 1 : prof*npcscores_p)) * &
         noise_r(k)
      k = k + 1
    ENDDO
  ENDDO

  pccomp_ad%pcscores = 0._jprb

  IF (LHOOK) CALL DR_HOOK('RTTOV_PCSCORES_AD', 1_jpim, ZHOOK_HANDLE)
END SUBROUTINE rttov_pcscores_ad
