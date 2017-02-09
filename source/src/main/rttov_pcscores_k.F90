!    Compute principal component scores
SUBROUTINE rttov_pcscores_k( &
            & opts,        &
            & chanprof,    &
            & chanprof_pc, &
            & pccomp,      &
            & pccomp_k,    &
            & coef_pccomp, &
            & total_k_pc)
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
!     1           02/12/2009   Oroginal version: Marco Matricardi. ECMWF.
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
       & rttov_coef_pccomp
  USE parkind1, ONLY : jprb
!INTF_OFF
  USE yomhook, ONLY : LHOOK, DR_HOOK
  USE parkind1, ONLY : jpim
!INTF_ON
  IMPLICIT NONE
!subroutine arguments:
  TYPE(rttov_options    ), INTENT(IN)    :: opts
  TYPE(rttov_chanprof   ), INTENT(IN)    :: chanprof   (:)     ! Channel indices
  TYPE(rttov_chanprof   ), INTENT(IN)    :: chanprof_pc(:)     ! Channel indices
  TYPE(rttov_pccomp     ), INTENT(INOUT) :: pccomp_k
  TYPE(rttov_pccomp     ), INTENT(INOUT) :: pccomp
  TYPE(rttov_coef_pccomp), INTENT(IN)    :: coef_pccomp
  REAL(KIND=jprb)        , INTENT(INOUT) :: total_k_pc(:, :, :)
!INTF_END
  INTEGER(KIND=jpim) :: j, chan, chan_pc, prof
  INTEGER(KIND=jpim) :: nprofiles                ! Number of profiles
  INTEGER(KIND=jpim) :: nchannels                ! Number of radiances computed (channels used * profiles)
  INTEGER(KIND=jpim) :: npcscores_p, nchannels_p, i

  REAL   (KIND=JPRB) :: noise_r(SIZE(coef_pccomp%noise_r))
  REAL   (KIND=JPRB) :: ZHOOK_HANDLE
!-----End of header-------------------------------------------------------------
!-------------------------------------------------------------------------------
!         1.   CALCULATE PRINCIPAL COMPONENT SCORES
!-------------------------------------------------------------------------------
  IF (LHOOK) CALL DR_HOOK('RTTOV_PCSCORES_K', 0_jpim, ZHOOK_HANDLE)

  npcscores_p = size(total_k_pc, 2)
  nprofiles   = size(total_k_pc, 3)
  nchannels   = size(chanprof)

  nchannels_p = nchannels / nprofiles 
  
  DO i = 1, nchannels_p
    chan = chanprof(i)%chan
    noise_r(i) = coef_pccomp%noise_r(chan)
  ENDDO

! DAR: Update for RTTOV 11.2
! DAR: RTTOV 11.2 - Improve pcscores performance by using transposed total_k_pc
! array facilitating improved use of vectorisation
  DO prof = 1, nprofiles
    DO j = 1, npcscores_p
      chan_pc = chanprof_pc(j)%chan

      total_k_pc(1:nchannels_p, j, prof) = &
        pccomp_k%pcscores((prof - 1) * npcscores_p + j) * &
        coef_pccomp%pcreg(opts%rt_ir%pc%ipcbnd, opts%rt_ir%pc%ipcreg)% &
        coefficients(1:nchannels_p, chan_pc) * &
        noise_r(1:nchannels_p)
    ENDDO
  ENDDO

  IF (LHOOK) CALL DR_HOOK('RTTOV_PCSCORES_K', 1_jpim, ZHOOK_HANDLE)
END SUBROUTINE rttov_pcscores_k
