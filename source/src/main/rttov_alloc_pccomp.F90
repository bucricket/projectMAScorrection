! Description:
!> @file
!!   Allocate/deallocate pccomp structure for PC-RTTOV.
!
!> @brief
!!   Allocate/deallocate pccomp structure for PC-RTTOV.
!!
!! @details
!!   The pccomp structure contains the output PC scores and
!!   reconstructed radiances, BTs for the PC-RTTOV direct model,
!!   the corresponding perturbations for the TL model, and the
!!   input gradients and perturbations for the AD and K models.
!!
!! @param[out]    err            status on exit
!! @param[in,out] pccomp         PC scores and radiances for PC-RTTOV
!! @param[in]     npcscores      number of PC scores to simulate
!! @param[in]     asw            1_jpim => allocate; 0_jpim => deallocate
!! @param[in]     init           set .TRUE. to initialise newly allocated structures, optional
!! @param[in]     nchannels_rec  total number of reconstructed radiances for all profiles, optional
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
!    Copyright 2015, EUMETSAT, All Rights Reserved.
!
SUBROUTINE rttov_alloc_pccomp( &
            & err,          &
            & pccomp,       &
            & npcscores,    &
            & asw,          &
            & init,         &
            & nchannels_rec)
!INTF_OFF
#include "throw.h"
!INTF_ON

  USE parkind1, ONLY : jpim, jplm
  USE rttov_types, ONLY : rttov_pccomp

  IMPLICIT NONE

  INTEGER(KIND=jpim), INTENT(OUT)             :: err
  TYPE(rttov_pccomp), INTENT(INOUT)           :: pccomp
  INTEGER(KIND=jpim), INTENT(IN)              :: npcscores
  INTEGER(KIND=jpim), INTENT(IN)              :: asw
  LOGICAL(KIND=jplm), OPTIONAL   , INTENT(IN) :: init
  INTEGER(KIND=jpim), OPTIONAL   , INTENT(IN) :: nchannels_rec
!INTF_END

#include "rttov_errorreport.interface"
#include "rttov_init_pccomp.interface"

  LOGICAL(KIND=jplm) :: init1
!- End of header --------------------------------------------------------

  TRY
  init1 = .FALSE.
  IF (PRESENT(init)) init1 = init

  IF (asw .EQ. 1) THEN
    NULLIFY (pccomp%pcscores, pccomp%bt_pccomp, pccomp%total_pccomp)
    ALLOCATE (pccomp%pcscores(npcscores), STAT = err)
    THROWM(err .NE. 0, "allocation of pcscores")
    IF (PRESENT(nchannels_rec)) THEN
      IF (nchannels_rec > 0) THEN
        ALLOCATE (pccomp%bt_pccomp(nchannels_rec), pccomp%total_pccomp(nchannels_rec), STAT = err)
        THROWM(err .NE. 0, "allocation of bt_pccomp, total_pccomp")
      ENDIF
    ENDIF
    IF (init1) CALL rttov_init_pccomp(pccomp)
  ENDIF

  IF (asw .EQ. 0) THEN
    IF (ASSOCIATED(pccomp%pcscores)) THEN
      DEALLOCATE (pccomp%pcscores, STAT = err)
      THROWM(err .NE. 0, "deallocation of pcscores")
    ENDIF
    IF (ASSOCIATED(pccomp%bt_pccomp)) THEN
      DEALLOCATE (pccomp%bt_pccomp, STAT = err)
      THROWM(err .NE. 0, "deallocation of bt_pccomp")
    ENDIF
    IF (ASSOCIATED(pccomp%total_pccomp)) THEN
      DEALLOCATE (pccomp%total_pccomp, STAT = err)
      THROWM(err .NE. 0, "deallocation of total_pccomp")
    ENDIF
    NULLIFY (pccomp%pcscores, pccomp%bt_pccomp, pccomp%total_pccomp)
  ENDIF
  CATCH
END SUBROUTINE 
