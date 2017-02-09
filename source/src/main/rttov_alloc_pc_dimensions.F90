SUBROUTINE rttov_alloc_pc_dimensions( &
            & err,          &
            & opts,         &
            & npcscores,    &
            & nprofiles,    &
            & chanprof_in,  &
            & chanprof_pc,  &
            & asw,          &
            & channels_rec)
! Description:
!   Allocates/deallocates channel structures for PC-RTTOV
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
!INTF_OFF
#include "throw.h"
!INTF_ON
  USE parkind1, ONLY : jpim
  USE rttov_types, ONLY : rttov_chanprof, rttov_options
  IMPLICIT NONE
  INTEGER(KIND=jpim)  , INTENT(OUT)           :: err 
  TYPE(rttov_options ), INTENT(IN)            :: opts
  INTEGER(KIND=jpim)  , INTENT(IN)            :: npcscores
  INTEGER(KIND=jpim)  , INTENT(IN)            :: nprofiles
  TYPE(rttov_chanprof), POINTER               :: chanprof_pc(:)
  TYPE(rttov_chanprof), POINTER               :: chanprof_in(:)
  INTEGER(KIND=jpim)  , INTENT(IN)            :: asw
  INTEGER(KIND=jpim)  , INTENT(IN) , OPTIONAL :: channels_rec(:)
!INTF_END

#include "rttov_errorreport.interface"

  INTEGER(KIND=jpim) :: ichan        , iprof, ichan1
  INTEGER(KIND=jpim) :: nchannels_rec

  TRY

  IF (asw .EQ. 1) THEN
    NULLIFY (chanprof_pc, chanprof_in)
    ALLOCATE (chanprof_pc(npcscores), STAT = err)
    THROWM(err.ne.0,"Allocation of chanprof_pc failed")
    ichan1 = 1
    DO iprof = 1, nprofiles
      DO ichan = 1, npcscores/nprofiles
        chanprof_pc(ichan1)%chan = ichan
        chanprof_pc(ichan1)%prof = iprof
        ichan1                   = ichan1 + 1_jpim
      ENDDO
    ENDDO
    IF ( opts%rt_ir%pc%addradrec ) THEN
      IF (Present(channels_rec)) THEN
        nchannels_rec = size(channels_rec)
        ALLOCATE (chanprof_in(nprofiles * nchannels_rec), STAT = err)
        THROWM(err.ne.0,"Allocation of chanprof_in failed")
        ichan1 = 1
        DO iprof = 1, nprofiles
          DO ichan = 1, nchannels_rec
            chanprof_in(ichan1)%chan = channels_rec(ichan)
            chanprof_in(ichan1)%prof = iprof
            ichan1                   = ichan1 + 1_jpim
          ENDDO
        ENDDO
      ELSE
        err = errorstatus_fatal
        THROWM(err.ne.0,"channels_rec missing with opts%rt_ir%pc%addradrec true")
      ENDIF
    ENDIF
  ENDIF
  IF (asw .EQ. 0) THEN
    IF (Associated(chanprof_in)) THEN
      DEALLOCATE (chanprof_in, STAT = err)
      THROWM(err.ne.0,"DeAllocation of chanprof_in failed")
    ENDIF
    DEALLOCATE (chanprof_pc, STAT = err)
    THROWM(err.ne.0,"DeAllocation of chanprof_pc failed")
    NULLIFY (chanprof_in, chanprof_pc)
  ENDIF

  CATCH

END SUBROUTINE

