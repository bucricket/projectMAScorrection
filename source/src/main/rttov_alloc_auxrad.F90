SUBROUTINE rttov_alloc_auxrad( &
            & err,       &
            & auxrad,    &
            & nlevels,   &
            & nchannels, &
            & asw)
! Description:
!   Allocates/deallocates the auxiliary radiance structure
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
  USE rttov_types, ONLY : radiance_aux
  INTEGER(KIND=jpim), INTENT(OUT)   :: err
  TYPE(radiance_aux), INTENT(INOUT) :: auxrad
  INTEGER(KIND=jpim), INTENT(IN)    :: nlevels
  INTEGER(KIND=jpim), INTENT(IN)    :: nchannels
  INTEGER(KIND=jpim), INTENT(IN)    :: asw
!INTF_END

#include "rttov_errorreport.interface"

  TRY
  IF (asw .EQ. 1) THEN
    ALLOCATE ( &
      auxrad%air_t_eff(nlevels, nchannels), &
      auxrad%surf_t_eff(nchannels), &
      auxrad%skin_t_eff(nchannels), &
      auxrad%cosmic_t_eff(nchannels), STAT = err)

    ALLOCATE ( &
      auxrad%air(nlevels, nchannels), &
      auxrad%surfair(nchannels), &
      auxrad%skin(nchannels), &
      auxrad%cosmic(nchannels), STAT = err)
    THROWM(err.NE.0,"Allocation of auxrad failed")
  ENDIF
  IF (asw .EQ. 0) THEN
    DEALLOCATE (auxrad%air_t_eff, auxrad%surf_t_eff, auxrad%skin_t_eff, &
      auxrad%cosmic_t_eff, STAT = err)

    DEALLOCATE (auxrad%air, auxrad%surfair, auxrad%skin, auxrad%cosmic, &
      STAT = err)
    THROWM(err.NE.0,"DeAllocation of auxrad failed")
    NULLIFY (auxrad%air, auxrad%surfair, auxrad%skin, auxrad%cosmic)
  ENDIF
  CATCH
END SUBROUTINE 
