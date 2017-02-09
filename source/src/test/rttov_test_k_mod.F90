MODULE rttov_test_k_mod
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
! Current Code Owner: SAF NWP
!
! History:
! Version   Date     Comment
! -------   ----     -------
! 2010/03 Code cleaning Pascal Brunel, Philippe Marguinaud
!
#include "throw.h"

Use rttov_chain, Only : chain, pchain, chain_profile_type, &
                        free_chain, advance_chain, get_pointer_chain, &
                        zero_chain
Use parkind1, Only : jpim, jprb, jplm
Use rttov_types, Only : profile_type

Implicit None

Private
Public :: make_chain_profile, free_chain_profile, assign_chain_profile

#include "rttov_errorreport.interface"

CONTAINS

  SUBROUTINE make_chain_profile (err, chain_profiles, c, t, profiles, lzero)
!
  INTEGER(KIND=jpim), INTENT(OUT) :: err
  TYPE(chain),        POINTER     :: chain_profiles(:)
  TYPE(pchain),       POINTER     :: c(:)
  CHARACTER(LEN=*),   INTENT(IN)  :: t
  TYPE(profile_type), INTENT(IN), TARGET :: profiles(:)
  LOGICAL(KIND=jplm), INTENT(IN)  :: lzero
!
  INTEGER(KIND=jpim) :: iprof, nprof
!
TRY
    nprof = size(profiles)

    ALLOCATE (chain_profiles(nprof), c(nprof), STAT = err)
    THROWM(err.ne.0,"Allocation of chain_profiles failed")

    DO iprof = 1, nprof

      CALL chain_profile_type( &
            & err,                    &
            & chain_profiles(iprof),  &
            & Trim(t),                &
            & a0 = profiles(iprof))
      THROW(err.ne.0)

      c(iprof)%p => chain_profiles(iprof)
      
      IF (lzero) &
      CALL zero_chain(c(iprof)%p)

    ENDDO

CATCH
  END SUBROUTINE

  SUBROUTINE free_chain_profile (err, chain_profiles, c)
!
  INTEGER(KIND=jpim), INTENT(OUT) :: err
  TYPE(chain),        POINTER     :: chain_profiles(:)
  TYPE(pchain),       POINTER     :: c(:)
!
  INTEGER(KIND=jpim) :: iprof, nprof
!
TRY
  IF (associated(chain_profiles)) THEN
    nprof = size(chain_profiles)
    DO iprof = 1, nprof
      CALL free_chain(chain_profiles(iprof))
    ENDDO
    DEALLOCATE (chain_profiles, c, STAT = err)
    THROWM(err.ne.0,"DeAllocation of chain_profiles failed")
  ENDIF
CATCH
  END SUBROUTINE


  SUBROUTINE assign_chain_profile (c, w)
!
  TYPE(pchain), INTENT(INOUT) :: c(:)
  REAL(KIND=jprb), INTENT(IN) :: w(:)
!
  INTEGER(KIND=JPIM) :: nprof, iprof
  REAL(KIND=JPRB), POINTER :: x
!
    nprof = size(w)
    DO iprof = 1, nprof
      CALL get_pointer_chain(c(iprof)%p, x)
      x = w(iprof)
      CALL advance_chain(c(iprof)%p)
    ENDDO

  END SUBROUTINE

END MODULE
