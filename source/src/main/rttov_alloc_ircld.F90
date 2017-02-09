!
SUBROUTINE rttov_alloc_ircld( &
            & err,       &
            & opts,      &
            & nprofiles, &
            & irclds,    &
            & nlayers,   &
            & asw,       &
            & init,      &
            & direct)
! Description:
!   Allocation/deallocation of a ircld structure
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
!    Copyright 2007, EUMETSAT, All Rights Reserved.
!
!INTF_OFF
#include "throw.h"
!INTF_ON
  USE rttov_types, ONLY : ircld_type, rttov_options
  USE parkind1, ONLY : jpim, jplm
!INTF_OFF
  USE rttov_const, ONLY : ncldtyp
  USE parkind1, ONLY : jprb
!INTF_ON
  IMPLICIT NONE

  INTEGER(KIND=jpim),  INTENT(OUT)          :: err      ! return code
  TYPE(rttov_options), INTENT(IN)           :: opts
  INTEGER(KIND=jpim),  INTENT(IN)           :: nprofiles
  INTEGER(KIND=jpim),  INTENT(IN)           :: nlayers
  TYPE(ircld_type),    INTENT(INOUT)        :: irclds
  INTEGER(KIND=jpim),  INTENT(IN)           :: asw      ! 1=allocate, 0=deallocate
  LOGICAL(KIND=jplm),  INTENT(IN), OPTIONAL :: init
  LOGICAL(KIND=jplm),  INTENT(IN), OPTIONAL :: direct
!INTF_END
#include "rttov_errorreport.interface"
#include "rttov_init_ircld.interface"

  LOGICAL(KIND=jplm) :: init1, direct1
!- End of header --------------------------------------------------------
  TRY
  init1 = .FALSE.
  IF (PRESENT(init)) init1 = init

  direct1 = .FALSE.
  IF (PRESENT(direct)) direct1 = direct

  IF (asw .EQ. 1) THEN
    CALL nullify_struct()

    ! This is used as a lookup: for each channel/layer/stream it tells us
    ! whether the layer is cloudy or non-cloudy
    IF (direct1) THEN
      IF (opts%rt_ir%addclouds) THEN
        ALLOCATE (irclds%icldarr(0:2*nlayers, nlayers, nprofiles), STAT = err)
        THROWM(err .NE. 0, "allocation of irclds % icldarr ")
      ELSE
        ALLOCATE (irclds%icldarr(0:0, nlayers, nprofiles), STAT = err)
        THROWM(err .NE. 0, "allocation of irclds % icldarr ")
      ENDIF
    ENDIF

    IF (opts%rt_ir%addaerosl .AND. .NOT. opts%rt_ir%user_aer_opt_param) THEN
      ALLOCATE (irclds%tave(nlayers, nprofiles),     &
                irclds%wmixave(nlayers, nprofiles),  &
                irclds%xpresave(nlayers, nprofiles), &
                irclds%esw(nlayers, nprofiles),      &
                irclds%esi(nlayers, nprofiles),      &
                irclds%ppv(nlayers, nprofiles), STAT = err)
      THROWM(err .NE. 0, "allocation of irclds")
    ENDIF

    IF (opts%rt_ir%addclouds) THEN
      IF (.NOT. opts%rt_ir%user_cld_opt_param) THEN
        ALLOCATE (irclds%cldtyp(ncldtyp, nlayers, nprofiles), STAT = err)
        THROWM(err .NE. 0, "allocation of irclds % cldtyp ")
      ENDIF
      ALLOCATE (irclds%xstr(2*nlayers, nprofiles), STAT = err)
      THROWM(err .NE. 0, "allocation of irclds % xstr ")
      ALLOCATE (irclds%cldcfr(nlayers, nprofiles), STAT = err)
      THROWM(err .NE. 0, "allocation of irclds % cldcfr ")
      ALLOCATE (irclds%maxcov(nlayers, nprofiles), STAT = err)
      THROWM(err .NE. 0, "allocation of irclds % maxcov ")
      ALLOCATE (irclds%xstrmax(nlayers, nprofiles), STAT = err)
      THROWM(err .NE. 0, "allocation of irclds % xstrmax ")
      ALLOCATE (irclds%xstrmin(nlayers, nprofiles), STAT = err)
      THROWM(err .NE. 0, "allocation of irclds % xstrmin ")
      ALLOCATE (irclds%a(2*nlayers, nprofiles), STAT = err)
      THROWM(err .NE. 0, "allocation of irclds % a ")

      IF (direct1) THEN
        ALLOCATE (irclds%xstrref1(2*nlayers, 2*nlayers, nprofiles), STAT = err)
        THROWM(err .NE. 0, "allocation of irclds % xstrref1 ")
        ALLOCATE (irclds%xstrref2(2*nlayers, 2*nlayers, nprofiles), STAT = err)
        THROWM(err .NE. 0, "allocation of irclds % xstrref2 ")
        ALLOCATE (irclds%xstrref(2*nlayers, nprofiles), STAT = err)
        THROWM(err .NE. 0, "allocation of irclds % xstrref ")
        ALLOCATE (irclds%xstrminref(nlayers, nprofiles), STAT = err)
        THROWM(err .NE. 0, "allocation of irclds % xstrminref ")
        ALLOCATE (irclds%ntotref(nlayers, nprofiles), STAT = err)
        THROWM(err .NE. 0, "allocation of irclds % ntotref ")
        ALLOCATE (irclds%indexstr(2*nlayers, nprofiles), STAT = err)
        THROWM(err .NE. 0, "allocation of irclds % indexstr ")
        ALLOCATE (irclds%icount1ref(2*nlayers, nprofiles), STAT = err)
        THROWM(err .NE. 0, "allocation of irclds % icount1ref ")
        ALLOCATE (irclds%iloopin(2*nlayers, nprofiles), STAT = err)
        THROWM(err .NE. 0, "allocation of irclds % iloopin ")
        ALLOCATE (irclds%flag(2*nlayers, nprofiles), STAT = err)
        THROWM(err .NE. 0, "allocation of irclds % flag ")
        ALLOCATE (irclds%iflag(2*nlayers, nprofiles), STAT = err)
        THROWM(err .NE. 0, "allocation of irclds % iflag")
      ENDIF

      IF (direct1) THEN
        ALLOCATE (irclds%nstreamref(nprofiles), &
                  irclds%iloop(nprofiles),      &
                  irclds%icount(nprofiles),     &
                  irclds%icounstr(nprofiles),   &
                  irclds%icount1(nprofiles), STAT = err)
        THROWM(err .NE. 0, "allocation of irclds")
      ENDIF
    ENDIF

    ALLOCATE (irclds%xstrclr(nprofiles), STAT = err)
    THROWM(err .NE. 0, "allocation of irclds % xstrclr")
    irclds%xstrclr = 1._jprb

    IF (direct1) THEN
      ALLOCATE (irclds%nstream(nprofiles), STAT = err)
      THROWM(err .NE. 0, "allocation of irclds % nstream")
      irclds%nstream = 0_jpim
    ENDIF

    IF (init1) THEN
      IF (direct1 .AND. opts%rt_ir%addclouds) THEN
        irclds%nstreamref = 0_jpim
        irclds%iloop      = 0_jpim
        irclds%icount     = 0_jpim
        irclds%icounstr   = 0_jpim
        irclds%icount1    = 0_jpim
      ENDIF
      CALL rttov_init_ircld(irclds)
    ENDIF
  ENDIF

  IF (asw .EQ. 0) THEN
    IF (ASSOCIATED(irclds%tave)) THEN
      DEALLOCATE (irclds%tave, STAT = err)
      THROWM(err .NE. 0, "deallocation of irclds % tave")
    ENDIF
    IF (ASSOCIATED(irclds%wmixave)) THEN
      DEALLOCATE (irclds%wmixave, STAT = err)
      THROWM(err .NE. 0, "deallocation of irclds % wmixave")
    ENDIF
    IF (ASSOCIATED(irclds%xpresave)) THEN
      DEALLOCATE (irclds%xpresave, STAT = err)
      THROWM(err .NE. 0, "deallocation of irclds % xpresave")
    ENDIF
    IF (ASSOCIATED(irclds%esw)) THEN
      DEALLOCATE (irclds%esw, STAT = err)
      THROWM(err .NE. 0, "deallocation of irclds % esw")
    ENDIF
    IF (ASSOCIATED(irclds%esi)) THEN
      DEALLOCATE (irclds%esi, STAT = err)
      THROWM(err .NE. 0, "deallocation of irclds % esi")
    ENDIF
    IF (ASSOCIATED(irclds%ppv)) THEN
      DEALLOCATE (irclds%ppv, STAT = err)
      THROWM(err .NE. 0, "deallocation of irclds % ppv")
    ENDIF
    IF (ASSOCIATED(irclds%icldarr)) THEN
      DEALLOCATE (irclds%icldarr, STAT = err)
      THROWM(err .NE. 0, "deallocation of irclds % icldarr")
    ENDIF
    IF (ASSOCIATED(irclds%xstrref1)) THEN
      DEALLOCATE (irclds%xstrref1, STAT = err)
      THROWM(err .NE. 0, "deallocation of irclds % xstrref1")
    ENDIF
    IF (ASSOCIATED(irclds%xstrref2)) THEN
      DEALLOCATE (irclds%xstrref2, STAT = err)
      THROWM(err .NE. 0, "deallocation of irclds % xstrref2")
    ENDIF
    IF (ASSOCIATED(irclds%cldtyp)) THEN
      DEALLOCATE (irclds%cldtyp, STAT = err)
      THROWM(err .NE. 0, "deallocation of irclds % cldtyp")
    ENDIF
    IF (ASSOCIATED(irclds%xstr)) THEN
      DEALLOCATE (irclds%xstr, STAT = err)
      THROWM(err .NE. 0, "deallocation of irclds % xstr")
    ENDIF
    IF (ASSOCIATED(irclds%xstrminref)) THEN
      DEALLOCATE (irclds%xstrminref, STAT = err)
      THROWM(err .NE. 0, "deallocation of irclds % xstrminref")
    ENDIF
    IF (ASSOCIATED(irclds%xstrref)) THEN
      DEALLOCATE (irclds%xstrref, STAT = err)
      THROWM(err .NE. 0, "deallocation of irclds % xstrref")
    ENDIF
    IF (ASSOCIATED(irclds%cldcfr)) THEN
      DEALLOCATE (irclds%cldcfr, STAT = err)
      THROWM(err .NE. 0, "deallocation of irclds % cldcfr")
    ENDIF
    IF (ASSOCIATED(irclds%maxcov)) THEN
      DEALLOCATE (irclds%maxcov, STAT = err)
      THROWM(err .NE. 0, "deallocation of irclds % maxcov")
    ENDIF
    IF (ASSOCIATED(irclds%xstrmax)) THEN
      DEALLOCATE (irclds%xstrmax, STAT = err)
      THROWM(err .NE. 0, "deallocation of irclds % xstrmax")
    ENDIF
    IF (ASSOCIATED(irclds%xstrmin)) THEN
      DEALLOCATE (irclds%xstrmin, STAT = err)
      THROWM(err .NE. 0, "deallocation of irclds % xstrmin")
    ENDIF
    IF (ASSOCIATED(irclds%a)) THEN
      DEALLOCATE (irclds%a, STAT = err)
      THROWM(err .NE. 0, "deallocation of irclds % a")
    ENDIF
    IF (ASSOCIATED(irclds%ntotref)) THEN
      DEALLOCATE (irclds%ntotref, STAT = err)
      THROWM(err .NE. 0, "deallocation of irclds % ntotref")
    ENDIF
    IF (ASSOCIATED(irclds%indexstr)) THEN
      DEALLOCATE (irclds%indexstr, STAT = err)
      THROWM(err .NE. 0, "deallocation of irclds % indexstr")
    ENDIF
    IF (ASSOCIATED(irclds%icount1ref)) THEN
      DEALLOCATE (irclds%icount1ref, STAT = err)
      THROWM(err .NE. 0, "deallocation of irclds % icount1ref")
    ENDIF
    IF (ASSOCIATED(irclds%iloopin)) THEN
      DEALLOCATE (irclds%iloopin, STAT = err)
      THROWM(err .NE. 0, "deallocation of irclds % iloopin")
    ENDIF
    IF (ASSOCIATED(irclds%flag)) THEN
      DEALLOCATE (irclds%flag, STAT = err)
      THROWM(err .NE. 0, "deallocation of irclds % flag")
    ENDIF
    IF (ASSOCIATED(irclds%iflag)) THEN
      DEALLOCATE (irclds%iflag, STAT = err)
      THROWM(err .NE. 0, "deallocation of irclds % iflag")
    ENDIF

    IF (ASSOCIATED(irclds%nstream)) THEN
      DEALLOCATE (irclds%nstream, STAT = err)
      THROWM(err .NE. 0, "deallocation of irclds % nstream")
    ENDIF
    IF (ASSOCIATED(irclds%nstreamref)) THEN
      DEALLOCATE (irclds%nstreamref, STAT = err)
      THROWM(err .NE. 0, "deallocation of irclds % nstreamref")
    ENDIF
    IF (ASSOCIATED(irclds%iloop)) THEN
      DEALLOCATE (irclds%iloop, STAT = err)
      THROWM(err .NE. 0, "deallocation of irclds % iloop")
    ENDIF
    IF (ASSOCIATED(irclds%icount)) THEN
      DEALLOCATE (irclds%icount, STAT = err)
      THROWM(err .NE. 0, "deallocation of irclds % icount")
    ENDIF
    IF (ASSOCIATED(irclds%icounstr)) THEN
      DEALLOCATE (irclds%icounstr, STAT = err)
      THROWM(err .NE. 0, "deallocation of irclds % icounstr")
    ENDIF
    IF (ASSOCIATED(irclds%icount1)) THEN
      DEALLOCATE (irclds%icount1, STAT = err)
      THROWM(err .NE. 0, "deallocation of irclds % icount1")
    ENDIF
    IF (ASSOCIATED(irclds%xstrclr)) THEN
      DEALLOCATE (irclds%xstrclr, STAT = err)
      THROWM(err .NE. 0, "deallocation of irclds % xstrclr")
    ENDIF

    CALL nullify_struct()

  ENDIF
  CATCH

CONTAINS

  SUBROUTINE nullify_struct()
    NULLIFY (irclds%tave)
    NULLIFY (irclds%wmixave)
    NULLIFY (irclds%xpresave)
    NULLIFY (irclds%esw)
    NULLIFY (irclds%esi)
    NULLIFY (irclds%ppv)
    NULLIFY (irclds%icldarr)
    NULLIFY (irclds%xstrref1)
    NULLIFY (irclds%xstrref2)
    NULLIFY (irclds%cldtyp)
    NULLIFY (irclds%xstr)
    NULLIFY (irclds%xstrminref)
    NULLIFY (irclds%xstrref)
    NULLIFY (irclds%cldcfr)
    NULLIFY (irclds%maxcov)
    NULLIFY (irclds%xstrmax)
    NULLIFY (irclds%xstrmin)
    NULLIFY (irclds%a)
    NULLIFY (irclds%ntotref)
    NULLIFY (irclds%indexstr)
    NULLIFY (irclds%icount1ref)
    NULLIFY (irclds%iloopin)
    NULLIFY (irclds%flag)
    NULLIFY (irclds%iflag)
    NULLIFY (irclds%nstream)
    NULLIFY (irclds%nstreamref)
    NULLIFY (irclds%iloop)
    NULLIFY (irclds%icount)
    NULLIFY (irclds%icounstr)
    NULLIFY (irclds%icount1)
    NULLIFY (irclds%xstrclr)
  END SUBROUTINE nullify_struct

END SUBROUTINE rttov_alloc_ircld
