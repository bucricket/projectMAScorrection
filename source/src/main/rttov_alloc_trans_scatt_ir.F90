!
SUBROUTINE rttov_alloc_trans_scatt_ir( &
            & err,                   &
            & opts,                  &
            & transmission_scatt_ir, &
            & nchannels,             &
            & ncldtyp,               &
            & nlayers,               &
            & asw,                   &
            & nstreams,              &
            & stream,                &
            & init)
! Description:
! allocation/deallocation of a transmission_scatt_ir structure
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
! Method:
!
! Current Code Owner: SAF NWP
!
! History:
! Version   Date     Comment
! -------   ----     -------
!  1.0       22/20/2002  Creation
!  1.1       03/11/2009  Transmittances / optical depths on levels (A Geer)
!  1.2       02/12/2009  Introduced new variables for the mixed cloud scheme (Marco Matricardi)
!
! 2010/03 Code cleaning Pascal Brunel, Philippe Marguinaud
!
!
! Code Description:
!   Language:           Fortran 90.
!   Software Standards: "European Standards for Writing and
!     Documenting Exchangeable Fortran 90 Code".
!
! Declarations:
! Modules used:
! Imported Parameters:
!INTF_OFF
#include "throw.h"
!INTF_ON
! Imported Type Definitions:
  USE rttov_types, ONLY : transmission_scatt_ir_type, rttov_options
  USE parkind1, ONLY : jpim, jplm
!INTF_OFF
  USE rttov_const, ONLY : naer_max
!INTF_ON
  IMPLICIT NONE
! subroutine arguments
! scalar arguments with intent(out):
  INTEGER(KIND=jpim)              , INTENT(OUT)             :: err                  ! return code
  TYPE(rttov_options)             , INTENT(IN)              :: opts
  INTEGER(KIND=jpim)              , INTENT(IN)              :: nlayers
  INTEGER(KIND=jpim)              , INTENT(IN)              :: nchannels
  INTEGER(KIND=jpim)              , INTENT(IN)              :: ncldtyp
  TYPE(transmission_scatt_ir_type), INTENT(INOUT)           :: transmission_scatt_ir
  INTEGER(KIND=jpim)              , INTENT(IN)              :: asw                  ! 1=allocate, 0=deallocate
  INTEGER(KIND=jpim)              , INTENT(IN)   , OPTIONAL :: nstreams
  LOGICAL(KIND=jplm)              , INTENT(IN)   , OPTIONAL :: stream
  LOGICAL(KIND=jplm)              , INTENT(IN)   , OPTIONAL :: init
!INTF_END
#include "rttov_errorreport.interface"
#include "rttov_init_trans_scatt_ir.interface"
! Local Arrays and Scalars:
  INTEGER(KIND=jpim) :: nlevels, su, naer
  LOGICAL(KIND=jplm) :: init1, stream1
!- End of header --------------------------------------------------------
  TRY
  nlevels = nlayers + 1
  init1   = .FALSE.
  stream1 = .FALSE.
  IF (PRESENT(init)) init1   = init
  IF (PRESENT(stream)) stream1 = stream
  IF (stream1 .AND. (.NOT. PRESENT(nstreams))) err = errorstatus_fatal
  THROWM( ERR .NE. 0, "nstreams argument missing" )

  IF (asw .EQ. 1) THEN

    CALL nullify_struct()

    IF (stream1) THEN
      ALLOCATE (transmission_scatt_ir%OPDPAC(nlevels, 0:nstreams, nchannels), STAT = ERR)
      THROWM( ERR .NE. 0, "mem allocation error")

      IF (opts%rt_ir%addsolar) THEN
        ALLOCATE (transmission_scatt_ir%OPDPACSUN(nlevels, 0:nstreams, nchannels), STAT = ERR)
        THROWM( ERR .NE. 0, "mem allocation error")
      ENDIF
    ELSE
      IF (opts%rt_ir%addclouds) THEN
        su = 1
      ELSE
        su = 0
      ENDIF

      ALLOCATE (transmission_scatt_ir%OPDPACL(0:su, nlayers, nchannels), STAT = ERR)
      THROWM( ERR .NE. 0, "mem allocation error")

      ALLOCATE (transmission_scatt_ir%OPDPAAER(nlayers, nchannels), STAT = ERR)
      THROWM( ERR .NE. 0, "mem allocation error")
      ALLOCATE (transmission_scatt_ir%OPDPSAER(nlayers, nchannels), STAT = ERR)
      THROWM( ERR .NE. 0, "mem allocation error")

      IF (opts%rt_ir%addsolar) THEN
        ALLOCATE (transmission_scatt_ir%OPDPACLSUN(0:su, nlayers, nchannels), STAT = ERR)
        THROWM( ERR .NE. 0, "mem allocation error")
        ALLOCATE (transmission_scatt_ir%AZPHACUP(0:su, nlayers, nchannels), STAT = ERR)
        THROWM( ERR .NE. 0, "mem allocation error")
        ALLOCATE (transmission_scatt_ir%AZPHACDO(0:su, nlayers, nchannels), STAT = ERR)
        THROWM( ERR .NE. 0, "mem allocation error")
        ALLOCATE (transmission_scatt_ir%OPDPABS(0:su, nlayers, nchannels), STAT = ERR)
        THROWM( ERR .NE. 0, "mem allocation error")
        ALLOCATE (transmission_scatt_ir%OPDPSCA(0:su, nlayers, nchannels), STAT = ERR)
        THROWM( ERR .NE. 0, "mem allocation error")
        ALLOCATE (transmission_scatt_ir%OPDPEXT(0:su, nlayers, nchannels), STAT = ERR)
        THROWM( ERR .NE. 0, "mem allocation error")
        ALLOCATE (transmission_scatt_ir%SSA(0:su, nlayers, nchannels), STAT = ERR)
        THROWM( ERR .NE. 0, "mem allocation error")

        ALLOCATE (transmission_scatt_ir%AZPHAERUPA(nlayers, nchannels), STAT = ERR)
        THROWM( ERR .NE. 0, "mem allocation error")
        ALLOCATE (transmission_scatt_ir%AZPHAERDOA(nlayers, nchannels), STAT = ERR)
        THROWM( ERR .NE. 0, "mem allocation error")
      ENDIF

      IF (opts%rt_ir%addaerosl) THEN
        IF (opts%rt_ir%user_aer_opt_param) THEN
          naer = 1
        ELSE
          naer = naer_max
        ENDIF
        ALLOCATE (transmission_scatt_ir%GPARAERA(nlayers, nchannels), STAT = ERR)
        THROWM( ERR .NE. 0, "mem allocation error")
        ALLOCATE (transmission_scatt_ir%GPARAER(nlayers, nchannels), STAT = ERR)
        THROWM( ERR .NE. 0, "mem allocation error")
        ALLOCATE (transmission_scatt_ir%OPDPAERLA(nlayers, nchannels), STAT = ERR)
        THROWM( ERR .NE. 0, "mem allocation error")

        IF (opts%rt_ir%addsolar) THEN
          ALLOCATE (transmission_scatt_ir%PHASINTUPREF(naer, nlayers, nchannels), STAT = ERR)
          THROWM( ERR .NE. 0, "mem allocation error")
          ALLOCATE (transmission_scatt_ir%PHASINTDOREF(naer, nlayers, nchannels), STAT = ERR)
          THROWM( ERR .NE. 0, "mem allocation error")
          ALLOCATE (transmission_scatt_ir%AZPHAERUP(nlayers, nchannels), STAT = ERR)
          THROWM( ERR .NE. 0, "mem allocation error")
          ALLOCATE (transmission_scatt_ir%AZPHAERDO(nlayers, nchannels), STAT = ERR)
          THROWM( ERR .NE. 0, "mem allocation error")
        ENDIF
      ENDIF

      IF (opts%rt_ir%addclouds) THEN

        ALLOCATE (transmission_scatt_ir%OPDPA(nlayers, nchannels), STAT = ERR)
        THROWM( ERR .NE. 0, "mem allocation error")
        ALLOCATE (transmission_scatt_ir%OPDPS(nlayers, nchannels), STAT = ERR)
        THROWM( ERR .NE. 0, "mem allocation error")
        ALLOCATE (transmission_scatt_ir%GPAR(nlayers, nchannels), STAT = ERR)
        THROWM( ERR .NE. 0, "mem allocation error")
        ALLOCATE (transmission_scatt_ir%GPARTOT(nlayers, nchannels), STAT = ERR)
        THROWM( ERR .NE. 0, "mem allocation error")
        ALLOCATE (transmission_scatt_ir%OPDPCLDLA(nlayers, nchannels), STAT = ERR)
        THROWM( ERR .NE. 0, "mem allocation error")

        IF (opts%rt_ir%addsolar) THEN
          ALLOCATE (transmission_scatt_ir%AZPHUP(nlayers, nchannels), STAT = ERR)
          THROWM( ERR .NE. 0, "mem allocation error")
          ALLOCATE (transmission_scatt_ir%AZPHDO(nlayers, nchannels), STAT = ERR)
          THROWM( ERR .NE. 0, "mem allocation error")
          ALLOCATE (transmission_scatt_ir%AZPHUPTOT(nlayers, nchannels), STAT = ERR)
          THROWM( ERR .NE. 0, "mem allocation error")
          ALLOCATE (transmission_scatt_ir%AZPHDOTOT(nlayers, nchannels), STAT = ERR)
          THROWM( ERR .NE. 0, "mem allocation error")
          ALLOCATE (transmission_scatt_ir%AZPHUPCLS(ncldtyp, nlayers, nchannels), STAT = ERR)
          THROWM( ERR .NE. 0, "mem allocation error")
          ALLOCATE (transmission_scatt_ir%AZPHDOCLS(ncldtyp, nlayers, nchannels), STAT = ERR)
          THROWM( ERR .NE. 0, "mem allocation error")
        ENDIF

        IF (.NOT. opts%rt_ir%user_cld_opt_param) THEN
          ALLOCATE (transmission_scatt_ir%OPDPACLS(ncldtyp, nlayers, nchannels), STAT = ERR)
          THROWM( ERR .NE. 0, "mem allocation error")
          ALLOCATE (transmission_scatt_ir%OPDPSCLS(ncldtyp, nlayers, nchannels), STAT = ERR)
          THROWM( ERR .NE. 0, "mem allocation error")
          ALLOCATE (transmission_scatt_ir%GPARCLS(ncldtyp, nlayers, nchannels), STAT = ERR)
          THROWM( ERR .NE. 0, "mem allocation error")
        ENDIF
      ENDIF
    ENDIF

    IF (init1) CALL rttov_init_trans_scatt_ir(transmission_scatt_ir)

  ENDIF

  IF (asw .EQ. 0) THEN
    IF (ASSOCIATED(transmission_scatt_ir%OPDPAAER)) DEALLOCATE (transmission_scatt_ir%OPDPAAER, STAT = ERR)
    THROWM( ERR .NE. 0, "mem deallocation error")
    IF (ASSOCIATED(transmission_scatt_ir%OPDPSAER)) DEALLOCATE (transmission_scatt_ir%OPDPSAER, STAT = ERR)
    THROWM( ERR .NE. 0, "mem deallocation error")
    IF (ASSOCIATED(transmission_scatt_ir%GPARAERA)) DEALLOCATE (transmission_scatt_ir%GPARAERA, STAT = ERR)
    THROWM( ERR .NE. 0, "mem deallocation error")
    IF (ASSOCIATED(transmission_scatt_ir%PHASINTUPREF)) DEALLOCATE (transmission_scatt_ir%PHASINTUPREF, STAT = ERR)
    THROWM( ERR .NE. 0, "mem deallocation error")
    IF (ASSOCIATED(transmission_scatt_ir%PHASINTDOREF)) DEALLOCATE (transmission_scatt_ir%PHASINTDOREF, STAT = ERR)
    THROWM( ERR .NE. 0, "mem deallocation error")
    IF (ASSOCIATED(transmission_scatt_ir%AZPHAERUP)) DEALLOCATE (transmission_scatt_ir%AZPHAERUP, STAT = ERR)
    THROWM( ERR .NE. 0, "mem deallocation error")
    IF (ASSOCIATED(transmission_scatt_ir%AZPHAERDO)) DEALLOCATE (transmission_scatt_ir%AZPHAERDO, STAT = ERR)
    THROWM( ERR .NE. 0, "mem deallocation error")
    IF (ASSOCIATED(transmission_scatt_ir%GPARAER)) DEALLOCATE (transmission_scatt_ir%GPARAER, STAT = ERR)
    THROWM( ERR .NE. 0, "mem deallocation error")
    IF (ASSOCIATED(transmission_scatt_ir%AZPHAERUPA)) DEALLOCATE (transmission_scatt_ir%AZPHAERUPA, STAT = ERR)
    THROWM( ERR .NE. 0, "mem deallocation error")
    IF (ASSOCIATED(transmission_scatt_ir%AZPHAERDOA)) DEALLOCATE (transmission_scatt_ir%AZPHAERDOA, STAT = ERR)
    THROWM( ERR .NE. 0, "mem deallocation error")
    IF (ASSOCIATED(transmission_scatt_ir%OPDPAERLA)) DEALLOCATE (transmission_scatt_ir%OPDPAERLA, STAT = ERR)
    THROWM( ERR .NE. 0, "mem deallocation error")
    IF (ASSOCIATED(transmission_scatt_ir%AZPHUP)) DEALLOCATE (transmission_scatt_ir%AZPHUP, STAT = ERR)
    THROWM( ERR .NE. 0, "mem deallocation error")
    IF (ASSOCIATED(transmission_scatt_ir%AZPHDO)) DEALLOCATE (transmission_scatt_ir%AZPHDO, STAT = ERR)
    THROWM( ERR .NE. 0, "mem deallocation error")
    IF (ASSOCIATED(transmission_scatt_ir%OPDPA)) DEALLOCATE (transmission_scatt_ir%OPDPA, STAT = ERR)
    THROWM( ERR .NE. 0, "mem deallocation error")
    IF (ASSOCIATED(transmission_scatt_ir%OPDPS)) DEALLOCATE (transmission_scatt_ir%OPDPS, STAT = ERR)
    THROWM( ERR .NE. 0, "mem deallocation error")
    IF (ASSOCIATED(transmission_scatt_ir%GPAR)) DEALLOCATE (transmission_scatt_ir%GPAR, STAT = ERR)
    THROWM( ERR .NE. 0, "mem deallocation error")
    IF (ASSOCIATED(transmission_scatt_ir%OPDPCLDLA)) DEALLOCATE (transmission_scatt_ir%OPDPCLDLA, STAT = ERR)
    THROWM( ERR .NE. 0, "mem deallocation error")
    IF (ASSOCIATED(transmission_scatt_ir%AZPHACUP)) DEALLOCATE (transmission_scatt_ir%AZPHACUP, STAT = ERR)
    THROWM( ERR .NE. 0, "mem deallocation error")
    IF (ASSOCIATED(transmission_scatt_ir%AZPHACDO)) DEALLOCATE (transmission_scatt_ir%AZPHACDO, STAT = ERR)
    THROWM( ERR .NE. 0, "mem deallocation error")
    IF (ASSOCIATED(transmission_scatt_ir%OPDPABS)) DEALLOCATE (transmission_scatt_ir%OPDPABS, STAT = ERR)
    THROWM( ERR .NE. 0, "mem deallocation error")
    IF (ASSOCIATED(transmission_scatt_ir%OPDPSCA)) DEALLOCATE (transmission_scatt_ir%OPDPSCA, STAT = ERR)
    THROWM( ERR .NE. 0, "mem deallocation error")
    IF (ASSOCIATED(transmission_scatt_ir%OPDPEXT)) DEALLOCATE (transmission_scatt_ir%OPDPEXT, STAT = ERR)
    THROWM( ERR .NE. 0, "mem deallocation error")
    IF (ASSOCIATED(transmission_scatt_ir%OPDPACLSUN)) DEALLOCATE (transmission_scatt_ir%OPDPACLSUN, STAT = ERR)
    THROWM( ERR .NE. 0, "mem deallocation error")
    IF (ASSOCIATED(transmission_scatt_ir%OPDPACSUN)) DEALLOCATE (transmission_scatt_ir%OPDPACSUN, STAT = ERR)
    THROWM( ERR .NE. 0, "mem deallocation error")
    IF (ASSOCIATED(transmission_scatt_ir%OPDPACL)) DEALLOCATE (transmission_scatt_ir%OPDPACL, STAT = ERR)
    THROWM( ERR .NE. 0, "mem deallocation error")
    IF (ASSOCIATED(transmission_scatt_ir%OPDPAC)) DEALLOCATE (transmission_scatt_ir%OPDPAC, STAT = ERR)
    THROWM( ERR .NE. 0, "mem deallocation error")
    IF (ASSOCIATED(transmission_scatt_ir%SSA)) DEALLOCATE (transmission_scatt_ir%SSA, STAT = ERR)
    THROWM( ERR .NE. 0, "mem deallocation error")
    IF (ASSOCIATED(transmission_scatt_ir%OPDPACLS)) DEALLOCATE (transmission_scatt_ir%OPDPACLS, STAT = ERR)
    THROWM( ERR .NE. 0, "mem deallocation error")
    IF (ASSOCIATED(transmission_scatt_ir%OPDPSCLS)) DEALLOCATE (transmission_scatt_ir%OPDPSCLS, STAT = ERR)
    THROWM( ERR .NE. 0, "mem deallocation error")
    IF (ASSOCIATED(transmission_scatt_ir%GPARCLS)) DEALLOCATE (transmission_scatt_ir%GPARCLS, STAT = ERR)
    THROWM( ERR .NE. 0, "mem deallocation error")
    IF (ASSOCIATED(transmission_scatt_ir%GPARTOT)) DEALLOCATE (transmission_scatt_ir%GPARTOT, STAT = ERR)
    THROWM( ERR .NE. 0, "mem deallocation error")
    IF (ASSOCIATED(transmission_scatt_ir%AZPHUPCLS)) DEALLOCATE (transmission_scatt_ir%AZPHUPCLS, STAT = ERR)
    THROWM( ERR .NE. 0, "mem deallocation error")
    IF (ASSOCIATED(transmission_scatt_ir%AZPHDOCLS)) DEALLOCATE (transmission_scatt_ir%AZPHDOCLS, STAT = ERR)
    THROWM( ERR .NE. 0, "mem deallocation error")
    IF (ASSOCIATED(transmission_scatt_ir%AZPHUPTOT)) DEALLOCATE (transmission_scatt_ir%AZPHUPTOT, STAT = ERR)
    THROWM( ERR .NE. 0, "mem deallocation error")
    IF (ASSOCIATED(transmission_scatt_ir%AZPHDOTOT)) DEALLOCATE (transmission_scatt_ir%AZPHDOTOT, STAT = ERR)
    THROWM( ERR .NE. 0, "mem deallocation error")

    CALL nullify_struct()

  ENDIF
  CATCH

CONTAINS

  SUBROUTINE nullify_struct()
    NULLIFY (transmission_scatt_ir%OPDPAAER)
    NULLIFY (transmission_scatt_ir%OPDPSAER)
    NULLIFY (transmission_scatt_ir%GPARAERA)
    NULLIFY (transmission_scatt_ir%PHASINTUPREF)
    NULLIFY (transmission_scatt_ir%PHASINTDOREF)
    NULLIFY (transmission_scatt_ir%AZPHAERUP)
    NULLIFY (transmission_scatt_ir%AZPHAERDO)
    NULLIFY (transmission_scatt_ir%GPARAER)
    NULLIFY (transmission_scatt_ir%AZPHAERUPA)
    NULLIFY (transmission_scatt_ir%AZPHAERDOA)
    NULLIFY (transmission_scatt_ir%OPDPAERLA)
    NULLIFY (transmission_scatt_ir%AZPHUP)
    NULLIFY (transmission_scatt_ir%AZPHUPCLS)
    NULLIFY (transmission_scatt_ir%AZPHDO)
    NULLIFY (transmission_scatt_ir%AZPHDOCLS)
    NULLIFY (transmission_scatt_ir%AZPHUPTOT)
    NULLIFY (transmission_scatt_ir%AZPHDOTOT)
    NULLIFY (transmission_scatt_ir%OPDPA)
    NULLIFY (transmission_scatt_ir%OPDPACLS)
    NULLIFY (transmission_scatt_ir%OPDPS)
    NULLIFY (transmission_scatt_ir%OPDPSCLS)
    NULLIFY (transmission_scatt_ir%GPAR)
    NULLIFY (transmission_scatt_ir%GPARTOT)
    NULLIFY (transmission_scatt_ir%GPARCLS)
    NULLIFY (transmission_scatt_ir%OPDPCLDLA)
    NULLIFY (transmission_scatt_ir%AZPHACUP)
    NULLIFY (transmission_scatt_ir%AZPHACDO)
    NULLIFY (transmission_scatt_ir%OPDPABS)
    NULLIFY (transmission_scatt_ir%OPDPSCA)
    NULLIFY (transmission_scatt_ir%OPDPEXT)
    NULLIFY (transmission_scatt_ir%OPDPACLSUN)
    NULLIFY (transmission_scatt_ir%OPDPACSUN)
    NULLIFY (transmission_scatt_ir%OPDPACL)
    NULLIFY (transmission_scatt_ir%OPDPAC)
    NULLIFY (transmission_scatt_ir%SSA)
  END SUBROUTINE nullify_struct

END SUBROUTINE rttov_alloc_trans_scatt_ir
