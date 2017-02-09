!
SUBROUTINE rttov_alloc_raytracing( &
            & ERR,          &
            & nraytracings, &
            & addsolar,     &
            & raytracings,  &
            & nlevels,      &
            & asw,          &
            & init)
! Description:
! allocation/deallocation of a raytracing structure
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
!  1.0       01/12/2002  New F90 code
!  1.1       11/10/2007  Add addclouds, addaerosl, init logicals
!                        nullify unused pointers P.Marguinaud
!  1.2       02/12/2009  Pathsat, Patsun, Patheff and Ltick are allocated
!                        on number of layers (Marco Matricardi).
! Code Description:
!   Language:           Fortran 90.
!   Software Standards: "European Standards for Writing and
!     Documenting Exchangeable Fortran 90 Code".
!
! 2010/03 Code cleaning Pascal Brunel, Philippe Marguinaud
!
!
! Declarations:
! Modules used:
! Imported Parameters:
!INTF_OFF
#include "throw.h"
!INTF_ON
! Imported Type Definitions:
  USE rttov_types, ONLY : raytracing_Type
  USE parkind1, ONLY : jpim, jplm
  IMPLICIT NONE
! subroutine arguments
! scalar arguments with intent(out):
  INTEGER(KIND=jpim)   , INTENT(OUT)             :: ERR         ! return code
  INTEGER(KIND=jpim)   , INTENT(IN)              :: nraytracings
  LOGICAL(KIND=jplm)   , INTENT(IN)              :: addsolar
  INTEGER(KIND=jpim)   , INTENT(IN)              :: nlevels
  TYPE(raytracing_Type), INTENT(INOUT)           :: raytracings
  INTEGER(KIND=jpim)   , INTENT(IN)              :: asw         ! 1=allocate, 0=deallocate
  LOGICAL(KIND=jplm)   , INTENT(IN)   , OPTIONAL :: init
!INTF_END
#include "rttov_errorreport.interface"
#include "rttov_init_raytracing.interface"
! Local Arrays and Scalars:
  INTEGER(KIND=jpim) :: nlayers
  LOGICAL(KIND=jplm) :: init1
!- End of header --------------------------------------------------------
  TRY
  nlayers = nlevels - 1
  init1   = .FALSE.
  IF (Present(init)) init1 = init
!  Allocate section
  IF (asw .EQ. 1) THEN
    ALLOCATE (                                                                                        &
      & raytracings%HGPL(nlevels, nraytracings), raytracings%DMAIR_R(nlevels, nraytracings),          &
      & raytracings%R(nlevels, nraytracings), raytracings%REFRACTIVITY(nlevels, nraytracings),        &
      & raytracings%R_R(nlevels, nraytracings), raytracings%Z_R(nlayers, nraytracings),               &
      & raytracings%RATOESAT(nlayers, nraytracings), raytracings%ZASAT(nlayers, nraytracings),        &
      & raytracings%INT(nlevels, nraytracings), raytracings%ZTEMP(nlevels, nraytracings),             &
      & raytracings%PPW(nlevels, nraytracings), raytracings%DISPCO2(nlevels, nraytracings),           &
      & raytracings%PATHSAT(nlayers, nraytracings), raytracings%PATHSAT_rsqrt(nlayers, nraytracings), &
      & raytracings%PATHSAT_sqrt(nlayers, nraytracings), raytracings%LTICK(nlayers, nraytracings),    &
      & raytracings%CO2_CM(nraytracings), STAT = err)
    THROWM(err.ne.0,"Allocation of raytracing failed")

    IF (addsolar) THEN
      ALLOCATE (                                       &
        & raytracings%RATOESUN(nlayers, nraytracings), &
        & raytracings%ZASUN(nlayers, nraytracings),    &
        & raytracings%PATHSUN(nlayers, nraytracings),  &
        & raytracings%PATHEFF(nlayers, nraytracings), STAT = err)
      THROWM(err.ne.0,"Allocation of solar raytracing failed")
    ENDIF

    IF (init1) THEN
      CALL rttov_init_raytracing(addsolar, raytracings)
    ENDIF

  ENDIF
  IF (asw .EQ. 0) THEN
    DEALLOCATE (                                        &
      & raytracings%HGPL, raytracings%DMAIR_R,          &
      & raytracings%R, raytracings%REFRACTIVITY,        &
      & raytracings%R_R, raytracings%Z_R,               &
      & raytracings%RATOESAT, raytracings%ZASAT,        &
      & raytracings%INT, raytracings%ZTEMP,             &
      & raytracings%PPW, raytracings%DISPCO2,           &
      & raytracings%PATHSAT, raytracings%PATHSAT_rsqrt, &
      & raytracings%PATHSAT_sqrt, raytracings%LTICK,    &
      & raytracings%CO2_CM,  STAT = err)
    THROWM(err.ne.0,"Deallocation of raytracing failed")

    IF (addsolar) THEN
      DEALLOCATE (              &
        & raytracings%RATOESUN, &
        & raytracings%ZASUN,    &
        & raytracings%PATHSUN,  &
        & raytracings%PATHEFF, STAT = err)
      THROWM(err.ne.0,"Deallocation of solar raytracing failed")
    ENDIF

    NULLIFY (raytracings%HGPL)
    NULLIFY (raytracings%DMAIR_R)
    NULLIFY (raytracings%R)
    NULLIFY (raytracings%REFRACTIVITY)
    NULLIFY (raytracings%R_R)
    NULLIFY (raytracings%Z_R)
    NULLIFY (raytracings%RATOESUN)
    NULLIFY (raytracings%RATOESAT)
    NULLIFY (raytracings%ZASUN)
    NULLIFY (raytracings%ZASAT)
    NULLIFY (raytracings%INT)
    NULLIFY (raytracings%ZTEMP)
    NULLIFY (raytracings%PPW)
    NULLIFY (raytracings%DISPCO2)
    NULLIFY (raytracings%PATHSAT)
    NULLIFY (raytracings%PATHSAT_SQRT)
    NULLIFY (raytracings%PATHSAT_RSQRT)
    NULLIFY (raytracings%PATHSUN)
    NULLIFY (raytracings%PATHEFF)
    NULLIFY (raytracings%LTICK)
    NULLIFY (raytracings%CO2_CM)
  ENDIF
  CATCH
END SUBROUTINE rttov_alloc_raytracing
