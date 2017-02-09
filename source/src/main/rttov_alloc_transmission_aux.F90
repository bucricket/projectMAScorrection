!
SUBROUTINE rttov_alloc_transmission_aux( &
            & err,              &
            & opts,             &
            & transmission_aux, &
            & nlayers,          &
            & nchannels,        &
            & asw,              &
            & nstreams,         &
            & init,             &
            & direct)
! Description:
! allocation/deallocation of a profile structure
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
!  1.2       03/11/2009  Transmittances / optical depths on levels (A Geer)
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
  USE rttov_types, ONLY : rttov_options, transmission_type_aux
  USE parkind1, ONLY : jpim, jplm
  IMPLICIT NONE
! subroutine arguments
! scalar arguments with intent(out):
  INTEGER(KIND=jpim)         , INTENT(OUT)          :: err             ! return code
  TYPE(rttov_options)        , INTENT(IN)           :: opts
  TYPE(transmission_type_aux), INTENT(INOUT)        :: transmission_aux
  INTEGER(KIND=jpim)         , INTENT(IN)           :: nlayers
  INTEGER(KIND=jpim)         , INTENT(IN)           :: nchannels
  INTEGER(KIND=jpim)         , INTENT(IN)           :: asw             ! 1=allocate, 0=deallocate
  INTEGER(KIND=jpim)         , INTENT(IN)           :: nstreams
  LOGICAL(KIND=jplm)         , INTENT(IN), OPTIONAL :: init
  LOGICAL(KIND=jplm)         , INTENT(IN), OPTIONAL :: direct
!INTF_END
#include "rttov_errorreport.interface"
#include "rttov_init_transmission_aux.interface"
! Local Arrays and Scalars:
  INTEGER(KIND=jpim) :: nlevels, su
  LOGICAL(KIND=jplm) :: init1, direct1
!- End of header --------------------------------------------------------
  TRY
  nlevels = nlayers + 1
  init1   = .FALSE.
  IF (PRESENT(init)) init1 = init
  direct1   = .FALSE.
  IF (PRESENT(direct)) direct1 = direct

  IF (opts%rt_ir%addclouds) THEN
    su = 1
  ELSE
    su = 0
  ENDIF

  IF (asw .EQ. 1) THEN
    CALL nullify_struct()

    IF (direct1) THEN
      ALLOCATE (transmission_aux%fac1(nlevels, 0:nstreams, nchannels), STAT = ERR)
      THROWM( ERR .NE. 0 , "allocation of transmission_aux % fac1")
      ALLOCATE (transmission_aux%surf_fac(0:nstreams, nchannels), STAT = ERR)
      THROWM( ERR .NE. 0 , "allocation of transmission_aux % surf_fac")
    ENDIF

    ! Thermal path1
    ALLOCATE (transmission_aux%thermal_path1)
    ALLOCATE (transmission_aux%thermal_path1%tau_level(nlevels, 0:nstreams, nchannels), STAT = ERR)
    THROWM( ERR .NE. 0 , "allocation of transmission_aux % thermal_path1 % tau_level")
    ALLOCATE (transmission_aux%thermal_path1%od_singlelayer(0:su, nlayers, nchannels), STAT = ERR)
    THROWM( ERR .NE. 0 , "allocation of transmission_aux % thermal_path1 % od_singlelayer")
    ALLOCATE (transmission_aux%thermal_path1%od_sfrac(0:nstreams, nchannels), STAT = ERR)
    THROWM( ERR .NE. 0 , "allocation of transmission_aux % thermal_path1 % od_sfrac")
    ALLOCATE (transmission_aux%thermal_path1%tau_surf(0:nstreams, nchannels), STAT = ERR)
    THROWM( ERR .NE. 0 , "allocation of transmission_aux % thermal_path1 % tau_surf")
    ALLOCATE (transmission_aux%thermal_path1%tau_surf_ac(0:nstreams, nchannels), STAT = ERR)
    THROWM( ERR .NE. 0 , "allocation of transmission_aux % thermal_path1 % tau_surf_ac")

    IF (direct1) THEN
      ALLOCATE (transmission_aux%thermal_path1%fac2(nlevels, 0:nstreams, nchannels), STAT = ERR)
      THROWM( ERR .NE. 0 , "allocation of transmission_aux % thermal_path1 % fac2")
      ALLOCATE (transmission_aux%thermal_path1%tau_level_r(nlevels, 0:nstreams, nchannels), STAT = ERR)
      THROWM( ERR .NE. 0 , "allocation of transmission_aux % thermal_path1 % tau_level_r")
      ALLOCATE (transmission_aux%thermal_path1%od_singlelayer_r(0:su, nlayers, nchannels), STAT = ERR)
      THROWM( ERR .NE. 0 , "allocation of transmission_aux % thermal_path1 % od_singlelayer_r")
      ALLOCATE (transmission_aux%thermal_path1%tau_surf_r(0:nstreams, nchannels), STAT = ERR)
      THROWM( ERR .NE. 0 , "allocation of transmission_aux % thermal_path1 % tau_surf_r")
      ALLOCATE (transmission_aux%thermal_path1%od_sfrac_r(0:nstreams, nchannels), STAT = ERR)
      THROWM( ERR .NE. 0 , "allocation of transmission_aux % thermal_path1 % od_sfrac_r")
    ELSE
      NULLIFY(transmission_aux%thermal_path1%fac2)
      NULLIFY(transmission_aux%thermal_path1%tau_level_r)
      NULLIFY(transmission_aux%thermal_path1%od_singlelayer_r)
      NULLIFY(transmission_aux%thermal_path1%tau_surf_r)
      NULLIFY(transmission_aux%thermal_path1%od_sfrac_r)
    ENDIF

    IF(opts%rt_all%do_lambertian .OR. opts%rt_mw%do_lambertian .OR. opts%rt_ir%do_lambertian) THEN
      ALLOCATE (transmission_aux%thermal_path1%tau_level_p(nlevels, 0:nstreams, nchannels), STAT = ERR)
      THROWM( ERR .NE. 0 , "allocation of transmission_aux % thermal_path1 % tau_level_p")
      ALLOCATE (transmission_aux%thermal_path1%tau_level_p_r(nlevels, 0:nstreams, nchannels), STAT = ERR)
      THROWM( ERR .NE. 0 , "allocation of transmission_aux % thermal_path1 % tau_level_p_r")
      ALLOCATE (transmission_aux%thermal_path1%tau_surf_p(0:nstreams, nchannels), STAT = ERR)
      THROWM( ERR .NE. 0 , "allocation of transmission_aux % thermal_path1 % tau_surf_p")
      ALLOCATE (transmission_aux%thermal_path1%tau_surf_p_r(0:nstreams, nchannels), STAT = ERR)
      THROWM( ERR .NE. 0 , "allocation of transmission_aux % thermal_path1 % tau_surf_p_r")
    ELSE
      NULLIFY(transmission_aux%thermal_path1%tau_level_p)
      NULLIFY(transmission_aux%thermal_path1%tau_level_p_r)
      NULLIFY(transmission_aux%thermal_path1%tau_surf_p)
      NULLIFY(transmission_aux%thermal_path1%tau_surf_p_r)
    ENDIF

    IF (opts%rt_ir%addsolar) THEN
      ! Solar path2
      ALLOCATE (transmission_aux%solar_path2)
      ALLOCATE (transmission_aux%solar_path2%tau_level(nlevels, 0:nstreams, nchannels), STAT = ERR)
      THROWM( ERR .NE. 0 , "allocation of transmission_aux % solar_path2 % tau_level")
      ALLOCATE (transmission_aux%solar_path2%od_singlelayer(0:su, nlayers, nchannels), STAT = ERR)
      THROWM( ERR .NE. 0 , "allocation of transmission_aux % solar_path2 % od_singlelayer")
      ALLOCATE (transmission_aux%solar_path2%od_sfrac(0:nstreams, nchannels), STAT = ERR)
      THROWM( ERR .NE. 0 , "allocation of transmission_aux % solar_path2 % od_sfrac")
      ALLOCATE (transmission_aux%solar_path2%tau_surf(0:nstreams, nchannels), STAT = ERR)
      THROWM( ERR .NE. 0 , "allocation of transmission_aux % solar_path2 % tau_surf")
      ALLOCATE (transmission_aux%solar_path2%od_frac_ac(0:nstreams, nchannels), STAT = ERR)
      THROWM( ERR .NE. 0 , "allocation of transmission_aux % solar_path2 % od_frac_ac")
      ALLOCATE (transmission_aux%solar_path2%tau_surf_ac(0:nstreams, nchannels), STAT = ERR)
      THROWM( ERR .NE. 0 , "allocation of transmission_aux % solar_path2 % tau_surf_ac")

      ! Solar path1
      ALLOCATE (transmission_aux%solar_path1)
      ALLOCATE (transmission_aux%solar_path1%tau_level(nlevels, 0:nstreams, nchannels), STAT = ERR)
      THROWM( ERR .NE. 0 , "allocation of transmission_aux % solar_path1 % tau_level")
      ALLOCATE (transmission_aux%solar_path1%tau_surf(0:nstreams, nchannels), STAT = ERR)
      THROWM( ERR .NE. 0 , "allocation of transmission_aux % solar_path1 % tau_surf")

      NULLIFY(transmission_aux%solar_path1%od_singlelayer)
      NULLIFY(transmission_aux%solar_path1%od_sfrac)
      NULLIFY(transmission_aux%solar_path1%tau_surf_ac)
      NULLIFY(transmission_aux%solar_path1%fac2)
      NULLIFY(transmission_aux%solar_path1%tau_level_r)
      NULLIFY(transmission_aux%solar_path1%tau_surf_r)

      IF (opts%rt_ir%addclouds .OR. opts%rt_ir%addaerosl) THEN
        ALLOCATE (transmission_aux%solar_path1%od_singlelayer(0:su, nlayers, nchannels), STAT = ERR)
        THROWM( ERR .NE. 0 , "allocation of transmission_aux % solar_path1 % od_singlelayer")
        ALLOCATE (transmission_aux%solar_path1%od_sfrac(0:nstreams, nchannels), STAT = ERR)
        THROWM( ERR .NE. 0 , "allocation of transmission_aux % solar_path1 % od_sfrac")
        ALLOCATE (transmission_aux%solar_path1%tau_surf_ac(0:nstreams, nchannels), STAT = ERR)
        THROWM( ERR .NE. 0 , "allocation of transmission_aux % solar_path1 % tau_surf_ac")

        IF (direct1) THEN
          ALLOCATE (transmission_aux%solar_path1%fac2(nlevels, 0:nstreams, nchannels), STAT = ERR)
          THROWM( ERR .NE. 0 , "allocation of transmission_aux % solar_path1 % fac2")
          ALLOCATE (transmission_aux%solar_path1%tau_level_r(nlevels, 0:nstreams, nchannels), STAT = ERR)
          THROWM( ERR .NE. 0 , "allocation of transmission_aux % solar_path1 % tau_level_r")
          ALLOCATE (transmission_aux%solar_path1%tau_surf_r(0:nstreams, nchannels), STAT = ERR)
          THROWM( ERR .NE. 0 , "allocation of transmission_aux % solar_path1 % tau_surf_r")
        ENDIF

      ENDIF
    ENDIF

    IF (init1) CALL rttov_init_transmission_aux(opts, transmission_aux)
  ENDIF

  IF (asw .EQ. 0) THEN
    IF (ASSOCIATED(transmission_aux%fac1)) DEALLOCATE (transmission_aux%fac1, STAT = ERR)
    THROWM( ERR .NE. 0 , "deallocation of transmission_aux%fac1")
    IF (ASSOCIATED(transmission_aux%surf_fac)) DEALLOCATE (transmission_aux%surf_fac, STAT = ERR)
    THROWM( ERR .NE. 0 , "deallocation of transmission_aux%surf_fac")

    IF(ASSOCIATED(transmission_aux%thermal_path1)) THEN
      IF (ASSOCIATED(transmission_aux%thermal_path1%tau_surf_p)) &
          DEALLOCATE (transmission_aux%thermal_path1%tau_surf_p, STAT = ERR)
      THROWM( ERR .NE. 0 , "allocation of transmission_aux % thermal_path1 % tau_surf_p")
      IF (ASSOCIATED(transmission_aux%thermal_path1%tau_surf_p_r)) &
          DEALLOCATE (transmission_aux%thermal_path1%tau_surf_p_r, STAT = ERR)
      THROWM( ERR .NE. 0 , "allocation of transmission_aux % thermal_path1 % tau_surf_p_r")
      IF (ASSOCIATED(transmission_aux%thermal_path1%tau_level_p)) &
          DEALLOCATE (transmission_aux%thermal_path1%tau_level_p, STAT = ERR)
      THROWM( ERR .NE. 0 , "allocation of transmission_aux % thermal_path1 % tau_level_p")
      IF (ASSOCIATED(transmission_aux%thermal_path1%tau_level_p_r)) &
          DEALLOCATE (transmission_aux%thermal_path1%tau_level_p_r, STAT = ERR)
      THROWM( ERR .NE. 0 , "allocation of transmission_aux % thermal_path1 % tau_level_p_r")

      IF (ASSOCIATED(transmission_aux%thermal_path1%tau_level)) &
          DEALLOCATE (transmission_aux%thermal_path1%tau_level, STAT = ERR)
      THROWM( ERR .NE. 0 , "deallocation of transmission_aux%thermal_path1%tau_level")
      IF (ASSOCIATED(transmission_aux%thermal_path1%od_singlelayer)) &
          DEALLOCATE (transmission_aux%thermal_path1%od_singlelayer, STAT = ERR)
      THROWM( ERR .NE. 0 , "deallocation of transmission_aux%thermal_path1%od_singlelayer")
      IF (ASSOCIATED(transmission_aux%thermal_path1%od_sfrac)) &
          DEALLOCATE (transmission_aux%thermal_path1%od_sfrac, STAT = ERR)
      THROWM( ERR .NE. 0 , "deallocation of transmission_aux%thermal_path1%od_sfrac")
      IF (ASSOCIATED(transmission_aux%thermal_path1%tau_surf)) &
          DEALLOCATE (transmission_aux%thermal_path1%tau_surf, STAT = ERR)
      THROWM( ERR .NE. 0 , "deallocation of transmission_aux%thermal_path1%tau_surf")
      IF (ASSOCIATED(transmission_aux%thermal_path1%tau_surf_ac)) &
          DEALLOCATE (transmission_aux%thermal_path1%tau_surf_ac, STAT = ERR)
      THROWM( ERR .NE. 0 , "deallocation of transmission_aux%thermal_path1%tau_surf_ac")

      IF (ASSOCIATED(transmission_aux%thermal_path1%fac2)) &
          DEALLOCATE (transmission_aux%thermal_path1%fac2, STAT = ERR)
      THROWM( ERR .NE. 0 , "deallocation of transmission_aux%thermal_path1%fac2")
      IF (ASSOCIATED(transmission_aux%thermal_path1%tau_level_r)) &
          DEALLOCATE (transmission_aux%thermal_path1%tau_level_r, STAT = ERR)
      THROWM( ERR .NE. 0 , "deallocation of transmission_aux%thermal_path1%tau_level_r")
      IF (ASSOCIATED(transmission_aux%thermal_path1%tau_surf_r)) &
          DEALLOCATE (transmission_aux%thermal_path1%tau_surf_r, STAT = ERR)
      THROWM( ERR .NE. 0 , "deallocation of transmission_aux%thermal_path1%tau_surf_r")
      IF (ASSOCIATED(transmission_aux%thermal_path1%od_singlelayer_r)) &
          DEALLOCATE (transmission_aux%thermal_path1%od_singlelayer_r, STAT = ERR)
      THROWM( ERR .NE. 0 , "deallocation of transmission_aux%thermal_path1%od_singlelayer_r")
      IF (ASSOCIATED(transmission_aux%thermal_path1%od_sfrac_r)) &
          DEALLOCATE (transmission_aux%thermal_path1%od_sfrac_r, STAT = ERR)
      THROWM( ERR .NE. 0 , "deallocation of transmission_aux%thermal_path1%od_sfrac_r")

      DEALLOCATE (transmission_aux%thermal_path1, STAT = ERR)
      THROWM( ERR .NE. 0 , "deallocation of transmission_aux%thermal_path1")
    ENDIF

    IF (ASSOCIATED(transmission_aux%solar_path2)) THEN
      IF (ASSOCIATED(transmission_aux%solar_path2%tau_level)) &
          DEALLOCATE (transmission_aux%solar_path2%tau_level, STAT = ERR)
      THROWM( ERR .NE. 0 , "deallocation of transmission_aux%solar_path2%tau_level")
      IF (ASSOCIATED(transmission_aux%solar_path2%od_singlelayer)) &
          DEALLOCATE (transmission_aux%solar_path2%od_singlelayer, STAT = ERR)
      THROWM( ERR .NE. 0 , "deallocation of transmission_aux%solar_path2%od_singlelayer")
      IF (ASSOCIATED(transmission_aux%solar_path2%od_sfrac)) &
          DEALLOCATE (transmission_aux%solar_path2%od_sfrac, STAT = ERR)
      THROWM( ERR .NE. 0 , "deallocation of transmission_aux%solar_path2%od_sfrac")
      IF (ASSOCIATED(transmission_aux%solar_path2%tau_surf)) &
          DEALLOCATE (transmission_aux%solar_path2%tau_surf, STAT = ERR)
      THROWM( ERR .NE. 0 , "deallocation of transmission_aux%solar_path2%tau_surf")
      IF (ASSOCIATED(transmission_aux%solar_path2%od_frac_ac)) &
          DEALLOCATE (transmission_aux%solar_path2%od_frac_ac, STAT = ERR)
      THROWM( ERR .NE. 0 , "deallocation of transmission_aux%solar_path2%od_frac_ac ")
      IF (ASSOCIATED(transmission_aux%solar_path2%tau_surf_ac)) &
          DEALLOCATE (transmission_aux%solar_path2%tau_surf_ac, STAT = ERR)
      THROWM( ERR .NE. 0 , "deallocation of transmission_aux%solar_path2%tau_surf_ac")

      DEALLOCATE (transmission_aux%solar_path2, STAT = ERR)
      THROWM( ERR .NE. 0 , "deallocation of transmission_aux%solar_path2")
    ENDIF

    IF (ASSOCIATED(transmission_aux%solar_path1)) THEN
      IF (ASSOCIATED(transmission_aux%solar_path1%od_singlelayer)) &
          DEALLOCATE (transmission_aux%solar_path1%od_singlelayer, STAT = ERR)
      THROWM( ERR .NE. 0 , "deallocation of transmission_aux%solar_path1%od_singlelayer")
      IF (ASSOCIATED(transmission_aux%solar_path1%od_sfrac)) &
          DEALLOCATE (transmission_aux%solar_path1%od_sfrac, STAT = ERR)
      THROWM( ERR .NE. 0 , "deallocation of transmission_aux%solar_path1%od_sfrac")
      IF (ASSOCIATED(transmission_aux%solar_path1%tau_surf_ac)) &
          DEALLOCATE (transmission_aux%solar_path1%tau_surf_ac, STAT = ERR)
      THROWM( ERR .NE. 0 , "deallocation of transmission_aux%solar_path1%tau_surf_ac")

      IF (ASSOCIATED(transmission_aux%solar_path1%fac2)) &
          DEALLOCATE (transmission_aux%solar_path1%fac2, STAT = ERR)
      THROWM( ERR .NE. 0 , "deallocation of transmission_aux%solar_path1%fac2")
      IF (ASSOCIATED(transmission_aux%solar_path1%tau_level_r)) &
          DEALLOCATE (transmission_aux%solar_path1%tau_level_r, STAT = ERR)
      THROWM( ERR .NE. 0 , "deallocation of transmission_aux%solar_path1%tau_level_r")
      IF (ASSOCIATED(transmission_aux%solar_path1%tau_surf_r)) &
          DEALLOCATE (transmission_aux%solar_path1%tau_surf_r, STAT = ERR)
      THROWM( ERR .NE. 0 , "deallocation of transmission_aux%solar_path1%tau_surf_r")

      IF (ASSOCIATED(transmission_aux%solar_path1%tau_level)) &
          DEALLOCATE (transmission_aux%solar_path1%tau_level, STAT = ERR)
      THROWM( ERR .NE. 0 , "deallocation of transmission_aux%solar_path1%tau_level")
      IF (ASSOCIATED(transmission_aux%solar_path1%tau_surf)) &
          DEALLOCATE (transmission_aux%solar_path1%tau_surf, STAT = ERR)
      THROWM( ERR .NE. 0 , "deallocation of transmission_aux%solar_path1%tau_surf")

      DEALLOCATE (transmission_aux%solar_path1, STAT = ERR)
      THROWM( ERR .NE. 0 , "deallocation of transmission_aux%solar_path1")
    ENDIF

    CALL nullify_struct()
  ENDIF
  CATCH

CONTAINS

  SUBROUTINE nullify_struct()
    NULLIFY (transmission_aux%fac1)
    NULLIFY (transmission_aux%surf_fac)
    NULLIFY (transmission_aux%thermal_path1)
    NULLIFY (transmission_aux%solar_path2)
    NULLIFY (transmission_aux%solar_path1)
  END SUBROUTINE nullify_struct

END SUBROUTINE rttov_alloc_transmission_aux
