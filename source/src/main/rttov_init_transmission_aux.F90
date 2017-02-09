!
SUBROUTINE rttov_init_transmission_aux(opts, transmission_aux)
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
! Imported Type Definitions:
  USE rttov_types, ONLY : rttov_options, transmission_Type_aux
!INTF_OFF
  USE parkind1, ONLY : jprb
!INTF_ON
  IMPLICIT NONE
! subroutine arguments
! scalar arguments with intent(out):
  TYPE(rttov_options),         INTENT(IN)    :: opts
  TYPE(transmission_type_aux), INTENT(INOUT) :: transmission_aux
!INTF_END
! Local Arrays and Scalars:
!- End of header --------------------------------------------------------
  
  ! Don't need to initialise these:
!   IF (ASSOCIATED (transmission_aux%fac1))            transmission_aux%fac1            = 0._jprb

  ! Thermal path1
  IF (ASSOCIATED (transmission_aux%thermal_path1%tau_surf))       transmission_aux%thermal_path1%tau_surf       = 0._jprb
  IF (ASSOCIATED (transmission_aux%thermal_path1%tau_level))      transmission_aux%thermal_path1%tau_level      = 0._jprb
  IF (ASSOCIATED (transmission_aux%thermal_path1%od_singlelayer)) transmission_aux%thermal_path1%od_singlelayer = 0._jprb
  IF (ASSOCIATED (transmission_aux%thermal_path1%od_sfrac))       transmission_aux%thermal_path1%od_sfrac       = 0._jprb
  IF (ASSOCIATED (transmission_aux%thermal_path1%tau_surf_ac))    transmission_aux%thermal_path1%tau_surf_ac    = 0._jprb

  IF (opts%rt_ir%addsolar) THEN
    ! Solar path2
    IF (ASSOCIATED (transmission_aux%solar_path2%tau_surf))       transmission_aux%solar_path2%tau_surf       = 0._jprb
    IF (ASSOCIATED (transmission_aux%solar_path2%tau_level))      transmission_aux%solar_path2%tau_level      = 0._jprb
    IF (ASSOCIATED (transmission_aux%solar_path2%od_singlelayer)) transmission_aux%solar_path2%od_singlelayer = 0._jprb
    IF (ASSOCIATED (transmission_aux%solar_path2%od_sfrac))       transmission_aux%solar_path2%od_sfrac       = 0._jprb
    IF (ASSOCIATED (transmission_aux%solar_path2%od_frac_ac))     transmission_aux%solar_path2%od_frac_ac     = 0._jprb
    IF (ASSOCIATED (transmission_aux%solar_path2%tau_surf_ac))    transmission_aux%solar_path2%tau_surf_ac    = 0._jprb

    ! Solar path1
    IF (ASSOCIATED (transmission_aux%solar_path1%tau_level))      transmission_aux%solar_path1%tau_level      = 0._jprb
    IF (ASSOCIATED (transmission_aux%solar_path1%tau_surf))       transmission_aux%solar_path1%tau_surf       = 0._jprb
    IF (opts%rt_ir%addclouds .OR. opts%rt_ir%addaerosl) THEN
      IF (ASSOCIATED (transmission_aux%solar_path1%tau_surf_r))     transmission_aux%solar_path1%tau_surf_r     = 0._jprb
      IF (ASSOCIATED (transmission_aux%solar_path1%tau_level_r))    transmission_aux%solar_path1%tau_level_r    = 0._jprb
      IF (ASSOCIATED (transmission_aux%solar_path1%od_singlelayer)) transmission_aux%solar_path1%od_singlelayer = 0._jprb
      IF (ASSOCIATED (transmission_aux%solar_path1%od_sfrac))       transmission_aux%solar_path1%od_sfrac       = 0._jprb
      IF (ASSOCIATED (transmission_aux%solar_path1%tau_surf_ac))    transmission_aux%solar_path1%tau_surf_ac    = 0._jprb
    ENDIF
  ENDIF

END SUBROUTINE 
