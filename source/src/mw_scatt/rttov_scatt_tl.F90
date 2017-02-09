! Description:
!> @file
!!   Runs RTTOV-SCATT tangent linear (TL) model
!
!> @brief
!!   Runs RTTOV-SCATT tangent linear (TL) model
!!
!! @details
!!   Given input profile and cloud/hydrometeor profile structures and
!!   corresponding profile perturbations, computes the resulting radiance
!!   perturbation from the tangent linear (TL) to the direct model
!!   evaluated at the input profile.
!!
!!   The same profiles structure as standard RTTOV is used to input the
!!   pressure, temperature, water vapour and surface variables. The additional
!!   cloud and hydrometeor profiles are input via the cld_profiles structure:
!!   these profiles specify the the constituent amounts on the "full" pressure
!!   levels defined in the profiles structure and they apply to a domain
!!   bounded by "half" pressure levels defined in the cld_profiles%ph(:) array.
!!   The half-pressure levels lie between the full pressure levels with the
!!   top-most half-pressure typically being zero (space) and the last half-
!!   pressure being the surface pressure, giving (nlevels+1) half-pressures
!!   in total.
!!
!!   RTTOV-SCATT allows a limited number of options to be set in the opts_scatt
!!   argument which configure the simulations.
!!
!!   The chanprof and frequencies arguments should be populated by calling the
!!   rttov_scatt_setupindex subroutine. Both should be allocated with size equal
!!   to the total number of channels being simulated.
!!
!!
!!   The methodology and validation is described in the following:
!!
!!   Bauer, P., 2002: Microwave radiative transfer modeling in clouds and precipitation.
!!     Part I: Model description.
!!     NWP SAF Report No. NWPSAF-EC-TR-005, 21 pp.
!!
!!   Moreau, E., P. Bauer and F. Chevallier, 2002: Microwave radiative transfer modeling in clouds and precipitation.
!!     Part II: Model evaluation.
!!     NWP SAF Report No. NWPSAF-EC-TR-006, 27 pp.
!!
!!   Chevallier, F., and P. Bauer, 2003:
!!     Model rain and clouds over oceans:comparison with SSM/I observations. Mon. Wea. Rev., 131, 1240-1255.
!!
!!   Smith, E. A., P. Bauer, F. S. Marzano, C. D. Kummerow, D. McKague, A. Mugnai, G. Panegrossi, 2002:
!!     Intercomparison of microwave radiative transfer models for precipitating clouds.
!!     IEEE Trans. Geosci. Remote Sens. 40, 541-549.
!!
!!   Bauer, P., Moreau, E., Chevallier, F., O'Keeffe, U., 2006:
!!     Multiple-scattering microwave radiative transfer for data assimilation applications
!!     Quart. J. R. Meteorol. Soc. 132, 1259-1281
!!
!!   Geer, A.J., Bauer, P. and O'Dell, C.W., 2009:
!!     A Revised Cloud Overlap Scheme for Fast Microwave Radiative Transfer in Rain and Cloud.
!!     Journal of Applied Meteorology and Climatology. 48, 2257-2270
!!
!! @param[out]    errorstatus      status on exit
!! @param[in]     opts_scatt       RTTOV-SCATT options to configure the simulations
!! @param[in]     nlevels          number of levels in profiles structure
!! @param[in]     chanprof         specifies channels and profiles to simulate
!! @param[in]     frequencies      frequency indices for each channel
!! @param[in]     profiles         input atmospheric profiles and surface variables
!! @param[in]     cld_profiles     input cloud and hydrometeor profiles
!! @param[in]     coef_rttov       RTTOV coefficients structure
!! @param[in]     coef_scatt       RTTOV-SCATT Mietable structure
!! @param[in]     calcemiss        flags for internal RTTOV surface emissivity calculation
!! @param[in,out] emissivity       input/output surface emissivities
!! @param[in]     profiles_tl      input atmospheric profile and surface variable perturbations
!! @param[in]     cld_profiles_tl  input cloud and hydrometeor profile perturbations
!! @param[in,out] emissivity_tl    input/output surface emissivity perturbations
!! @param[in,out] radiance         output radiances and corresponding BTs
!! @param[in,out] radiance_tl      output radiance and BT perturbations
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
Subroutine rttov_scatt_tl( &
     & errorstatus,        &! out
     & opts_scatt,         &! in
     & nlevels,            &! in
     & chanprof,           &! in
     & frequencies,        &! in
     & profiles,           &! in  
     & cld_profiles,       &! in
     & coef_rttov,         &! in
     & coef_scatt,         &! in
     & calcemiss,          &! in
     & emissivity,         &! inout
     & profiles_tl,        &! in
     & cld_profiles_tl,    &! in
     & emissivity_tl,      &! inout
     & radiance,           &! inout
     & radiance_tl)         ! inout

  ! History:
  ! Version   Date     Comment
  ! -------   ----     -------
  !   1.0    09/2002      Initial version     (F. Chevallier)
  !   1.1    05/2003      RTTOV7.3 compatible (F. Chevallier)
  !   1.2    03/2004      Added polarimetry   (R. Saunders)
  !   1.3    08/2004      Polarimetry fixes   (U. O'Keeffe)
  !   1.4    11/2004      Clean-up            (P. Bauer)
  !   1.5    07/2005      Polarimetry fixes   (U. O'Keeffe)
  !   1.6    11/2005      Add errorstatus to iniscatt arguments and use a temporary
  !                       radiance type for the calcpolarisation call (J Cameron)
  !   1.7    11/2007      RTTOV9 version      (A. Geer)
  !   1.8    07/2008      Clear sky speed-ups (A. Geer)
  !   1.9    10/2008      Revised cloud partitioning (A. Geer)
  !   1.10   11/2009      User may supply the average cloud fraction (A. Geer)
  !   1.11   11/2009      Adapted for RTTOV10 (A. Geer)
  !   1.12   04/2010      Tidied up after code cleaning (A. Geer)
  !   1.13   02/2010      Revised cloud partitioning is now the default (A. Geer)
  !   1.14   08/2010      Fix for polarimetric channels until they can be 
  !                       hanled properly by RTTOV_SCATT (W. Bell)
  !   1.15   06/2015      Allow foam_fraction input (L-F Meunier)
  !
  Use rttov_types, Only :    &
       & rttov_coefs          ,&
       & rttov_scatt_coef     ,&
       & profile_Type         ,&
       & profile_cloud_Type   ,&
       & radiance_Type        ,&
       & rttov_chanprof       ,&
       & rttov_emissivity     ,&
       & rttov_options_scatt

  Use parkind1, Only : jpim, jplm
!INTF_OFF
  Use parkind1, Only : jprb

  Use rttov_types, Only :    &
       & geometry_Type        ,&
       & profile_scatt_aux    ,&
       & transmission_Type    ,&
       & rttov_options

  Use rttov_const, Only :   &
       & errorstatus_success ,&
       & errorstatus_fatal, &
       & ccthres          ,&
       & sensor_id_po

  USE YOMHOOK, ONLY: LHOOK , DR_HOOK
!INTF_ON
  Implicit None


!* Subroutine arguments:
  Type(rttov_options_scatt), Intent(in) :: opts_scatt              ! RTTOV-SCATT options
  Integer (Kind=jpim), Intent (in)      :: nlevels                 ! Number of levels
  Type(rttov_chanprof),Intent (in)      :: chanprof(:)             ! Indices
  Type (profile_Type), Intent (in)      :: profiles(:)             ! Atmospheric profiles
  Integer (Kind=jpim), Intent (in)      :: frequencies (size(chanprof)) ! Frequency indices
  Integer (Kind=jpim), Intent (out)     :: errorstatus                  ! Error return flag

  Logical (Kind=jplm),    Intent (in)    :: calcemiss     (size(chanprof))  ! Switch for emissivity calculation
  Type(rttov_emissivity), Intent (inout) :: emissivity    (size(chanprof))  ! Surface emissivity
  Type(rttov_emissivity), Intent (inout) :: emissivity_tl (size(chanprof))  ! Surface emissivity
  
  Type (profile_Type),        Intent (in)    :: profiles_tl     (size(profiles))
  Type (rttov_coefs),         Intent (in)    :: coef_rttov                  ! RTTOV Coefficients
  Type (rttov_scatt_coef),    Intent (in)    :: coef_scatt                  ! RTTOV_SCATT Coefficients
  Type (profile_cloud_Type),  Intent (in)    :: cld_profiles    (size(profiles)) ! Cloud profiles 
  Type (profile_cloud_Type),  Intent (in)    :: cld_profiles_tl (size(profiles))   
  Type (radiance_Type),       Intent (inout) :: radiance                    ! Radiances
  Type (radiance_Type),       Intent (inout) :: radiance_tl          

!INTF_END

#include "rttov_tl.interface"
#include "rttov_iniscatt_tl.interface"
#include "rttov_eddington_tl.interface"
#include "rttov_errorreport.interface"
 
  Integer (Kind=jpim), target :: sa__mclayer    (size(chanprof))
  Integer (Kind=jpim), target :: sa_tl__mclayer (size(chanprof))
   
  Real (Kind=jprb), target :: t__tau_total     (size(chanprof))
  Real (Kind=jprb), target :: t__tau_levels    (nlevels,size(chanprof))
  Real (Kind=jprb), target :: t_tl__tau_total  (size(chanprof))
  Real (Kind=jprb), target :: t_tl__tau_levels (nlevels,size(chanprof))
  
  Real (Kind=jprb), target :: sa__cfrac   (size(profiles))  
  Real (Kind=jprb), target :: sa__ems_bnd (size(chanprof))
  Real (Kind=jprb), target :: sa__ref_bnd (size(chanprof))
  Real (Kind=jprb), target :: sa__ems_cld (size(chanprof))
  Real (Kind=jprb), target :: sa__ref_cld (size(chanprof))
  
  Real (Kind=jprb), target :: sa__tbd (size(profiles),nlevels+1)
  
  Real (Kind=jprb), target :: sa__delta  (size(chanprof),nlevels)
  Real (Kind=jprb), target :: sa__tau    (size(chanprof),nlevels)
  Real (Kind=jprb), target :: sa__ext    (size(chanprof),nlevels)
  Real (Kind=jprb), target :: sa__ssa    (size(chanprof),nlevels)
  Real (Kind=jprb), target :: sa__asm    (size(chanprof),nlevels)
  Real (Kind=jprb), target :: sa__lambda (size(chanprof),nlevels)
  Real (Kind=jprb), target :: sa__h      (size(chanprof),nlevels)
  
  Real (Kind=jprb), target :: sa__b0     (size(profiles),nlevels)
  Real (Kind=jprb), target :: sa__b1     (size(profiles),nlevels)
  Real (Kind=jprb), target :: sa__bn     (size(profiles),nlevels)
  Real (Kind=jprb), target :: sa__dz     (size(profiles),nlevels)
  Real (Kind=jprb), target :: sa__clw    (size(profiles),nlevels)
  Real (Kind=jprb), target :: sa__ciw    (size(profiles),nlevels)
  Real (Kind=jprb), target :: sa__totalice (size(profiles),nlevels)
  Real (Kind=jprb), target :: sa__rain   (size(profiles),nlevels)
  Real (Kind=jprb), target :: sa__sp     (size(profiles),nlevels)

  Real (Kind=jprb), target :: sa_tl__cfrac   (size(profiles))
  Real (Kind=jprb), target :: sa_tl__ems_bnd (size(chanprof))
  Real (Kind=jprb), target :: sa_tl__ref_bnd (size(chanprof))
  Real (Kind=jprb), target :: sa_tl__ems_cld (size(chanprof))
  Real (Kind=jprb), target :: sa_tl__ref_cld (size(chanprof))
  
  Real (Kind=jprb), target :: sa_tl__tbd (size(profiles),nlevels+1)
  
  Real (Kind=jprb), target :: sa_tl__delta  (size(chanprof),nlevels)
  Real (Kind=jprb), target :: sa_tl__tau    (size(chanprof),nlevels)
  Real (Kind=jprb), target :: sa_tl__ext    (size(chanprof),nlevels)
  Real (Kind=jprb), target :: sa_tl__ssa    (size(chanprof),nlevels)
  Real (Kind=jprb), target :: sa_tl__asm    (size(chanprof),nlevels)
  Real (Kind=jprb), target :: sa_tl__lambda (size(chanprof),nlevels)
  Real (Kind=jprb), target :: sa_tl__h      (size(chanprof),nlevels)
  
  Real (Kind=jprb), target :: sa_tl__b0     (size(profiles),nlevels)
  Real (Kind=jprb), target :: sa_tl__b1     (size(profiles),nlevels)
  Real (Kind=jprb), target :: sa_tl__bn     (size(profiles),nlevels)
  Real (Kind=jprb), target :: sa_tl__dz     (size(profiles),nlevels)
  Real (Kind=jprb), target :: sa_tl__clw    (size(profiles),nlevels)
  Real (Kind=jprb), target :: sa_tl__ciw    (size(profiles),nlevels)
  Real (Kind=jprb), target :: sa_tl__totalice (size(profiles),nlevels)
  Real (Kind=jprb), target :: sa_tl__rain   (size(profiles),nlevels)
  Real (Kind=jprb), target :: sa_tl__sp     (size(profiles),nlevels)

!* Local variables:
  Integer (Kind=jpim) :: nprofiles, nchannels
  Logical (Kind=jplm) :: lpolarimetric(size(chanprof))
  Integer (Kind=jpim) :: iprof, ichan, ichan_act 
  Real    (Kind=jprb) :: cld_bt        (size(chanprof))            
  Real    (Kind=jprb) :: cld_bt_tl     (size(chanprof))  

  Type (transmission_Type) :: transmission, transmission_tl
  Type (geometry_Type)     :: angles (size(profiles))
  Type (profile_scatt_aux) :: scatt_aux, scatt_aux_tl

  Type(rttov_options)    :: opts

  Character (len=80) :: errMessage
  Character (len=15) :: NameOfRoutine = 'rttov_scatt_tl '
  
  Real (KIND=JPRB) :: ZHOOK_HANDLE        

  !- End of header --------------------------------------------------------
  
  if (lhook) call dr_hook('RTTOV_SCATT_TL',0_jpim,zhook_handle)

  nprofiles = size(profiles)
  nchannels = size(chanprof)

  errorstatus = errorstatus_success

  transmission % tau_total  => t__tau_total
  transmission % tau_levels => t__tau_levels

  transmission_tl % tau_total  => t_tl__tau_total
  transmission_tl % tau_levels => t_tl__tau_levels

  scatt_aux % cfrac    => sa__cfrac
  scatt_aux % ems_bnd  => sa__ems_bnd
  scatt_aux % ref_bnd  => sa__ref_bnd
  scatt_aux % ems_cld  => sa__ems_cld
  scatt_aux % ref_cld  => sa__ref_cld
  scatt_aux % tbd      => sa__tbd
  scatt_aux % mclayer  => sa__mclayer
  scatt_aux % delta    => sa__delta
  scatt_aux % tau      => sa__tau
  scatt_aux % ext      => sa__ext
  scatt_aux % ssa      => sa__ssa
  scatt_aux % asm      => sa__asm
  scatt_aux % lambda   => sa__lambda
  scatt_aux % h        => sa__h
  scatt_aux % b0       => sa__b0
  scatt_aux % b1       => sa__b1
  scatt_aux % bn       => sa__bn
  scatt_aux % dz       => sa__dz
  scatt_aux % clw      => sa__clw
  scatt_aux % ciw      => sa__ciw
  scatt_aux % totalice => sa__totalice
  scatt_aux % rain     => sa__rain
  scatt_aux % sp       => sa__sp

  scatt_aux_tl % cfrac    => sa_tl__cfrac
  scatt_aux_tl % ems_bnd  => sa_tl__ems_bnd
  scatt_aux_tl % ref_bnd  => sa_tl__ref_bnd
  scatt_aux_tl % ems_cld  => sa_tl__ems_cld
  scatt_aux_tl % ref_cld  => sa_tl__ref_cld
  scatt_aux_tl % tbd      => sa_tl__tbd
  scatt_aux_tl % mclayer  => sa_tl__mclayer
  scatt_aux_tl % delta    => sa_tl__delta
  scatt_aux_tl % tau      => sa_tl__tau
  scatt_aux_tl % ext      => sa_tl__ext
  scatt_aux_tl % ssa      => sa_tl__ssa
  scatt_aux_tl % asm      => sa_tl__asm
  scatt_aux_tl % lambda   => sa_tl__lambda
  scatt_aux_tl % h        => sa_tl__h
  scatt_aux_tl % b0       => sa_tl__b0
  scatt_aux_tl % b1       => sa_tl__b1
  scatt_aux_tl % bn       => sa_tl__bn
  scatt_aux_tl % dz       => sa_tl__dz
  scatt_aux_tl % clw      => sa_tl__clw
  scatt_aux_tl % ciw      => sa_tl__ciw
  scatt_aux_tl % totalice => sa_tl__totalice
  scatt_aux_tl % rain     => sa_tl__rain
  scatt_aux_tl % sp       => sa_tl__sp

  ! Check inputs
  ! ------------
  Do iprof = 1, nprofiles
    If (  profiles(iprof) % s2m % p /= cld_profiles(iprof) % ph(nlevels+1)  ) Then
      errorstatus = errorstatus_fatal
      Write( errMessage, '( "Surface pressure and lowest half level should be identical")' )
      Call rttov_errorreport (errorstatus_fatal, errMessage, NameOfRoutine)
      IF (LHOOK) CALL DR_HOOK('RTTOV_SCATT',1_jpim,ZHOOK_HANDLE)
      Return
    End If
  End Do

  ! Identify polarimetric channels for fix to use only the clear sky part of RT calculation 
  do ichan = 1, nchannels
    ichan_act = chanprof(ichan)%chan
    lpolarimetric(ichan) = ( (coef_rttov % coef % id_sensor == sensor_id_po) &
                   & .and.   (coef_rttov % coef % fastem_polar(ichan_act) + 1_jpim >= 6_jpim) )
  enddo

!*         1.   Gas absorption

  ! Profiles will be interpolated from model/RTTOV-SCATT levels to 
  ! RTTOV coefficient levels within RTTOV itself.
  opts%interpolation%addinterp   = .true.
  opts%rt_ir%addclouds           = .false.
  opts%rt_ir%addsolar            = .false.
  opts%rt_ir%addaerosl           = .false.
  opts%rt_ir%pc%addpc            = .false.
  opts%rt_ir%pc%addradrec        = .false.

  opts%rt_mw%fastem_version           = opts_scatt%fastem_version
  opts%rt_mw%supply_foam_fraction     = opts_scatt%supply_foam_fraction
  opts%rt_all%use_q2m                 = opts_scatt%use_q2m
  opts%config                         = opts_scatt%config
  opts%interpolation%interp_mode      = opts_scatt%interp_mode
  opts%interpolation%reg_limit_extrap = opts_scatt%reg_limit_extrap

  Call rttov_tl(             &
     & errorstatus,          &! out
     & chanprof,             &! in
     & opts,                 &! in
     & profiles,             &! in
     & profiles_tl,          &! in
     & coef_rttov,           &! in
     & transmission,         &! inout
     & transmission_tl,      &! inout
     & radiance,             &! inout
     & radiance_tl,          &! inout
     & calcemiss,            &! in
     & emissivity,           &! inout
     & emissivity_tl)         ! inout

  If ( errorstatus == errorstatus_fatal ) Then
     Write( errMessage, '( "error in rttov_tl")' )
     Call rttov_errorreport (errorstatus_fatal, errMessage, NameOfRoutine)
     IF (LHOOK) CALL DR_HOOK('RTTOV_SCATT_TL',1_jpim,ZHOOK_HANDLE)
     Return
  End If
 
  scatt_aux_tl % ems_cld (:) = emissivity_tl (:) % emis_in
  scatt_aux    % ems_cld (:) = emissivity    (:) % emis_in
  scatt_aux_tl % ref_cld (:) = -1.0_JPRB * emissivity_tl (:) % emis_in
  scatt_aux    % ref_cld (:) =  1.0_JPRB - emissivity    (:) % emis_in

!*  2.   Initialisations for Eddington
  Call rttov_iniscatt_tl(        &
        & errorstatus,           &! out
        & opts,                  &! in
        & nlevels,               &! in
        & nchannels,             &! in
        & nprofiles,             &! in
        & chanprof,              &! in
        & frequencies,           &! in
        & profiles,              &! in
        & profiles_tl,           &! in
        & cld_profiles,          &! in
        & cld_profiles_tl,       &! in
        & coef_rttov%coef,       &! in
        & coef_scatt,            &! in
        & transmission,          &! in
        & transmission_tl,       &! in
        & calcemiss,             &! in
        & opts_scatt%lusercfrac, &! in
        & angles,                &! out
        & scatt_aux,             &! inout
        & scatt_aux_tl)           ! inout

  If ( errorstatus == errorstatus_fatal ) Then
     Write( errMessage, '( "error in rttov_iniscatt_tl")' )
     Call rttov_errorreport (errorstatus_fatal, errMessage, NameOfRoutine)
     IF (LHOOK) CALL DR_HOOK('RTTOV_SCATT_TL',1_jpim,ZHOOK_HANDLE)
     Return
  End If

!* 3.   Eddington (in temperature space)
  Call rttov_eddington_tl(   &
        & nlevels,        &! in
        & nchannels,         &! in
        & nprofiles,         &! in
        & chanprof%prof,     &! in
        & angles,            &! in
        & profiles,          &! in  
        & profiles_tl,       &! in  
        & scatt_aux,         &! in
        & scatt_aux_tl,      &! in
        & cld_bt,            &! out  
        & cld_bt_tl)          ! out   

!*  4.   Combine clear and cloudy parts

  Do ichan = 1, nchannels
    iprof = chanprof (ichan) % prof
    if (scatt_aux % cfrac (iprof) > ccthres .and. .not. lpolarimetric(ichan)) then      
     
      radiance_tl % bt (ichan) = &
       &   scatt_aux % cfrac (iprof)                      * cld_bt_tl (ichan)  & 
       & + (cld_bt (ichan) - radiance % bt_clear (ichan)) * scatt_aux_tl % cfrac (iprof) &
       & + (1.0_JPRB - scatt_aux % cfrac (iprof))         * radiance_tl % bt_clear (ichan) 

      radiance % bt (ichan) = cld_bt (ichan) * scatt_aux % cfrac (iprof)   & 
       & + radiance % bt_clear (ichan) * (1.0_JPRB - scatt_aux % cfrac (iprof)) 

    else

      radiance_tl % bt (ichan) = radiance_tl % bt_clear (ichan) 
      radiance    % bt (ichan) = radiance    % bt_clear (ichan) 

    endif

  End Do
   
  if (lhook) call dr_hook('RTTOV_SCATT_TL',1_jpim,zhook_handle) 
  
End Subroutine rttov_scatt_tl
