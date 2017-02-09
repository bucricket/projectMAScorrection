! Description:
!> @file
!!   Runs RTTOV-SCATT adjoint (AD) or Jacobian (K) models
!
!> @brief
!!   Runs RTTOV-SCATT adjoint (AD) or Jacobian (K) models
!!
!! @details
!!   The same interface is used for the adjoint and Jacobian models.
!!
!!   If the profiles_ad(:) array is the same size as the profiles(:)
!!   array the adjoint model is run. Given the gradient of a scalar
!!   function with respect to channel BTs, this computes the gradient
!!   of the same scalar function with respect to the profile variables.
!!   This is the transpose of the tangent linear model.
!!
!!   Alternatively if the profiles_ad(:) structure is the same size
!!   as the chanprof(:) array the Jacobian model is run. This calculates
!!   the derivative of the satellite-seen BTs with respect to each profile
!!   variable in profiles_ad and cld_profiles_ad for each channel.
!!
!!   In both cases cld_profiles_ad(:) must be the same size as profiles_ad(:).
!!
!!   For RTTOV-SCATT the input radiance perturbation is always specified in BT
!!   in radiance_ad%bt(:).
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
!! @param[in,out] profiles_ad      output gradient wrt atmospheric profile and surface variables
!! @param[in,out] cld_profiles_ad  output gradient wrt cloud and hydrometeor profile
!! @param[in,out] emissivity_ad    output gradient wrt surface emissivities
!! @param[in,out] radiance         output radiances and corresponding BTs
!! @param[in,out] radiance_ad      input perturbations wrt BTs
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
Subroutine rttov_scatt_ad( & 
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
     & profiles_ad,        &! inout
     & cld_profiles_ad,    &! inout
     & emissivity_ad,      &! inout
     & radiance,           &! inout
     & radiance_ad)         ! inout
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
  !   1.15   11/2012      RTTOV11; remove "old cloud fraction" approach (A. Geer)
  !   1.16   06/2015      Allow foam_fraction input (L-F Meunier)

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
       & errorstatus_success, &
       & errorstatus_fatal,   &
       & adk_adjoint,         &
       & adk_k,               &
       & ccthres,             &
       & sensor_id_po

  USE YOMHOOK, ONLY: LHOOK , DR_HOOK
!INTF_ON
  Implicit None


!* Subroutine arguments:
  Type(rttov_options_scatt), Intent(in) :: opts_scatt                   ! RTTOV-SCATT options
  Integer (Kind=jpim), Intent (in)      :: nlevels                      ! Number of levels
  Type(rttov_chanprof),Intent (in)      :: chanprof(:)                  ! Indices
  Type (profile_Type), Intent (in)      :: profiles(:)                  ! Atmospheric profiles
  Type (profile_Type), Intent (inout)   :: profiles_ad(:)    
  Integer (Kind=jpim), Intent (in)      :: frequencies (size(chanprof)) ! Frequency indices
  Integer (Kind=jpim), Intent (out)     :: errorstatus                  ! Error return flag

  Logical (Kind=jplm),    Intent (in)    :: calcemiss     (size(chanprof))  ! Switch for emissivity calculation
  Type(rttov_emissivity), Intent (inout) :: emissivity    (size(chanprof))  ! Surface emissivity
  Type(rttov_emissivity), Intent (inout) :: emissivity_ad (size(chanprof))  ! Surface emissivity

  Type (rttov_coefs),         Intent (in)    :: coef_rttov                     ! RTTOV Coefficients
  Type (rttov_scatt_coef),    Intent (in)    :: coef_scatt                     ! RTTOV_SCATT Coefficients
  Type (profile_cloud_Type),  Intent (in)    :: cld_profiles    (size(profiles))    ! Cloud profiles 
  Type (profile_cloud_Type),  Intent (inout) :: cld_profiles_ad (size(profiles_ad))   
  Type (radiance_Type),       Intent (inout) :: radiance         ! Radiances
  Type (radiance_Type),       Intent (inout) :: radiance_ad          

!INTF_END

#include "rttov_direct.interface"
#include "rttov_iniscatt.interface"
#include "rttov_eddington.interface"
#include "rttov_ad.interface"
#include "rttov_k.interface"
#include "rttov_iniscatt_ad.interface"
#include "rttov_eddington_ad.interface"
#include "rttov_errorreport.interface"
   
  Integer (Kind=jpim), target :: sa__mclayer    (size(chanprof))
  Integer (Kind=jpim), target :: sa_ad__mclayer (size(chanprof))
      
  Real (Kind=jprb), target :: t__tau_total     (size(chanprof))
  Real (Kind=jprb), target :: t__tau_levels    (nlevels,size(chanprof))
  Real (Kind=jprb), target :: t_ad__tau_total  (size(chanprof))
  Real (Kind=jprb), target :: t_ad__tau_levels (nlevels,size(chanprof))
  
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

  Real (Kind=jprb), target :: sa_ad__cfrac   (size(profiles_ad))
  Real (Kind=jprb), target :: sa_ad__ems_bnd (size(chanprof))
  Real (Kind=jprb), target :: sa_ad__ref_bnd (size(chanprof))
  Real (Kind=jprb), target :: sa_ad__ems_cld (size(chanprof))
  Real (Kind=jprb), target :: sa_ad__ref_cld (size(chanprof))
  
  Real (Kind=jprb), target :: sa_ad__tbd (size(profiles_ad),nlevels+1)
  
  Real (Kind=jprb), target :: sa_ad__delta  (size(chanprof),nlevels)
  Real (Kind=jprb), target :: sa_ad__tau    (size(chanprof),nlevels)
  Real (Kind=jprb), target :: sa_ad__ext    (size(chanprof),nlevels)
  Real (Kind=jprb), target :: sa_ad__ssa    (size(chanprof),nlevels)
  Real (Kind=jprb), target :: sa_ad__asm    (size(chanprof),nlevels)
  Real (Kind=jprb), target :: sa_ad__lambda (size(chanprof),nlevels)
  Real (Kind=jprb), target :: sa_ad__h      (size(chanprof),nlevels)
  
  Real (Kind=jprb), target :: sa_ad__b0     (size(profiles_ad),nlevels)
  Real (Kind=jprb), target :: sa_ad__b1     (size(profiles_ad),nlevels)
  Real (Kind=jprb), target :: sa_ad__bn     (size(profiles_ad),nlevels)
  Real (Kind=jprb), target :: sa_ad__dz     (size(profiles_ad),nlevels)
  Real (Kind=jprb), target :: sa_ad__clw    (size(profiles_ad),nlevels)
  Real (Kind=jprb), target :: sa_ad__ciw    (size(profiles_ad),nlevels)
  Real (Kind=jprb), target :: sa_ad__totalice (size(profiles_ad),nlevels)
  Real (Kind=jprb), target :: sa_ad__rain   (size(profiles_ad),nlevels)
  Real (Kind=jprb), target :: sa_ad__sp     (size(profiles_ad),nlevels)

!* Local variables:
  Integer (Kind=jpim) :: nprofiles, nprofilesad, nchannels
  Logical (Kind=jplm) :: lpolarimetric(size(chanprof))
  Integer (Kind=jpim) :: iprof, ichan, ichan_act
  Integer (Kind=jpim) :: iprofad, adk  

  Real (Kind=jprb), Dimension (size(chanprof)) :: cld_bt, cld_bt_ad

  Type (transmission_Type)   :: transmission, transmission_ad
  Type (geometry_Type)       :: angles (size(profiles))
  Type (profile_scatt_aux)   :: scatt_aux, scatt_aux_ad   

  Type(rttov_options)    :: opts

  Character (len=80) :: errMessage
  Character (len=15) :: NameOfRoutine = 'rttov_scatt_ad '

  REAL(KIND=JPRB) :: ZHOOK_HANDLE
  
  !- End of header --------------------------------------------------------
  
  if (lhook) call dr_hook('RTTOV_SCATT_AD',0_jpim,zhook_handle)        

  nprofiles   = size(profiles)
  nprofilesad = size(profiles_ad)
  nchannels   = size(chanprof)

  errorstatus = errorstatus_success

  if (nprofilesad == nprofiles) then 
     adk = adk_adjoint   ! Adjoint mode
  else if (nprofilesad == nchannels) then
     adk = adk_k         ! K mode
  else
     Write( errMessage, '( "incorrect number of profiles in adjoint/K")' )
     Call rttov_errorreport (errorstatus_fatal, errMessage, NameOfRoutine)
     IF (LHOOK) CALL DR_HOOK('RTTOV_SCATT_AD',1_jpim,ZHOOK_HANDLE)
     Return
  endif 

  transmission % tau_total  => t__tau_total
  transmission % tau_levels => t__tau_levels

  transmission_ad % tau_total  => t_ad__tau_total
  transmission_ad % tau_levels => t_ad__tau_levels

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

  scatt_aux_ad % cfrac    => sa_ad__cfrac
  scatt_aux_ad % ems_bnd  => sa_ad__ems_bnd
  scatt_aux_ad % ref_bnd  => sa_ad__ref_bnd
  scatt_aux_ad % ems_cld  => sa_ad__ems_cld
  scatt_aux_ad % ref_cld  => sa_ad__ref_cld
  scatt_aux_ad % tbd      => sa_ad__tbd
  scatt_aux_ad % mclayer  => sa_ad__mclayer
  scatt_aux_ad % delta    => sa_ad__delta
  scatt_aux_ad % tau      => sa_ad__tau
  scatt_aux_ad % ext      => sa_ad__ext
  scatt_aux_ad % ssa      => sa_ad__ssa
  scatt_aux_ad % asm      => sa_ad__asm
  scatt_aux_ad % lambda   => sa_ad__lambda
  scatt_aux_ad % h        => sa_ad__h
  scatt_aux_ad % b0       => sa_ad__b0
  scatt_aux_ad % b1       => sa_ad__b1
  scatt_aux_ad % bn       => sa_ad__bn
  scatt_aux_ad % dz       => sa_ad__dz
  scatt_aux_ad % clw      => sa_ad__clw
  scatt_aux_ad % ciw      => sa_ad__ciw
  scatt_aux_ad % totalice => sa_ad__totalice
  scatt_aux_ad % rain     => sa_ad__rain
  scatt_aux_ad % sp       => sa_ad__sp

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
  opts%rt_all%switchrad          = .true.   ! Required, in AD only, and for some historical reason

  opts%rt_mw%fastem_version           = opts_scatt%fastem_version
  opts%rt_mw%supply_foam_fraction     = opts_scatt%supply_foam_fraction
  opts%rt_all%use_q2m                 = opts_scatt%use_q2m
  opts%config                         = opts_scatt%config
  opts%interpolation%interp_mode      = opts_scatt%interp_mode
  opts%interpolation%reg_limit_extrap = opts_scatt%reg_limit_extrap

  Call rttov_direct(               &!
         & errorstatus,            &! out
         & chanprof,               &! in
         & opts,                   &! in
         & profiles,               &! in
         & coef_rttov,             &! in
         & transmission,           &! inout
         & radiance,               &! inout
         & calcemis   = calcemiss, &! in
         & emissivity = emissivity) ! inout

  If ( errorstatus == errorstatus_fatal ) Then
     Write( errMessage, '( "error in rttov_direct")' )
     Call rttov_errorreport (errorstatus_fatal, errMessage, NameOfRoutine)
     IF (LHOOK) CALL DR_HOOK('RTTOV_SCATT_AD',1_jpim,ZHOOK_HANDLE)
     Return
  End If

  scatt_aux % ems_cld (:) = emissivity (:) % emis_in
  scatt_aux % ref_cld (:) = 1.0_JPRB - emissivity (:) % emis_in
  
!*  2.   Initialisations for Eddington
  Call rttov_iniscatt(           &
        & errorstatus,           &! out
        & opts,                  &! in
        & nlevels,               &! in
        & nchannels,             &! in
        & nprofiles,             &! in
        & chanprof,              &! in
        & frequencies,           &! in
        & profiles,              &! in  
        & cld_profiles,          &! in
        & coef_rttov%coef,       &! in  
        & coef_scatt,            &! in  
        & transmission,          &! in
        & calcemiss,             &! in
        & opts_scatt%lusercfrac, &! in
        & angles,                &! out
        & scatt_aux   )           ! inout   

  If ( errorstatus == errorstatus_fatal ) Then
     Write( errMessage, '( "error in rttov_iniscatt")' )
     Call rttov_errorreport (errorstatus_fatal, errMessage, NameOfRoutine)
     IF (LHOOK) CALL DR_HOOK('RTTOV_SCATT_AD',1_jpim,ZHOOK_HANDLE)
     Return
  End If

!* 3.   Eddington (in temperature space)
   Call rttov_eddington(     &
        & nlevels,           &! in
        & nchannels,         &! in
        & nprofiles,         &! in
        & chanprof%prof,     &! in
        & angles,            &! in
        & profiles,          &! in  
        & scatt_aux,         &! in
        & cld_bt)             ! out   

!*  4.   Combine clear and cloudy parts
  Do ichan = 1, nchannels
    iprof = chanprof(ichan) % prof

    if (scatt_aux % cfrac (iprof) > ccthres .and. .not. lpolarimetric(ichan)) then 

      radiance % bt (ichan) = cld_bt (ichan) * scatt_aux % cfrac (iprof) & 
        & + radiance % bt_clear (ichan) * (1.0_JPRB - scatt_aux % cfrac (iprof))

    else

      radiance % bt (ichan) = radiance % bt_clear (ichan)

    endif
  End Do
    
!* ADJOINT PART
  
  scatt_aux_ad % cfrac   (:)   = 0.0_JPRB
  scatt_aux_ad % ems_bnd (:)   = 0.0_JPRB
  scatt_aux_ad % ref_bnd (:)   = 0.0_JPRB
  scatt_aux_ad % ems_cld (:)   = 0.0_JPRB
  scatt_aux_ad % ref_cld (:)   = 0.0_JPRB
  scatt_aux_ad % tbd     (:,:) = 0.0_JPRB
  scatt_aux_ad % delta   (:,:) = 0.0_JPRB
  scatt_aux_ad % tau     (:,:) = 0.0_JPRB
  scatt_aux_ad % ext     (:,:) = 0.0_JPRB
  scatt_aux_ad % ssa     (:,:) = 0.0_JPRB
  scatt_aux_ad % asm     (:,:) = 0.0_JPRB
  scatt_aux_ad % lambda  (:,:) = 0.0_JPRB
  scatt_aux_ad % h       (:,:) = 0.0_JPRB
  scatt_aux_ad % b0      (:,:) = 0.0_JPRB
  scatt_aux_ad % b1      (:,:) = 0.0_JPRB
  scatt_aux_ad % bn      (:,:) = 0.0_JPRB
  scatt_aux_ad % dz      (:,:) = 0.0_JPRB
  scatt_aux_ad % clw     (:,:) = 0.0_JPRB
  scatt_aux_ad % ciw     (:,:) = 0.0_JPRB
  scatt_aux_ad % totalice (:,:) = 0.0_JPRB
  scatt_aux_ad % rain    (:,:) = 0.0_JPRB
  scatt_aux_ad % sp      (:,:) = 0.0_JPRB

  transmission_ad % tau_total       (:)   = 0.0_JPRB
  transmission_ad % tau_levels      (:,:) = 0.0_JPRB 

  emissivity_ad(:)%emis_out=0.0_JPRB 

  cld_bt_ad (:) = 0.0_JPRB

!* 4.   Combine clear and cloudy parts  
  Do ichan = 1, nchannels
    iprof = chanprof (ichan) % prof

    if (adk == adk_adjoint) then
      iprofad = iprof
    else if (adk == adk_k) then
      iprofad = ichan
    endif

    if (scatt_aux % cfrac (iprof) > ccthres .and. .not. lpolarimetric(ichan)) then 

      cld_bt_ad (ichan) = cld_bt_ad (ichan) + &
       & scatt_aux % cfrac (iprof) * radiance_ad % bt (ichan) 
      
      scatt_aux_ad % cfrac (iprofad) = scatt_aux_ad % cfrac (iprofad) + &
       & (cld_bt (ichan) - radiance % bt_clear (ichan)) * radiance_ad % bt (ichan) 

      radiance_ad  % bt_clear (ichan) = radiance_ad % bt_clear (ichan) + &
       & (1.0_JPRB - scatt_aux % cfrac (iprof)) * radiance_ad % bt (ichan)

      radiance_ad % bt (ichan) = 0.0_JPRB

    else

      radiance_ad % bt_clear (ichan) = radiance_ad % bt_clear (ichan) + &
       & radiance_ad % bt (ichan)

      radiance_ad % bt (ichan) = 0.0_JPRB

    endif

  End Do

!* 3.   Eddington (in temperature space)
  Call rttov_eddington_ad(   &
        & nlevels,           &! in
        & nchannels,         &! in
        & nprofiles,         &! in
        & nprofilesad,       &! in
        & chanprof % prof,   &! in
        & angles,            &! in
        & profiles,          &! in  
        & profiles_ad,       &! inout  
        & scatt_aux,         &! in
        & scatt_aux_ad,      &! inout
        & cld_bt,            &! out  
        & cld_bt_ad)          ! inout   

!*  2.   Initialisations for Eddington  
  Call rttov_iniscatt_ad(        &
        & errorstatus,           &! out
        & opts,                  &! in
        & nlevels,               &! in
        & nchannels,             &! in
        & nprofiles,             &! in
        & nprofilesad,           &! in
        & chanprof,              &! in
        & frequencies,           &! in
        & profiles,              &! in
        & profiles_ad,           &! inout
        & cld_profiles,          &! in
        & cld_profiles_ad,       &! inout
        & coef_rttov%coef,       &! in
        & coef_scatt,            &! in
        & transmission,          &! in
        & transmission_ad,       &! inout
        & calcemiss,             &! in
        & opts_scatt%lusercfrac, &! in
        & angles,                &! out
        & scatt_aux,             &! inout
        & scatt_aux_ad)           ! inout

  If ( errorstatus == errorstatus_fatal ) Then
     Write( errMessage, '( "error in rttov_iniscatt_ad")' )
     Call rttov_errorreport (errorstatus_fatal, errMessage, NameOfRoutine)
     IF (LHOOK) CALL DR_HOOK('RTTOV_SCATT_AD',1_jpim,ZHOOK_HANDLE)
     Return
  End If
  
  emissivity_ad(:)%emis_in = emissivity_ad(:)%emis_in - scatt_aux_ad % ref_cld (:)
  scatt_aux_ad % ref_cld (:) = 0.0_JPRB      
  emissivity_ad(:)%emis_in = emissivity_ad(:)%emis_in + scatt_aux_ad % ems_cld (:)
  scatt_aux_ad % ems_cld (:) = 0.0_JPRB   
  
!*         1.   Gas absorption

  if (adk == adk_adjoint) then  
    Call rttov_ad(        &
       & errorstatus,     &! out
       & chanprof,        &! in
       & opts,            &
       & profiles,        &! in
       & profiles_ad,     &! inout
       & coef_rttov,      &! in
       & transmission,    &! inout
       & transmission_ad, &! inout
       & radiance,        &! inout
       & radiance_ad,     &! inout 
       & calcemiss,       &! in
       & emissivity,      &! inout
       & emissivity_ad)    ! inout
  else if (adk == adk_k) then 
    Call rttov_k(        &
       & errorstatus,     &! out
       & chanprof,        &! in
       & opts,            &! in
       & profiles,        &! in
       & profiles_ad,     &! inout
       & coef_rttov,      &! in
       & transmission,    &! inout
       & transmission_ad, &! inout
       & radiance,        &! inout
       & radiance_ad,     &! inout 
       & calcemiss,       &! in
       & emissivity,      &! inout
       & emissivity_ad)    ! inout
  endif

  If ( errorstatus == errorstatus_fatal ) Then
     Write( errMessage, '( "error in rttov_ad")' )
     Call rttov_errorreport (errorstatus_fatal, errMessage, NameOfRoutine)
     IF (LHOOK) CALL DR_HOOK('RTTOV_SCATT_AD',1_jpim,ZHOOK_HANDLE)
     Return
  End If

!*  4b.   Combine clear and cloudy parts - again - because rttov_ad/k has
!*        overwritten the bt array with clear-sky values.
  Do ichan = 1, nchannels
    iprof = chanprof (ichan) % prof

    if (scatt_aux % cfrac (iprof) > ccthres .and. .not. lpolarimetric(ichan)) then 

      radiance % bt (ichan) = cld_bt (ichan) * scatt_aux % cfrac (iprof) &
       & + radiance % bt_clear (ichan) * (1.0_JPRB - scatt_aux % cfrac (iprof))

    else

      radiance % bt (ichan) = radiance % bt_clear (ichan)

    endif

  End Do

  if (lhook) call dr_hook('RTTOV_SCATT_AD',1_jpim,zhook_handle)

End Subroutine rttov_scatt_ad
