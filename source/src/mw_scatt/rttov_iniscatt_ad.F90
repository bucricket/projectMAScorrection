!
Subroutine rttov_iniscatt_ad (&
      & errorstatus,       &! out
      & opts,              &! in
      & nlevels,           &! in
      & nchannels,         &! in
      & nprofiles,         &! in
      & nprofilesad,       &! in
      & chanprof,          &! in
      & frequencies,       &! in
      & profiles,          &! in  
      & profiles_ad,       &! inout
      & cld_profiles,      &! in 
      & cld_profiles_ad,   &! inout 
      & coef_rttov,        &! in
      & coef_scatt,        &! in
      & transmission,      &! in
      & transmission_ad,   &! inout
      & calcemiss,         &! in
      & usercfrac,         &! in
      & angles,            &! out
      & scatt_aux,         &! inout
      & scatt_aux_ad)       ! inout 

  !
  ! Description:
  ! AD of routine to
  ! Calculate some variables related to the input precipitation profile
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
  !    Copyright 2002, EUMETSAT, All Rights Reserved.
  !
  ! Method:
  !
  ! Current Code Owner: SAF NWP
  !
  ! History:
  ! Version   Date     Comment
  ! -------   ----     -------
  !   1.0    09/2002      Initial version     (F. Chevallier)
  !   1.1    05/2003      RTTOV7.3 compatible (F. Chevallier)
  !   1.2    03/2004      Added polarimetry   (R. Saunders)
  !   1.3    08/2004      Polarimetry fixes   (U. O'Keefe)
  !   1.4    11/2004      Clean-up            (P. Bauer)
  !   1.5    10/2005      Fixes for rttov8 indexing   (U. O'Keeffe)
  !   1.6    11/2005      Limit lines to 132 characters
  !                       add errorstatus to arguments
  !                       change stop to return (J Cameron)
  !   1.7    09/2006      Add if loop to stop use of iccmax index
  !                       if = 0 (A. Doherty)
  !   1.8    11/2007      RTTOV9 / cleanup (A. Geer)   
  !   1.9    03/2008      Revised cloud partitioning (A. Geer)  
  !   1.10   03/2009      Safety check on cloud fraction (A. Geer) 
  !   1.11   11/2009      RTTOV transmittances / optical depths come on levels (A Geer)
  !   1.12   11/2012      Remove "old cloud" option; new subroutine for hydrometeors (A Geer)
  !
  ! Code Description:
  !   Language:           Fortran 90.
  !   Software Standards: "European Standards for Writing and
  !   Documenting Exchangeable Fortran 90 Code".
  !
  ! Declarations:
  ! Modules used:
  ! Imported Type Definitions:

  Use rttov_types, Only :     &
       & rttov_coef          ,&
       & rttov_scatt_coef    ,&
       & transmission_type   ,&
       & geometry_Type       ,&
       & profile_scatt_aux   ,&
       & profile_Type        ,&
       & profile_cloud_Type  ,&
       & rttov_chanprof      ,&
       & rttov_options

  Use parkind1, Only : jpim, jplm
!INTF_OFF
  Use rttov_types, Only :     &
       & transmission_type_aux

  Use rttov_const, Only:      &
       & errorstatus_success, &
       & errorstatus_fatal,   &
       & gravity,             &
       & pressure_top,        &
       & rgp,                 &
       & rm,                  &
       & adk_adjoint,         &
       & adk_k

  Use parkind1, Only : jprb, jplm

  USE YOMHOOK, ONLY: LHOOK , DR_HOOK
!INTF_ON

  Implicit None

!* Subroutine arguments:
  Type(rttov_options), Intent (in) :: opts       ! RTTOV options
  Integer (Kind=jpim), Intent (in) :: nlevels ! Number of levels
  Integer (Kind=jpim), Intent (in) :: nprofiles  ! Number of profiles
  Integer (Kind=jpim), Intent (in) :: nprofilesad  ! Number of profiles in adjoint
  Integer (Kind=jpim), Intent (in) :: nchannels  ! Number of channels*profiles=radiances
  Integer (Kind=jpim), Intent(out) :: errorstatus             ! Error return code
  Type(rttov_chanprof), Intent(in) :: chanprof    (nchannels) ! Channel and profile indices
  Integer (Kind=jpim), Intent (in) :: frequencies (nchannels) ! Frequency indices
  Logical (Kind=jplm), Intent (in) :: calcemiss   (nchannels) ! Emissivity flags
  Logical (Kind=jplm), Intent (in) :: usercfrac               ! User has supplied cloud fraction

  Type (profile_Type),        Intent (in)    :: profiles        (nprofiles)   ! Atmospheric profiles
  Type (profile_Type),        Intent (inout) :: profiles_ad     (nprofilesad) ! Atmospheric profiles
  Type (rttov_coef),          Intent (in)    :: coef_rttov                    ! RTTOV Coefficients
  Type (rttov_scatt_coef),    Intent (in)    :: coef_scatt                    ! RTTOV_SCATT Coefficients
  Type (profile_cloud_Type),  Intent (in)    :: cld_profiles    (nprofiles)   ! Cloud profiles
  Type (profile_cloud_Type),  Intent (inout) :: cld_profiles_ad (nprofilesad) ! Cloud profiles
  Type (transmission_Type),   Intent (in)    :: transmission                  ! Transmittances and optical depths
  Type (transmission_Type),   Intent (inout) :: transmission_ad               ! Transmittances and optical depths
  Type (geometry_Type),       Intent (out)   :: angles          (nprofiles)   ! Zenith angles
  Type (profile_scatt_aux),   Intent (inout) :: scatt_aux                     ! Auxiliary profile variables
  Type (profile_scatt_aux),   Intent (inout) :: scatt_aux_ad                  ! Auxiliary profile variables

!INTF_END
 
!* Local variables
  Integer (Kind=jpim) :: ilayer, iprof, ichan
  Integer (Kind=jpim) :: iprofad, adk
  Real    (Kind=jprb) :: p1, p2, pm, p1_ad, p2_ad, pm_ad, dp2dz

  Real    (Kind=jprb), Dimension (nprofiles,nlevels)   :: presf   ! Pressure levels [hPa]
  Real    (Kind=jprb), Dimension (nprofiles,nlevels+1) :: presfh  ! Half-level pressure levels [hPa]
  Real    (Kind=jprb), Dimension (nchannels,nlevels)   :: od
  Real    (Kind=jprb), Dimension (nprofiles,nlevels)   :: dztop  ! Thickness of top half of level [km]
  Real    (Kind=jprb), Dimension (nprofiles,nlevels)   :: dzbot  ! Thickness of bottom half of level [km]
  Real    (Kind=jprb), Dimension (2:nlevels)           :: dzr    ! Thickness of RTTOV level [km]
  Real    (Kind=jprb), Dimension (nchannels,nlevels+1) :: od_rttov ! Single layer optical depths on RTTOV levels

  Real    (Kind=jprb), Dimension (nchannels)           :: zod_up_cld   ! Optical depth from top of the atmosphere 
  Real    (Kind=jprb), Dimension (nprofilesad,nlevels)   :: presf_ad     ! Pressure levels [hPa]
  Real    (Kind=jprb), Dimension (nprofilesad,nlevels+1) :: presfh_ad    ! Half-level pressure levels [hPa]
  Real    (Kind=jprb), Dimension (nchannels,nlevels)   :: od_ad
  Real    (Kind=jprb), Dimension (nprofilesad,nlevels)   :: dztop_ad 
  Real    (Kind=jprb), Dimension (nprofilesad,nlevels)   :: dzbot_ad  
  Real    (Kind=jprb), Dimension (2:nlevels)           :: dzr_ad    
  Real    (Kind=jprb), Dimension (nchannels,nlevels+1) :: od_rttov_ad 
  Real    (Kind=jprb), Dimension (nchannels)           :: zod_up_cld_ad    ! Optical depth from top of the atmosphere 
   
  Real    (Kind=jprb), Dimension (nchannels,nlevels)   :: ext_0
  Real    (Kind=jprb), Dimension (nchannels,nlevels)   :: ext_1, ssa_1, asm_1
  Real    (Kind=jprb), Dimension (nchannels,nlevels)   :: ext_2, ssa_2, asm_2
  Real    (Kind=jprb), Dimension (nchannels,nlevels)   :: ext_3, ssa_3, asm_3
 
  Type (transmission_Type_aux) :: transmissioncld     ! Clear+cloud transmittances with cloud
  Type (transmission_Type_aux) :: transmissioncld_ad  ! Clear+cloud transmittances with cloud

  Character (len=80) :: errMessage
  Character (len=18) :: NameOfRoutine = 'rttov_iniscatt_ad '

  REAL(KIND=JPRB) :: ZHOOK_HANDLE

#include "rttov_mieproc.interface"
#include "rttov_iniedd.interface"
#include "rttov_hydro.interface"
#include "rttov_calcemis_mw.interface"
#include "rttov_mieproc_ad.interface"
#include "rttov_iniedd_ad.interface"
#include "rttov_hydro_ad.interface"
#include "rttov_calcemis_mw_ad.interface"
#include "rttov_setgeometry.interface"
#include "rttov_errorreport.interface"
#include "rttov_calcemis_mw_k.interface"

  !- End of header --------------------------------------------------------

  if (lhook) call dr_hook('RTTOV_INISCATT_AD',0_jpim,zhook_handle)

  if (nprofilesad == nprofiles) then 
    adk = adk_adjoint   ! Adjoint mode
  else if (nprofilesad == nchannels) then
    adk = adk_k         ! K mode
  endif 

  ! NB this only truly needed if using FASTEM 3 
  allocate (transmissioncld % thermal_path1)
  allocate (transmissioncld % thermal_path1 % tau_surf (0:0,nchannels))
  allocate (transmissioncld_ad % thermal_path1)
  allocate (transmissioncld_ad % thermal_path1 % tau_surf (0:0,nchannels))

  errorstatus = errorstatus_success

  dp2dz = -1.0E-03_JPRB * rgp / gravity / rm 

  scatt_aux % ext (:,:) = 0.0_JPRB
  scatt_aux % ssa (:,:) = 0.0_JPRB
  scatt_aux % asm (:,:) = 0.0_JPRB

!* Security on user-defined pressures
  Do iprof = 1, nprofiles
     Do ilayer = 1, nlevels
        If (profiles (iprof) % p (ilayer) >= pressure_top) Then
            presf (iprof,ilayer) = profiles (iprof) % p (ilayer)
       else
            presf (iprof,ilayer) = pressure_top
       Endif
     Enddo
     Do ilayer = 1, nlevels + 1
        If (cld_profiles (iprof) % ph (ilayer) >= pressure_top ) Then
            presfh (iprof,ilayer) = cld_profiles (iprof) % ph (ilayer)
        else
            presfh (iprof,ilayer) = pressure_top
        Endif
     Enddo
  Enddo

!* Set up geometric variables
  Call rttov_setgeometry ( &
    & opts,            & ! in
    & .FALSE._jplm,    & ! in (inc. solar radiation)
    & profiles,        & ! in
    & coef=coef_rttov, & ! in
    & angles=angles)     ! out

!* Compute temperature at layer boundaries (K)
  Do iprof = 1, nprofiles
     scatt_aux % tbd (iprof,nlevels+1) = profiles (iprof) % s2m % t
     scatt_aux % tbd (iprof,1)         = profiles (iprof) % t(1)
  Enddo

  Do ilayer = 1, nlevels-1
     Do iprof = 1, nprofiles     
        p1 = presf  (iprof,ilayer+1)
        p2 = presf  (iprof,ilayer  )
        pm = presfh (iprof,ilayer+1)

        scatt_aux % tbd (iprof,ilayer+1) =  profiles (iprof) % t (ilayer+1) &
                                       & + (profiles (iprof) % t (ilayer)   & 
                                       & -  profiles (iprof) % t (ilayer+1)) &
                                       & / log(p2/p1) * log(pm/p1)        
     Enddo
  Enddo

!* Nadir heights (km)
  Do ilayer = nlevels, 1, -1
     Do iprof = 1, nprofiles
        p1 = presfh (iprof,ilayer+1)
        p2 = presfh (iprof,ilayer  )
        pm = presf  (iprof,ilayer  )

        If (p1 <= p2) then
           errorstatus = errorstatus_fatal
           Write( errMessage, '( "iniscatt : problem with user-defined pressure layering")' )
           Call rttov_errorreport (errorstatus, errMessage, NameOfRoutine)
           if (lhook) call dr_hook('RTTOV_INISCATT',1_jpim,zhook_handle)
           Return
        End If

        scatt_aux % dz (iprof,ilayer) = dp2dz * Log(p2/p1) * profiles (iprof) % t (ilayer)

        dzbot (iprof,ilayer) = dp2dz * Log(pm/p1) * profiles (iprof) % t (ilayer)
        dztop (iprof,ilayer) = dp2dz * Log(p2/pm) * profiles (iprof) % t (ilayer)

     Enddo
  Enddo

!* Get single-layer optical depths (at nadir and in hPa-1) and put onto model half levels
  Do ichan = 1, nchannels
     iprof = chanprof(ichan) % prof

! Top RTTOV level to space    
     od_rttov (ichan,1)         = -1.0_jprb * log( transmission % tau_levels (1,ichan) )     
     Do ilayer = 2, nlevels 
        od_rttov (ichan,ilayer) = log( transmission % tau_levels (ilayer-1,ichan) ) &
                              & - log( transmission % tau_levels (ilayer,ichan) )
     Enddo
! Surface to bottom RTTOV (full pressure) level
     od_rttov (ichan,nlevels+1) = log( transmission % tau_levels (nlevels,ichan) ) &
                              & - log( transmission % tau_total (ichan) )

! RTTOV layer thickness
     Do ilayer = 2, nlevels  
       dzr(ilayer) = dzbot(iprof,ilayer-1) + dztop(iprof,ilayer)
     Enddo

! Re-allocate optical depths between half pressure levels        
     od (ichan,1)         = od_rttov(ichan,1) &
                        & + od_rttov(ichan,2) * dzbot(iprof,1) / dzr(2)     
     Do ilayer = 2, nlevels - 1  
        od (ichan,ilayer) = od_rttov(ichan,ilayer)   * dztop(iprof,ilayer) / dzr(ilayer) &
                        & + od_rttov(ichan,ilayer+1) * dzbot(iprof,ilayer) / dzr(ilayer+1) 
     Enddo
     od (ichan,nlevels)   = od_rttov(ichan,nlevels)  * dztop(iprof,nlevels) / dzr(nlevels) &
                        & + od_rttov(ichan,nlevels+1) 
      
  Enddo

  !* Optical depths in km-1 and at nadir
  Do ilayer = 1, nlevels
     Do ichan = 1, nchannels
        iprof = chanprof(ichan) % prof
     
        scatt_aux % ext (ichan,ilayer) = od (ichan,ilayer) &
          & / scatt_aux % dz (iprof,ilayer) * angles (iprof) % coszen 
     
        ext_0 (ichan,ilayer) = scatt_aux % ext (ichan,ilayer)     

        if (scatt_aux % ext (ichan,ilayer) < 1.0E-10_JPRB) scatt_aux % ext (ichan,ilayer) = 1.0E-10_JPRB
     Enddo

 Enddo

!* Store clear-sky absorption/scattering parameters
  ext_1 (:,:) = scatt_aux % ext (:,:)
  ssa_1 (:,:) = scatt_aux % ssa (:,:)
  asm_1 (:,:) = scatt_aux % asm (:,:)

!* Creation of the hydrometeor profiles
  Call rttov_hydro (        &
       & nlevels,           &! in
       & nprofiles,         &! in   
       & usercfrac,         &! in
       & presf,             &! in 
       & profiles,          &! in
       & cld_profiles,      &! in
       & coef_scatt,        &! in
       & scatt_aux)          ! inout 

!* Cloud/rain absorption/scattering parameters
  Call rttov_mieproc (      &
       & nlevels,           &! in
       & nchannels,         &! in
       & nprofiles,         &! in
       & frequencies,       &! in
       & chanprof%prof,     &! in
       & profiles,          &! in
       & coef_scatt,        &! in
       & scatt_aux)          ! inout 
       
!* Store clear+cloud+rain absorption/scattering parameters
  ext_2 (:,:) = scatt_aux % ext (:,:)
  ssa_2 (:,:) = scatt_aux % ssa (:,:)
  asm_2 (:,:) = scatt_aux % asm (:,:)

!* Scattering parameters for Eddington RT
  Call rttov_iniedd(        &
       & nlevels,           &! in
       & nchannels ,        &! in
       & nprofiles ,        &! in
       & chanprof%prof,     &! in
       & angles    ,        &! in
       & scatt_aux)          ! inout 

!* Store delta-scaled clear+cloud+rain absorption/scattering parameters
  ext_3 (:,:) = scatt_aux % ext (:,:)
  ssa_3 (:,:) = scatt_aux % ssa (:,:)
  asm_3 (:,:) = scatt_aux % asm (:,:)

!* Surface emissivities
  zod_up_cld (:) = 0.0_JPRB
  
  Do ichan = 1, nchannels
     iprof = chanprof(ichan) % prof
     
     Do ilayer = 1, nlevels     
        zod_up_cld (ichan) = zod_up_cld (ichan) + scatt_aux % ext (ichan,ilayer) * scatt_aux % dz (iprof,ilayer) 
     Enddo
     if (zod_up_cld (ichan) >= 30.0_JPRB) zod_up_cld (ichan) = 30.0_JPRB
     transmissioncld % thermal_path1 % tau_surf (0,ichan) = Exp(-1.0_JPRB * zod_up_cld (ichan) / angles (iprof) % coszen)
  Enddo
  
  Call rttov_calcemis_mw(      &
       & opts,                 &! in
       & profiles,             &! in
       & angles,               &! in
       & coef_rttov,           &! in
       & chanprof,            &! in
       & transmissioncld,      &! in
       & calcemiss,            &! in
       & scatt_aux % ems_cld,  &! inout
       & scatt_aux % ref_cld,  &! out
       & errorstatus          ) ! inout 

!* Hemispheric emissivity (= Fastem's effective emissivity)
  scatt_aux % ems_bnd (:) = scatt_aux % ems_cld (:)
  scatt_aux % ref_bnd (:) = scatt_aux % ref_cld (:)

!* ADJOINT PART
!* Hemispheric emissivity (= Fastem's effective emissivity)
  scatt_aux_ad % ems_cld (:) = scatt_aux_ad % ems_cld (:) + scatt_aux_ad % ems_bnd (:)
  scatt_aux_ad % ems_bnd (:) = 0.0_JPRB
  
  scatt_aux_ad % ref_cld (:) = scatt_aux_ad % ref_cld (:) + scatt_aux_ad % ref_bnd (:)  
  scatt_aux_ad % ref_bnd (:) = 0.0_JPRB
  
  transmissioncld_ad % thermal_path1 % tau_surf (0,:) = 0.0_JPRB

  if (adk == adk_adjoint) then   
    Call rttov_calcemis_mw_ad(           &
            & opts,                      &! in
            & profiles,                  &! in
            & profiles_ad,               &! inout
            & angles,                    &! in
            & coef_rttov,                &! in
            & chanprof,                 &! in
            & transmissioncld    ,       &! in
            & transmissioncld_ad,        &! in
            & calcemiss,                 &! in
            & scatt_aux_ad % ems_cld,    &! inout
            & scatt_aux_ad % ref_cld)     ! inout 
  else if (adk == adk_k) then 
    Call rttov_calcemis_mw_k(            &
            & opts,                      &! in
            & profiles,                  &! in
            & profiles_ad,               &! inout
            & angles,                    &! in
            & coef_rttov,                &! in
            & chanprof,                 &! in
            & transmissioncld    ,       &! in
            & transmissioncld_ad,        &! in
            & calcemiss,                 &! in
            & scatt_aux_ad % ems_cld,    &! inout
            & scatt_aux_ad % ref_cld)     ! inout 
  endif
 
  zod_up_cld_ad (:) = 0.0_JPRB
  
  Do ichan = 1, nchannels
     iprof = chanprof(ichan) % prof
     
     zod_up_cld_ad (ichan) = zod_up_cld_ad (ichan) - transmissioncld_ad % thermal_path1 % tau_surf (0,ichan) &
                         & * transmissioncld % thermal_path1 % tau_surf (0,ichan) / angles (iprof) % coszen
     transmissioncld_ad % thermal_path1 % tau_surf (0,ichan) = 0.0_JPRB

     if (zod_up_cld (ichan) == 30.0_JPRB) zod_up_cld_ad (ichan) = 0.0_JPRB

     Do ilayer = 1, nlevels
        iprof = chanprof(ichan) % prof
        if (adk == adk_adjoint) then
          iprofad = iprof  
        else if (adk == adk_k) then
          iprofad = ichan  
        endif
   
        scatt_aux_ad % ext (ichan,ilayer) = scatt_aux_ad % ext (ichan,ilayer) & 
                                        & + scatt_aux % dz  (iprof,ilayer) * zod_up_cld_ad (ichan)
        scatt_aux_ad % dz  (iprofad,ilayer) = scatt_aux_ad % dz  (iprofad,ilayer) & 
                                        & + scatt_aux % ext (ichan,ilayer) * zod_up_cld_ad (ichan)
     Enddo
  Enddo
  zod_up_cld_ad (:) = 0.0_JPRB

  scatt_aux % ext (:,:) = ext_2 (:,:) 
  scatt_aux % ssa (:,:) = ssa_2 (:,:) 
  scatt_aux % asm (:,:) = asm_2 (:,:) 

!* Scattering parameters for Eddington RT
  Call rttov_iniedd_ad(     &
       & nlevels,           &! in
       & nchannels ,        &! in
       & nprofiles ,        &! in
       & nprofilesad,       &! in
       & chanprof%prof,     &! in
       & angles    ,        &! in
       & scatt_aux ,        &! inout
       & scatt_aux_ad)       ! inout 
       
!* Cloud/rain absorption/scattering parameters
  scatt_aux % ext (:,:) = ext_1 (:,:) 
  scatt_aux % ssa (:,:) = ssa_1 (:,:) 
  scatt_aux % asm (:,:) = asm_1 (:,:) 

  Call rttov_mieproc_ad (&
       & nlevels,           &! in
       & nchannels,         &! in
       & nprofiles,         &! in
       & nprofilesad,       &! in
       & frequencies,       &! in
       & chanprof%prof,     &! in
       & profiles,          &! in
       & coef_scatt,        &! in
       & scatt_aux,         &! inout
       & scatt_aux_ad)       ! inout 

  presfh_ad (:,:) = 0.0_JPRB
  presf_ad  (:,:) = 0.0_JPRB

  !* Creation of the hydrometeor profiles
  Call rttov_hydro_ad (        &
       & nlevels,           &! in
       & nprofiles,         &! in   
       & nprofilesad,       &! in  
       & nchannels,         &! in  
       & usercfrac,         &! in
       & chanprof,          &! in
       & presf,             &! in 
       & presf_ad,          &! inout  
       & profiles,          &! in
       & profiles_ad,       &! inout 
       & cld_profiles,      &! in
       & cld_profiles_ad,   &! inout 
       & coef_scatt,        &! in
       & scatt_aux,         &! inout 
       & scatt_aux_ad)       ! inout 
       
  !* Optical depths in km-1 and at nadir
  Do ilayer = 1,nlevels
     Do ichan = 1, nchannels
        iprof = chanprof(ichan) % prof
        if (adk == adk_adjoint) then
          iprofad = iprof  
        else if (adk == adk_k) then
          iprofad = ichan  
        endif
     
        If (ext_0 (ichan,ilayer) < 1.0E-10_JPRB) scatt_aux_ad % ext (ichan,ilayer) = 0.0_JPRB

        od_ad (ichan,ilayer) = scatt_aux_ad % ext (ichan,ilayer) &
          & / scatt_aux % dz (iprof,ilayer) * angles (iprof) % coszen   
        scatt_aux_ad % dz (iprofad,ilayer) = scatt_aux_ad % dz (iprofad,ilayer) &
          & - od (ichan,ilayer) * angles (iprof) % coszen &
          & / (scatt_aux % dz (iprof,ilayer) * scatt_aux % dz (iprof,ilayer)) &
          & * scatt_aux_ad % ext (ichan,ilayer) 
        scatt_aux_ad % ext (ichan,ilayer) = 0.0_JPRB
     Enddo
  Enddo

  dzbot_ad(:,:)    = 0.0_JPRB
  dztop_ad(:,:)    = 0.0_JPRB
  od_rttov_ad(:,:) = 0.0_JPRB

  Do ichan = 1, nchannels

     iprof = chanprof(ichan) % prof
     if (adk == adk_adjoint) then
        iprofad = iprof  
     else if (adk == adk_k) then
        iprofad = ichan  
     endif

     dzr_ad(:)         = 0.0_JPRB

     Do ilayer = 2, nlevels  
       dzr(ilayer) = dzbot(iprof,ilayer-1) + dztop(iprof,ilayer)
     Enddo

! Re-allocate optical depths between half pressure levels 
     od_rttov_ad(ichan,1) = od_rttov_ad(ichan,1) + od_ad (ichan,1)  
     od_rttov_ad(ichan,2) = od_rttov_ad(ichan,2) + od_ad (ichan,1) &
                        & * dzbot(iprof,1) / dzr(2) 
     dzbot_ad(iprofad,1)  = dzbot_ad(iprofad,1)  + od_ad (ichan,1) &
                        & * od_rttov(ichan,2) / dzr(2) 
     dzr_ad(2)            = dzr_ad(2) - od_ad (ichan,1) * od_rttov(ichan,2) &
                        & * dzbot(iprof,1) / dzr(2)**2   
   
     Do ilayer = 2, nlevels - 1  
        od_rttov_ad(ichan,ilayer)   = od_rttov_ad(ichan,ilayer)                     & 
          & + od_ad (ichan,ilayer) * dztop(iprof,ilayer) / dzr(ilayer)
        dztop_ad(iprofad,ilayer)    = dztop_ad(iprofad,ilayer)                      &  
          & + od_ad (ichan,ilayer) * od_rttov(ichan,ilayer) / dzr(ilayer)
        dzr_ad(ilayer)              = dzr_ad(ilayer)                                &  
          & - od_ad (ichan,ilayer) * od_rttov(ichan,ilayer) * dztop(iprof,ilayer)   &
          & / dzr(ilayer)**2
        od_rttov_ad(ichan,ilayer+1) = od_rttov_ad(ichan,ilayer+1)                   &
          & + od_ad (ichan,ilayer) * dzbot(iprof,ilayer) / dzr(ilayer+1)
        dzbot_ad(iprofad,ilayer)    = dzbot_ad(iprofad,ilayer)                      &
          & + od_ad (ichan,ilayer) * od_rttov(ichan,ilayer+1) / dzr(ilayer+1)     
        dzr_ad(ilayer+1)            = dzr_ad(ilayer+1)                              &         
          & - od_ad (ichan,ilayer) * od_rttov(ichan,ilayer+1) * dzbot(iprof,ilayer) &
          & / dzr(ilayer+1)**2
     Enddo

     od_rttov_ad(ichan,nlevels)   = od_rttov_ad(ichan,nlevels)                    &
       & + od_ad (ichan,nlevels) * dztop(iprof,nlevels) / dzr(nlevels)
     dztop_ad(iprofad,nlevels)    = dztop_ad(iprofad,nlevels)                     &
       & + od_ad (ichan,nlevels) * od_rttov(ichan,nlevels) / dzr(nlevels)
     dzr_ad(nlevels)              = dzr_ad(nlevels)                               &             
       & - od_ad (ichan,nlevels) * od_rttov(ichan,nlevels) * dztop(iprof,nlevels) &
       & / dzr(nlevels)**2
     od_rttov_ad(ichan,nlevels+1) = od_rttov_ad(ichan,nlevels+1) + od_ad (ichan,nlevels)

     od_ad (ichan,:) = 0.0_JPRB

! RTTOV layer thickness
     Do ilayer = 2, nlevels  
       dzbot_ad(iprofad,ilayer-1) = dzbot_ad(iprofad,ilayer-1) + dzr_ad(ilayer)
       dztop_ad(iprofad,ilayer)   = dztop_ad(iprofad,ilayer)   + dzr_ad(ilayer)
       dzr_ad(ilayer)           = 0.0_JPRB
     Enddo

     transmission_ad % tau_levels (nlevels,ichan) =     &
          & + transmission_ad % tau_levels (nlevels,ichan) &
          & + od_rttov_ad(ichan,nlevels+1) / transmission % tau_levels (nlevels,ichan) 

     transmission_ad % tau_total (ichan) =     &
          & + transmission_ad % tau_total (ichan) &
          & - od_rttov_ad(ichan,nlevels+1) / transmission % tau_total (ichan) 
     
     Do ilayer = nlevels, 2, -1

        transmission_ad % tau_levels (ilayer-1,ichan) =     &
          & + transmission_ad % tau_levels (ilayer-1,ichan) &
          & + od_rttov_ad(ichan,ilayer) / transmission % tau_levels (ilayer-1,ichan) 

        transmission_ad % tau_levels (ilayer,ichan) =     &
          & + transmission_ad % tau_levels (ilayer,ichan) &
          & - od_rttov_ad(ichan,ilayer) / transmission % tau_levels (ilayer,ichan) 
       
     Enddo

     transmission_ad % tau_levels (1,ichan) = transmission_ad % tau_levels (1,ichan) &
                   & - od_rttov_ad (ichan,1) / transmission % tau_levels (1,ichan)

  Enddo

  od_rttov_ad (:,:) = 0.0_JPRB  
 
!* Nadir heights (km)
  Do ilayer = 1, nlevels
     Do iprofad = 1, nprofilesad

        if (adk == adk_adjoint) then
          iprof = iprofad  
        else if (adk == adk_k) then
          iprof = chanprof(iprofad) % prof  
        endif

        p1 = presfh (iprof,ilayer+1)
        p2 = presfh (iprof,ilayer  )
        pm = presf  (iprof,ilayer  )

        p1_ad = 0.0_JPRB
        p2_ad = 0.0_JPRB
        pm_ad = 0.0_JPRB

        p2_ad = p2_ad + dp2dz / p2 * profiles (iprof) % t (ilayer) * dztop_ad (iprofad,ilayer)
        pm_ad = pm_ad - dp2dz / pm * profiles (iprof) % t (ilayer) * dztop_ad (iprofad,ilayer)
        profiles_ad (iprofad) % t (ilayer) = profiles_ad (iprofad) % t (ilayer) & 
                                       & + dp2dz * Log(p2/pm) * dztop_ad (iprofad,ilayer) 
        dztop_ad (iprofad,ilayer) = 0.0_JPRB

        pm_ad = pm_ad + dp2dz / pm * profiles (iprof) % t (ilayer) * dzbot_ad (iprofad,ilayer)
        p1_ad = p1_ad - dp2dz / p1 * profiles (iprof) % t (ilayer) * dzbot_ad (iprofad,ilayer)
        profiles_ad (iprofad) % t (ilayer) = profiles_ad (iprofad) % t (ilayer) & 
                                       & + dp2dz * Log(pm/p1) * dzbot_ad (iprofad,ilayer) 
        dzbot_ad (iprofad,ilayer) = 0.0_JPRB

        p2_ad = p2_ad + dp2dz / p2 * profiles (iprof) % t (ilayer) * scatt_aux_ad % dz (iprofad,ilayer)
        p1_ad = p1_ad - dp2dz / p1 * profiles (iprof) % t (ilayer) * scatt_aux_ad % dz (iprofad,ilayer)
        profiles_ad (iprofad) % t (ilayer) = profiles_ad (iprofad) % t (ilayer) & 
                                       & + dp2dz * Log(p2/p1) * scatt_aux_ad % dz (iprofad,ilayer) 
        scatt_aux_ad % dz (iprofad,ilayer) = 0.0_JPRB

        presf_ad  (iprofad,ilayer)   = presf_ad  (iprofad,ilayer)   + pm_ad
        presfh_ad (iprofad,ilayer)   = presfh_ad (iprofad,ilayer)   + p2_ad
        presfh_ad (iprofad,ilayer+1) = presfh_ad (iprofad,ilayer+1) + p1_ad
     Enddo
  Enddo

!* Temperature at layer boundaries (K)
  Do ilayer = nlevels - 1, 1, -1
     Do iprofad = 1, nprofilesad

        if (adk == adk_adjoint) then
           iprof = iprofad  
        else if (adk == adk_k) then
           iprof = chanprof(iprofad) % prof  
        endif

        p1 = presf  (iprof,ilayer+1)
        p2 = presf  (iprof,ilayer  )
        pm = presfh (iprof,ilayer+1)

        profiles_ad (iprofad) % t (ilayer+1) = profiles_ad (iprofad) % t (ilayer+1) &
         & + scatt_aux_ad % tbd (iprofad,ilayer+1)                                  &
         & - 1.0_JPRB / Log(p2/p1) * Log(pm/p1) * scatt_aux_ad % tbd (iprofad,ilayer+1) 
        profiles_ad (iprofad) % t (ilayer)   = profiles_ad (iprofad) % t (ilayer)   &
         & + 1.0_JPRB / Log(p2/p1) * Log(pm/p1) * scatt_aux_ad % tbd (iprofad,ilayer+1) 
     
        p1_ad = (profiles (iprof) % t (ilayer) - profiles (iprof) % t (ilayer+1)) &
            & / (Log(p2/p1) * Log(p2/p1)) / p1 * Log(pm/p1) * scatt_aux_ad % tbd (iprofad,ilayer+1) &
            & - (profiles (iprof) % t (ilayer) - profiles (iprof) % t (ilayer+1)) &
            & /  Log(p2/p1) / p1 * scatt_aux_ad % tbd (iprofad,ilayer+1)     
        p2_ad = -1.0_JPRB &
            & * (profiles (iprof) % t (ilayer) - profiles (iprof) % t (ilayer+1)) &
            & *  Log(pm/p1) / (Log(p2/p1) * Log(p2/p1)) / p2 * scatt_aux_ad % tbd (iprofad,ilayer+1) 
        pm_ad = (profiles (iprof) % t (ilayer) - profiles (iprof) % t (ilayer+1)) &
            & /  Log(p2/p1) / pm * scatt_aux_ad % tbd (iprofad,ilayer+1) 
        scatt_aux_ad % tbd (iprofad,ilayer+1) = 0.0_JPRB

        presf_ad  (iprofad,ilayer+1) = presf_ad  (iprofad,ilayer+1) + p1_ad
        presf_ad  (iprofad,ilayer  ) = presf_ad  (iprofad,ilayer  ) + p2_ad
        presfh_ad (iprofad,ilayer+1) = presfh_ad (iprofad,ilayer+1) + pm_ad
     Enddo
  Enddo

  Do iprofad = 1, nprofilesad
     profiles_ad (iprofad) % s2m % t = profiles_ad (iprofad) % s2m % t + scatt_aux_ad % tbd (iprofad,nlevels+1) 
     profiles_ad (iprofad) % t (1)   = profiles_ad (iprofad) % t (1)   + scatt_aux_ad % tbd (iprofad,1) 
  Enddo
  scatt_aux_ad % tbd (:,:) = 0.0_JPRB

!* Security on user-defined pressures
  Do iprofad = 1, nprofilesad

     if (adk == adk_adjoint) then
        iprof = iprofad  
     else if (adk == adk_k) then
        iprof = chanprof(iprofad) % prof 
     endif

     Do ilayer = 1, nlevels
        If (profiles    (iprof) % p (ilayer) >= pressure_top) &
         &  profiles_ad (iprofad) % p (ilayer) = profiles_ad (iprofad) % p (ilayer) + presf_ad (iprofad,ilayer) 
        presf_ad  (iprofad,ilayer) = 0.0_JPRB
     Enddo
     Do ilayer = 1, nlevels + 1
        If (cld_profiles    (iprof) % ph (ilayer) >= pressure_top) &
         &  cld_profiles_ad (iprofad) % ph (ilayer) = cld_profiles_ad (iprofad) % ph (ilayer) + presfh_ad (iprofad,ilayer) 
        presfh_ad (iprofad,ilayer) = 0.0_JPRB
     Enddo
  Enddo

  scatt_aux % ext (:,:) = ext_3 (:,:) 
  scatt_aux % ssa (:,:) = ssa_3 (:,:) 
  scatt_aux % asm (:,:) = asm_3 (:,:) 

!* Deallocate
  deallocate (transmissioncld    % thermal_path1 % tau_surf)
  deallocate (transmissioncld    % thermal_path1)
  deallocate (transmissioncld_ad % thermal_path1 % tau_surf)
  deallocate (transmissioncld_ad % thermal_path1)

  if (lhook) call dr_hook('RTTOV_INISCATT_AD',1_jpim,zhook_handle)
 
End Subroutine rttov_iniscatt_ad
