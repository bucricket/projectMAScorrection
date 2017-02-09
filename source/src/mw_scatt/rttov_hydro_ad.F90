!
Subroutine rttov_hydro_ad (&
      & nlevels,           &! in
      & nprofiles,         &! in
      & nprofilesad,       &! in
      & nchannels,         &! in
      & usercfrac,         &! in
      & chanprof,          &! in
      & presf,             &! in
      & presfad,           &! inout
      & profiles,          &! in  
      & profiles_ad,       &! inout  
      & cld_profiles,      &! in 
      & cld_profiles_ad,   &! inout 
      & coef_scatt,        &! in
      & scatt_aux,         &! inout 
      & scatt_aux_ad)       ! inout 

  !
  ! Description:
  ! Initialises the hydrometeor profile in the cloudy sub-column.
  !
  ! Copyright:
  !    This software was developed within the context of
  !    the EUMETSAT Satellite Application Facility on
  !    Numerical Weather Prediction (NWP SAF), under the
  !    Cooperation Agreement dated 25 November 1998, between
  !    EUMETSAT and the Met Office, UK, by one or more painiscattrtners
  !    within the NWP SAF. The partners in the NWP SAF are
  !    the Met Office, ECMWF, KNMI and MeteoFrance.
  !
  !    Copyright 2012, EUMETSAT, All Rights Reserved.
  !
  ! Method:
  ! 
  ! Current Code Owner: SAF NWP
  !
  ! History:
  ! Version   Date     Comment
  ! -------   ----     -------
  !   1.0    09/2012   Initial version     (A Geer)
  !
  ! Code Description:
  !   Language:           Fortran 90.
  !   Software Standards: "European Standards for Writing and
  !     Documenting Exchangeable Fortran 90 Code".
  ! 
  ! Declarations:
  ! Modules used:
  ! Imported Type Definitions:
  Use rttov_types, Only :     &
       & rttov_scatt_coef    ,&
       & profile_scatt_aux   ,&
       & profile_Type        ,&
       & profile_cloud_Type  ,&
       & rttov_chanprof 

  Use parkind1, Only : jpim, jplm, jprb

!INTF_OFF
  Use rttov_const, Only:      &
       & rgp,                 &
       & rm,                  &
       & rho_rain,            &
       & rho_snow,            &
       & ccthres,             &
       & adk_adjoint,         &
       & adk_k

  USE YOMHOOK, ONLY: LHOOK , DR_HOOK
!INTF_ON
  Implicit None
  
!* Subroutine arguments:
  Integer (Kind=jpim), Intent (in) :: nlevels   ! Number of levels
  Integer (Kind=jpim), Intent (in) :: nprofiles ! Number of profiles
  Integer (Kind=jpim), Intent (in) :: nprofilesad  ! Number of profiles in adjoint
  Integer (Kind=jpim), Intent (in) :: nchannels  ! Number of channels
  Logical (Kind=jplm), Intent (in) :: usercfrac               ! User has supplied cloud fraction
  Type(rttov_chanprof), Intent(in) :: chanprof    (nchannels) ! Channel and profile indices
  Real (Kind=jprb),           Intent (in)    :: presf(nprofiles,nlevels) ! Pressure levels [hPa]
  Real (Kind=jprb),           Intent (inout) :: presfad(nprofilesad,nlevels) ! Pressure levels [hPa]
  Type (profile_Type),        Intent (in)    :: profiles (nprofiles)     ! Atmospheric profiles
  Type (profile_Type),        Intent (inout) :: profiles_ad (nprofilesad)     ! Atmospheric profiles
  Type (rttov_scatt_coef),    Intent (in)    :: coef_scatt               ! RTTOV_SCATT Coefficients
  Type (profile_cloud_Type),  Intent (in)    :: cld_profiles (nprofiles) ! Cloud profiles 
  Type (profile_cloud_Type),  Intent (inout) :: cld_profiles_ad (nprofilesad) ! Cloud profiles AD
  Type (profile_scatt_aux),   Intent (inout) :: scatt_aux                ! Auxiliary profile variables
  Type (profile_scatt_aux),   Intent (inout) :: scatt_aux_ad             ! Auxiliary profile variables
!INTF_END

!* Local variables
  Integer (Kind=jpim) :: ilayer, iprof
  Integer (Kind=jpim) :: iprofad, adk
  Real    (Kind=jprb) :: de2mr, mmr_to_density(nprofiles,nlevels), mmr_to_density_ad
  Real    (Kind=jprb) :: cfrac(nprofiles)
  Real    (Kind=jprb) :: hydro_weights(nprofiles,nlevels), hydro_column(nprofiles) ! Total hydrometeor amounts [g m^-2]
  Real    (Kind=jprb) :: hydro_weights_ad(nlevels), hydro_column_ad  
  Real    (Kind=jprb), Dimension (nprofiles,nlevels)   :: clw_scale, ciw_scale, &
   & totalice_scale, rain_scale, sp_scale
  Real    (Kind=jprb), Dimension (nprofiles,nlevels)   :: clw_precf, ciw_precf, &
   & rain_precf, sp_precf, totalice_precf

  Character (len=15) :: NameOfRoutine = 'rttov_hydro_ad '
  
  REAL(KIND=JPRB) :: ZHOOK_HANDLE
  
  !- End of header --------------------------------------------------------

  if (lhook) call dr_hook('RTTOV_HYDRO_AD',0_jpim,zhook_handle)

  if (nprofilesad == nprofiles) then 
    adk = adk_adjoint   ! Adjoint mode
  else if (nprofilesad == nchannels) then
    adk = adk_k         ! K mode
  endif 

  de2mr =  1.0E+05_JPRB * rm / rgp

  !* Initialise cloud and rain properties of the cloudy/rainy column
  do ilayer=1,nlevels
    do iprof = 1, nprofiles  
      scatt_aux % clw  (iprof,ilayer) = cld_profiles (iprof) % clw  (ilayer) 
      scatt_aux % rain (iprof,ilayer) = cld_profiles (iprof) % rain (ilayer) 
      if ( cld_profiles (iprof) % use_totalice ) then
        scatt_aux % totalice (iprof,ilayer) = cld_profiles (iprof) % totalice (ilayer) 
        scatt_aux % ciw      (iprof,ilayer) = 0.0_JPRB
        scatt_aux % sp       (iprof,ilayer) = 0.0_JPRB
      else
        scatt_aux % totalice (iprof,ilayer) = 0.0_JPRB
        scatt_aux % ciw      (iprof,ilayer) = cld_profiles (iprof) % ciw  (ilayer) 
        scatt_aux % sp       (iprof,ilayer) = cld_profiles (iprof) % sp   (ilayer) 
      endif 
    enddo
  enddo

  !* Save for adjoint
  clw_scale     (:,:) = scatt_aux % clw      (:,:)     
  ciw_scale     (:,:) = scatt_aux % ciw      (:,:)
  totalice_scale(:,:) = scatt_aux % totalice (:,:) 

  !* Change units
  do ilayer = 1, nlevels
     do iprof = 1, nprofiles

        !* Condensate from g/g to g/m^3
        mmr_to_density(iprof,ilayer) = presf (iprof,ilayer) * de2mr / profiles (iprof) % t (ilayer) 
        scatt_aux % clw      (iprof,ilayer) = scatt_aux % clw      (iprof,ilayer) * mmr_to_density(iprof,ilayer)
        scatt_aux % ciw      (iprof,ilayer) = scatt_aux % ciw      (iprof,ilayer) * mmr_to_density(iprof,ilayer)
        scatt_aux % totalice (iprof,ilayer) = scatt_aux % totalice (iprof,ilayer) * mmr_to_density(iprof,ilayer)

        if( cld_profiles (iprof) % mmr_snowrain ) then

          !* Save for adjoint
          rain_scale (iprof,ilayer) = scatt_aux % rain (iprof,ilayer) 
          sp_scale   (iprof,ilayer) = scatt_aux % sp   (iprof,ilayer) 

          scatt_aux % rain (iprof,ilayer) = scatt_aux % rain (iprof,ilayer) * mmr_to_density(iprof,ilayer)  
          scatt_aux % sp   (iprof,ilayer) = scatt_aux % sp   (iprof,ilayer) * mmr_to_density(iprof,ilayer)  

        else

          !* Rates from kg/m^2/s to g/m^3
          scatt_aux % rain (iprof,ilayer) = scatt_aux % rain (iprof,ilayer) / rho_rain
          scatt_aux % sp   (iprof,ilayer) = scatt_aux % sp   (iprof,ilayer) / rho_snow

          scatt_aux % rain (iprof,ilayer) = scatt_aux % rain (iprof,ilayer) * 3600.0_JPRB 
          scatt_aux % sp   (iprof,ilayer) = scatt_aux % sp   (iprof,ilayer) * 3600.0_JPRB

          !* Save for adjoint
          rain_scale (iprof,ilayer) = scatt_aux % rain (iprof,ilayer) 
          sp_scale   (iprof,ilayer) = scatt_aux % sp   (iprof,ilayer) 
 
          if (scatt_aux % rain (iprof,ilayer) > 0.0_JPRB) &
            & scatt_aux % rain (iprof,ilayer) = (scatt_aux % rain (iprof,ilayer) * &
            & coef_scatt % conv_rain (1))**(coef_scatt % conv_rain (2)) 
          if (scatt_aux % sp   (iprof,ilayer) > 0.0_JPRB) &
            & scatt_aux % sp   (iprof,ilayer) = (scatt_aux % sp   (iprof,ilayer) * &
            & coef_scatt % conv_sp   (1))**(coef_scatt % conv_sp   (2))

        endif
     enddo
  enddo

  !* Save for adjoint
  clw_precf     (:,:) = scatt_aux % clw      (:,:)     
  ciw_precf     (:,:) = scatt_aux % ciw      (:,:)
  rain_precf    (:,:) = scatt_aux % rain     (:,:) 
  sp_precf      (:,:) = scatt_aux % sp       (:,:) 
  totalice_precf(:,:) = scatt_aux % totalice (:,:) 

  scatt_aux % cfrac (:)   = 0.0_JPRB

  do iprof = 1, nprofiles

    if( usercfrac ) then

      !* User-supplied cloud fraction
      scatt_aux % cfrac (iprof) = cld_profiles (iprof) % cfrac

    else

      !* Calculate a hydrometeor-weighted average cloudy sky fraction 

      !* Partial column of hydrometeors in g/m^2
      hydro_weights(iprof,:) = &
        & ( scatt_aux % rain (iprof,:) + scatt_aux % sp (iprof,:) &
        &   + scatt_aux % ciw (iprof,:)  + scatt_aux % clw (iprof,:)  &
        &   + scatt_aux % totalice (iprof,:) ) &
        & * scatt_aux % dz (iprof,:)
      hydro_column(iprof) = sum( hydro_weights(iprof,:) )

      !* Weighted mean cloud fraction
      if (hydro_column(iprof) > 1e-10_JPRB) then
        scatt_aux % cfrac (iprof) = &
          & sum(hydro_weights(iprof,:) * cld_profiles (iprof) % cc (:)) / hydro_column(iprof)
        cfrac(iprof) = scatt_aux % cfrac (iprof)
        if ( scatt_aux % cfrac (iprof) < 0.0_JPRB ) scatt_aux % cfrac (iprof) = 0.0_JPRB
        if ( scatt_aux % cfrac (iprof) > 1.0_JPRB ) scatt_aux % cfrac (iprof) = 1.0_JPRB
      else
        scatt_aux % cfrac (iprof) = 0.0_JPRB
      endif

    endif

    !* Partition all cloud and rain into the cloudy column
    if (scatt_aux % cfrac (iprof) > ccthres) Then
        scatt_aux % clw  (iprof,:) = scatt_aux % clw  (iprof,:) / scatt_aux % cfrac (iprof)
        scatt_aux % ciw  (iprof,:) = scatt_aux % ciw  (iprof,:) / scatt_aux % cfrac (iprof)
        scatt_aux % totalice  (iprof,:) = scatt_aux % totalice  (iprof,:) / scatt_aux % cfrac (iprof)
        scatt_aux % rain (iprof,:) = scatt_aux % rain (iprof,:) / scatt_aux % cfrac (iprof)
        scatt_aux % sp   (iprof,:) = scatt_aux % sp   (iprof,:) / scatt_aux % cfrac (iprof)
    else
        scatt_aux % clw  (iprof,:) = 0.0_JPRB
        scatt_aux % ciw  (iprof,:) = 0.0_JPRB
        scatt_aux % totalice  (iprof,:) = 0.0_JPRB
        scatt_aux % rain (iprof,:) = 0.0_JPRB
        scatt_aux % sp   (iprof,:) = 0.0_JPRB
    endif
  enddo

  ! ADJOINT PART

  ! Calculate a hydrometeor-weighted average cloudy sky fraction 
  Do iprofad = 1, nprofilesad

    if (adk == adk_adjoint) then
      iprof = iprofad
    else if (adk == adk_k) then
      iprof = chanprof(iprofad) % prof
    endif

    !* Partition all cloud and rain into the cloudy column
    If (scatt_aux % cfrac (iprof) > ccthres) Then
      scatt_aux_ad % cfrac (iprofad) = scatt_aux_ad % cfrac (iprofad)        &
        & - sum ( scatt_aux_ad % clw (iprofad,:)  * clw_precf  (iprof,:)   &
        & +       scatt_aux_ad % ciw (iprofad,:)  * ciw_precf  (iprof,:)   &
        & +       scatt_aux_ad % totalice (iprofad,:)  * totalice_precf  (iprof,:)   &
        & +       scatt_aux_ad % rain (iprofad,:) * rain_precf (iprof,:)   &
        & +       scatt_aux_ad % sp (iprofad,:)   * sp_precf   (iprof,:) ) &
        & / ( scatt_aux % cfrac (iprof) ** 2 )

      scatt_aux_ad % clw  (iprofad,:) = scatt_aux_ad % clw  (iprofad,:) / scatt_aux % cfrac (iprof)
      scatt_aux_ad % ciw  (iprofad,:) = scatt_aux_ad % ciw  (iprofad,:) / scatt_aux % cfrac (iprof)
      scatt_aux_ad % totalice (iprofad,:) = scatt_aux_ad %totalice (iprofad,:) / scatt_aux % cfrac (iprof)
      scatt_aux_ad % rain (iprofad,:) = scatt_aux_ad % rain (iprofad,:) / scatt_aux % cfrac (iprof) 
      scatt_aux_ad % sp   (iprofad,:) = scatt_aux_ad % sp   (iprofad,:) / scatt_aux % cfrac (iprof) 

    Else
      scatt_aux_ad % clw  (iprofad,:) = 0.0_JPRB
      scatt_aux_ad % ciw  (iprofad,:) = 0.0_JPRB
      scatt_aux_ad % totalice  (iprofad,:) = 0.0_JPRB
      scatt_aux_ad % rain (iprofad,:) = 0.0_JPRB
      scatt_aux_ad % sp   (iprofad,:) = 0.0_JPRB
    Endif

    if( usercfrac ) then

      !* User-supplied cloud fraction
      cld_profiles_ad (iprofad) % cfrac = scatt_aux_ad % cfrac (iprofad)  
      scatt_aux_ad % cfrac (iprofad) = 0.0_JPRB

    else

      !* Weighted mean cloud fraction
      if (hydro_column(iprof) > 1e-10_JPRB) then

        if ( cfrac(iprof) < 0.0_JPRB ) scatt_aux_ad % cfrac (iprofad) = 0.0_JPRB
        if ( cfrac(iprof) > 1.0_JPRB ) scatt_aux_ad % cfrac (iprofad) = 0.0_JPRB

        cld_profiles_ad (iprofad) % cc (:) = cld_profiles_ad (iprofad) % cc (:) &
          & + scatt_aux_ad % cfrac (iprofad) * hydro_weights(iprof,:) / hydro_column(iprof)

        hydro_weights_ad(:) = scatt_aux_ad % cfrac (iprofad) * cld_profiles (iprof) % cc (:) &
          & / hydro_column(iprof)

        hydro_column_ad = -1.0_JPRB * scatt_aux_ad % cfrac (iprofad) &
          & * sum(hydro_weights(iprof,:) * cld_profiles (iprof) % cc (:)) &
          & / ( hydro_column(iprof) ** 2 )

      else
        hydro_weights_ad(:) = 0.0_JPRB
        hydro_column_ad     = 0.0_JPRB
      endif

      scatt_aux_ad % cfrac (iprofad) = 0.0_JPRB

      !* Partial column of hydrometeors in g/m^2

      hydro_weights_ad(:) = hydro_weights_ad(:) + hydro_column_ad
      hydro_column_ad = 0.0_JPRB

      scatt_aux_ad % rain (iprofad,:) = scatt_aux_ad % rain (iprofad,:) &
                & + hydro_weights_ad(:) * scatt_aux % dz (iprof,:)
      scatt_aux_ad % sp (iprofad,:)   = scatt_aux_ad % sp (iprofad,:)   &
                & + hydro_weights_ad(:) * scatt_aux % dz (iprof,:)
      scatt_aux_ad % ciw (iprofad,:)  = scatt_aux_ad % ciw (iprofad,:)  &
                & + hydro_weights_ad(:) * scatt_aux % dz (iprof,:)
      scatt_aux_ad % totalice (iprofad,:)  = scatt_aux_ad % totalice (iprofad,:)  &
                & + hydro_weights_ad(:) * scatt_aux % dz (iprof,:)
      scatt_aux_ad % clw (iprofad,:)  = scatt_aux_ad % clw (iprofad,:)  &
                & + hydro_weights_ad(:) * scatt_aux % dz (iprof,:)

      scatt_aux_ad % dz (iprofad,:) = scatt_aux_ad % dz (iprofad,:) &
        & + hydro_weights_ad(:) * ( clw_precf (iprof,:) + ciw_precf (iprof,:) &
        & + rain_precf (iprof,:) + sp_precf (iprof,:) + totalice_precf (iprof,:))

      hydro_weights_ad(:) = 0.0_JPRB

    endif

  Enddo

  do ilayer = 1,nlevels
     do iprofad = 1, nprofilesad

        if (adk == adk_adjoint) then
          iprof = iprofad  
        else if (adk == adk_k) then
          iprof = chanprof(iprofad) % prof  
        endif

        mmr_to_density_ad = 0.0_JPRB
   
        if( cld_profiles (iprof) % mmr_snowrain ) then

          mmr_to_density_ad = mmr_to_density_ad + rain_scale (iprof,ilayer) * scatt_aux_ad % rain (iprofad,ilayer)
          mmr_to_density_ad = mmr_to_density_ad + sp_scale   (iprof,ilayer) * scatt_aux_ad % sp   (iprofad,ilayer)
          scatt_aux_ad % rain (iprofad,ilayer) = scatt_aux_ad % rain (iprofad,ilayer) * mmr_to_density(iprof,ilayer)  
          scatt_aux_ad % sp   (iprofad,ilayer) = scatt_aux_ad % sp   (iprofad,ilayer) * mmr_to_density(iprof,ilayer)  

        else

          !* Rates from kg/m^2/s to g/m^3
          if (sp_scale   (iprof,ilayer) > 0.0_JPRB) then 
            scatt_aux_ad % sp   (iprofad,ilayer) = scatt_aux_ad % sp   (iprofad,ilayer) &
              & * (coef_scatt % conv_sp (2)) * (sp_scale (iprof,ilayer)**(coef_scatt % conv_sp (2) - 1.0_JPRB)) &
              & * (coef_scatt % conv_sp (1))**(coef_scatt % conv_sp (2)) 
          else
            scatt_aux_ad % sp   (iprofad,ilayer) = 0.0_JPRB
          endif

          if (rain_scale (iprof,ilayer) > 0.0_JPRB) then 
            scatt_aux_ad % rain (iprofad,ilayer) = scatt_aux_ad % rain (iprofad,ilayer) &
              & * (coef_scatt % conv_rain (2)) * (rain_scale (iprof,ilayer)**(coef_scatt % conv_rain (2) - 1.0_JPRB)) &
              & * (coef_scatt % conv_rain (1))**(coef_scatt % conv_rain (2)) 
          else
              scatt_aux_ad % rain (iprofad,ilayer) = 0.0_JPRB
          endif

          scatt_aux_ad % sp   (iprofad,ilayer) = scatt_aux_ad % sp   (iprofad,ilayer) * 3600.0_JPRB
          scatt_aux_ad % rain (iprofad,ilayer) = scatt_aux_ad % rain (iprofad,ilayer) * 3600.0_JPRB 

          scatt_aux_ad % sp   (iprofad,ilayer) = scatt_aux_ad % sp   (iprofad,ilayer) / rho_snow
          scatt_aux_ad % rain (iprofad,ilayer) = scatt_aux_ad % rain (iprofad,ilayer) / rho_rain

        endif

        mmr_to_density_ad = mmr_to_density_ad + totalice_scale (iprof,ilayer) * scatt_aux_ad % totalice (iprofad,ilayer)
        mmr_to_density_ad = mmr_to_density_ad + ciw_scale      (iprof,ilayer) * scatt_aux_ad % ciw      (iprofad,ilayer)
        mmr_to_density_ad = mmr_to_density_ad + clw_scale      (iprof,ilayer) * scatt_aux_ad % clw      (iprofad,ilayer)

        scatt_aux_ad % totalice (iprofad,ilayer) = scatt_aux_ad % totalice (iprofad,ilayer) * mmr_to_density(iprof,ilayer) 
        scatt_aux_ad % ciw      (iprofad,ilayer) = scatt_aux_ad % ciw      (iprofad,ilayer) * mmr_to_density(iprof,ilayer) 
        scatt_aux_ad % clw      (iprofad,ilayer) = scatt_aux_ad % clw      (iprofad,ilayer) * mmr_to_density(iprof,ilayer)  

        presfad (iprofad,ilayer) = presfad (iprofad,ilayer) &
                               & + de2mr / profiles (iprof) % t (ilayer) * mmr_to_density_ad   
        profiles_ad (iprofad) % t (ilayer) = profiles_ad (iprofad) % t (ilayer) & 
                                         & - presf (iprof,ilayer) * de2mr &
                                         & / (profiles (iprof) % t (ilayer) * profiles (iprof) % t (ilayer)) &
                                         & * mmr_to_density_ad   

        mmr_to_density_ad = 0.0_JPRB

     enddo
  enddo

  !* Initialise cloud and rain properties of the cloudy/rainy column
  do ilayer=1,nlevels
    do iprofad = 1, nprofilesad 

      if (adk == adk_adjoint) then
        iprof = iprofad  
      else if (adk == adk_k) then
        iprof = chanprof(iprofad) % prof  
      endif

      cld_profiles_ad (iprofad) % clw  (ilayer) = cld_profiles_ad (iprofad) % clw  (ilayer) & 
                                            & + scatt_aux_ad % clw  (iprofad,ilayer)
      cld_profiles_ad (iprofad) % rain (ilayer) = cld_profiles_ad (iprofad) % rain (ilayer) & 
                                            & + scatt_aux_ad % rain (iprofad,ilayer)

      if ( cld_profiles (iprof) % use_totalice ) then
        cld_profiles_ad (iprofad) % totalice (ilayer) = cld_profiles_ad (iprofad) % totalice (ilayer) & 
                                                  & + scatt_aux_ad % totalice (iprofad,ilayer)
      else
        cld_profiles_ad (iprofad) % ciw (ilayer) = cld_profiles_ad (iprofad) % ciw  (ilayer) & 
                                             & + scatt_aux_ad % ciw  (iprofad,ilayer)
        cld_profiles_ad (iprofad) % sp  (ilayer) = cld_profiles_ad (iprofad) % sp   (ilayer) & 
                                             & + scatt_aux_ad % sp   (iprofad,ilayer)
      endif 

      scatt_aux_ad % totalice (iprofad,ilayer) = 0.0_JPRB
      scatt_aux_ad % ciw      (iprofad,ilayer) = 0.0_JPRB 
      scatt_aux_ad % sp       (iprofad,ilayer) = 0.0_JPRB
      scatt_aux_ad % clw      (iprofad,ilayer) = 0.0_JPRB
      scatt_aux_ad % rain     (iprofad,ilayer) = 0.0_JPRB

    enddo
  enddo

  if (lhook) call dr_hook('RTTOV_HYDRO_AD',1_jpim,zhook_handle)

end subroutine rttov_hydro_ad
