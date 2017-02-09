!
Subroutine rttov_hydro_tl (&
      & nlevels,           &! in
      & nprofiles,         &! in
      & usercfrac,         &! in
      & presf,             &! in
      & presftl,           &! in
      & profiles,          &! in  
      & profiles_tl,       &! in  
      & cld_profiles,      &! in 
      & cld_profiles_tl,   &! in 
      & coef_scatt,        &! in
      & scatt_aux,         &! inout 
      & scatt_aux_tl)       ! inout 

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
       & profile_cloud_Type  

  Use parkind1, Only : jpim, jplm, jprb

!INTF_OFF
  Use rttov_const, Only:      &
       & rgp,                 &
       & rm,                  &
       & rho_rain,            &
       & rho_snow,            &
       & ccthres 

  USE YOMHOOK, ONLY: LHOOK , DR_HOOK
!INTF_ON
  Implicit None
  
!* Subroutine arguments:
  Integer (Kind=jpim), Intent (in) :: nlevels   ! Number of levels
  Integer (Kind=jpim), Intent (in) :: nprofiles ! Number of profiles
  Logical (Kind=jplm), Intent (in) :: usercfrac               ! User has supplied cloud fraction
  Real (Kind=jprb),           Intent (in)    :: presf(nprofiles,nlevels) ! Pressure levels [hPa]
  Real (Kind=jprb),           Intent (in)    :: presftl(nprofiles,nlevels) ! Pressure levels [hPa]
  Type (profile_Type),        Intent (in)    :: profiles (nprofiles)     ! Atmospheric profiles
  Type (profile_Type),        Intent (in)    :: profiles_tl (nprofiles)  ! Atmospheric profiles
  Type (rttov_scatt_coef),    Intent (in)    :: coef_scatt               ! RTTOV_SCATT Coefficients
  Type (profile_cloud_Type),  Intent (in)    :: cld_profiles (nprofiles) ! Cloud profiles 
  Type (profile_cloud_Type),  Intent (in)    :: cld_profiles_tl (nprofiles) ! Cloud profiles TL
  Type (profile_scatt_aux),   Intent (inout) :: scatt_aux                ! Auxiliary profile variables
  Type (profile_scatt_aux),   Intent (inout) :: scatt_aux_tl             ! Auxiliary profile variables
!INTF_END

!* Local variables
  Integer (Kind=jpim) :: ilayer, iprof
  Real    (Kind=jprb) :: de2mr, mmr_to_density, mmr_to_density_tl
  Real    (Kind=jprb) :: hydro_weights(nlevels), hydro_column    ! Total hydrometeor amounts [g m^-2]
  Real    (Kind=jprb) :: hydro_weights_tl(nlevels), hydro_column_tl  

  Character (len=15) :: NameOfRoutine = 'rttov_hydro_tl '
  
  REAL(KIND=JPRB) :: ZHOOK_HANDLE
  
  !- End of header --------------------------------------------------------

  if (lhook) call dr_hook('RTTOV_HYDRO_TL',0_jpim,zhook_handle)

  de2mr =  1.0E+05_JPRB * rm / rgp

  !* Initialise cloud and rain properties of the cloudy/rainy column
  do ilayer=1,nlevels
    do iprof = 1, nprofiles  
      scatt_aux    % clw  (iprof,ilayer) = cld_profiles    (iprof) % clw  (ilayer) 
      scatt_aux    % rain (iprof,ilayer) = cld_profiles    (iprof) % rain (ilayer) 
      scatt_aux_tl % clw  (iprof,ilayer) = cld_profiles_tl (iprof) % clw  (ilayer) 
      scatt_aux_tl % rain (iprof,ilayer) = cld_profiles_tl (iprof) % rain (ilayer) 
      if ( cld_profiles (iprof) % use_totalice ) then
        scatt_aux    % totalice (iprof,ilayer) = cld_profiles    (iprof) % totalice (ilayer) 
        scatt_aux    % ciw      (iprof,ilayer) = 0.0_JPRB
        scatt_aux    % sp       (iprof,ilayer) = 0.0_JPRB
        scatt_aux_tl % totalice (iprof,ilayer) = cld_profiles_tl (iprof) % totalice (ilayer) 
        scatt_aux_tl % ciw      (iprof,ilayer) = 0.0_JPRB
        scatt_aux_tl % sp       (iprof,ilayer) = 0.0_JPRB
      else
        scatt_aux    % totalice (iprof,ilayer) = 0.0_JPRB
        scatt_aux    % ciw      (iprof,ilayer) = cld_profiles    (iprof) % ciw  (ilayer) 
        scatt_aux    % sp       (iprof,ilayer) = cld_profiles    (iprof) % sp   (ilayer) 
        scatt_aux_tl % totalice (iprof,ilayer) = 0.0_JPRB
        scatt_aux_tl % ciw      (iprof,ilayer) = cld_profiles_tl (iprof) % ciw  (ilayer) 
        scatt_aux_tl % sp       (iprof,ilayer) = cld_profiles_tl (iprof) % sp   (ilayer) 
      endif 
    enddo
  enddo

  !* Change units
  do ilayer = 1, nlevels
     do iprof = 1, nprofiles

        !* Condensate from g/g to g/m^3
        mmr_to_density    = presf   (iprof,ilayer) * de2mr / profiles    (iprof) % t (ilayer) 
        mmr_to_density_tl = presftl (iprof,ilayer) * de2mr / profiles    (iprof) % t (ilayer) &
                        & - presf   (iprof,ilayer) * de2mr * profiles_tl (iprof) % t (ilayer) &
                        & / (profiles (iprof) % t (ilayer) * profiles (iprof) % t (ilayer))
        
        scatt_aux_tl % clw (iprof,ilayer) = scatt_aux_tl % clw (iprof,ilayer) * mmr_to_density &
                                        & + scatt_aux    % clw (iprof,ilayer) * mmr_to_density_tl
        scatt_aux    % clw (iprof,ilayer) = scatt_aux    % clw (iprof,ilayer) * mmr_to_density

        scatt_aux_tl % ciw (iprof,ilayer) = scatt_aux_tl % ciw (iprof,ilayer) * mmr_to_density &
                                        & + scatt_aux    % ciw (iprof,ilayer) * mmr_to_density_tl
        scatt_aux    % ciw (iprof,ilayer) = scatt_aux    % ciw (iprof,ilayer) * mmr_to_density

        scatt_aux_tl % totalice (iprof,ilayer) = scatt_aux_tl % totalice (iprof,ilayer) * mmr_to_density &
                                             & + scatt_aux    % totalice (iprof,ilayer) * mmr_to_density_tl
        scatt_aux    % totalice (iprof,ilayer) = scatt_aux    % totalice (iprof,ilayer) * mmr_to_density

        if( cld_profiles (iprof) % mmr_snowrain ) then

          scatt_aux_tl % rain (iprof,ilayer) = scatt_aux_tl % rain (iprof,ilayer) * mmr_to_density &
                                           & + scatt_aux    % rain (iprof,ilayer) * mmr_to_density_tl
          scatt_aux    % rain (iprof,ilayer) = scatt_aux    % rain (iprof,ilayer) * mmr_to_density

          scatt_aux_tl % sp (iprof,ilayer) = scatt_aux_tl % sp (iprof,ilayer) * mmr_to_density &
                                         & + scatt_aux    % sp (iprof,ilayer) * mmr_to_density_tl
          scatt_aux    % sp (iprof,ilayer) = scatt_aux    % sp (iprof,ilayer) * mmr_to_density 

        else

          !* Rates from kg/m^2/s to g/m^3
          scatt_aux_tl % rain (iprof,ilayer) = scatt_aux_tl % rain (iprof,ilayer) / rho_rain
          scatt_aux    % rain (iprof,ilayer) = scatt_aux    % rain (iprof,ilayer) / rho_rain
          scatt_aux_tl % sp   (iprof,ilayer) = scatt_aux_tl % sp   (iprof,ilayer) / rho_snow
          scatt_aux    % sp   (iprof,ilayer) = scatt_aux    % sp   (iprof,ilayer) / rho_snow

          scatt_aux_tl % rain (iprof,ilayer) = scatt_aux_tl % rain (iprof,ilayer) * 3600.0_JPRB 
          scatt_aux    % rain (iprof,ilayer) = scatt_aux    % rain (iprof,ilayer) * 3600.0_JPRB 
          scatt_aux_tl % sp   (iprof,ilayer) = scatt_aux_tl % sp   (iprof,ilayer) * 3600.0_JPRB
          scatt_aux    % sp   (iprof,ilayer) = scatt_aux    % sp   (iprof,ilayer) * 3600.0_JPRB

          if (scatt_aux    % rain (iprof,ilayer) > 0.0_JPRB) then
            scatt_aux_tl % rain (iprof,ilayer) =  scatt_aux_tl % rain (iprof,ilayer) &
              & * (coef_scatt % conv_rain (2)) * (scatt_aux % rain (iprof,ilayer)**(coef_scatt % conv_rain (2) - 1.0_JPRB)) &
               & * (coef_scatt % conv_rain (1))**(coef_scatt % conv_rain (2)) 
            scatt_aux    % rain (iprof,ilayer) = (scatt_aux    % rain (iprof,ilayer) &
                                          & *  coef_scatt % conv_rain (1))**(coef_scatt % conv_rain (2)) 
          else
            scatt_aux_tl % rain (iprof,ilayer) = 0.0_JPRB
          endif     
          if (scatt_aux    % sp   (iprof,ilayer) > 0.0_JPRB) then
            scatt_aux_tl % sp   (iprof,ilayer) = scatt_aux_tl % sp    (iprof,ilayer) &
              & * (coef_scatt % conv_sp   (2)) * (scatt_aux % sp   (iprof,ilayer)**(coef_scatt % conv_sp   (2) - 1.0_JPRB)) &
              & * (coef_scatt % conv_sp   (1))**(coef_scatt % conv_sp   (2)) 
            scatt_aux    % sp   (iprof,ilayer) = (scatt_aux   % sp    (iprof,ilayer) &
              & *  coef_scatt % conv_sp   (1))**(coef_scatt%conv_sp     (2))
          else
            scatt_aux_tl % sp   (iprof,ilayer) = 0.0_JPRB
          endif   

        endif

     enddo
  enddo

  scatt_aux    % cfrac (:) = 0.0_JPRB
  scatt_aux_tl % cfrac (:) = 0.0_JPRB

  do iprof = 1, nprofiles

    if( usercfrac ) then

      !* User-supplied cloud fraction
      scatt_aux    % cfrac (iprof) = cld_profiles    (iprof) % cfrac
      scatt_aux_tl % cfrac (iprof) = cld_profiles_tl (iprof) % cfrac

    else

      !* Calculate a hydrometeor-weighted average cloudy sky fraction 

      !* Partial column of hydrometeors in g/m^2
      hydro_weights(:) = &
        & ( scatt_aux % rain (iprof,:) + scatt_aux % sp (iprof,:) &
        & + scatt_aux % ciw (iprof,:)  + scatt_aux % clw (iprof,:) &
        & + scatt_aux % totalice (iprof,:) ) &
        & * scatt_aux % dz (iprof,:)
      hydro_column = sum(hydro_weights)

      hydro_weights_tl(:) = &
        & ( scatt_aux_tl % rain (iprof,:) + scatt_aux_tl % sp (iprof,:) &
        & + scatt_aux_tl % ciw (iprof,:)  + scatt_aux_tl % clw (iprof,:) &
        & + scatt_aux_tl % totalice (iprof,:) ) &
        & * scatt_aux % dz (iprof,:) &
        & + ( scatt_aux % rain (iprof,:) + scatt_aux % sp (iprof,:) &
        & + scatt_aux % ciw (iprof,:)  + scatt_aux % clw (iprof,:) &
        & + scatt_aux % totalice (iprof,:) ) &
        & * scatt_aux_tl % dz (iprof,:) 
      hydro_column_tl = sum(hydro_weights_tl)

      !* Weighted mean cloud fraction    
      if (hydro_column > 1e-10_JPRB) then

        scatt_aux % cfrac (iprof) = &
          & sum(hydro_weights(:) * cld_profiles (iprof) % cc (:)) / hydro_column 

        scatt_aux_tl % cfrac (iprof) = &
          & + sum(hydro_weights_tl(:) * cld_profiles (iprof) % cc (:)) / hydro_column &
          & + sum(hydro_weights(:) * cld_profiles_tl (iprof) % cc (:)) / hydro_column &
          & - hydro_column_tl * sum(hydro_weights(:) * cld_profiles (iprof) % cc (:)) &
          & / ( hydro_column**2 )

        if ( scatt_aux % cfrac (iprof) < 0.0_JPRB ) then
          scatt_aux    % cfrac (iprof) = 0.0_JPRB
          scatt_aux_tl % cfrac (iprof) = 0.0_JPRB
        endif

        if ( scatt_aux % cfrac (iprof) > 1.0_JPRB ) then
          scatt_aux    % cfrac (iprof) = 1.0_JPRB
          scatt_aux_tl % cfrac (iprof) = 0.0_JPRB
        endif

      else
        scatt_aux % cfrac (iprof)    = 0.0_JPRB
        scatt_aux_tl % cfrac (iprof) = 0.0_JPRB
      endif

    endif

    !* Partition all cloud and rain into the cloudy column
    If (scatt_aux % cfrac (iprof) > ccthres) Then

      scatt_aux_tl % clw  (iprof,:) = &
        &   scatt_aux_tl % clw  (iprof,:) / scatt_aux % cfrac (iprof) &
        & - scatt_aux_tl % cfrac (iprof) * scatt_aux % clw  (iprof,:) &
        & / (scatt_aux % cfrac (iprof)**2)
      scatt_aux_tl % ciw  (iprof,:) = &
        &   scatt_aux_tl % ciw  (iprof,:) / scatt_aux % cfrac (iprof) &
        & - scatt_aux_tl % cfrac (iprof) * scatt_aux % ciw  (iprof,:) &
        & / (scatt_aux % cfrac (iprof)**2)
      scatt_aux_tl % totalice  (iprof,:) = &
        &   scatt_aux_tl % totalice  (iprof,:) / scatt_aux % cfrac (iprof) &
        & - scatt_aux_tl % cfrac (iprof) * scatt_aux % totalice  (iprof,:) &
        & / (scatt_aux % cfrac (iprof)**2)
      scatt_aux_tl % rain  (iprof,:) = &
        &   scatt_aux_tl % rain  (iprof,:) / scatt_aux % cfrac (iprof) &
        & - scatt_aux_tl % cfrac (iprof) * scatt_aux % rain  (iprof,:) &
        & / (scatt_aux % cfrac (iprof)**2)
      scatt_aux_tl % sp  (iprof,:) = &
        &   scatt_aux_tl % sp  (iprof,:) / scatt_aux % cfrac (iprof) &
        & - scatt_aux_tl % cfrac (iprof) * scatt_aux % sp  (iprof,:) &
        & / (scatt_aux % cfrac (iprof)**2)

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
      scatt_aux_tl % clw  (iprof,:) = 0.0_JPRB
      scatt_aux_tl % ciw  (iprof,:) = 0.0_JPRB
      scatt_aux_tl % totalice  (iprof,:) = 0.0_JPRB
      scatt_aux_tl % rain (iprof,:) = 0.0_JPRB
      scatt_aux_tl % sp   (iprof,:) = 0.0_JPRB
    Endif
  Enddo

  if (lhook) call dr_hook('RTTOV_HYDRO_TL',1_jpim,zhook_handle)

end subroutine rttov_hydro_tl
