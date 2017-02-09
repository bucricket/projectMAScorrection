Subroutine rttov_iniedd_tl (&        
     & nlevels,       &! in
     & nchannels ,    &! in
     & nprofiles ,    &! in
     & lprofiles ,    &! in
     & angles    ,    &! in
     & scatt_aux,     &! inout
     & scatt_aux_tl)   ! inout 

  ! Description:
  ! to compute variables specific to Eddington approximation to RT
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
  ! - Bauer, P., 2002: Microwave radiative transfer modeling in clouds and precipitation.
  !     Part I: Model description.
  !     NWP SAF Report No. NWPSAF-EC-TR-005, 27 pp.
  ! - Chevallier, F., and P. Bauer, 2003:
  !     Model rain and clouds over oceans: comparison with SSM/I observations.
  !     Mon. Wea. Rev., 131, 1240-1255.
  ! - Moreau, E., P. Bauer and F. Chevallier, 2002: Microwave radiative transfer modeling in clouds and precipitation.
  !     Part II: Model evaluation.
  !     NWP SAF Report No. NWPSAF-EC-TR-006, 27 pp.
  !
  ! Current Code Owner: SAF NWP
  !
  ! History:
  ! Version   Date     Comment
  ! -------   ----     -------
  !  1.0       09/2002   Initial version     (P. Bauer, E. Moreau)
  !  1.1       05/2003   RTTOV7.3 compatible (F. Chevallier)
  !  1.2       03/2004   Added polarimetry   (R. Saunders)
  !  1.3       08/2004   Polarimetry fixes   (U. O'Keefe)
  !  1.4       11/2004   Clean-up            (P. Bauer)
  !  1.5       11/2007   RTTOV9 version      (A. Geer)
  !  1.6       07/2008   Consistent mininimum SSA (A. Geer)
  !  1.7       07/2008   Clear sky speed-ups (A. Geer)
  !  1.8       03/2010   Optimisation + don't delta scale outside Eddington (A. Geer)
  !
  ! Code Description:
  !   Language:           Fortran 90.
  !   Software Standards: "European Standards for Writing and
  !   Documenting Exchangeable Fortran 90 Code".
  !
  ! Declarations:
  ! Modules used:
  ! Imported Type Definitions:

  Use rttov_types, Only :    &
       & geometry_Type        ,&
       & profile_scatt_aux 

  Use parkind1, Only : jpim
!INTF_OFF
  Use rttov_const, Only : min_ssa, ccthres

  Use parkind1, Only : jprb

  USE YOMHOOK, ONLY: LHOOK , DR_HOOK
!INTF_ON
  Implicit None

!* Subroutine arguments:
  Integer (Kind=jpim), Intent (in) :: nlevels               ! Number of levels
  Integer (Kind=jpim), Intent (in) :: nprofiles             ! Number of profiles
  Integer (Kind=jpim), Intent (in) :: nchannels             ! Number of radiances
  Integer (Kind=jpim), Intent (in) :: lprofiles (nchannels) ! Profile indices

  Type (geometry_Type),     Intent (in)    :: angles (nprofiles) ! Zenith angles
  Type (profile_scatt_aux), Intent (inout) :: scatt_aux          ! Auxiliary profile variables
  Type (profile_scatt_aux), Intent (inout) :: scatt_aux_tl       ! Auxiliary profile variables

!INTF_END

!* Local variables
  Integer(Kind=jpim) :: ilayer, iprof, ichan, ii

  REAL(KIND=JPRB) :: ZHOOK_HANDLE

  !- End of header --------------------------------------------------------

  if (lhook) call dr_hook('RTTOV_INIEDD_TL',0_jpim,zhook_handle)

  scatt_aux_tl % delta   = 0.0_JPRB
  scatt_aux_tl % lambda  = 0.0_JPRB
  scatt_aux_tl % h       = 0.0_JPRB
  scatt_aux_tl % tau     = 0.0_JPRB
  
  scatt_aux    % delta   = 0.0_JPRB
  scatt_aux    % lambda  = 0.0_JPRB
  scatt_aux    % h       = 0.0_JPRB
  scatt_aux    % tau     = 1.0_JPRB
  
  scatt_aux    % mclayer = 0

!* Layer interface temperatures, lapse rates
  do ilayer = 1, nlevels
   do ii=1,nprofiles
     scatt_aux_tl % b0 (ii,ilayer) =  scatt_aux_tl % tbd (ii,ilayer+1)
     scatt_aux    % b0 (ii,ilayer) =  scatt_aux    % tbd (ii,ilayer+1)
     scatt_aux_tl % bn (ii,ilayer) =  scatt_aux_tl % tbd (ii,ilayer  )
     scatt_aux    % bn (ii,ilayer) =  scatt_aux    % tbd (ii,ilayer  )
     scatt_aux_tl % b1 (ii,ilayer) = (scatt_aux_tl % bn  (ii,ilayer)  -  scatt_aux_tl % b0 (ii,ilayer)) &
                                & /  scatt_aux    % dz  (ii,ilayer)  - (scatt_aux    % bn (ii,ilayer)  &
                                & -  scatt_aux    % b0  (ii,ilayer)) *  scatt_aux_tl % dz (ii,ilayer)  &
                                & /  scatt_aux    % dz  (ii,ilayer)  /  scatt_aux    % dz (ii,ilayer)
     scatt_aux    % b1 (ii,ilayer) = (scatt_aux    % bn  (ii,ilayer)  -  scatt_aux    % b0 (ii,ilayer)) &
                               &  /  scatt_aux    % dz  (ii,ilayer)
   end do
  end do

  do ichan = 1, nchannels
!* Cloud top level index
    scatt_aux % mclayer (ichan) = nlevels+1 
    do ilayer = 1, nlevels 
      if (scatt_aux % ssa (ichan,ilayer) > min_ssa ) then 
        scatt_aux % mclayer (ichan) = ilayer
        exit
      endif
    end do
    if (scatt_aux % mclayer (ichan) > nlevels-2 .and. scatt_aux % mclayer (ichan) /= nlevels+1) &
      & scatt_aux % mclayer (ichan) = nlevels-2 !* DGBF imposes minimum number of layers in Eddington
  end do 

!* Delta-scaling
  do ilayer = 1, nlevels
    do ichan = 1, nchannels
      iprof = lprofiles (ichan)
      if ( scatt_aux % cfrac(iprof) > ccthres ) then 

        if (ilayer >= scatt_aux % mclayer (ichan)) then   
          scatt_aux_tl % ext  (ichan,ilayer) = (1.0_JPRB - scatt_aux % ssa (ichan,ilayer) * scatt_aux % asm (ichan,ilayer) &
                                   & * scatt_aux % asm (ichan,ilayer)) * scatt_aux_tl % ext  (ichan,ilayer) &
                                   & + scatt_aux % ext (ichan,ilayer) * scatt_aux % asm (ichan,ilayer) &
                                   & * (-1.0_JPRB * scatt_aux_tl % ssa  (ichan,ilayer) * scatt_aux % asm (ichan,ilayer) &
                                   & -2.0_JPRB * scatt_aux % ssa (ichan,ilayer) * scatt_aux_tl % asm  (ichan,ilayer)) 
          scatt_aux    % ext  (ichan,ilayer) = (1.0_JPRB - scatt_aux % ssa (ichan,ilayer) * scatt_aux % asm (ichan,ilayer) &
                                   & * scatt_aux % asm (ichan,ilayer)) * scatt_aux % ext (ichan,ilayer) 
          scatt_aux_tl % ssa  (ichan,ilayer) = (1.0_JPRB /  (1.0_JPRB - scatt_aux % asm (ichan,ilayer) &
                                   & * scatt_aux % asm (ichan,ilayer) &
                                   & * scatt_aux % ssa (ichan,ilayer))**2.0_JPRB ) * ( (1.0_JPRB - scatt_aux % asm (ichan,ilayer) &
                                   & * scatt_aux % asm (ichan,ilayer) * scatt_aux % ssa (ichan,ilayer))         &
                                   & * (( 1.0_JPRB - scatt_aux % asm (ichan,ilayer) * scatt_aux % asm (ichan,ilayer)) *&
                                   & scatt_aux_tl % ssa (ichan,ilayer) &
                                   & +  (-2.0_JPRB * scatt_aux % asm (ichan,ilayer) * scatt_aux_tl % asm (ichan,ilayer) * &
                                   &scatt_aux % ssa (ichan,ilayer)) ) &
                                   & - (1.0_JPRB - scatt_aux % asm (ichan,ilayer) * scatt_aux % asm (ichan,ilayer)) * &
                                   &scatt_aux % ssa (ichan,ilayer) &
                                   & * (  (-1.0_JPRB * scatt_aux % asm (ichan,ilayer) * scatt_aux % asm (ichan,ilayer) * &
                                   &scatt_aux_tl % ssa  (ichan,ilayer)) &
                                   & +  (-2.0_JPRB * scatt_aux % asm (ichan,ilayer) * scatt_aux_tl % asm  (ichan,ilayer) * &
                                   &scatt_aux % ssa (ichan,ilayer)) ) ) 
          scatt_aux    % ssa  (ichan,ilayer) = (1.0_JPRB - scatt_aux % asm (ichan,ilayer) * scatt_aux % asm (ichan,ilayer)) &
                                   & * scatt_aux % ssa (ichan,ilayer) / (1.0_JPRB - scatt_aux % asm (ichan,ilayer)  &
                                   & * scatt_aux % asm (ichan,ilayer) * scatt_aux % ssa (ichan,ilayer)) 
          scatt_aux_tl % asm  (ichan,ilayer) = scatt_aux_tl % asm (ichan,ilayer) / (1.0_JPRB + scatt_aux % asm (ichan,ilayer)) &
                                   & - scatt_aux % asm (ichan,ilayer) * scatt_aux_tl % asm (ichan,ilayer) / &
                                   &(1.0_JPRB + scatt_aux % asm (ichan,ilayer)) &
                                   & / (1.0_JPRB + scatt_aux % asm (ichan,ilayer)) 
          scatt_aux    % asm  (ichan,ilayer) = scatt_aux % asm (ichan,ilayer) / (1.0_JPRB + scatt_aux % asm (ichan,ilayer))

          scatt_aux_tl % lambda (ichan,ilayer) = (1.0_JPRB / ( 2._JPRB * sqrt (3.0_JPRB * scatt_aux % ext (ichan,ilayer) *&
                                     &scatt_aux % ext (ichan,ilayer)  &
                                     & * (1.0_JPRB - scatt_aux % ssa (ichan,ilayer)) &
                                     & * (1.0_JPRB - scatt_aux % ssa (ichan,ilayer) * scatt_aux % asm (ichan,ilayer))))) &
                                     & * ( 6.0_JPRB * scatt_aux % ext (ichan,ilayer) * scatt_aux_tl % ext (ichan,ilayer) &
                                     & * (1.0_JPRB - scatt_aux % ssa (ichan,ilayer)) * &
                                     &(1.0_JPRB - scatt_aux % ssa (ichan,ilayer) * scatt_aux % asm (ichan,ilayer)) &
                                     & - 3.0_JPRB * scatt_aux % ext (ichan,ilayer) * scatt_aux % ext(ichan,ilayer) &
                                     & * scatt_aux_tl % ssa (ichan,ilayer) *  &
                                     &(1.0_JPRB - scatt_aux % ssa (ichan,ilayer) * scatt_aux % asm (ichan,ilayer)) &
                                     & - 3.0_JPRB * scatt_aux % ext (ichan,ilayer) * scatt_aux % ext (ichan,ilayer) * &
                                     &(1.0_JPRB - scatt_aux % ssa (ichan,ilayer)) &
                                     & * (scatt_aux % ssa (ichan,ilayer) * scatt_aux_tl % asm (ichan,ilayer)+ &
                                     &scatt_aux_tl % ssa (ichan,ilayer) &
                                     & * scatt_aux % asm (ichan,ilayer)) ) 

          scatt_aux    % lambda (ichan,ilayer) = sqrt (3.0_JPRB * scatt_aux % ext (ichan,ilayer) * scatt_aux % ext (ichan,ilayer) &
                                     & * (1.0_JPRB - scatt_aux % ssa (ichan,ilayer)) &
                                     & * (1.0_JPRB - scatt_aux % ssa (ichan,ilayer) * scatt_aux % asm (ichan,ilayer))) 

          scatt_aux_tl % h      (ichan,ilayer) = 1.5_JPRB * scatt_aux_tl % ext (ichan,ilayer) &
                                     & * (1.0_JPRB - scatt_aux % ssa (ichan,ilayer) * scatt_aux % asm (ichan,ilayer)) &
                                     & -1.5_JPRB * scatt_aux % ext (ichan,ilayer) &
                                     & * (scatt_aux_tl % ssa (ichan,ilayer) * scatt_aux % asm (ichan,ilayer) + &
                                     &scatt_aux % ssa (ichan,ilayer) * scatt_aux_tl % asm (ichan,ilayer)) 
          scatt_aux    % h      (ichan,ilayer) = 1.5_JPRB * scatt_aux % ext (ichan,ilayer) &
                                     & * (1.0_JPRB - scatt_aux % ssa (ichan,ilayer) * scatt_aux % asm (ichan,ilayer)) 

          if (scatt_aux    % h (ichan,ilayer) < 0.00001_JPRB) then
              scatt_aux_tl % h (ichan,ilayer) = 0.0_JPRB
              scatt_aux    % h (ichan,ilayer) = 0.00001_JPRB
          endif

        endif     

        scatt_aux_tl % delta (ichan,ilayer) = (scatt_aux_tl % ext (ichan,ilayer) * scatt_aux    % dz (iprof,ilayer) &
                                  & +  scatt_aux    % ext (ichan,ilayer) * scatt_aux_tl % dz (iprof,ilayer)) / &
                                  &angles (iprof) % coszen 
        scatt_aux    % delta (ichan,ilayer) = (scatt_aux    % ext (ichan,ilayer) * scatt_aux    % dz (iprof,ilayer)) /&
         & angles (iprof) % coszen

        if (scatt_aux    % delta (ichan,ilayer) >= 30.0_JPRB) then
            scatt_aux    % delta (ichan,ilayer)  = 30.0_JPRB
            scatt_aux_tl % delta (ichan,ilayer)  =  0.0_JPRB
        endif
  
        scatt_aux    % tau    (ichan,ilayer) =  1.0_JPRB / exp (scatt_aux % delta (ichan,ilayer))
        scatt_aux_tl % tau    (ichan,ilayer) = -1.0_JPRB * scatt_aux_tl % delta (ichan,ilayer) * scatt_aux % tau    (ichan,ilayer)

      endif
    end do
  end do
 
  if (lhook) call dr_hook('RTTOV_INIEDD_TL',1_jpim,zhook_handle)
  
End subroutine rttov_iniedd_tl
