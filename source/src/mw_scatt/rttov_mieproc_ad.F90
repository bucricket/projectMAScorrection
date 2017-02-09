!
Subroutine rttov_mieproc_ad (&        
     & nlevels,           &! in
     & nchannels,         &! in
     & nprofiles,         &! in
     & nprofilesad,       &! in
     & frequencies,       &! in
     & lprofiles,         &! in
     & profiles,          &! in
     & coef_scatt,        &! in
     & scatt_aux,         &! inout
     & scatt_aux_ad)       ! inout 
  !
  ! Description:
  ! Calculates scattering parameters from Mie tables
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
  ! - Bauer, P., 2002: Microwave radiative transfer modeling in clouds and 
  !     precipitation. Part I: Model description.
  !     NWP SAF Report No. NWPSAF-EC-TR-005, 27 pp.
  ! - Chevallier, F., and P. Bauer, 2003:
  !     Model rain and clouds over oceans: comparison with SSM/I observations.
  !     Mon. Wea. Rev., 131, 1240-1255.
  ! - Moreau, E., P. Bauer and F. Chevallier, 2002: Microwave radiative transfer 
  !     modeling in clouds and precipitation. Part II: Model evaluation.
  !     NWP SAF Report No. NWPSAF-EC-TR-006, 27 pp.
  !
  ! Current Code Owner: SAF NWP
  !
  ! History:
  ! Version   Date     Comment
  ! -------   ----     -------
  !  1.0    09/2002      Initial version         (E. Moreau)
  !  1.1    05/2003      RTTOV7.3 compatible     (F. Chevallier)
  !  1.3    03/2004      Polarimetry code added  (R. Saunders)
  !  1.4    11/2004      Clean-up                (P. Bauer)
  !  1.5    10/2006      Introduce interpolation to zero for scatt coeffs 
  !                      for small LWC (U.O'Keeffe). Cleaner version (A. Geer)
  !  1.6    11/2007      RTTOV9 versions         (A. Geer)
  !  1.7    06/2008      Fix 2 rarely occurring adjoint bugs (A. Geer)
  !  1.8    06/2008      Performance enhancements (D. Salmond)
  !  1.9    06/2008      Fix minor bug in small LWC interpolation (A. Geer)
  !  1.10   07/2008      Clear sky speed-ups     (A. Geer)
  !  1.11   03/2009      Prevent extrapolation beyond table bounds (A. Geer)
  !  1.12   03/2010      Optimisation (A. Geer)
  !
  ! Code Description:
  !   Language:           Fortran 90.
  !   Software Standards: "European Standards for Writing and
  !   Documenting Exchangeable Fortran 90 Code".
  !
  !
  ! Declarations:
  ! Modules used:
  ! Imported Type Definitions:

  Use rttov_types, Only :    &
       & profile_scatt_aux    ,&
       & profile_Type   ,&
       & rttov_scatt_coef

  Use parkind1, Only : jpim
!INTF_OFF
  Use rttov_const, Only: adk_adjoint, adk_k, ccthres

  Use parkind1, Only : jprb

  USE YOMHOOK, ONLY: LHOOK , DR_HOOK
!INTF_ON
  Implicit None

!* Subroutine arguments:
  Integer (Kind=jpim), Intent (in) :: nlevels                 ! Number oflevels
  Integer (Kind=jpim), Intent (in) :: nchannels               ! Number of channels*profiles=radiances
  Integer (Kind=jpim), Intent (in) :: nprofiles               ! Number of profiles
  Integer (Kind=jpim), Intent (in) :: nprofilesad             ! Number of profiles in adjoint variables
  Integer (Kind=jpim), Intent (in) :: frequencies (nchannels) ! Frequency indices
  Integer (Kind=jpim), Intent (in) :: lprofiles   (nchannels) ! Profile indices

  Type (profile_Type),      Intent (in)    :: profiles    (nprofiles)  
  Type (rttov_scatt_coef),  Intent (in)    :: coef_scatt               ! RTTOV_SCATT Coefficients
  Type (profile_scatt_aux), Intent (inout) :: scatt_aux                ! Auxiliary profile variables
  Type (profile_scatt_aux), Intent (inout) :: scatt_aux_ad             ! Auxiliary profile variables

!INTF_END

!* Local variables:
  Integer (Kind=jpim) :: iwc(coef_scatt % nhydro), itemp(coef_scatt % nhydro), itype, ichan, iprof, ifreq, ilayer
  Integer (Kind=jpim) :: ifreq_last, iprof_last
  Integer (Kind=jpim) :: adk, iprofad
  Real    (Kind=jprb) :: wc(coef_scatt % nhydro)   , temp   ,  kp   , ap   , gp   , s_k   , s_a   , s_g   , zln10
  Real    (Kind=jprb) :: kpp,app,gpp
  Real    (Kind=jprb) :: wc_ad,  kp_ad, ap_ad, gp_ad    
  Real    (Kind=jprb) :: cont(coef_scatt % nhydro), cont_ad(coef_scatt % nhydro), cont_min, offset_t(coef_scatt % nhydro)
  
  Real    (Kind=jprb) :: ext, ssa, asm

  REAL (KIND=JPRB) :: ZHOOK_HANDLE

  !- End of header --------------------------------------------------------      

  if (lhook) call dr_hook('RTTOV_MIEPROC_AD',0_jpim,zhook_handle)

  if (nprofilesad == nprofiles) then 
    adk = adk_adjoint   ! Adjoint mode
  else if (nprofilesad == nchannels) then
    adk = adk_k         ! K mode
  endif 

  zln10 = log (10.0_JPRB)
  cont_min = 10.0_JPRB ** ( (1.0_JPRB + coef_scatt % offset_water) / coef_scatt % scale_water )

  offset_t(1) = coef_scatt % offset_temp_rain
  offset_t(2) = coef_scatt % offset_temp_sp
  offset_t(3) = coef_scatt % offset_temp_liq
  offset_t(4) = coef_scatt % offset_temp_ice
  offset_t(5) = coef_scatt % offset_temp_totalice

  !* Loops over channels, levels, hydrometeor types
  nlayer_loop1: do ilayer = 1, nlevels

    ifreq_last = -1
    iprof_last = -1

    nchan_loop1: do ichan = 1, nchannels

      iprof = lprofiles   (ichan)
      ifreq = frequencies (ichan)  

      if( scatt_aux % cfrac(iprof) > ccthres ) then   

        if(iprof /= iprof_last) then

          cont(1)  = scatt_aux % rain (iprof,ilayer)
          cont(2)  = scatt_aux % sp  (iprof,ilayer) 
          cont(3)  = scatt_aux % clw (iprof,ilayer) 
          cont(4)  = scatt_aux % ciw (iprof,ilayer)
          cont(5)  = scatt_aux % totalice (iprof,ilayer)

          ntype_loop1: do itype = 1, coef_scatt % nhydro

            !* Nearest index for Mie-table: LWC/IWC
            if (cont(itype) >= cont_min) then 
              wc(itype) = coef_scatt % scale_water * log10 (cont(itype)) - coef_scatt % offset_water 
            else if (cont(itype) < cont_min .and. cont(itype) > 0.0_JPRB) then
              wc(itype) = cont(itype) / cont_min
            else
              wc(itype) = 0.0_JPRB
            endif

            iwc(itype) = floor (wc(itype))
            if (iwc(itype) > coef_scatt % mwc - 1) then
              ! Prevent extrapolation 
              iwc(itype) = coef_scatt % mwc - 1
              wc(itype)  = coef_scatt % mwc
            endif

            !* Nearest index for Mie-table: T (w/o melting layer)
            temp = profiles (iprof) % t (ilayer) - offset_t(itype)

            itemp(itype) = anint (temp)
            if (itemp(itype) <                      1) itemp(itype) = 1
            if (itemp(itype) > coef_scatt % mtemp - 1) itemp(itype) = coef_scatt % mtemp - 1

          enddo ntype_loop1
        endif

        if(ifreq /= ifreq_last .or. iprof /= iprof_last) then

          kpp=0.0_JPRB
          app=0.0_JPRB
          gpp=0.0_JPRB

          ntype_loop2: do itype = 1, coef_scatt % nhydro

            if (iwc(itype) >= 1) then

              s_k = coef_scatt % ext (ifreq,itype,itemp(itype),iwc(itype)+1) &
                & - coef_scatt % ext (ifreq,itype,itemp(itype),iwc(itype))
              s_a = coef_scatt % ssa (ifreq,itype,itemp(itype),iwc(itype)+1) &
                & - coef_scatt % ssa (ifreq,itype,itemp(itype),iwc(itype))
              s_g = coef_scatt % asp (ifreq,itype,itemp(itype),iwc(itype)+1) &
                & - coef_scatt % asp (ifreq,itype,itemp(itype),iwc(itype))

              kp  = coef_scatt % ext (ifreq,itype,itemp(itype),iwc(itype)) + s_k * (wc(itype) - iwc(itype))
              ap  = coef_scatt % ssa (ifreq,itype,itemp(itype),iwc(itype)) + s_a * (wc(itype) - iwc(itype))
              gp  = coef_scatt % asp (ifreq,itype,itemp(itype),iwc(itype)) + s_g * (wc(itype) - iwc(itype))

            else

              ! For small water contents, interpolate linearly to zero 
              s_k   = coef_scatt % ext (ifreq,itype,itemp(itype),1) 
              s_a   = coef_scatt % ssa (ifreq,itype,itemp(itype),1) 
              s_g   = coef_scatt % asp (ifreq,itype,itemp(itype),1) 

              kp  = max( s_k * wc(itype), 1E-10_JPRB )              
              if (wc(itype) > 1E-10_JPRB) then
                ap = s_a * wc(itype)  
                gp = s_g * wc(itype)  
              else
                ap = 0.0_JPRB
                gp = 0.0_JPRB
              endif             

            endif
            kpp=kpp+kp
            app=app+kp * ap
            gpp=gpp+kp * ap * gp
          enddo ntype_loop2
        endif
     
        ifreq_last=ifreq
        iprof_last=iprof
        scatt_aux % ext (ichan,ilayer) = scatt_aux % ext (ichan,ilayer) + kpp
        scatt_aux % ssa (ichan,ilayer) = scatt_aux % ssa (ichan,ilayer) + app
        scatt_aux % asm (ichan,ilayer) = scatt_aux % asm (ichan,ilayer) + gpp

      endif
    enddo nchan_loop1
  enddo nlayer_loop1
  
! ext (:,:) = scatt_aux % ext (:,:)
! ssa (:,:) = scatt_aux % ssa (:,:)
! asm (:,:) = scatt_aux % asm (:,:)

  do ilayer = 1, nlevels
   do ichan = 1, nchannels
     ext  = scatt_aux % ext (ichan,ilayer)
     ssa  = scatt_aux % ssa (ichan,ilayer)
     asm  = scatt_aux % asm (ichan,ilayer)

!* ADJOINT PART
     if (ext  >= 20.0_JPRB) then
        scatt_aux % ext (ichan,ilayer)  = 20.0_JPRB
        scatt_aux_ad % ext (ichan,ilayer) = 0.0_JPRB
     endif

     if (ssa  > 0.0_JPRB ) then
        scatt_aux % ssa (ichan,ilayer) =  ssa /  ext 
        scatt_aux_ad % ext (ichan,ilayer) = scatt_aux_ad % ext (ichan,ilayer) &
          & - scatt_aux_ad % ssa (ichan,ilayer) * ssa   &
          & / (ext * ext ) 
        scatt_aux_ad % ssa (ichan,ilayer) = scatt_aux_ad % ssa (ichan,ilayer) / ext 
     endif

     if (asm  > 0.0_JPRB ) then
        scatt_aux % asm (ichan,ilayer) = asm /  ssa 
        scatt_aux_ad % ssa (ichan,ilayer) = scatt_aux_ad % ssa (ichan,ilayer) &
          & - scatt_aux_ad % asm (ichan,ilayer) * asm        &
          & / (ssa * ssa )
        scatt_aux_ad     % asm (ichan,ilayer) = scatt_aux_ad  % asm (ichan,ilayer) / ssa 
     endif
     
   enddo
  enddo
  
!* Loops over channels, levels, hydrometeor types
  nlayer_loop2: do ilayer = nlevels, 1, -1

    ifreq_last = -1
    iprof_last = -1

    nchan_loop2: do ichan = 1, nchannels

      iprof = lprofiles   (ichan)    
      if (adk == adk_adjoint) then
        iprofad = iprof  
      else if (adk == adk_k) then
        iprofad = ichan  
      endif
      ifreq = frequencies (ichan)  

      if (scatt_aux % cfrac (iprof) > ccthres) then 

        if(iprof /= iprof_last) then

          cont(1)  = scatt_aux % rain (iprof,ilayer)
          cont(2)  = scatt_aux % sp  (iprof,ilayer) 
          cont(3)  = scatt_aux % clw (iprof,ilayer) 
          cont(4)  = scatt_aux % ciw (iprof,ilayer)
          cont(5)  = scatt_aux % totalice (iprof,ilayer)

          ntype_loop3: do itype = coef_scatt % nhydro, 1, -1

            !* Nearest index for Mie-table: LWC/IWC
            if (cont(itype) >= cont_min) then 
              wc(itype) = coef_scatt % scale_water * log10 (cont(itype)) - coef_scatt % offset_water 
            else if (cont(itype) < cont_min .and. cont(itype) > 0.0_JPRB) then
              wc(itype) = cont(itype) / cont_min
            else
              wc(itype) = 0.0_JPRB
            endif

            iwc(itype) = floor (wc(itype))
            if (iwc(itype) > coef_scatt % mwc - 1) then
              ! Prevent extrapolation 
              iwc(itype) = coef_scatt % mwc - 1
              wc(itype)  = coef_scatt % mwc
            endif

            !* Nearest index for Mie-table: T (w/o melting layer)
            temp = profiles (iprof) % t (ilayer) - offset_t(itype)

            itemp(itype) = anint (temp)
            if (itemp(itype) <                      1) itemp(itype) = 1
            if (itemp(itype) > coef_scatt % mtemp - 1) itemp(itype) = coef_scatt % mtemp - 1

          enddo ntype_loop3

        endif

          ntype_loop4: do itype = coef_scatt % nhydro, 1, -1

            if (iwc(itype) >= 1) then
              s_k = coef_scatt % ext (ifreq,itype,itemp(itype),iwc(itype)+1) &
                & - coef_scatt % ext (ifreq,itype,itemp(itype),iwc(itype))
              s_a = coef_scatt % ssa (ifreq,itype,itemp(itype),iwc(itype)+1) &
                & - coef_scatt % ssa (ifreq,itype,itemp(itype),iwc(itype))
              s_g = coef_scatt % asp (ifreq,itype,itemp(itype),iwc(itype)+1) &
                & - coef_scatt % asp (ifreq,itype,itemp(itype),iwc(itype))

              kp  = coef_scatt % ext (ifreq,itype,itemp(itype),iwc(itype)) + s_k * (wc(itype) - iwc(itype))
              ap  = coef_scatt % ssa (ifreq,itype,itemp(itype),iwc(itype)) + s_a * (wc(itype) - iwc(itype))
              gp  = coef_scatt % asp (ifreq,itype,itemp(itype),iwc(itype)) + s_g * (wc(itype) - iwc(itype))
           else
              ! For small water contents, interpolate linearly to zero 
              s_k   = coef_scatt % ext (ifreq,itype,itemp(itype),1) 
              s_a   = coef_scatt % ssa (ifreq,itype,itemp(itype),1) 
              s_g   = coef_scatt % asp (ifreq,itype,itemp(itype),1) 

              kp  = max( s_k * wc(itype), 1E-10_JPRB )              
              if (wc(itype) > 1E-10_JPRB) then
                ap = s_a * wc(itype)  
                gp = s_g * wc(itype)   
              else
                ap = 0.0_JPRB
                gp = 0.0_JPRB
              endif              

           endif

           kp_ad   = 0.0_JPRB
           ap_ad   = 0.0_JPRB
           gp_ad   = 0.0_JPRB
           wc_ad   = 0.0_JPRB

           kp_ad = kp_ad + ap * gp * scatt_aux_ad % asm (ichan,ilayer)
           ap_ad = ap_ad + kp * gp * scatt_aux_ad % asm (ichan,ilayer)
           gp_ad = gp_ad + kp * ap * scatt_aux_ad % asm (ichan,ilayer)

           kp_ad = kp_ad +      ap * scatt_aux_ad % ssa (ichan,ilayer) 
           ap_ad = ap_ad +      kp * scatt_aux_ad % ssa (ichan,ilayer)

           kp_ad = kp_ad +           scatt_aux_ad % ext (ichan,ilayer)

           if (iwc(itype) >= 1) then

             wc_ad = wc_ad + s_g * gp_ad
             gp_ad = 0.0_JPRB
             wc_ad = wc_ad + s_a * ap_ad
             ap_ad = 0.0_JPRB
             wc_ad = wc_ad + s_k * kp_ad
             kp_ad = 0.0_JPRB

           else

             if (wc(itype) > 1E-10_JPRB) then
               wc_ad = wc_ad + s_g * gp_ad
               gp_ad = 0.0_JPRB
               wc_ad = wc_ad + s_a * ap_ad
               ap_ad = 0.0_JPRB
             endif
             if( kp > 1E-10_JPRB) wc_ad = wc_ad + s_k * kp_ad
             kp_ad = 0.0_JPRB

           endif

           cont_ad(itype) = 0.0_JPRB
           if (cont(itype) >= cont_min) then 
             cont_ad(itype) = cont_ad(itype) + coef_scatt % scale_water * wc_ad / (zln10 * cont(itype)) 
           else if (cont(itype) < cont_min .and. cont(itype) > 0.0_JPRB) then
             cont_ad(itype) = cont_ad(itype) + wc_ad / cont_min 
  !          else
  !            cont_ad(itype) = cont_ad(itype) + 0.0_JPRB
           endif

           wc_ad = 0.0_JPRB

         enddo ntype_loop4

       if (scatt_aux    % rain (iprof,ilayer) > 0.0_JPRB) &
        & scatt_aux_ad % rain (iprofad,ilayer) = scatt_aux_ad % rain (iprofad,ilayer) + cont_ad(1)
       if (scatt_aux    % sp   (iprof,ilayer) > 0.0_JPRB) &
        & scatt_aux_ad % sp   (iprofad,ilayer) = scatt_aux_ad % sp   (iprofad,ilayer) + cont_ad(2)
       if (scatt_aux    % clw  (iprof,ilayer) > 0.0_JPRB) &
        & scatt_aux_ad % clw  (iprofad,ilayer) = scatt_aux_ad % clw  (iprofad,ilayer) + cont_ad(3)
       if (scatt_aux    % ciw  (iprof,ilayer) > 0.0_JPRB ) &
        & scatt_aux_ad % ciw  (iprofad,ilayer) = scatt_aux_ad % ciw  (iprofad,ilayer) + cont_ad(4)
       if (scatt_aux    % totalice  (iprof,ilayer) > 0.0_JPRB ) &
        & scatt_aux_ad % totalice  (iprofad,ilayer) = scatt_aux_ad % totalice (iprofad,ilayer) + cont_ad(5)

       ifreq_last=ifreq
       iprof_last=iprof

      endif 
    enddo nchan_loop2
  enddo nlayer_loop2

  if (lhook) call dr_hook('RTTOV_MIEPROC_AD',1_jpim,zhook_handle)

End subroutine rttov_mieproc_ad
