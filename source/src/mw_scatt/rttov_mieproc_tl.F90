!
Subroutine rttov_mieproc_tl (&        
     & nlevels,           &! in
     & nchannels,         &! in
     & nprofiles,         &! in
     & frequencies,       &! in
     & lprofiles,         &! in
     & profiles,          &! in
     & coef_scatt,        &! in
     & scatt_aux,         &! inout
     & scatt_aux_tl)       ! inout 
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
  !  1.2       11/2004   Clean-up            (P. Bauer)
  !  1.3       10/2006   Introduce interpolation to zero for scatt coeffs 
  !                      for small LWC (U.O'Keeffe). Cleaner version (A. Geer)
  !  1.4       11/2007   RTTOV9 version      (A. Geer)
  !  1.5       06/2008   Performance enhancements (D. Salmond)
  !  1.6       06/2008   Fix minor bug in small LWC interpolation (A. Geer)
  !  1.7       07/2008   Clear sky speed-ups (A. Geer)
  !  1.8       03/2009   Prevent extrapolation beyond table bounds (A. Geer)
  !  1.9       03/2010   Optimisation (A. Geer)
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

  Use rttov_types, Only :      &
       & profile_scatt_aux    ,&
       & profile_Type   ,&
       & rttov_scatt_coef

  Use parkind1, Only : jpim
!INTF_OFF
  Use rttov_const, Only : ccthres

  Use parkind1, Only : jprb

  USE YOMHOOK, ONLY: LHOOK , DR_HOOK
!INTF_ON
  Implicit None

!* Subroutine arguments:
  Integer (Kind=jpim), Intent (in) :: nlevels                 ! Number of  levels
  Integer (Kind=jpim), Intent (in) :: nchannels               ! Number of channels*profiles=radiances
  Integer (Kind=jpim), Intent (in) :: nprofiles               ! Number of profiles
  Integer (Kind=jpim), Intent (in) :: frequencies (nchannels) ! Frequency indices
  Integer (Kind=jpim), Intent (in) :: lprofiles   (nchannels) ! Profile indices

  Type (profile_Type),      Intent (in)    :: profiles    (nprofiles) ! profiles 
  Type (rttov_scatt_coef),  Intent (in)    :: coef_scatt              ! RTTOV_SCATT Coefficients
  Type (profile_scatt_aux), Intent (inout) :: scatt_aux               ! Auxiliary profile variables
  Type (profile_scatt_aux), Intent (inout) :: scatt_aux_tl            ! Auxiliary profile variables

!INTF_END

!* Local variables:
  Integer (Kind=jpim) :: iwc(coef_scatt % nhydro), itemp(coef_scatt % nhydro), itype, ichan, ifreq, iprof, ilayer
  Integer (Kind=jpim) :: ifreq_last, iprof_last
  Real    (Kind=jprb) :: wc(coef_scatt % nhydro)   , temp   ,  kp   , ap   , gp   , s_k   , s_a   , s_g   , zln10     
  Real    (Kind=jprb) :: kpp,app,gpp,kpp_tl,app_tl,gpp_tl
  Real    (Kind=jprb) :: wc_tl(coef_scatt % nhydro), kp_tl, ap_tl, gp_tl
  Real    (Kind=jprb) :: cont(coef_scatt % nhydro), cont_tl(coef_scatt % nhydro), cont_min, offset_t(coef_scatt % nhydro)
  
  REAL(KIND=JPRB) :: ZHOOK_HANDLE

  !- End of header --------------------------------------------------------
  
  if (lhook) call dr_hook('RTTOV_MIEPROC_TL',0_jpim,zhook_handle)

  zln10 = log (10.0_JPRB)
  cont_min = 10.0_JPRB ** ( (1.0_JPRB + coef_scatt % offset_water) / coef_scatt % scale_water )

  offset_t(1) = coef_scatt % offset_temp_rain
  offset_t(2) = coef_scatt % offset_temp_sp
  offset_t(3) = coef_scatt % offset_temp_liq
  offset_t(4) = coef_scatt % offset_temp_ice
  offset_t(5) = coef_scatt % offset_temp_totalice

!* Loops over channels, levels, hydrometeor types
  nlayer_loop: do ilayer = 1, nlevels

    ifreq_last = -1
    iprof_last = -1

    nchan_loop: do ichan = 1, nchannels

      iprof = lprofiles   (ichan)     
      ifreq = frequencies (ichan)  

      if( scatt_aux % cfrac(iprof) > ccthres ) then   

        if(iprof /= iprof_last ) then

          cont(1)     = scatt_aux % rain (iprof,ilayer)
          cont_tl(1)  = scatt_aux_tl % rain (iprof,ilayer)
          cont(2)     = scatt_aux % sp  (iprof,ilayer) 
          cont_tl(2)  = scatt_aux_tl % sp (iprof,ilayer)
          cont(3)     = scatt_aux % clw (iprof,ilayer) 
          cont_tl(3)  = scatt_aux_tl % clw (iprof,ilayer)
          cont(4)     = scatt_aux % ciw (iprof,ilayer)
          cont_tl(4)  = scatt_aux_tl % ciw (iprof,ilayer)
          cont(5)     = scatt_aux % totalice (iprof,ilayer)
          cont_tl(5)  = scatt_aux_tl % totalice (iprof,ilayer)

          ntype_loop1: do itype = 1, coef_scatt % nhydro

            !* nearest index for Mie-table: LWC/IWC
            if (cont(itype) >= cont_min) then 
              wc(itype)    = coef_scatt % scale_water * log10 (cont(itype)) - coef_scatt % offset_water
              wc_tl(itype) = coef_scatt % scale_water * cont_tl(itype) / (zln10 * cont(itype))   
            else if (cont(itype) < cont_min .and. cont(itype) > 0.0_JPRB) then
              wc(itype)    = cont(itype) / cont_min
              wc_tl(itype) = cont_tl(itype) / cont_min
            else
              wc(itype)    = 0.0_JPRB
              wc_tl(itype) = 0.0_JPRB
            endif

            iwc(itype) = floor (wc(itype))
            if (iwc(itype) > coef_scatt % mwc - 1) then
              ! Prevent extrapolation 
              iwc(itype) = coef_scatt % mwc - 1
              wc(itype)  = coef_scatt % mwc
            endif

            !* nearest index for Mie-table: T (w/o melting layer)
            temp = profiles (iprof) % t (ilayer) - offset_t(itype)

            itemp(itype) = anint (temp)
            if (itemp(itype) <                      1) itemp(itype) = 1
            if (itemp(itype) > coef_scatt % mtemp - 1) itemp(itype) = coef_scatt % mtemp - 1

          enddo ntype_loop1

        endif

        if(ifreq /= ifreq_last .or. iprof /= iprof_last ) then

          kpp=0.0_JPRB
          app=0.0_JPRB
          gpp=0.0_JPRB
          kpp_tl=0.0_JPRB
          app_tl=0.0_JPRB
          gpp_tl=0.0_JPRB

          ntype_loop: do itype = 1, coef_scatt % nhydro 

            if (iwc(itype) >= 1) then

              s_k   = coef_scatt % ext (ifreq,itype,itemp(itype),iwc(itype)+1) &
                  & - coef_scatt % ext (ifreq,itype,itemp(itype),iwc(itype))
              s_a   = coef_scatt % ssa (ifreq,itype,itemp(itype),iwc(itype)+1) &
                  & - coef_scatt % ssa (ifreq,itype,itemp(itype),iwc(itype))
              s_g   = coef_scatt % asp (ifreq,itype,itemp(itype),iwc(itype)+1) &
                  & - coef_scatt % asp (ifreq,itype,itemp(itype),iwc(itype))

              kp    = coef_scatt % ext (ifreq,itype,itemp(itype),iwc(itype)) + s_k * (wc(itype) - iwc(itype))
              kp_tl = s_k * wc_tl(itype)

              ap    = coef_scatt % ssa (ifreq,itype,itemp(itype),iwc(itype)) + s_a * (wc(itype) - iwc(itype))
              ap_tl = s_a * wc_tl(itype)

              gp    = coef_scatt % asp (ifreq,itype,itemp(itype),iwc(itype)) + s_g * (wc(itype) - iwc(itype))
              gp_tl = s_g * wc_tl(itype)       
                                   
            else

              ! For small water contents, interpolate linearly to zero 
              s_k   = coef_scatt % ext (ifreq,itype,itemp(itype),1) 
              s_a   = coef_scatt % ssa (ifreq,itype,itemp(itype),1) 
              s_g   = coef_scatt % asp (ifreq,itype,itemp(itype),1) 

              kp  = max( s_k * wc(itype), 1E-10_JPRB )              
              if (kp > 1E-10_JPRB) then
                kp_tl = s_k * wc_tl(itype) 
              else 
                kp_tl = 0.0_JPRB  
              endif

              if (wc(itype) > 1E-10_JPRB) then
                ap    = s_a * wc(itype)
                ap_tl = s_a * wc_tl(itype)
                gp    = s_g * wc(itype) 
                gp_tl = s_g * wc_tl(itype)   
              else
                ap    = 0.0_JPRB 
                ap_tl = 0.0_JPRB
                gp    = 0.0_JPRB 
                gp_tl = 0.0_JPRB
              endif

            endif
            kpp=kpp+kp
            app=app+kp*ap
            gpp=gpp+kp*ap*gp
            kpp_tl=kpp_tl+kp_tl
            app_tl=app_tl+kp_tl*ap+kp*ap_tl
            gpp_tl=gpp_tl+kp_tl*ap*gp+kp*ap_tl*gp+kp*ap*gp_tl

          enddo ntype_loop
        endif

        ifreq_last=ifreq
        iprof_last=iprof
        scatt_aux    % ext (ichan,ilayer) = scatt_aux    % ext (ichan,ilayer) + kpp
        scatt_aux_tl % ext (ichan,ilayer) = scatt_aux_tl % ext (ichan,ilayer) + kpp_tl

        scatt_aux    % ssa (ichan,ilayer) = scatt_aux    % ssa (ichan,ilayer) + app
        scatt_aux_tl % ssa (ichan,ilayer) = scatt_aux_tl % ssa (ichan,ilayer) + app_tl

        scatt_aux    % asm (ichan,ilayer) = scatt_aux    % asm (ichan,ilayer) + gpp
        scatt_aux_tl % asm (ichan,ilayer) = scatt_aux_tl % asm (ichan,ilayer) + gpp_tl

      endif
    enddo nchan_loop
  enddo nlayer_loop

  do ilayer = 1, nlevels
  do ichan = 1, nchannels
    if (scatt_aux    % asm (ichan,ilayer) > 0.0_JPRB) then
    scatt_aux_tl % asm (ichan,ilayer) = scatt_aux_tl % asm (ichan,ilayer) / scatt_aux % ssa (ichan,ilayer) &
                        & - scatt_aux_tl % ssa (ichan,ilayer) * scatt_aux % asm (ichan,ilayer) &
                        & / scatt_aux    % ssa (ichan,ilayer) / scatt_aux % ssa (ichan,ilayer) 
    scatt_aux    % asm (ichan,ilayer) = scatt_aux    % asm (ichan,ilayer) / scatt_aux % ssa (ichan,ilayer)
     endif
    if (scatt_aux    % ssa (ichan,ilayer) > 0.0_JPRB) then
    scatt_aux_tl % ssa (ichan,ilayer) = (scatt_aux_tl % ssa (ichan,ilayer) * scatt_aux    % ext (ichan,ilayer) &
                        & -  scatt_aux    % ssa (ichan,ilayer) * scatt_aux_tl % ext (ichan,ilayer)) &
                        & / (scatt_aux    % ext (ichan,ilayer) * scatt_aux    % ext (ichan,ilayer)) 
    scatt_aux    % ssa (ichan,ilayer) =  scatt_aux    % ssa (ichan,ilayer) / scatt_aux    % ext (ichan,ilayer)
     endif
    if (scatt_aux % ext (ichan,ilayer) >= 20.0_JPRB) then
         scatt_aux_tl % ext (ichan,ilayer) =  0.0_JPRB
         scatt_aux    % ext (ichan,ilayer) = 20.0_JPRB
     endif
  enddo
  enddo

  if (lhook) call dr_hook('RTTOV_MIEPROC_TL',1_jpim,zhook_handle)

End subroutine rttov_mieproc_tl
