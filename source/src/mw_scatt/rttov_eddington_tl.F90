!
Subroutine rttov_eddington_tl ( &
     & nlevels,        &! in
     & nchannels,         &! in
     & nprofiles,         &! in
     & lprofiles,         &! in
     & angles,            &! in
     & profiles,          &! in
     & profiles_tl,       &! in
     & scatt_aux,         &! in
     & scatt_aux_tl,      &! in
     & cld_bt,            &! out
     & cld_bt_tl)          ! out 

  ! Description:
  ! to compute Eddington approximation to RT
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
  !  1.0       09/2002   Initial version     (P. Bauer, E. Moreau)
  !  1.1       05/2003   RTTOV7.3 compatible (F. Chevallier)
  !  1.2       03/2004   Added polarimetry   (R. Saunders)
  !  1.3       08/2004   Polarimetry fixes   (U. O'Keefe)
  !  1.4       11/2004   Clean-up            (P. Bauer)
  !  1.5       11/2007   RTTOV-9 version     (A. Geer)
  !  1.6       07/2008   Speed up and make consistent minimum SSA (A. Geer)
  !  1.7       07/2008   Clear sky speed-ups (A. Geer)
  !  1.8       06/2009   Bug fixed in index (T. Wilhelmsson)
  !  1.9       03/2010   Use mclayer rather than min_ssa (A. Geer)
  !
  ! Code Description:
  !   Language:           Fortran 90.
  !   Software Standards: "European Standards for Writing and
  !   Documenting Exchangeable Fortran 90 Code".
  !
  ! Declarations:
  ! Modules used:
  ! Imported Type Definitions:

  Use rttov_types, Only :      &
       & geometry_Type        ,&
       & profile_Type         ,&
       & profile_scatt_aux

  Use parkind1, Only : jpim     ,jprb
!INTF_OFF
  Use rttov_const, Only:       &
       & tcosmic, ccthres

  USE YOMHOOK, ONLY: LHOOK , DR_HOOK
!INTF_ON
  Implicit None

!* Subroutine arguments:
  Integer (Kind=jpim), Intent (in) :: nlevels            ! Number of NWP-levels
  Integer (Kind=jpim), Intent (in) :: nprofiles             ! Number of profiles
  Integer (Kind=jpim), Intent (in) :: nchannels             ! Number of channels*profiles=radiances
  Integer (Kind=jpim), Intent (in) :: lprofiles (nchannels) ! Profile indices
 
  Type (geometry_Type),       Intent (in)    :: angles       (nprofiles) ! Zenith angles
  Type (profile_Type),        Intent (in)    :: profiles     (nprofiles) ! Profiles 
  Type (profile_Type),        Intent (in)    :: profiles_tl  (nprofiles) ! Profiles 
  Type (profile_scatt_aux),   Intent (in)    :: scatt_aux                ! Auxiliary profile variables for RTTOV_SCATT
  Type (profile_scatt_aux),   Intent (in)    :: scatt_aux_tl             ! Auxiliary profile variables for RTTOV_SCATT
  Real (Kind=jprb),           Intent (out)   :: cld_bt       (nchannels) ! Radiances
  Real (Kind=jprb),           Intent (out)   :: cld_bt_tl    (nchannels) ! Radiances

!INTF_END
 
!* Local variables
  Real (Kind=jprb), Dimension (nchannels,nlevels) :: dp       ! D+ for boundary conditions
  Real (Kind=jprb), Dimension (nchannels,nlevels) :: dm       ! D- for boundary conditions
  Real (Kind=jprb), Dimension (nchannels,nlevels) :: j_up     ! Upward radiance source terms
  Real (Kind=jprb), Dimension (nchannels,nlevels) :: j_do     ! Downward radiance source terms
  Real (Kind=jprb), Dimension (nchannels,nlevels) :: dp_tl    ! D+ for boundary conditions
  Real (Kind=jprb), Dimension (nchannels,nlevels) :: dm_tl    ! D- for boundary conditions
  Real (Kind=jprb), Dimension (nchannels,nlevels) :: j_up_tl  ! Upward radiance source terms
  Real (Kind=jprb), Dimension (nchannels,nlevels) :: j_do_tl  ! Downward radiance source terms

  Real (Kind=jprb), Dimension (nchannels) :: irad_do, ftop       ! Downward radiances 
  Real (Kind=jprb), Dimension (nchannels) :: irad_up             ! Upward radiances   
  Real (Kind=jprb), Dimension (nchannels) :: irad_sfc            ! Inward radiances at surface
  Real (Kind=jprb), Dimension (nchannels) :: irad_space          ! Inward radiances from space
  Real (Kind=jprb), Dimension (nchannels) :: tau_t               ! Transmittances integrated over all levels
  Real (Kind=jprb), Dimension (nchannels) :: irad_do_tl, ftop_tl ! Downward radiances 
  Real (Kind=jprb), Dimension (nchannels) :: irad_up_tl          ! Upward radiances 
  Real (Kind=jprb), Dimension (nchannels) :: irad_sfc_tl         ! Inward radiances at surface
  Real (Kind=jprb), Dimension (nchannels) :: irad_space_tl       ! Inward radiances from space
  Real (Kind=jprb), Dimension (nchannels) :: tau_t_tl            ! Transmittances integrated over all levels

  Real (Kind=jprb), Dimension (nchannels,nlevels) :: j_part  ! Part of the linear in tau
  Real (Kind=jprb), Dimension (nchannels,nlevels) :: j_part_tl  

  Integer(Kind=jpim)  :: ilayer, jlayer, iprof, ichan        
  
  REAL(KIND=JPRB) :: ZHOOK_HANDLE

  !- End of header --------------------------------------------------------

#include "rttov_boundaryconditions_tl.interface"
#include "rttov_integratesource_tl.interface"

  IF (LHOOK) CALL DR_HOOK('RTTOV_EDDINGTON_TL',0_jpim,ZHOOK_HANDLE)

  j_up    (:,:) = 0.0_JPRB
  j_up_tl (:,:) = 0.0_JPRB
  j_do    (:,:) = 0.0_JPRB
  j_do_tl (:,:) = 0.0_JPRB

  !* Channels * Profiles      
  do ichan = 1, nchannels
    iprof = lprofiles (ichan)
     
    !* Top/bottom
    irad_sfc_tl (ichan) = scatt_aux_tl % ems_cld (ichan) * profiles    (iprof) % skin % t &
                       & + scatt_aux    % ems_cld (ichan) * profiles_tl (iprof) % skin % t
    irad_sfc    (ichan) = scatt_aux    % ems_cld (ichan) * profiles    (iprof) % skin % t

    irad_space_tl (ichan) = 0.0_JPRB
    irad_space    (ichan) = tcosmic

    if (scatt_aux % cfrac (iprof) > ccthres) then

      !* Clear-sky source terms (but only actually used if SSA is neglible)
      do ilayer = 1, nlevels
        if (ilayer < scatt_aux % mclayer(ichan)) then
         if ( scatt_aux % delta (ichan,ilayer) > 1E-8 ) then

           j_part (ichan,ilayer) = &
            &   (scatt_aux % b0 (iprof,ilayer) - scatt_aux % bn (iprof,ilayer)) &
            & * (1.0_JPRB - scatt_aux % tau (ichan,ilayer))                     &
            & / scatt_aux % delta (ichan,ilayer)
           j_part_tl(ichan,ilayer) = &
            &   (scatt_aux_tl % b0 (iprof,ilayer) - scatt_aux_tl % bn (iprof,ilayer)) & 
            & * (1.0_JPRB - scatt_aux % tau (ichan,ilayer))                           &
            & / scatt_aux % delta (ichan,ilayer)                                      &
            & - scatt_aux_tl % tau (ichan,ilayer) / scatt_aux % delta (ichan,ilayer)  &    
            & * (scatt_aux % b0 (iprof,ilayer) - scatt_aux % bn (iprof,ilayer))       &  
            & - scatt_aux_tl % delta (ichan,ilayer) * j_part (ichan,ilayer)           &  
            & / scatt_aux % delta (ichan,ilayer)

           j_up (ichan,ilayer)    = scatt_aux % bn (iprof,ilayer)              &
            & - scatt_aux % b0 (iprof,ilayer) * scatt_aux % tau (ichan,ilayer) &
            & + j_part(ichan,ilayer)
           j_up_tl (ichan,ilayer) = scatt_aux_tl % bn (iprof,ilayer)              &
            & - scatt_aux_tl % b0 (iprof,ilayer) * scatt_aux % tau (ichan,ilayer) &
            & - scatt_aux % b0 (iprof,ilayer) * scatt_aux_tl % tau (ichan,ilayer) &
            & + j_part_tl(ichan,ilayer)

           j_do (ichan,ilayer)    = scatt_aux % b0 (iprof,ilayer)              &
            & - scatt_aux % bn (iprof,ilayer) * scatt_aux % tau (ichan,ilayer) &
            & - j_part(ichan,ilayer)
           j_do_tl (ichan,ilayer) = scatt_aux_tl % b0 (iprof,ilayer)              &
            & - scatt_aux_tl % bn (iprof,ilayer) * scatt_aux % tau (ichan,ilayer) &
            & - scatt_aux % bn (iprof,ilayer) * scatt_aux_tl % tau (ichan,ilayer) &
            & - j_part_tl(ichan,ilayer)

         else

           ! Linear in tau formula is computationally unstable for small delta - avoid
           j_up_tl (ichan,ilayer) = scatt_aux_tl % b0 (iprof,ilayer) &
             & * (1.0_JPRB - scatt_aux % tau (ichan,ilayer)) &
             & - scatt_aux    % b0 (iprof,ilayer) * scatt_aux_tl % tau (ichan,ilayer)
           j_up    (ichan,ilayer) = scatt_aux    % b0 (iprof,ilayer) &
             & * (1.0_JPRB - scatt_aux % tau (ichan,ilayer)) 

           j_do_tl (ichan,ilayer) = j_up_tl (ichan,ilayer)
           j_do    (ichan,ilayer) = j_up    (ichan,ilayer)

         endif
       endif
      enddo

      !* Downward radiance at cloud top
      irad_do_tl (ichan) = irad_space_tl (ichan)
      irad_do    (ichan) = irad_space    (ichan)
  
      Do ilayer = 1, scatt_aux % mclayer (ichan) - 1
        irad_do_tl (ichan) = irad_do_tl (ichan) * scatt_aux % tau (ichan,ilayer) &
         & + irad_do (ichan) * scatt_aux_tl % tau (ichan,ilayer) + j_do_tl (ichan,ilayer)
        irad_do    (ichan) = irad_do (ichan) * scatt_aux % tau (ichan,ilayer) &
         & + j_do (ichan,ilayer)
      End do
     
      ftop    (ichan) = irad_do    (ichan)
      ftop_tl (ichan) = irad_do_tl (ichan)
    endif
  End do 

  !* Get D+, D- from boundary conditions
  Call rttov_boundaryconditions_tl (&
    & nlevels,          &! in
    & nchannels,        &! in
    & nprofiles,        &! in
    & lprofiles,        &! in
    & scatt_aux,        &! in
    & scatt_aux_tl,     &! in
    & profiles,         &! in
    & profiles_tl ,     &! in
    & ftop,             &! in
    & ftop_tl,          &! in
    & dp,               &! out
    & dp_tl,            &! out
    & dm,               &! out
    & dm_tl)             ! out 

  !* Integrate radiance source terms 
  Call rttov_integratesource_tl (&
    & nlevels,          &! in
    & nchannels,        &! in
    & nprofiles,        &! in
    & lprofiles,        &! in
    & angles,           &! in
    & scatt_aux,        &! in
    & scatt_aux_tl,     &! in
    & dp,               &! in
    & dp_tl,            &! in
    & dm,               &! in
    & dm_tl,            &! in
    & j_do,             &! inout
    & j_do_tl,          &! inout
    & j_up,             &! inout
    & j_up_tl)           ! inout 

  !* Integrate downward radiances/transmittance
  irad_do_tl (:) = 0.0_JPRB
  irad_do    (:) = irad_space (:)

  irad_up_tl (:) = irad_sfc_tl (:)
  irad_up    (:) = irad_sfc    (:)

  tau_t_tl   (:) = 0.0_JPRB
  tau_t      (:) = 1.0_JPRB

  Do ilayer = 1, nlevels
     jlayer = nlevels + 1 - ilayer
  
     irad_do_tl (:) = irad_do_tl (:) * scatt_aux    % tau (:,ilayer) &
                  & + irad_do    (:) * scatt_aux_tl % tau (:,ilayer) + j_do_tl (:,ilayer)
     irad_do    (:) = irad_do    (:) * scatt_aux    % tau (:,ilayer) + j_do    (:,ilayer)

     irad_up_tl (:) = irad_up_tl (:) * scatt_aux    % tau (:,jlayer) &
                  & + irad_up    (:) * scatt_aux_tl % tau (:,jlayer) + j_up_tl (:,jlayer)
     irad_up    (:) = irad_up    (:) * scatt_aux    % tau (:,jlayer) + j_up    (:,jlayer)

     tau_t_tl   (:) = tau_t_tl   (:) * scatt_aux    % tau (:,jlayer) &
                  & + tau_t      (:) * scatt_aux_tl % tau (:,jlayer)
     tau_t      (:) = tau_t      (:) * scatt_aux    % tau (:,jlayer)
  Enddo
                           
  cld_bt_tl (:) = irad_up_tl (:) + (scatt_aux_tl % ref_cld (:) * irad_do    (:) &
                               & +  scatt_aux    % ref_cld (:) * irad_do_tl (:)) * tau_t    (:) &
                               & + (scatt_aux    % ref_cld (:) * irad_do    (:)) * tau_t_tl (:)
  cld_bt    (:) = irad_up    (:) +  scatt_aux    % ref_cld (:) * irad_do    (:)  * tau_t    (:)
  
  IF (LHOOK) CALL DR_HOOK('RTTOV_EDDINGTON_TL',1_jpim,ZHOOK_HANDLE)

End Subroutine rttov_eddington_tl
