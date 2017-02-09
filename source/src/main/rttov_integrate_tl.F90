Subroutine rttov_integrate_tl( &
   addcosmic, opts, maxnstreams, chanprof,      &! in
   emissivity, emissivity_tl,    &! in
   reflectance, reflectance_tl,    &! in
   refl_norm,        &! in
   thermrefl, thermrefl_tl,  &! in
   do_lambertian,    &! in
   thermal,          &! in
   dothermal,        &! in
   solar,            &! in
   dosolar,          &! in
   solar_spectrum,   &! in
   transmission_aux, transmission_aux_tl,  &! in
   transmission_scatt_ir, transmission_scatt_ir_tl,&
   profiles, profiles_tl, &! in
   profiles_dry, profiles_dry_tl, &! in
   aux_prof, aux_prof_tl,      &! in
   coef,     &! in
   raytracing, raytracing_tl,    &! in
   ircld, ircld_tl,         &! in
   rad , &! in
   auxrad ,          &! in
   auxrad_stream, auxrad_stream_tl, &! in
   rad_tl           ) ! inout
!
! Description:
! To perform TL of integration of radiative transfer equation
! in rttov suite, calculating radiances and brightness temperature.
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
! Eyre J.R. 1991 A fast radiative transfer model for satellite sounding
! systems.  ECMWF Research Dept. Tech. Memo. 176 (available from the
! librarian at ECMWF).
!
! Saunders R.W., M. Matricardi and P. Brunel 1999 An Improved Fast Radiative
! Transfer Model for Assimilation of Satellite Radiance Observations.
! QJRMS, 125, 1407-1425.
!
! Matricardi, M. 2003 RTIASI-4, a new version of the ECMWF fast radiative
! transfer model for the infrared atmospheric sounding interferometer.
! ECMWF Research Dept. Tech. Memo. 425 (available from the librarian at ECMWF)
!
! Current Code Owner: SAF NWP
!
! History:
! Version   Date     Comment
! -------   ----     -------
!          25/06/91.    Original code.  J.R.EYRE    *ECMWF*
!          21/08/00.    Emissivity and reflectivity handled separately. Steve English
!          31/01/01.    More cloud computations. F. Chevallier
!          23/03/01     New coef. format, new channel numbers (P. Brunel)
!          31/01/01.    More cloud computations. F. Chevallier
!          28/09/01     Cosmic background temp added G.Deblonde
!          18/01/2002   Thread safe (D.Salmond)
!  1.0     01/12/2002   New F90 code with structures (P Brunel A Smith)
!  1.1     02/01/2003   Added comments (R Saunders)
!  1.2     06/05/2003   Init rad%downcld to 0 in section 1 (P  Brunel)
!  1.3     26/09/2003   Modified to allow for multiple polarisations (S English)
!  1.4     28/02/2005   Improved vectorisation (D Dent)
!  1.5     29/03/2005   Add end of header comment (J. Cameron)
!  1.6     01/06/2005   Marco Matricardi (ECMWF):
!             --        IASI capability added.
!             --        Linear in tau approximation for RT equation introduced.
!             --        Solar radiation introduced for IASI and AIRS.
!  1.7     03/01/2007   Corrected bugs in'old code' R Saunders
!  1.8     08/02/2007   Removed polarisation index (R Saunders)
!  1.9     11/03/2007   Reintroduced overcast radiance (R Saunders)
!  1.10    15/11/2007   Changed skin to surfair for overcast radiances(RSaunders)
!  1.11    27/11/2007   Optimised for NEC/IBM (D Salmond)
!  1.12    21/12/2007   Added polarimetric option (R. Saunders)
!  1.13    15/07/2009   User defined ToA. Layers distinct from levels. Top-layer
!                       brought into layer looping to shorten code (P.Rayer)
!  1.14    03/11/2009   Transmittances on levels (A Geer)
!  1.15    02/12/2009   Introduced principal component capability. Pathsat, Pathsun and
!                       related quantities are now layer arrays (Marco Matricardi).
!  1.16    05/07/2010   Remove addsolar flag from profiles structure (J Hocking)
!  1.17    14/10/2010   Remove rt8_mode (J Hocking)
!  1.18    14/12/2010   Use traj0_sta%solar array to flag channels for which solar calculations
!                       should be performed (J Hocking)
!  2.0     14/12/2011   Re-written (D Rundle)
!
! 2010/03 Code cleaning Pascal Brunel, Philippe Marguinaud
!
!
!
! Code Description:
!   Language:           Fortran 90.
!   Software Standards: "European Standards for Writing and
!     Documenting Exchangeable Fortran 90 Code".
!

  use rttov_types, Only : rttov_chanprof, rttov_coef, rttov_options, profile_type, profile_aux, &
       transmission_type_aux, transmission_scatt_ir_type, radiance_type, ircld_type, raytracing_type, &
       radiance_aux, rttov_emissivity, rttov_reflectance

  use parkind1, Only : jpim, jprb, jplm
!INTF_OFF

  use rttov_const, Only : sensor_id_po, pi_r, deg2rad, &
       overcast_albedo_wvn, overcast_albedo1, overcast_albedo2, &
       gas_unit_compatibility

  use yomhook, Only : LHOOK, DR_HOOK
!INTF_ON

  Implicit None

!subroutine arguments:

  logical(jplm),   intent(in)                      :: addcosmic   ! switch for adding cosmic background
  type(rttov_options),  intent(in)                 :: opts        ! options structure
  integer(jpim),   intent(in)                      :: maxnstreams                !
  type(rttov_chanprof), intent(in)                 :: chanprof(:)     ! Channel indices
  type(profile_type),   intent(in)                 :: profiles(:)     ! Profiles
  type(profile_type),   intent(in)                 :: profiles_dry(size(profiles)) ! Gas profiles (ppmv dry)
  type(rttov_emissivity),  intent(in), optional    :: emissivity(size(chanprof))   ! surface emissivity
  type(rttov_reflectance), intent(in), optional    :: reflectance(size(chanprof))  ! surface solar reflectance
  real(jprb),      intent(in)                      :: refl_norm(size(chanprof))
  real(jprb),      intent(in)                      :: thermrefl(size(chanprof))    ! surface thermal reflectance
  logical(jplm),   intent(in)                      :: do_lambertian(:)
  logical(jplm),   intent(in)                      :: thermal(size(chanprof))
  logical(jplm),   intent(in)                      :: dothermal
  logical(jplm),   intent(in)                      :: solar(size(chanprof))
  logical(jplm),   intent(in)                      :: dosolar
  real(jprb),      intent(in)                      :: solar_spectrum(size(chanprof))
  type(ircld_type),     intent(in)                 :: ircld
  type(raytracing_type),intent(in)                 :: raytracing
  type(transmission_type_aux), intent(in)          :: transmission_aux            ! transmittances and single-layer od
  type(transmission_scatt_ir_type), intent(in)     :: transmission_scatt_ir
  type(profile_aux) ,   intent(in)                 :: aux_prof ! auxillary profiles info.
  type(rttov_coef),     intent(in)                 :: coef
  type(radiance_aux),   intent(in)                 :: auxrad_stream
  type(radiance_type),  intent(in)                 :: rad    ! radiances (mw/cm-1/ster/sq.m) and BTs
  type(radiance_aux),   intent(in)                 :: auxrad ! auxillary radiances

  Type(rttov_emissivity),  Intent(inout), Optional :: emissivity_tl(size(chanprof))
  Type(rttov_reflectance), Intent(inout), Optional :: reflectance_tl(size(chanprof))
  Real(jprb),              Intent(in)              :: thermrefl_tl(size(chanprof))
  Type(profile_Type),      Intent(in)              :: profiles_tl(size(profiles))
  Type(profile_Type),      Intent(in)              :: profiles_dry_tl(size(profiles))
  Type(ircld_type),        Intent(in)              :: ircld_tl
  Type(raytracing_type),   Intent(in)              :: raytracing_tl
  Type(transmission_Type_aux), Intent(in)          :: transmission_aux_tl
  type(transmission_scatt_ir_type) ,intent(in)     :: transmission_scatt_ir_tl
  Type(profile_aux) ,  Intent(in)                  :: aux_prof_tl
  Type(radiance_Type), Intent(inout)               :: rad_tl ! in because of mem allocation
  Type(radiance_aux),  Intent(inout)               :: auxrad_stream_tl

!INTF_END

#include "rttov_calcrad_tl.interface"

!local variables: 
  integer(jpim) :: i, lev, ist, isti, lay, nchannels, nlayers, nlevels, narray(3), iprof, prof, chan, nstreams ! counter variables
  logical(jplm) :: keyradonly ! flag to indicate only calculate primary radiance outputs
  real(jprb)    :: cfraction(size(chanprof)), cfraction_tl(size(chanprof)), pfraction(size(chanprof)) ! cloud fraction
  logical(jplm) :: sateqsun(profiles(1)%nlayers,size(profiles(:))) ! True where the solar zenith angle equal to observation angle
  integer(jpim) :: pol_id(size(chanprof))       ! polarisation index

  Real(jprb)    :: rad_air_tl(profiles(1)%nlevels, size(chanprof)), rad_surfair_tl(size(chanprof)), rad_skin_tl(size(chanprof))

  integer(jpim) :: iv2lay(size(chanprof)), iv2lev(size(chanprof)), iv3lay(size(chanprof)), iv3lev(size(chanprof))
  real(jprb)    :: refl, refl_norm_scat

  REAL(JPRB)    :: ZHOOK_HANDLE

!- End of header ------------------------------------------------------

if (LHOOK) CALL DR_HOOK('RTTOV_INTEGRATE_TL',0_jpim,ZHOOK_HANDLE)

#define tau_surf_p_r_tl transmission_aux_tl%thermal_path1%Tau_surf_p_r(ist,i)
#define tau_surf_p_tl transmission_aux_tl%thermal_path1%Tau_surf_p(ist,i)
#define tau_surf_tl transmission_aux_tl%thermal_path1%Tau_surf(ist,i)

#define tau_layer_p_r_tl transmission_aux_tl%thermal_path1%Tau_level_p_r(lay,ist,i)
#define tau_layer_p_tl transmission_aux_tl%thermal_path1%Tau_level_p(lay,ist,i)
#define tau_layer_tl transmission_aux_tl%thermal_path1%Tau_level(lay,ist,i)

#define tau_level_p_r_tl transmission_aux_tl%thermal_path1%Tau_level_p_r(lay+1,ist,i)
#define tau_level_p_tl transmission_aux_tl%thermal_path1%Tau_level_p(lay+1,ist,i)
#define tau_level_tl transmission_aux_tl%thermal_path1%Tau_level(lay+1,ist,i)

#define daux_lay_tl (rad_air_tl(lay+1,i) - rad_air_tl(lay,i))
#define dtau_lay_tl (tau_layer_tl - tau_level_tl)
#define od_singlelayer_tl transmission_aux_tl%thermal_path1%Od_singlelayer(isti,lay,i)

#define tau_surf_p_r transmission_aux%thermal_path1%Tau_surf_p_r(ist,i)
#define tau_surf_p transmission_aux%thermal_path1%Tau_surf_p(ist,i)
#define tau_surf_r transmission_aux%thermal_path1%Tau_surf_r(ist,i)
#define tau_surf transmission_aux%thermal_path1%Tau_surf(ist,i)

#define tau_layer_p_r transmission_aux%thermal_path1%Tau_level_p_r(lay,ist,i)
#define tau_layer_p transmission_aux%thermal_path1%Tau_level_p(lay,ist,i)
#define tau_layer_r transmission_aux%thermal_path1%Tau_level_r(lay,ist,i)
#define tau_layer transmission_aux%thermal_path1%Tau_level(lay,ist,i)

#define tau_level_p_r transmission_aux%thermal_path1%Tau_level_p_r(lay+1,ist,i)
#define tau_level_p transmission_aux%thermal_path1%Tau_level_p(lay+1,ist,i)
#define tau_level_r transmission_aux%thermal_path1%Tau_level_r(lay+1,ist,i)
#define tau_level transmission_aux%thermal_path1%Tau_level(lay+1,ist,i)

#define daux_lay (auxrad%air(lay+1,i) - auxrad%air(lay,i))
#define dtau_lay (tau_layer - tau_level)
#define od_singlelayer_r transmission_aux%thermal_path1%Od_singlelayer_r(isti,lay,i)
#define od_singlelayer transmission_aux%thermal_path1%Od_singlelayer(isti,lay,i)

!---------------------------
!0. Initialise useful variables
!---------------------------
nchannels = size(chanprof)
nlayers = profiles(1)%nlayers
nlevels = nlayers + 1

! DAR: narray is a more compact way of passing nchannels, nlevels and maxnstreams to the internal subroutines
! Again module variables would be a more elegant way of solving this problem and would results in less duplicated code
narray(1) = nlevels; narray(2) = maxnstreams; narray(3) = nchannels

keyradonly = opts%rt_ir%addaerosl .or. opts%rt_ir%pc%addpc

do i = 1, nchannels
  prof = chanprof(i)%prof
  cfraction(i) = aux_prof%s(prof)%cfraction
  cfraction_tl(i) = aux_prof_tl%s(prof)%cfraction
  pfraction(i) = aux_prof%s(prof)%pfraction_surf
enddo

if (dosolar) then
  do iprof = 1, size(profiles)
    do lay = 1, nlayers
      if(raytracing%pathsat(lay,iprof) == raytracing%pathsun(lay,iprof)) then
        sateqsun(lay,iprof) = .true.
      else
        sateqsun(lay,iprof) = .false.
      endif
    enddo
  enddo
endif

do i = 1, nchannels
  prof = chanprof(i)%prof
! case-1: surf lies above lev=nlevels
  iv3lev(i) = aux_prof%s(prof)%nearestlev_surf - 1   ! lowest lev above surf
! case-2: surf lies below lev=nlevels
  if (pfraction(i) < 0.0_JPRB) iv3lev(i) = iv3lev(i) + 1  ! iv3lev=iv2lev=lowest lev above surf

  iv2lev(i) = aux_prof%s(prof)%nearestlev_surf       ! highest lev below surf
  iv2lay(i) = iv2lev(i) - 1                          ! same layer as that numbered by  iv2 in RTTOV-9
  iv3lay(i) = iv3lev(i) - 1                          ! same layer as that numbered by  iv3 in RTTOV-9
enddo

if(coef%id_sensor == sensor_id_po) then
  do i = 1, nchannels
    chan = chanprof(i)%chan
    pol_id(i) = coef%fastem_polar(chan) + 1_jpim
  enddo
else
  pol_id(:) = 0_jpim
endif

if (keyradonly) rad_tl%overcast = 0.0_jprb
rad_tl%clear = 0.0_jprb
auxrad_stream_tl%cloudy = 0.0_jprb
if (dosolar) then
  auxrad_stream_tl%up_solar = 0.0_jprb
  auxrad_stream_tl%meanrad_up_solar = 0.0_jprb
  auxrad_stream_tl%down_solar = 0.0_jprb
  auxrad_stream_tl%meanrad_down_solar = 0.0_jprb
endif

!----------------------------
!1. calculate layer radiances
!----------------------------
if (dothermal) Call rttov_calcrad_tl(chanprof, &! in
                                     profiles, profiles_tl,  &! in
                                     coef,         &! in
                                     auxrad, & !in
                                     rad_skin_tl, rad_surfair_tl, rad_air_tl) ! out

!-------------------------------------
!2. calculate atmospheric contribution
!-------------------------------------
if (dothermal) call calc_atmospheric_radiance_tl(thermal, transmission_aux, transmission_aux_tl, auxrad, ircld, &
                                                 chanprof, narray, rad_air_tl, auxrad_stream_tl)

!---Scattering of the solar beam--------------------------------------------------
if((opts%rt_ir%addaerosl .OR. opts%rt_ir%addclouds) .AND. dosolar) &
    call  solar_scattering_air_tl(transmission_aux, transmission_aux_tl, ircld, raytracing, raytracing_tl, &
                                       chanprof, narray, transmission_scatt_ir_tl, &
                                       solar, solar_spectrum, reflectance%refl_out, sateqsun, &
                                       auxrad_stream, auxrad_stream_tl)

!-------------------------------------------------------------------------------
!2a calculate near-surface layer contribution
!-------------------------------------------------------------------------------
if (dothermal) call calc_near_surf_contribution_tl(thermal, transmission_aux, transmission_aux_tl, auxrad, &
                                                   ircld, chanprof, narray, iv3lay, iv3lev, pol_id,& !in
                                                   rad_surfair_tl, rad_air_tl, auxrad_stream_tl) !inout

!---Scattering of the solar beam---------------------------
if((opts%rt_ir%addaerosl .OR. opts%rt_ir%addclouds) .AND. dosolar) &
    call solar_scattering_near_surf_tl(transmission_aux, transmission_aux_tl, &
                                       ircld, raytracing, raytracing_tl, transmission_scatt_ir_tl, &
                                       solar, solar_spectrum, reflectance%refl_out, sateqsun, iv3lay, iv3lev, &
                                       chanprof, narray, pfraction, auxrad_stream, auxrad_stream_tl)

!-------------------------------------------------------------------------------
!2b calculate clear-sky Rayleigh scattering contribution
!-------------------------------------------------------------------------------
if (dosolar) &
  call solar_rayleigh_tl(narray, iv3lev, opts, coef, solar, solar_spectrum, chanprof, &
                         raytracing, raytracing_tl, ircld, profiles, profiles_tl, &
                         profiles_dry, profiles_dry_tl, &
                         transmission_aux, transmission_aux_tl, auxrad_stream_tl)

do i = 1, nchannels
  prof = chanprof(i)%prof
  nstreams = ircld%nstream(prof)
  ist = 0_jpim

  !-------------------------------------------------------------------------------
  ! Add the thermal and solar atmospheric contributions to the clear and cloudy streams
  !-------------------------------------------------------------------------------
  if (thermal(i)) then
    if (do_lambertian(i)) then
      rad_tl%clear(i) = rad_tl%clear(i) + auxrad_stream_tl%meanrad_up(ist,i) + &
        auxrad_stream_tl%meanrad_down(ist,i) * (thermrefl(i) * &
        tau_surf_p * tau_surf) + &
        auxrad_stream%meanrad_down(ist,i) * &
        (thermrefl_tl(i) * tau_surf_p * tau_surf + &
        thermrefl(i) * (tau_surf_p_tl * tau_surf + &
                        tau_surf_p * tau_surf_tl))
    else
      rad_tl%clear(i) = rad_tl%clear(i) + auxrad_stream_tl%meanrad_up(ist,i) + &
                        tau_surf * &
                        (thermrefl(i) * &
                          (auxrad_stream_tl%meanrad_down(ist,i) * tau_surf + &
                          2._JPRB * tau_surf_tl * auxrad_stream%meanrad_down(ist,i)) + &
                        thermrefl_tl(i) * auxrad_stream%meanrad_down(ist,i) * tau_surf)
    endif
  endif

  if (solar(i)) then
    refl_norm_scat = COS(profiles(prof)%zenangle * deg2rad)
    rad_tl%clear(i) = rad_tl%clear(i) + &
                      auxrad_stream_tl%meanrad_up_solar(ist,i) + &
                      transmission_aux%solar_path1%Tau_surf(ist,i) * refl_norm_scat * &
                      (reflectance(i)%refl_out * &
                        (auxrad_stream_tl%meanrad_down_solar(ist,i) * &
                        transmission_aux%solar_path1%Tau_surf(ist,i) + &
                        2._JPRB * auxrad_stream%meanrad_down_solar(ist,i) * &
                        transmission_aux_tl%solar_path1%Tau_surf(ist,i)) + &
                      reflectance_tl(i)%refl_out * &
                      auxrad_stream%meanrad_down_solar(ist,i) * &
                      transmission_aux%solar_path1%Tau_surf(ist,i))
  endif

  if (thermal(i)) auxrad_stream_tl%up(iv2lay(i),ist,i) = auxrad_stream_tl%meanrad_up(ist,i)
  if (solar(i)) auxrad_stream_tl%up_solar(iv2lay(i),ist,i) = auxrad_stream_tl%meanrad_up_solar(ist,i)

  if (thermal(i)) then
    do ist = 1, nstreams
      if (do_lambertian(i)) then
        auxrad_stream_tl%cloudy(ist,i) = &
            auxrad_stream_tl%meanrad_up(ist,i) + &
            auxrad_stream_tl%meanrad_down(ist,i) * (thermrefl(i) * &
            tau_surf_p * tau_surf) + &
            auxrad_stream%meanrad_down(ist,i) * &
            (thermrefl_tl(i) * tau_surf_p * tau_surf + &
            thermrefl(i) * (tau_surf_p_tl * tau_surf + &
                            tau_surf_p * tau_surf_tl))
      else
        auxrad_stream_tl%cloudy(ist,i) = &
            auxrad_stream_tl%meanrad_up(ist,i) + &
            tau_surf * (thermrefl(i) * &
            (auxrad_stream_tl%meanrad_down(ist,i) * tau_surf + &
            auxrad_stream%meanrad_down(ist,i) * &
            tau_surf_tl * 2._JPRB) + &
            thermrefl_tl(i) * auxrad_stream%meanrad_down(ist,i) * &
            tau_surf)
      endif
    enddo
  endif

  if (solar(i)) then
    do ist = 1, nstreams
      auxrad_stream_tl%cloudy(ist,i) = auxrad_stream_tl%cloudy(ist,i) + &
                        auxrad_stream_tl%meanrad_up_solar(ist,i) + &
                        transmission_aux%solar_path1%Tau_surf(ist,i) * refl_norm_scat * &
                        (reflectance(i)%refl_out * &
                          (auxrad_stream_tl%meanrad_down_solar(ist,i) * &
                          transmission_aux%solar_path1%Tau_surf(ist,i) + &
                          2._JPRB * auxrad_stream%meanrad_down_solar(ist,i) * &
                          transmission_aux_tl%solar_path1%Tau_surf(ist,i)) + &
                        reflectance_tl(i)%refl_out * &
                        auxrad_stream%meanrad_down_solar(ist,i) * &
                        transmission_aux%solar_path1%Tau_surf(ist,i))
    enddo
  endif

enddo


!-----------------------
!3. calculate surface emission contribution
!-----------------------
if (dothermal) then
  do i = 1, nchannels
    if (thermal(i)) then
    ! clear sky radiance without reflection term
      ist = 0_jpim

      rad_tl%clear(i) = rad_tl%clear(i) + &
                        auxrad%skin(i) * (emissivity_tl(i)%emis_out * tau_surf + &
                                          emissivity(i)%emis_out * tau_surf_tl) + &
                        rad_skin_tl(i) *  emissivity(i)%emis_out * tau_surf

      prof = chanprof(i)%prof
      do ist=1, ircld%nstream(prof)
        auxrad_stream_tl%cloudy(ist,i) = auxrad_stream_tl%cloudy(ist,i) + auxrad%skin(i) * &
                                        (emissivity_tl(i)%emis_out * tau_surf     + &
                                          emissivity(i)%emis_out * tau_surf_tl) + &
                                          rad_skin_tl(i) * emissivity(i)%emis_out * &
                                          tau_surf
      enddo
    endif
  enddo
endif

!--------------------------------
!4. Add solar surface contribution
!--------------------------------

if (dosolar) &
  call solar_surface_contribution_tl(transmission_aux, transmission_aux_tl, ircld, &
                                     reflectance%refl_out, reflectance_tl%refl_out, refl_norm, &
                                     chanprof, solar, solar_spectrum, narray, auxrad_stream_tl, rad_tl)

!--------------------------------
!5. cosmic temperature correction
!--------------------------------

!calculate planck function corresponding to tcosmic = 2.7k
!deblonde tcosmic for microwave sensors only
if (addcosmic) then
  do i = 1, nchannels
    prof = chanprof(i)%prof
    nstreams = ircld%nstream(prof)

    ist = 0_jpim

    if (do_lambertian(i)) then
      rad_tl%clear(i) = rad_tl%clear(i) + &
        auxrad%cosmic(i) * ( &
        thermrefl_tl(i) * tau_surf_p * tau_surf + &
        thermrefl(i) * ( &
        tau_surf_p_tl * tau_surf + &
        tau_surf_p * tau_surf_tl))
    else
      rad_tl%clear(i) = rad_tl%clear(i) + &
                        auxrad%cosmic(i) * tau_surf * &
                       (thermrefl_tl(i) * tau_surf + &
                        thermrefl(i) * 2.0_jprb * tau_surf_tl)
    endif
  enddo
endif

!---------------------------------------------------
!6. Calculate overcast radiances
!---------------------------------------------------

! Overcast radiances only calculated for non-scattering case

if (.not. keyradonly) then
  ist = 0_jpim
!cdir nodep
  do i = 1, nchannels
    if (thermal(i)) then
      do lay = 1, nlayers
          lev = lay + 1
          ! overcast radiances at given cloud top
          rad_tl%overcast(lay,i) = auxrad_stream_tl%up(lay,ist,i) + &
                                  rad_air_tl(lev,i) * tau_level + &
                                  auxrad%air(lev,i) * tau_level_tl
      enddo
    else if (solar(i)) then
      prof = chanprof(i)%prof
      chan = chanprof(i)%chan
      if (reflectance(i)%refl_cloud_top > 0) then
        refl = reflectance(i)%refl_cloud_top
      else
        if (coef%ff_cwn(chan) > overcast_albedo_wvn) then
          refl = overcast_albedo1 * pi_r
        else
          refl = overcast_albedo2 * pi_r
        endif
      endif
      do lay = 1, nlayers
        lev = lay + 1
        ! overcast radiances at given cloud top
        rad_tl%overcast(lay,i) = solar_spectrum(i) * refl / raytracing%pathsun(lay, prof) * &
                (transmission_aux_tl%solar_path2%Tau_level(lev,ist,i) - &
                transmission_aux%solar_path2%Tau_level(lev,ist,i) * &
                raytracing_tl%pathsun(lay, prof) / raytracing%pathsun(lay, prof)) + &
                auxrad_stream_tl%up_solar(lay,ist,i) + &
                refl / raytracing%pathsat(lay, prof) * &
                transmission_aux%solar_path1%Tau_level(lev,ist,i) * &
                (auxrad_stream_tl%down_solar(lay,ist,i) * &
                transmission_aux%solar_path1%Tau_level(lev,ist,i) - &
                auxrad_stream%down_solar(lay,ist,i) * &
                raytracing_tl%pathsat(lay, prof) / raytracing%pathsat(lay, prof) * &
                transmission_aux%solar_path1%Tau_level(lev,ist,i) + &
                auxrad_stream%down_solar(lay,ist,i) * &
                2_jpim * transmission_aux_tl%solar_path1%Tau_level(lev,ist,i))
      
      enddo
    else
      rad_tl%overcast(:,i) = 0._jprb
    endif
  enddo

  if (dothermal) then
    ! Add surface component to overcast radiances
    do i = 1, nchannels
      if (thermal(i)) then
        lay = iv2lay(i)

        rad_tl%overcast(lay,i) = auxrad_stream_tl%up(lay,ist,i) + &
                                  tau_surf_tl * auxrad%surfair(i) + &
                                  tau_surf    * rad_surfair_tl(i)
      endif
    enddo
  endif
endif

!---------------------------
! 7. Calculate total radiance
!---------------------------
! The simple cloudy scheme is not applied to aerosol-affected radiances
if (opts%rt_ir%addaerosl .OR. opts%rt_ir%addclouds) then
!----------------------------------------
! Calculate complex cloudy radiances
!----------------------------------------
   rad_tl%cloudy = 0._jprb ! nchannels
   do i = 1, nchannels
      prof = chanprof(i)%prof

      do ist = 1, ircld%nstream(prof)
         rad_tl%cloudy(i) = rad_tl%cloudy(i) + & 
                  auxrad_stream_tl%cloudy(ist,i) * (ircld%xstr(ist+1,prof) - ircld%xstr(ist,prof)) + &
                  auxrad_stream%cloudy(ist,i) * (ircld_tl%xstr(ist+1,prof) - ircld_tl%xstr(ist,prof))
      enddo

      rad_tl%cloudy(i)= rad_tl%cloudy(i) + &
                        rad_tl%clear(i) * ircld%XSTRCLR(prof) + &
                        rad%clear(i)    * ircld_tl%XSTRCLR(prof)
   enddo

   rad_tl%total(1:nchannels) = rad_tl%cloudy(1:nchannels)
else
!---------------------------
! Calculate total radiance (clear case/simple cloud)
!---------------------------

  if(opts%rt_ir%pc%addpc) then
    rad_tl%total(1:nchannels) = rad_tl%clear(1:nchannels)
  else

    ! Interpolate to given cloud-top pressures
    do i = 1, nchannels
      prof = chanprof(i)%prof
      lay = aux_prof%s(prof)%nearestlev_ctp - 1

      rad_tl%cloudy(i) = rad_tl%overcast(lay,i) + &
                   (rad_tl%overcast(lay-1,i) - rad_tl%overcast(lay,i)) * aux_prof%s(prof)%pfraction_ctp + &
                   (rad%overcast(lay-1,i)    - rad%overcast(lay,i)   ) * aux_prof_tl%s(prof)%pfraction_ctp
    enddo

    rad_tl%total(1:nchannels) = rad_tl%clear(1:nchannels) + &
                                cfraction(1:nchannels)    * &
                                (rad_tl%cloudy(1:nchannels) - rad_tl%clear(1:nchannels) ) + &
                                cfraction_tl(1:nchannels) * (rad%cloudy(1:nchannels) - rad%clear(1:nchannels))
  endif
endif

if (LHOOK) CALL DR_HOOK('RTTOV_INTEGRATE_TL',1_jpim,ZHOOK_HANDLE)

contains

subroutine calc_atmospheric_radiance_tl(thermal, transmission_aux, transmission_aux_tl, auxrad, ircld, &
                                        chanprof, narray, rad_air_tl, auxrad_stream_tl)

  use rttov_types, Only : transmission_type_aux, radiance_aux, ircld_type, rttov_chanprof
  use parkind1, only : jpim, jprb, jplm

  Implicit None

  logical(jplm), intent(in)                     :: thermal(:)
  type(transmission_type_aux), intent(in)       :: transmission_aux, transmission_aux_tl
  type(radiance_aux), intent(in)                :: auxrad
  type(ircld_type),     intent(in)              :: ircld
  type(rttov_chanprof), intent(in)              :: chanprof(:)
  integer(jpim), intent(in)                     :: narray(:)
  Real(jprb), intent(in)                        :: rad_air_tl(:,:)
  type(radiance_aux), intent(inout)             :: auxrad_stream_tl

  integer(jpim) :: nlevels, nlayers, nstreams, nchannels, prof
  integer(jpim) :: i, ist, lay

  Real(jprb) :: up_laym1, down_laym1
  Real(jprb) :: tau_lev_r_tl, tau_lev_rm1_tl, od_singlelayer_r_tl

  Real(jprb) :: b1_2, b2_2, b3_2
  Real(jprb) :: b1_tl, b2_tl, b3_tl

!unpack narray
  nlevels = narray(1)
  nlayers = nlevels - 1_jpim
  nchannels = narray(3)

  if(transmission_aux%anynegtau .gt. 0.0_jprb) then

     do i = 1, nchannels
      if (thermal(i)) then
        prof = chanprof(i)%prof
        nstreams = ircld%nstream(prof)
        do ist = 0, nstreams
          do lay = 1, nlayers
              isti = ircld%icldarr(ist,lay,prof)

              B1_2 = auxrad%air(lay,i) * dtau_lay
              B2_2 = daux_lay * tau_level
              B3_2 = daux_lay * dtau_lay * od_singlelayer_r

              od_singlelayer_r_tl = -od_singlelayer_tl * od_singlelayer_r**2
              tau_lev_r_tl   = -tau_level_tl * tau_level_r**2
              tau_lev_rm1_tl = -tau_layer_tl * tau_layer_r**2

              B1_TL = rad_air_tl(lay,i) * dtau_lay + &
                      auxrad%air(lay,i) * dtau_lay_tl
              B2_TL = daux_lay_tl * tau_level + &
                      daux_lay * tau_level_tl
              B3_TL = od_singlelayer_r * &
                      (daux_lay_tl * dtau_lay + &
                      daux_lay * dtau_lay_tl) + &
                      od_singlelayer_r_tl * &
                      (daux_lay * dtau_lay)

              auxrad_stream_tl%up(lay,ist,i)   = transmission_aux%fac1(lay,ist,i) * (B1_TL - B2_TL + B3_TL)
              auxrad_stream_tl%down(lay,ist,i) = transmission_aux%fac1(lay,ist,i) * &
                                    transmission_aux%thermal_path1%fac2(lay+1,ist,i) * &
                                  (tau_level_r * ( &
                                  (B1_TL - B3_TL) * (tau_layer_r) + &
                                  (B2_TL * tau_level_r + &
                                    B2_2 * (2._jprb * tau_lev_r_tl))) + &
                                  (B1_2 - B3_2) * (tau_lev_rm1_tl * tau_level_r + &
                                    tau_lev_r_tl * tau_layer_r))
          enddo
        enddo
      endif
    enddo

! DAR: This is slow. Hopefully it'll never run.
    do  i = 1,nchannels
      if (thermal(i)) then
        prof = chanprof(i)%prof
        nstreams = ircld%nstream(prof)
        do ist = 0, nstreams
            do lay = 1,nlayers
              if(tau_level < 0.0_jprb .and. tau_layer >= 0._jprb) then 
                  auxrad_stream_tl%up(lay, ist, i) = 0.5_JPRB  * &
                                      ((auxrad%air(lay+1,i) + auxrad%air(lay,i)) * &
                                      (tau_layer_tl - tau_level_tl) + &
                                      (rad_air_tl(lay+1,i) + rad_air_tl(lay,i)) * &
                                      (tau_layer - tau_level))
                  auxrad_stream_tl%down(lay, ist, i) = 0.0_jprb
              endif
            enddo
        enddo
      endif
    enddo

!Cumulative sum
  do i = 1,nchannels
    if (thermal(i)) then
      prof = chanprof(i)%prof
      nstreams = ircld%nstream(prof)
      do ist = 0, nstreams
         do lay = 2, nlayers
            auxrad_stream_tl%up(lay,ist,i) = auxrad_stream_tl%up(lay,ist,i) + &
                                             auxrad_stream_tl%up(lay-1,ist,i)
            auxrad_stream_tl%down(lay,ist,i) = auxrad_stream_tl%down(lay,ist,i) + &
                                               auxrad_stream_tl%down(lay-1,ist,i)
         enddo
      enddo
    endif
  enddo
else ! DAR: fast code (on IBM) when no special cases
  do i = 1, nchannels
    if (thermal(i)) then
      prof = chanprof(i)%prof
      do ist = 0, ircld%nstream(prof)
         up_laym1 = 0._jprb
         down_laym1 = 0._jprb
         do lay = 1, nlayers
            isti = ircld%icldarr(ist,lay,prof)

            B1_2 = auxrad%air(lay,i) * dtau_lay
            B2_2 = daux_lay * tau_level
            B3_2 = daux_lay * dtau_lay * od_singlelayer_r

            od_singlelayer_r_tl = -od_singlelayer_tl * od_singlelayer_r**2
            tau_lev_r_tl   = -tau_level_tl * tau_level_r**2
            tau_lev_rm1_tl = -tau_layer_tl * tau_layer_r**2

            B1_TL = rad_air_tl(lay,i) * dtau_lay + &
                    auxrad%air(lay,i) * dtau_lay_tl
            B2_TL = daux_lay_tl * tau_level  + &
                    daux_lay * tau_level_tl
            B3_TL = od_singlelayer_r * &
                    (daux_lay_tl * dtau_lay + &
                    daux_lay * dtau_lay_tl) + &
                    od_singlelayer_r_tl * (daux_lay * dtau_lay)

            auxrad_stream_tl%up(lay,ist,i) = up_laym1 + transmission_aux%fac1(lay,ist,i) * &
                    (B1_TL - B2_TL + B3_TL)

            if (do_lambertian(i)) then
              auxrad_stream_tl%down(lay,ist,i) = down_laym1 + transmission_aux%fac1(lay,ist,i) * &
                transmission_aux%thermal_path1%fac2(lay+1,ist,i) * 0.5_jprb * &
                ((rad_air_tl(lay,i) + rad_air_tl(lay+1,i)) * &
                ((tau_layer_p) - (tau_level_p)) * &
                ((tau_level_p_r) * (tau_layer_p_r)) + &
                (auxrad%air(lay,i) + auxrad%air(lay+1,i)) * (&
                ((tau_layer_p_tl) - (tau_level_p_tl)) * &
                ((tau_level_p_r) * (tau_layer_p_r)) + &
                ((tau_layer_p) - (tau_level_p)) * ( &
                ((tau_level_p_r_tl) * (tau_layer_p_r)) + &
                (tau_level_p_r) * (tau_layer_p_r_tl))))
            else
              auxrad_stream_tl%down(lay,ist,i) = down_laym1 + transmission_aux%fac1(lay,ist,i) * &
                      transmission_aux%thermal_path1%fac2(lay+1,ist,i) * &
                      (tau_level_r * ( &
                      (B1_TL - B3_TL) * (tau_layer_r) + &
                      (B2_TL * tau_level_r + &
                      B2_2 * (2._jprb * tau_lev_r_tl))) + &
                      (B1_2 - B3_2)  * (tau_lev_rm1_tl * tau_level_r + &
                      tau_lev_r_tl * tau_layer_r))
            endif

            up_laym1   = auxrad_stream_tl%up(lay,ist,i)
            down_laym1 = auxrad_stream_tl%down(lay,ist,i)
         enddo
      enddo
    endif
  enddo
endif

do i = 1, nchannels
  if (thermal(i)) then
    prof = chanprof(i)%prof
    do ist = 0, ircld%nstream(prof)
      do lay = 1, nlayers
        if(auxrad_stream%down_ref(lay,ist,i) < 0._jprb) then
          auxrad_stream_tl%down(lay,ist,i) = 0.0_jprb
        endif
      enddo
    enddo
  endif
enddo

end subroutine calc_atmospheric_radiance_tl

! DAR: This is the next biggest consumer of CPU time and should be looked at next. The trouble seems to be that you have to change a
!      lot of non-consecutive data. Maybe the RTTOV 9 comments can be removed
subroutine calc_near_surf_contribution_tl(thermal, transmission_aux, transmission_aux_tl, auxrad, ircld, & !in
                                       chanprof, narray, iv3lay, iv3lev, pol_id, & !in
                                       rad_surfair_tl, rad_air_tl, auxrad_stream_tl) ! inout

  use rttov_const, only : min_od
  use rttov_types, Only : rttov_chanprof, ircld_type, transmission_type_aux, radiance_aux
  use parkind1, only : jpim, jprb, jplm

  Implicit None

  logical(jplm), intent(in)               :: thermal(:)
  type(transmission_type_aux), intent(in) :: transmission_aux, transmission_aux_tl
  type(radiance_aux), intent(in)          :: auxrad
  type(ircld_type), intent(in)            :: ircld
  type(rttov_chanprof), intent(in)        :: chanprof(:)
  integer(jpim), intent(in)               :: narray(:)
  integer(jpim), intent(in)               :: iv3lay(:)
  integer(jpim), intent(in)               :: iv3lev(:)
  integer(jpim), intent(in)               :: pol_id(:)
  real(jprb), intent(in)                  :: rad_surfair_tl(:), rad_air_tl(:,:)
  type(radiance_aux), intent(inout)       :: auxrad_stream_tl

  integer(jpim) :: nlayers, nlevels, nstreams, nchannels
  integer(jpim) :: prof

  integer(jpim) :: i, lev, lay, ist

  real(jprb) ::  rad_tmp, rad_tmp_tl

#define B1_3 (auxrad%air(lev,i) * (tau_level  - tau_surf))
#define B2_3 ((auxrad%surfair(i) - auxrad%air(lev,i)) * tau_surf)
#define B3_3 ((auxrad%surfair(i) - auxrad%air(lev,i)) * (tau_level - tau_surf) * (transmission_aux%thermal_path1%od_sfrac_r(ist,i)))
#define B1_TL1 (rad_air_tl(lev,i) * (tau_level  - tau_surf))
#define B1_TL2 (auxrad%air(lev,i) * (tau_level_tl - tau_surf_tl))
#define B2_TL1 (tau_surf * (rad_surfair_tl(i) - rad_air_tl(lev,i)))
#define B2_TL2 (tau_surf_tl * (auxrad%surfair(i) - auxrad%air(lev,i)))
#define B3_TL1 ((auxrad%surfair(i) - auxrad%air(lev,i)) * (tau_level_tl - tau_surf_tl))
#define B3_TL2 ((rad_surfair_tl(i) - rad_air_tl(lev,i)) * (tau_level - tau_surf))
#define B3_TL3 (-B3_3 * transmission_aux_tl%thermal_path1%od_sfrac(ist,i))
#define B3_TL4 (transmission_aux%thermal_path1%od_sfrac_r(ist,i))

#define B1_L_TL ((tau_level_p_tl)  - (tau_surf_p_tl))
#define B2_L_TL (rad_surfair_tl(i) + rad_air_tl(lev,i))
#define B3_L_TL1 ((tau_level_p_r_tl) * (tau_surf_p_r))
#define B3_L_TL2 ((tau_level_p_r) * (tau_surf_p_r_tl))
#define B1_L ((tau_level_p)  - (tau_surf_p))
#define B2_L (auxrad%surfair(i) + auxrad%air(lev,i))
#define B3_L ((tau_level_p_r) * (tau_surf_p_r))


! z=x/y -=> z' = x'-zy'

! unpack
  nlevels = narray(1)
  nlayers = nlevels - 1_jpim
  nstreams = narray(2)
  nchannels = narray(3)

  Do i = 1, nchannels
    if (thermal(i)) then
      prof = chanprof(i)%prof
      nstreams = ircld%nstream(prof)

      lay = iv3lay(i)
      lev = iv3lev(i)

      do ist = 0, nstreams !ircld%nstreams(prof)

          if(tau_surf < 0._JPRB) then ! DAR 0 < min_tau
            if(tau_level >= 0._JPRB) then

              auxrad_stream_tl%meanrad_up(ist,i) = 0.5_JPRB * (((rad_surfair_tl(i) + rad_air_tl(lev,i)) * &
                                    (tau_level - tau_surf)) + &
                                  ((auxrad%surfair(i)     + auxrad%air(lev,i)) * &
                                    (tau_level_tl - tau_surf_tl)))

              if (pol_id(i) >= 6_jpim ) auxrad_stream_tl%meanrad_up(ist,i) = 0.0_jprb
            endif

            auxrad_stream_tl%meanrad_down(ist,i) = auxrad_stream_tl%down(lay,ist,i)
            auxrad_stream_tl%meanrad_up(ist,i)   = auxrad_stream_tl%up(lay,ist,i) + &
                                                   auxrad_stream_tl%meanrad_up(ist,i)

          else !Tau_level<0 .or. transmission_aux%thermal_path1%tau_surf(ist,i) >= 0

            if(transmission_aux%thermal_path1%od_sfrac(ist,i) < min_od .OR. &
                ((tau_level - tau_surf) < min_od)) THEN

              ! small optical depth or optical depth change set radiance to zero

              auxrad_stream_tl%meanrad_down(ist,i) = auxrad_stream_tl%down(lay,ist,i)
              auxrad_stream_tl%meanrad_up(ist,i)   = auxrad_stream_tl%up(lay,ist,i)

            else
              ! compute linear in tau radiances
              ! PGF90 compatibility for b3_tl macros
              auxrad_stream_tl%meanrad_up(ist,i) = &
              B1_TL1 + &
              B1_TL2 - &
              (B2_TL1 + &
              B2_TL2) + &
(B3_TL1 + & 
              B3_TL2 + &
B3_TL3) * &
              B3_TL4

              if (do_lambertian(i)) then
                rad_tmp_tl = &
                  transmission_aux%surf_fac(ist,i) * 0.5_jprb * (&
                  B1_L_TL * &
                  B2_L * &
                  B3_L + &
                  B1_L * (&
                  B2_L_TL * &
                  B3_L + &
                  B2_L * &
                  (B3_L_TL1 + &
                   B3_L_TL2)))

              else
                rad_tmp = transmission_aux%surf_fac(ist,i) * &
                (B1_3 - &
                B3_3) * &
                (tau_level_r * tau_surf_r) + &
                B2_3 * (tau_surf_r)**2_jpim

                rad_tmp_tl = tau_surf_r * ( &
                    (tau_level * ((B1_TL1 + &
                    B1_TL2) - &
(B3_TL1 + &
                    B3_TL2 + &
B3_TL3) * &
                    B3_TL4) - &
                    (B1_3 - &
                    B3_3) * &
                    tau_level_tl) * &
                    tau_level_r**2_jpim + &
                    ((B2_TL1 + &
                    B2_TL2) * &
                    tau_surf - B2_3 * tau_surf_tl) * &
                    tau_surf_r**2_jpim - &
                    rad_tmp * tau_surf_tl)
                  ! f = 1/z * (x/y + w/z) => f' = 1/z * ((yx'-xy')/y^2 + (wz'-zw')/z^2 - f)

              endif

              if (pol_id(i) >= 6_jpim ) auxrad_stream_tl%meanrad_up(ist,i) = 0.0_jprb

              auxrad_stream_tl%meanrad_down(ist,i) = auxrad_stream_tl%down(lay,ist,i) + rad_tmp_tl
              auxrad_stream_tl%meanrad_up(ist,i)   = auxrad_stream_tl%up(lay,ist,i)   + &
                                                       auxrad_stream_tl%meanrad_up(ist,i)

            endif
          endif
          if (pol_id(i) >= 6_jpim) auxrad_stream_tl%meanrad_up(ist,i) = 0.0_jprb
      enddo
    endif
  enddo

end subroutine calc_near_surf_contribution_tl

subroutine solar_scattering_air_tl(transmission_aux, transmission_aux_tl, ircld, raytracing, raytracing_tl, &
                                   chanprof, narray, transmission_scatt_ir_tl, &
                                   solar, solar_spectrum, refl, sateqsun, auxrad_stream, auxrad_stream_tl)

  use rttov_const, only : z4pi_r
  use rttov_types, Only : transmission_type_aux, ircld_type, radiance_aux, raytracing_type, &
                          rttov_chanprof, transmission_scatt_ir_type
  use parkind1, only : jpim, jplm, jprb

  Implicit None

  type(transmission_type_aux), intent(in)      :: transmission_aux, transmission_aux_tl
  type(ircld_type),     intent(in)             :: ircld
  type(raytracing_type), intent(in)            :: raytracing, raytracing_tl
  type(radiance_aux), intent(in)               :: auxrad_stream
  type(rttov_chanprof), intent(in)             :: chanprof(:)
  type(transmission_scatt_ir_type), intent(in) :: transmission_scatt_ir_tl
  logical(jplm), intent(in)                    :: solar(:)
  real(jprb), intent(in)                       :: solar_spectrum(:)
  real(jprb), intent(in)                       :: refl(:)
  logical(jplm), intent(in)                    :: sateqsun(:,:)
  integer(jpim), intent(in)                    :: narray(:)
  type(radiance_aux), intent(inout)            :: auxrad_stream_tl

!local variables
  integer(jpim) :: chan, prof
  integer(jpim) :: nlevels, nlayers, nstreams, maxnstreams, nchannels
  integer(jpim) :: i, ist, lay

  real(jprb)    :: temp_tl(1:narray(1)-1_jpim,0:narray(2))

  Real(jprb) :: fac1_2_tl(1:narray(1)-1_jpim,0:narray(2)), fac2_2_tl(1:narray(1)-1_jpim,0:narray(2)), &
                fac3_2_tl(1:narray(1)-1_jpim), fac4_2_tl(1:narray(1)-1_jpim,0:narray(2)), &
                fac5_2_tl(1:narray(1)-1_jpim,0:narray(2)), fac6_2_tl(1:narray(1)-1_jpim,0:narray(2)), &
                fac7_2_tl(1:narray(1)-1_jpim)

! macro definitions for different factors

#define fac1_2 auxrad_stream%Fac1_2(lay,ist,i)
#define fac2_2 auxrad_stream%Fac2_2(lay,ist,i)
#define fac3_2 auxrad_stream%Fac3_2(lay,ist,i)
#define fac4_2 auxrad_stream%Fac4_2(lay,ist,i)
#define fac5_2 auxrad_stream%Fac5_2(lay,ist,i)
#define fac6_2 auxrad_stream%Fac6_2(lay,ist,i)
#define fac7_2 auxrad_stream%Fac7_2(lay,ist,i)
#define dfac54_2 (fac5_2 - fac4_2)

#define fac1_2_tl Fac1_2_tl(lay,ist)
#define fac2_2_tl Fac2_2_tl(lay,ist)
#define fac3_2_tl Fac3_2_tl(lay)
#define fac4_2_tl Fac4_2_tl(lay,ist)
#define fac5_2_tl Fac5_2_tl(lay,ist)
#define fac6_2_tl Fac6_2_tl(lay,ist)
#define fac7_2_tl Fac7_2_tl(lay)
#define dfac54_2_tl (fac5_2_tl - fac4_2_tl)

#define tausun2_layer_tl transmission_aux_tl%solar_path2%Tau_level(lay,ist,i)
#define tausun2_level_tl transmission_aux_tl%solar_path2%Tau_level(lay+1,ist,i)
#define tausun2_layer transmission_aux%solar_path2%Tau_level(lay,ist,i)
#define tausun2_level transmission_aux%solar_path2%Tau_level(lay+1,ist,i)

#define tausun1_layer_tl transmission_aux_tl%solar_path1%Tau_level(lay,ist,i)
#define tausun1_level_tl transmission_aux_tl%solar_path1%Tau_level(lay+1,ist,i)
#define tausun1_layer transmission_aux%solar_path1%Tau_level(lay,ist,i)
#define tausun1_level transmission_aux%solar_path1%Tau_level(lay+1,ist,i)
#define tausun1_level_r transmission_aux%solar_path1%Tau_level_r(lay+1,ist,i)
#define tausun1_surf transmission_aux%solar_path1%Tau_surf(ist,i)
#define tausun1_surf_r transmission_aux%solar_path1%Tau_surf_r(ist,i)
#define tausun1_surf_tl transmission_aux_tl%solar_path1%Tau_surf(ist,i)

! unpack
  nlevels = narray(1)
  nlayers = nlevels - 1_jpim
  maxnstreams = narray(2)
  nchannels = narray(3)

  do i = 1, nchannels
    if(solar(i)) then

      chan = chanprof(i)%chan
      prof = chanprof(i)%prof
      nstreams = ircld%nstream(prof)

      do ist = 0, nstreams
        do lay = 1, nlayers
          isti = ircld%icldarr(ist,lay,prof)

          fac6_2_tl = solar_spectrum(i) * z4pi_r * &
                      transmission_scatt_ir_tl%azphacdo(isti,lay,i)

          fac1_2_tl = solar_spectrum(i) * z4pi_r * &
                      transmission_scatt_ir_tl%azphacup(isti,lay,i)

          fac2_2_tl = transmission_scatt_ir_tl%ssa(isti,lay,i)

          fac4_2_tl = -transmission_aux_tl%solar_path2%Od_singlelayer(isti,lay,i) * fac4_2
          fac5_2_tl = -transmission_aux_tl%solar_path1%Od_singlelayer(isti,lay,i) * fac5_2
        enddo
      enddo

      do lay = 1, nlayers
        fac3_2_tl = raytracing_tl%pathsat(lay,prof) / raytracing%patheff(lay,prof) - &
                    raytracing%pathsat(lay,prof) * raytracing_tl%patheff(lay,prof) / &
                    raytracing%patheff(lay,prof)**2_jpim
      enddo

  !----------------Upward single scattering of the solar beam-----------------------
      do ist = 0, nstreams
        do lay = 1, nlayers
          auxrad_stream_tl%up_solar(lay,ist,i) = (tausun2_layer - tausun2_level) * &
                (fac1_2_tl * fac2_2 * fac3_2 + & 
                 fac2_2_tl * fac3_2 * fac1_2 + & 
                 fac3_2_tl * fac1_2 * fac2_2) + &
                (fac1_2 * fac2_2 * fac3_2 * &
                (tausun2_layer_tl - tausun2_level_tl))
        enddo
      enddo

      do ist = 0, nstreams
        do lay = 2, nlayers
          auxrad_stream_tl%up_solar(lay,ist,i) = auxrad_stream_tl%up_solar(lay,ist,i) + &
                                                 auxrad_stream_tl%up_solar(lay-1,ist,i)
        enddo
      enddo

  !-------------------Downward single scattering of the solar beam------------------

      if (refl(i) > 0._jprb) then

        do lay = 1, nlayers
          if(.not. sateqsun(lay,prof)) then
            fac7_2_tl = (raytracing%pathsun(lay,prof) * raytracing_tl%pathsat(lay,prof) - &
                         raytracing_tl%pathsun(lay,prof) * raytracing%pathsat(lay,prof)) / &
                        (raytracing%pathsun(lay,prof) - raytracing%pathsat(lay,prof))**2_jpim
          endif
        enddo

        do ist = 0, nstreams
          do lay = 1, nlayers
! nested multiplication of derivative
            if(.not. sateqsun(lay,prof)) then
              temp_tl(lay,ist) = &
                transmission_aux%solar_path1%fac2(lay+1,ist,i) * &
                (1.0_jprb / (tausun1_layer * tausun1_level)) * (&
                  fac6_2_tl * &
                  fac2_2 * &
                  fac7_2 * &
                  dfac54_2 * &
                  tausun2_layer + &
                  fac6_2 * (&
                    fac2_2_tl * &
                    fac7_2 * &
                    dfac54_2 * &
                    tausun2_layer + &
                    fac2_2 * (&
                      fac7_2_tl * &
                      dfac54_2 * &
                      tausun2_layer + &
                      fac7_2 * (&
                        dfac54_2_tl * &
                        tausun2_layer + &
                        dfac54_2 * (&
                          (tausun2_layer_tl - &
                          (tausun2_layer/(tausun1_layer * tausun1_level))*&
                          (tausun1_layer_tl * tausun1_level + &
                          tausun1_layer * tausun1_level_tl)))))))
            else
              isti = ircld%icldarr(ist,lay,prof)
              temp_tl(lay,ist) = &
                transmission_aux%solar_path1%fac2(lay+1,ist,i) * &
                (1.0_jprb / (tausun1_layer * tausun1_level)) * (&
                  fac6_2_tl * &
                  fac2_2 * &
                  fac4_2 * &
                  transmission_aux%solar_path2%Od_singlelayer(isti,lay,i) * &
                  tausun2_layer + &
                  fac6_2 * (&
                    fac2_2_tl * &
                    fac4_2 * &
                    transmission_aux%solar_path2%Od_singlelayer(isti,lay,i) * &
                    tausun2_layer + &
                    fac2_2 * (&
                      fac4_2_tl * &
                      transmission_aux%solar_path2%Od_singlelayer(isti,lay,i) * &
                      tausun2_layer + &
                      fac4_2 * (&
                        transmission_aux_tl%solar_path2%Od_singlelayer(isti,lay,i) * &
                        tausun2_layer + &
                        dfac54_2 * (&
                          (tausun2_layer_tl - &
                          (tausun2_layer/(tausun1_layer * tausun1_level))*&
                          (tausun1_layer_tl * tausun1_level + &
                          tausun1_layer * tausun1_level_tl)))))))
            endif
          enddo
        enddo

        do ist = 0, nstreams
          lay = 1_jpim
          if(auxrad_stream%down_ref_solar(lay,ist,i) < 0._jprb) then
            auxrad_stream_tl%down_solar(lay,ist,i) = 0._jprb
          else
            auxrad_stream_tl%down_solar(lay,ist,i) = temp_tl(lay,ist)
          endif

          do lay = 2, nlayers
            temp_tl(lay,ist) = temp_tl(lay-1,ist) + temp_tl(lay,ist)
            if(auxrad_stream%down_ref_solar(lay,ist,i) < 0._jprb) then
              auxrad_stream_tl%down_solar(lay,ist,i) = 0.0_jprb
            else
              auxrad_stream_tl%down_solar(lay,ist,i) = temp_tl(lay,ist)
            endif
          enddo
        enddo

      endif ! refl(i) > 0.

    endif
  enddo

end subroutine solar_scattering_air_tl

subroutine solar_scattering_near_surf_tl(transmission_aux, transmission_aux_tl, &
                                         ircld, raytracing, raytracing_tl, transmission_scatt_ir_tl, &
                                         solar, solar_spectrum, refl, sateqsun, iv3lay, iv3lev, chanprof, narray, &
                                         pfraction, auxrad_stream, auxrad_stream_tl)

  use rttov_const, only : z4pi_r
  use rttov_types, Only : transmission_type_aux, ircld_type, radiance_aux, &
                          rttov_chanprof, transmission_scatt_ir_type, raytracing_type
  use parkind1, only : jpim, jprb, jplm 

  Implicit None

  type(transmission_type_aux), intent(in)      :: transmission_aux, transmission_aux_tl
  type(ircld_type),     intent(in)             :: ircld
  type(raytracing_type), intent(in)            :: raytracing, raytracing_tl
  type(transmission_scatt_ir_type), intent(in) :: transmission_scatt_ir_tl
  logical(jplm), intent(in)                    :: solar(:)
  real(jprb), intent(in)                       :: solar_spectrum(:)
  real(jprb), intent(in)                       :: refl(:)
  logical(jplm), intent(in)                    :: sateqsun(:,:)
  integer(jpim), intent(in)                    :: iv3lay(:), iv3lev(:)
  type(rttov_chanprof), intent(in)             :: chanprof(:)
  integer(jpim), intent(in)                    :: narray(:)
  real(jprb), intent(in)                       :: pfraction(:)

  type(radiance_aux), intent(in)               :: auxrad_stream
  type(radiance_aux), intent(inout)            :: auxrad_stream_tl

  integer(jpim) :: prof, chan
  integer(jpim) :: nlevels, nlayers, nstreams, nchannels
  integer(jpim) :: i, ist, lay, lev, lay1

  Real(jprb) :: fac1_3_tl(0:narray(2)), fac2_3_tl(0:narray(2)), fac3_3_tl, fac4_3_tl(0:narray(2)), &
                fac5_3_tl(0:narray(2)), fac6_3_tl(0:narray(2)), fac7_3_tl(0:narray(2))

#define fac1_3 auxrad_stream%Fac1_3(ist,i)
#define fac2_3 auxrad_stream%Fac2_3(ist,i)
#define fac3_3 auxrad_stream%Fac3_3(ist,i)
#define fac4_3 auxrad_stream%Fac4_3(ist,i)
#define fac5_3 auxrad_stream%Fac5_3(ist,i)
#define fac6_3 auxrad_stream%Fac6_3(ist,i)
#define fac7_3 auxrad_stream%Fac7_3(ist,i)
#define dfac54_3 (fac5_3 - fac4_3)

#define dfac54_3_tl (fac5_3_tl(ist) - fac4_3_tl(ist))

! unpack
  nlevels = narray(1)
  nlayers = nlevels - 1_jpim
  nstreams = narray(2)
  nchannels = narray(3)  

  do i=1, nchannels
     if(solar(i)) then
        prof = chanprof(i)%prof
        chan = chanprof(i)%chan
        nstreams = ircld%nstream(prof)
        lay = iv3lay(i)
        lev = iv3lev(i)

        if(pfraction(i) < 0.0_JPRB) then
           lay1 = lay
        else
           lay1 = lay + 1_jpim
        endif

        do ist = 0, nstreams
          isti = ircld%icldarr(ist,lay,prof)
          fac4_3_tl(ist) = -transmission_aux_tl%solar_path2%od_sfrac(ist,i) * fac4_3 !(exp(-x))'=-exp(-x)
          fac5_3_tl(ist) = -transmission_aux_tl%solar_path1%od_sfrac(ist,i) * fac5_3

          fac1_3_tl(ist) = solar_spectrum(i) * z4pi_r * transmission_scatt_ir_tl%azphacup(isti,lay,i)

          fac2_3_tl(ist) = transmission_scatt_ir_tl%ssa(isti,lay,i)
        enddo

        fac3_3_tl = raytracing_tl%pathsat(lay,prof) / raytracing%patheff(lay,prof) - &
                    raytracing%pathsat(lay,prof) * raytracing_tl%patheff(lay,prof) / &
                    raytracing%patheff(lay,prof)**2_jpim

!--------------Upward single scattering of the solar beam-------------------------
        do ist = 0, nstreams       
           auxrad_stream_tl%meanrad_up_solar(ist,i) = fac1_3_tl(ist) * fac2_3 * fac3_3 * &
                                 (tausun2_level - transmission_aux%solar_path2%Tau_surf(ist,i)) + &
                                  fac1_3 * (fac2_3_tl(ist) * fac3_3 * &
                                 (tausun2_level - transmission_aux%solar_path2%Tau_surf(ist,i)) + &
                                  fac2_3 * (fac3_3_tl  * &
                                 (tausun2_level - transmission_aux%solar_path2%Tau_surf(ist,i)) + &
                                  fac3_3 * (&
                                 (tausun2_level_tl - transmission_aux_tl%solar_path2%Tau_surf(ist,i)))))
        enddo

        do ist = 0, nstreams
          auxrad_stream_tl%meanrad_up_solar(ist,i)   = auxrad_stream_tl%up_solar(lay,ist,i)   + &
                                                       auxrad_stream_tl%meanrad_up_solar(ist,i)
        enddo

!--------------Downward single scattering of the solar beam-----------------------

        if (refl(i) > 0._jprb) then

          if(.not. sateqsun(lay1,prof))then
             do ist = 0, nstreams
                isti = ircld%icldarr(ist,lay,prof)
                fac6_3_tl(ist) = solar_spectrum(i) * z4pi_r * transmission_scatt_ir_tl%azphacdo(isti,lay,i)

                fac7_3_tl(ist) = (raytracing%pathsun(lay1,prof)    * raytracing_tl%pathsat(lay1,prof)  - &
                                  raytracing_tl%pathsun(lay1,prof) * raytracing%pathsat(lay1,prof)   ) / &
                                 (raytracing%pathsun(lay1,prof) - raytracing%pathsat(lay1,prof))**2

                auxrad_stream_tl%meanrad_down_solar(ist,i) = transmission_aux%solar_path1%fac2(lev,ist,i) * &
                            (tausun1_level_r * tausun1_surf_r) * (&
                             fac6_3_tl(ist) * fac2_3 * fac7_3 * dfac54_3 * tausun2_level + &
                             fac6_3 * (&
                                      fac2_3_tl(ist) * fac7_3 * dfac54_3 * tausun2_level + &
                                              fac2_3 * (&
                                               fac7_3_tl(ist) * dfac54_3 * tausun2_level + &
                                                       fac7_3 * (&
                                                             dfac54_3_tl * tausun2_level + &
                                                                dfac54_3 * (&
                            (tausun2_level_tl - &
                            (tausun2_level * (tausun1_level_r * tausun1_surf_r))*&
                            (tausun1_level_tl * tausun1_surf + &
                             tausun1_level * tausun1_surf_tl)))))))
             enddo
          else
             do ist = 0, nstreams
                auxrad_stream_tl%meanrad_down_solar(ist,i) = transmission_aux%solar_path1%fac2(lev,ist,i) * &
                            (tausun1_level_r * tausun1_surf_r) * (&
                             fac6_3_tl(ist) * fac2_3 * fac4_3 * transmission_aux%solar_path2%od_sfrac(ist,i) * &
                             tausun2_level + &
                               fac6_3 * (&
                                 fac2_3_tl(ist) * fac4_3 * transmission_aux%solar_path2%od_sfrac(ist,i) * &
                                 tausun2_level + &
                                   fac2_3 * (&
                                     fac4_3_tl(ist) * transmission_aux%solar_path2%od_sfrac(ist,i) * &
                                     tausun2_level + &
                                       fac4_3 * (&
                                         transmission_aux_tl%solar_path2%od_sfrac(ist,i) * &
                                         tausun2_level + &
                                         transmission_aux%solar_path2%od_sfrac(ist,i) * (&
                            (tausun2_level_tl - &
                            (tausun2_level/(tausun1_level * tausun1_surf)) * &
                            (tausun1_level_tl * tausun1_surf + &
                            tausun1_level * tausun1_surf_tl)))))))

             enddo
          endif

          do ist = 0, nstreams
            auxrad_stream_tl%meanrad_down_solar(ist,i) = auxrad_stream_tl%down_solar(lay,ist,i) + &
                                                         auxrad_stream_tl%meanrad_down_solar(ist,i)

            if(auxrad_stream%meanrad_down_solar(ist,i) < 0._jprb)then
              auxrad_stream_tl%meanrad_down_solar(ist,i) = 0._jprb
            endif
          enddo

       endif ! refl(i) > 0.

     endif

  enddo

end subroutine solar_scattering_near_surf_tl

subroutine solar_rayleigh_tl(narray, iv3lev, opts, coef, solar, solar_spectrum, chanprof, &
                             raytracing, raytracing_tl, ircld, profiles, profiles_tl, &
                             profiles_dry, profiles_dry_tl, &
                             transmission_aux, transmission_aux_tl, auxrad_stream_tl)

  use rttov_types, Only : rttov_coef, rttov_chanprof, raytracing_type, profile_type, transmission_type_aux, &
                          radiance_aux, rttov_options, ircld_type

  use parkind1, only : jprb, jplm, jpim

  use rttov_const, Only : z4pi_r, deg2rad, gravity, na, Mh2o, Mair, &
                          ray_scs_wlm, ray_scs_a1, ray_scs_b1, ray_scs_c1, ray_scs_d1, &
                          ray_scs_a2, ray_scs_b2, ray_scs_c2, ray_scs_d2, ray_min_wvn

  implicit none

  integer(jpim), intent(in)                  :: narray(:)
  integer(jpim), intent(in)                  :: iv3lev(:)
  type(rttov_options), intent(in)            :: opts
  type(rttov_coef), intent(in)               :: coef
  logical(jplm), intent(in)                  :: solar(:)
  real(jprb), intent(in)                     :: solar_spectrum(:)
  type(rttov_chanprof), intent(in)           :: chanprof(:)
  type(raytracing_type), intent(in)          :: raytracing
  type(raytracing_type), intent(in)          :: raytracing_tl
  type(ircld_type), intent(in)               :: ircld
  type(profile_type), intent(in)             :: profiles(:)
  type(profile_type), intent(in)             :: profiles_tl(:)
  type(profile_type), intent(in)             :: profiles_dry(:)
  type(profile_type), intent(in)             :: profiles_dry_tl(:)
  type(transmission_type_aux), intent(in)    :: transmission_aux
  type(transmission_type_aux), intent(in)    :: transmission_aux_tl
  type(radiance_aux), intent(inout)          :: auxrad_stream_tl

  real(kind=jprb) :: wlm, ss_param, v_h2o, v_h2o_dry, Mwet
  real(kind=jprb) :: cossat, cossol, cosazi, cosscata_term1, cosscata_term2, cosscata
  real(kind=jprb) :: ray_phase, solar_src, solar_src_updn
  real(kind=jprb) :: v_h2o_tl, v_h2o_dry_tl, Mwet_tl, plev_tl, plev1_tl, dp_tl
  real(kind=jprb) :: cosscata_term1_tl, cosscata_term2_tl, cosscata_tl
  real(kind=jprb) :: ray_phase_tl, solar_src_tl, solar_src_updn_tl

  integer(kind=jpim) :: ist, nchannels, nlevels, chan, prof, i, lev, lay
  real(kind=jprb) :: rayrad_up_tl(0:narray(1)-1,0:narray(2)), rayrad_dn_tl(0:narray(1)-1,0:narray(2))

  nchannels = narray(3)
  nlevels = narray(1)

  ! The scattering cross-section calculation is based on Bucholzt '95.
  ! The phase function could be made to account for depolarisation.
  ! The solar geometry is approximated using the sun-surface path values,
  ! and similarly for the sun-level-satellite transmittances.

  ! Rayleigh scattering currently only included for channels less than 2um (ray_min_wvn).

  do i = 1, nchannels
    prof = chanprof(i)%prof
    chan = chanprof(i)%chan

    if(solar(i) .and. coef%ff_cwn(chan) > ray_min_wvn) then

      ! Calculate layer-independent scattering parameter
      wlm = 10000.0_jprb / coef%ff_cwn(chan)    ! Wavelength in microns

      if (wlm < ray_scs_wlm) then
        ss_param = ray_scs_a1 * wlm ** (ray_scs_b1 + ray_scs_c1 * wlm + ray_scs_d1/wlm)
      else
        ss_param = ray_scs_a2 * wlm ** (ray_scs_b2 + ray_scs_c2 * wlm + ray_scs_d2/wlm)
      endif
      ss_param = ss_param * 0.01_jprb ** 2_jpim * na * z4pi_r / gravity

      rayrad_up_tl(0,:) = 0._jprb
      rayrad_dn_tl(0,:) = 0._jprb

      ! Sum contributions from atmospheric layers
      do lev = 2, nlevels
        lay = lev - 1_jpim

        ! Layer H2O by volume as fraction:
        v_h2o_dry = 0.5_jprb * (profiles_dry(prof)%q(lev-1) + profiles_dry(prof)%q(lev)) * 1.E-06_jprb
        v_h2o_dry_tl = 0.5_jprb * (profiles_dry_tl(prof)%q(lev-1) + profiles_dry_tl(prof)%q(lev)) * 1.E-06_jprb
        if (profiles(prof)%gas_units /= gas_unit_compatibility) then
          v_h2o = v_h2o_dry / (1._jprb + v_h2o_dry)
          v_h2o_tl = v_h2o_dry_tl * v_h2o * ( 1._jprb - v_h2o / (1._jprb + v_h2o_dry)) / v_h2o_dry
        else
          v_h2o = v_h2o_dry
          v_h2o_tl = v_h2o_dry_tl
        endif

        ! Average molar weight of wet air for the layer (kg)
        Mwet = ((1._jprb - v_h2o) * Mair + v_h2o * Mh2o) * 1.E-3_jprb
        Mwet_tl = v_h2o_tl * (Mh2o - Mair) * 1.E-3_jprb

        ! cosine of scattering angle - raytracing%zasat/zasun contain the sine of the angles
        ! NB there are no perturbations around the profile angles, but the raytracing angles
        !    are dependent on atmospheric parameters.
        cossat = 1._jprb - raytracing%zasat(lay, prof) * raytracing%zasat(lay, prof)
        cossol = 1._jprb - raytracing%zasun(lay, prof) * raytracing%zasun(lay, prof)
        cosscata_term1 = sqrt(cossat * cossol)
        cosazi = cos((profiles(prof)%azangle - profiles(prof)%sunazangle)*deg2rad)
        cosscata_term2 = raytracing%zasat(lay, prof) * raytracing%zasun(lay, prof) * cosazi
        solar_src = solar_spectrum(i) * ss_param * &
              (profiles(prof)%p(lev) - profiles(prof)%p(lev-1)) * 100.0_jprb * &
              raytracing%pathsat(lay, prof) / Mwet

        cosscata = - cosscata_term1 - cosscata_term2
        ray_phase = 0.75_jprb * (1.0_jprb + cosscata * cosscata)
        solar_src_updn = solar_src * ray_phase

        cosscata_term1_tl = -1._jprb / cosscata_term1 * &
                            (raytracing_tl%zasat(lay, prof) * raytracing%zasat(lay, prof) * cossol + &
                             raytracing_tl%zasun(lay, prof) * raytracing%zasun(lay, prof) * cossat)
        cosscata_term2_tl = raytracing_tl%zasat(lay, prof) * raytracing%zasun(lay, prof) * cosazi + &
                            raytracing_tl%zasun(lay, prof) * raytracing%zasat(lay, prof) * cosazi

        if (opts%interpolation%lgradp) then
          plev_tl = profiles_tl(prof)%p(lev)
          plev1_tl = profiles_tl(prof)%p(lev-1)
        else
          plev_tl = 0._jprb
          plev1_tl = 0._jprb
        endif

        solar_src_tl = solar_spectrum(i) * ss_param * &
                 ((plev_tl - plev1_tl) * 100.0_jprb * &
                 raytracing%pathsat(lay, prof) / Mwet + &
                 (profiles(prof)%p(lev) - profiles(prof)%p(lev-1)) * 100.0_jprb * &
                 raytracing_tl%pathsat(lay, prof) / Mwet + &
                 (profiles(prof)%p(lev) - profiles(prof)%p(lev-1)) * 100.0_jprb * &
                 raytracing%pathsat(lay, prof) * (-1._jprb) * Mwet_tl / (Mwet * Mwet))

        cosscata_tl = - cosscata_term1_tl - cosscata_term2_tl
        ray_phase_tl = 2._jprb * 0.75_jprb * cosscata_tl * cosscata
        solar_src_updn_tl = solar_src_tl * ray_phase + solar_src * ray_phase_tl

        do ist = 0, ircld%nstream(prof)
          rayrad_up_tl(lay,ist) = rayrad_up_tl(lay-1,ist) + &
                                  transmission_aux_tl%solar_path2%Tau_level(lev-1,ist,i) * solar_src_updn + &
                                  transmission_aux%solar_path2%Tau_level(lev-1,ist,i) * solar_src_updn_tl

          rayrad_dn_tl(lay,ist) = rayrad_dn_tl(lay-1,ist) + &
                                  (solar_src_updn_tl * transmission_aux%solar_path2%Tau_level(lev-1,ist,i) + &
                                  solar_src_updn * transmission_aux_tl%solar_path2%Tau_level(lev-1,ist,i) - &
                                  solar_src_updn * transmission_aux%solar_path2%Tau_level(lev-1,ist,i) * &
                                  3_jpim * transmission_aux_tl%solar_path1%Tau_level(lev-1,ist,i) / &
                                  transmission_aux%solar_path1%Tau_level(lev-1,ist,i)) / &
                                  transmission_aux%solar_path1%Tau_level(lev-1,ist,i) ** 3_jpim
        enddo
      enddo

      ! Add Rayleigh contributions to radiance totals
      do ist = 0, ircld%nstream(prof)
        auxrad_stream_tl%up_solar(:,ist,i) = auxrad_stream_tl%up_solar(:,ist,i) + &
                                             rayrad_up_tl(1:nlayers,ist)

        auxrad_stream_tl%meanrad_up_solar(ist,i) = auxrad_stream_tl%meanrad_up_solar(ist,i) + &
                                                   rayrad_up_tl(iv3lay(i),ist)

        auxrad_stream_tl%down_solar(:,ist,i) = auxrad_stream_tl%down_solar(:,ist,i) + &
                                               rayrad_dn_tl(1:nlayers,ist)

        auxrad_stream_tl%meanrad_down_solar(ist,i) = auxrad_stream_tl%meanrad_down_solar(ist,i) + &
                                                     rayrad_dn_tl(iv3lay(i),ist)
      enddo

      ! Calculate the contribution from the part-layer above the surface

      ! Layer H2O by volume as fraction:
      if (opts%rt_all%use_q2m) then
        v_h2o_dry = 0.5_jprb * (profiles_dry(prof)%q(iv3lev(i)) + profiles_dry(prof)%s2m%q) * 1.E-06_jprb
        v_h2o_dry_tl = 0.5_jprb * (profiles_dry_tl(prof)%q(iv3lev(i)) + profiles_dry_tl(prof)%s2m%q) * 1.E-06_jprb
      else
        v_h2o_dry = profiles_dry(prof)%q(iv3lev(i)) * 1.E-06_jprb
        v_h2o_dry_tl = profiles_dry_tl(prof)%q(iv3lev(i)) * 1.E-06_jprb
      endif
      if (profiles(prof)%gas_units /= gas_unit_compatibility) then
        v_h2o = v_h2o_dry / (1._jprb + v_h2o_dry)
        v_h2o_tl = v_h2o_dry_tl * v_h2o * ( 1._jprb - v_h2o / (1._jprb + v_h2o_dry)) / v_h2o_dry
      else
        v_h2o = v_h2o_dry
        v_h2o_tl = v_h2o_dry_tl
      endif

      ! Average molar weight of wet air for the layer (kg):
      Mwet = ((1._jprb - v_h2o) * Mair + v_h2o * Mh2o) * 1.E-3_jprb
      Mwet_tl = v_h2o_tl * (Mh2o - Mair) * 1.E-3_jprb

      ! cosine of scattering angle
      ! NB there are no pertubations around profile angles
      cossat = 1._jprb - raytracing%zasat(iv2lay(i), prof) * raytracing%zasat(iv2lay(i), prof)
      cossol = 1._jprb - raytracing%zasun(iv2lay(i), prof) * raytracing%zasun(iv2lay(i), prof)
      cosscata_term1 = sqrt(cossat * cossol)
      cosazi = cos((profiles(prof)%azangle - profiles(prof)%sunazangle)*deg2rad)
      cosscata_term2 = raytracing%zasat(iv2lay(i), prof) * raytracing%zasun(iv2lay(i), prof) * cosazi
      solar_src = solar_spectrum(i) * & ! mW m^-2 (cm^-1)^-1
            abs(profiles(prof)%s2m%p - profiles(prof)%p(iv3lev(i))) * 100.0_jprb * &  ! convert hPa to Pa
            raytracing%pathsat(iv2lay(i), prof) * ss_param / Mwet

      cosscata = - cosscata_term1 - cosscata_term2
      ray_phase = 0.75_jprb * (1.0_jprb + cosscata * cosscata)
      solar_src_updn = solar_src * ray_phase

      cosscata_term1_tl = -1._jprb / cosscata_term1 * &
                      (raytracing_tl%zasat(iv2lay(i), prof) * raytracing%zasat(iv2lay(i), prof) * cossol + &
                      raytracing_tl%zasun(iv2lay(i), prof) * raytracing%zasun(iv2lay(i), prof) * cossat)

      cosscata_term2_tl = raytracing_tl%zasat(iv2lay(i), prof) * raytracing%zasun(iv2lay(i), prof) * cosazi + &
                          raytracing_tl%zasun(iv2lay(i), prof) * raytracing%zasat(iv2lay(i), prof) * cosazi

      if (opts%interpolation%lgradp) then
        plev_tl = profiles_tl(prof)%p(iv3lev(i))
      else
        plev_tl = 0._jprb
      endif

      if (profiles(prof)%s2m%p >= profiles(prof)%p(iv3lev(i))) then
        dp_tl = profiles_tl(prof)%s2m%p - plev_tl
      else
        dp_tl = plev_tl - profiles_tl(prof)%s2m%p
      endif

      solar_src_tl = solar_spectrum(i) * ss_param * &
                     (dp_tl * 100.0_jprb * &
                     raytracing%pathsat(iv2lay(i), prof) / Mwet + &
                     abs(profiles(prof)%s2m%p - profiles(prof)%p(iv3lev(i))) * 100.0_jprb * &
                     raytracing_tl%pathsat(iv2lay(i), prof) / Mwet + &
                     abs(profiles(prof)%s2m%p - profiles(prof)%p(iv3lev(i))) * 100.0_jprb * &
                     raytracing%pathsat(iv2lay(i), prof) * (-1._jprb) * Mwet_tl / (Mwet * Mwet))

      cosscata_tl = - cosscata_term1_tl - cosscata_term2_tl
      ray_phase_tl = 2._jprb * 0.75_jprb * cosscata_tl * cosscata
      solar_src_updn_tl = solar_src_tl * ray_phase + solar_src * ray_phase_tl

      ! Add near-surface contributions to the total radiances
      do ist = 0, ircld%nstream(prof)
        auxrad_stream_tl%meanrad_up_solar(ist,i) = auxrad_stream_tl%meanrad_up_solar(ist,i) + &
                                  transmission_aux_tl%solar_path2%Tau_level(iv3lev(i),ist,i) * solar_src_updn + &
                                  transmission_aux%solar_path2%Tau_level(iv3lev(i),ist,i) * solar_src_updn_tl

        auxrad_stream_tl%meanrad_down_solar(ist,i) = auxrad_stream_tl%meanrad_down_solar(ist,i) + &
                                  (solar_src_updn_tl * transmission_aux%solar_path2%Tau_level(iv3lev(i),ist,i) + &
                                  solar_src_updn * transmission_aux_tl%solar_path2%Tau_level(iv3lev(i),ist,i) - &
                                  solar_src_updn * transmission_aux%solar_path2%Tau_level(iv3lev(i),ist,i) * &
                                  3_jpim * transmission_aux_tl%solar_path1%Tau_level(iv3lev(i),ist,i) / &
                                  transmission_aux%solar_path1%Tau_level(iv3lev(i),ist,i)) / &
                                  transmission_aux%solar_path1%Tau_level(iv3lev(i),ist,i) ** 3_jpim
      enddo

    endif
  enddo

end subroutine solar_rayleigh_tl

subroutine solar_surface_contribution_tl(transmission_aux, transmission_aux_tl, ircld, &
                                             reflectance, reflectance_tl, refl_norm, chanprof, solar, &
                                             solar_spectrum, narray, auxrad_stream_tl, rad_tl)

  use rttov_types, only : rttov_chanprof, ircld_type, transmission_type_aux, radiance_type, radiance_aux
  use parkind1, only : jpim, jprb, jplm 

  Implicit None

  type(transmission_type_aux), intent(in)    :: transmission_aux, transmission_aux_tl
  type(ircld_type),     intent(in)           :: ircld
  real(jprb), intent(in)                     :: reflectance(:), reflectance_tl(:) ! surface reflectance
  real(jprb), intent(in)                     :: refl_norm(:)
  type(rttov_chanprof), intent(in)           :: chanprof(:)
  logical(jplm), intent(in)                  :: solar(:)
  real(jprb), intent(in)                     :: solar_spectrum(:)
  integer(jpim), intent(in)                  :: narray(:)

  type(radiance_type), intent(inout)         :: rad_tl
  type(radiance_aux), intent(inout)          :: auxrad_stream_tl

  integer(jpim) :: prof
  integer(jpim) :: nlevels, nlayers, nstreams, nchannels
  integer(jpim) :: i, ist
  real(jprb)    :: temp_tl(0:narray(2)) ! ist

! unpack
  nlevels = narray(1)
  nlayers = nlevels - 1_jpim
  nchannels = narray(3)  

  Do i = 1, nchannels
     if(solar(i)) then
        prof = chanprof(i)%prof
        nstreams = ircld%nstream(prof)

        do ist = 0, nstreams
          temp_tl(ist) = solar_spectrum(i) * refl_norm(i) * &
                      (reflectance_tl(i) * transmission_aux%solar_path2%Tau_surf(ist,i) + &
                      reflectance(i) * transmission_aux_tl%solar_path2%Tau_surf(ist,i))
        enddo

        rad_tl%clear(i) = rad_tl%clear(i) + temp_tl(0)
        auxrad_stream_tl%cloudy(1:nstreams,i) = auxrad_stream_tl%cloudy(1:nstreams,i) + temp_tl(1:nstreams)
     endif
  enddo
end subroutine solar_surface_contribution_tl

end subroutine rttov_integrate_tl
