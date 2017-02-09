subroutine rttov_integrate(addcosmic, opts, maxnstreams, chanprof, &! in
                           emissivity, reflectance, refl_norm, thermrefl, do_lambertian, &! in
                           thermal, dothermal, solar, dosolar, solar_spectrum, &!in
                           transmission_aux, transmission_scatt_ir, &!in
                           profiles, profiles_dry, aux_prof, coef, raytracing, ircld, &! in
                           rad, rad2, auxrad, auxrad_stream)      ! inout

! Description:
! To perform integration of radiative transfer equation
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
! ECMWF Research Dept. Tech. Memo. 425 (available from the librarian at ECMWF).
!
! Matricardi, M. 2005 The inclusion of aerosols and clouds in RTIASI,the
! ECMWF radiative transfer model for the infrared atmospheric sounding
! interferometer.
! ECMWF Research Dept. Tech. Memo. 474 (available from the librarian at ECMWF).
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
!  1.4     03/09/2004   Mods. for Vectorisation (D Salmond ECMWF & BCarruthers, Cray)
!  1.5     28/02/2005   Further mods to vectorisation (D Dent)
!  1.6     01/06/2005   Marco Matricardi (ECMWF):
!             --        IASI capability added.
!             --        Linear in tau approximation for RT equation introduced.
!             --        Solar radiation introduced for IASI and AIRS.
!  1.7     01/06/2006   Marco Matricardi (ECMWF):
!                       Multiple scattering in the infrared introduced for
!                       water clouds,cirrus clouds and aerosols.
!  1.8     26/01/2007   Removed polarisation (R Saunders)
!  1.9     08/03/2007   Reintroduced overcast cloud (R Saunders)
!  1.10    13/07/2007   Added extra variables requested by P. Watts(RSaunders)
!  1.11    15/11/2007   Changed skin to surfair for overcast radiance calc (R Saunders)
!  1.12    27/11/2007   Optimised for NEC/IBM (D. Salmond)
!  1.13    21/12/2007   Added polarimetric option (R. Saunders)
!  1.14    27/02/2009   Profile levels to include ToA. Distinguish between
!                       layer arrays and level arrays - size, index labels,
!                       looping (P. Rayer)
!  1.15    17/04/2009   Model-top brought into calculation of source function (P.Rayer)
!  1.16    29/06/2009   Top-layer brought into layer looping to shorten code (P.Rayer)
!  1.17    03/11/2009   Transmittances on levels (A Geer)
!  1.18    02/12/2009   Introduced principal component capability. Pathsat, Pathsun and
!                       related quantities are now layer arrays (Marco Matricardi).
!  1.19    05/07/2010   Remove addsolar flag from profiles structure (J Hocking)
!  1.20    14/10/2010   Remove rt8_mode (J Hocking)
!  1.21    14/12/2010   Use traj0_sta%sun array to flag channels for which solar calculations
!                       should be performed (J Hocking)
!  2.0     14/12/2011   Re-written (D Rundle)
!  2.01       03/2012   Code cleaning Pascal Brunel, Philippe Marguinaud
!  2.02    08/11/2012   Added option of MW lambertian approximation (R Saunders)
!!
!
! Code Description:
!   Language:           Fortran 90.
!   Software Standards: "European Standards for Writing and
!     Documenting Exchangeable Fortran 90 Code".

  use parkind1, Only : jpim, jprb, jplm

  use rttov_types, Only : rttov_chanprof, rttov_coef, rttov_options, profile_type, profile_aux, &
                          transmission_type_aux, transmission_scatt_ir_type, radiance_type, radiance2_type, &
                          ircld_type, raytracing_type, radiance_aux, rttov_emissivity, rttov_reflectance

!INTF_OFF
  use rttov_const, Only : sensor_id_po, min_od, min_tau, pi_r, z4pi_r, &
                          deg2rad, gravity, na, Mh2o, Mair, gas_unit_compatibility, &
                          overcast_albedo_wvn, overcast_albedo1, overcast_albedo2, &
                          ray_scs_wlm, ray_scs_a1, ray_scs_b1, ray_scs_c1, &
                          ray_scs_d1, ray_scs_a2, ray_scs_b2, ray_scs_c2, ray_scs_d2, ray_min_wvn

  use yomhook, Only : LHOOK, DR_HOOK
!INTF_ON

  Implicit None

!subroutine arguments:

  logical(jplm),   intent(in)                   :: addcosmic   ! switch for adding cosmic background
  type(rttov_options),  intent(in)              :: opts        ! options structure
  integer(jpim),   intent(in)                   :: maxnstreams                !
  type(rttov_chanprof), intent(in)              :: chanprof(:)     ! Channel indices
  type(profile_type),   intent(in)              :: profiles(:)     ! Profiles
  type(profile_type),   intent(in)              :: profiles_dry(:) ! Gas profiles (ppmv dry)
  type(rttov_emissivity),  intent(in), optional :: emissivity(size(chanprof))   ! surface emissivity
  type(rttov_reflectance), intent(in), optional :: reflectance(size(chanprof))  ! surface solar reflectance
  real(jprb),      intent(in)                   :: refl_norm(size(chanprof))
  real(jprb),      intent(in)                   :: thermrefl(size(chanprof))    ! surface thermal reflectance
  logical(jplm),   intent(in)                   :: do_lambertian(size(chanprof))
  logical(jplm),   intent(in)                   :: thermal(size(chanprof))
  logical(jplm),   intent(in)                   :: dothermal
  logical(jplm),   intent(in)                   :: solar(size(chanprof))
  logical(jplm),   intent(in)                   :: dosolar
  real(jprb),      intent(in)                   :: solar_spectrum(size(chanprof)) ! ToA solar irradiance
  type(ircld_type),     intent(in)              :: ircld
  type(raytracing_type),intent(in)              :: raytracing
  type(transmission_type_aux), intent(in)       :: transmission_aux               ! transmittances and single-layer od
  type(transmission_scatt_ir_type), intent(in)  :: transmission_scatt_ir
  type(profile_aux) ,   intent(in)              :: aux_prof ! auxillary profiles info.
  type(rttov_coef),     intent(in)              :: coef
  type(radiance_aux),   intent(inout)           :: auxrad_stream
  type(radiance_type),  intent(inout)           :: rad    ! radiances (mw/cm-1/ster/sq.m) and BTs
  type(radiance2_type), intent(inout), optional :: rad2
  type(radiance_aux),   intent(inout)           :: auxrad ! auxillary radiances

!INTF_END

!module variables: 

  integer(jpim)              :: i, lev, ist, isti, lay, nchannels, nlayers, nlevels, iprof, totnstreams ! counter variables
  real(jprb)                 :: refl, refl_norm_scat
  logical(jplm)              :: keyradonly ! flag to indicate only calculate primary radiance outputs

  integer(jpim), allocatable :: nstreams(:)
  integer(jpim), allocatable :: iv2lev(:), iv2lay(:), iv3lev(:), iv3lay(:)
  integer(jpim), allocatable :: pol_id(:)       ! polarisation index  

  real(jprb), allocatable    :: cfraction(:), pfraction(:)        ! cloud fraction
  logical(jplm), allocatable :: sateqsun(:,:) ! True where the solar zenith angle equal to observation angle

! Define macros for commonly used variables
#define prof chanprof(i)%Prof
#define chan chanprof(i)%Chan

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
#define fac1 transmission_aux%Fac1(lay,ist,i)
#define fac2_thermal_path1 transmission_aux%thermal_path1%fac2(lay+1,ist,i)
#define fac2_solar_path1 transmission_aux%solar_path1%fac2(lay+1,ist,i)
#define surf_fac transmission_aux%Surf_fac(ist,i)

#include "rttov_calcrad.interface"

  REAL(JPRB) :: ZHOOK_HANDLE

!- End of header ------------------------------------------------------

if (LHOOK) CALL DR_HOOK('RTTOV_INTEGRATE',0_jpim,ZHOOK_HANDLE)

!-------------------------------------------------------------------------------
!0. Initialise useful variables
!-------------------------------------------------------------------------------

call init_rttov_integrate_mod_vars(chanprof, profiles, ircld, aux_prof, coef, raytracing)

totnstreams = maxnstreams

keyradonly = opts%rt_ir%addaerosl .or. opts%rt_ir%pc%addpc

if (present(rad2)) then
  rad2%down = 0.0_jprb
  if (keyradonly) then
    rad2%up = 0.0_jprb
    rad2%upclear = 0.0_jprb
    rad2%dnclear = 0.0_jprb
    rad2%refldnclear = 0.0_jprb
    rad2%surf = 0.0_jprb
  endif
endif

if (keyradonly) rad%overcast = 0.0_jprb

rad%clear = 0.0_jprb
auxrad_stream%cloudy = 0.0_jprb
if (dosolar) then
  auxrad_stream%up_solar = 0.0_jprb
  auxrad_stream%meanrad_up_solar = 0.0_jprb
  auxrad_stream%down_solar = 0.0_jprb
  auxrad_stream%meanrad_down_solar = 0.0_jprb
endif

!-------------------------------------------------------------------------------
!1. calculate layer radiances
!-------------------------------------------------------------------------------

if (dothermal) call rttov_calcrad(addcosmic, chanprof, profiles,  coef, &!in
                   auxrad) !out

!-------------------------------------------------------------------------------
!2. calculate atmospheric contribution from layers
!-------------------------------------------------------------------------------
 
if (dothermal) call calc_atmospheric_radiance(thermal, transmission_aux, auxrad, &! in
                               auxrad_stream) ! out

!---Scattering of the solar beam--------------------------------------------------
if((opts%rt_ir%addaerosl .OR. opts%rt_ir%addclouds) .AND. dosolar) &
     call  solar_scattering_air(transmission_aux, solar, solar_spectrum, reflectance%refl_out, &
                                raytracing, chanprof, transmission_scatt_ir, auxrad_stream)

!-------------------------------------------------------------------------------
!2a calculate near-surface layer contribution
!-------------------------------------------------------------------------------

if (dothermal) call calc_near_surf_contribution(thermal, transmission_aux, auxrad, auxrad_stream)

!---Scattering of the solar beam---------------------------
if((opts%rt_ir%addaerosl .OR. opts%rt_ir%addclouds) .AND. dosolar) &
     call solar_scattering_near_surf(transmission_aux, solar, reflectance%refl_out, chanprof, auxrad_stream)

!-------------------------------------------------------------------------------
!2b calculate clear-sky Rayleigh scattering contribution
!-------------------------------------------------------------------------------

if (dosolar) &
  call solar_rayleigh(coef, opts, solar, solar_spectrum, chanprof, raytracing, profiles, &
                      profiles_dry, transmission_aux, auxrad_stream)

!-------------------------------------------------------------------------------
! Add thermal and solar atmospheric contributions to the clear and cloudy streams
!-------------------------------------------------------------------------------
do i = 1, nchannels
  ist = 0_jpim
  if (thermal(i))then
    if( do_lambertian(i) )then   !Lambertian reflected downwelling option
      rad%clear(i) = rad%clear(i) + auxrad_stream%meanrad_up(ist,i) + &
                  auxrad_stream%meanrad_down(ist,i) * thermrefl(i) * &
                  tau_surf_p * tau_surf
    else                        !Specular reflected downwelling option
      rad%clear(i) = rad%clear(i) + auxrad_stream%meanrad_up(ist,i) + &
                  auxrad_stream%meanrad_down(ist,i) * thermrefl(i) * tau_surf**2_jpim
    endif
  endif

  if (solar(i)) then
    ! For surface-reflected component, use input BRDF - but this isn't really very good (especially for sunglint).
    ! Normalisation is cos(sat_zen_angle) for the downward scattered beam
    refl_norm_scat = COS(profiles(prof)%zenangle * deg2rad)

    rad%clear(i) = rad%clear(i) + auxrad_stream%meanrad_up_solar(ist,i) + &
                   auxrad_stream%meanrad_down_solar(ist,i) * &
                   reflectance(i)%refl_out * refl_norm_scat * &
                   transmission_aux%solar_path1%Tau_surf(ist,i)**2_jpim
  endif

  ! These lines replace the upward radiances from the level at the bottom of the layer
  ! containing the surface with the calculated surface->ToA radiance (only for stream 0)

  if (thermal(i)) auxrad_stream%up(iv2lay(i),ist,i) = auxrad_stream%meanrad_up(ist,i)
  if (solar(i)) auxrad_stream%up_solar(iv2lay(i),ist,i) = auxrad_stream%meanrad_up_solar(ist,i)

  ! Same as above, but for each cloud stream
  if (thermal(i)) then
    do ist = 1, nstreams(i)
      if( do_lambertian(i) )then
        auxrad_stream%cloudy(ist,i) = &
                            auxrad_stream%cloudy(ist,i) + auxrad_stream%meanrad_up(ist,i) + &
                            auxrad_stream%meanrad_down(ist,i) * thermrefl(i) * &
                            tau_surf_p * tau_surf
      else
        auxrad_stream%cloudy(ist,i) = &
                            auxrad_stream%cloudy(ist,i) + auxrad_stream%meanrad_up(ist,i) + &
                            auxrad_stream%meanrad_down(ist,i) * thermrefl(i) * tau_surf**2_jpim
      endif
    enddo
  endif

  if (solar(i)) then
    do ist = 1, nstreams(i)
      auxrad_stream%cloudy(ist,i) = auxrad_stream%cloudy(ist,i) + &
                                    auxrad_stream%meanrad_up_solar(ist,i) + &
                                    auxrad_stream%meanrad_down_solar(ist,i) * &
                                    reflectance(i)%refl_out * refl_norm_scat * &
                                    transmission_aux%solar_path1%Tau_surf(ist,i)**2_jpim
    enddo
  endif
enddo


!-------------------------------------------------------------------------------
!3. calculate surface emission contribution
!-------------------------------------------------------------------------------
if (dothermal) then
!cdir nodep
  do i = 1, nchannels
    if (thermal(i)) then
      ist = 0_jpim
      ! clear sky radiance without reflection term
      rad%clear(i) = rad%clear(i) + auxrad%skin(i) * emissivity(i)%emis_out * tau_surf

!cdir nodep
      do ist = 1, nstreams(i)
        auxrad_stream%cloudy(ist,i) = auxrad_stream%cloudy(ist,i) + &
                                      auxrad%skin(i) * emissivity(i)%emis_out * tau_surf
      enddo
    endif
  enddo
endif


!-------------------------------------------------------------------------------
!4. Add solar surface reflected contribution
!-------------------------------------------------------------------------------

if (dosolar) &
  call solar_surface_contribution(transmission_aux, solar, solar_spectrum, reflectance%refl_out, refl_norm, &
                                  auxrad_stream, rad)

!-------------------------------------------------------------------------------
!5. cosmic temperature correction
!-------------------------------------------------------------------------------
if (addcosmic) then
  ist = 0
  do i = 1, nchannels
    if (do_lambertian(i)) then
       rad%clear(i) = rad%clear(i) + thermrefl(i) * &
                      tau_surf_p * &
                      tau_surf * &
                      auxrad%cosmic(i)
    else
       rad%clear(i) = rad%clear(i) + thermrefl(i) * &
                      transmission_aux%thermal_path1%Tau_surf(ist,i)**2_jpim * &
                      auxrad%cosmic(i)
    endif
  enddo
endif

!-------------------------------------------------------------------------------
! Calculate secondary radiances
!-------------------------------------------------------------------------------

! These outputs only include thermal contributions (no solar) and are only calculated
! for non-scattering and non-PC simulations if the rad2 parameter is present

if (present(rad2) .and. dothermal .and. .not. keyradonly) then

!cdir nodep
  do i = 1, nchannels
    if (thermal(i)) then

      ist = 0_jpim

      ! Note the surface layer value in auxrad_stream%up has been modified above
      rad2%up(:,i) = auxrad_stream%up(:,ist,i)

      do lay = 2, nlayers
        if (tau_level > min_tau ) Then
            rad2%down(lay,i) = auxrad_stream%down(lay,ist,i) * tau_level
        else
            rad2%down(lay,i) = rad2%down(lay-1,i)
        endif
        rad2%surf(lay,i) = auxrad%air(lay+1,i)
      enddo

      if (tau_surf > min_tau ) Then
        rad2%down(iv2lay(i),i) = auxrad_stream%meanrad_down(ist,i) * tau_surf
      else
        rad2%down(iv2lay(i),i) = rad2%down(iv3lay(i),i)
      endif
      rad2%surf(iv2lay(i),i) = auxrad%skin(i)

      ! clear-sky upwelling radiance (without surface reflected term)
      rad2%upclear(i) = auxrad_stream%meanrad_up(ist,i) + &
                        auxrad%skin(i) * emissivity(i)%emis_out * tau_surf


      ! clear-sky downwelling radiance
      rad2%dnclear(i) = auxrad_stream%meanrad_down(ist,i) * tau_surf

      ! reflected clear-sky downwelling radiance
      if( do_lambertian(i) )then   !Lambertian reflected downwelling option
        rad2%refldnclear(i) = rad2%dnclear(i) * thermrefl(i) * tau_surf_p
      else                        ! Specular reflection
        rad2%refldnclear(i) = rad2%dnclear(i) * thermrefl(i) * tau_surf
      endif
    endif
  enddo

endif


!-------------------------------------------------------------------------------
!6. calculate overcast radiances
!-------------------------------------------------------------------------------

! Overcast radiances only calculated for non-scattering, non-PC simulations

if (.not. keyradonly) then

  ! rad%overcast includes solar contribution ONLY for pure-solar channels because assumed
  !   cloud emissivity is 1.0 (i.e. no reflection) for all thermal channels
  ist = 0_jpim
  do i = 1, nchannels
    if (thermal(i)) then
      do lay = 1, nlayers
        lev = lay + 1
        ! overcast radiances at given cloud top
        rad%overcast(lay,i) = auxrad_stream%up(lay,ist,i) + auxrad%air(lev,i) * tau_level
      enddo
    else if (solar(i)) then
      ! Very crude model: assumes clouds are Lambertian reflectors with fixed albedo
      ! Use input cloud top BRDF if user has supplied it, otherwise use default BRDF
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
        rad%overcast(lay,i) = solar_spectrum(i) * refl / raytracing%pathsun(lay, prof) * &
                              transmission_aux%solar_path2%Tau_level(lev,ist,i) + &
                              auxrad_stream%up_solar(lay,ist,i) + &
                              auxrad_stream%down_solar(lay,ist,i) * refl / raytracing%pathsat(lay, prof) * &
                              transmission_aux%solar_path1%Tau_level(lev,ist,i) ** 2_jpim
      enddo
    else
      rad%overcast(:,i) = 0._jprb
    endif
  enddo

  if (dothermal) then
    ! Add surface component to thermal overcast radiances
    do i = 1, nchannels
      if (thermal(i)) then
        lay = iv2lay(i)
        rad%overcast(lay,i) = auxrad_stream%up(lay,ist,i) + tau_surf * auxrad%surfair(i)
      endif
    enddo
  endif
endif

!-------------------------------------------------------------------------------
! 7. Calculate total radiance
!-------------------------------------------------------------------------------
! The simple cloudy scheme is not applied to aerosol-affected radiances
if (opts%rt_ir%addaerosl .OR. opts%rt_ir%addclouds) then
!---------------------------------------------------
! Calculate complex cloudy radiances
!---------------------------------------------------
   rad%cloudy = 0._jprb
   do i = 1, nchannels
      do ist = 1, nstreams(i)
         rad%cloudy(i) = rad%cloudy(i) + &
                         auxrad_stream%cloudy(ist,i) * &
                         (ircld%xstr(ist+1,prof) - ircld%xstr(ist,prof))
      enddo
      rad%cloudy(i) = rad%cloudy(i) + rad%clear(i) * ircld%XSTRCLR(prof)
   enddo

   rad%total(1:nchannels) = rad%cloudy(1:nchannels)
else
!---------------------------------------------------
! Calculate total radiance (clear case/simple cloud)
!---------------------------------------------------

  if(opts%rt_ir%pc%addpc) then
    rad%total(1:nchannels) = rad%clear(1:nchannels)
  else

    ! Interpolate to given cloud-top pressures
    do i = 1, nchannels
      lay = aux_prof%s(prof)%nearestlev_ctp - 1
      rad%cloudy(i) = rad%overcast(lay,i) * &
                      (1.0_JPRB - aux_prof%s(prof)%pfraction_ctp) + &
                      rad%overcast(lay-1,i) * (aux_prof%s(prof)%pfraction_ctp)
    enddo

    rad%total(1:nchannels) = rad%clear(1:nchannels) + &
                             cfraction(1:nchannels) * &
                            (rad%cloudy(1:nchannels) - rad%clear(1:nchannels) )
  endif
endif


! deallocate module arrays - DAR: can move these to a more useful array if they are commonly used?
deallocate(nstreams, iv2lev, iv2lay, iv3lev, iv3lay, pol_id, cfraction, pfraction, sateqsun)

if (LHOOK) CALL DR_HOOK('RTTOV_INTEGRATE',1_jpim,ZHOOK_HANDLE)

contains

! DAR: This subroutine calculates the individual layer clear-sky radiances and then does a cumulative sum to determine the 
!      radiance observed from the top of atmosphere to a particular layer.
! DAR: As expected, this routine consumes most of the time (loops over channels, streams and levels)
!      and has been most heavily optimised (for IBM only so far, Intel shows neutral impact - will look into this)
!      I have taken out the code that switches the order of the loops on the NEC and will let MF test this impact       
!      See subroutine comments for more details of individual changes.

subroutine calc_atmospheric_radiance(thermal, transmission_aux, auxrad, auxrad_stream)

  Implicit None

  logical(jplm), intent(in)                     :: thermal(:)
  type(transmission_type_aux), intent(in)       :: transmission_aux
  type(radiance_aux), intent(in)                :: auxrad
  type(radiance_aux), intent(inout)             :: auxrad_stream


! DAR: fac and fac2 contain real 1 or 0 depending on whether a calculation should be performed or not.
!      IBM performance was suffering as a result of doing lots of branching (mispredicts?) so these arrays are populated in advance
!      and the calculation is performed regardless.

! DAR: I've removed a lot of the 'temporary' variables that were eventually summed together because they were taking up a signficant
!      amount of time on the IBM. Now the big sum is done with all the variables in full because the 16 prefetch streams can handle 
!      this (POWER6) - maybe will have to split this up for POWER7 (only 12). Will do more testing on Intel with vtune to check 
!      performance impact.

! DAR: I've also removed od_singlelayer_r and added it to transmission_aux so it can be reused the TL/AD/K code when ready 
!      and doesn't have to be recalculated - This should be faster...
  
! DAR: Slow code path - may be faster on intel?

! If negative transmittance then do this block

  if(transmission_aux%anynegtau .gt. 0.0_jprb) then
   do i = 1, nchannels
      if (thermal(i)) then
        do ist = 0, nstreams(i)
          do lay = 1, nlayers
            isti = ircld%icldarr(ist,lay,prof)
! DAR: As explained above, every difference is calculated
!      explicitly as this is apparently faster on the IBM. Also, changed indices on auxrad_stream%up/down to be consistent with
!      every other variable.
          auxrad_stream%up(lay,ist,i) = fac1 * &
                (dtau_lay * &
                (auxrad%air(lay,i) + &
                daux_lay * &
                od_singlelayer_r) - &
                daux_lay * &
                tau_level)

          auxrad_stream%down(lay,ist,i) = &
                fac1 * &
                fac2_thermal_path1 * &
                ((dtau_lay * &
                (auxrad%air(1,i) - &
                daux_lay * &
                od_singlelayer_r) * &
                (tau_level_r * &
                tau_layer_r)) + &
                daux_lay * &
                tau_level_r)
          enddo
        enddo
      endif
   enddo

! DAR: This is slow. Hopefully it'll never run and can be deleted
   do  i = 1, nchannels
      if (thermal(i)) then
        do ist = 0, nstreams(i)
          do lay = 1,nlayers
              if(tau_level < 0.0_jprb .and. &
                tau_layer >= 0._jprb) then
                auxrad_stream%up(lay,ist,i) = (0.5_JPRB * &
                    (auxrad%air(lay+1,i) + &
                    auxrad%air(lay,i))) * &
                    dtau_lay
                auxrad_stream%down(lay,ist, i) = 0.0_jprb
              endif
          enddo
        enddo
      endif
   enddo

   do i = 1,nchannels
      if (thermal(i)) then
        do ist = 0, nstreams(i)
          do lay = 2, nlayers
              auxrad_stream%up(lay,ist,i) = auxrad_stream%up(lay, ist,i) + auxrad_stream%up(lay-1,ist,i)
              auxrad_stream%down(lay,ist,i) = auxrad_stream%down(lay,ist,i) + auxrad_stream%down(lay-1,ist,i)
          enddo
        enddo
      endif
   enddo
   
else ! DAR: fast code (on IBM) when no special cases. Note cumulative sum is done at same time as it's significantly quicker than 
   !      doing it later. This code is still very hard to read though.
! Positive transmittance so do this block
   do i = 1, nchannels
      if (thermal(i)) then
        do ist = 0, nstreams(i)
          lay = 1
          isti = ircld%icldarr(ist,lay,prof)
          auxrad_stream%up(lay,ist,i) = &
              fac1 * &
              (dtau_lay * &
              (auxrad%air(lay,i) + &
              daux_lay * &
              od_singlelayer_r) - &
              daux_lay * &
              tau_level)
         if (do_lambertian(i)) then        !Lambertian reflected downwelling option no linear-in-tau
           auxrad_stream%down(lay,ist,i) = &
             fac1 * &
             fac2_thermal_path1 * &
             0.5_jprb * (auxrad%air(lay,i) + auxrad%air(lay+1,i)) * &
             ((tau_layer_p) - (tau_level_p)) *  &
             ((tau_level_p_r) * (tau_layer_p_r)) 
         else                               !Specular reflected downwelling with linear-in-tau
          auxrad_stream%down(lay,ist,i) = &
              fac1 * &
              fac2_thermal_path1 * &
              ((dtau_lay * &
              (auxrad%air(lay,i) - &
              daux_lay * &
              od_singlelayer_r) * &
              (tau_level_r * &
              tau_layer_r)) + &
              daux_lay * &
              tau_level_r)          
          endif

          DO lay = 2, nlayers
            isti = ircld%icldarr(ist,lay,prof)
            auxrad_stream%up(lay,ist,i) = auxrad_stream%up(lay-1, ist, i) + &
              fac1 * &
              (dtau_lay * &
              (auxrad%air(lay,i) + &
              daux_lay * &
              od_singlelayer_r) - &
              daux_lay * &
              tau_level)
            IF (do_lambertian(i)) THEN        !Lambertian reflected downwelling option no linear-in-tau
              auxrad_stream%down(lay,ist,i) = auxrad_stream%down(lay-1, ist, i) + &
                fac1 * &
                fac2_thermal_path1 * &
                0.5_jprb * (auxrad%air(lay,i) + auxrad%air(lay+1,i)) * &
                ((tau_layer_p) - (tau_level_p)) *  &
                ((tau_level_p_r) * (tau_layer_p_r))
            ELSE                                 !Specular reflected downwelling with linear-in-tau
              auxrad_stream%down(lay,ist,i) = &
                auxrad_stream%down(lay-1, ist, i) + &
                fac1 * &
                fac2_thermal_path1 * &
                ((dtau_lay * &
                (auxrad%air(lay,i) - &
                daux_lay * &
                od_singlelayer_r) * &
                (tau_level_r * &
                tau_layer_r)) + &
                daux_lay * &
                tau_level_r)
           endif
          enddo
        enddo
      endif
   enddo
endif

do i=1,nchannels
  if(thermal(i)) then
    do ist = 0, nstreams(i)
       do lay = 1, nlayers
         auxrad_stream%down_ref(lay,ist,i) = auxrad_stream%down(lay,ist,i)
         auxrad_stream%down(lay,ist,i) = max(auxrad_stream%down_ref(lay,ist,i), 0.0_jprb)
       enddo
    enddo
  end if
enddo

end subroutine calc_atmospheric_radiance

! DAR: This is the next biggest consumer of CPU time and should be looked at next. The trouble seems to be that you have to change a
!      lot of non-consecutive data. Maybe the RTTOV 9 comments can be removed

subroutine calc_near_surf_contribution(thermal, transmission_aux, auxrad, & !in
                                       auxrad_stream) !inout
  Implicit None

  logical(jplm), intent(in)               :: thermal(:)
  type(transmission_type_aux), intent(in) :: transmission_aux
  type(radiance_aux), intent(in)          :: auxrad
  type(radiance_aux),   intent(inout)     :: auxrad_stream

#define B1_3  auxrad%air(lev,i) * (tau_level  - tau_surf)
#define B2_3 (auxrad%surfair(i) - auxrad%air(lev,i)) * tau_surf
#define B3_3 (auxrad%surfair(i) - auxrad%air(lev,i)) * (tau_level - tau_surf) * (transmission_aux%thermal_path1%od_sfrac_r(ist,i))
#define B1_L ((tau_level_p)  - (tau_surf_p))
#define B2_L (auxrad%surfair(i) + auxrad%air(lev,i))
#define B3_L ((tau_level_p_r) * (tau_surf_p_r))

  Do i = 1, nchannels
    if (thermal(i)) then
      lay = iv3lay(i)
      lev = iv3lev(i)

      do ist = 0, nstreams(i)

          if(tau_surf < 0._JPRB) then 
            auxrad_stream%meanrad_down(ist,i) = 0.0_JPRB

            if(tau_level >= 0._JPRB) then
                auxrad_stream%meanrad_up(ist,i) = &
                  (0.5_JPRB * &
                  (auxrad%surfair(i) + &
                  auxrad%air(lev,i))) * &
                  (tau_level - &
                  tau_surf)
            endif

          else !tau_surf >= 0

            if(transmission_aux%thermal_path1%od_sfrac(ist,i) < min_od .OR. &
              ((tau_level - tau_surf) < min_od)) THEN
                ! small optical depth or optical depth change set radiance to zero
                auxrad_stream%meanrad_down(ist,i) = 0.0_JPRB
                auxrad_stream%meanrad_up(ist,i) = 0.0_JPRB
            else
                ! compute linear in tau radiances
                auxrad_stream%meanrad_up(ist,i) = &
                  B1_3 - &
                  B2_3 + &
                  B3_3
                 if (do_lambertian(i)) then           !Lambertian reflected downwelling option
                   auxrad_stream%meanrad_down(ist,i) = &
                                surf_fac * 0.5_jprb * &
                                B2_L * &
                                B1_L * &
                                B3_L

                 else                                ! Specular reflection option
                   auxrad_stream%meanrad_down(ist,i) = &
                                surf_fac * & 
                                (B1_3 - &
                                B3_3) * &
                                tau_level_r * tau_surf_r + &
                                B2_3 * tau_surf_r**2
                 endif
            endif
          endif

          if (pol_id(i) >= 6_jpim ) then
            auxrad_stream%meanrad_up(ist,i) = 0.0_jprb
          else
            auxrad_stream%meanrad_up(ist,i) = auxrad_stream%meanrad_up(ist,i) + auxrad_stream%up(lay,ist,i)
          endif

          auxrad_stream%meanrad_down(ist,i) = auxrad_stream%meanrad_down(ist,i) + auxrad_stream%down(lay,ist,i)
          auxrad_stream%meanrad_down(ist,i) = max(auxrad_stream%meanrad_down(ist,i), 0._jprb)
      enddo
    endif
  enddo

end subroutine calc_near_surf_contribution

subroutine solar_scattering_air(transmission_aux, solar, solar_spectrum, refl, raytracing, chanprof, &
                         transmission_scatt_ir, auxrad_stream)
  Implicit None

  type(transmission_type_aux), intent(in)       :: transmission_aux
  type(raytracing_type), intent(in)             :: raytracing
  type(radiance_aux), intent(inout)             :: auxrad_stream
  type(rttov_chanprof), intent(in)              :: chanprof(:)
  type(transmission_scatt_ir_type), intent(in)  :: transmission_scatt_ir

  logical(jplm), intent(in)                     :: solar(:)
  real(jprb), intent(in)                        :: solar_spectrum(:)
  real(jprb), intent(in)                        :: refl(:)
  real(jprb)                                    :: temp(nlayers,0:totnstreams)

! The solar_path1 and solar_path2 quantities are completely consistent (i.e.
! they are based on the same optical depth regression, namely the solar one).
! It is important to note that solar_path2%tau_level is on the sun-surface-satellite path
! but solar_path2%od_single_layer is on the sun-surface path.

#define fac1_2 auxrad_stream%Fac1_2(lay,ist,i)
#define fac2_2 auxrad_stream%Fac2_2(lay,ist,i)
#define fac3_2 auxrad_stream%Fac3_2(lay,ist,i)
#define fac4_2 auxrad_stream%Fac4_2(lay,ist,i)
#define fac5_2 auxrad_stream%Fac5_2(lay,ist,i)
#define fac6_2 auxrad_stream%Fac6_2(lay,ist,i)
#define fac7_2 auxrad_stream%Fac7_2(lay,ist,i)
#define dfac54_2 (fac5_2 - fac4_2)
#define tausun_layer transmission_aux%solar_path2%Tau_level(lay,ist,i)

  do i = 1, nchannels
    if (solar(i)) then

      do ist = 0, nstreams(i)
        do lay = 1, nlayers
          isti = ircld%icldarr(ist,lay,prof)
          auxrad_stream%Fac6_2(lay,ist,i) = solar_spectrum(i) * z4pi_r * &
                                          transmission_scatt_ir%azphacdo(isti,lay,i)
          auxrad_stream%Fac1_2(lay,ist,i) = solar_spectrum(i) * z4pi_r * &
                                          transmission_scatt_ir%azphacup(isti,lay,i)
          auxrad_stream%Fac2_2(lay,ist,i) = transmission_scatt_ir%ssa(isti,lay,i)
        enddo
      enddo

      do lay = 1, nlayers
        auxrad_stream%Fac3_2(lay,:,i) = raytracing%pathsat(lay,prof) / raytracing%patheff(lay,prof)
      enddo

      do ist = 0, nstreams(i)
        do lay = 1, nlayers
          isti = ircld%icldarr(ist,lay,prof)
          auxrad_stream%Fac4_2(lay,ist,i) = exp(-transmission_aux%solar_path2%Od_singlelayer(isti,lay,i))
          auxrad_stream%Fac5_2(lay,ist,i) = exp(-transmission_aux%solar_path1%Od_singlelayer(isti,lay,i))
        enddo
      enddo

      do lay = 1, nlayers
         if(.not. sateqsun(lay,prof)) then
            auxrad_stream%Fac7_2(lay,:,i) = raytracing%pathsat(lay,prof) / &
                                           (raytracing%pathsun(lay,prof) - raytracing%pathsat(lay,prof))
         endif
      enddo

  !----------------Upward single scattering of the solar beam-----------------------
  !
  ! auxrad_stream has already been layer integrated so must add integrated (not single layer) flux to each layer
        do ist=0, nstreams(i)
          temp(:,ist) = (auxrad_stream%Fac1_2(:,ist,i) * &
                        auxrad_stream%Fac2_2(:,ist,i) * &
                        auxrad_stream%Fac3_2(:,ist,i) * &
                        (transmission_aux%solar_path2%Tau_level(1:nlevels-1,ist,i) - &
                        transmission_aux%solar_path2%Tau_level(2:nlevels,ist,i)))
        enddo

        do ist=0, nstreams(i)
           do lay=2, nlayers
              temp(lay,ist) = temp(lay,ist) + temp(lay-1,ist) 
           enddo
        enddo

        do ist=0, nstreams(i)
          auxrad_stream%up_solar(:,ist,i) = temp(:,ist)
        enddo

  !-------------------Downward single scattering of the solar beam------------------
      ! Downward beam isn't necessary if the surface is dark
      if (refl(i) > 0._jprb) then

        do ist = 0, nstreams(i)
          do lay = 1, nlayers

            temp(lay,ist) = fac2_solar_path1 * &
                            fac6_2 * &
                            fac2_2 * &
                            tausun_layer * &
                            transmission_aux%solar_path1%Tau_level_r(lay,ist,i) * &
                            transmission_aux%solar_path1%Tau_level_r(lay+1,ist,i)

            if (.not. sateqsun(lay,prof)) then
              temp(lay,ist) = temp(lay,ist) * &
                              fac7_2 * &
                              dfac54_2
            else
              isti = ircld%icldarr(ist,lay,prof)
              temp(lay,ist) = temp(lay,ist) * &
                              fac4_2 * &
                              transmission_aux%solar_path2%Od_singlelayer(isti,lay,i)
            endif
          enddo
        enddo

        do ist = 0, nstreams(i)
          lay = 1_jpim
          auxrad_stream%down_ref_solar(lay,ist,i) = temp(lay,ist)
          auxrad_stream%down_solar(lay,ist,i) = max(auxrad_stream%down_ref_solar(lay,ist,i), 0.0_jprb)
          do lay = 2, nlayers
               temp(lay,ist) = temp(lay-1,ist) + temp(lay,ist)
               auxrad_stream%down_ref_solar(lay,ist,i) = temp(lay,ist)
               auxrad_stream%down_solar(lay,ist,i) = max(auxrad_stream%down_ref_solar(lay,ist,i), 0.0_jprb)
          enddo
        enddo

      endif ! refl(i) > 0.

    endif !solar(i)
  enddo

end subroutine solar_scattering_air

subroutine solar_scattering_near_surf(transmission_aux, solar, refl, chanprof, &
                                          auxrad_stream)
  Implicit None
  
  type(transmission_type_aux), intent(in)    :: transmission_aux
  logical(jplm), intent(in)                  :: solar(:)
  real(jprb), intent(in)                     :: refl(:)
  type(rttov_chanprof), intent(in)           :: chanprof(:)

  type(radiance_aux), intent(inout)          :: auxrad_stream

  integer(jpim) :: lay1
  real(jprb) :: temp(0:totnstreams)

#define fac1_3 auxrad_stream%Fac1_3(ist,i)
#define fac2_3 auxrad_stream%Fac2_3(ist,i)
#define fac3_3 auxrad_stream%Fac3_3(ist,i)
#define fac4_3 auxrad_stream%Fac4_3(ist,i)
#define fac5_3 auxrad_stream%Fac5_3(ist,i)
#define fac6_3 auxrad_stream%Fac6_3(ist,i)
#define fac7_3 auxrad_stream%Fac7_3(ist,i)
#define dfac54_3 (fac5_3 - fac4_3)
#define dtausun_surf (transmission_aux%solar_path2%Tau_level(lev,ist,i) - transmission_aux%solar_path2%Tau_surf(ist,i))
#define tausun_level transmission_aux%solar_path2%Tau_level(lev,ist,i) 

  do i=1, nchannels
     if(solar(i)) then

        ! lay is the layer above the one containing the surface
        ! lev is the nearest layer above the surface

        lay = iv3lay(i)
        lev = iv3lev(i)

        ! lay1 is the layer containing the surface or the bottom layer
        !   if the surface lies below the bottom of the profile

        if (pfraction(i) < 0._jprb) then
          lay1 = lay
        else
          lay1 = lay + 1
        endif

        auxrad_stream%Fac4_3(0:nstreams(i),i) = exp(-transmission_aux%solar_path2%od_sfrac(0:nstreams(i),i))
        auxrad_stream%Fac5_3(0:nstreams(i),i) = exp(-transmission_aux%solar_path1%od_sfrac(0:nstreams(i),i))

        do ist = 0, nstreams(i)

           fac1_3 = fac1_2
           fac2_3 = fac2_2
           fac3_3 = fac3_2

          !--------------Upward single scattering of the solar beam-------------------------
           auxrad_stream%meanrad_up_solar(ist,i) = &
             fac1_3 * &
             fac2_3 * &
             fac3_3 * &
             dtausun_surf

        enddo

        do ist = 0, nstreams(i)
          auxrad_stream%meanrad_up_solar(ist,i) = auxrad_stream%meanrad_up_solar(ist,i) + &
                                                  auxrad_stream%up_solar(lay,ist,i)
        enddo

          !--------------Downward single scattering of the solar beam-----------------------

        if (refl(i) > 0._jprb) then

          do ist = 0, nstreams(i)
             fac6_3 = fac6_2
             temp(ist) = fac2_solar_path1 * &
                         fac6_3 * &
                         fac2_3 * &
                         tausun_level * &
                         (transmission_aux%solar_path1%Tau_level_r(lay+1,ist,i) * &
                         transmission_aux%solar_path1%Tau_surf_r(ist,i))
          enddo

          if(.not. sateqsun(lay1,prof))then
             do ist = 0, nstreams(i)
                auxrad_stream%Fac7_3(ist,i) = auxrad_stream%Fac7_2(lay1,ist,i)

                auxrad_stream%meanrad_down_solar(ist,i) = &
                                                    temp(ist) * &
                                                    fac7_3 * &
                                                    dfac54_3
             enddo
          else
             do ist = 0, nstreams(i)
                auxrad_stream%meanrad_down_solar(ist,i) = &
                                    temp(ist) * fac4_3 * &
                                    transmission_aux%solar_path2%od_sfrac(ist,i)
             enddo
          endif


          do ist = 0, nstreams(i)
            auxrad_stream%meanrad_down_solar(ist,i) = auxrad_stream%meanrad_down_solar(ist,i) + &
                                                      auxrad_stream%down_solar(lay,ist,i)

            auxrad_stream%meanrad_down_solar(ist,i) = max(auxrad_stream%meanrad_down_solar(ist,i), 0._jprb)
          enddo

        endif ! refl(i) > 0.

     endif
  end do

end subroutine solar_scattering_near_surf

subroutine solar_rayleigh(coef, opts, solar, solar_spectrum, chanprof, raytracing, profiles, &
                          profiles_dry, transmission_aux, auxrad_stream)

  implicit none

  type(rttov_coef), intent(in)               :: coef
  type(rttov_options), intent(in)            :: opts
  logical(jplm), intent(in)                  :: solar(:)
  real(jprb), intent(in)                     :: solar_spectrum(:)
  type(rttov_chanprof), intent(in)           :: chanprof(:)
  type(raytracing_type), intent(in)          :: raytracing
  type(profile_type), intent(in)             :: profiles(:)
  type(profile_type), intent(in)             :: profiles_dry(:)
  type(transmission_type_aux), intent(in)    :: transmission_aux
  type(radiance_aux), intent(inout)          :: auxrad_stream

  real(kind=jprb) :: wlm, ss_param, v_h2o, Mwet, cosscata_term1, cosscata_term2, cosscata
  real(kind=jprb) :: ray_phase, solar_src, solar_src_updn
  real(kind=jprb) :: rayrad_up(0:nlayers,0:maxnstreams), rayrad_dn(0:nlayers,0:maxnstreams)

  ! The scattering cross-section calculation is based on Bucholzt '95.
  ! The phase function could be made to account for depolarisation.
  ! The solar geometry is approximated using the sun-surface path values,
  ! and similarly for the sun-level-satellite transmittances.

  ! Rayleigh scattering currently only included for channels less than 2um (ray_min_wvn).

  do i = 1, nchannels
    if (solar(i) .and. coef%ff_cwn(chan) > ray_min_wvn) then

      ! Calculate layer-independent scattering parameter
      wlm = 10000.0_jprb / coef%ff_cwn(chan)    ! Wavelength in microns

      if (wlm < ray_scs_wlm) then
        ss_param = ray_scs_a1 * wlm ** (ray_scs_b1 + ray_scs_c1 * wlm + ray_scs_d1/wlm)
      else
        ss_param = ray_scs_a2 * wlm ** (ray_scs_b2 + ray_scs_c2 * wlm + ray_scs_d2/wlm)
      endif
      ss_param = ss_param * 0.01_jprb ** 2_jpim * na * z4pi_r / gravity

      rayrad_up(0,:) = 0._jprb
      rayrad_dn(0,:) = 0._jprb

      ! Sum contributions from atmospheric layers
      do lev = 2, nlevels
        lay = lev - 1

        ! Layer H2O by volume as fraction:
        v_h2o = 0.5_jprb * (profiles_dry(prof)%q(lev-1) + profiles_dry(prof)%q(lev)) * 1.E-06_jprb
        ! Do nothing in compatibility mode: using ppmv dry where ppmv wet should be used.
        if (profiles(prof)%gas_units /= gas_unit_compatibility) then
          v_h2o = v_h2o / (1._jprb + v_h2o)
        endif

        ! Average molar weight of wet air for the layer (kg)
        Mwet = ((1._jprb - v_h2o) * Mair + v_h2o * Mh2o) * 1.E-3_jprb 

        ! cosine of scattering angle - raytracing%zasat/zasun contain the sine of the angles
        cosscata_term1 = sqrt((1._jprb - raytracing%zasat(lay, prof) * raytracing%zasat(lay, prof)) * &
                         (1._jprb - raytracing%zasun(lay, prof) * raytracing%zasun(lay, prof)))
        cosscata_term2 = raytracing%zasat(lay, prof) * raytracing%zasun(lay, prof) * &
                         cos((profiles(prof)%azangle - profiles(prof)%sunazangle)*deg2rad)

        solar_src = solar_spectrum(i) * & ! mW m^-2 (cm^-1)^-1
              (profiles(prof)%p(lev) - profiles(prof)%p(lev-1)) * 100.0_jprb * & ! convert hPa to Pa
              raytracing%pathsat(lay, prof) * ss_param / Mwet

        cosscata = - cosscata_term1 - cosscata_term2
        ray_phase = 0.75_jprb * (1.0_jprb + cosscata * cosscata)
        solar_src_updn = solar_src * ray_phase

        ! Phase function symmetry means upwelling and downwelling radiances are the same
        do ist = 0, nstreams(i)
          rayrad_up(lay,ist) = rayrad_up(lay-1,ist) + solar_src_updn * &
                               transmission_aux%solar_path2%Tau_level(lev-1,ist,i)

          rayrad_dn(lay,ist) = rayrad_dn(lay-1,ist) + solar_src_updn * &
                               transmission_aux%solar_path2%Tau_level(lev-1,ist,i) / &
                               transmission_aux%solar_path1%Tau_level(lev-1,ist,i) ** 3_jpim
        enddo
      enddo

      ! Add Rayleigh contributions to radiance totals
      do ist = 0, nstreams(i)
        auxrad_stream%up_solar(:,ist,i) = auxrad_stream%up_solar(:,ist,i) + &
                                          rayrad_up(1:nlayers,ist)

        auxrad_stream%meanrad_up_solar(ist,i) = auxrad_stream%meanrad_up_solar(ist,i) + &
                                                rayrad_up(iv3lay(i),ist)

        auxrad_stream%down_solar(:,ist,i) = auxrad_stream%down_solar(:,ist,i) + &
                                            rayrad_dn(1:nlayers,ist)

        auxrad_stream%meanrad_down_solar(ist,i) = auxrad_stream%meanrad_down_solar(ist,i) + &
                                                  rayrad_dn(iv3lay(i),ist)
      enddo

      ! Calculate the contribution from the part-layer above the surface

      ! Layer H2O by volume as fraction:
      if (opts%rt_all%use_q2m) then
        v_h2o = 0.5_jprb * (profiles_dry(prof)%q(iv3lev(i)) + profiles_dry(prof)%s2m%q) * 1.E-06_jprb
      else
        v_h2o = profiles_dry(prof)%q(iv3lev(i)) * 1.E-06_jprb
      endif
      ! Do nothing in compatibility mode: using ppmv dry where ppmv wet should be used.
      if (profiles(prof)%gas_units /= gas_unit_compatibility) then
        v_h2o = v_h2o / (1._jprb + v_h2o)
      endif

      ! Average molar weight of wet air for the layer (kg):
      Mwet = ((1._jprb - v_h2o) * Mair + v_h2o * Mh2o) * 1.E-3_jprb

      ! cosine of scattering angle
      cosscata_term1 = sqrt((1._jprb - raytracing%zasat(iv2lay(i), prof) * raytracing%zasat(iv2lay(i), prof)) * &
                        (1._jprb - raytracing%zasun(iv2lay(i), prof) * raytracing%zasun(iv2lay(i), prof)))
      cosscata_term2 = raytracing%zasat(iv2lay(i), prof) * raytracing%zasun(iv2lay(i), prof) * &
                        cos((profiles(prof)%azangle - profiles(prof)%sunazangle)*deg2rad)
      solar_src = solar_spectrum(i) * & ! mW m^-2 (cm^-1)^-1
            abs(profiles(prof)%s2m%p - profiles(prof)%p(iv3lev(i))) * 100.0_jprb * &  ! convert hPa to Pa
            raytracing%pathsat(iv2lay(i), prof) * ss_param / Mwet

      cosscata = - cosscata_term1 - cosscata_term2
      ray_phase = 0.75_jprb * (1.0_jprb + cosscata * cosscata)
      solar_src_updn = solar_src * ray_phase

      ! Add near-surface contributions to the total radiances
      do ist = 0, nstreams(i)
        auxrad_stream%meanrad_up_solar(ist,i) = auxrad_stream%meanrad_up_solar(ist,i) + solar_src_updn * &
                                                transmission_aux%solar_path2%Tau_level(iv3lev(i),ist,i)

        auxrad_stream%meanrad_down_solar(ist,i) = auxrad_stream%meanrad_down_solar(ist,i) + solar_src_updn * &
                                                  transmission_aux%solar_path2%Tau_level(iv3lev(i),ist,i) / &
                                                  transmission_aux%solar_path1%Tau_level(iv3lev(i),ist,i) ** 3_jpim
      enddo

    endif

  enddo

end subroutine solar_rayleigh

subroutine solar_surface_contribution(transmission_aux, solar, solar_spectrum, reflectance, &
                                          refl_norm, auxrad_stream, rad)
  Implicit None

  type(transmission_type_aux), intent(in)    :: transmission_aux
  logical(jplm), intent(in)                  :: solar(:)
  real(jprb), intent(in)                     :: solar_spectrum(:)
  real(jprb), intent(in)                     :: reflectance(:) ! surface reflectance
  real(jprb), intent(in)                     :: refl_norm(:)
  type(radiance_aux), intent(inout)          :: auxrad_stream
  type(radiance_type), intent(inout)         :: rad

  real(jprb)    :: temp(0:totnstreams) ! ist

  Do i = 1, nchannels
     if(solar(i)) then

        temp(0:nstreams(i)) = solar_spectrum(i) * reflectance(i) * refl_norm(i) * &
                              transmission_aux%solar_path2%Tau_surf(0:nstreams(i),i)

        rad%clear(i) = rad%clear(i) + temp(0)
        auxrad_stream%cloudy(1:nstreams(i),i)= auxrad_stream%cloudy(1:nstreams(i),i) + temp(1:nstreams(i))

     endif
  enddo
end subroutine solar_surface_contribution

subroutine init_rttov_integrate_mod_vars(chanprof, profiles, ircld, aux_prof, coef, raytracing)

  type(rttov_chanprof), intent(in)              :: chanprof(:)     ! Channel indices
  type(profile_type),   intent(in)              :: profiles(:) ! Profiles
  type(ircld_type),     intent(in)              :: ircld
  type(raytracing_type),intent(in)              :: raytracing
  type(profile_aux) ,   intent(in)              :: aux_prof ! auxillary profiles info.
  type(rttov_coef),     intent(in)              :: coef

!store calculation size
nchannels = size(chanprof)
nlayers = profiles(1)%nlayers
nlevels = nlayers + 1

! allocate arrays
! DAR: Will this be a problem for malloc locking for parallel users?
allocate(nstreams(nchannels), iv2lev(nchannels), iv2lay(nchannels), iv3lev(nchannels), iv3lay(nchannels), pol_id(nchannels), &
         cfraction(nchannels), pfraction(nchannels), sateqsun(nlayers,size(profiles)))

! populate module arrays
Do i = 1, nchannels
   nstreams(i) = ircld%nstream(prof)
enddo

do i = 1, nchannels
   cfraction(i) = aux_prof%s(prof)%cfraction
   pfraction(i) = aux_prof%s(prof)%pfraction_surf
enddo

Do i = 1, nchannels

! case-1: surf lies above lev=nlevels
   iv3lev(i) = aux_prof%s(prof)%nearestlev_surf - 1   ! lowest lev above surf
! case-2: surf lies below lev=nlevels
   if (pfraction(i) < 0.0_JPRB) iv3lev(i) = iv3lev(i) + 1  ! iv3lev=iv2lev=lowest lev above surf

   iv2lev(i) = aux_prof%s(prof)%nearestlev_surf       ! highest lev below surf
   iv2lay(i) = iv2lev(i) - 1                          ! same layer as that numbered by  iv2 in RTTOV-9
   iv3lay(i) = iv3lev(i) - 1                          ! same layer as that numbered by  iv3 in RTTOV-9
End Do

if(coef%id_sensor == sensor_id_po) then
   do i = 1, nchannels
      pol_id(i) = coef%fastem_polar(chan) + 1_jpim
   enddo
Else
   pol_id(:) = 0_jpim
Endif 

if (dosolar) then
  do iprof = 1, size(profiles)
     do lay = 1, nlayers
        sateqsun(lay,iprof) = .false.
        if(raytracing%pathsat(lay,iprof) == raytracing%pathsun(lay,iprof)) sateqsun(lay,iprof) = .true.
     enddo
  enddo
endif

end subroutine init_rttov_integrate_mod_vars

end subroutine rttov_integrate

