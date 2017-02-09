Subroutine rttov_integrate_ad( &
   addcosmic, opts, maxnstreams, chanprof,                       &
   emissivity,                   emissivity_ad,                  &
   reflectance,                  reflectance_ad,                 &
   refl_norm,                                                    &
   thermrefl,                    thermrefl_ad,                   &
   do_lambertian,                                                &
   thermal,                                                      &
   dothermal,                                                    &
   solar,                                                        &
   dosolar,                                                      &
   solar_spectrum,                                               &
   transmission_aux,             transmission_aux_ad,            &
   transmission_scatt_ir,        transmission_scatt_ir_ad,       &
   profiles,                     profiles_ad,                    &
   profiles_dry,                 profiles_dry_ad,                &
   aux_prof,                     aux_prof_ad,                    &
   coef,                                                         &
   raytracing,                   raytracing_ad,                  &
   ircld,                        ircld_ad,                       &
   rad,                                                          &
   auxrad,                                                       &
                                 auxrad_stream,                  &
                                 auxrad_stream_ad,               &
                                 rad_ad)

  ! in,                           !inout
!
! Description:
! To perform AD of integration of radiative transfer equation
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
! ECMWF Research Dept. Tech. Memo. 425 (available from the librarian at ECMWF)  !
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
!  1.4     10/05/2004   Fixed bug in call to rttov_calcpolarisation (D.Salmond)
!  1.5     06/09/2004   Mods. for Vectorisation (D Salmond ECMWF & B  Carruthers, Cray)
!  1.6     28/02/2005   More improvements to vectorisation (D Dent)
!  1.7     29/03/2005   Add end of header comment (J. Cameron)
!  1.8     03/03/2006   Marco Matricardi (ECMWF):
!             --        IASI capability added.
!             --        Linear in tau approximation for RT equation introduced.
!             --        Solar radiation introduced for IASI and AIRS.
!  1.9     02/01/2007   Corrected some bugs in the old code (R Saunders)
!  1.10    08/02/2007   Removed polarisation indeces (R Saunders)
!  1.11    11/03/2007   Reintroduced overcast radiance (R Saunders)
!  1.12    27/11/2007   Optimised for IBM/NEC (D Salmond)
!  1.13    17/11/2007   Allowed adjoint/K sensitivity to clear BT/radiance (A Geer)
!  1.14    21/12/2007   Added polarimetric option (R. Saunders)
!  1.15    15/08/2009   User defined ToA. Layers distinct from levels. Top-layer
!                       brought into layer looping to shorten code (P.Rayer)
!  1.16    03/11/2009   Transmittances on levels (A Geer)
!  1.17    02/12/2009   Introduced principal component capability. Pathsat, Pathsun and
!                       related quantities are now layer arrays (Marco Matricardi).
!  1.18    05/07/2010   Remove addsolar flag from profiles structure (J Hocking)
!  1.19    14/10/2010   Remove rt8_mode (J Hocking)
!  1.20    14/12/2010   Use traj0_sta%solar array to flag channels for which solar calculations
!                       should be performed (J Hocking)
!  2.0     14/12/2011   Re-written (D Rundle)
!
! 2010/03 Code cleaning Pascal Brunel, Philippe Marguinaud
!
!
! Code Description:
!   Language:           Fortran 90.
!   Software Standards: "European Standards for Writing and
!     Documenting Exchangeable Fortran 90 Code".

  USE rttov_types, ONLY : &
    rttov_chanprof, rttov_coef, profile_Type, profile_aux, transmission_type_aux, transmission_scatt_ir_type, &
    radiance_Type, rttov_options, ircld_type, raytracing_type, radiance_aux, rttov_emissivity, &
    rttov_reflectance

  USE parkind1, ONLY : jpim, jprb, jplm

!INTF_OFF
  USE rttov_const, ONLY : sensor_id_po, pi_r, deg2rad, &
    overcast_albedo_wvn, overcast_albedo1, overcast_albedo2, &
    gas_unit_compatibility

  Use yomhook, Only : LHOOK, DR_HOOK
!INTF_ON

  IMPLICIT NONE

!subroutine arguments:
  LOGICAL(jplm), INTENT(in)                        :: addcosmic
  TYPE(rttov_options), INTENT(in)                  :: opts
  INTEGER(jpim), INTENT(in)                        :: maxnstreams
  TYPE(rttov_chanprof), INTENT(in)                 :: chanprof(:)
  TYPE(profile_Type), INTENT(in)                   :: profiles(:)
  TYPE(profile_Type), INTENT(in)                   :: profiles_dry(SIZE(profiles))
  TYPE(rttov_emissivity),  INTENT(in), OPTIONAL    :: emissivity(SIZE(chanprof))
  TYPE(rttov_reflectance), INTENT(in), OPTIONAL    :: reflectance(SIZE(chanprof))
  REAL(jprb), INTENT(in)                           :: refl_norm(SIZE(chanprof))
  REAL(jprb), INTENT(in)                           :: thermrefl(SIZE(chanprof))
  LOGICAL(jplm), INTENT(in)                        :: do_lambertian(SIZE(chanprof))
  LOGICAL(jplm), INTENT(in)                        :: thermal(SIZE(chanprof))
  LOGICAL(jplm), INTENT(in)                        :: dothermal
  LOGICAL(jplm), INTENT(in)                        :: solar(SIZE(chanprof))
  LOGICAL(jplm), INTENT(in)                        :: dosolar
  REAL(jprb), INTENT(in)                           :: solar_spectrum(SIZE(chanprof))
  TYPE(rttov_coef), INTENT(in)                     :: coef
  TYPE(profile_aux) , INTENT(in)                   :: aux_prof
  TYPE(transmission_Type_aux), INTENT(in)          :: transmission_aux
  TYPE(transmission_scatt_ir_type), INTENT(in)     :: transmission_scatt_ir
  TYPE(ircld_type), INTENT(in)                     :: ircld
  TYPE(raytracing_type), INTENT(in)                :: raytracing
  TYPE(radiance_Type), INTENT(in)                  :: rad
  TYPE(radiance_aux), INTENT(in)                   :: auxrad
  TYPE(radiance_aux), INTENT(in)                   :: auxrad_stream
  
  TYPE(radiance_aux), INTENT(inout)                :: auxrad_stream_ad
  TYPE(rttov_emissivity),  INTENT(inout), OPTIONAL :: emissivity_ad(SIZE(chanprof))
  TYPE(rttov_reflectance), INTENT(inout), OPTIONAL :: reflectance_ad(SIZE(chanprof))
  REAL(jprb), INTENT(inout)                        :: thermrefl_ad(SIZE(chanprof))
  TYPE(profile_Type), INTENT(inout)                :: profiles_ad(SIZE(profiles))
  TYPE(profile_Type), INTENT(inout)                :: profiles_dry_ad(SIZE(profiles))
  TYPE(profile_aux), INTENT(inout)                 :: aux_prof_ad
  TYPE(transmission_Type_aux), INTENT(inout)       :: transmission_aux_ad
  TYPE(transmission_scatt_ir_type), INTENT(inout)  :: transmission_scatt_ir_ad
  TYPE(ircld_type), INTENT(inout)                  :: ircld_ad
  TYPE(raytracing_type), INTENT(inout)             :: raytracing_ad
  TYPE(radiance_Type), INTENT(inout)               :: rad_ad
  
!INTF_END
  
#include "rttov_calcrad_ad.interface"

!local variables:
  REAL(jprb)    :: cfraction(SIZE(chanprof))
  REAL(jprb)    :: cfraction_ad(SIZE(chanprof))
  
  REAL(jprb)    :: rad_air_ad(profiles(1)%nlevels,SIZE(chanprof))
  REAL(jprb)    :: rad_surfair_ad(SIZE(chanprof))
  REAL(jprb)    :: rad_skin_ad(SIZE(chanprof))
  
  REAL(jprb)    :: pfraction(SIZE(chanprof))        ! cloud fraction
  LOGICAL(jplm) :: sateqsun(profiles(1)%nlayers,SIZE(profiles(:))) ! True where the solar zenith angle equal to observation angle
  INTEGER(jpim) :: pol_id(SIZE(chanprof))       ! polarisation index
  
  INTEGER(jpim) :: iv2lay(SIZE(chanprof)), iv2lev(SIZE(chanprof)), iv3lay(SIZE(chanprof)), iv3lev(SIZE(chanprof))
  
  INTEGER(jpim) :: i, lev, lay, iprof, ist, isti, nlayers, nlevels, nchannels
  
  REAL(jprb)    :: p, refl, refl_norm_scat
  
  LOGICAL(jplm) :: keyradonly ! flag to indicate only calculate primary radiance outputs
  
  INTEGER(jpim) :: narray(3)
  
  REAL(JPRB) :: ZHOOK_HANDLE
  
  IF (LHOOK) CALL DR_HOOK('RTTOV_INTEGRATE_AD',0_jpim,ZHOOK_HANDLE)
  !- End of header --------------------------------------------------------
  
#define prof chanprof(i)%Prof
#define chan chanprof(i)%Chan

!X.  initialisation of local variables

  cfraction_ad(:)   = 0._JPRB
  rad_surfair_ad(:) = 0._JPRB
  rad_skin_ad(:)    = 0._JPRB
  rad_air_ad(:,:)   = 0._jprb
   
  ! auxrad_stream_ad is initialised in rttov_ad

!---------------------------
!0. Initialise useful variables
!---------------------------

  nchannels = SIZE(chanprof)
  nlayers = profiles(1)%nlayers
  nlevels = nlayers + 1

  DO i = 1, nchannels
    IF(do_lambertian(i)) THEN
      transmission_aux_ad%thermal_path1%tau_level_p(:,:,i)   = 0._jprb
      transmission_aux_ad%thermal_path1%tau_level_p_r(:,:,i) = 0._jprb
      transmission_aux_ad%thermal_path1%tau_surf_p_r(:,i)    = 0._jprb
      transmission_aux_ad%thermal_path1%tau_surf_p(:,i)      = 0._jprb
    ENDIF
  ENDDO
  
! DAR: narray is a more compact way of passing nchannels, nlevels and maxnstreams to the internal subroutines
! Again module variables would be a more elegant way of solving this problem and would results in less duplicated code
narray(1) = nlevels; narray(2) = maxnstreams; narray(3) = nchannels

keyradonly = opts%rt_ir%addaerosl .or. opts%rt_ir%pc%addpc

do i = 1, nchannels
  cfraction(i) = aux_prof%s(prof)%cfraction
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
    pol_id(i) = coef%fastem_polar(chan) + 1_jpim
  enddo
else
  pol_id(:) = 0_jpim
endif


!---------------------------
! 7. Calculate total radiance
!---------------------------
if (opts%rt_ir%addaerosl .OR. opts%rt_ir%addclouds) then
   rad_ad%cloudy(:) = rad_ad%total(:)
!----------------------------------------
! Calculate complex cloudy radiances
!----------------------------------------
   do i = 1, nchannels
      rad_ad%clear(i)        = rad_ad%cloudy(i) * ircld%XSTRCLR(prof)
      ircld_ad%XSTRCLR(prof) = ircld_ad%XSTRCLR(prof) + rad_ad%cloudy(i) * rad%clear(i)

      do ist = ircld%nstream(prof), 1, -1 !reverse loop nec. due to data dependence
         auxrad_stream_ad%cloudy(ist,i) = auxrad_stream_ad%cloudy(ist,i) + &
                                          rad_ad%cloudy(i) * (ircld%xstr(ist+1,prof) - ircld%xstr(ist,prof))

         ircld_ad%xstr(ist+1,prof)      = ircld_ad%xstr(ist+1,prof) + &
                                          rad_ad%cloudy(i) * auxrad_stream%cloudy(ist,i)

         ircld_ad%xstr(ist,prof)        = ircld_ad%xstr(ist,prof) - &
                                          rad_ad%cloudy(i) * auxrad_stream%cloudy(ist,i)
      enddo
   enddo
else
!---------------------------
! Calculate total radiance (clear case/simple cloud)
!---------------------------
  if (opts%rt_ir%pc%addpc) then
    rad_ad%clear(:)  = rad_ad%total(:)
  else
    rad_ad%clear(1:nchannels)  = rad_ad%clear(1:nchannels) + (1.0_jprb - cfraction(:)) * rad_ad%total(1:nchannels)
    rad_ad%cloudy(1:nchannels) = cfraction(:) * rad_ad%total(1:nchannels)
    cfraction_ad(:)  = (rad%cloudy(1:nchannels) - rad%clear(1:nchannels)) * rad_ad%total(1:nchannels)

    ! Interpolate to given cloud-top pressures
    do i = 1, nchannels
      lay = aux_prof%s(prof)%nearestlev_ctp - 1
      p   = aux_prof%s(prof)%pfraction_ctp

      rad_ad%overcast(lay,i)            = & !rad_ad%overcast(lay,i) + & !first use of rad_ad%overcast
                                         (1._jprb - p) * rad_ad%cloudy(i)

      rad_ad%overcast(lay-1,i)          = & !rad_ad%overcast(lay-1,i) + & !first use of rad_ad%overcast
                                          p * rad_ad%cloudy(i)

      aux_prof_ad%s(prof)%pfraction_ctp = aux_prof_ad%s(prof)%pfraction_ctp + &
                                         (rad%overcast(lay-1,i) - rad%overcast(lay,i)) * rad_ad%cloudy(i)
    enddo

  endif
! rad_ad%cloudy done
endif


!---------------------------------------------------
!6. Calculate overcast radiances
!---------------------------------------------------
if (.not. keyradonly) then
  ist = 0_jpim

!cdir nodep
  do i = 1, nchannels
    if (thermal(i)) then
      do lay = 1, nlayers
        ! overcast radiances at given cloud top

        ! Make exception for iv2lay because up(iv2lay) was replaced by meanrad_up in direct and
        ! overcast(iv2lay) was replaced by surface contribution
        if (lay == iv2lay(i)) then
          auxrad_stream_ad%meanrad_up(ist,i) = auxrad_stream_ad%meanrad_up(ist,i) + rad_ad%overcast(lay,i)
          rad_surfair_ad(i) = rad_surfair_ad(i) + rad_ad%overcast(lay,i) * &
                              transmission_aux%thermal_path1%tau_surf(ist,i)
          transmission_aux_ad%thermal_path1%tau_surf(ist,i) = &
                              transmission_aux_ad%thermal_path1%tau_surf(ist,i) + &
                              rad_ad%overcast(lay,i) * auxrad%surfair(i)
        else
          auxrad_stream_ad%up(lay,ist,i) = auxrad_stream_ad%up(lay,ist,i) + rad_ad%overcast(lay,i)
          rad_air_ad(lay+1,i) = &!rad_air_ad(lay+1,i) + &
                              rad_ad%overcast(lay,i) * transmission_aux%thermal_path1%tau_level(lay+1,ist,i)
          transmission_aux_ad%thermal_path1%tau_level(lay+1,ist,i) = &
                              transmission_aux_ad%thermal_path1%tau_level(lay+1,ist,i) + &
                              auxrad%air(lay+1,i) * rad_ad%overcast(lay,i)
        endif
      enddo
    else if (solar(i)) then
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
        transmission_aux_ad%solar_path2%Tau_level(lev,ist,i) = &
                transmission_aux_ad%solar_path2%Tau_level(lev,ist,i) + &
                solar_spectrum(i) * refl / raytracing%pathsun(lay, prof) * &
                rad_ad%overcast(lay,i)
        raytracing_ad%pathsun(lay, prof) = raytracing_ad%pathsun(lay, prof) - &
                solar_spectrum(i) * refl / raytracing%pathsun(lay, prof)**2_jpim * &
                transmission_aux%solar_path2%Tau_level(lev,ist,i) * rad_ad%overcast(lay,i)

        transmission_aux_ad%solar_path1%Tau_level(lev,ist,i) = &
                transmission_aux_ad%solar_path1%Tau_level(lev,ist,i) + &
                refl / raytracing%pathsat(lay, prof) * &
                auxrad_stream%down_solar(lay,ist,i) * &
                transmission_aux%solar_path1%Tau_level(lev,ist,i) * &
                2_jpim * rad_ad%overcast(lay,i)                
        raytracing_ad%pathsat(lay, prof) = &
                raytracing_ad%pathsat(lay, prof) - &
                refl / raytracing%pathsat(lay, prof) * &
                transmission_aux%solar_path1%Tau_level(lev,ist,i)**2_jpim * &
                auxrad_stream%down_solar(lay,ist,i) / raytracing%pathsat(lay, prof) * &
                rad_ad%overcast(lay,i)
        auxrad_stream_ad%down_solar(lay,ist,i) = &
                auxrad_stream_ad%down_solar(lay,ist,i) + &
                refl / raytracing%pathsat(lay, prof) * &
                transmission_aux%solar_path1%Tau_level(lev,ist,i)**2_jpim * &
                rad_ad%overcast(lay,i)
                
        ! Make exception for iv2lay because up_solar(iv2lay) was replaced by meanrad_up_solar in direct
        if (lay == iv2lay(i)) then
          auxrad_stream_ad%meanrad_up_solar(ist,i) = &
                 auxrad_stream_ad%meanrad_up_solar(ist,i) + rad_ad%overcast(lay,i)
        else
          auxrad_stream_ad%up_solar(lay,ist,i) = &
                 auxrad_stream_ad%up_solar(lay,ist,i) + rad_ad%overcast(lay,i)
        endif
      enddo
    else
      rad_ad%overcast(:,i) = 0._jprb
    endif
  enddo
  ! rad_ad%overcast done
  ! rad_ad%clear done

  ! aux_prof_ad%cfraction write out now
  do i = 1, nchannels
    aux_prof_ad%s(prof)%cfraction = aux_prof_ad%s(prof)%cfraction + cfraction_ad(i)
  enddo
endif

!--------------------------------
!5. cosmic temperature correction
!--------------------------------

!calculate planck function corresponding to tcosmic=2.7k
!deblonde tcosmic for microwave sensors only

if (addcosmic) then
  do i = 1, nchannels
    ist = 0_jpim

    if (do_lambertian(i)) then
      thermrefl_ad(i) = thermrefl_ad(i) + &
        auxrad%cosmic(i) * &
        transmission_aux%thermal_path1%tau_surf_p(ist,i) * &
        transmission_aux%thermal_path1%tau_surf(ist,i) * &
        rad_ad%clear(i)

      transmission_aux_ad%thermal_path1%tau_surf(ist,i) = &! transmission_aux_ad%thermal_path1%tau_surf(ist,i) + &
        auxrad%cosmic(i) * thermrefl(i) * &
        transmission_aux%thermal_path1%tau_surf_p(ist,i) * &
        rad_ad%clear(i)

      transmission_aux_ad%thermal_path1%tau_surf_p(ist,i) = &
        transmission_aux_ad%thermal_path1%tau_surf_p(ist,i) + &
        auxrad%cosmic(i) * thermrefl(i) * &
        transmission_aux%thermal_path1%tau_surf(ist,i) * &
        rad_ad%clear(i)
    ELSE
      thermrefl_ad(i) = thermrefl_ad(i) + &
                           rad_ad%clear(i) * auxrad%cosmic(i) * &
                           transmission_aux%thermal_path1%tau_surf(ist,i)**2_jpim

      transmission_aux_ad%thermal_path1%tau_surf(ist,i) = &! transmission_aux_ad%thermal_path1%tau_surf(ist,i) +
                           rad_ad%clear(i) * 2.0_jprb * thermrefl(i) * auxrad%cosmic(i) * &
                           transmission_aux%thermal_path1%tau_surf(ist,i)
    endif
  enddo
endif

!----------------------------------------------------------------------------------------
! 4. solar surface contribution
!----------------------------------------------------------------------------------------
if (dosolar) &
  call solar_surface_contribution_ad(transmission_aux, transmission_aux_ad, &
                                     ircld, reflectance%refl_out, reflectance_ad%refl_out, refl_norm, &
                                     chanprof, solar, solar_spectrum, narray, &
                                     auxrad_stream_ad, rad_ad)

!-----------------------
!3. surface emission contribution
!-----------------------
if (dothermal) then
  do i = 1, nchannels
    if (thermal(i)) then
      ist = 0_jpim

      emissivity_ad(i)%emis_out = emissivity_ad(i)%emis_out + &
                        rad_ad%clear(i) * auxrad%skin(i) * transmission_aux%thermal_path1%tau_surf(ist,i)

      transmission_aux_ad%thermal_path1%tau_surf(ist,i) = transmission_aux_ad%thermal_path1%tau_surf(ist,i) + &
                                            rad_ad%clear(i) * auxrad%skin(i) * emissivity(i)%emis_out

      rad_skin_ad(i) = rad_skin_ad(i) + &
                      rad_ad%clear(i) * emissivity(i)%emis_out * transmission_aux%thermal_path1%tau_surf(ist,i)

      do ist = 1, ircld%nstream(prof) !rev loop not necessary
        emissivity_ad(i)%emis_out = emissivity_ad(i)%emis_out + &
                            auxrad_stream_ad%cloudy(ist,i) * auxrad%skin(i) * &
                            transmission_aux%thermal_path1%tau_surf(ist,i)

        transmission_aux_ad%thermal_path1%tau_surf(ist,i) = transmission_aux_ad%thermal_path1%tau_surf(ist,i) + &
                                              auxrad_stream_ad%cloudy(ist,i) * auxrad%skin(i) * emissivity(i)%emis_out

        rad_skin_ad(i) = rad_skin_ad(i) + &
                          auxrad_stream_ad%cloudy(ist,i) * emissivity(i)%emis_out * &
                          transmission_aux%thermal_path1%tau_surf(ist,i)
      enddo
    endif
  enddo
endif

!-------------------------------------
!2. calculate atmospheric contribution
!-------------------------------------
!cdir nodep
do i = 1, nchannels
  ist = 0_jpim

  if (thermal(i)) then
    auxrad_stream_ad%meanrad_up(ist,i) = auxrad_stream_ad%meanrad_up(ist,i) + rad_ad%clear(i)

    if (do_lambertian(i)) then
      auxrad_stream_ad%meanrad_down(ist,i) = &!auxrad_stream_ad%meanrad_down(ist,i) + &
        thermrefl(i) * transmission_aux%thermal_path1%tau_surf_p(ist,i) * &
        transmission_aux%thermal_path1%tau_surf(ist,i) * rad_ad%clear(i)

      thermrefl_ad(i) = thermrefl_ad(i) + &
        auxrad_stream%meanrad_down(ist,i) * &
        transmission_aux%thermal_path1%tau_surf_p(ist,i) * &
        transmission_aux%thermal_path1%tau_surf(ist,i) * rad_ad%clear(i)

      transmission_aux_ad%thermal_path1%tau_surf_p(ist,i) = &
        transmission_aux_ad%thermal_path1%tau_surf_p(ist,i) + &
        auxrad_stream%meanrad_down(ist,i) * thermrefl(i) * &
        transmission_aux%thermal_path1%tau_surf(ist,i) * rad_ad%clear(i)

      transmission_aux_ad%thermal_path1%tau_surf(ist,i) = &
        transmission_aux_ad%thermal_path1%tau_surf(ist,i) + &
        auxrad_stream%meanrad_down(ist,i) * thermrefl(i) * &
        transmission_aux%thermal_path1%tau_surf_p(ist,i) * rad_ad%clear(i)
    else

      auxrad_stream_ad%meanrad_down(ist,i) = & ! auxrad_stream_ad%meanrad_down(ist,i) + &
                                rad_ad%clear(i) * thermrefl(i) * transmission_aux%thermal_path1%tau_surf(ist,i)**2

      transmission_aux_ad%thermal_path1%tau_surf(ist,i) = transmission_aux_ad%thermal_path1%tau_surf(ist,i) + &
                                            rad_ad%clear(i) * 2._JPRB * thermrefl(i) * &
                                            transmission_aux%thermal_path1%tau_surf(ist,i) * &
                                            auxrad_stream%meanrad_down(ist,i)

      thermrefl_ad(i) = thermrefl_ad(i) + &
                            rad_ad%clear(i) * transmission_aux%thermal_path1%tau_surf(ist,i)**2_jpim * &
                            auxrad_stream%meanrad_down(ist,i)
    endif
  endif

  if (solar(i)) then
    auxrad_stream_ad%meanrad_up_solar(ist,i) = auxrad_stream_ad%meanrad_up_solar(ist,i) + rad_ad%clear(i)

    refl_norm_scat = COS(profiles(prof)%zenangle * deg2rad)

    reflectance_ad(i)%refl_out = reflectance_ad(i)%refl_out + &
              transmission_aux%solar_path1%Tau_surf(ist,i)**2_jpim * refl_norm_scat * &
              auxrad_stream%meanrad_down_solar(ist,i) * rad_ad%clear(i)

    auxrad_stream_ad%meanrad_down_solar(ist,i) = &! meanrad_down_solar_ad(ist,i) + &
              transmission_aux%solar_path1%Tau_surf(ist,i)**2_jpim * &
              reflectance(i)%refl_out * refl_norm_scat * rad_ad%clear(i)

    transmission_aux_ad%solar_path1%Tau_surf(ist,i) = transmission_aux_ad%solar_path1%Tau_surf(ist,i) + &
              transmission_aux%solar_path1%Tau_surf(ist,i) * &
              reflectance(i)%refl_out * refl_norm_scat * 2._JPRB * rad_ad%clear(i) * &
              auxrad_stream%meanrad_down_solar(ist,i)
  endif

  if (thermal(i)) then
    do ist = 1, ircld%nstream(prof)

      ! add upward and downward parts
      auxrad_stream_ad%meanrad_up(ist,i) = auxrad_stream_ad%meanrad_up(ist,i) + &
                            auxrad_stream_ad%cloudy(ist,i)

      if (do_lambertian(i)) then
        auxrad_stream_ad%meanrad_down(ist,i) = &!auxrad_stream_ad%meanrad_down(ist,i) + &
          thermrefl(i) * transmission_aux%thermal_path1%tau_surf_p(ist,i) * &
          transmission_aux%thermal_path1%tau_surf(ist,i) * auxrad_stream_ad%cloudy(ist,i)

        thermrefl_ad(i) = thermrefl_ad(i) + &
          auxrad_stream%meanrad_down(ist,i) * &
          transmission_aux%thermal_path1%tau_surf_p(ist,i) * &
          transmission_aux%thermal_path1%tau_surf(ist,i) * auxrad_stream_ad%cloudy(ist,i)

        transmission_aux_ad%thermal_path1%tau_surf_p(ist,i) = &
          transmission_aux_ad%thermal_path1%tau_surf_p(ist,i) + &
          auxrad_stream%meanrad_down(ist,i) * thermrefl(i) * &
          transmission_aux%thermal_path1%tau_surf(ist,i) * auxrad_stream_ad%cloudy(ist,i)

        transmission_aux_ad%thermal_path1%tau_surf(ist,i) = &
          transmission_aux_ad%thermal_path1%tau_surf(ist,i) + &
          auxrad_stream%meanrad_down(ist,i) * thermrefl(i) * &
          transmission_aux%thermal_path1%tau_surf_p(ist,i) * auxrad_stream_ad%cloudy(ist,i)
      else
        auxrad_stream_ad%meanrad_down(ist,i) = & ! auxrad_stream_ad%meanrad_down(ist,i) + &
                                auxrad_stream_ad%cloudy(ist,i) * &
                                transmission_aux%thermal_path1%tau_surf(ist,i)**2 * thermrefl(i)

        transmission_aux_ad%thermal_path1%tau_surf(ist,i) = transmission_aux_ad%thermal_path1%tau_surf(ist,i) + &
                                              auxrad_stream_ad%cloudy(ist,i) * &
                                              2._JPRB * transmission_aux%thermal_path1%tau_surf(ist,i) * &
                                              thermrefl(i) * &
                                              auxrad_stream%meanrad_down(ist,i)

        thermrefl_ad(i) = thermrefl_ad(i) + &
                            auxrad_stream_ad%cloudy(ist,i) * transmission_aux%thermal_path1%tau_surf(ist,i)**2 &
                            * auxrad_stream%meanrad_down(ist,i)
      endif
    enddo
  endif

  if (solar(i)) then
    do ist = 1, ircld%nstream(prof)
      auxrad_stream_ad%meanrad_up_solar(ist,i) = auxrad_stream_ad%meanrad_up_solar(ist,i) + &
                                   auxrad_stream_ad%cloudy(ist,i)

      reflectance_ad(i)%refl_out = reflectance_ad(i)%refl_out + &
                transmission_aux%solar_path1%Tau_surf(ist,i)**2_jpim * refl_norm_scat * &
                auxrad_stream%meanrad_down_solar(ist,i) * auxrad_stream_ad%cloudy(ist,i)

      auxrad_stream_ad%meanrad_down_solar(ist,i) = &! meanrad_down_solar_ad(ist,i) + &
                transmission_aux%solar_path1%Tau_surf(ist,i)**2_jpim * &
                reflectance(i)%refl_out * refl_norm_scat * auxrad_stream_ad%cloudy(ist,i)

      transmission_aux_ad%solar_path1%Tau_surf(ist,i) = transmission_aux_ad%solar_path1%Tau_surf(ist,i) + &
                transmission_aux%solar_path1%Tau_surf(ist,i) * &
                reflectance(i)%refl_out * refl_norm_scat * 2._JPRB * auxrad_stream_ad%cloudy(ist,i) * &
                auxrad_stream%meanrad_down_solar(ist,i)
    enddo
  endif

enddo

!-------------------------------------------------------------------------------
!2b calculate clear-sky Rayleigh scattering contribution
!-------------------------------------------------------------------------------
if (dosolar) &
  call solar_rayleigh_ad(narray, iv3lev, opts, coef, solar, solar_spectrum, chanprof, &
                         raytracing, raytracing_ad, ircld, profiles, profiles_ad, &
                         profiles_dry, profiles_dry_ad, &
                         transmission_aux, transmission_aux_ad, auxrad_stream_ad)


!-------------------------------------------------------------------------------
!2a calculate near-surface layer contribution
!-------------------------------------------------------------------------------
if((opts%rt_ir%addaerosl .OR. opts%rt_ir%addclouds) .and. dosolar) &
     call solar_scattering_near_surf_ad(transmission_aux, transmission_aux_ad, &
                                        ircld, raytracing, raytracing_ad, transmission_scatt_ir_ad, &
                                        solar, solar_spectrum, reflectance%refl_out, sateqsun, iv3lay, iv3lev, &
                                        chanprof, narray, pfraction, auxrad_stream, auxrad_stream_ad)

if (dothermal) call calc_near_surf_contribution_ad(thermal, transmission_aux, transmission_aux_ad, auxrad, &
                       ircld, chanprof, narray, iv3lay, iv3lev, pol_id, rad_surfair_ad, rad_air_ad, auxrad_stream_ad)

!-------------------------------------
!2. calculate atmospheric contribution
!-------------------------------------
if((opts%rt_ir%addaerosl .OR. opts%rt_ir%addclouds) .and. dosolar) &
     call solar_scattering_air_ad(transmission_aux, transmission_aux_ad, ircld, raytracing, raytracing_ad, &
                                  transmission_scatt_ir_ad, solar, solar_spectrum, reflectance%refl_out, &
                                  sateqsun, chanprof, narray, auxrad_stream, auxrad_stream_ad)

if (dothermal) call calc_atmospheric_radiance_ad(thermal, transmission_aux, transmission_aux_ad, &
                                                 auxrad, ircld, chanprof, narray, rad_air_ad, auxrad_stream_ad)

rad_ad%clear(:) = 0._jprb

!----------------------------
!1. calculate layer radiances
!----------------------------
IF (dothermal) CALL rttov_calcrad_ad(chanprof, profiles, profiles_ad, coef, auxrad, &
                                     rad_skin_ad, rad_surfair_ad, rad_air_ad)

IF (LHOOK) CALL DR_HOOK('RTTOV_INTEGRATE_AD',1_jpim,ZHOOK_HANDLE)

contains

SUBROUTINE calc_atmospheric_radiance_ad(thermal, transmission_aux, transmission_aux_ad, &
                                        auxrad, ircld, chanprof, narray, rad_air_ad, auxrad_stream_ad)

  use rttov_types, Only : transmission_type_aux, radiance_aux, ircld_type, rttov_chanprof
  use rttov_const, only : min_tau
  use parkind1, only : jpim, jprb, jplm

  Implicit None

  logical(jplm), intent(in)                     :: thermal(:)
  type(transmission_type_aux), intent(in)       :: transmission_aux
  type(radiance_aux), intent(in)                :: auxrad
  type(ircld_type),     intent(in)              :: ircld
  type(rttov_chanprof), intent(in)              :: chanprof(:)
  integer(jpim), intent(in)                     :: narray(:)

  Real(jprb), intent(inout)                     :: rad_air_ad(:,:)
  type(transmission_type_aux), intent(inout)    :: transmission_aux_ad
  type(radiance_aux), intent(inout)             :: auxrad_stream_ad

  integer(jpim) :: nlevels, nlayers, nchannels
  integer(jpim) :: i, ist, lay, lev, levm1

  Real(jprb) :: tau_prod
  Real(jprb) :: b1_ad, b2_ad, b3_ad
  Real(jprb) :: temp
  REAL(jprb) :: z(3)

  !unpack narray
  nlevels = narray(1)
  nlayers = nlevels - 1_jpim
  nchannels = narray(3)

#define tau_lev transmission_aux%thermal_path1%Tau_level(lay+1,ist,i)
#define tau_levm1 transmission_aux%thermal_path1%Tau_level(lay,ist,i)
#define tau_lev_r transmission_aux%thermal_path1%Tau_level_r(lay+1,ist,i)
#define tau_levm1_r transmission_aux%thermal_path1%Tau_level_r(lay,ist,i)
#define tausun_lev transmission_aux%solar_path2%Tau_level(lay+1,ist,i)
#define tausun_levm1 transmission_aux%solar_path2%Tau_level(lay,ist,i)

#define dtau_lay (tau_levm1 - tau_lev)
#define daux_lay (auxrad%air(lev,i) - auxrad%air(levm1,i))
#define B1_2 auxrad%air(levm1,i) * dtau_lay
#define B2_2 daux_lay * tau_lev
#define B3_2 daux_lay * dtau_lay * transmission_aux%thermal_path1%od_singlelayer_r(isti,lay,i)

  DO i = 1, nchannels
    IF (thermal(i)) THEN
      DO ist = 0, ircld%nstream(prof)

        do lay = nlayers, 2, -1
          if(auxrad_stream%down_ref(lay,ist,i) < 0._jprb) &
             auxrad_stream_ad%down(lay,ist,i) = 0._jprb

          auxrad_stream_ad%up(lay-1,ist,i) = auxrad_stream_ad%up(lay-1,ist,i) + auxrad_stream_ad%up(lay,ist,i)
          auxrad_stream_ad%down(lay-1,ist,i) = auxrad_stream_ad%down(lay-1,ist,i) + auxrad_stream_ad%down(lay,ist,i)
        enddo

        lay = 1
        if(auxrad_stream%down_ref(lay,ist,i) < 0._jprb) then
           auxrad_stream_ad%down(lay,ist,i) = 0._jprb
        endif

        DO lay = 1, nlayers
           levm1 = lay
           lev = lay + 1
           isti = ircld%icldarr(ist,lay,prof)
           IF(transmission_aux%anynegtau > 0._jprb) THEN

              if(lay > 1 .and. tau_levm1 >= 0.and. &
                 tau_lev < 0) then

                 If (tau_lev > min_tau) Then
                    tau_prod = tau_lev_r * tau_levm1_r
                    auxrad_stream_ad%up(lay,ist,i) = auxrad_stream_ad%up(lay,ist,i) + &
                         auxrad_stream_ad%down(lay,ist,i) * tau_prod

                    transmission_aux_ad%thermal_path1%tau_level(lev,ist,i) = &
                         transmission_aux_ad%thermal_path1%tau_level(lev,ist,i) - &
                         auxrad_stream_ad%down(lay,ist,i) * &
                         tau_levm1 * (tau_prod*tau_prod) * &
                         0.5_JPRB * (auxrad%air(lev,i) + auxrad%air(levm1,i)) * &
                         (tau_levm1 - &
                         tau_lev )

                    transmission_aux_ad%thermal_path1%tau_level(levm1,ist,i) = &
                         transmission_aux_ad%thermal_path1%tau_level(levm1,ist,i) - &
                         auxrad_stream_ad%down(lay,ist,i) * &
                         tau_lev * (tau_prod*tau_prod) * &
                         0.5_JPRB * (auxrad%air(lev,i) + auxrad%air(levm1,i)) * &
                         (tau_levm1 - &
                         tau_lev)
                 End If

                 rad_air_ad(lev,i) = rad_air_ad(lev,i) + &
                      auxrad_stream_ad%up(lay,ist,i) * 0.5_JPRB * &
                      (tau_levm1 - &
                      tau_lev)

                 rad_air_ad(levm1,i) = rad_air_ad(levm1,i) + &
                      auxrad_stream_ad%up(lay,ist,i) * 0.5_JPRB * &
                      (tau_levm1 - &
                      tau_lev)

                 transmission_aux_ad%thermal_path1%tau_level(levm1,ist,i) = &
                      transmission_aux_ad%thermal_path1%tau_level(levm1,ist,i) + &
                      auxrad_stream_ad%up(lay,ist,i) * 0.5_JPRB * (auxrad%air(lev,i) + auxrad%air(levm1,i))

                 transmission_aux_ad%thermal_path1%tau_level(lev,ist,i) = &
                      transmission_aux_ad%thermal_path1%tau_level(lev,ist,i) - &
                      auxrad_stream_ad%up(lay,ist,i) * 0.5_JPRB * (auxrad%air(lev,i) + auxrad%air(levm1,i))
              else
                 temp = transmission_aux%fac1(lay,ist,i) * &
                        transmission_aux%thermal_path1%fac2(lay+1,ist,i) * auxrad_stream_ad%down(lay,ist,i)

                 B1_AD = temp * (tau_levm1_r * tau_lev_r)
                 B2_AD = temp * (tau_lev_r)**2
                 B3_AD = -temp * (tau_levm1_r * tau_lev_r)

                 B1_AD = B1_AD + auxrad_stream_ad%up(lay,ist,i)
                 B2_AD = B2_AD - auxrad_stream_ad%up(lay,ist,i)
                 B3_AD = B3_AD + auxrad_stream_ad%up(lay,ist,i)

                 transmission_aux_ad%thermal_path1%tau_level(levm1,ist,i) = &
                      transmission_aux_ad%thermal_path1%tau_level(levm1,ist,i) + &
                      transmission_aux%fac1(lay,ist,i) * (&
                      B1_AD * auxrad%air(levm1,i) - &
                      temp * &
                      (B1_2 - &
                      B3_2) * &
                      (tau_lev_r * &
                      tau_levm1_r**2) + &
                      B3_AD * &
                      daux_lay * &
                      transmission_aux%thermal_path1%od_singlelayer_r(isti,lay,i))

                 transmission_aux_ad%thermal_path1%tau_level(lev,ist,i) = &
                      transmission_aux_ad%thermal_path1%tau_level(lev,ist,i) + &
                      transmission_aux%fac1(lay,ist,i) * (&
                      -temp * &
                      (B1_2 - &
                      B3_2) * &
                      (tau_lev_r**2 * &
                      tau_levm1_r) - &
                      B3_AD * &
                      daux_lay * &
                      transmission_aux%thermal_path1%od_singlelayer_r(isti,lay,i) - &
                      temp * &
                      2.0_jprb * B2_2 * &
                      (tau_lev_r**3) - &
                      B1_AD * &
                      auxrad%air(levm1,i) + &
                      B2_AD * &
                      daux_lay)

                 transmission_aux_ad%thermal_path1%od_singlelayer(isti,lay,i) = &
                      transmission_aux_ad%thermal_path1%od_singlelayer(isti,lay,i) - &
                      transmission_aux%fac1(lay,ist,i) * &
                      B3_AD * &
                      B3_2 * &
                      transmission_aux%thermal_path1%od_singlelayer_r(isti,lay,i)

                 rad_air_ad(levm1,i) = rad_air_ad(levm1,i) + transmission_aux%fac1(lay,ist,i) * (&
                      -B3_AD * &
                      dtau_lay * &
                      transmission_aux%thermal_path1%od_singlelayer_r(isti,lay,i) + &
                      B1_AD * &
                      dtau_lay - &
                      B2_AD * &
                      tau_lev)

                 rad_air_ad(lev,i) = rad_air_ad(lev,i) + transmission_aux%fac1(lay,ist,i) * (&
                      B3_AD * &
                      dtau_lay * &
                      transmission_aux%thermal_path1%od_singlelayer_r(isti,lay,i) + &
                      B2_AD * &
                      tau_lev)
             endif
           ELSE
             IF (transmission_aux%fac1(lay,ist,i) > 0._jprb) THEN           
               IF (do_lambertian(i)) THEN
                 temp  = 0._jprb
                 B1_AD = 0._jprb
                 B2_AD = 0._jprb
                 B3_AD = 0._jprb
                 IF(transmission_aux%fac1(lay,ist,i) * transmission_aux%thermal_path1%fac2(lay+1,ist,i) > 0._jprb) THEN
                   z(1) = 0.5_jprb * &
                     (transmission_aux%thermal_path1%tau_level_p(lay,ist,i) - &
                      transmission_aux%thermal_path1%tau_level_p(lev,ist,i)) * &
                      transmission_aux%thermal_path1%tau_level_p_r(lev,ist,i) * &
                      transmission_aux%thermal_path1%tau_level_p_r(lay,ist,i) * &
                      auxrad_stream_ad%down(lay,ist,i)
                   
                   rad_air_ad(lay,i)  = rad_air_ad(lay,i) + z(1)
                   rad_air_ad(lay+1,i) = rad_air_ad(lay+1,i) + z(1)
                   
                   z(2) = 0.5_jprb * (auxrad%air(lay,i) + auxrad%air(lay+1,i)) * &
                     (transmission_aux%thermal_path1%tau_level_p_r(lev,ist,i) * &
                      transmission_aux%thermal_path1%tau_level_p_r(lay,ist,i)) * &
                      auxrad_stream_ad%down(lay,ist,i)
                   
                   transmission_aux_ad%thermal_path1%tau_level_p(lay,ist,i) = &
                     transmission_aux_ad%thermal_path1%tau_level_p(lay,ist,i) + z(2)

                   transmission_aux_ad%thermal_path1%tau_level_p(lev,ist,i) = &
                     transmission_aux_ad%thermal_path1%tau_level_p(lev,ist,i) - z(2)
                  
                   z(3) = 0.5_jprb * (auxrad%air(lay,i) + auxrad%air(lay+1,i)) * &
                     (transmission_aux%thermal_path1%tau_level_p(lay,ist,i) - &
                      transmission_aux%thermal_path1%tau_level_p(lev,ist,i)) * &
                      auxrad_stream_ad%down(lay,ist,i)
                   
                   transmission_aux_ad%thermal_path1%tau_level_p_r(lev,ist,i) = &
                     transmission_aux_ad%thermal_path1%tau_level_p_r(lev,ist,i) + &
                     z(3) * transmission_aux%thermal_path1%tau_level_p_r(lay,ist,i)
                   
                   transmission_aux_ad%thermal_path1%tau_level_p_r(lay,ist,i) = &
                     transmission_aux_ad%thermal_path1%tau_level_p_r(lay,ist,i) + &
                     z(3) * transmission_aux%thermal_path1%tau_level_p_r(lev,ist,i)                   
                 ENDIF
               ELSE
                 temp = transmission_aux%thermal_path1%fac2(lay+1,ist,i) * &
                        auxrad_stream_ad%down(lay,ist,i)
                 B1_AD = temp  * (tau_levm1_r * tau_lev_r)
                 B2_AD = temp * (tau_lev_r)**2
                 B3_AD = -temp * (tau_levm1_r * tau_lev_r)
               ENDIF

               B1_AD = B1_AD + auxrad_stream_ad%up(lay,ist,i)
               B2_AD = B2_AD - auxrad_stream_ad%up(lay,ist,i)
               B3_AD = B3_AD + auxrad_stream_ad%up(lay,ist,i)
               
               transmission_aux_ad%thermal_path1%tau_level(levm1,ist,i) = &
                 transmission_aux_ad%thermal_path1%tau_level(levm1,ist,i) + &
                 B1_AD * auxrad%air(levm1,i) - &
                 temp * &
                 (B1_2 - &
                 B3_2) * &
                 (tau_lev_r * &
                 tau_levm1_r**2) + &
                 B3_AD * daux_lay * transmission_aux%thermal_path1%od_singlelayer_r(isti,lay,i)
               
               transmission_aux_ad%thermal_path1%tau_level(lev,ist,i) = &
                 transmission_aux_ad%thermal_path1%tau_level(lev,ist,i) - &
                 temp * &
                 (B1_2 - &
                 B3_2) * &
                 (tau_lev_r**2 * tau_levm1_r) - &
                 B3_AD * daux_lay * &
                 transmission_aux%thermal_path1%od_singlelayer_r(isti,lay,i) - &
                 temp * 2.0_jprb * B2_2 * &
                 (tau_lev_r**3) - &
                 B1_AD * auxrad%air(levm1,i) + &
                 B2_AD * daux_lay
               
               transmission_aux_ad%thermal_path1%od_singlelayer(isti,lay,i) = &
                 transmission_aux_ad%thermal_path1%od_singlelayer(isti,lay,i) - &
                 B3_AD * &
                 B3_2 * &
                 transmission_aux%thermal_path1%od_singlelayer_r(isti,lay,i)
               
               rad_air_ad(levm1,i) = rad_air_ad(levm1,i) - &
                 B3_AD * dtau_lay * transmission_aux%thermal_path1%od_singlelayer_r(isti,lay,i) + &
                 B1_AD * dtau_lay - &
                 B2_AD * tau_lev
               
               rad_air_ad(lev,i) = rad_air_ad(lev,i) + &
                 B3_AD * dtau_lay * transmission_aux%thermal_path1%od_singlelayer_r(isti,lay,i) + &
                 B2_AD * tau_lev
             ENDIF
           ENDIF
         ENDDO
       END DO
     ENDIF
   END DO
 END SUBROUTINE calc_atmospheric_radiance_ad

subroutine calc_near_surf_contribution_ad(thermal, transmission_aux, transmission_aux_ad, auxrad, ircld, &
                                          chanprof, narray, iv3lay, iv3lev, pol_id, &
                                          rad_surfair_ad, rad_air_ad, auxrad_stream_ad)

  use rttov_const, only : min_tau, min_od
  use rttov_types, Only : rttov_chanprof, ircld_type, transmission_type_aux, radiance_aux
  use parkind1, only : jpim, jprb, jplm

  Implicit None

  logical(jplm), intent(in)                  :: thermal(:)
  type(transmission_type_aux), intent(in)    :: transmission_aux
  type(radiance_aux), intent(in)             :: auxrad
  type(ircld_type), intent(in)               :: ircld
  type(rttov_chanprof), intent(in)           :: chanprof(:)
  integer(jpim), intent(in)                  :: narray(:)
  integer(jpim), intent(in)                  :: iv3lay(:)
  integer(jpim), intent(in)                  :: iv3lev(:)
  integer(jpim), intent(in)                  :: pol_id(:)

  type(transmission_type_aux), intent(inout) :: transmission_aux_ad
  real(jprb), intent(inout)                  :: rad_surfair_ad(:), rad_air_ad(:,:)
  type(radiance_aux), intent(inout)          :: auxrad_stream_ad

!local
  integer(jpim) :: nlayers, nlevels, nstreams, nchannels

  real(jprb)    :: b1_3, b2_3, b3_3
  real(jprb)    :: b1_ad, b2_ad, b3_ad

  integer(jpim) :: i, lev, lay, ist
  REAL(jprb)    :: z(3)

! unpack
  nlevels = narray(1)
  nlayers = nlevels - 1_jpim
  nstreams = narray(2)
  nchannels = narray(3)

  Do i = 1, nchannels
    if (thermal(i)) then
      lay = iv3lay(i)
      lev = iv3lev(i)

      do ist = 0, ircld%nstream(prof) !rev loop not necessary
        ! assume there is no atmospheric source term for 3rd/4th stokes vector elements
        if (pol_id(i) >= 6_jpim) auxrad_stream_ad%meanrad_up(ist,i) = 0.0_jprb

        auxrad_stream_ad%up(lay,ist,i) = auxrad_stream_ad%up(lay,ist,i) + &
                            auxrad_stream_ad%meanrad_up(ist,i)

        auxrad_stream_ad%down(lay,ist,i) = auxrad_stream_ad%down(lay,ist,i) + &
                              auxrad_stream_ad%meanrad_down(ist,i)

        if(transmission_aux%thermal_path1%tau_surf(ist,i) < 0._JPRB) then
           if(transmission_aux%thermal_path1%tau_level(lev,ist,i) >= 0._JPRB) then ! nearly always true?

              rad_surfair_ad(i) = rad_surfair_ad(i) + &
                                  0.5_JPRB * auxrad_stream_ad%meanrad_up(ist,i) * &
                                 (transmission_aux%thermal_path1%tau_level(lev,ist,i) - &
                                  transmission_aux%thermal_path1%tau_surf(ist,i))

              rad_air_ad(lev,i) = rad_air_ad(lev,i) + &
                                  0.5_JPRB * auxrad_stream_ad%meanrad_up(ist,i)* &
                                 (transmission_aux%thermal_path1%tau_level(lev,ist,i) - &
                                  transmission_aux%thermal_path1%tau_surf(ist,i))

              transmission_aux_ad%thermal_path1%tau_level(lev,ist,i) = &
                                  transmission_aux_ad%thermal_path1%tau_level(lev,ist,i) + &
                                  0.5_JPRB * auxrad_stream_ad%meanrad_up(ist,i) * &
                                 (auxrad%surfair(i) + auxrad%air(lev,i))

              transmission_aux_ad%thermal_path1%tau_surf(ist,i) = &
                                  transmission_aux_ad%thermal_path1%tau_surf(ist,i) - &
                                  0.5_JPRB * auxrad_stream_ad%meanrad_up(ist,i) * &
                                 (auxrad%surfair(i) + auxrad%air(lev,i))
           endif
        else
           if(transmission_aux%thermal_path1%od_sfrac(ist,i) < min_od .OR. &
              ((transmission_aux%thermal_path1%tau_level(lev,ist,i) - &
                transmission_aux%thermal_path1%tau_surf(ist,i)) < min_od)) THEN
           else
              IF (transmission_aux%thermal_path1%tau_surf(ist,i) > min_tau) THEN
                if (do_lambertian(i)) then
                  z(1) = 0.5_jprb * (transmission_aux%thermal_path1%tau_level_p(lev,ist,i) - &
                                     transmission_aux%thermal_path1%tau_surf_p(ist,i)) * &
                                     transmission_aux%thermal_path1%tau_level_p_r(lev,ist,i) * &
                                     transmission_aux%thermal_path1%tau_surf_p_r(ist,i) * &
                                     auxrad_stream_ad%meanrad_down(ist,i)
                   
                  rad_surfair_ad(i) = rad_surfair_ad(i) + z(1)
                  rad_air_ad(lev,i) = rad_air_ad(lev,i) + z(1)
                   
                  z(2) = 0.5_jprb * (auxrad%surfair(i) + auxrad%air(lev,i)) * &
                    (transmission_aux%thermal_path1%tau_level_p_r(lev,ist,i) * &
                     transmission_aux%thermal_path1%tau_surf_p_r(ist,i)) * &
                    auxrad_stream_ad%meanrad_down(ist,i)
                   
                  transmission_aux_ad%thermal_path1%tau_surf_p(ist,i) = &
                    transmission_aux_ad%thermal_path1%tau_surf_p(ist,i) - z(2)

                  transmission_aux_ad%thermal_path1%tau_level_p(lev,ist,i) = &
                    transmission_aux_ad%thermal_path1%tau_level_p(lev,ist,i) + &
                    z(2)
                   
                  z(3) = 0.5_jprb * (auxrad%surfair(i) + auxrad%air(lev,i)) * &
                    (transmission_aux%thermal_path1%tau_level_p(lev,ist,i) - &
                     transmission_aux%thermal_path1%tau_surf_p(ist,i)) * &
                    auxrad_stream_ad%meanrad_down(ist,i)
                   
                  transmission_aux_ad%thermal_path1%tau_level_p_r(lev,ist,i) = &
                    transmission_aux_ad%thermal_path1%tau_level_p_r(lev,ist,i) + &
                    z(3) * transmission_aux%thermal_path1%tau_surf_p_r(ist,i)

                  transmission_aux_ad%thermal_path1%tau_surf_p_r(ist,i) = &
                    transmission_aux_ad%thermal_path1%tau_surf_p_r(ist,i) + &
                    z(3) * transmission_aux%thermal_path1%tau_level_p_r(lev,ist,i)

                  B1_AD = 0.0_JPRB
                  B2_AD = 0.0_JPRB
                  B3_AD = 0.0_JPRB

                else
  !direct calc
                   B1_3 = auxrad%air(lev,i) * &
                         (transmission_aux%thermal_path1%tau_level(lev,ist,i) - &
                         transmission_aux%thermal_path1%tau_surf(ist,i))

                   B2_3 = (auxrad%surfair(i) - auxrad%air(lev,i)) * transmission_aux%thermal_path1%tau_surf(ist,i)

                   B3_3 = (auxrad%surfair(i) - auxrad%air(lev,i)) * &
                          (transmission_aux%thermal_path1%tau_level(lev,ist,i) - &
                          transmission_aux%thermal_path1%tau_surf(ist,i)) * &
                          (transmission_aux%thermal_path1%od_sfrac_r(ist,i))
                   !first use B1,b2,b3c
                   B1_AD = auxrad_stream_ad%meanrad_down(ist,i) * &
                           (transmission_aux%thermal_path1%tau_level_r(lev,ist,i) * &
                           transmission_aux%thermal_path1%tau_surf_r(ist,i))
                   B3_AD = -auxrad_stream_ad%meanrad_down(ist,i) * &
                           (transmission_aux%thermal_path1%tau_level_r(lev,ist,i) * &
                           transmission_aux%thermal_path1%tau_surf_r(ist,i))
                   B2_AD = auxrad_stream_ad%meanrad_down(ist,i) * &
                           (transmission_aux%thermal_path1%tau_surf_r(ist,i))**2

                   transmission_aux_ad%thermal_path1%tau_surf(ist,i) = &
                         transmission_aux_ad%thermal_path1%tau_surf(ist,i) - &
                         auxrad_stream_ad%meanrad_down(ist,i) * 2._jprb * B2_3 * &
                         transmission_aux%thermal_path1%tau_surf_r(ist,i)**3_jpim

                   transmission_aux_ad%thermal_path1%tau_surf(ist,i) = &
                         transmission_aux_ad%thermal_path1%tau_surf(ist,i) - &
                         auxrad_stream_ad%meanrad_down(ist,i) * (B1_3 - B3_3) * &
                         transmission_aux%thermal_path1%tau_level_r(lev,ist,i) * &
                         transmission_aux%thermal_path1%tau_surf_r(ist,i)**2_jpim

                   transmission_aux_ad%thermal_path1%tau_level(lev,ist,i) = &
                         transmission_aux_ad%thermal_path1%tau_level(lev,ist,i) - &
                         auxrad_stream_ad%meanrad_down(ist,i) * (B1_3 - B3_3) * &
                         transmission_aux%thermal_path1%tau_level_r(lev,ist,i)**2_jpim * &
                         transmission_aux%thermal_path1%tau_surf_r(ist,i)
                endif
              ELSE
                 B1_AD = 0.0_JPRB
                 B2_AD = 0.0_JPRB
                 B3_AD = 0.0_JPRB
              ENDIF

              B1_AD = B1_AD + auxrad_stream_ad%meanrad_up(ist,i)
              B2_AD = B2_AD - auxrad_stream_ad%meanrad_up(ist,i)
              B3_AD = B3_AD + auxrad_stream_ad%meanrad_up(ist,i)

              transmission_aux_ad%thermal_path1%od_sfrac(ist,i) = &
                                                    transmission_aux_ad%thermal_path1%od_sfrac(ist,i) - &
                                                    B3_AD * (auxrad%surfair(i)-auxrad%air(lev,i)) * &
                                                   (transmission_aux%thermal_path1%tau_level(lev,ist,i) - &
                                                    transmission_aux%thermal_path1%tau_surf(ist,i)) * &
                                                   (transmission_aux%thermal_path1%od_sfrac_r(ist,i)**2)

              transmission_aux_ad%thermal_path1%tau_surf(ist,i) = &
                                                    transmission_aux_ad%thermal_path1%tau_surf(ist,i) - &
                                                    B3_AD * (auxrad%surfair(i)-auxrad%air(lev,i)) * &
                                                    transmission_aux%thermal_path1%od_sfrac_r(ist,i) + &
                                                    B2_AD * (auxrad%surfair(i)-auxrad%air(lev,i)) - &
                                                    B1_AD * auxrad%air(lev,i)

              transmission_aux_ad%thermal_path1%tau_level(lev,ist,i) = &
                                                    transmission_aux_ad%thermal_path1%tau_level(lev,ist,i) + &
                                                    B3_AD * (auxrad%surfair(i) - auxrad%air(lev,i)) * &
                                                    transmission_aux%thermal_path1%od_sfrac_r(ist,i) + &
                                                    B1_AD * auxrad%air(lev,i)

              rad_air_ad(lev,i) = rad_air_ad(lev,i) - &
                   B3_AD * (transmission_aux%thermal_path1%tau_level(lev,ist,i) - &
                   transmission_aux%thermal_path1%tau_surf(ist,i)) * &
                   transmission_aux%thermal_path1%od_sfrac_r(ist,i) - &
                   B2_AD * transmission_aux%thermal_path1%tau_surf(ist,i) + &
                   B1_AD * (transmission_aux%thermal_path1%tau_level(lev,ist,i) - &
                   transmission_aux%thermal_path1%tau_surf(ist,i))

              rad_surfair_ad(i) = rad_surfair_ad(i) + &
                   B3_AD * (transmission_aux%thermal_path1%tau_level(lev,ist,i) - &
                   transmission_aux%thermal_path1%tau_surf(ist,i)) * &
                   transmission_aux%thermal_path1%od_sfrac_r(ist,i) + &
                   B2_AD * transmission_aux%thermal_path1%tau_surf(ist,i)
           endif
        endif
      Enddo
    endif
  End do
end subroutine calc_near_surf_contribution_ad

subroutine solar_scattering_air_ad(transmission_aux, transmission_aux_ad, ircld, raytracing, raytracing_ad, &
                                   transmission_scatt_ir_ad, solar, solar_spectrum, refl, sateqsun, &
                                   chanprof, narray, auxrad_stream, auxrad_stream_ad)

  use rttov_const, only : min_tau, z4pi_r
  use rttov_types, Only : transmission_type_aux, ircld_type, radiance_aux, &
                          rttov_chanprof, transmission_scatt_ir_type, raytracing_type
  use parkind1, only : jpim, jprb, jplm 

  Implicit None
  
  type(transmission_type_aux), intent(in)         :: transmission_aux
  type(ircld_type),     intent(in)                :: ircld
  type(raytracing_type), intent(in)               :: raytracing
  logical(jplm), intent(in)                       :: solar(:)
  real(jprb), intent(in)                          :: solar_spectrum(:)
  real(jprb), intent(in)                          :: refl(:)
  logical(jplm), intent(in)                       :: sateqsun(:,:)
  type(rttov_chanprof), intent(in)                :: chanprof(:)
  integer(jpim), intent(in)                       :: narray(:)
  type(radiance_aux), intent(in)                  :: auxrad_stream

  type(transmission_type_aux), intent(inout)      :: transmission_aux_ad
  type(raytracing_type), intent(inout)            :: raytracing_ad
  type(transmission_scatt_ir_type), intent(inout) :: transmission_scatt_ir_ad
  type(radiance_aux), intent(inout)               :: auxrad_stream_ad

  integer(jpim) :: nlevels, nlayers, nstreams, nchannels
  integer(jpim) :: i, ist, lay, lev, levm1
  
  Real(jprb) :: fac_2_ad(7)
  Real(jprb) :: temp(3,narray(1)-1)

! unpack
  nlevels = narray(1)
  nlayers = nlevels - 1_jpim
  nstreams = narray(2)
  nchannels = narray(3)

  do i = 1, nchannels
     if(solar(i)) then

        do ist = 0, ircld%nstream(prof) !rev loop not necessary

          if (refl(i) > 0._jprb) then

            do lay = nlayers, 2, -1
              if(auxrad_stream%down_ref_solar(lay,ist,i) < 0._jprb) &
                 auxrad_stream_ad%down_solar(lay,ist,i) = 0._jprb

              auxrad_stream_ad%down_solar(lay-1,ist,i) = auxrad_stream_ad%down_solar(lay-1,ist,i) + &
                                                         auxrad_stream_ad%down_solar(lay,ist,i) * &
                                                         transmission_aux%solar_path1%fac2(lay+1,ist,i)
            enddo


            temp(1,1:nlayers) = auxrad_stream%fac6_2(1:nlayers,ist,i) * auxrad_stream%fac2_2(1:nlayers,ist,i) * &
                                auxrad_stream_ad%down_solar(1:nlayers,ist,i)
            temp(2,1:nlayers) = transmission_aux%solar_path1%tau_level_r(1:nlayers,ist,i) * &
                              transmission_aux%solar_path1%tau_level_r(2:nlayers+1,ist,i)
           endif

           do lay = nlayers, 2, -1
             auxrad_stream_ad%up_solar(lay-1,ist,i) = auxrad_stream_ad%up_solar(lay-1,ist,i) + &
                                                      auxrad_stream_ad%up_solar(lay,ist,i)
           enddo

           do lay = 1, nlayers !rev loop not necessary

              levm1 = lay
              lev   = lay + 1
              isti = ircld%icldarr(ist,lay,prof)

              fac_2_ad(:) = 0.0_jprb ! DAR: if errors then may need to zero everything. Not just 2 and 6

              if (refl(i) > 0._jprb) then

                if(lay == 1 .or. transmission_aux%solar_path1%Tau_level(lay+1,ist,i) > min_tau) then
                   !-------------------Downward single scattering of the solar beam--------------------------
                   if(.not. sateqsun(lay,prof)) then
                      temp(3,lay) = (auxrad_stream%fac5_2(lay,ist,i) - auxrad_stream%fac4_2(lay,ist,i)) * &
                                     auxrad_stream%fac7_2(lay,ist,i)

                      if(lay > 1_jpim) then
                         transmission_aux_ad%solar_path1%tau_level(levm1,ist,i) = &
                              transmission_aux_ad%solar_path1%tau_level(levm1,ist,i) - &
                              temp(1,lay) * temp(2,lay)**2 * temp(3,lay) * tausun_levm1 * &
                              transmission_aux%solar_path1%Tau_level(lay+1,ist,i)

                         transmission_aux_ad%solar_path2%tau_level(levm1,ist,i) = &
                              transmission_aux_ad%solar_path2%tau_level(levm1,ist,i) + &
                              temp(1,lay) * temp(2,lay) * temp(3,lay)
                      endif

                      fac_2_ad(4) = &!fac4_2_ad
                                  - temp(1,lay) * temp(2,lay) * &
                                  auxrad_stream%fac7_2(lay,ist,i) * transmission_aux%solar_path2%tau_level(lay,ist,i)

                      fac_2_ad(5) = &!fac5_2_ad + &
                                  temp(1,lay) * temp(2,lay) * &
                                  auxrad_stream%fac7_2(lay,ist,i) * transmission_aux%solar_path2%tau_level(lay,ist,i)

                      fac_2_ad(7) = &!fac7_2_ad + &
                                  temp(1,lay) * temp(2,lay) * &
                                 (auxrad_stream%fac5_2(lay,ist,i) - auxrad_stream%fac4_2(lay,ist,i)) * &
                                  transmission_aux%solar_path2%tau_level(lay,ist,i)

                      raytracing_ad%pathsat(lay,prof) = raytracing_ad%pathsat(lay,prof) + &
                           fac_2_ad(7) * raytracing%pathsun(lay,prof) / &
                          (raytracing%pathsun(lay,prof) - raytracing%pathsat(lay,prof))**2

                      raytracing_ad%pathsun(lay,prof) = raytracing_ad%pathsun(lay,prof) - &
                           fac_2_ad(7) * raytracing%pathsat(lay,prof) / &
                          (raytracing%pathsun(lay,prof) - raytracing%pathsat(lay,prof))**2
                   else
                      temp(3,lay) = auxrad_stream%fac4_2(lay,ist,i) * &
                                    transmission_aux%solar_path2%od_singlelayer(isti,lay,i)

                      if(lay > 1) then
  !                       transmission_aux_ad%solar_path1%tau_level(levm1,ist,i) = &
  !                            transmission_aux_ad%solar_path1%tau_level(levm1,ist,i) - &
  !                            temp(1,lay) * temp(2,lay) * temp(3,lay)
  !DAR 3108 think this should be temp22 * temp3 * tau_lev as above from integrate_tl

                         transmission_aux_ad%solar_path1%tau_level(levm1,ist,i) = &
                              transmission_aux_ad%solar_path1%tau_level(levm1,ist,i) - &
                              temp(1,lay) * temp(2,lay)**2_jpim * temp(3,lay) * &
                              transmission_aux%solar_path1%Tau_level(lay+1,ist,i)

                         transmission_aux_ad%solar_path2%tau_level(levm1,ist,i) = &
                              transmission_aux_ad%solar_path2%tau_level(levm1,ist,i) + &
                              temp(1,lay) * temp(2,lay) * temp(3,lay)
                      endif

                      transmission_aux_ad%solar_path2%od_singlelayer(isti,lay,i) = &
                           transmission_aux_ad%solar_path2%od_singlelayer(isti,lay,i) + &
                           temp(1,lay) * temp(2,lay) * auxrad_stream%fac4_2(lay,ist,i) * tausun_levm1

                      fac_2_ad(4) = &!fac4_2_ad + &
                                  temp(1,lay) * temp(2,lay) * &
                                  transmission_aux%solar_path2%tau_level(lay,ist,i) * &
                                  transmission_aux%solar_path2%od_singlelayer(isti,lay,i)
                   endif

                   transmission_aux_ad%solar_path1%tau_level(lev,ist,i) = &
                        transmission_aux_ad%solar_path1%tau_level(lev,ist,i) - &
                        temp(1,lay) * temp(3,lay) * temp(2,lay)**2_jpim * tausun_levm1 * &
                        transmission_aux%solar_path1%Tau_level(lay,ist,i)

                   fac_2_ad(2) = &!fac2_2_ad + &
                               auxrad_stream_ad%down_solar(lay,ist,i) * auxrad_stream%fac6_2(lay,ist,i) * &
                               temp(2,lay) * temp(3,lay) * tausun_levm1

                   fac_2_ad(6) = &!fac6_2_ad + &
                               auxrad_stream_ad%down_solar(lay,ist,i) * auxrad_stream%fac2_2(lay,ist,i) * &
                               temp(2,lay) * temp(3,lay) * tausun_levm1
                else
                   auxrad_stream_ad%down_solar(lay,ist,i) = 0._JPRB
                endif  ! min_tau

              endif ! refl(i) > 0.

       !----------------Upward single scattering of the solar beam-------------------------------

              isti = ircld%icldarr(ist,lay,prof)

              fac_2_ad(1) = &!fac1_2_ad + &
                          auxrad_stream_ad%up_solar(lay,ist,i) * &
                          auxrad_stream%fac2_2(lay,ist,i) * auxrad_stream%fac3_2(lay,ist,i) * &
                         (tausun_levm1 - tausun_lev)

              fac_2_ad(2) = fac_2_ad(2) + &
                          auxrad_stream_ad%up_solar(lay,ist,i) * &
                          auxrad_stream%fac1_2(lay,ist,i) * auxrad_stream%fac3_2(lay,ist,i) * &
                         (tausun_levm1 - tausun_lev)

              fac_2_ad(3) = &!fac3_2_ad + &
                          auxrad_stream_ad%up_solar(lay,ist,i) * &
                          auxrad_stream%fac1_2(lay,ist,i) * auxrad_stream%fac2_2(lay,ist,i) * &
                         (tausun_levm1 - tausun_lev)

              if(lay > 1) then
                 transmission_aux_ad%solar_path2%tau_level(levm1,ist,i) = &
                      transmission_aux_ad%solar_path2%tau_level(levm1,ist,i) + &
                      auxrad_stream_ad%up_solar(lay,ist,i) * &
                      auxrad_stream%fac1_2(lay,ist,i) * auxrad_stream%fac2_2(lay,ist,i) * &
                      auxrad_stream%fac3_2(lay,ist,i)
              endif

              transmission_aux_ad%solar_path2%tau_level(lev,ist,i) = &
                  transmission_aux_ad%solar_path2%tau_level(lev,ist,i) - &
                  auxrad_stream_ad%up_solar(lay,ist,i) * &
                  auxrad_stream%fac1_2(lay,ist,i) * auxrad_stream%fac2_2(lay,ist,i) * auxrad_stream%fac3_2(lay,ist,i)

              transmission_aux_ad%solar_path1%od_singlelayer(isti,lay,i) = &
                  transmission_aux_ad%solar_path1%od_singlelayer(isti,lay,i) - &
                  fac_2_ad(5) * auxrad_stream%fac5_2(lay,ist,i) !exp(od_singlelayer)

              transmission_aux_ad%solar_path2%od_singlelayer(isti,lay,i) = &
                   transmission_aux_ad%solar_path2%od_singlelayer(isti,lay,i) - &
                   fac_2_ad(4) * auxrad_stream%fac4_2(lay,ist,i) !exp(od_sunsinglelayer)

              raytracing_ad%pathsat(lay,prof) = raytracing_ad%pathsat(lay,prof) + &
                   fac_2_ad(3) / raytracing%patheff(lay,prof)

              raytracing_ad%patheff(lay,prof) = raytracing_ad%patheff(lay,prof) - &
                   fac_2_ad(3) * raytracing%pathsat(lay,prof) / raytracing%patheff(lay,prof) ** 2

             transmission_scatt_ir_ad%ssa(isti,lay,i) = &
                  transmission_scatt_ir_ad%ssa(isti,lay,i) + fac_2_ad(2)

             transmission_scatt_ir_ad%azphacup(isti,lay,i) = &
                  transmission_scatt_ir_ad%azphacup(isti,lay,i) + &
                  fac_2_ad(1) * solar_spectrum(i) * z4pi_r

              transmission_scatt_ir_ad%azphacdo(isti,lay,i) = &
                  transmission_scatt_ir_ad%azphacdo(isti,lay,i) + &
                  fac_2_ad(6) * solar_spectrum(i) * z4pi_r
          enddo
       enddo
    endif
  enddo
end subroutine solar_scattering_air_ad

subroutine solar_scattering_near_surf_ad(transmission_aux, transmission_aux_ad, &
                                         ircld, raytracing, raytracing_ad, transmission_scatt_ir_ad, &
                                         solar, solar_spectrum, refl, sateqsun, iv3lay, iv3lev, chanprof, narray, &
                                         pfraction, auxrad_stream, auxrad_stream_ad)

  use rttov_const, only : min_tau, z4pi_r
  use rttov_types, Only : transmission_type_aux, ircld_type, radiance_aux, &
                          rttov_chanprof, transmission_scatt_ir_type, raytracing_type
  use parkind1, only : jpim, jprb, jplm 

  Implicit None

  type(transmission_type_aux), intent(in)         :: transmission_aux
  type(ircld_type),     intent(in)                :: ircld
  type(raytracing_type), intent(in)               :: raytracing
  logical(jplm), intent(in)                       :: solar(:)
  real(jprb), intent(in)                          :: solar_spectrum(:)
  real(jprb), intent(in)                          :: refl(:)
  logical(jplm), intent(in)                       :: sateqsun(:,:)
  integer(jpim), intent(in)                       :: iv3lay(:), iv3lev(:)
  type(rttov_chanprof), intent(in)                :: chanprof(:)
  integer(jpim), intent(in)                       :: narray(:)
  real(jprb), intent(in)                          :: pfraction(:)
  type(radiance_aux), intent(in)                  :: auxrad_stream

  type(transmission_type_aux), intent(inout)      :: transmission_aux_ad
  type(raytracing_type), intent(inout)            :: raytracing_ad
  type(transmission_scatt_ir_type), intent(inout) :: transmission_scatt_ir_ad
  type(radiance_aux), intent(inout)               :: auxrad_stream_ad

  integer(jpim) :: nlevels, nlayers, nstreams, nchannels
  integer(jpim) :: i, ist, lay, lev, lay1

  Real(jprb) :: fac_3_ad(7)
  Real(jprb) :: temp(4,0:narray(2))

! unpack
  nlevels = narray(1)
  nlayers = nlevels - 1_jpim
  nstreams = narray(2)
  nchannels = narray(3)

  Do i = 1, nchannels
     if(solar(i))then
        lay = iv3lay(i)
        lev = iv3lev(i)

        if(pfraction(i) < 0.0_JPRB )then
           lay1 = iv3lay(i)
        else
           lay1 = iv3lay(i) + 1
        endif

       do ist = 0, ircld%nstream(prof) !rev loop not necessary

         if (refl(i) > 0._jprb) then

            if(auxrad_stream%meanrad_down_solar(ist,i) < 0._jprb) then
               auxrad_stream_ad%meanrad_down_solar(ist,i) = 0._jprb
            endif

            auxrad_stream_ad%down_solar(lay,ist,i) = auxrad_stream_ad%down_solar(lay,ist,i) + &
                                                     auxrad_stream_ad%meanrad_down_solar(ist,i)

            fac_3_ad(:) = 0.0_jprb
            temp(1,ist) = auxrad_stream%fac6_3(ist,i) * auxrad_stream%fac2_3(ist,i) * &
                          auxrad_stream_ad%meanrad_down_solar(ist,i)
            temp(2,ist) = transmission_aux%solar_path1%tau_level_r(lev,ist,i) * &
                          transmission_aux%solar_path1%tau_surf_r(ist,i)
            ! By doing the multiplication in the following way we avoid overflows
            temp(3,ist) = (transmission_aux%solar_path1%fac2(lev,ist,i) * temp(2,ist)) * temp(2,ist)

            if(.not. sateqsun(lay,prof)) then
               temp(4,ist) = (auxrad_stream%fac5_3(ist,i) - auxrad_stream%fac4_3(ist,i)) * auxrad_stream%fac7_3(ist,i)

               fac_3_ad(4) = &!fac4_3_ad &
                           - temp(1,ist) * temp(2,ist) * auxrad_stream%fac7_3(ist,i) * &
                             transmission_aux%solar_path2%tau_level(lev,ist,i)

               fac_3_ad(5) = &!fac5_3_ad + &
                           + temp(1,ist) * temp(2,ist) * auxrad_stream%fac7_3(ist,i) * &
                             transmission_aux%solar_path2%tau_level(lev,ist,i)

               fac_3_ad(7) = &!fac7_3_ad + &
                           temp(1,ist) * temp(2,ist) * &
                          (auxrad_stream%fac5_3(ist,i) - auxrad_stream%fac4_3(ist,i)) * &
                           transmission_aux%solar_path2%tau_level(lev,ist,i)

               raytracing_ad%pathsat(lay1,prof) = raytracing_ad%pathsat(lay1,prof) + &
                                                  fac_3_ad(7) * raytracing%pathsun(lay1,prof) / &
                                                 (raytracing%pathsun(lay1,prof) - raytracing%pathsat(lay1,prof))**2

               raytracing_ad%pathsun(lay1,prof) = raytracing_ad%pathsun(lay1,prof) - &
                                                  fac_3_ad(7) * raytracing%pathsat(lay1,prof) / &
                                                 (raytracing%pathsun(lay1,prof) - raytracing%pathsat(lay1,prof))**2
            else
               temp(4,ist) = auxrad_stream%fac4_3(ist,i) * transmission_aux%solar_path2%od_sfrac(ist,i)

               transmission_aux_ad%solar_path2%od_sfrac(ist,i) = transmission_aux_ad%solar_path2%od_sfrac(ist,i) + &
                                                        temp(1,ist) * temp(2,ist) * &
                                                        auxrad_stream%fac4_3(ist,i) * &
                                                        transmission_aux%solar_path2%tau_level(lev,ist,i)

               fac_3_ad(4) = &!fac4_3_ad + &
                           temp(1,ist) * temp(2,ist) * &
                           transmission_aux%solar_path2%tau_level(lev,ist,i) * &
                           transmission_aux%solar_path2%od_sfrac(ist,i)
            endif

            if(transmission_aux%solar_path1%tau_level(lev,ist,i) > min_tau) then

               transmission_aux_ad%solar_path1%tau_surf(ist,i) = transmission_aux_ad%solar_path1%tau_surf(ist,i) - &
                                                     temp(1,ist) * temp(3,ist) * temp(4,ist) * &
                                                     transmission_aux%solar_path2%tau_level(lev,ist,i) * &
                                                     transmission_aux%solar_path1%tau_level(lev,ist,i)

               transmission_aux_ad%solar_path1%tau_level(lev,ist,i) = &
                                                          transmission_aux_ad%solar_path1%tau_level(lev,ist,i) - &
                                                          temp(1,ist) * temp(3,ist) * temp(4,ist) * &
                                                          transmission_aux%solar_path2%tau_level(lev,ist,i) * &
                                                          transmission_aux%solar_path1%tau_surf(ist,i)

               transmission_aux_ad%solar_path2%tau_level(lev,ist,i) = &
                                                          transmission_aux_ad%solar_path2%tau_level(lev,ist,i) + &
                                                          temp(1,ist) * temp(2,ist) * temp(4,ist)

               fac_3_ad(2) = &!fac2_3_ad + &
                           auxrad_stream_ad%meanrad_down_solar(ist,i) * auxrad_stream%fac6_3(ist,i) * &
                           temp(2,ist) * temp(4,ist) * transmission_aux%solar_path2%tau_level(lev,ist,i)

               fac_3_ad(6) = &!fac6_3_ad + &
                           auxrad_stream_ad%meanrad_down_solar(ist,i) * auxrad_stream%fac2_3(ist,i) * &
                           temp(2,ist) * temp(4,ist) * transmission_aux%solar_path2%tau_level(lev,ist,i)
            endif

          else

            fac_3_ad(2) = 0._jprb
            fac_3_ad(4) = 0._jprb
            fac_3_ad(5) = 0._jprb
            fac_3_ad(6) = 0._jprb

          endif ! refl(i) > 0.

       !--------------Upward single scattering of the solar beam---------------------------------

          isti = ircld%icldarr(ist,lay,prof)

          auxrad_stream_ad%up_solar(lay,ist,i) = auxrad_stream_ad%up_solar(lay,ist,i) + &
                                                 auxrad_stream_ad%meanrad_up_solar(ist,i)

          transmission_aux_ad%solar_path2%tau_level(lev,ist,i) = &
                     transmission_aux_ad%solar_path2%tau_level(lev,ist,i) + &
                     auxrad_stream_ad%meanrad_up_solar(ist,i) * &
                     auxrad_stream%fac1_3(ist,i) * auxrad_stream%fac2_3(ist,i) * auxrad_stream%fac3_3(ist,i)

          transmission_aux_ad%solar_path2%tau_surf(ist,i) = &
                     transmission_aux_ad%solar_path2%tau_surf(ist,i) - &
                     auxrad_stream_ad%meanrad_up_solar(ist,i) * &
                     auxrad_stream%fac1_3(ist,i) * auxrad_stream%fac2_3(ist,i) * auxrad_stream%fac3_3(ist,i)

          fac_3_ad(3) = &!fac3_3_ad + &
                      auxrad_stream_ad%meanrad_up_solar(ist,i) * &
                      auxrad_stream%fac1_3(ist,i) * auxrad_stream%fac2_3(ist,i) * &
                     (transmission_aux%solar_path2%tau_level(lev,ist,i) - &
                      transmission_aux%solar_path2%tau_surf(ist,i))

          fac_3_ad(2) = fac_3_ad(2) + &
                      auxrad_stream_ad%meanrad_up_solar(ist,i) * &
                      auxrad_stream%fac1_3(ist,i) * auxrad_stream%fac3_3(ist,i) * &
                     (transmission_aux%solar_path2%tau_level(lev,ist,i) - &
                      transmission_aux%solar_path2%tau_surf(ist,i))

          fac_3_ad(1) = &!fac1_3_ad + &
                      auxrad_stream_ad%meanrad_up_solar(ist,i) * &
                      auxrad_stream%fac2_3(ist,i) * auxrad_stream%fac3_3(ist,i) * &
                     (transmission_aux%solar_path2%tau_level(lev,ist,i) - &
                      transmission_aux%solar_path2%tau_surf(ist,i))

          transmission_scatt_ir_ad%azphacdo(isti,lay,i) = &
               transmission_scatt_ir_ad%azphacdo(isti,lay,i) + &
               fac_3_ad(6) * solar_spectrum(i) * z4pi_r

          transmission_aux_ad%solar_path1%od_sfrac(ist,i) = transmission_aux_ad%solar_path1%od_sfrac(ist,i) - &
               fac_3_ad(5) * auxrad_stream%fac5_3(ist,i)

          transmission_aux_ad%solar_path2%od_sfrac(ist,i) = transmission_aux_ad%solar_path2%od_sfrac(ist,i) - &
               fac_3_ad(4) * auxrad_stream%fac4_3(ist,i)

          raytracing_ad%pathsat(lay,prof) = raytracing_ad%pathsat(lay,prof) + &
               fac_3_ad(3) / raytracing%patheff(lay,prof)

          raytracing_ad%patheff(lay,prof) = raytracing_ad%patheff(lay,prof) - &
               fac_3_ad(3) * raytracing%pathsat(lay,prof) / raytracing%patheff(lay,prof) ** 2

          transmission_scatt_ir_ad%ssa(isti,lay,i) = transmission_scatt_ir_ad%ssa(isti,lay,i) + &
               fac_3_ad(2)

          transmission_scatt_ir_ad%azphacup(isti,lay,i) = &
               transmission_scatt_ir_ad%azphacup(isti,lay,i) + &
               fac_3_ad(1) * solar_spectrum(i) * z4pi_r

       enddo
    endif
  enddo
end subroutine solar_scattering_near_surf_ad

subroutine solar_rayleigh_ad(narray, iv3lev, opts, coef, solar, solar_spectrum, chanprof, &
                             raytracing, raytracing_ad, ircld, profiles, profiles_ad, &
                             profiles_dry, profiles_dry_ad, &
                             transmission_aux, transmission_aux_ad, auxrad_stream_ad)

  use rttov_types, Only : rttov_coef, rttov_chanprof, raytracing_type, profile_type, transmission_type_aux, &
                          radiance_aux, rttov_options, ircld_type

  use rttov_const, Only : z4pi_r, deg2rad, gravity, na, Mh2o, Mair, &
                          ray_scs_wlm, ray_scs_a1, ray_scs_b1, ray_scs_c1, ray_scs_d1, &
                          ray_scs_a2, ray_scs_b2, ray_scs_c2, ray_scs_d2, ray_min_wvn

  use parkind1, only : jprb, jplm, jpim

  implicit none

  integer(jpim), intent(in)                  :: narray(:)
  integer(jpim), intent(in)                  :: iv3lev(:)
  type(rttov_options), intent(in)            :: opts
  type(rttov_coef), intent(in)               :: coef
  logical(jplm), intent(in)                  :: solar(:)
  real(jprb), intent(in)                     :: solar_spectrum(:)
  type(rttov_chanprof), intent(in)           :: chanprof(:)
  type(raytracing_type), intent(in)          :: raytracing
  type(raytracing_type), intent(inout)       :: raytracing_ad
  type(ircld_type), intent(in)               :: ircld
  type(profile_type), intent(in)             :: profiles(:)
  type(profile_type), intent(inout)          :: profiles_ad(:)
  type(profile_type), intent(in)             :: profiles_dry(:)
  type(profile_type), intent(inout)          :: profiles_dry_ad(:)
  type(transmission_type_aux), intent(in)    :: transmission_aux
  type(transmission_type_aux), intent(inout) :: transmission_aux_ad
  type(radiance_aux), intent(inout)          :: auxrad_stream_ad

  real(kind=jprb) :: wlm, ss_param, v_h2o, v_h2o_dry, Mwet
  real(kind=jprb) :: cossat, cossol, cosazi, cosscata_term1, cosscata_term2, cosscata
  real(kind=jprb) :: ray_phase, solar_src, solar_src_updn
  real(kind=jprb) :: v_h2o_ad, v_h2o_dry_ad, Mwet_ad
  real(kind=jprb) :: cosscata_term1_ad, cosscata_term2_ad, cosscata_ad
  real(kind=jprb) :: ray_phase_ad, solar_src_ad, solar_src_updn_ad

  integer(kind=jpim) :: ist, nchannels, nlevels, i, lev, lay
  real(kind=jprb) :: rayrad_up_ad(0:narray(1)-1,0:narray(2)), rayrad_dn_ad(0:narray(1)-1,0:narray(2))

  nchannels = narray(3)
  nlevels = narray(1)

  do i = 1, nchannels

    if(solar(i) .and. coef%ff_cwn(chan) > ray_min_wvn) then

      ! Calculate layer-independent scattering parameter
      wlm = 10000.0_jprb / coef%ff_cwn(chan)    ! Wavelength in microns

      if (wlm < ray_scs_wlm) then
        ss_param = ray_scs_a1 * wlm ** (ray_scs_b1 + ray_scs_c1 * wlm + ray_scs_d1/wlm)
      else
        ss_param = ray_scs_a2 * wlm ** (ray_scs_b2 + ray_scs_c2 * wlm + ray_scs_d2/wlm)
      endif
      ss_param = ss_param * 0.01_jprb ** 2_jpim * na * z4pi_r / gravity

      rayrad_up_ad(:,:) = 0._jprb
      rayrad_dn_ad(:,:) = 0._jprb

      ! Add the contribution from the part-layer above the surface
      if (opts%rt_all%use_q2m) then
        v_h2o_dry = 0.5_jprb * (profiles_dry(prof)%q(iv3lev(i)) + profiles_dry(prof)%s2m%q) * 1.E-06_jprb
      else
        v_h2o_dry = profiles_dry(prof)%q(iv3lev(i)) * 1.E-06_jprb
      endif
      if (profiles(prof)%gas_units /= gas_unit_compatibility) then
        v_h2o = v_h2o_dry / (1._jprb + v_h2o_dry)
      else
        v_h2o = v_h2o_dry
      endif
      Mwet = ((1._jprb  - v_h2o) * Mair + v_h2o * Mh2o) * 1.E-3_jprb

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

      cosscata_term1_ad = 0._jprb
      cosscata_term2_ad = 0._jprb
      solar_src_updn_ad = 0._jprb

      do ist = 0, ircld%nstream(prof)
        solar_src_updn_ad = solar_src_updn_ad + &
                          auxrad_stream_ad%meanrad_down_solar(ist,i) * &
                          transmission_aux%solar_path2%Tau_level(iv3lev(i),ist,i) / &
                          transmission_aux%solar_path1%Tau_level(iv3lev(i),ist,i) ** 3_jpim

        transmission_aux_ad%solar_path2%Tau_level(iv3lev(i),ist,i) = &
                          transmission_aux_ad%solar_path2%Tau_level(iv3lev(i),ist,i) + &
                          auxrad_stream_ad%meanrad_down_solar(ist,i) * &
                          solar_src_updn / transmission_aux%solar_path1%Tau_level(iv3lev(i),ist,i) ** 3_jpim

        transmission_aux_ad%solar_path1%Tau_level(iv3lev(i),ist,i) = &
                          transmission_aux_ad%solar_path1%Tau_level(iv3lev(i),ist,i) - &
                          auxrad_stream_ad%meanrad_down_solar(ist,i) * &
                          solar_src_updn * 3_jpim * transmission_aux%solar_path2%Tau_level(iv3lev(i),ist,i) / &
                          transmission_aux%solar_path1%Tau_level(iv3lev(i),ist,i) ** 4_jpim

        solar_src_updn_ad = solar_src_updn_ad + auxrad_stream_ad%meanrad_up_solar(ist,i) * &
                      transmission_aux%solar_path2%Tau_level(iv3lev(i),ist,i)

        transmission_aux_ad%solar_path2%Tau_level(iv3lev(i),ist,i) = &
                      transmission_aux_ad%solar_path2%Tau_level(iv3lev(i),ist,i) + &
                      solar_src_updn * auxrad_stream_ad%meanrad_up_solar(ist,i)
      enddo

      ray_phase_ad = solar_src * solar_src_updn_ad  ! First use of ray_phase_ad, no accumulation
      solar_src_ad = solar_src_updn_ad * ray_phase  ! First use of solar_src_ad
      cosscata_ad = 2._jprb * 0.75_jprb * ray_phase_ad * cosscata ! First use of cossacata_ad
      cosscata_term1_ad = cosscata_term1_ad - cosscata_ad
      cosscata_term2_ad = cosscata_term2_ad - cosscata_ad

      ! First use of Mwet_ad, no accumulation
      Mwet_ad = solar_src_ad * solar_spectrum(i) * ss_param * &
                raytracing%pathsat(iv2lay(i), prof) * &
                abs(profiles(prof)%s2m%p - profiles(prof)%p(iv3lev(i))) * 100.0_jprb * &
                (-1._jprb) / (Mwet * Mwet)

      raytracing_ad%pathsat(iv2lay(i), prof) = raytracing_ad%pathsat(iv2lay(i), prof) + &
                solar_src_ad * solar_spectrum(i) * ss_param / Mwet * &
                abs(profiles(prof)%s2m%p - profiles(prof)%p(iv3lev(i))) * 100.0_jprb

      ! Note sign of accumulation in the two pressure terms
      if (profiles(prof)%s2m%p >= profiles(prof)%p(iv3lev(i))) then
        profiles_ad(prof)%s2m%p = profiles_ad(prof)%s2m%p + &
                  solar_src_ad * solar_spectrum(i) * ss_param * &
                  raytracing%pathsat(iv2lay(i), prof) * 100.0_jprb / Mwet

        if (opts%interpolation%lgradp) then
          profiles_ad(prof)%p(iv3lev(i)) = profiles_ad(prof)%p(iv3lev(i)) - &
                    solar_src_ad * solar_spectrum(i) * ss_param * &
                    raytracing%pathsat(iv2lay(i), prof) * 100.0_jprb / Mwet
        endif
      else
        profiles_ad(prof)%s2m%p = profiles_ad(prof)%s2m%p - &
                  solar_src_ad * solar_spectrum(i) * ss_param * &
                  raytracing%pathsat(iv2lay(i), prof) * 100.0_jprb / Mwet

        if (opts%interpolation%lgradp) then
          profiles_ad(prof)%p(iv3lev(i)) = profiles_ad(prof)%p(iv3lev(i)) + &
                    solar_src_ad * solar_spectrum(i) * ss_param * &
                    raytracing%pathsat(iv2lay(i), prof) * 100.0_jprb / Mwet
        endif
      endif

      raytracing_ad%zasat(iv2lay(i), prof) = raytracing_ad%zasat(iv2lay(i), prof) + &
                      cosscata_term2_ad * raytracing%zasun(iv2lay(i), prof) * cosazi
      raytracing_ad%zasun(iv2lay(i), prof) = raytracing_ad%zasun(iv2lay(i), prof) + &
                      cosscata_term2_ad * raytracing%zasat(iv2lay(i), prof) * cosazi

      raytracing_ad%zasat(iv2lay(i), prof) = raytracing_ad%zasat(iv2lay(i), prof) - &
                      cosscata_term1_ad * raytracing%zasat(iv2lay(i), prof) * cossol / cosscata_term1
      raytracing_ad%zasun(iv2lay(i), prof) = raytracing_ad%zasun(iv2lay(i), prof) - &
                      cosscata_term1_ad * raytracing%zasun(iv2lay(i), prof) * cossat / cosscata_term1

      ! First use of v_h2o_ad, no accumulation
      v_h2o_ad = Mwet_ad * (Mh2o - Mair) * 1.E-3_jprb

      ! First use of v_h2o_dry_ad, no accumulation
      if (profiles(prof)%gas_units /= gas_unit_compatibility) then
        v_h2o_dry_ad = v_h2o_ad * v_h2o * ( 1._jprb - v_h2o / (1._jprb + v_h2o_dry)) / v_h2o_dry
      else
        v_h2o_dry_ad = v_h2o_ad
      endif

      if (opts%rt_all%use_q2m) then
        profiles_dry_ad(prof)%q(iv3lev(i)) = profiles_dry_ad(prof)%q(iv3lev(i)) + &
                          0.5_jprb * v_h2o_dry_ad * 1.E-06_jprb
        profiles_dry_ad(prof)%s2m%q = profiles_dry_ad(prof)%s2m%q + &
                          0.5_jprb * v_h2o_dry_ad * 1.E-06_jprb
      else
        profiles_dry_ad(prof)%q(iv3lev(i)) = profiles_dry_ad(prof)%q(iv3lev(i)) + &
                          v_h2o_dry_ad * 1.E-06_jprb
      endif

      ! Atmospheric contribution
      do ist = 0, ircld%nstream(prof)
        rayrad_dn_ad(iv3lay(i),ist) = rayrad_dn_ad(iv3lay(i),ist) + auxrad_stream_ad%meanrad_down_solar(ist,i)

        rayrad_dn_ad(1:nlayers,ist) = rayrad_dn_ad(1:nlayers,ist) + auxrad_stream_ad%down_solar(:,ist,i)

        rayrad_up_ad(iv3lay(i),ist) = rayrad_up_ad(iv3lay(i),ist) + auxrad_stream_ad%meanrad_up_solar(ist,i)

        rayrad_up_ad(1:nlayers,ist) = rayrad_up_ad(1:nlayers,ist) + auxrad_stream_ad%up_solar(:,ist,i)
      enddo

      do lev = nlevels, 2, -1
        lay = lev - 1_jpim

        v_h2o_dry = 0.5_jprb * (profiles_dry(prof)%q(lev-1) + profiles_dry(prof)%q(lev)) * 1.E-06_jprb
        if (profiles(prof)%gas_units /= gas_unit_compatibility) then
          v_h2o = v_h2o_dry / (1._jprb + v_h2o_dry)
        else
          v_h2o = v_h2o_dry
        endif
        Mwet = ((1._jprb  - v_h2o) * Mair + v_h2o * Mh2o) * 1.E-3_jprb

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

        cosscata_term1_ad = 0._jprb
        cosscata_term2_ad = 0._jprb
        solar_src_updn_ad = 0._jprb

        do ist = 0, ircld%nstream(prof)
          solar_src_updn_ad = solar_src_updn_ad + rayrad_dn_ad(lay,ist) * &
                            transmission_aux%solar_path2%Tau_level(lev-1,ist,i) / &
                            transmission_aux%solar_path1%Tau_level(lev-1,ist,i) ** 3_jpim

          transmission_aux_ad%solar_path2%Tau_level(lev-1,ist,i) = &
                            transmission_aux_ad%solar_path2%Tau_level(lev-1,ist,i) + &
                            rayrad_dn_ad(lay,ist) * solar_src_updn / &
                            transmission_aux%solar_path1%Tau_level(lev-1,ist,i) ** 3_jpim

          transmission_aux_ad%solar_path1%Tau_level(lev-1,ist,i) = &
                            transmission_aux_ad%solar_path1%Tau_level(lev-1,ist,i) - &
                            rayrad_dn_ad(lay,ist) * solar_src_updn * &
                            3_jpim * transmission_aux%solar_path2%Tau_level(lev-1,ist,i) / &
                            transmission_aux%solar_path1%Tau_level(lev-1,ist,i) ** 4_jpim

          rayrad_dn_ad(lay-1,ist) = rayrad_dn_ad(lay-1,ist) + rayrad_dn_ad(lay,ist)

          solar_src_updn_ad = solar_src_updn_ad + rayrad_up_ad(lay,ist) * &
                        transmission_aux%solar_path2%Tau_level(lev-1,ist,i)

          transmission_aux_ad%solar_path2%Tau_level(lev-1,ist,i) = &
                        transmission_aux_ad%solar_path2%Tau_level(lev-1,ist,i) + &
                        solar_src_updn * rayrad_up_ad(lay,ist)

          rayrad_up_ad(lay-1,ist) = rayrad_up_ad(lay-1,ist) + rayrad_up_ad(lay,ist)
        enddo

        ray_phase_ad = solar_src * solar_src_updn_ad  ! First use of ray_phase_ad, no accumulation
        solar_src_ad = solar_src_updn_ad * ray_phase  ! First use of solar_src_ad
        cosscata_ad = 2._jprb * 0.75_jprb * ray_phase_ad * cosscata ! First use of cossacata_ad
        cosscata_term1_ad = cosscata_term1_ad - cosscata_ad
        cosscata_term2_ad = cosscata_term2_ad - cosscata_ad

        raytracing_ad%pathsat(lay, prof) = raytracing_ad%pathsat(lay, prof) + &
                 solar_src_ad * solar_spectrum(i) * ss_param * &
                 (profiles(prof)%p(lev) - profiles(prof)%p(lev-1)) * 100.0_jprb / Mwet

        ! First use of Mwet_ad, no accumulation
        Mwet_ad = solar_src_ad * solar_spectrum(i) * ss_param * &
                 (profiles(prof)%p(lev) - profiles(prof)%p(lev-1)) * 100.0_jprb * &
                 raytracing%pathsat(lay, prof) * (-1._jprb) / (Mwet * Mwet)

        if (opts%interpolation%lgradp) then
          profiles_ad(prof)%p(lev) = profiles_ad(prof)%p(lev) + &
                          solar_src_ad * solar_spectrum(i) * ss_param * &
                          100.0_jprb * raytracing%pathsat(lay, prof) / Mwet
          profiles_ad(prof)%p(lev-1) = profiles_ad(prof)%p(lev-1) - &
                          solar_src_ad * solar_spectrum(i) * ss_param * &
                          100.0_jprb * raytracing%pathsat(lay, prof) / Mwet
        endif

        raytracing_ad%zasat(lay, prof) = raytracing_ad%zasat(lay, prof) + &
                            cosscata_term2_ad * raytracing%zasun(lay, prof) * cosazi
        raytracing_ad%zasun(lay, prof) = raytracing_ad%zasun(lay, prof) + &
                            cosscata_term2_ad * raytracing%zasat(lay, prof) * cosazi

        raytracing_ad%zasat(lay, prof) = raytracing_ad%zasat(lay, prof) - &
                            cosscata_term1_ad * raytracing%zasat(lay, prof) * cossol / cosscata_term1
        raytracing_ad%zasun(lay, prof) = raytracing_ad%zasun(lay, prof) - &
                            cosscata_term1_ad * raytracing%zasun(lay, prof) * cossat / cosscata_term1

        v_h2o_ad = Mwet_ad * (Mh2o - Mair) * 1.E-3_jprb
        if (profiles(prof)%gas_units /= gas_unit_compatibility) then
          v_h2o_dry_ad = v_h2o_ad * v_h2o * ( 1._jprb - v_h2o / (1._jprb + v_h2o_dry)) / v_h2o_dry
        else
          v_h2o_dry_ad = v_h2o_ad
        endif

        profiles_dry_ad(prof)%q(lev-1) = profiles_dry_ad(prof)%q(lev-1) + &
                          0.5_jprb * v_h2o_dry_ad * 1.E-06_jprb
        profiles_dry_ad(prof)%q(lev) = profiles_dry_ad(prof)%q(lev) + &
                          0.5_jprb * v_h2o_dry_ad * 1.E-06_jprb
      enddo

    endif
  enddo

end subroutine solar_rayleigh_ad

subroutine solar_surface_contribution_ad(transmission_aux, transmission_aux_ad, ircld, &
                                             reflectance, reflectance_ad, refl_norm, chanprof, solar, &
                                             solar_spectrum, narray, auxrad_stream_ad, rad_ad)

  use rttov_types, only : rttov_chanprof, ircld_type, transmission_type_aux, radiance_type, radiance_aux
  use parkind1, only : jpim, jprb, jplm 


  Implicit None

  type(transmission_type_aux), intent(in)       :: transmission_aux
  type(ircld_type),     intent(in)              :: ircld
  real(jprb), intent(in)                        :: reflectance(:)
  real(jprb), intent(in)                        :: refl_norm(:)
  type(rttov_chanprof), intent(in)              :: chanprof(:)
  logical(jplm), intent(in)                     :: solar(:)
  Real(jprb), Intent(in)                        :: solar_spectrum(:)
  integer(jpim), intent(in)                     :: narray(:)

  type(radiance_type), intent(in)               :: rad_ad
  type(radiance_aux), intent(in)                :: auxrad_stream_ad

  type(transmission_type_aux), intent(inout)    :: transmission_aux_ad
  real(jprb), intent(inout)                     :: reflectance_ad(:)

  integer(jpim) :: nlevels, nlayers, nstreams, nchannels, i, ist
  real(jprb)    :: temp_ad(0:narray(2)) ! ist

! unpack
  nlevels = narray(1)
  nlayers = nlevels - 1_jpim
  nchannels = narray(3)  

  Do i = 1, nchannels
     if(solar(i)) then
        nstreams = ircld%nstream(prof)

        temp_ad(0) = rad_ad%clear(i)
        temp_ad(1:nstreams) = auxrad_stream_ad%cloudy(1:nstreams,i)

        do ist = 0, nstreams
          reflectance_ad(i) = reflectance_ad(i) + solar_spectrum(i) * refl_norm(i) * &
                              temp_ad(ist) * transmission_aux%solar_path2%Tau_surf(ist,i)

          transmission_aux_ad%solar_path2%Tau_surf(ist,i) = transmission_aux_ad%solar_path2%Tau_surf(ist,i) + &
                              solar_spectrum(i) * refl_norm(i) * temp_ad(ist) * reflectance(i)
        enddo

     endif
  enddo
end subroutine solar_surface_contribution_ad

End Subroutine rttov_integrate_ad
