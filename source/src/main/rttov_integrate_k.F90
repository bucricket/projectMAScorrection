SUBROUTINE rttov_integrate_k( &
  addcosmic, opts, maxnstreams, chanprof,                      &
  emissivity,                   emissivity_k,                  &
  reflectance,                  reflectance_k,                 &
  refl_norm,                                                   &
  thermrefl,                    thermrefl_k,                   &
  do_lambertian,                                               &
  thermal,                                                     &
  dothermal,                                                   &
  solar,                                                       &
  dosolar,                                                     &
  solar_spectrum,                                              &
  transmission_aux,             transmission_aux_k,            &
  transmission_scatt_ir,        transmission_scatt_ir_k,       &
  profiles,                     profiles_k,                    &
  profiles_dry,                 profiles_dry_k,                &
  aux_prof,                     aux_prof_k,                    &
  coef,                                                        &
  raytracing,                   raytracing_k,                  &
  ircld,                        ircld_k,                       &
  rad,                                                         &
  auxrad,                                                      &
  auxrad_stream,                auxrad_stream_k,               &
                                rad_k)
  
  ! in,                           !inout
  !
  ! Description:
  ! To perform K of integration of radiative transfer equation
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
  ! Matricardi, M. 2005 The inclusion of aerosols and clouds in RTIASI, the ECMWF
  ! fast radiative transfer model for the infrared atmospheric sounding interferometer.
  ! ECMWF Research Dept. Tech. Memo. 474 (available from the librarian at ECMWF)  !
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
  !  1.9     01.07/2006   Marco Matricardi (ECMWF):
  !             --        A parameterization of multiple scattering for aerosols
  !                       and clouds has been included in RTTOV.
  !             --        The "stream method" has been introduced to compute
  !                       radiances in presence of partially cloudy layers
  !             --        The single scattering of the solar beam has been
  !                       introduced.
  !          30/08/2006   New version of this routine created based on AD
  !  1.10    08/01/2007   Corrected some bugs in the 'old' code (R Saunders)
  !  1.11    09/02/2007   Removed polarisation index
  !  1.12    13/03/2007   Reintroduced overcast radiance (R Saunders)
  !  1.13    27/11/2007   Optimised for IBM/NEC (D Salmond)
  !  1.14    17/11/2007   Allowed adjoint/K sensitivity to clear BT/radiance (A Geer)
  !  1.15    21/12/2007   Added polarimetric option (R. Saunders)
  !  1.16    15/09/2009   User defined ToA. Layers distinct from levels. Top-layer
  !                       brought into layer looping to shorten code (P.Rayer)
  !  1.17    03/11/2009   Transmittances on levels (A Geer)
  !  1.18    02/12/2009   Introduced principal component capability. Pathsat, Pathsun and
  !                       related quantities are now layer arrays (Marco Matricardi).
  !  1.19    05/07/2010   Remove addsolar flag from profiles structure (J Hocking)
  !  1.20    14/10/2010   Remove rt8_mode (J Hocking)
  !  1.21    14/12/2010   Use traj0_sta%solar array to flag channels for which solar calculations
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
  !
  
  USE rttov_types, ONLY : &
    rttov_chanprof, rttov_coef, profile_type, profile_aux, transmission_type_aux, transmission_scatt_ir_type, &
    radiance_type, rttov_options, ircld_type, raytracing_type, radiance_aux, rttov_emissivity, &
    rttov_reflectance
  
  USE parkind1, ONLY : jpim, jprb, jplm
  
!INTF_OFF
  USE rttov_const, ONLY : sensor_id_po, pi_r, deg2rad, &
    overcast_albedo_wvn, overcast_albedo1, overcast_albedo2, &
    gas_unit_compatibility
  
  USE yomhook, ONLY : LHOOK, DR_HOOK
!INTF_ON
  
  IMPLICIT NONE
  
  !subroutine arguments:
  LOGICAL(jplm), INTENT(in)                        :: addcosmic
  TYPE(rttov_options), INTENT(in)                  :: opts
  INTEGER(jpim), INTENT(in)                        :: maxnstreams
  TYPE(rttov_chanprof), INTENT(in)                 :: chanprof(:)
  TYPE(profile_type), INTENT(in)                   :: profiles(:)
  TYPE(profile_type), INTENT(in)                   :: profiles_dry(SIZE(profiles))
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
  TYPE(transmission_type_aux), INTENT(in)          :: transmission_aux
  TYPE(transmission_scatt_ir_type), INTENT(in)     :: transmission_scatt_ir
  TYPE(ircld_type), INTENT(in)                     :: ircld
  TYPE(raytracing_type), INTENT(in)                :: raytracing
  TYPE(radiance_type), INTENT(in)                  :: rad
  TYPE(radiance_aux), INTENT(in)                   :: auxrad
  TYPE(radiance_aux), INTENT(in)                   :: auxrad_stream
  
  TYPE(radiance_aux), INTENT(inout)                :: auxrad_stream_k
  TYPE(rttov_emissivity),  INTENT(inout), OPTIONAL :: emissivity_k(SIZE(chanprof))
  TYPE(rttov_reflectance), INTENT(inout), OPTIONAL :: reflectance_k(SIZE(chanprof))
  REAL(jprb), INTENT(inout)                        :: thermrefl_k(SIZE(chanprof))
  TYPE(profile_type), INTENT(inout)                :: profiles_k(SIZE(chanprof))
  TYPE(profile_type), INTENT(inout)                :: profiles_dry_k(SIZE(chanprof))
  TYPE(profile_aux), INTENT(inout)                 :: aux_prof_k
  TYPE(transmission_type_aux), INTENT(inout)       :: transmission_aux_k
  TYPE(transmission_scatt_ir_type), INTENT(inout)  :: transmission_scatt_ir_k
  TYPE(ircld_type), INTENT(inout)                  :: ircld_k
  TYPE(raytracing_type), INTENT(inout)             :: raytracing_k
  TYPE(radiance_type), INTENT(inout)               :: rad_k
  
!INTF_END
  
#include "rttov_calcrad_k.interface"
  
  !local variables:
  REAL(jprb)    :: cfraction(SIZE(chanprof))
  REAL(jprb)    :: cfraction_k(SIZE(chanprof))
  
  REAL(jprb)    :: rad_air_k(profiles(1)%nlevels,SIZE(chanprof))
  REAL(jprb)    :: rad_surfair_k(SIZE(chanprof))
  REAL(jprb)    :: rad_skin_k(SIZE(chanprof))
  
  REAL(jprb)    :: pfraction(SIZE(chanprof))        ! cloud fraction
  LOGICAL(jplm) :: sateqsun(profiles(1)%nlayers,SIZE(profiles(:))) ! True where the solar zenith angle equal to observation angle
  INTEGER(jpim) :: pol_id(SIZE(chanprof))       ! polarisation index
  
  INTEGER(jpim) :: iv2lay(SIZE(chanprof)), iv2lev(SIZE(chanprof)), iv3lay(SIZE(chanprof)), iv3lev(SIZE(chanprof))
  
  INTEGER(jpim) :: i, lev, lay, iprof, ist, isti, nlayers, nlevels, nchannels
  
  REAL(jprb)    :: p, refl, refl_norm_scat
  
  LOGICAL(jplm) :: keyradonly ! flag to indicate only calculate primary radiance outputs
  
  INTEGER(jpim) :: narray(3)
  
  REAL(JPRB) :: ZHOOK_HANDLE
  
  IF (LHOOK) CALL DR_HOOK('RTTOV_INTEGRATE_K',0_jpim,ZHOOK_HANDLE)
  !- End of header --------------------------------------------------------
  
#define prof chanprof(i)%Prof
#define chan chanprof(i)%Chan
  
  !X.  initialisation of local variables
  
  rad_surfair_k(:) = 0._JPRB
  rad_skin_k(:)    = 0._JPRB
  rad_air_k(:,:)   = 0._jprb
  
  ! auxrad_stream_k is initialised in rttov_k
  
  !---------------------------
  !0. Initialise useful variables
  !---------------------------
  nchannels = SIZE(chanprof)
  nlayers = profiles(1)%nlayers
  nlevels = nlayers + 1
  
  DO i = 1, nchannels
    IF(do_lambertian(i)) THEN
      transmission_aux_k%thermal_path1%tau_level_p(:,:,i)   = 0._jprb
      transmission_aux_k%thermal_path1%tau_level_p_r(:,:,i) = 0._jprb
      transmission_aux_k%thermal_path1%tau_surf_p_r(:,i)    = 0._jprb
      transmission_aux_k%thermal_path1%tau_surf_p(:,i)      = 0._jprb
    ENDIF
  ENDDO
  
  ! DAR: narray is a more compact way of passing nchannels, nlevels and maxnstreams to the internal subroutines
  ! Again module variables would be a more elegant way of solving this problem and would results in less duplicated code
  narray(1) = nlevels; narray(2) = maxnstreams; narray(3) = nchannels
  
  keyradonly = opts%rt_ir%addaerosl .OR. opts%rt_ir%pc%addpc
  
  DO i = 1, nchannels
    cfraction(i) = aux_prof%s(prof)%cfraction
    pfraction(i) = aux_prof%s(prof)%pfraction_surf
  ENDDO

  IF (dosolar) THEN
    DO iprof = 1, SIZE(profiles)
      DO lay = 1, nlayers
        IF(raytracing%pathsat(lay,iprof) == raytracing%pathsun(lay,iprof)) THEN
          sateqsun(lay,iprof) = .TRUE.
        ELSE
          sateqsun(lay,iprof) = .FALSE.
        ENDIF
      ENDDO
    ENDDO
  ENDIF

  DO i = 1, nchannels
    ! case-1: surf lies above lev=nlevels
    iv3lev(i) = aux_prof%s(prof)%nearestlev_surf - 1   ! lowest lev above surf
    ! case-2: surf lies below lev=nlevels
    IF (pfraction(i) < 0.0_JPRB) iv3lev(i) = iv3lev(i) + 1  ! iv3lev=iv2lev=lowest lev above surf
    
    iv2lev(i) = aux_prof%s(prof)%nearestlev_surf       ! highest lev below surf
    iv2lay(i) = iv2lev(i) - 1                          ! same layer as that numbered by  iv2 in RTTOV-9
    iv3lay(i) = iv3lev(i) - 1                          ! same layer as that numbered by  iv3 in RTTOV-9
  ENDDO
  
  IF (coef%id_sensor == sensor_id_po) THEN
    DO i = 1, nchannels
      pol_id(i) = coef%fastem_polar(chan) + 1_jpim
    ENDDO
  ELSE
    pol_id(:) = 0_jpim
  ENDIF
  
  
  !---------------------------
  ! 7. Calculate total radiance
  !---------------------------
  cfraction_k(1:nchannels) = 0._JPRB ! DAR: This line not in integrate_ad
  
  IF (opts%rt_ir%addaerosl .OR. opts%rt_ir%addclouds) THEN
    rad_k%cloudy(:) = rad_k%total(:)
    !----------------------------------------
    ! Calculate complex cloudy radiances
    !----------------------------------------
    DO i = 1, nchannels
      rad_k%clear(i)     = rad_k%cloudy(i) * ircld%XSTRCLR(prof)
      ircld_k%XSTRCLR(i) = ircld_k%XSTRCLR(i) + rad_k%cloudy(i) * rad%clear(i)
      
      DO ist = ircld%nstream(prof), 1, -1 !reverse loop nec. due to data dependence
        auxrad_stream_k%cloudy(ist,i) = auxrad_stream_k%cloudy(ist,i) + &
          rad_k%cloudy(i) * (ircld%xstr(ist+1,prof) - ircld%xstr(ist,prof))
        
        ircld_k%xstr(ist+1,i)         = ircld_k%xstr(ist+1,i) + &
          rad_k%cloudy(i) * auxrad_stream%cloudy(ist,i)
        
        ircld_k%xstr(ist,i)           = ircld_k%xstr(ist,i) - &
          rad_k%cloudy(i) * auxrad_stream%cloudy(ist,i)
      ENDDO
    ENDDO
  ELSE
    !---------------------------
    ! Calculate total radiance (clear case/simple cloud)
    !---------------------------
    IF (opts%rt_ir%pc%addpc) THEN
      rad_k%clear(:)  = rad_k%total(:)
    ELSE
      rad_k%clear(1:nchannels)  = rad_k%clear(1:nchannels) + (1.0_jprb - cfraction(:)) * rad_k%total(1:nchannels)
      rad_k%cloudy(1:nchannels) = cfraction(:) * rad_k%total(1:nchannels)
      cfraction_k(:)  = (rad%cloudy(1:nchannels) - rad%clear(1:nchannels)) * rad_k%total(1:nchannels)
      
      ! Interpolate to given cloud-top pressures
      DO i = 1, nchannels
        lay = aux_prof%s(prof)%nearestlev_ctp - 1
        p   = aux_prof%s(prof)%pfraction_ctp
        
        rad_k%overcast(lay,i)         = & !rad_k%overcast(lay,i) + & !first use of rad_k%overcast
          (1._jprb - p) * rad_k%cloudy(i)
        
        rad_k%overcast(lay-1,i)       = & !rad_k%overcast(lay-1,i) + & !first use of rad_k%overcast
          p * rad_k%cloudy(i)
        
        aux_prof_k%s(i)%pfraction_ctp = aux_prof_k%s(i)%pfraction_ctp + &
          (rad%overcast(lay-1,i) - rad%overcast(lay,i)) * rad_k%cloudy(i)
      ENDDO
      
    ENDIF
    ! rad_k%cloudy done
  ENDIF
  
  !---------------------------------------------------
  ! 6. Calculate overcast radiances
  !---------------------------------------------------
  IF (.NOT. keyradonly) THEN
    ist = 0_jpim
    
    !cdir nodep
    DO i = 1, nchannels
      IF (thermal(i)) THEN
        DO lay = 1, nlayers
          ! overcast radiances at given cloud top
          
          ! Make exception for iv2lay because up(iv2lay) was replaced by meanrad_up in direct and
          ! overcast(iv2lay) was replaced by surface contribution
          IF (lay == iv2lay(i)) THEN
            auxrad_stream_k%meanrad_up(ist,i) = auxrad_stream_k%meanrad_up(ist,i) + rad_k%overcast(lay,i)
            rad_surfair_k(i) = rad_surfair_k(i) + rad_k%overcast(lay,i) * &
              transmission_aux%thermal_path1%tau_surf(ist,i)
            transmission_aux_k%thermal_path1%tau_surf(ist,i) = &
              transmission_aux_k%thermal_path1%tau_surf(ist,i) + &
              rad_k%overcast(lay,i) * auxrad%surfair(i)
          ELSE
            auxrad_stream_k%up(lay,ist,i) = auxrad_stream_k%up(lay,ist,i) + rad_k%overcast(lay,i)
            rad_air_k(lay+1,i) = &!rad_air_k(lay+1,i) + &
              rad_k%overcast(lay,i) * transmission_aux%thermal_path1%tau_level(lay+1,ist,i)
            transmission_aux_k%thermal_path1%tau_level(lay+1,ist,i) = &
              transmission_aux_k%thermal_path1%tau_level(lay+1,ist,i) + &
              auxrad%air(lay+1,i) * rad_k%overcast(lay,i)
          ENDIF
        ENDDO
      ELSE IF (solar(i)) THEN
        IF (reflectance(i)%refl_cloud_top > 0) THEN
          refl = reflectance(i)%refl_cloud_top
        ELSE
          IF (coef%ff_cwn(chan) > overcast_albedo_wvn) THEN
            refl = overcast_albedo1 * pi_r
          ELSE
            refl = overcast_albedo2 * pi_r
          ENDIF
        ENDIF
        DO lay = 1, nlayers
          lev = lay + 1
          ! overcast radiances at given cloud top
          transmission_aux_k%solar_path2%Tau_level(lev,ist,i) = &
            transmission_aux_k%solar_path2%Tau_level(lev,ist,i) + &
            solar_spectrum(i) * refl / raytracing%pathsun(lay, prof) * &
            rad_k%overcast(lay,i)
          raytracing_k%pathsun(lay, i) = raytracing_k%pathsun(lay, i) - &
            solar_spectrum(i) * refl / raytracing%pathsun(lay, prof)**2_jpim * &
            transmission_aux%solar_path2%Tau_level(lev,ist,i) * rad_k%overcast(lay,i)
          
          transmission_aux_k%solar_path1%Tau_level(lev,ist,i) = &
            transmission_aux_k%solar_path1%Tau_level(lev,ist,i) + &
            refl / raytracing%pathsat(lay, prof) * &
            auxrad_stream%down_solar(lay,ist,i) * &
            transmission_aux%solar_path1%Tau_level(lev,ist,i) * &
            2_jpim * rad_k%overcast(lay,i)
          raytracing_k%pathsat(lay, i) = &
            raytracing_k%pathsat(lay, i) - &
            refl / raytracing%pathsat(lay, prof) * &
            transmission_aux%solar_path1%Tau_level(lev,ist,i)**2_jpim * &
            auxrad_stream%down_solar(lay,ist,i) / raytracing%pathsat(lay, prof) * &
            rad_k%overcast(lay,i)
          auxrad_stream_k%down_solar(lay,ist,i) = &
            auxrad_stream_k%down_solar(lay,ist,i) + &
            refl / raytracing%pathsat(lay, prof) * &
            transmission_aux%solar_path1%Tau_level(lev,ist,i)**2_jpim * &
            rad_k%overcast(lay,i)
          
          ! Make exception for iv2lay because up_solar(iv2lay) was replaced by meanrad_up_solar in direct
          IF (lay == iv2lay(i)) THEN
            auxrad_stream_k%meanrad_up_solar(ist,i) = &
              auxrad_stream_k%meanrad_up_solar(ist,i) + rad_k%overcast(lay,i)
          ELSE
            auxrad_stream_k%up_solar(lay,ist,i) = &
              auxrad_stream_k%up_solar(lay,ist,i) + rad_k%overcast(lay,i)
          ENDIF
        ENDDO
      ELSE
        rad_k%overcast(:,i) = 0._jprb
      ENDIF
    ENDDO
    ! rad_k%overcast done
    ! rad_k%clear done
    
    ! N.b. this is where rttov_integrate_k differs from rttov_integrate_ad -
    ! AD code is per profile and k code is over nchannels (hence i)
    DO i = 1, nchannels
      aux_prof_k%s(i)%cfraction = aux_prof_k%s(i)%cfraction + cfraction_k(i)
    ENDDO
  ENDIF
  
  !--------------------------------
  !5. cosmic temperature correction
  !--------------------------------
  
  !calculate planck function corresponding to tcosmic=2.7k
  !deblonde tcosmic for microwave sensors only
  
  IF (addcosmic) THEN
    DO i = 1, nchannels
      ist = 0_jpim
      
      IF (do_lambertian(i)) THEN
        thermrefl_k(i) = thermrefl_k(i) + &
          auxrad%cosmic(i) * &
          transmission_aux%thermal_path1%tau_surf_p(ist,i) * &
          transmission_aux%thermal_path1%tau_surf(ist,i) * &
          rad_k%clear(i)
        
        transmission_aux_k%thermal_path1%tau_surf(ist,i) = &! transmission_aux_k%thermal_path1%tau_surf(ist,i) + &
          auxrad%cosmic(i) * thermrefl(i) * &
          transmission_aux%thermal_path1%tau_surf_p(ist,i) * &
          rad_k%clear(i)
        
        transmission_aux_k%thermal_path1%tau_surf_p(ist,i) = &
          transmission_aux_k%thermal_path1%tau_surf_p(ist,i) + &
          auxrad%cosmic(i) * thermrefl(i) * &
          transmission_aux%thermal_path1%tau_surf(ist,i) * &
          rad_k%clear(i)
      ELSE
        
        thermrefl_k(i) = thermrefl_k(i) + &
          rad_k%clear(i) * auxrad%cosmic(i) * &
          transmission_aux%thermal_path1%tau_surf(ist,i)**2_jpim
        
        transmission_aux_k%thermal_path1%tau_surf(ist,i) = & !transmission_aux_k%thermal_path1%tau_surf(ist,i) +
          rad_k%clear(i) * 2.0_jprb * thermrefl(i) * auxrad%cosmic(i) * &
          transmission_aux%thermal_path1%tau_surf(ist,i)
      ENDIF
    ENDDO
  ENDIF
  
  !----------------------------------------------------------------------------------------
  ! 4. solar surface contribution
  !----------------------------------------------------------------------------------------
  IF (dosolar) &
    CALL solar_surface_contribution_k(transmission_aux, transmission_aux_k, &
    ircld, reflectance%refl_out, reflectance_k%refl_out, refl_norm, &
    chanprof, solar, solar_spectrum, &
    narray, auxrad_stream_k, rad_k)
  
  !-----------------------
  !3. surface emission contribution
  !-----------------------
  IF (dothermal) THEN
    DO i = 1, nchannels
      IF (thermal(i)) THEN
        ist = 0_jpim
        
        emissivity_k(i)%emis_out = emissivity_k(i)%emis_out + &
          rad_k%clear(i) * auxrad%skin(i) * transmission_aux%thermal_path1%tau_surf(ist,i)
        
        transmission_aux_k%thermal_path1%tau_surf(ist,i) = transmission_aux_k%thermal_path1%tau_surf(ist,i) + &
          rad_k%clear(i) * auxrad%skin(i) * emissivity(i)%emis_out
        
        rad_skin_k(i) = rad_skin_k(i) + &
          rad_k%clear(i) * emissivity(i)%emis_out * transmission_aux%thermal_path1%tau_surf(ist,i)
        
        DO ist = 1, ircld%nstream(prof) !rev loop not necessary
          emissivity_k(i)%emis_out = emissivity_k(i)%emis_out + &
            auxrad_stream_k%cloudy(ist,i) * auxrad%skin(i) * &
            transmission_aux%thermal_path1%tau_surf(ist,i)
          
          transmission_aux_k%thermal_path1%tau_surf(ist,i) = transmission_aux_k%thermal_path1%tau_surf(ist,i) + &
            auxrad_stream_k%cloudy(ist,i) * auxrad%skin(i) * emissivity(i)%emis_out
          
          rad_skin_k(i) = rad_skin_k(i) + &
            auxrad_stream_k%cloudy(ist,i) * emissivity(i)%emis_out * &
            transmission_aux%thermal_path1%tau_surf(ist,i)
        ENDDO
      ENDIF
    ENDDO
  ENDIF
  
  
  !-------------------------------------
  !2. calculate atmospheric contribution
  !-------------------------------------
  !cdir nodep
  DO i = 1, nchannels
    ist = 0_jpim
    
    IF (thermal(i)) THEN
      auxrad_stream_k%meanrad_up(ist,i) = auxrad_stream_k%meanrad_up(ist,i) + rad_k%clear(i)
      
      IF (do_lambertian(i)) THEN
        auxrad_stream_k%meanrad_down(ist,i) = &!auxrad_stream_k%meanrad_down(ist,i) + &
          thermrefl(i) * transmission_aux%thermal_path1%tau_surf_p(ist,i) * &
          transmission_aux%thermal_path1%tau_surf(ist,i) * rad_k%clear(i)
        
        thermrefl_k(i) = thermrefl_k(i) + &
          auxrad_stream%meanrad_down(ist,i) * &
          transmission_aux%thermal_path1%tau_surf_p(ist,i) * &
          transmission_aux%thermal_path1%tau_surf(ist,i) * rad_k%clear(i)
        
        transmission_aux_k%thermal_path1%tau_surf_p(ist,i) = &
          transmission_aux_k%thermal_path1%tau_surf_p(ist,i) + &
          auxrad_stream%meanrad_down(ist,i) * thermrefl(i) * &
          transmission_aux%thermal_path1%tau_surf(ist,i) * rad_k%clear(i)
        
        transmission_aux_k%thermal_path1%tau_surf(ist,i) = &
          transmission_aux_k%thermal_path1%tau_surf(ist,i) + &
          auxrad_stream%meanrad_down(ist,i) * thermrefl(i) * &
          transmission_aux%thermal_path1%tau_surf_p(ist,i) * rad_k%clear(i)
      ELSE
        
        auxrad_stream_k%meanrad_down(ist,i) = & ! &auxrad_stream_k%meanrad_down(ist,i) + &
          rad_k%clear(i) * thermrefl(i) * transmission_aux%thermal_path1%tau_surf(ist,i)**2
        
        transmission_aux_k%thermal_path1%tau_surf(ist,i) = transmission_aux_k%thermal_path1%tau_surf(ist,i) + &
          rad_k%clear(i) * 2._JPRB * thermrefl(i) * &
          transmission_aux%thermal_path1%tau_surf(ist,i) * &
          auxrad_stream%meanrad_down(ist,i)
        
        thermrefl_k(i) = thermrefl_k(i) + &
          rad_k%clear(i) * transmission_aux%thermal_path1%tau_surf(ist,i)**2 * &
          auxrad_stream%meanrad_down(ist,i)
      ENDIF
      
    ENDIF
    
    IF (solar(i)) THEN
      auxrad_stream_k%meanrad_up_solar(ist,i) = auxrad_stream_k%meanrad_up_solar(ist,i) + rad_k%clear(i)
      
      refl_norm_scat = COS(profiles(prof)%zenangle * deg2rad)
      
      reflectance_k(i)%refl_out = reflectance_k(i)%refl_out + &
        transmission_aux%solar_path1%Tau_surf(ist,i)**2_jpim * refl_norm_scat * &
        auxrad_stream%meanrad_down_solar(ist,i) * rad_k%clear(i)
      
      auxrad_stream_k%meanrad_down_solar(ist,i) = &! meanrad_down_solar_k(ist,i) + &
        transmission_aux%solar_path1%Tau_surf(ist,i)**2_jpim * &
        reflectance(i)%refl_out * refl_norm_scat * rad_k%clear(i)
      
      transmission_aux_k%solar_path1%Tau_surf(ist,i) = transmission_aux_k%solar_path1%Tau_surf(ist,i) + &
        transmission_aux%solar_path1%Tau_surf(ist,i) * &
        reflectance(i)%refl_out * refl_norm_scat * 2._JPRB * rad_k%clear(i) * &
        auxrad_stream%meanrad_down_solar(ist,i)
    ENDIF
    
    IF (thermal(i)) THEN
      DO ist = 1, ircld%nstream(prof)
        ! add upward and downward parts
        auxrad_stream_k%meanrad_up(ist,i) = auxrad_stream_k%meanrad_up(ist,i) + &
          auxrad_stream_k%cloudy(ist,i)
        
        IF (do_lambertian(i)) THEN
          auxrad_stream_k%meanrad_down(ist,i) = &!auxrad_stream_k%meanrad_down(ist,i) + &
            thermrefl(i) * transmission_aux%thermal_path1%tau_surf_p(ist,i) * &
            transmission_aux%thermal_path1%tau_surf(ist,i) * auxrad_stream_k%cloudy(ist,i)

          thermrefl_k(i) = thermrefl_k(i) + &
            auxrad_stream%meanrad_down(ist,i) * &
            transmission_aux%thermal_path1%tau_surf_p(ist,i) * &
            transmission_aux%thermal_path1%tau_surf(ist,i) * auxrad_stream_k%cloudy(ist,i)

          transmission_aux_k%thermal_path1%tau_surf_p(ist,i) = &
            transmission_aux_k%thermal_path1%tau_surf_p(ist,i) + &
            auxrad_stream%meanrad_down(ist,i) * thermrefl(i) * &
            transmission_aux%thermal_path1%tau_surf(ist,i) * auxrad_stream_k%cloudy(ist,i)

          transmission_aux_k%thermal_path1%tau_surf(ist,i) = &
            transmission_aux_k%thermal_path1%tau_surf(ist,i) + &
            auxrad_stream%meanrad_down(ist,i) * thermrefl(i) * &
            transmission_aux%thermal_path1%tau_surf_p(ist,i) * auxrad_stream_k%cloudy(ist,i)
        ELSE
          auxrad_stream_k%meanrad_down(ist,i) = & ! auxrad_stream_k%meanrad_down(ist,i) + &
            auxrad_stream_k%cloudy(ist,i) * transmission_aux%thermal_path1%tau_surf(ist,i)**2 * &
            thermrefl(i)

          transmission_aux_k%thermal_path1%tau_surf(ist,i) = transmission_aux_k%thermal_path1%tau_surf(ist,i) + &
            auxrad_stream_k%cloudy(ist,i) * &
            2._JPRB * transmission_aux%thermal_path1%tau_surf(ist,i) * thermrefl(i) * &
            auxrad_stream%meanrad_down(ist,i)

          thermrefl_k(i) = thermrefl_k(i) + &
            auxrad_stream_k%cloudy(ist,i) * &
            transmission_aux%thermal_path1%tau_surf(ist,i)**2 * &
            auxrad_stream%meanrad_down(ist,i)
        ENDIF
      ENDDO
    ENDIF
      
    IF (solar(i)) THEN
      DO ist = 1, ircld%nstream(prof)

        auxrad_stream_k%meanrad_up_solar(ist,i) = auxrad_stream_k%meanrad_up_solar(ist,i) + &
          auxrad_stream_k%cloudy(ist,i)
        
        reflectance_k(i)%refl_out = reflectance_k(i)%refl_out + &
          transmission_aux%solar_path1%Tau_surf(ist,i)**2_jpim * refl_norm_scat * &
          auxrad_stream%meanrad_down_solar(ist,i) * auxrad_stream_k%cloudy(ist,i)
        
        auxrad_stream_k%meanrad_down_solar(ist,i) = &! meanrad_down_solar_k(ist,i) + &
          transmission_aux%solar_path1%Tau_surf(ist,i)**2_jpim * &
          reflectance(i)%refl_out * refl_norm_scat * auxrad_stream_k%cloudy(ist,i)
        
        transmission_aux_k%solar_path1%Tau_surf(ist,i) = transmission_aux_k%solar_path1%Tau_surf(ist,i) + &
          transmission_aux%solar_path1%Tau_surf(ist,i) * &
          reflectance(i)%refl_out * refl_norm_scat * 2._JPRB * auxrad_stream_k%cloudy(ist,i) * &
          auxrad_stream%meanrad_down_solar(ist,i)
      ENDDO
    ENDIF
  ENDDO
  
  !-------------------------------------------------------------------------------
  !2b calculate clear-sky Rayleigh scattering contribution
  !-------------------------------------------------------------------------------
  IF (dosolar) &
    CALL solar_rayleigh_k(narray, iv3lev, opts, coef, solar, solar_spectrum, chanprof, &
    raytracing, raytracing_k, ircld, profiles, profiles_k, &
    profiles_dry, profiles_dry_k, &
    transmission_aux, transmission_aux_k, auxrad_stream_k)
  
  !-------------------------------------------------------------------------------
  !2a calculate near-surface layer contribution
  !-------------------------------------------------------------------------------
  IF((opts%rt_ir%addaerosl .OR. opts%rt_ir%addclouds) .AND. dosolar) &
    CALL solar_scattering_near_surf_k(transmission_aux, transmission_aux_k, ircld, raytracing, raytracing_k, &
    transmission_scatt_ir_k, solar, solar_spectrum, reflectance%refl_out, &
    sateqsun, iv3lay, iv3lev, chanprof, narray, pfraction, auxrad_stream, &
    auxrad_stream_k)
  
  IF (dothermal) CALL calc_near_surf_contribution_k(thermal, transmission_aux, transmission_aux_k, auxrad, ircld, &
    chanprof, narray, iv3lay, iv3lev, pol_id, rad_surfair_k, rad_air_k, auxrad_stream_k)
  
  !-------------------------------------
  !2. calculate atmospheric contribution
  !-------------------------------------
  IF((opts%rt_ir%addaerosl .OR. opts%rt_ir%addclouds) .AND. dosolar) &
    CALL solar_scattering_air_k(transmission_aux, transmission_aux_k, ircld, raytracing, raytracing_k, &
    transmission_scatt_ir_k, solar, solar_spectrum, reflectance%refl_out, &
    sateqsun, chanprof, narray, auxrad_stream, auxrad_stream_k)
  
  IF (dothermal) CALL calc_atmospheric_radiance_k(thermal, transmission_aux, transmission_aux_k, &
    auxrad, ircld, chanprof, narray, rad_air_k, auxrad_stream_k)
  
  rad_k%clear(:) = 0._jprb
  
  !----------------------------
  !1. calculate layer radiances
  !----------------------------
  
  IF (dothermal) CALL rttov_calcrad_k(chanprof, profiles, profiles_k, coef, auxrad, &
    rad_skin_k, rad_surfair_k, rad_air_k)
  
  IF (LHOOK) CALL DR_HOOK('RTTOV_INTEGRATE_K',1_jpim,ZHOOK_HANDLE)
  
CONTAINS
  
  
  SUBROUTINE calc_atmospheric_radiance_k(thermal, transmission_aux, transmission_aux_k, &
    auxrad, ircld, chanprof, narray, rad_air_k, auxrad_stream_k)
    
    USE rttov_types, ONLY : transmission_type_aux, radiance_aux, ircld_type, rttov_chanprof
    USE rttov_const, ONLY : min_tau
    USE parkind1, ONLY : jpim, jprb, jplm
    
    IMPLICIT NONE
    
    LOGICAL(jplm), INTENT(in)                     :: thermal(:)
    TYPE(transmission_type_aux), INTENT(in)       :: transmission_aux
    TYPE(radiance_aux), INTENT(in)                :: auxrad
    TYPE(ircld_type),     INTENT(in)              :: ircld
    TYPE(rttov_chanprof), INTENT(in)              :: chanprof(:)
    INTEGER(jpim), INTENT(in)                     :: narray(:)
    
    REAL(jprb), INTENT(inout)                     :: rad_air_k(:,:)
    TYPE(transmission_type_aux), INTENT(inout)    :: transmission_aux_k
    TYPE(radiance_aux), INTENT(inout)             :: auxrad_stream_k
    
    INTEGER(jpim) :: nlevels, nlayers, nchannels
    INTEGER(jpim) :: i, ist, lay, lev, levm1
    
    REAL(jprb) :: tau_prod
    REAL(jprb) :: b1_k, b2_k, b3_k
    REAL(jprb) :: temp
    REAL(jprb) :: z(3), ztemp(6)
    
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
          
          DO lay = nlayers, 2, -1
            IF(auxrad_stream%down_ref(lay,ist,i) < 0._jprb) &
              auxrad_stream_k%down(lay,ist,i) = 0._jprb
            
            auxrad_stream_k%up(lay-1,ist,i) = auxrad_stream_k%up(lay-1,ist,i) + auxrad_stream_k%up(lay,ist,i)
            auxrad_stream_k%down(lay-1,ist,i) = auxrad_stream_k%down(lay-1,ist,i) + auxrad_stream_k%down(lay,ist,i)
          ENDDO
          
          lay = 1
          IF(auxrad_stream%down_ref(lay,ist,i) < 0._jprb) THEN
            auxrad_stream_k%down(lay,ist,i) = 0._jprb
          ENDIF
          
          DO lay = 1, nlayers
            levm1 = lay
            lev = lay+1
            isti = ircld%icldarr(ist,lay,prof)
            IF(transmission_aux%anynegtau > 0._jprb) THEN
              
              IF(lay > 1 .AND. tau_levm1 >= 0.AND. tau_lev < 0) THEN
                
                IF (tau_lev > min_tau) THEN
                  tau_prod = tau_lev_r * tau_levm1_r
                  auxrad_stream_k%up(lay,ist,i) = auxrad_stream_k%up(lay,ist,i) + &
                    auxrad_stream_k%down(lay,ist,i) * tau_prod
                  
                  transmission_aux_k%thermal_path1%tau_level(lev,ist,i) = &
                    transmission_aux_k%thermal_path1%tau_level(lev,ist,i) - &
                    auxrad_stream_k%down(lay,ist,i) * &
                    tau_levm1 * (tau_prod*tau_prod) * &
                    0.5_JPRB * (auxrad%air(lev,i) + auxrad%air(levm1,i)) * (tau_levm1 - tau_lev )
                  
                  transmission_aux_k%thermal_path1%tau_level(levm1,ist,i) = &
                    transmission_aux_k%thermal_path1%tau_level(levm1,ist,i) - &
                    auxrad_stream_k%down(lay,ist,i) * &
                    tau_lev * (tau_prod*tau_prod) * &
                    0.5_JPRB * (auxrad%air(lev,i) + auxrad%air(levm1,i)) * (tau_levm1 - tau_lev)
                END IF
                
                rad_air_k(lev,i) = rad_air_k(lev,i) + &
                  auxrad_stream_k%up(lay,ist,i) * 0.5_JPRB * (tau_levm1 - tau_lev)
                
                rad_air_k(levm1,i) = rad_air_k(levm1,i) + &
                  auxrad_stream_k%up(lay,ist,i) * 0.5_JPRB * (tau_levm1 - tau_lev)
                
                transmission_aux_k%thermal_path1%tau_level(levm1,ist,i) = &
                  transmission_aux_k%thermal_path1%tau_level(levm1,ist,i) + &
                  auxrad_stream_k%up(lay,ist,i) * 0.5_JPRB * (auxrad%air(lev,i) + auxrad%air(levm1,i))
                
                transmission_aux_k%thermal_path1%tau_level(lev,ist,i) = &
                  transmission_aux_k%thermal_path1%tau_level(lev,ist,i) - &
                  auxrad_stream_k%up(lay,ist,i) * 0.5_JPRB * (auxrad%air(lev,i) + auxrad%air(levm1,i))
              ELSE
                temp = transmission_aux%fac1(lay,ist,i) * &
                  transmission_aux%thermal_path1%fac2(lay+1,ist,i) * auxrad_stream_k%down(lay,ist,i)
                
                B1_K = temp * (tau_levm1_r * tau_lev_r)
                B2_K = temp * (tau_lev_r)**2
                B3_K = -temp * (tau_levm1_r * tau_lev_r)
                
                B1_K = B1_K + auxrad_stream_k%up(lay,ist,i)
                B2_K = B2_K - auxrad_stream_k%up(lay,ist,i)
                B3_K = B3_K + auxrad_stream_k%up(lay,ist,i)
                
                transmission_aux_k%thermal_path1%tau_level(levm1,ist,i) = &
                  transmission_aux_k%thermal_path1%tau_level(levm1,ist,i) + &
                  transmission_aux%fac1(lay,ist,i) * (&
                  B1_K * auxrad%air(levm1,i) - &
                  temp * &
                  (B1_2 - &
                  B3_2) * &
                  (tau_lev_r * tau_levm1_r**2) + &
                  B3_K * daux_lay * transmission_aux%thermal_path1%od_singlelayer_r(isti,lay,i))
                
                transmission_aux_k%thermal_path1%tau_level(lev,ist,i) = &
                  transmission_aux_k%thermal_path1%tau_level(lev,ist,i) + &
                  transmission_aux%fac1(lay,ist,i) * (&
                  -temp * &
                  (B1_2 - &
                  B3_2) * &
                  (tau_lev_r**2 * tau_levm1_r) - &
                  B3_K * daux_lay * transmission_aux%thermal_path1%od_singlelayer_r(isti,lay,i) - &
                  temp * 2.0_jprb * B2_2 * (tau_lev_r**3) - &
                  B1_K * auxrad%air(levm1,i) + &
                  B2_K * daux_lay)
                
                transmission_aux_k%thermal_path1%od_singlelayer(isti,lay,i) = &
                  transmission_aux_k%thermal_path1%od_singlelayer(isti,lay,i) - &
                  transmission_aux%fac1(lay,ist,i) * B3_K * &
                  B3_2 * &
                  transmission_aux%thermal_path1%od_singlelayer_r(isti,lay,i)
                
                rad_air_k(levm1,i) = rad_air_k(levm1,i) + transmission_aux%fac1(lay,ist,i) * (&
                  -B3_K * dtau_lay * transmission_aux%thermal_path1%od_singlelayer_r(isti,lay,i) + &
                  B1_K * dtau_lay - &
                  B2_K * tau_lev)
                
                rad_air_k(lev,i) = rad_air_k(lev,i) + transmission_aux%fac1(lay,ist,i) * (&
                  B3_K * dtau_lay * transmission_aux%thermal_path1%od_singlelayer_r(isti,lay,i) + &
                  B2_K * tau_lev)
              ENDIF
            ELSE
              IF (transmission_aux%fac1(lay,ist,i) > 0._jprb) THEN
                IF (do_lambertian(i)) THEN
                  temp = 0._jprb
                  B1_K = 0._jprb
                  B2_K = 0._jprb
                  B3_K = 0._jprb
                  
                  IF(transmission_aux%fac1(lay,ist,i) * transmission_aux%thermal_path1%fac2(lay+1,ist,i) > 0._jprb) THEN
                    z(1) = 0.5_jprb * &
                      (transmission_aux%thermal_path1%tau_level_p(lay,ist,i) - &
                      transmission_aux%thermal_path1%tau_level_p(lev,ist,i)) * &
                      transmission_aux%thermal_path1%tau_level_p_r(lev,ist,i) * &
                      transmission_aux%thermal_path1%tau_level_p_r(lay,ist,i) * &
                      auxrad_stream_k%down(lay,ist,i)
                    
                    rad_air_k(lay,i)  = rad_air_k(lay,i) + z(1)
                    rad_air_k(lay+1,i) = rad_air_k(lay+1,i) + z(1)
                    
                    z(2) = 0.5_jprb * (auxrad%air(lay,i) + auxrad%air(lay+1,i)) * &
                      (transmission_aux%thermal_path1%tau_level_p_r(lev,ist,i) * &
                      transmission_aux%thermal_path1%tau_level_p_r(lay,ist,i)) * &
                      auxrad_stream_k%down(lay,ist,i)
                    
                    transmission_aux_k%thermal_path1%tau_level_p(lay,ist,i) = &
                      transmission_aux_k%thermal_path1%tau_level_p(lay,ist,i) + z(2)
                    
                    transmission_aux_k%thermal_path1%tau_level_p(lev,ist,i) = &
                      transmission_aux_k%thermal_path1%tau_level_p(lev,ist,i) - z(2)
                    
                    z(3) = 0.5_jprb * (auxrad%air(lay,i) + auxrad%air(lay+1,i)) * &
                      (transmission_aux%thermal_path1%tau_level_p(lay,ist,i) - &
                      transmission_aux%thermal_path1%tau_level_p(lev,ist,i)) * &
                      auxrad_stream_k%down(lay,ist,i)
                    
                    transmission_aux_k%thermal_path1%tau_level_p_r(lev,ist,i) = &
                      transmission_aux_k%thermal_path1%tau_level_p_r(lev,ist,i) + &
                      z(3) * transmission_aux%thermal_path1%tau_level_p_r(lay,ist,i)
                    
                    transmission_aux_k%thermal_path1%tau_level_p_r(lay,ist,i) = &
                      transmission_aux_k%thermal_path1%tau_level_p_r(lay,ist,i) + &
                      z(3) * transmission_aux%thermal_path1%tau_level_p_r(lev,ist,i)
                  ENDIF
                ELSE
                  temp = transmission_aux%thermal_path1%fac2(lay+1,ist,i) * auxrad_stream_k%down(lay,ist,i) * &
                    (tau_lev_r * tau_levm1_r)
                ENDIF

                B1_K = temp + auxrad_stream_k%up(lay,ist,i)
                B2_K = temp * tau_lev_r * tau_levm1 - auxrad_stream_k%up(lay,ist,i)
                B3_K = -temp + auxrad_stream_k%up(lay,ist,i)
                
                ztemp(1) = B3_K * transmission_aux%thermal_path1%od_singlelayer_r(isti,lay,i)
                ztemp(2) = ztemp(1) * daux_lay
                ztemp(3) = ztemp(1) * dtau_lay
                ztemp(4) = temp * &
                  (B1_2 - &
                   B3_2)
                ztemp(5) = B1_K * auxrad%air(levm1,i)
                ztemp(6) = B2_K * tau_lev
                
                 transmission_aux_k%thermal_path1%tau_level(levm1,ist,i) = &
                   transmission_aux_k%thermal_path1%tau_level(levm1,ist,i) + &
                   ztemp(5) - &
                   ztemp(4) * tau_levm1_r + &
                   ztemp(2)
                               
                 transmission_aux_k%thermal_path1%tau_level(lev,ist,i) = &
                   transmission_aux_k%thermal_path1%tau_level(lev,ist,i) - &
                   ztemp(4) * tau_lev_r - &
                   ztemp(2) - &
                   temp * 2.0_jprb * B2_2 * (tau_lev_r**2 * tau_levm1) - &
                   ztemp(5) + &
                   B2_K * daux_lay
                
                 transmission_aux_k%thermal_path1%od_singlelayer(isti,lay,i) = &
                   transmission_aux_k%thermal_path1%od_singlelayer(isti,lay,i) - &
                   B3_2 * ztemp(1)

                rad_air_k(levm1,i) = rad_air_k(levm1,i) - &
                  ztemp(3) + &
                  B1_K * dtau_lay - &
                  ztemp(6)
                
                rad_air_k(lev,i) = rad_air_k(lev,i) + &
                  ztemp(3) + &
                  ztemp(6)
              ENDIF
            ENDIF
          ENDDO
        END DO
      ENDIF
    END DO
  END SUBROUTINE calc_atmospheric_radiance_k
  
  SUBROUTINE calc_near_surf_contribution_k(thermal, transmission_aux, transmission_aux_k, auxrad, ircld, chanprof, &
    narray, iv3lay, iv3lev, pol_id, rad_surfair_k, rad_air_k, auxrad_stream_k)
    
    USE rttov_const, ONLY : min_tau, min_od
    USE rttov_types, ONLY : rttov_chanprof, ircld_type, transmission_type_aux, radiance_aux
    USE parkind1, ONLY : jpim, jprb, jplm
    
    IMPLICIT NONE
    
    LOGICAL(jplm), INTENT(in)                  :: thermal(:)
    TYPE(transmission_type_aux), INTENT(in)    :: transmission_aux
    TYPE(radiance_aux), INTENT(in)             :: auxrad
    TYPE(ircld_type), INTENT(in)               :: ircld
    TYPE(rttov_chanprof), INTENT(in)           :: chanprof(:)
    INTEGER(jpim), INTENT(in)                  :: narray(:)
    INTEGER(jpim), INTENT(in)                  :: iv3lay(:)
    INTEGER(jpim), INTENT(in)                  :: iv3lev(:)
    INTEGER(jpim), INTENT(in)                  :: pol_id(:)
    
    TYPE(transmission_type_aux), INTENT(inout) :: transmission_aux_k
    REAL(jprb), INTENT(inout)                  :: rad_surfair_k(:), rad_air_k(:,:)
    TYPE(radiance_aux), INTENT(inout)          :: auxrad_stream_k
    
    !local
    INTEGER(jpim) :: nlayers, nlevels, nstreams, nchannels
    
    REAL(jprb)    :: b1_3, b2_3, b3_3
    REAL(jprb)    :: b1_k, b2_k, b3_k
    
    INTEGER(jpim) :: i, lev, lay, ist
    REAL(jprb) :: z(3)
    
    ! unpack
    nlevels = narray(1)
    nlayers = nlevels - 1_jpim
    nstreams = narray(2)
    nchannels = narray(3)
    
    !cdir nodep
    DO i = 1, nchannels
      IF (thermal(i)) THEN
        lay = iv3lay(i)
        lev = iv3lev(i)
        
        DO ist = 0, ircld%nstream(prof) !rev loop not necessary
          ! assume there is no atmospheric source term for 3rd/4th stokes vector elements
          IF (pol_id(i) >= 6_jpim) auxrad_stream_k%meanrad_up(ist,i) = 0.0_jprb
          
          auxrad_stream_k%up(lay,ist,i) = auxrad_stream_k%up(lay,ist,i) + &
            auxrad_stream_k%meanrad_up(ist,i)
          
          auxrad_stream_k%down(lay,ist,i) = auxrad_stream_k%down(lay,ist,i) + &
            auxrad_stream_k%meanrad_down(ist,i)
          
          IF(transmission_aux%thermal_path1%tau_surf(ist,i) < 0._JPRB) THEN
            IF(transmission_aux%thermal_path1%tau_level(lev,ist,i) >= 0._JPRB) THEN ! nearly always true?
              
              rad_surfair_k(i) = rad_surfair_k(i) + &
                0.5_JPRB * auxrad_stream_k%meanrad_up(ist,i) * &
                (transmission_aux%thermal_path1%tau_level(lev,ist,i) - &
                transmission_aux%thermal_path1%tau_surf(ist,i))
              
              rad_air_k(lev,i) = rad_air_k(lev,i) + &
                0.5_JPRB * auxrad_stream_k%meanrad_up(ist,i)* &
                (transmission_aux%thermal_path1%tau_level(lev,ist,i) - &
                transmission_aux%thermal_path1%tau_surf(ist,i))
              
              transmission_aux_k%thermal_path1%tau_level(lev,ist,i) = &
                transmission_aux_k%thermal_path1%tau_level(lev,ist,i) + &
                0.5_JPRB * auxrad_stream_k%meanrad_up(ist,i) * &
                (auxrad%surfair(i) + auxrad%air(lev,i))
              
              transmission_aux_k%thermal_path1%tau_surf(ist,i) = &
                transmission_aux_k%thermal_path1%tau_surf(ist,i) - &
                0.5_JPRB * auxrad_stream_k%meanrad_up(ist,i) * &
                (auxrad%surfair(i) + auxrad%air(lev,i))
            ENDIF
          ELSE
            IF(transmission_aux%thermal_path1%od_sfrac(ist,i) < min_od .OR. &
              ((transmission_aux%thermal_path1%tau_level(lev,ist,i) - &
              transmission_aux%thermal_path1%tau_surf(ist,i)) < min_od)) THEN
            ELSE
              IF (transmission_aux%thermal_path1%tau_surf(ist,i) > min_tau) THEN
                IF (do_lambertian(i)) THEN
                  z(1) = 0.5_jprb * (transmission_aux%thermal_path1%tau_level_p(lev,ist,i) - &
                    transmission_aux%thermal_path1%tau_surf_p(ist,i)) * &
                    transmission_aux%thermal_path1%tau_level_p_r(lev,ist,i) * &
                    transmission_aux%thermal_path1%tau_surf_p_r(ist,i) * &
                    auxrad_stream_k%meanrad_down(ist,i)
                  
                  rad_surfair_k(i) = rad_surfair_k(i) + z(1)
                  rad_air_k(lev,i) = rad_air_k(lev,i) + z(1)
                  
                  z(2) = 0.5_jprb * (auxrad%surfair(i) + auxrad%air(lev,i)) * &
                    (transmission_aux%thermal_path1%tau_level_p_r(lev,ist,i) * &
                    transmission_aux%thermal_path1%tau_surf_p_r(ist,i)) * &
                    auxrad_stream_k%meanrad_down(ist,i)
                  
                  transmission_aux_k%thermal_path1%tau_surf_p(ist,i) = &
                    transmission_aux_k%thermal_path1%tau_surf_p(ist,i) - z(2)
                  
                  transmission_aux_k%thermal_path1%tau_level_p(lev,ist,i) = &
                    transmission_aux_k%thermal_path1%tau_level_p(lev,ist,i) + &
                    z(2)
                  
                  z(3) = 0.5_jprb * (auxrad%surfair(i) + auxrad%air(lev,i)) * &
                    (transmission_aux%thermal_path1%tau_level_p(lev,ist,i) - &
                    transmission_aux%thermal_path1%tau_surf_p(ist,i)) * &
                    auxrad_stream_k%meanrad_down(ist,i)
                  
                  transmission_aux_k%thermal_path1%tau_level_p_r(lev,ist,i) = &
                    transmission_aux_k%thermal_path1%tau_level_p_r(lev,ist,i) + &
                    z(3) * transmission_aux%thermal_path1%tau_surf_p_r(ist,i)
                  
                  transmission_aux_k%thermal_path1%tau_surf_p_r(ist,i) = &
                    transmission_aux_k%thermal_path1%tau_surf_p_r(ist,i) + &
                    z(3) * transmission_aux%thermal_path1%tau_level_p_r(lev,ist,i)
                  
                  B1_K = 0.0_JPRB
                  B2_K = 0.0_JPRB
                  B3_K = 0.0_JPRB
                ELSE
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
                  B1_K = auxrad_stream_k%meanrad_down(ist,i) * &
                    (transmission_aux%thermal_path1%tau_level_r(lev,ist,i) * &
                    transmission_aux%thermal_path1%tau_surf_r(ist,i))
                  B3_K = -auxrad_stream_k%meanrad_down(ist,i) * &
                    (transmission_aux%thermal_path1%tau_level_r(lev,ist,i) * &
                    transmission_aux%thermal_path1%tau_surf_r(ist,i))
                  B2_K = auxrad_stream_k%meanrad_down(ist,i) * &
                    (transmission_aux%thermal_path1%tau_surf_r(ist,i))**2
                  
                  transmission_aux_k%thermal_path1%tau_surf(ist,i) = &
                    transmission_aux_k%thermal_path1%tau_surf(ist,i) - &
                    auxrad_stream_k%meanrad_down(ist,i) * 2._jprb * B2_3 * &
                    transmission_aux%thermal_path1%tau_surf_r(ist,i)**3_jpim
                  
                  transmission_aux_k%thermal_path1%tau_surf(ist,i) = &
                    transmission_aux_k%thermal_path1%tau_surf(ist,i) - &
                    auxrad_stream_k%meanrad_down(ist,i) * (B1_3 - B3_3) * &
                    transmission_aux%thermal_path1%tau_level_r(lev,ist,i) * &
                    transmission_aux%thermal_path1%tau_surf_r(ist,i)**2_jpim
                  
                  transmission_aux_k%thermal_path1%tau_level(lev,ist,i) = &
                    transmission_aux_k%thermal_path1%tau_level(lev,ist,i) - &
                    auxrad_stream_k%meanrad_down(ist,i) * (B1_3 - B3_3) * &
                    transmission_aux%thermal_path1%tau_level_r(lev,ist,i)**2_jpim * &
                    transmission_aux%thermal_path1%tau_surf_r(ist,i)
                ENDIF
              ELSE
                B1_K = 0.0_JPRB
                B2_K = 0.0_JPRB
                B3_K = 0.0_JPRB
              ENDIF
              
              B1_K = B1_K + auxrad_stream_k%meanrad_up(ist,i)
              B2_K = B2_K - auxrad_stream_k%meanrad_up(ist,i)
              B3_K = B3_K + auxrad_stream_k%meanrad_up(ist,i)
              
              transmission_aux_k%thermal_path1%od_sfrac(ist,i) = &
                transmission_aux_k%thermal_path1%od_sfrac(ist,i) - &
                B3_K * (auxrad%surfair(i)-auxrad%air(lev,i)) * &
                (transmission_aux%thermal_path1%tau_level(lev,ist,i) - &
                transmission_aux%thermal_path1%tau_surf(ist,i)) * &
                (transmission_aux%thermal_path1%od_sfrac_r(ist,i)**2)
              
              transmission_aux_k%thermal_path1%tau_surf(ist,i) = &
                transmission_aux_k%thermal_path1%tau_surf(ist,i) - &
                B3_K * (auxrad%surfair(i)-auxrad%air(lev,i)) * &
                transmission_aux%thermal_path1%od_sfrac_r(ist,i) + &
                B2_K * (auxrad%surfair(i)-auxrad%air(lev,i)) - &
                B1_K * auxrad%air(lev,i)
              
              transmission_aux_k%thermal_path1%tau_level(lev,ist,i) = &
                transmission_aux_k%thermal_path1%tau_level(lev,ist,i) + &
                B3_K * (auxrad%surfair(i) - auxrad%air(lev,i)) * &
                transmission_aux%thermal_path1%od_sfrac_r(ist,i) + &
                B1_K * auxrad%air(lev,i)
              
              rad_air_k(lev,i) = rad_air_k(lev,i) - &
                B3_K * (transmission_aux%thermal_path1%tau_level(lev,ist,i) - &
                transmission_aux%thermal_path1%tau_surf(ist,i)) * &
                transmission_aux%thermal_path1%od_sfrac_r(ist,i) - &
                B2_K * transmission_aux%thermal_path1%tau_surf(ist,i) + &
                B1_K * (transmission_aux%thermal_path1%tau_level(lev,ist,i) - &
                transmission_aux%thermal_path1%tau_surf(ist,i))
              
              rad_surfair_k(i) = rad_surfair_k(i) + &
                B3_K * (transmission_aux%thermal_path1%tau_level(lev,ist,i) - &
                transmission_aux%thermal_path1%tau_surf(ist,i)) * &
                transmission_aux%thermal_path1%od_sfrac_r(ist,i) + &
                B2_K * transmission_aux%thermal_path1%tau_surf(ist,i)
            ENDIF
          ENDIF
        ENDDO
      ENDIF
    END DO
  END SUBROUTINE calc_near_surf_contribution_k
  
  SUBROUTINE solar_scattering_air_k(transmission_aux, transmission_aux_k, ircld, raytracing, raytracing_k, &
    transmission_scatt_ir_k, solar, solar_spectrum, refl, sateqsun, &
    chanprof, narray, auxrad_stream, auxrad_stream_k)
    
    USE rttov_const, ONLY : min_tau, z4pi_r
    USE rttov_types, ONLY : transmission_type_aux, ircld_type, radiance_aux, &
      rttov_chanprof, transmission_scatt_ir_type, raytracing_type
    USE parkind1, ONLY : jpim, jprb, jplm
    
    IMPLICIT NONE
    
    TYPE(transmission_type_aux), INTENT(in)         :: transmission_aux
    TYPE(ircld_type),     INTENT(in)                :: ircld
    TYPE(raytracing_type), INTENT(in)               :: raytracing
    LOGICAL(jplm), INTENT(in)                       :: solar(:)
    REAL(jprb), INTENT(in)                          :: solar_spectrum(:)
    REAL(jprb), INTENT(in)                          :: refl(:)
    LOGICAL(jplm), INTENT(in)                       :: sateqsun(:,:)
    TYPE(rttov_chanprof), INTENT(in)                :: chanprof(:)
    INTEGER(jpim), INTENT(in)                       :: narray(:)
    TYPE(radiance_aux), INTENT(in)                  :: auxrad_stream
    
    TYPE(transmission_type_aux), INTENT(inout)      :: transmission_aux_k
    TYPE(raytracing_type), INTENT(inout)            :: raytracing_k
    TYPE(transmission_scatt_ir_type), INTENT(inout) :: transmission_scatt_ir_k
    TYPE(radiance_aux), INTENT(inout)               :: auxrad_stream_k
    
    INTEGER(jpim) :: nlevels, nlayers, nstreams, nchannels
    INTEGER(jpim) :: i, ist, lay, lev, levm1
    
    REAL(jprb) :: fac_2_k(7)
    REAL(jprb) :: temp(3,narray(1)-1)
    
    ! unpack
    nlevels = narray(1)
    nlayers = nlevels - 1_jpim
    nstreams = narray(2)
    nchannels = narray(3)
    
    !cdir nodep
    DO i = 1, nchannels
      IF(solar(i)) THEN
        
        DO ist = 0, ircld%nstream(prof)!rev loop not necessary
          
          IF (refl(i) > 0._jprb) THEN
            
            DO lay = nlayers, 2, -1
              IF(auxrad_stream%down_ref_solar(lay,ist,i) < 0._jprb) &
                auxrad_stream_k%down_solar(lay,ist,i) = 0._jprb
              
              auxrad_stream_k%down_solar(lay-1,ist,i) = auxrad_stream_k%down_solar(lay-1,ist,i) + &
                auxrad_stream_k%down_solar(lay,ist,i) * transmission_aux%solar_path1%fac2(lay+1,ist,i)
            ENDDO
            
            
            temp(1,1:nlayers) = auxrad_stream%fac6_2(1:nlayers,ist,i) * auxrad_stream%fac2_2(1:nlayers,ist,i) * &
              auxrad_stream_k%down_solar(1:nlayers,ist,i)
            temp(2,1:nlayers) = transmission_aux%solar_path1%tau_level_r(1:nlayers,ist,i) * &
              transmission_aux%solar_path1%tau_level_r(2:nlayers+1,ist,i)
          ENDIF
          
          DO lay = nlayers, 2, -1
            auxrad_stream_k%up_solar(lay-1,ist,i) = auxrad_stream_k%up_solar(lay-1,ist,i) + &
              auxrad_stream_k%up_solar(lay,ist,i)
          ENDDO
          
          DO lay = 1, nlayers !rev loop not necessary
            
            levm1 = lay
            lev   = lay + 1
            isti = ircld%icldarr(ist,lay,prof)

            fac_2_k(:) = 0.0_jprb ! DAR: if errors then may need to zero everything. Not just 2 and 6
            
            IF (refl(i) > 0._jprb) THEN
              
              IF(lay == 1 .OR. transmission_aux%solar_path1%Tau_level(lay+1,ist,i) > min_tau) THEN
                !-------------------Downward single scattering of the solar beam--------------------------
                IF(.NOT. sateqsun(lay,prof)) THEN
                  temp(3,lay) = (auxrad_stream%fac5_2(lay,ist,i) - auxrad_stream%fac4_2(lay,ist,i)) * &
                    auxrad_stream%fac7_2(lay,ist,i)
                  
                  IF(lay > 1_jpim) THEN
                    transmission_aux_k%solar_path1%tau_level(levm1,ist,i) = &
                      transmission_aux_k%solar_path1%tau_level(levm1,ist,i) - &
                      temp(1,lay) * temp(2,lay)**2 * temp(3,lay) * tausun_levm1 * &
                      transmission_aux%solar_path1%Tau_level(lay+1,ist,i)
                    
                    transmission_aux_k%solar_path2%tau_level(levm1,ist,i) = &
                      transmission_aux_k%solar_path2%tau_level(levm1,ist,i) + &
                      temp(1,lay) * temp(2,lay) * temp(3,lay)
                  ENDIF
                  
                  fac_2_k(4) = &!fac4_2_k
                    - temp(1,lay) * temp(2,lay) * &
                    auxrad_stream%fac7_2(lay,ist,i) * transmission_aux%solar_path2%tau_level(lay,ist,i)
                  
                  fac_2_k(5) = &!fac5_2_k + &
                    temp(1,lay) * temp(2,lay) * &
                    auxrad_stream%fac7_2(lay,ist,i) * transmission_aux%solar_path2%tau_level(lay,ist,i)
                  
                  fac_2_k(7) = &!fac7_2_k + &
                    temp(1,lay) * temp(2,lay) * &
                    (auxrad_stream%fac5_2(lay,ist,i) - auxrad_stream%fac4_2(lay,ist,i)) * &
                    transmission_aux%solar_path2%tau_level(lay,ist,i)
                  
                  raytracing_k%pathsat(lay,i) = raytracing_k%pathsat(lay,i) + &
                    fac_2_k(7) * raytracing%pathsun(lay,prof) / &
                    (raytracing%pathsun(lay,prof) - raytracing%pathsat(lay,prof))**2
                  
                  raytracing_k%pathsun(lay,i) = raytracing_k%pathsun(lay,i) - &
                    fac_2_k(7) * raytracing%pathsat(lay,prof) / &
                    (raytracing%pathsun(lay,prof) - raytracing%pathsat(lay,prof))**2
                  
                ELSE
                  temp(3,lay) = auxrad_stream%fac4_2(lay,ist,i) * &
                    transmission_aux%solar_path2%od_singlelayer(isti,lay,i)
                  
                  IF(lay > 1) THEN
                    !                       transmission_aux_k%solar_path1%tau_level(levm1,ist,i) = &
                    !                            transmission_aux_k%solar_path1%tau_level(levm1,ist,i) - &
                    !                            temp(1,lay) * temp(2,lay) * temp(3,lay)
                    !DAR 3108 think this should be temp22 * temp3 * tau_lev as above from integrate_tl
                    transmission_aux_k%solar_path1%tau_level(levm1,ist,i) = &
                      transmission_aux_k%solar_path1%tau_level(levm1,ist,i) - &
                      temp(1,lay) * temp(2,lay)**2_jpim * temp(3,lay) * &
                      transmission_aux%solar_path1%Tau_level(lay+1,ist,i)
                    
                    transmission_aux_k%solar_path2%tau_level(levm1,ist,i) = &
                      transmission_aux_k%solar_path2%tau_level(levm1,ist,i) + &
                      temp(1,lay) * temp(2,lay) * temp(3,lay)
                  ENDIF
                  
                  transmission_aux_k%solar_path2%od_singlelayer(isti,lay,i) = &
                    transmission_aux_k%solar_path2%od_singlelayer(isti,lay,i) + &
                    temp(1,lay) * temp(2,lay) * auxrad_stream%fac4_2(lay,ist,i) * tausun_levm1
                  
                  fac_2_k(4) = &!fac4_2_k + &
                    temp(1,lay) * temp(2,lay) * &
                    transmission_aux%solar_path2%tau_level(lay,ist,i) * &
                    transmission_aux%solar_path2%od_singlelayer(isti,lay,i)
                ENDIF
                
                transmission_aux_k%solar_path1%tau_level(lev,ist,i) = &
                  transmission_aux_k%solar_path1%tau_level(lev,ist,i) - &
                  temp(1,lay) * temp(3,lay) * temp(2,lay)**2_jpim * tausun_levm1 * &
                  transmission_aux%solar_path1%Tau_level(lay,ist,i)
                
                fac_2_k(2) = &!fac2_2_k + &
                  auxrad_stream_k%down_solar(lay,ist,i) * auxrad_stream%fac6_2(lay,ist,i) * &
                  temp(2,lay) * temp(3,lay) * tausun_levm1
                
                fac_2_k(6) = &!fac6_2_k + &
                  auxrad_stream_k%down_solar(lay,ist,i) * auxrad_stream%fac2_2(lay,ist,i) * &
                  temp(2,lay) * temp(3,lay) * tausun_levm1
              ELSE
                auxrad_stream_k%down_solar(lay,ist,i) = 0._JPRB
              ENDIF  ! min_tau
              
            ENDIF ! refl(i) > 0.
            
            !----------------Upward single scattering of the solar beam-------------------------------

            isti = ircld%icldarr(ist,lay,prof)

            fac_2_k(1) = &!fac1_2_k + &
              auxrad_stream_k%up_solar(lay,ist,i) * &
              auxrad_stream%fac2_2(lay,ist,i) * auxrad_stream%fac3_2(lay,ist,i) * &
              (tausun_levm1 - tausun_lev)
            
            fac_2_k(2) = fac_2_k(2) + &
              auxrad_stream_k%up_solar(lay,ist,i) * &
              auxrad_stream%fac1_2(lay,ist,i) * auxrad_stream%fac3_2(lay,ist,i) * &
              (tausun_levm1 - tausun_lev)
            
            fac_2_k(3) = &!fac3_2_k + &
              auxrad_stream_k%up_solar(lay,ist,i) * &
              auxrad_stream%fac1_2(lay,ist,i) * auxrad_stream%fac2_2(lay,ist,i) * &
              (tausun_levm1 - tausun_lev)
            
            IF(lay > 1) THEN
              transmission_aux_k%solar_path2%tau_level(levm1,ist,i) = &
                transmission_aux_k%solar_path2%tau_level(levm1,ist,i) + &
                auxrad_stream_k%up_solar(lay,ist,i) * &
                auxrad_stream%fac1_2(lay,ist,i) * auxrad_stream%fac2_2(lay,ist,i) * &
                auxrad_stream%fac3_2(lay,ist,i)
            ENDIF
            
            transmission_aux_k%solar_path2%tau_level(lev,ist,i) = &
              transmission_aux_k%solar_path2%tau_level(lev,ist,i) - &
              auxrad_stream_k%up_solar(lay,ist,i) * &
              auxrad_stream%fac1_2(lay,ist,i) * auxrad_stream%fac2_2(lay,ist,i) * &
              auxrad_stream%fac3_2(lay,ist,i)
            
            transmission_aux_k%solar_path1%od_singlelayer(isti,lay,i) = &
              transmission_aux_k%solar_path1%od_singlelayer(isti,lay,i) - &
              fac_2_k(5) * auxrad_stream%fac5_2(lay,ist,i) !exp(od_singlelayer)
            
            transmission_aux_k%solar_path2%od_singlelayer(isti,lay,i) = &
              transmission_aux_k%solar_path2%od_singlelayer(isti,lay,i) - &
              fac_2_k(4) * auxrad_stream%fac4_2(lay,ist,i) !exp(od_sunsinglelayer)

            raytracing_k%pathsat(lay,i) = raytracing_k%pathsat(lay,i) + &
                 fac_2_k(3) / raytracing%patheff(lay,prof)

            raytracing_k%patheff(lay,i) = raytracing_k%patheff(lay,i) - &
                 fac_2_k(3) * raytracing%pathsat(lay,prof) / raytracing%patheff(lay,prof) ** 2

            transmission_scatt_ir_k%ssa(isti,lay,i) = &
              transmission_scatt_ir_k%ssa(isti,lay,i) + fac_2_k(2)
            
            transmission_scatt_ir_k%azphacup(isti,lay,i) = &
              transmission_scatt_ir_k%azphacup(isti,lay,i) + &
              fac_2_k(1) * solar_spectrum(i) * z4pi_r
            
            transmission_scatt_ir_k%azphacdo(isti,lay,i) = &
              transmission_scatt_ir_k%azphacdo(isti,lay,i) + &
              fac_2_k(6) * solar_spectrum(i) * z4pi_r
          ENDDO
        ENDDO
      ENDIF
    ENDDO
  END SUBROUTINE solar_scattering_air_k
  
  SUBROUTINE solar_scattering_near_surf_k(transmission_aux, transmission_aux_k, ircld, raytracing, raytracing_k, &
    transmission_scatt_ir_k, solar, solar_spectrum, refl, sateqsun, &
    iv3lay, iv3lev, chanprof, narray, pfraction, auxrad_stream, auxrad_stream_k)
    
    USE rttov_const, ONLY : min_tau, z4pi_r
    USE rttov_types, ONLY : transmission_type_aux, ircld_type, radiance_aux, &
      rttov_chanprof, transmission_scatt_ir_type, raytracing_type
    USE parkind1, ONLY : jpim, jprb, jplm
    
    IMPLICIT NONE
    
    TYPE(transmission_type_aux), INTENT(in)         :: transmission_aux
    TYPE(ircld_type),     INTENT(in)                :: ircld
    TYPE(raytracing_type), INTENT(in)               :: raytracing
    LOGICAL(jplm), INTENT(in)                       :: solar(:)
    REAL(jprb), INTENT(in)                          :: solar_spectrum(:)
    REAL(jprb), INTENT(in)                          :: refl(:)
    LOGICAL(jplm), INTENT(in)                       :: sateqsun(:,:)
    INTEGER(jpim), INTENT(in)                       :: iv3lay(:), iv3lev(:)
    TYPE(rttov_chanprof), INTENT(in)                :: chanprof(:)
    INTEGER(jpim), INTENT(in)                       :: narray(:)
    REAL(jprb), INTENT(in)                          :: pfraction(:)
    TYPE(radiance_aux), INTENT(in)                  :: auxrad_stream
    
    TYPE(transmission_type_aux), INTENT(inout)      :: transmission_aux_k
    TYPE(raytracing_type), INTENT(inout)            :: raytracing_k
    TYPE(transmission_scatt_ir_type), INTENT(inout) :: transmission_scatt_ir_k
    TYPE(radiance_aux), INTENT(inout)               :: auxrad_stream_k
    
    INTEGER(jpim) :: nlevels, nlayers, nstreams, nchannels
    INTEGER(jpim) :: i, ist, lay, lev, lay1
    
    REAL(jprb) :: fac_3_k(7)
    REAL(jprb) :: temp(4,0:narray(2))
    
    ! unpack
    nlevels = narray(1)
    nlayers = nlevels - 1_jpim
    nstreams = narray(2)
    nchannels = narray(3)
    
    DO i = 1, nchannels
      IF(solar(i))THEN
        lay = iv3lay(i)
        lev = iv3lev(i)
        
        IF(pfraction(i) < 0.0_JPRB )THEN
          lay1 = iv3lay(i)
        ELSE
          lay1 = iv3lay(i) + 1
        ENDIF
        
        DO ist = 0, ircld%nstream(prof) !rev loop not necessary
          
          IF (refl(i) > 0._jprb) THEN
            
            IF(auxrad_stream%meanrad_down_solar(ist,i) < 0._jprb) THEN
              auxrad_stream_k%meanrad_down_solar(ist,i) = 0._jprb
            ENDIF
            
            auxrad_stream_k%down_solar(lay,ist,i) = auxrad_stream_k%down_solar(lay,ist,i) + &
              auxrad_stream_k%meanrad_down_solar(ist,i)
            
            fac_3_k(:) = 0.0_jprb
            temp(1,ist) = auxrad_stream%fac6_3(ist,i) * auxrad_stream%fac2_3(ist,i) * &
              auxrad_stream_k%meanrad_down_solar(ist,i)
            temp(2,ist) = transmission_aux%solar_path1%tau_level_r(lev,ist,i) * &
              transmission_aux%solar_path1%tau_surf_r(ist,i)
            ! By doing the multiplication in the following way we avoid overflows
            temp(3,ist) = (transmission_aux%solar_path1%fac2(lev,ist,i) * temp(2,ist)) * temp(2,ist)
            
            IF(.NOT. sateqsun(lay,prof)) THEN
              temp(4,ist) = (auxrad_stream%fac5_3(ist,i) - auxrad_stream%fac4_3(ist,i)) * auxrad_stream%fac7_3(ist,i)
              
              fac_3_k(4) = &!fac4_3_k &
                - temp(1,ist) * temp(2,ist) * auxrad_stream%fac7_3(ist,i) * &
                transmission_aux%solar_path2%tau_level(lev,ist,i)
              
              fac_3_k(5) = &!fac5_3_k + &
                + temp(1,ist) * temp(2,ist) * auxrad_stream%fac7_3(ist,i) * &
                transmission_aux%solar_path2%tau_level(lev,ist,i)
              
              fac_3_k(7) = &!fac7_3_k + &
                temp(1,ist) * temp(2,ist) * &
                (auxrad_stream%fac5_3(ist,i) - auxrad_stream%fac4_3(ist,i)) * &
                transmission_aux%solar_path2%tau_level(lev,ist,i)
              
              raytracing_k%pathsat(lay1,i) = raytracing_k%pathsat(lay1,i) + &
                fac_3_k(7) * raytracing%pathsun(lay1,prof) / &
                (raytracing%pathsun(lay1,prof) - raytracing%pathsat(lay1,prof))**2
              
              raytracing_k%pathsun(lay1,i) = raytracing_k%pathsun(lay1,i) - &
                fac_3_k(7) * raytracing%pathsat(lay1,prof) / &
                (raytracing%pathsun(lay1,prof) - raytracing%pathsat(lay1,prof))**2
            ELSE
              temp(4,ist) = auxrad_stream%fac4_3(ist,i) * transmission_aux%solar_path2%od_sfrac(ist,i)
              
              transmission_aux_k%solar_path2%od_sfrac(ist,i) = transmission_aux_k%solar_path2%od_sfrac(ist,i) + &
                temp(1,ist) * temp(2,ist) * &
                auxrad_stream%fac4_3(ist,i) * &
                transmission_aux%solar_path2%tau_level(lev,ist,i)
              
              fac_3_k(4) = &!fac4_3_k + &
                temp(1,ist) * temp(2,ist) * &
                transmission_aux%solar_path2%tau_level(lev,ist,i) * &
                transmission_aux%solar_path2%od_sfrac(ist,i)
            ENDIF
            
            IF(transmission_aux%solar_path1%tau_level(lev,ist,i) > min_tau) THEN
              
              transmission_aux_k%solar_path1%tau_surf(ist,i) = transmission_aux_k%solar_path1%tau_surf(ist,i) - &
                temp(1,ist) * temp(3,ist) * temp(4,ist) * &
                transmission_aux%solar_path2%tau_level(lev,ist,i) * &
                transmission_aux%solar_path1%tau_level(lev,ist,i)
              
              transmission_aux_k%solar_path1%tau_level(lev,ist,i) = &
                transmission_aux_k%solar_path1%tau_level(lev,ist,i) - &
                temp(1,ist) * temp(3,ist) * temp(4,ist) * &
                transmission_aux%solar_path2%tau_level(lev,ist,i) * &
                transmission_aux%solar_path1%tau_surf(ist,i)
              
              transmission_aux_k%solar_path2%tau_level(lev,ist,i) = &
                transmission_aux_k%solar_path2%tau_level(lev,ist,i) + &
                temp(1,ist) * temp(2,ist) * temp(4,ist)
              
              fac_3_k(2) = &!fac2_3_k + &
                auxrad_stream_k%meanrad_down_solar(ist,i) * auxrad_stream%fac6_3(ist,i) * &
                temp(2,ist) * temp(4,ist) * transmission_aux%solar_path2%tau_level(lev,ist,i)
              
              fac_3_k(6) = &!fac6_3_k + &
                auxrad_stream_k%meanrad_down_solar(ist,i) * auxrad_stream%fac2_3(ist,i) * &
                temp(2,ist) * temp(4,ist) * transmission_aux%solar_path2%tau_level(lev,ist,i)
            ENDIF
            
          ELSE
            
            fac_3_k(2) = 0._jprb
            fac_3_k(4) = 0._jprb
            fac_3_k(5) = 0._jprb
            fac_3_k(6) = 0._jprb
            
          ENDIF ! refl(i) > 0.
          
          !--------------Upward single scattering of the solar beam---------------------------------
          
          isti = ircld%icldarr(ist,lay,prof)

          auxrad_stream_k%up_solar(lay,ist,i) = auxrad_stream_k%up_solar(lay,ist,i) + &
            auxrad_stream_k%meanrad_up_solar(ist,i)
          
          transmission_aux_k%solar_path2%tau_level(lev,ist,i) = &
            transmission_aux_k%solar_path2%tau_level(lev,ist,i) + &
            auxrad_stream_k%meanrad_up_solar(ist,i) * &
            auxrad_stream%fac1_3(ist,i) * auxrad_stream%fac2_3(ist,i) * auxrad_stream%fac3_3(ist,i)
          
          transmission_aux_k%solar_path2%tau_surf(ist,i) = &
            transmission_aux_k%solar_path2%tau_surf(ist,i) - &
            auxrad_stream_k%meanrad_up_solar(ist,i) * &
            auxrad_stream%fac1_3(ist,i) * auxrad_stream%fac2_3(ist,i) * auxrad_stream%fac3_3(ist,i)
          
          fac_3_k(3) = &!fac3_3_k + &
            auxrad_stream_k%meanrad_up_solar(ist,i) * &
            auxrad_stream%fac1_3(ist,i) * auxrad_stream%fac2_3(ist,i) * &
            (transmission_aux%solar_path2%tau_level(lev,ist,i) - &
            transmission_aux%solar_path2%tau_surf(ist,i))
          
          fac_3_k(2) = fac_3_k(2) + &
            auxrad_stream_k%meanrad_up_solar(ist,i) * &
            auxrad_stream%fac1_3(ist,i) * auxrad_stream%fac3_3(ist,i) * &
            (transmission_aux%solar_path2%tau_level(lev,ist,i) - &
            transmission_aux%solar_path2%tau_surf(ist,i))
          
          fac_3_k(1) = &!fac1_3_k + &
            auxrad_stream_k%meanrad_up_solar(ist,i) * &
            auxrad_stream%fac2_3(ist,i) * auxrad_stream%fac3_3(ist,i) * &
            (transmission_aux%solar_path2%tau_level(lev,ist,i) - &
            transmission_aux%solar_path2%tau_surf(ist,i))
          
          transmission_scatt_ir_k%azphacdo(isti,lay,i) = &
            transmission_scatt_ir_k%azphacdo(isti,lay,i) + &
            fac_3_k(6) * solar_spectrum(i) * z4pi_r
          
          transmission_aux_k%solar_path1%od_sfrac(ist,i) = &
            transmission_aux_k%solar_path1%od_sfrac(ist,i) - &
            fac_3_k(5) * auxrad_stream%fac5_3(ist,i)
          
          transmission_aux_k%solar_path2%od_sfrac(ist,i) = &
            transmission_aux_k%solar_path2%od_sfrac(ist,i) - &
            fac_3_k(4) * auxrad_stream%fac4_3(ist,i)
          
          raytracing_k%pathsat(lay,i) = raytracing_k%pathsat(lay,i) + &
            fac_3_k(3) / raytracing%patheff(lay,prof)

          raytracing_k%patheff(lay,i) = raytracing_k%patheff(lay,i) - &
            fac_3_k(3) * raytracing%pathsat(lay,prof) / raytracing%patheff(lay,prof) ** 2

          transmission_scatt_ir_k%ssa(isti,lay,i) = transmission_scatt_ir_k%ssa(isti,lay,i) + &
            fac_3_k(2)
          
          transmission_scatt_ir_k%azphacup(isti,lay,i) = transmission_scatt_ir_k%azphacup(isti,lay,i) + &
            fac_3_k(1) * solar_spectrum(i) * z4pi_r
          
        ENDDO
      ENDIF
    ENDDO
  END SUBROUTINE solar_scattering_near_surf_k
  
  SUBROUTINE solar_rayleigh_k(narray, iv3lev, opts, coef, solar, solar_spectrum, chanprof, &
    raytracing, raytracing_k, ircld, profiles, profiles_k, &
    profiles_dry, profiles_dry_k, &
    transmission_aux, transmission_aux_k, auxrad_stream_k)
    
    USE rttov_types, ONLY : rttov_coef, rttov_chanprof, raytracing_type, profile_type, transmission_type_aux, &
      radiance_aux, rttov_options, ircld_type
    
    USE rttov_const, ONLY : z4pi_r, deg2rad, gravity, na, Mh2o, Mair, &
      ray_scs_wlm, ray_scs_a1, ray_scs_b1, ray_scs_c1, ray_scs_d1, &
      ray_scs_a2, ray_scs_b2, ray_scs_c2, ray_scs_d2, ray_min_wvn
    
    USE parkind1, ONLY : jprb, jplm, jpim
    
    IMPLICIT NONE
    
    INTEGER(jpim), INTENT(in)                  :: narray(:)
    INTEGER(jpim), INTENT(in)                  :: iv3lev(:)
    TYPE(rttov_options), INTENT(in)            :: opts
    TYPE(rttov_coef), INTENT(in)               :: coef
    LOGICAL(jplm), INTENT(in)                  :: solar(:)
    REAL(jprb), INTENT(in)                     :: solar_spectrum(:)
    TYPE(rttov_chanprof), INTENT(in)           :: chanprof(:)
    TYPE(raytracing_type), INTENT(in)          :: raytracing
    TYPE(raytracing_type), INTENT(inout)       :: raytracing_k
    TYPE(ircld_type), INTENT(in)               :: ircld
    TYPE(profile_type), INTENT(in)             :: profiles(:)
    TYPE(profile_type), INTENT(inout)          :: profiles_k(:)
    TYPE(profile_type), INTENT(in)             :: profiles_dry(:)
    TYPE(profile_type), INTENT(inout)          :: profiles_dry_k(:)
    TYPE(transmission_type_aux), INTENT(in)    :: transmission_aux
    TYPE(transmission_type_aux), INTENT(inout) :: transmission_aux_k
    TYPE(radiance_aux), INTENT(inout)          :: auxrad_stream_k
    
    REAL(kind=jprb) :: wlm, ss_param, v_h2o, v_h2o_dry, Mwet
    REAL(kind=jprb) :: cossat, cossol, cosazi, cosscata_term1, cosscata_term2, cosscata
    REAL(kind=jprb) :: ray_phase, solar_src, solar_src_updn
    REAL(kind=jprb) :: v_h2o_k, v_h2o_dry_k, Mwet_k
    REAL(kind=jprb) :: cosscata_term1_k, cosscata_term2_k, cosscata_k
    REAL(kind=jprb) :: ray_phase_k, solar_src_k, solar_src_updn_k
    
    INTEGER(kind=jpim) :: ist, nchannels, nlevels, i, lev, lay
    REAL(kind=jprb) :: rayrad_up_k(0:narray(1)-1,0:narray(2)), rayrad_dn_k(0:narray(1)-1,0:narray(2))
    
    nchannels = narray(3)
    nlevels = narray(1)
    
    DO i = 1, nchannels
      
      IF(solar(i) .AND. coef%ff_cwn(chan) > ray_min_wvn) THEN
        
        ! Calculate layer-independent scattering parameter
        wlm = 10000.0_jprb / coef%ff_cwn(chan)    ! Wavelength in microns
        
        IF (wlm < ray_scs_wlm) THEN
          ss_param = ray_scs_a1 * wlm ** (ray_scs_b1 + ray_scs_c1 * wlm + ray_scs_d1/wlm)
        ELSE
          ss_param = ray_scs_a2 * wlm ** (ray_scs_b2 + ray_scs_c2 * wlm + ray_scs_d2/wlm)
        ENDIF
        ss_param = ss_param * 0.01_jprb ** 2_jpim * na * z4pi_r / gravity
        
        rayrad_up_k(:,:) = 0._jprb
        rayrad_dn_k(:,:) = 0._jprb
        
        ! Add the contribution from the part-layer above the surface
        IF (opts%rt_all%use_q2m) THEN
          v_h2o_dry = 0.5_jprb * (profiles_dry(prof)%q(iv3lev(i)) + profiles_dry(prof)%s2m%q) * 1.E-06_jprb
        ELSE
          v_h2o_dry = profiles_dry(prof)%q(iv3lev(i)) * 1.E-06_jprb
        ENDIF
        IF (profiles(prof)%gas_units /= gas_unit_compatibility) THEN
          v_h2o = v_h2o_dry / (1._jprb + v_h2o_dry)
        ELSE
          v_h2o = v_h2o_dry
        ENDIF
        Mwet = ((1._jprb  - v_h2o) * Mair + v_h2o * Mh2o) * 1.E-3_jprb
        
        cossat = 1._jprb - raytracing%zasat(iv2lay(i), prof) * raytracing%zasat(iv2lay(i), prof)
        cossol = 1._jprb - raytracing%zasun(iv2lay(i), prof) * raytracing%zasun(iv2lay(i), prof)
        cosscata_term1 = SQRT(cossat * cossol)
        cosazi = COS((profiles(prof)%azangle - profiles(prof)%sunazangle)*deg2rad)
        cosscata_term2 = raytracing%zasat(iv2lay(i), prof) * raytracing%zasun(iv2lay(i), prof) * cosazi
        
        solar_src = solar_spectrum(i) * & ! mW m^-2 (cm^-1)^-1
          ABS(profiles(prof)%s2m%p - profiles(prof)%p(iv3lev(i))) * 100.0_jprb * &  ! convert hPa to Pa
          raytracing%pathsat(iv2lay(i), prof) * ss_param / Mwet
        cosscata = - cosscata_term1 - cosscata_term2
        ray_phase = 0.75_jprb * (1.0_jprb + cosscata * cosscata)
        solar_src_updn = solar_src * ray_phase
        
        cosscata_term1_k = 0._jprb
        cosscata_term2_k = 0._jprb
        solar_src_updn_k = 0._jprb
        
        DO ist = 0, ircld%nstream(prof)
          solar_src_updn_k = solar_src_updn_k + &
            auxrad_stream_k%meanrad_down_solar(ist,i) * &
            transmission_aux%solar_path2%Tau_level(iv3lev(i),ist,i) / &
            transmission_aux%solar_path1%Tau_level(iv3lev(i),ist,i) ** 3_jpim

          transmission_aux_k%solar_path2%Tau_level(iv3lev(i),ist,i) = &
            transmission_aux_k%solar_path2%Tau_level(iv3lev(i),ist,i) + &
            auxrad_stream_k%meanrad_down_solar(ist,i) * &
            solar_src_updn / transmission_aux%solar_path1%Tau_level(iv3lev(i),ist,i) ** 3_jpim

          transmission_aux_k%solar_path1%Tau_level(iv3lev(i),ist,i) = &
            transmission_aux_k%solar_path1%Tau_level(iv3lev(i),ist,i) - &
            auxrad_stream_k%meanrad_down_solar(ist,i) * &
            solar_src_updn * 3_jpim * transmission_aux%solar_path2%Tau_level(iv3lev(i),ist,i) / &
            transmission_aux%solar_path1%Tau_level(iv3lev(i),ist,i) ** 4_jpim

          solar_src_updn_k = solar_src_updn_k + auxrad_stream_k%meanrad_up_solar(ist,i) * &
            transmission_aux%solar_path2%Tau_level(iv3lev(i),ist,i)

          transmission_aux_k%solar_path2%Tau_level(iv3lev(i),ist,i) = &
            transmission_aux_k%solar_path2%Tau_level(iv3lev(i),ist,i) + &
            solar_src_updn * auxrad_stream_k%meanrad_up_solar(ist,i)
        ENDDO
        
        ray_phase_k = solar_src * solar_src_updn_k  ! First use of ray_phase_k, no accumulation
        solar_src_k = solar_src_updn_k * ray_phase  ! First use of solar_src_k
        cosscata_k = 2._jprb * 0.75_jprb * ray_phase_k * cosscata ! First use of cossacata_k
        cosscata_term1_k = cosscata_term1_k - cosscata_k
        cosscata_term2_k = cosscata_term2_k - cosscata_k
        
        ! First use of Mwet_k, no accumulation
        Mwet_k = solar_src_k * solar_spectrum(i) * ss_param * &
          raytracing%pathsat(iv2lay(i), prof) * &
          ABS(profiles(prof)%s2m%p - profiles(prof)%p(iv3lev(i))) * 100.0_jprb * &
          (-1._jprb) / (Mwet * Mwet)
        
        raytracing_k%pathsat(iv2lay(i), i) = raytracing_k%pathsat(iv2lay(i), i) + &
          solar_src_k * solar_spectrum(i) * ss_param / Mwet * &
          ABS(profiles(prof)%s2m%p - profiles(prof)%p(iv3lev(i))) * 100.0_jprb
        
        ! Note sign of accumulation in the two pressure terms
        IF (profiles(prof)%s2m%p >= profiles(prof)%p(iv3lev(i))) THEN
          profiles_k(i)%s2m%p = profiles_k(i)%s2m%p + &
            solar_src_k * solar_spectrum(i) * ss_param * &
            raytracing%pathsat(iv2lay(i), prof) * 100.0_jprb / Mwet
          
          IF (opts%interpolation%lgradp) THEN
            profiles_k(i)%p(iv3lev(i)) = profiles_k(i)%p(iv3lev(i)) - &
              solar_src_k * solar_spectrum(i) * ss_param * &
              raytracing%pathsat(iv2lay(i), prof) * 100.0_jprb / Mwet
          ENDIF
        ELSE
          profiles_k(i)%s2m%p = profiles_k(i)%s2m%p - &
            solar_src_k * solar_spectrum(i) * ss_param * &
            raytracing%pathsat(iv2lay(i), prof) * 100.0_jprb / Mwet
          
          IF (opts%interpolation%lgradp) THEN
            profiles_k(i)%p(iv3lev(i)) = profiles_k(i)%p(iv3lev(i)) + &
              solar_src_k * solar_spectrum(i) * ss_param * &
              raytracing%pathsat(iv2lay(i), prof) * 100.0_jprb / Mwet
          ENDIF
        ENDIF
        
        raytracing_k%zasat(iv2lay(i), i) = raytracing_k%zasat(iv2lay(i), i) + &
          cosscata_term2_k * raytracing%zasun(iv2lay(i), prof) * cosazi
        raytracing_k%zasun(iv2lay(i), i) = raytracing_k%zasun(iv2lay(i), i) + &
          cosscata_term2_k * raytracing%zasat(iv2lay(i), prof) * cosazi
        
        raytracing_k%zasat(iv2lay(i), i) = raytracing_k%zasat(iv2lay(i), i) - &
          cosscata_term1_k * raytracing%zasat(iv2lay(i), prof) * cossol / cosscata_term1
        raytracing_k%zasun(iv2lay(i), i) = raytracing_k%zasun(iv2lay(i), i) - &
          cosscata_term1_k * raytracing%zasun(iv2lay(i), prof) * cossat / cosscata_term1
        
        ! First use of v_h2o_k, no accumulation
        v_h2o_k = Mwet_k * (Mh2o - Mair) * 1.E-3_jprb

        IF (profiles(prof)%gas_units /= gas_unit_compatibility) THEN
          v_h2o_dry_k = v_h2o_k * v_h2o * ( 1._jprb - v_h2o / (1._jprb + v_h2o_dry)) / v_h2o_dry
        ELSE
          v_h2o_dry_k = v_h2o_k
        ENDIF

        IF (opts%rt_all%use_q2m) THEN
          profiles_dry_k(i)%q(iv3lev(i)) = profiles_dry_k(i)%q(iv3lev(i)) + &
            0.5_jprb * 1.E-06_jprb * v_h2o_dry_k
          profiles_dry_k(i)%s2m%q = profiles_dry_k(i)%s2m%q + &
            0.5_jprb * 1.E-06_jprb * v_h2o_dry_k
        ELSE
          profiles_dry_k(i)%q(iv3lev(i)) = profiles_dry_k(i)%q(iv3lev(i)) + &
            v_h2o_dry_k * 1.E-06_jprb
        ENDIF
        
        ! Atmospheric contribution
        DO ist = 0, ircld%nstream(prof)
          rayrad_dn_k(iv3lay(i),ist) = rayrad_dn_k(iv3lay(i),ist) + auxrad_stream_k%meanrad_down_solar(ist,i)

          rayrad_dn_k(1:nlayers,ist) = rayrad_dn_k(1:nlayers,ist) + auxrad_stream_k%down_solar(:,ist,i)

          rayrad_up_k(iv3lay(i),ist) = rayrad_up_k(iv3lay(i),ist) + auxrad_stream_k%meanrad_up_solar(ist,i)
          
          rayrad_up_k(1:nlayers,ist) = rayrad_up_k(1:nlayers,ist) + auxrad_stream_k%up_solar(:,ist,i)
        ENDDO
        
        DO lev = nlevels, 2, -1
          lay = lev - 1_jpim
          
          v_h2o_dry = 0.5_jprb * (profiles_dry(prof)%q(lev-1) + profiles_dry(prof)%q(lev)) * 1.E-06_jprb
          IF (profiles(prof)%gas_units /= gas_unit_compatibility) THEN
            v_h2o = v_h2o_dry / (1._jprb + v_h2o_dry)
          ELSE
            v_h2o = v_h2o_dry
          ENDIF
          Mwet = ((1._jprb  - v_h2o) * Mair + v_h2o * Mh2o) * 1.E-3_jprb
          
          cossat = 1._jprb - raytracing%zasat(lay, prof) * raytracing%zasat(lay, prof)
          cossol = 1._jprb - raytracing%zasun(lay, prof) * raytracing%zasun(lay, prof)
          cosscata_term1 = SQRT(cossat * cossol)
          cosazi = COS((profiles(prof)%azangle - profiles(prof)%sunazangle)*deg2rad)
          cosscata_term2 = raytracing%zasat(lay, prof) * raytracing%zasun(lay, prof) * cosazi
          
          solar_src = solar_spectrum(i) * ss_param * &
            (profiles(prof)%p(lev) - profiles(prof)%p(lev-1)) * 100.0_jprb * &
            raytracing%pathsat(lay, prof) / Mwet
          cosscata = - cosscata_term1 - cosscata_term2
          ray_phase = 0.75_jprb * (1.0_jprb + cosscata * cosscata)
          solar_src_updn = solar_src * ray_phase
          
          cosscata_term1_k = 0._jprb
          cosscata_term2_k = 0._jprb
          solar_src_updn_k = 0._jprb
          
          DO ist = 0, ircld%nstream(prof)
            solar_src_updn_k = solar_src_updn_k + rayrad_dn_k(lay,ist) * &
              transmission_aux%solar_path2%Tau_level(lev-1,ist,i) / &
              transmission_aux%solar_path1%Tau_level(lev-1,ist,i) ** 3_jpim

            transmission_aux_k%solar_path2%Tau_level(lev-1,ist,i) = &
              transmission_aux_k%solar_path2%Tau_level(lev-1,ist,i) + &
              rayrad_dn_k(lay,ist) * solar_src_updn / &
              transmission_aux%solar_path1%Tau_level(lev-1,ist,i) ** 3_jpim

            transmission_aux_k%solar_path1%Tau_level(lev-1,ist,i) = &
              transmission_aux_k%solar_path1%Tau_level(lev-1,ist,i) - &
              rayrad_dn_k(lay,ist) * solar_src_updn * &
              3_jpim * transmission_aux%solar_path2%Tau_level(lev-1,ist,i) / &
              transmission_aux%solar_path1%Tau_level(lev-1,ist,i) ** 4_jpim

            rayrad_dn_k(lay-1,ist) = rayrad_dn_k(lay-1,ist) + rayrad_dn_k(lay,ist)

            solar_src_updn_k = solar_src_updn_k + rayrad_up_k(lay,ist) * &
              transmission_aux%solar_path2%Tau_level(lev-1,ist,i)

            transmission_aux_k%solar_path2%Tau_level(lev-1,ist,i) = &
              transmission_aux_k%solar_path2%Tau_level(lev-1,ist,i) + &
              solar_src_updn * rayrad_up_k(lay,ist)

            rayrad_up_k(lay-1,ist) = rayrad_up_k(lay-1,ist) + rayrad_up_k(lay,ist)
          ENDDO
          
          ray_phase_k = solar_src * solar_src_updn_k  ! First use of ray_phase_k, no accumulation
          solar_src_k = solar_src_updn_k * ray_phase  ! First use of solar_src_k
          cosscata_k = 2._jprb * 0.75_jprb * ray_phase_k * cosscata ! First use of cossacata_k
          cosscata_term1_k = cosscata_term1_k - cosscata_k
          cosscata_term2_k = cosscata_term2_k - cosscata_k
          
          raytracing_k%pathsat(lay, i) = raytracing_k%pathsat(lay, i) + &
            solar_src_k * solar_spectrum(i) * ss_param * &
            (profiles(prof)%p(lev) - profiles(prof)%p(lev-1)) * 100.0_jprb / Mwet
          
          ! First use of Mwet_k, no accumulation
          Mwet_k = solar_src_k * solar_spectrum(i) * ss_param * &
            (profiles(prof)%p(lev) - profiles(prof)%p(lev-1)) * 100.0_jprb * &
            raytracing%pathsat(lay, prof) * (-1._jprb) / (Mwet * Mwet)
          
          IF (opts%interpolation%lgradp) THEN
            profiles_k(i)%p(lev) = profiles_k(i)%p(lev) + &
              solar_src_k * solar_spectrum(i) * ss_param * &
              100.0_jprb * raytracing%pathsat(lay, prof) / Mwet
            profiles_k(i)%p(lev-1) = profiles_k(i)%p(lev-1) - &
              solar_src_k * solar_spectrum(i) * ss_param * &
              100.0_jprb * raytracing%pathsat(lay, prof) / Mwet
          ENDIF
          
          raytracing_k%zasat(lay, i) = raytracing_k%zasat(lay, i) + &
            cosscata_term2_k * raytracing%zasun(lay, prof) * cosazi
          raytracing_k%zasun(lay, i) = raytracing_k%zasun(lay, i) + &
            cosscata_term2_k * raytracing%zasat(lay, prof) * cosazi
          
          raytracing_k%zasat(lay, i) = raytracing_k%zasat(lay, i) - &
            cosscata_term1_k * raytracing%zasat(lay, prof) * cossol / cosscata_term1
          raytracing_k%zasun(lay, i) = raytracing_k%zasun(lay, i) - &
            cosscata_term1_k * raytracing%zasun(lay, prof) * cossat / cosscata_term1
          
          v_h2o_k = Mwet_k * (Mh2o - Mair) * 1.E-3_jprb
          IF (profiles(prof)%gas_units /= gas_unit_compatibility) THEN
            v_h2o_dry_k = v_h2o_k * v_h2o * ( 1._jprb - v_h2o / (1._jprb + v_h2o_dry)) / v_h2o_dry
          ELSE
            v_h2o_dry_k = v_h2o_k
          ENDIF

          profiles_dry_k(i)%q(lev-1) = profiles_dry_k(i)%q(lev-1) + &
            0.5_jprb * v_h2o_dry_k * 1.E-06_jprb
          profiles_dry_k(i)%q(lev) = profiles_dry_k(i)%q(lev) + &
            0.5_jprb * v_h2o_dry_k * 1.E-06_jprb
        ENDDO
        
      ENDIF
    ENDDO
    
  END SUBROUTINE solar_rayleigh_k
  
  SUBROUTINE solar_surface_contribution_k(transmission_aux, transmission_aux_k, ircld, &
    reflectance, reflectance_k, refl_norm, chanprof, solar, &
    solar_spectrum, narray, auxrad_stream_k, rad_k)
    
    USE rttov_types, ONLY : rttov_chanprof, ircld_type, transmission_type_aux, radiance_type, radiance_aux
    USE parkind1, ONLY : jpim, jprb, jplm
    
    IMPLICIT NONE
    
    TYPE(transmission_type_aux), INTENT(in)       :: transmission_aux
    TYPE(ircld_type),     INTENT(in)              :: ircld
    REAL(jprb), INTENT(in)                        :: reflectance(:)
    REAL(jprb), INTENT(in)                        :: refl_norm(:)
    TYPE(rttov_chanprof), INTENT(in)              :: chanprof(:)
    LOGICAL(jplm), INTENT(in)                     :: solar(:)
    REAL(jprb), INTENT(in)                        :: solar_spectrum(:)
    INTEGER(jpim), INTENT(in)                     :: narray(:)
    
    TYPE(radiance_type), INTENT(in)               :: rad_k
    TYPE(radiance_aux), INTENT(in)                :: auxrad_stream_k
    
    TYPE(transmission_type_aux), INTENT(inout)    :: transmission_aux_k
    REAL(jprb), INTENT(inout)                     :: reflectance_k(:)
    
    INTEGER(jpim) :: nlevels, nlayers, nstreams, nchannels, i, ist
    REAL(jprb)    :: temp_k(0:narray(2)) ! ist
    
    ! unpack
    nlevels = narray(1)
    nlayers = nlevels - 1_jpim
    nchannels = narray(3)
    
    DO i = 1, nchannels
      IF(solar(i)) THEN
        nstreams = ircld%nstream(prof)
        
        temp_k(0) = rad_k%clear(i)
        temp_k(1:nstreams) = auxrad_stream_k%cloudy(1:nstreams,i)
        
        DO ist = 0, nstreams
          reflectance_k(i) = reflectance_k(i) + solar_spectrum(i) * refl_norm(i) * &
            temp_k(ist) * transmission_aux%solar_path2%Tau_surf(ist,i)
          
          transmission_aux_k%solar_path2%Tau_surf(ist,i) = transmission_aux_k%solar_path2%Tau_surf(ist,i) + &
            solar_spectrum(i) * refl_norm(i) * temp_k(ist) * reflectance(i)
        ENDDO
        
      ENDIF
    ENDDO
  END SUBROUTINE solar_surface_contribution_k
  
END SUBROUTINE rttov_integrate_k
