!     Calculate the secant of the zenith angle at the lower boundary
!     of each atmospheric layer.
SUBROUTINE rttov_locpat_k(opts, dosolar, chanprof, profiles, profiles_k, &
      profiles_dry, profiles_dry_k, aux, coef, &
      angles, raytracing, raytracing_k, do_pmc_calc)
!     Description:
!     The calculation of the secant of the zenith angle at the lower
!     boundary of each atmospheric layer is performed by evaluating
!     the height of the pressure levels.
!
!     Method.
!     The calculation of the height of the pressure levels is based
!     on the hydrostatic equation. To account for the presence of
!     water vapour, virtual temperature are used. The variation of
!     gravity with latitude is introduced using the international
!     gravity formula. The bending of rays as they traverse the atmosphere
!     is fully accounted for applying the Snell's law. For the computation
!     of the refractive index of air, an updated version of Edlen's formula
!     is used.
!
!     K.P. Birch and M.J.Downs:'An Updated Edlen Equation for the
!     refractive index of air'. Metrologia, 1993, 30, 155-162.
!
!     K.P. Birch and M.J.Downs:'Correction to the Updated Edlen Equation
!     for the refractive index of air'. Metrologia, 1994, 31, 315-316.
!
!     Owner:
!     EUMETSAT
!
!     History:
!     Version      Date        Comment
!     1            14/02/2003  Oroginal code. RTIASI.4. Marco Matricardi. ECMWF.
!     2            30/07/2004  RTIASI.5. Marco Matricardi. ECMWF.
!                              Cloud variables removed
!     3            27/02/2009  Profile levels to include ToA. All arrays in
!                              RAYTRACING are on levels except the layer arrays
!                              INT, DMAIR, LTICK. Change RATOESAT and LTICK to
!                              to allow agreement with RTTOV-9.
!     4            27/02/2009  Bug corrected in defining DISPCO2 (P. Rayer)
!     5            02/12/2009  The computation of the local zenith angle has been
!                              changed to take into account the extra top level.
!                              PATHSAT and PATHSUN are now layer arrays (Marco Matricardi).
!     6            05/07/2010  Remove addsolar flag from profiles structure (J Hocking)
!     7            15/10/2010  Ensure refraction assumed when addpc true (J Hocking)
!     8            10/01/2013  Add PMC shifts (P Rayer)
!
! 2010/03 Code cleaning Pascal Brunel, Philippe Marguinaud
!
!     Code description:
!       Language:              Fortran 90.
!       Software Standards:    "European Standards for Writing and Documenting
!                              Exchangeable Fortran 90 code".
!     Module used:
  USE rttov_types, ONLY : rttov_options, rttov_coef, profile_aux, profile_type, &
                          geometry_type, raytracing_type, rttov_chanprof
  USE parkind1, ONLY : jplm
!INTF_OFF
  USE rttov_const, ONLY : d1, d2, d3, d4, d5, dco2, ed1, ed2, ed3, ed4, &
                          ew1, ew2, htop, ctom, waver, rgc, mair, mh2o, flatt, &
                          eqrad, omega, grave, t0, p0, max_sol_zen, co2_conc, &
                          gas_unit_compatibility

  USE parkind1, ONLY : jpim, jprb
  USE yomhook, ONLY : LHOOK, DR_HOOK
  USE rttov_math_mod

!INTF_ON
  IMPLICIT NONE
  TYPE(rttov_options  ), INTENT(IN)    :: opts
  LOGICAL(jplm)        , INTENT(IN)    :: dosolar
  TYPE(RTTOV_CHANPROF ), INTENT(IN)    :: chanprof(:)
  TYPE(profile_type   ), INTENT(IN)    :: profiles(:)
  TYPE(profile_type   ), INTENT(INOUT) :: profiles_k(:)
  TYPE(profile_type   ), INTENT(IN)    :: profiles_dry(:)
  TYPE(profile_type   ), INTENT(INOUT) :: profiles_dry_k(:)
  TYPE(profile_aux    ), INTENT(IN)    :: aux
  TYPE(geometry_type  ), INTENT(IN)    :: angles(:)
  TYPE(rttov_coef     ), INTENT(IN)    :: coef
  TYPE(raytracing_type), INTENT(IN)    :: raytracing
  TYPE(raytracing_type), INTENT(INOUT) :: raytracing_k
  LOGICAL(jplm)        , INTENT(IN), OPTIONAL :: do_pmc_calc
!INTF_END

  INTEGER(kind=jpim) :: lev, lay, prof, i, iprof
  INTEGER(kind=jpim) :: nlevels, nlayers
  REAL   (kind=jprb) :: rearth(size(profiles))
  REAL   (kind=jprb) :: gravl(SIZE(profiles))   ! gravity at latitude lat [m/s^2]
  REAL   (kind=jprb) :: gravh(SIZE(profiles)), gravh_r(SIZE(profiles))
  REAL   (kind=jprb) :: eta, beta ! coefficients of the international gravity formula
  REAL   (kind=jprb) :: rlh(SIZE(profiles))
  REAL   (kind=jprb) :: dflat
  REAL   (kind=jprb) :: fac, c

  REAL   (kind=jprb) :: pres_k(profiles(1)%nlevels,SIZE(chanprof))
  REAL   (kind=jprb) :: temp_k(profiles(1)%nlevels,size(chanprof))
  REAL   (KIND=jprb) :: qwet(profiles(1)%nlevels)
  REAL   (KIND=jprb) :: qwet_k(profiles(1)%nlevels)

  INTEGER(kind=jpim) :: nchannels, nprofiles
  LOGICAL(KIND=jplm) :: do_pmc_calc1

  REAL   (kind=jprb) :: scale_eps

  INTEGER(jpim) :: map(SIZE(chanprof),2), prof_stat, last_prof
  REAL(jprb) :: ztemp(profiles(1)%nlevels, 4)

  REAL   (KIND=JPRB) :: ZHOOK_HANDLE

!-----end of header-------------------------------------------------------------
  IF (LHOOK) CALL DR_HOOK('RTTOV_LOCPAT_K', 0_jpim, ZHOOK_HANDLE)
  
  nchannels = SIZE(chanprof)
  nprofiles = SIZE(profiles)
  nlevels     = profiles(1)%nlevels
  nlayers     = nlevels - 1

  map(1,1) = chanprof(1)%prof
  map(1,2) = chanprof(1)%chan

  prof_stat = 1 ! assume profs are contiguous and monotonic
  DO i = 2, nchannels
    map(i,1) = chanprof(i)%prof
    map(i,2) = chanprof(i)%chan

    IF(map(i,1) < map(i-1,1)) THEN ! they are not.
      prof_stat = -1
    ENDIF
  ENDDO

  do_pmc_calc1 = .FALSE.
  IF (PRESENT(do_pmc_calc)) do_pmc_calc1 = do_pmc_calc

! DAR: Move constant calculations to top of locpat subroutine

! scale water profile in air density calculation
  scale_eps = 1.0_jprb - mh2o/mair

  dflat = (1.0_jprb - flatt) ** 2_jpim
  fac = (omega ** 2_jpim * (eqrad * 1000._jprb)) / (grave)
  beta = 5._jprb * fac / 2._jprb - flatt - 17._jprb * fac * flatt / 14._jprb
  eta = flatt * (5._jprb * fac - flatt) / 8._jprb

!-------1.  caculate the earth's radius at latitude lat assuming the-----------
!       earth is an ellipsoid of revolution                                    |
!------------------------------------------------------------------------------
  rearth(:) = SQRT(eqrad ** 2_jpim * dflat / &
               (angles(:)%sinlat ** 2_jpim + &
                dflat * angles(:)%coslat ** 2_jpim))

  rlh(:) = rearth(:) * 5e2_jprb ! DAR: used to be gravl / gravh

!-------2.  the value of earth's gravity at surface at latitude lat is---------
!       computed using the international gravity formula.                      |
!------------------------------------------------------------------------------
  gravl(:) = grave * &
    (1._jprb + beta * angles(:)%sinlat ** 2_jpim +      &
    eta * (2._jprb * angles(:)%sinlat * angles(:)%coslat) ** 2_jpim)

  gravh_r(:) = rearth(:) / (2.0e-3_jprb * gravl(:)) ! computational useful form

!-------3. calculate reciprocal of moist air density (for later calculation)
! density of dry air is adjusted to account for the presence of water vapour 
! by replacing temp with virtual temp

  c = 1000.0_jprb * rgc / (100.0_jprb * mair)

  IF (coef%pmc_shift .AND. do_pmc_calc1) CALL calc_ssu_co2_thickness_k()
!Start K
  temp_k = 0._jprb
  IF (opts%interpolation%LGRADP) pres_k = 0._jprb

! Extra bits

  IF (opts%rt_ir%addaerosl .OR. opts%rt_ir%addclouds) THEN
    DO i = 1, nchannels
      DO LAY = NLEVELS - 1, 1,  - 1
        RAYTRACING_K%HGPL(LAY, I)     = RAYTRACING_K%HGPL(LAY, I) + RAYTRACING_K%LTICK(LAY, I)
        RAYTRACING_K%HGPL(LAY + 1, I) = RAYTRACING_K%HGPL(LAY + 1, I) - RAYTRACING_K%LTICK(LAY, I)
      ENDDO
    ENDDO
  ENDIF

  DO i = 1, nchannels
    RAYTRACING_K%HGPL(2, I) = RAYTRACING_K%HGPL(2, I) + RAYTRACING_K%HGPL(1, I)
  ENDDO

!-------6.  Compute secant of the zenith angle at the lower boundary of each -
!       layer for the satellite-to-surface line of view and of each             |
!       layer for the sun-to-surface line of view (if possible and reqd.)       |

  IF(opts%rt_all%ADDREFRAC .OR. opts%rt_ir%pc%ADDPC) THEN
    raytracing_k%r = 0._jprb ! moved from init_raytracing
  ELSE
    raytracing_k%ppw = 0._jprb
  ENDIF

  IF (dosolar) THEN
    raytracing_k%pathsat = raytracing_k%pathsat + raytracing_k%patheff
    raytracing_k%pathsun = raytracing_k%pathsun + raytracing_k%patheff

    CALL calc_seczen_k(&
      raytracing%ratoesun, raytracing%zasun, raytracing%pathsun,&
      raytracing_k%ratoesun, raytracing_k%zasun, raytracing_k%pathsun,&
      do_sat = .FALSE._jplm) ! solar call
!  ELSE
!    raytracing_k%pathsun(:,:)  = 0._jprb
  ENDIF

  CALL calc_seczen_k(raytracing%ratoesat, raytracing%zasat, raytracing%pathsat,&
    raytracing_k%ratoesat, raytracing_k%zasat, raytracing_k%pathsat, &
    do_sat = .TRUE._jplm) ! satellite call

!-------5. compute atmospheric refractive index (if reqd.)---------------------
! compute the refractivity index for dry air at temperature tempp and pressure
! as a function of wave-number using an updated version of edlen equation
! not yet corrected for CO2 contribution - used in calc_seczen
  IF(opts%rt_all%ADDREFRAC .OR. opts%rt_ir%pc%ADDPC) CALL calc_refractivity_k()

!-------5. Compute height of pressure levels
  CALL calc_hgpl_k()

!-------4.  set up the height of the surface pressure level---------------------
! height of nearest level below surf counted from bottom of profile grid

  DO i = 1, nchannels
    prof = chanprof(i)%prof
    raytracing_k%hgpl(aux%s(prof)%nearestlev_surf, i) = 0._jprb
  ENDDO

  last_prof = -1
  DO i = 1, nchannels
    prof = chanprof(i)%prof
    IF (prof .NE. last_prof) THEN
      CALL reciprocal(&
        profiles(prof)%p(:) - raytracing%ppw(:,prof) * scale_eps, &
        ztemp(:,1))
      ztemp(:,2) = ztemp(:,1) * raytracing%dmair_r(:, prof)

      IF (profiles(prof)%gas_units == gas_unit_compatibility) THEN
        qwet(:) = profiles_dry(prof)%q(:)
      ELSE
        qwet(:) = profiles_dry(prof)%q(:) / (1._jprb + profiles_dry(prof)%q(:) * 1e-6_jprb)
      ENDIF

      last_prof = prof
    ENDIF

    temp_k(:,i) = temp_k(:,i) + c * raytracing_k%dmair_r(:, i) * ztemp(:,1)

    raytracing_k%ppw(:,i) = raytracing_k%ppw(:,i) + &
      raytracing_k%dmair_r(:, i) * scale_eps * ztemp(:,2)

    qwet_k(:) = &!qwet_k(:) + &
      1.0e-6_jprb * profiles(prof)%p(:) * raytracing_k%ppw(:,i)

! first calculate density for dry air (ideal gas)
    IF (opts%interpolation%LGRADP) THEN
!      CALL DIVIDE(-(raytracing%dmair_r(:, prof)) * &
!             raytracing_k%dmair_r(:, i), &
!             ztemp(:,1), &
!             pres_k(:, i), acc = .TRUE._jplm)

      pres_k(:, i) = pres_k(:, i) - raytracing_k%dmair_r(:, i) * ztemp(:,2)

      pres_k(:,i) = pres_k(:,i) + &
        1.0e-6_jprb * raytracing_k%ppw(:,i) * qwet(:)

      profiles_k(i)%p(:) = profiles_k(i)%p(:) + &
        pres_k(:,i)
    ENDIF

    qwet_k(2) = qwet_k(2) + qwet_k(1)
    qwet_k(1) = 0._jprb

    IF (profiles(prof)%gas_units == gas_unit_compatibility) THEN
      profiles_dry_k(i)%q(:) = profiles_dry_k(i)%q(:) + qwet_k(:)
    ELSE
      profiles_dry_k(i)%q(:) = profiles_dry_k(i)%q(:) + &
        qwet_k(:) * qwet(:) * (1._jprb - 1e-6_jprb * qwet(:)) / profiles_dry(prof)%q(:)
    ENDIF

    profiles_k(i)%t(2) = profiles_k(i)%t(2) + temp_k(1,i)
    temp_k(1,i) = 0._jprb

    profiles_k(i)%t(:) = profiles_k(i)%t(:) + temp_k(:,i)
  ENDDO

  IF (LHOOK) CALL DR_HOOK('RTTOV_LOCPAT_K', 1_jpim, ZHOOK_HANDLE)

CONTAINS 

  SUBROUTINE calc_refractivity_k()

    REAL(kind=jprb) :: dry_pp_h((profiles(1)%nlevels))
    REAL(kind=jprb) :: t_norm((profiles(1)%nlevels))
    REAL(kind=jprb) :: dry_pp_h_k((profiles(1)%nlevels))
    REAL(kind=jprb) :: U((profiles(1)%nlevels)), V_R((profiles(1)%nlevels))
    REAL(kind=jprb) :: DU_K((profiles(1)%nlevels)), DV_K((profiles(1)%nlevels))

    REAL(kind=jprb) :: DISP ! the value of the refractive index given by the dispersion equation.
    REAL(kind=jprb) :: DISP_MOIST
    REAL(kind=jprb) :: DISPCO2

! calculate constants

! dry air dispersion equation constant in updated Edlen eqn
    DISP = 1.0e-8_jprb * ( &
      d1 + d2 / (d3 - (waver * ctom) ** 2_jpim) + &
           d4 / (d5 - (waver * ctom) ** 2_jpim))

! moist air dispersion constant
    DISP_MOIST = htop * (ew1 - ew2 * (waver * ctom) ** 2_jpim) * 1.0e-10_jprb

    last_prof = -1
    DO i = 1, nchannels
! Store temporary FWD variables
      prof = chanprof(i)%prof
      IF(prof .ne. last_prof) THEN
        dry_pp_h(:) = htop * &
          (profiles(prof)%p(:) - raytracing%ppw(:, prof))
        t_norm(:) = profiles(prof)%t(:) - t0
        
        U = 1._jprb + &
          1.0e-8_jprb * (ed2 - ed3 * t_norm(:)) * dry_pp_h(:)
        CALL reciprocal(1._jprb + ed4 * t_norm(:), V_R(:))
        
        ztemp(:,1) = U(:) * V_R / ed1
        ztemp(:,2) = dry_pp_h(:) * V_R(:) / ed1
        ztemp(:,3) = -raytracing%refractivity(:, prof) * V_R(:)
        ztemp(:,4) = 1.0e-8_jprb * (ed2 - ed3 * t_norm(:))
        last_prof = prof
      ENDIF

!K
! compute moist contribution to refractivity
! and store final refractive index for atmospheric profile

      raytracing_k%ppw(:, i) = -&! raytracing_k%ppw(:, i) - &!
        (raytracing_k%r(:, i) * DISP_MOIST)

      ! apply CO2 correction
      IF (opts%rt_ir%co2_data) THEN
        raytracing_k%dispco2(:, i) = &! raytracing_k%dispco2(:, i) + &
             raytracing%refractivity(:, prof) * raytracing_k%r(:, i)

        raytracing_k%refractivity(:, i) = &!raytracing_k%refractivity(:, i) + &
          raytracing%dispco2(:, prof) * raytracing_k%r(:, i) 
      ELSE
        DISPCO2 = DISP * (1._jprb + dco2 * (CO2_CONC * 1e-6_jprb - 0.0003_jprb))
        raytracing_k%refractivity(:, i) = &! raytracing_k%refractivity(:, i) +
          raytracing_k%r(:, i) * DISPCO2 ! AD
      ENDIF

      dry_pp_h_k(:) = &!dry_pp_h_k(:) + &
        raytracing_k%refractivity(:, i) * ztemp(:,1)

      DU_k(:) = &!DU_k(:) + &
        raytracing_k%refractivity(:, i) * ztemp(:,2)

      DV_k(:) = &!DV_k(:) &
        raytracing_k%refractivity(:, i) * ztemp(:,3)

      temp_k(:,i) = temp_k(:,i) + &
        DV_k(:) * ed4 + & !combine two lines from TL     
        DU_k(:) * 1.0e-8_jprb * dry_pp_h(:) * (-ed3) 


      dry_pp_h_k(:) = dry_pp_h_k(:) + &
        DU_k(:) * ztemp(:,4)

      raytracing_k%ppw(:, i) = raytracing_k%ppw(:, i) - &
        dry_pp_h_k(:) * htop

      IF(opts%interpolation%LGRADP) THEN
        profiles_k(i)%p(:) = profiles_k(i)%p(:) + &
          htop * dry_pp_h_k(:)
      ENDIF
    ENDDO

    IF (opts%rt_ir%co2_data) THEN
     DO i = 1, nchannels
       profiles_dry_k(i)%co2(:) = profiles_dry_k(i)%co2(:) + &
         raytracing_k%dispco2(:, i) * disp * &
         dco2 * 1e-6_jprb
     ENDDO
   ELSE
!     DO i = 1, nchannels
!       profiles_dry_k(i)%co2(:) = 0._jprb
!     ENDDO
   ENDIF

 END SUBROUTINE calc_refractivity_k

 SUBROUTINE calc_hgpl_k()
   REAL(kind=jprb) :: dp((profiles(1)%nlevels))
   REAL(kind=jprb) :: dp_k((profiles(1)%nlevels),SIZE(chanprof))
!    REAL(jprb) :: ztemp_t(3,profiles(1)%nlevels)
   REAL(kind=jprb) :: ztemp((profiles(1)%nlevels),SIZE(profiles))

! Init   
!   raytracing_k%INT      = 0._jprb

   CALL INVSQRT(raytracing%ztemp(:, :), ztemp(:, :))

   last_prof = -1
   DO i = 1, nchannels
     prof = chanprof(i)%prof
     IF(prof .NE. last_prof) THEN     
!       dp(1:nlevels-1) = 0.5_jprb * &
!         (profiles(prof)%p(nlevels-1:1:-1) - profiles(prof)%p(nlevels:2:-1))
       DO lev = 1, nlevels - 1
         dp(lev) = 0.5_jprb * (profiles(prof)%p(lev) - &
                               profiles(prof)%p(lev + 1))
       ENDDO

!        DO lev = nlevels, aux%s(prof)%nearestlev_surf + 1, -1
!          ztemp_t(1,lev) = - 0.5_jprb * ztemp(lev,prof) * 1e-3_jprb
!          ztemp_t(2,lev) = 2._jprb * &
!                           (1.0e3_jprb * raytracing%hgpl(lev - 1, prof) - rlh(prof))
!          ztemp_t(3,lev) = raytracing%dmair_r(lev, prof) + &
!                           raytracing%dmair_r(lev - 1, prof)
!        ENDDO
!        DO lev = 1, aux%s(prof)%nearestlev_surf - 1
!          ztemp_t(1,lev) = - 0.5_jprb * ztemp(lev, prof) * 1e-3_jprb
!          ztemp_t(2,lev) = 2._jprb * &
!                           (1.0e3_jprb * raytracing%hgpl(lev + 1, prof) - rlh(prof))
!          ztemp_t(3,lev) = raytracing%dmair_r(lev + 1, prof) + &
!                           raytracing%dmair_r(lev, prof)
!        ENDDO
       last_prof = prof
     ENDIF

     DO lev = nlevels, aux%s(prof)%nearestlev_surf + 1, -1

       raytracing_k%ztemp(lev, i) = &!raytracing_k%ztemp(lev, i) -
         -0.5_jprb * ztemp(lev,prof) * 1e-3_jprb * raytracing_k%hgpl(lev, i)

       raytracing_k%INT(lev, i) = &!raytracing_k%INT(lev, i) - &
         -2._jprb * 100._jprb * gravh_r(prof) * raytracing_k%ztemp(lev, i)

       raytracing_k%hgpl(lev - 1, i) = raytracing_k%hgpl(lev - 1, i) + 1.0e3_jprb * &
         2._jprb * raytracing_k%ztemp(lev, i) * &
         (1.0e3_jprb * raytracing%hgpl(lev - 1, prof) - rlh(prof))

       raytracing_k%dmair_r(lev, i) = &
         raytracing_k%dmair_r(lev, i) + &
         raytracing_k%INT(lev, i) * dp(lev - 1)

       raytracing_k%dmair_r(lev - 1, i) = &
         raytracing_k%dmair_r(lev - 1, i) + &
         raytracing_k%INT(lev, i) * dp(lev - 1)

       IF(opts%interpolation%lgradp) THEN
         dp_k(lev - 1, i) = &!dp_k(lev - 1, i) + &
           raytracing_k%INT(lev, i) * &
           (raytracing%dmair_r(lev - 1, prof) + &
           raytracing%dmair_r(lev, prof))
       ENDIF
     ENDDO

     DO lev = 1, aux%s(prof)%nearestlev_surf - 1
       raytracing_k%ztemp(lev, i) = &!raytracing_k%ztemp(lev, i) -
         -0.5_jprb * ztemp(lev,prof) * 1e-3_jprb * raytracing_k%hgpl(lev, i)

       raytracing_k%INT(lev, i) = &!raytracing_k%INT(lev, i) - &
         -2._jprb * 100._jprb * gravh_r(prof) * raytracing_k%ztemp(lev, i)

       raytracing_k%hgpl(lev + 1, i) = raytracing_k%hgpl(lev + 1, i) + 1.0e3_jprb * &
         2._jprb * raytracing_k%ztemp(lev, i) * &
         (1.0e3_jprb * raytracing%hgpl(lev + 1, prof) - rlh(prof))

       raytracing_k%dmair_r(lev + 1, i) = &
         raytracing_k%dmair_r(lev + 1, i) - &
         raytracing_k%INT(lev, i) * dp(lev)

       raytracing_k%dmair_r(lev, i) = &
         raytracing_k%dmair_r(lev, i) - &
         raytracing_k%INT(lev, i) * dp(lev)

       IF(opts%interpolation%lgradp) THEN
         dp_k(lev, i) = &!dp_k(lev, i) - &
           -raytracing_k%INT(lev, i) * &
           (raytracing%dmair_r(lev, prof) + &
            raytracing%dmair_r(lev + 1, prof))
       ENDIF
     ENDDO

     IF(opts%interpolation%LGRADP) THEN
       DO lev = nlevels - 1, 1, -1
         pres_k(lev, i) = pres_k(lev, i) + 0.5_jprb * dp_k(lev, i)
         pres_k(lev + 1, i) = pres_k(lev + 1, i) - 0.5_jprb * dp_k(lev, i)
       ENDDO
     ENDIF

   ENDDO
 END SUBROUTINE calc_hgpl_k

!-------Compute secant of the zenith angle at the lower boundary of each -
!       layer for the sat/sun-to-surface line of view.                       |
  SUBROUTINE calc_seczen_k(ra,za,path,ra_k,za_k,path_k,do_sat)
    LOGICAL(jplm), INTENT(IN) :: do_sat ! do satellite secant or solar secant
    REAL(jprb), INTENT(IN) :: ra(:,:), za(:,:), path(:,:)
    REAL(jprb), INTENT(INOUT) :: ra_k(:,:), za_k(:,:), path_k(:,:)

    REAL(jprb) :: z(SIZE(ra(:,1)),SIZE(ra(1,:)))
    REAL(jprb) :: alt(SIZE(ra(1,:)))
    INTEGER(jpim) :: n

    IF(do_sat .AND. aux%on_coef_levels) THEN
      DO i = 1, nchannels
        prof = chanprof(i)%prof
        raytracing_k%pathsat(:,i) = raytracing_k%pathsat(:,i) + &
          raytracing_k%pathsat_sqrt(:,i) * raytracing%pathsat_rsqrt(:,prof)

        raytracing_k%pathsat_rsqrt(:,i) = &!raytracing_k%pathsat_rsqrt + &
          raytracing%pathsat(:,prof) * raytracing_k%pathsat_sqrt(:,i)
      ENDDO
      CALL INVSQRT_K(raytracing%pathsat_rsqrt, raytracing_k%pathsat, raytracing_k%pathsat_rsqrt, &
        acc=.TRUE._jplm, map = map, prof_stat = prof_stat)
    ENDIF
       
! calculate path from zenith angle using path = 1/sqrt(1-za**2)
    IF(opts%rt_ir%addsolar) THEN
      CALL SINTOSEC_K(za(1:nlayers,:), za_k(1:nlayers,:), path, path_k, acc = .TRUE._jplm, map = map, &
        prof_stat = prof_stat)
    ELSE
      CALL SINTOSEC_K(za(1:nlayers,:), za_k(1:nlayers,:), path, path_k, acc = .FALSE._jplm, map = map, &
        prof_stat = prof_stat) ! first use of za      
    ENDIF
!   path_k = 0._jprb; 

    IF(do_sat) THEN
      n = nlayers - 1
      DO iprof = 1, nprofiles
        alt(iprof) = rearth(iprof) + coef%fc_sat_height
        z(:, iprof) = alt(iprof) * raytracing%z_r(:, iprof)
      ENDDO
    ELSE ! do_sun
      n = 0
      DO iprof = 1, nprofiles
        IF(profiles(iprof)%sunzenangle >= 0.0 .AND. &
          profiles(iprof)%sunzenangle < max_sol_zen) THEN
          alt(iprof) = rearth(iprof)
          z(:, iprof) = alt(iprof) * raytracing%z_r(:, iprof)
        ENDIF
      ENDDO
    ENDIF

    IF(do_sat) THEN          
      DO i = 1, nchannels
        prof = chanprof(i)%prof
        ra_k(:,i) = &!ra_k(:,i) + 
                    za_k(1:nlayers,i) * angles(prof)%sinview
      ENDDO
    ELSE ! do_sun
      DO i = 1, nchannels
        prof = chanprof(i)%prof
        IF(profiles(prof)%sunzenangle >= 0.0 .AND. &
          profiles(prof)%sunzenangle < max_sol_zen) THEN

          ra_k(:,i) = &!ra_k(:,i) + 
                      za_k(1:nlayers,i) * angles(prof)%sinzen_sun
        ELSE
!          za_k(:,i) = 0.0_jprb
        ENDIF
      ENDDO
    ENDIF
!  za_k = 0._jprb;

    IF(opts%rt_all%ADDREFRAC .OR. opts%rt_ir%pc%ADDPC) THEN
      last_prof = -1
      DO i = 1, nchannels
        prof = chanprof(i)%prof
        IF(prof .NE. last_prof) THEN
          ztemp(1:nlayers,1) = z(1:nlayers,prof) * raytracing%r(nlevels - n, prof) * &
                       raytracing%r_r(2:nlevels,prof)**2_jpim
          ztemp(1:nlayers,2) = raytracing%r_r(2:nlevels,prof) * z(1:nlayers,prof)
          ztemp(1:nlayers,3) = raytracing%r_r(2:nlevels,prof) * &
            raytracing%r(nlevels - n, prof)
          last_prof = prof
        ENDIF

! the (nlevels - n) takes in to account the fact that 10.2 and earlier
! had the refraction calculation on unpacked, reversed levels.
! first use but called twice so keep accumulation

        raytracing_k%r(2:nlevels,i) = raytracing_k%r(2:nlevels,i) - & 
          ra_k(:,i) * ztemp(1:nlayers,1)

        raytracing_k%r(nlevels - n,i) = raytracing_k%r(nlevels - n, i) + &
          SUM(ztemp(1:nlayers,2) * ra_k(:,i))

        ra_k(:,i) = ra_k(:,i) * ztemp(1:nlayers,3)
      ENDDO
    ENDIF

    IF (do_sat) THEN
      DO i = 1, nchannels
        prof = chanprof(i)%prof
        raytracing_k%z_r(:,i) = &! raytracing_k%z_r(:,i) + &
          alt(prof) * ra_k(:,i)
      ENDDO

    ELSE ! do sun
      DO i = 1, nchannels
        prof = chanprof(i)%prof
        IF(profiles(prof)%sunzenangle >= 0.0 .AND. &
          profiles(prof)%sunzenangle < max_sol_zen) THEN

! note accumulation as first call is above
          raytracing_k%z_r(:,i) = &!raytracing_k%z_r(:,i) + & 
            alt(prof) * ra_k(:,i)
        ELSE
          raytracing_k%z_r(:,i) = 0._jprb
        ENDIF
      ENDDO
    ENDIF
!   ra_k = 0._jprb;

    CALL reciprocal_k(raytracing%z_r, raytracing_k%hgpl(2:nlayers+1, :), &
      raytracing_k%z_r, acc = .TRUE._jplm, map = map, prof_stat = prof_stat)

  END SUBROUTINE CALC_SECZEN_K
  
  SUBROUTINE calc_ssu_co2_thickness_k
    REAL(kind=jprb) :: pstd
    REAL(kind=jprb) :: tstd
    REAL(kind=jprb) :: pa(profiles(1)%nlevels-1),   pa_k(profiles(1)%nlevels-1)
    REAL(kind=jprb) :: dpa(profiles(1)%nlevels-1),  dpa_k(profiles(1)%nlevels-1)
    REAL(kind=jprb) :: rho(profiles(1)%nlevels-1),  rho_k(profiles(1)%nlevels-1)
    REAL(kind=jprb) :: temp(profiles(1)%nlevels-1), temp_k(profiles(1)%nlevels-1)
    REAL(kind=jprb) :: co2(profiles(1)%nlevels-1),  co2_k(profiles(1)%nlevels - 1)
    REAL(kind=jprb) :: dmair(profiles(1)%nlevels),  dmair_k(profiles(1)%nlevels)
    REAL(kind=jprb) :: dz, dz_k
    REAL(kind=jprb) :: dzstp_k

    pstd=p0*100._jprb
    tstd=t0
    
    gravh(:) = 1._jprb / gravh_r(:)

    DO i = 1, nchannels
      prof = chanprof(i)%prof
      CALL reciprocal(raytracing%dmair_r(:,prof), dmair(:))
      dmair_k(:) = 0._jprb
      DO lay = 1, nlayers
! fwd
        pa(lay)   = 100._jprb * (profiles(prof)%p(lay) + profiles(prof)%p(lay+1)) / 2._jprb
        dpa(lay)  = 100._jprb * (profiles(prof)%p(lay) - profiles(prof)%p(lay+1))
        rho(lay)  = (dmair(lay) + dmair(lay+1)) / 2._jprb
        temp(lay) = (profiles(prof)%t(lay) + profiles(prof)%t(lay+1)) / 2._jprb
        co2(lay)  = (profiles_dry(prof)%co2(lay) + profiles_dry(prof)%co2(lay+1)) / 2._jprb
! ad
        pa_k(lay)   = 0._jprb
        dpa_k(lay)  = 0._jprb
        rho_k(lay)  = 0._jprb
        temp_k(lay) = 0._jprb
        co2_k(lay)  = 0._jprb
      ENDDO
    
      DO lay = 1, nlayers
! fwd
        dz = -dpa(lay)/&
          ( (gravl(prof)-gravh(prof)*raytracing%hgpl(lay+1, prof)) * rho(lay) )
! ad
        dzstp_k = 0._jprb
        dzstp_k = dzstp_k + raytracing_k%co2_cm(i)
        dz_k = 0._jprb
        dz_k =  dz_k + &
          dzstp_k * pa(lay)/pstd * tstd/temp(lay)* co2(lay)*1e-6_jprb*100._jprb

        IF (opts%interpolation%lgradp) THEN
          pa_k(lay) = pa_k(lay) + &
            dzstp_k * dz / pstd * tstd/temp(lay)* co2(lay)*1e-6_jprb*100._jprb
!        ELSE
!          pa_k(lay) = 0._jprb
        ENDIF
        
        temp_k(lay) = temp_k(lay) - &
          dzstp_k * dz * pa(lay)/pstd * tstd/temp(lay)**2_jpim * &
          co2(lay)*1e-6_jprb*100._jprb

        co2_k(lay) = co2_k(lay) + &
          dzstp_k * dz * pa(lay)/pstd * tstd/temp(lay)*1e-6_jprb*100._jprb
!        dzstp_k = 0._jprb
        
        IF (opts%interpolation%lgradp) THEN
          dpa_k(lay) = dpa_k(lay) - dz_k / ( (gravl(prof)-gravh(prof)*raytracing%hgpl(lay+1, prof)) * rho(lay) )
!        ELSE
!          dpa_k(lay) = 0._jprb
        ENDIF

        raytracing_k%hgpl(lay+1, i) = raytracing_k%hgpl(lay+1, i) &
          - dz_k * dpa(lay) * gravh(prof)*rho(lay)/ &
          ( (gravl(prof)-gravh(prof)*raytracing%hgpl(lay+1, prof)) * rho(lay) )**2_jpim

        rho_k(lay) = rho_k(lay)                                              &
          + dz_k * dpa(lay)*(gravl(prof)-gravh(prof)*raytracing%hgpl(lay+1, prof)) / &
          ( (gravl(prof)-gravh(prof)*raytracing%hgpl(lay+1, prof)) * rho(lay) )**2_jpim
        
!        dz_k = 0._jprb
        
        profiles_dry_k(i)%co2(lay) = profiles_dry_k(i)%co2(lay) + &
          co2_k(lay)/2._jprb
        profiles_dry_k(i)%co2(lay+1) = profiles_dry_k(i)%co2(lay+1) + &
          co2_k(lay)/2._jprb
!        co2_k(lay) = 0._jprb

        profiles_k(i)%t(lay) = profiles_k(i)%t(lay) + &
          temp_k(lay)/2._jprb

        profiles_k(i)%t(lay+1) = profiles_k(i)%t(lay+1) + &
          temp_k(lay)/2._jprb
!        temp_k(lay) = 0._jprb
        
        dmair_k(lay) = dmair_k(lay) + rho_k(lay) / 2._jprb
        dmair_k(lay+1) = dmair_k(lay+1) + rho_k(lay) / 2._jprb
       
!        rho_k(lay) = 0._jprb
        
        IF (opts%interpolation%lgradp) THEN
          profiles_k(i)%p(lay) = profiles_k(i)%p(lay) + &
            dpa_k(lay)*100._jprb
          profiles_k(i)%p(lay+1) = profiles_k(i)%p(lay+1) - &
            dpa_k(lay)*100._jprb
!          dpa_k(lay) = 0._jprb
          profiles_k(i)%p(lay) = profiles_k(i)%p(lay) + &
            pa_k(lay)*100._jprb/2._jprb
          profiles_k(i)%p(lay+1) = profiles_k(i)%p(lay+1) + &
            pa_k(lay)*100._jprb/2._jprb
!          pa_k(lay) = 0._jprb
!        ELSE
!          profiles_k(i)%p(lay) = 0._jprb
!          profiles_k(i)%p(lay+1) = 0._jprb
!          dpa_k(lay) = 0._jprb
!          pa_k(lay) = 0._jprb
        ENDIF
      ENDDO
      CALL reciprocal_ad(dmair(:), raytracing_k%dmair_r(:,i), &
                         dmair_k(:), acc = .FALSE._jplm)
    END DO

  END SUBROUTINE calc_ssu_co2_thickness_k

END SUBROUTINE rttov_locpat_k
