!     Calculate the secant of the zenith angle at the lower boundary
!     of each atmospheric layer.
SUBROUTINE rttov_locpat_ad(opts, dosolar, profiles, profiles_ad, &
      profiles_dry, profiles_dry_ad, aux, coef, &
      angles, raytracing, raytracing_ad, do_pmc_calc)
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
                          geometry_type, raytracing_type
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
  TYPE(profile_type   ), INTENT(IN)    :: profiles(:)
  TYPE(profile_type   ), INTENT(INOUT) :: profiles_ad(:)
  TYPE(profile_type   ), INTENT(IN)    :: profiles_dry(:)
  TYPE(profile_type   ), INTENT(INOUT) :: profiles_dry_ad(:)
  TYPE(profile_aux    ), INTENT(IN)    :: aux
  TYPE(geometry_type  ), INTENT(IN)    :: angles(:)
  TYPE(rttov_coef     ), INTENT(IN)    :: coef
  TYPE(raytracing_type), INTENT(IN)    :: raytracing
  TYPE(raytracing_type), INTENT(INOUT) :: raytracing_ad
  LOGICAL(jplm)        , INTENT(IN), OPTIONAL :: do_pmc_calc
!INTF_END

  INTEGER(kind=jpim) :: iprof, lev, lay
  INTEGER(kind=jpim) :: nlevels, nlayers
  REAL   (kind=jprb) :: rearth(size(profiles))
  REAL   (kind=jprb) :: gravl(SIZE(profiles))   ! gravity at latitude lat [m/s^2]
  REAL   (kind=jprb) :: gravh(SIZE(profiles)), gravh_r(SIZE(profiles))
  REAL   (kind=jprb) :: eta, beta ! coefficients of the international gravity formula
  REAL   (kind=jprb) :: rlh(SIZE(profiles))
  REAL   (kind=jprb) :: dflat
  REAL   (kind=jprb) :: fac, c

  REAL   (kind=jprb) :: pres_ad(profiles(1)%nlevels,SIZE(profiles))
  REAL   (kind=jprb) :: temp_ad(profiles(1)%nlevels,SIZE(profiles))
  REAL   (KIND=jprb) :: qwet(profiles(1)%nlevels)
  REAL   (KIND=jprb) :: qwet_ad(profiles(1)%nlevels)

  INTEGER(kind=jpim) :: nprofiles
  LOGICAL(KIND=jplm) :: do_pmc_calc1

  REAL   (kind=jprb) :: scale_eps

  REAL   (KIND=JPRB) :: ZHOOK_HANDLE
!-----end of header-------------------------------------------------------------
  IF (LHOOK) CALL DR_HOOK('RTTOV_LOCPAT_AD', 0_jpim, ZHOOK_HANDLE)
  
  nprofiles = SIZE(profiles)
  nlevels     = profiles(1)%nlevels
  nlayers     = nlevels - 1

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

  c = 1000.0_jprb * rgc / (100.0_jprb * mair)

!Start AD

!-------7.  Compute the STP co2 thickness in cm

  IF (coef%pmc_shift .AND. do_pmc_calc1) CALL calc_ssu_co2_thickness_ad()

  temp_ad = 0._jprb
  !water_ad = 0._jprb
  IF (opts%interpolation%LGRADP) pres_ad = 0._jprb

! Extra bits
  IF (opts%rt_ir%addaerosl .OR. opts%rt_ir%addclouds) THEN
    DO iprof = 1, nprofiles
      DO LAY = NLEVELS - 1, 1,  - 1
        RAYTRACING_AD%HGPL(LAY, IPROF)     = RAYTRACING_AD%HGPL(LAY, IPROF) + RAYTRACING_AD%LTICK(LAY, IPROF)
        RAYTRACING_AD%HGPL(LAY + 1, IPROF) = RAYTRACING_AD%HGPL(LAY + 1, IPROF) - RAYTRACING_AD%LTICK(LAY, IPROF)
      ENDDO
    ENDDO
  ENDIF

  DO iprof = 1, NPROFILES
    RAYTRACING_AD%HGPL(2, IPROF) = RAYTRACING_AD%HGPL(2, IPROF) + RAYTRACING_AD%HGPL(1, IPROF)
  ENDDO
!-------6.  Compute secant of the zenith angle at the lower boundary of each -
!       layer for the satellite-to-surface line of view and of each             |
!       layer for the sun-to-surface line of view (if possible and reqd.)       |

  IF(opts%rt_all%ADDREFRAC .OR. opts%rt_ir%pc%ADDPC) THEN
    raytracing_ad%r = 0._jprb ! moved from init_raytracing
  ELSE
    raytracing_ad%ppw = 0._jprb
  ENDIF

  IF (dosolar) THEN
    raytracing_ad%pathsat = raytracing_ad%pathsat + raytracing_ad%patheff
    raytracing_ad%pathsun = raytracing_ad%pathsun + raytracing_ad%patheff

    CALL calc_seczen_ad(&
      raytracing%ratoesun, raytracing%zasun, raytracing%pathsun,&
      raytracing_ad%ratoesun, raytracing_ad%zasun, raytracing_ad%pathsun,&
      do_sat = .FALSE._jplm) ! solar call
!  ELSE
!    raytracing_ad%pathsun(:,:)  = 0._jprb
  ENDIF

  CALL calc_seczen_ad(raytracing%ratoesat, raytracing%zasat, raytracing%pathsat,&
    raytracing_ad%ratoesat, raytracing_ad%zasat, raytracing_ad%pathsat, &
    do_sat = .TRUE._jplm) ! satellite call


!-------5. compute atmospheric refractive index (if reqd.)---------------------
! compute the refractivity index for dry air at temperature tempp and pressure
! as a function of wave-number using an updated version of edlen equation
! not yet corrected for CO2 contribution - used in calc_seczen
  IF(opts%rt_all%ADDREFRAC .OR. opts%rt_ir%pc%ADDPC) CALL calc_refractivity_ad()

!-------5. Compute height of pressure levels
  CALL calc_hgpl_ad()
!-------4.  set up the height of the surface pressure level---------------------
! height of nearest level below surf counted from bottom of profile grid
  DO iprof = 1, nprofiles
    raytracing_ad%hgpl(aux%s(iprof)%nearestlev_surf, iprof) = 0._jprb
  ENDDO

  DO iprof = 1, nprofiles
! first calculate density for dry air (ideal gas)

    ! calculate water vapour in units of ppmv wet
    IF (profiles(iprof)%gas_units == gas_unit_compatibility) THEN
      qwet(:) = profiles_dry(iprof)%q(:)
    ELSE
      qwet(:) = profiles_dry(iprof)%q(:) / (1._jprb + profiles_dry(iprof)%q(:) * 1e-6_jprb)
    ENDIF

    CALL DIVIDE(c * raytracing_ad%dmair_r(:, iprof), &
      profiles(iprof)%p(:) - raytracing%ppw(:,iprof) * scale_eps, &
      temp_ad(:,iprof), acc = .TRUE._jplm)

    CALL DIVIDE(raytracing%dmair_r(:, iprof) * &
      scale_eps * raytracing_ad%dmair_r(:, iprof), &
      profiles(iprof)%p(:) - raytracing%ppw(:,iprof) * scale_eps, &
      raytracing_ad%ppw(:,iprof), acc = .TRUE._jplm)

      qwet_ad(:) = &!qwet_ad(:) + &
        1.0e-6_jprb * profiles(iprof)%p(:) * raytracing_ad%ppw(:,iprof)

    IF (opts%interpolation%LGRADP) THEN
      CALL DIVIDE(-raytracing%dmair_r(:, iprof) * &
        raytracing_ad%dmair_r(:, iprof), &
        profiles(iprof)%p(:) - raytracing%ppw(:,iprof) * scale_eps, &
        pres_ad(:, iprof), acc = .TRUE._jplm)

      pres_ad(:,iprof) = pres_ad(:,iprof) + &
        1.0e-6_jprb * raytracing_ad%ppw(:,iprof) * qwet(:)

      profiles_ad(iprof)%p(:) = profiles_ad(iprof)%p(:) + &
        pres_ad(:,iprof)
    ENDIF

    qwet_ad(2) = qwet_ad(2) + qwet_ad(1)
    qwet_ad(1) = 0._jprb

    IF (profiles(iprof)%gas_units == gas_unit_compatibility) THEN
      profiles_dry_ad(iprof)%q(:) = profiles_dry_ad(iprof)%q(:) + qwet_ad(:)
    ELSE
      profiles_dry_ad(iprof)%q(:) = profiles_dry_ad(iprof)%q(:) + &
        qwet_ad(:) * qwet(:) * (1._jprb - 1e-6_jprb * qwet(:)) / profiles_dry(iprof)%q(:)
    ENDIF

    profiles_ad(iprof)%t(2) = profiles_ad(iprof)%t(2) + temp_ad(1,iprof)

    temp_ad(1,iprof) = 0._jprb

    profiles_ad(iprof)%t(:) = profiles_ad(iprof)%t(:) + temp_ad(:,iprof)
  ENDDO

  IF (LHOOK) CALL DR_HOOK('RTTOV_LOCPAT_AD', 1_jpim, ZHOOK_HANDLE)
CONTAINS 

  SUBROUTINE calc_refractivity_ad()

    REAL(kind=jprb) :: dry_pp_h((profiles(1)%nlevels))
    REAL(kind=jprb) :: t_norm((profiles(1)%nlevels))
    REAL(kind=jprb) :: dry_pp_h_ad((profiles(1)%nlevels))
    REAL(kind=jprb) :: U((profiles(1)%nlevels)), V_R((profiles(1)%nlevels))
    REAL(kind=jprb) :: DU_AD((profiles(1)%nlevels)), DV_AD((profiles(1)%nlevels))

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

! compute and store (multiplicative) correction factor to account for presence 
! of CO2 in the air 
! assume constant CO2 using figure from 2005 (376ppmv) unless co2 data supplied

!AD
! store reciprocal of refractive index for performance - also reverse array index

    DO iprof = 1, nprofiles
! calculate useful temporary array to break refractivity calculation up 
! in to a better size

     ! calculate partial pressure of water vapour
      dry_pp_h(:) = htop * (profiles(iprof)%p(:) - raytracing%ppw(:, iprof))
      t_norm(:) = profiles(iprof)%t(:) - t0

      U = 1._jprb + 1.0e-8_jprb * (ed2 - ed3 * t_norm(:)) * dry_pp_h(:)
      CALL reciprocal(1._jprb + ed4 * t_norm(:), V_R)

!AD
! compute moist contribution to refractivity
! and store final refractive index for atmospheric profile
      raytracing_ad%ppw(:, iprof) = -&! raytracing_ad%ppw(:, iprof) - &!
        (raytracing_ad%r(:, iprof) * DISP_MOIST)

      ! apply CO2 correction
      IF (opts%rt_ir%co2_data) THEN
        raytracing_ad%dispco2(:, iprof) = &! raytracing_ad%dispco2(:, iprof) + &
             raytracing%refractivity(:, iprof) * raytracing_ad%r(:, iprof)

        raytracing_ad%refractivity(:, iprof) = &!raytracing_ad%refractivity(:, iprof) + &
          raytracing%dispco2(:, iprof) * raytracing_ad%r(:, iprof) 
      ELSE
        DISPCO2 = DISP * (1._jprb + dco2 * (CO2_CONC * 1e-6_jprb - 0.0003_jprb))
        raytracing_ad%refractivity(:, iprof) = &! raytracing_ad%refractivity(:, iprof) +
          raytracing_ad%r(:, iprof) * DISPCO2 ! AD
      ENDIF

      dry_pp_h_ad(:) = &!dry_pp_h_ad(:) + &
        raytracing_ad%refractivity(:, iprof) * (U * V_R)/ed1

      DU_ad(:) = &!DU_ad(:) + &
        dry_pp_h(:) * raytracing_ad%refractivity(:, iprof) * V_R(:) / ed1

      DV_ad(:) = &!DV_ad(:) &
        -raytracing%refractivity(:, iprof) * raytracing_ad%refractivity(:, iprof) * V_R(:)

      temp_ad(:,iprof) = temp_ad(:,iprof) + &
        ed4 * DV_ad(:) + & !combine two lines from TL     
!      temp_ad(:,iprof) = temp_ad(:,iprof) + &
        1.0e-8_jprb * dry_pp_h(:) * (-ed3) * DU_ad(:)

      dry_pp_h_ad(:) = dry_pp_h_ad(:) + &
        1.0e-8_jprb * (ed2 - ed3 * t_norm(:)) * DU_ad(:)

      raytracing_ad%ppw(:, iprof) = raytracing_ad%ppw(:, iprof) - &
        htop * dry_pp_h_ad(:)
     
      IF(opts%interpolation%LGRADP) THEN
        pres_ad(:,iprof) = pres_ad(:,iprof) + &
          htop * dry_pp_h_ad(:)

      ENDIF
    ENDDO

    IF (opts%rt_ir%co2_data) THEN
      DO iprof = 1, nprofiles
        profiles_dry_ad(iprof)%co2(:) = profiles_dry_ad(iprof)%co2(:) + &
          raytracing_ad%dispco2(:, iprof) * disp * &
          dco2 * 1e-6_jprb
      ENDDO
    ELSE
     !     DO iprof = 1, nprofiles
!       profiles_dry_ad(iprof)%co2(:) = 0._jprb
!     ENDDO
    ENDIF

 END SUBROUTINE calc_refractivity_ad

 SUBROUTINE calc_hgpl_ad()
   REAL(kind=jprb) :: dp((profiles(1)%nlevels))
   REAL(kind=jprb) :: dp_ad((profiles(1)%nlevels),SIZE(profiles))
   REAL(kind=jprb) :: ztemp((profiles(1)%nlevels),SIZE(profiles))

   CALL INVSQRT(raytracing%ztemp(:, :), ztemp(:, :))

   DO iprof = 1, nprofiles
!     dp(1:nlevels-1) = 0.5_jprb * &
!       (profiles(iprof)%p(nlevels-1:1:-1) - profiles(iprof)%p(nlevels:2:-1))
     DO lev = 1, nlevels - 1
       dp(lev) = 0.5_jprb * (profiles(iprof)%p(lev) - &
                            profiles(iprof)%p(lev + 1))
     ENDDO

     DO lev = nlevels, aux%s(iprof)%nearestlev_surf + 1, -1

       raytracing_ad%ztemp(lev, iprof) = &!raytracing_ad%ztemp(lev, iprof) -
         -0.5_jprb * ztemp(lev,iprof) * 1e-3_jprb * raytracing_ad%hgpl(lev, iprof)

       raytracing_ad%INT(lev, iprof) = &!raytracing_ad%INT(lev, iprof) - &
         -2._jprb * 100._jprb * gravh_r(iprof) * raytracing_ad%ztemp(lev, iprof)

       raytracing_ad%hgpl(lev - 1, iprof) = raytracing_ad%hgpl(lev - 1, iprof) + 1.0e3_jprb * &
         2._jprb * raytracing_ad%ztemp(lev, iprof) * &
         (1.0e3_jprb * raytracing%hgpl(lev - 1, iprof) - rlh(iprof))

       raytracing_ad%dmair_r(lev, iprof) = &
         raytracing_ad%dmair_r(lev, iprof) + &
         raytracing_ad%INT(lev, iprof) * dp(lev - 1)

       raytracing_ad%dmair_r(lev - 1, iprof) = &
         raytracing_ad%dmair_r(lev - 1, iprof) + &
         raytracing_ad%INT(lev, iprof) * dp(lev - 1)

       IF(opts%interpolation%lgradp) THEN
         dp_ad(lev - 1, iprof) = &!dp_ad(lev - 1, iprof) + &
           raytracing_ad%INT(lev, iprof) * &
           (raytracing%dmair_r(lev - 1, iprof) + &
           raytracing%dmair_r(lev, iprof))
       ENDIF
     ENDDO

     DO lev = 1, aux%s(iprof)%nearestlev_surf - 1
       raytracing_ad%ztemp(lev, iprof) = &!raytracing_ad%ztemp(lev, iprof) -
         -0.5_jprb * ztemp(lev,iprof) * 1e-3_jprb * raytracing_ad%hgpl(lev, iprof)

       raytracing_ad%INT(lev, iprof) = &!raytracing_ad%INT(lev, iprof) - &
         -2._jprb * 100._jprb * gravh_r(iprof) * raytracing_ad%ztemp(lev, iprof)

       raytracing_ad%hgpl(lev + 1, iprof) = raytracing_ad%hgpl(lev + 1, iprof) + 1.0e3_jprb * &
         2._jprb * raytracing_ad%ztemp(lev, iprof) * &
         (1.0e3_jprb * raytracing%hgpl(lev + 1, iprof) - rlh(iprof))

       raytracing_ad%dmair_r(lev + 1, iprof) = &
         raytracing_ad%dmair_r(lev + 1, iprof) - &
         raytracing_ad%INT(lev, iprof) * dp(lev)

       raytracing_ad%dmair_r(lev, iprof) = &
         raytracing_ad%dmair_r(lev, iprof) - &
         raytracing_ad%INT(lev, iprof) * dp(lev)

       IF(opts%interpolation%lgradp) THEN
         dp_ad(lev, iprof) = &!dp_ad(lev, iprof) - &
           -raytracing_ad%INT(lev, iprof) * &
           (raytracing%dmair_r(lev, iprof) + &
            raytracing%dmair_r(lev + 1, iprof))
       ENDIF
     ENDDO

     IF(opts%interpolation%LGRADP) THEN
       DO lev = nlevels - 1,1, -1
         pres_ad(lev, iprof) = pres_ad(lev, iprof) + 0.5_jprb * dp_ad(lev, iprof)
         pres_ad(lev + 1, iprof) = pres_ad(lev + 1, iprof) - 0.5_jprb * dp_ad(lev, iprof)
       ENDDO
     ENDIF

   ENDDO
 END SUBROUTINE calc_hgpl_ad

!-------Compute secant of the zenith angle at the lower boundary of each -
!       layer for the sat/sun-to-surface line of view.                       |
  SUBROUTINE calc_seczen_ad(ra,za,path,ra_ad,za_ad,path_ad,do_sat)
    LOGICAL(jplm), INTENT(IN) :: do_sat ! do satellite secant or solar secant
    REAL(jprb), INTENT(IN) :: ra(:,:), za(:,:), path(:,:)
    REAL(jprb), INTENT(INOUT) :: ra_ad(:,:), za_ad(:,:), path_ad(:,:)

    REAL(jprb) :: z(SIZE(ra(:,1)),SIZE(ra(1,:)))
    INTEGER(jpim) :: n

    IF(do_sat .AND. aux%on_coef_levels) THEN
      raytracing_ad%pathsat = raytracing_ad%pathsat + &
        raytracing_ad%pathsat_sqrt * raytracing%pathsat_rsqrt

      raytracing_ad%pathsat_rsqrt = &!raytracing_ad%pathsat_rsqrt + &
        raytracing%pathsat * raytracing_ad%pathsat_sqrt

      CALL INVSQRT_AD(raytracing%pathsat_rsqrt, raytracing_ad%pathsat, raytracing_ad%pathsat_rsqrt, acc=.TRUE._jplm)
    ENDIF

! calculate path from zenith angle using path = 1/sqrt(1-za**2)
    IF(opts%rt_ir%addsolar) THEN
      CALL SINTOSEC_AD(za(1:nlayers,:), za_ad(1:nlayers,:), path, path_ad, acc = .TRUE._jplm)
    ELSE
      CALL SINTOSEC_AD(za(1:nlayers,:), za_ad(1:nlayers,:), path, path_ad, acc = .FALSE._jplm) ! first use of za
    ENDIF

!   path_ad = 0._jprb; 

    DO iprof = 1, nprofiles
      IF(do_sat) THEN
        n = nlayers - 1

        z(:, iprof) = (rearth(iprof) + coef%fc_sat_height) * &
                      raytracing%z_r(:,iprof)

        ra_ad(:,iprof) = &!ra_ad(:,iprof) + 
                         za_ad(1:nlayers,iprof) * angles(iprof)%sinview
      ELSE ! do_sun
        IF(profiles(iprof)%sunzenangle >= 0.0 .AND. &
          profiles(iprof)%sunzenangle < max_sol_zen) THEN
          n = 0
          z(:,iprof) = rearth(iprof) * raytracing%z_r(:,iprof)
          ra_ad(:,iprof) = &!ra_ad(:,iprof) + 
                           za_ad(1:nlayers,iprof) * angles(iprof)%sinzen_sun
        ELSE
!          za_ad(:,iprof) = 0.0_jprb
        ENDIF
      ENDIF
    ENDDO
!  za_ad = 0._jprb;

    IF(opts%rt_all%ADDREFRAC .OR. opts%rt_ir%pc%ADDPC) THEN
      DO iprof = 1, nprofiles
! the (nlevels - n) takes in to account the fact that 10.2 and earlier
! had the refraction calculation on unpacked, reversed levels.
        raytracing_ad%r(2:nlevels,iprof) = raytracing_ad%r(2:nlevels,iprof) - & !first use
          z(:,iprof) * raytracing%r(nlevels - n, iprof) * raytracing%r_r(2:nlevels,iprof)**2_jpim * ra_ad(:,iprof)

! DAR - this is how you replace a loop with a vector in the AD
        raytracing_ad%r(nlevels - n,iprof) = raytracing_ad%r(nlevels - n, iprof) + &
          SUM(raytracing%r_r(2:nlevels,iprof) * z(:,iprof) * ra_ad(:,iprof)) 

        ra_ad(:,iprof) = &
          raytracing%r_r(2:nlevels,iprof) * raytracing%r(nlevels - n, iprof) * &
          ra_ad(:,iprof)
      ENDDO
    ENDIF

    IF (do_sat) THEN
      DO iprof = 1, nprofiles
        raytracing_ad%z_r(:,iprof) = &!raytracing_ad%z_r(:,iprof) + &
          (rearth(iprof) + coef%fc_sat_height) * ra_ad(:,iprof)
      ENDDO
    ELSE ! do sun
      DO iprof = 1, nprofiles
        IF(profiles(iprof)%sunzenangle >= 0.0 .AND. &
          profiles(iprof)%sunzenangle < max_sol_zen) THEN
          raytracing_ad%z_r(1:nlayers,iprof) = &!raytracing_ad%z_r(1:nlayers,iprof) + & 
            rearth(iprof) * ra_ad(1:nlayers,iprof)
        ENDIF
      ENDDO
    ENDIF
!   ra_ad = 0._jprb;

    CALL reciprocal_ad(raytracing%z_r,raytracing_ad%hgpl(2:nlayers+1, :), &
         raytracing_ad%z_r, acc = .TRUE._jplm)
  END SUBROUTINE CALC_SECZEN_AD
  !------------------------------------------------------------------------------
!-------11  Compute the STP co2 thickess in cm
!------------------------------------------------------------------------------
! loops go top to bottom
      
  SUBROUTINE calc_ssu_co2_thickness_ad
    REAL(kind=jprb) :: pstd
    REAL(kind=jprb) :: tstd
    REAL(kind=jprb) :: pa(profiles(1)%nlevels-1),   pa_ad(profiles(1)%nlevels-1)
    REAL(kind=jprb) :: dpa(profiles(1)%nlevels-1),  dpa_ad(profiles(1)%nlevels-1)
    REAL(kind=jprb) :: rho(profiles(1)%nlevels-1),  rho_ad(profiles(1)%nlevels-1)
    REAL(kind=jprb) :: temp(profiles(1)%nlevels-1), temp_ad(profiles(1)%nlevels-1)
    REAL(kind=jprb) :: co2(profiles(1)%nlevels-1),  co2_ad(profiles(1)%nlevels - 1)
    REAL(kind=jprb) :: dmair(profiles(1)%nlevels),  dmair_ad(profiles(1)%nlevels)
    REAL(kind=jprb) :: dz, dz_ad
    REAL(kind=jprb) :: dzstp_ad

    pstd=p0*100._jprb
    tstd=t0
    
    gravh(:) = 1._jprb / gravh_r(:)

    DO iprof = 1, nprofiles
      CALL reciprocal(raytracing%dmair_r(:,iprof), dmair(:))
      dmair_ad(:) = 0._jprb
      DO lay = 1, nlayers
! fwd
        pa(lay)   = 100._jprb * (profiles(iprof)%p(lay) + profiles(iprof)%p(lay+1)) / 2._jprb
        dpa(lay)  = 100._jprb * (profiles(iprof)%p(lay) - profiles(iprof)%p(lay+1))
        rho(lay)  = (dmair(lay) + dmair(lay+1)) / 2._jprb
        temp(lay) = (profiles(iprof)%t(lay) + profiles(iprof)%t(lay+1)) / 2._jprb
        co2(lay)  = (profiles_dry(iprof)%co2(lay) + profiles_dry(iprof)%co2(lay+1)) / 2._jprb
! ad
        pa_ad(lay)   = 0._jprb
        dpa_ad(lay)  = 0._jprb
        rho_ad(lay)  = 0._jprb
        temp_ad(lay) = 0._jprb
        co2_ad(lay)  = 0._jprb
      ENDDO
    
      DO lay = 1, nlayers
! fwd
        dz = -dpa(lay)/&
          ( (gravl(iprof)-gravh(iprof)*raytracing%hgpl(lay+1, iprof)) * rho(lay) )
! ad
        dzstp_ad = 0._jprb
        dzstp_ad = dzstp_ad + raytracing_ad%co2_cm(iprof)
        dz_ad = 0._jprb
        dz_ad =  dz_ad + &
          dzstp_ad * pa(lay)/pstd * tstd/temp(lay)* co2(lay)*1e-6_jprb*100._jprb

        IF (opts%interpolation%lgradp) THEN
          pa_ad(lay) = pa_ad(lay) + &
            dzstp_ad * dz / pstd * tstd/temp(lay)* co2(lay)*1e-6_jprb*100._jprb
!        ELSE
!          pa_ad(lay) = 0._jprb
        ENDIF
        
        temp_ad(lay) = temp_ad(lay) - &
          dzstp_ad * dz * pa(lay)/pstd * tstd/temp(lay)**2_jpim * &
          co2(lay)*1e-6_jprb*100._jprb

        co2_ad(lay) = co2_ad(lay) + &
          dzstp_ad * dz * pa(lay)/pstd * tstd/temp(lay)*1e-6_jprb*100._jprb
!        dzstp_ad = 0._jprb
        
        IF (opts%interpolation%lgradp) THEN
          dpa_ad(lay) = dpa_ad(lay) - dz_ad / ( (gravl(iprof)-gravh(iprof)*raytracing%hgpl(lay+1, iprof)) * rho(lay) )
!        ELSE
!          dpa_ad(lay) = 0._jprb
        ENDIF

        raytracing_ad%hgpl(lay+1, iprof) = raytracing_ad%hgpl(lay+1, iprof) &
          - dz_ad * dpa(lay) * gravh(iprof)*rho(lay)/ &
          ( (gravl(iprof)-gravh(iprof)*raytracing%hgpl(lay+1, iprof)) * rho(lay) )**2_jpim

        rho_ad(lay) = rho_ad(lay)                                              &
          + dz_ad * dpa(lay)*(gravl(iprof)-gravh(iprof)*raytracing%hgpl(lay+1, iprof)) / &
          ( (gravl(iprof)-gravh(iprof)*raytracing%hgpl(lay+1, iprof)) * rho(lay) )**2_jpim
        
!        dz_ad = 0._jprb
        
        profiles_dry_ad(iprof)%co2(lay) = profiles_dry_ad(iprof)%co2(lay) + &
          co2_ad(lay)/2._jprb
        profiles_dry_ad(iprof)%co2(lay+1) = profiles_dry_ad(iprof)%co2(lay+1) + &
          co2_ad(lay)/2._jprb
!        co2_ad(lay) = 0._jprb

        profiles_ad(iprof)%t(lay) = profiles_ad(iprof)%t(lay) + &
          temp_ad(lay)/2._jprb

        profiles_ad(iprof)%t(lay+1) = profiles_ad(iprof)%t(lay+1) + &
          temp_ad(lay)/2._jprb
!        temp_ad(lay) = 0._jprb
        
        dmair_ad(lay) = dmair_ad(lay) + rho_ad(lay) / 2._jprb
        dmair_ad(lay+1) = dmair_ad(lay+1) + rho_ad(lay) / 2._jprb
       
!        rho_ad(lay) = 0._jprb
        
        IF (opts%interpolation%lgradp) THEN
          profiles_ad(iprof)%p(lay) = profiles_ad(iprof)%p(lay) + &
            dpa_ad(lay)*100._jprb
          profiles_ad(iprof)%p(lay+1) = profiles_ad(iprof)%p(lay+1) - &
            dpa_ad(lay)*100._jprb
!          dpa_ad(lay) = 0._jprb
          profiles_ad(iprof)%p(lay) = profiles_ad(iprof)%p(lay) + &
            pa_ad(lay)*100._jprb/2._jprb
          profiles_ad(iprof)%p(lay+1) = profiles_ad(iprof)%p(lay+1) + &
            pa_ad(lay)*100._jprb/2._jprb
!          pa_ad(lay) = 0._jprb
!        ELSE
!          profiles_ad(iprof)%p(lay) = 0._jprb
!          profiles_ad(iprof)%p(lay+1) = 0._jprb
!          dpa_ad(lay) = 0._jprb
!          pa_ad(lay) = 0._jprb
        ENDIF
      ENDDO
      CALL reciprocal_ad(dmair(:), raytracing_ad%dmair_r(:,iprof), &
                         dmair_ad(:), acc = .FALSE._jplm)
    END DO

  END SUBROUTINE calc_ssu_co2_thickness_ad

END SUBROUTINE rttov_locpat_ad
