!     Calculate the secant of the zenith angle at the lower boundary
!     of each atmospheric layer.
SUBROUTINE rttov_locpat_tl(opts, dosolar, profiles, profiles_tl, &
      profiles_dry, profiles_dry_tl, aux, coef, &
      angles, raytracing, raytracing_tl, do_pmc_calc)
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
  USE rttov_const, ONLY : d1, d2, d3, d4, d5, dco2, ed1, ed2, ed3, ed4,&
                          ew1, ew2, htop, ctom, waver, rgc, mair, mh2o, flatt,  &
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
  TYPE(profile_type   ), INTENT(IN)    :: profiles_tl(:)
  TYPE(profile_type   ), INTENT(IN)    :: profiles_dry(:)
  TYPE(profile_type   ), INTENT(IN)    :: profiles_dry_tl(:)
  TYPE(profile_aux    ), INTENT(IN)    :: aux
  TYPE(geometry_type  ), INTENT(IN)    :: angles(:)
  TYPE(rttov_coef     ), INTENT(IN)    :: coef
  TYPE(raytracing_type), INTENT(IN)    :: raytracing
  TYPE(raytracing_type), INTENT(INOUT) :: raytracing_tl
  LOGICAL(jplm)        , INTENT(IN), OPTIONAL :: do_pmc_calc
!INTF_END

  INTEGER(kind=jpim) :: iprof, lev, lay
  INTEGER(kind=jpim) :: nlevels , nlayers
  REAL   (kind=jprb) :: rearth(size(profiles))
  REAL   (kind=jprb) :: gravl(SIZE(profiles))   ! gravity at latitude lat [m/s^2]
  REAL   (kind=jprb) :: gravh(SIZE(profiles)), gravh_r(SIZE(profiles))
  REAL   (kind=jprb) :: eta, beta ! coefficients of the international gravity formula
  REAL   (kind=jprb) :: rlh(SIZE(profiles))

  REAL   (kind=jprb) :: dflat
  REAL   (kind=jprb) :: fac, c

  REAL   (kind=jprb) :: pres_tl(profiles(1)%nlevels,SIZE(profiles))
  REAL   (kind=jprb) :: temp_tl(profiles(1)%nlevels,SIZE(profiles))
  REAL   (KIND=jprb) :: qwet(profiles(1)%nlevels)
  REAL   (KIND=jprb) :: qwet_tl(profiles(1)%nlevels)

  INTEGER(kind=jpim) :: nprofiles
  LOGICAL(KIND=jplm) :: do_pmc_calc1

  REAL   (kind=jprb) :: scale_eps

  REAL   (KIND=JPRB) :: ZHOOK_HANDLE

!-----end of header-------------------------------------------------------------
  IF (LHOOK) CALL DR_HOOK('RTTOV_LOCPAT_TL', 0_jpim, ZHOOK_HANDLE)
  
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

!-------3. calculate reciprocal of moist air density (for later calculation)
! density of dry air is adjusted to account for the presence of water vapour 
! by replacing temp with virtual temp

  c = 1000.0_jprb * rgc / (100.0_jprb * mair)

  DO iprof = 1, nprofiles
! first calculate density for dry air (ideal gas)

    temp_tl(:,iprof) = profiles_tl(iprof)%t(:)
    temp_tl(1,iprof) = profiles_tl(iprof)%t(2)

    ! calculate water vapour in units of ppmv wet
    IF (profiles(iprof)%gas_units == gas_unit_compatibility) THEN
      qwet(:) = profiles_dry(iprof)%q(:)
      qwet_tl(:) = profiles_dry_tl(iprof)%q(:)
    ELSE
      qwet(:) = profiles_dry(iprof)%q(:) / (1._jprb + profiles_dry(iprof)%q(:) * 1e-6_jprb)
      qwet_tl(:) = profiles_dry_tl(iprof)%q(:) * qwet(:) * &
                         (1._jprb - 1e-6_jprb * qwet(:)) / profiles_dry(iprof)%q(:)
    ENDIF
    qwet_tl(1) = qwet_tl(2) ! As for temp_tl above

    IF (opts%interpolation%LGRADP) THEN
      pres_tl(:,iprof) = profiles_tl(iprof)%p(:)
      raytracing_tl%ppw(:,iprof) = 1e-6_jprb * &
           (pres_tl(:,iprof) * qwet(:) + &
            profiles(iprof)%p(:) * qwet_tl(:))

      CALL DIVIDE(c * temp_tl(:,iprof) - raytracing%dmair_r(:, iprof) * &
           (pres_tl(:,iprof) - scale_eps * raytracing_tl%ppw(:,iprof)), &
           profiles(iprof)%p(:) - raytracing%ppw(:,iprof) * scale_eps, &
           raytracing_tl%dmair_r(:, iprof), acc = .FALSE._jplm)

    ELSE
      raytracing_tl%ppw(:,iprof) = 1e-6_jprb * &
           profiles(iprof)%p(:) * qwet_tl(:)

      CALL DIVIDE(c * temp_tl(:,iprof) + raytracing%dmair_r(:, iprof) * &
           scale_eps * raytracing_tl%ppw(:,iprof), &
           profiles(iprof)%p(:) - raytracing%ppw(:,iprof) * scale_eps, &
           raytracing_tl%dmair_r(:, iprof), acc = .FALSE._jplm)
    ENDIF
  ENDDO

!-------4.  set up the height of the surface pressure level---------------------
! height of nearest level below surf counted from bottom of profile grid
  DO iprof = 1, nprofiles
    raytracing_tl%hgpl(aux%s(iprof)%nearestlev_surf, iprof) = 0._jprb
  ENDDO

!-------5. Compute height of pressure levels
  CALL calc_hgpl_tl()

!-------5. compute atmospheric refractive index (if reqd.)---------------------
! compute the refractivity index for dry air at temperature tempp and pressure
! as a function of wave-number using an updated version of edlen equation
! not yet corrected for CO2 contribution - used in calc_seczen
  IF(opts%rt_all%ADDREFRAC .OR. opts%rt_ir%pc%ADDPC) CALL calc_refractivity_tl()

!-------6.  Compute secant of the zenith angle at the lower boundary of each -
!       layer for the satellite-to-surface line of view and of each             |
!       layer for the sun-to-surface line of view (if possible and reqd.)       |
  CALL calc_seczen_tl(&
    raytracing%ratoesat, raytracing%zasat, raytracing%pathsat, &
    raytracing_tl%ratoesat, raytracing_tl%zasat, raytracing_tl%pathsat, &
    do_sat=.TRUE._jplm) ! satellite call

  IF (dosolar) THEN
    CALL calc_seczen_tl(&
      raytracing%ratoesun, raytracing%zasun, raytracing%pathsun,&
      raytracing_tl%ratoesun, raytracing_tl%zasun, raytracing_tl%pathsun,&
      do_sat=.FALSE._jplm) ! solar call
    raytracing_tl%patheff = raytracing_tl%pathsat + raytracing_tl%pathsun
  ENDIF

! Extra bits
  DO iprof = 1, nprofiles
    raytracing_tl%hgpl(1, iprof) = raytracing_tl%hgpl(2, iprof)
    IF (opts%rt_ir%addaerosl .OR. opts%rt_ir%addclouds) THEN
      DO lay = 1, nlayers ! layer thickness
        raytracing_tl%ltick(lay, iprof) = raytracing_tl%hgpl(lay, iprof) - &
          raytracing_tl%hgpl(lay + 1, iprof)
      ENDDO
    ENDIF
  ENDDO

!-------7.  Compute the STP co2 thickness in cm
  IF (coef%pmc_shift .AND. do_pmc_calc1) CALL calc_ssu_co2_thickness_tl()

  IF (LHOOK) CALL DR_HOOK('RTTOV_LOCPAT_TL', 1_jpim, ZHOOK_HANDLE)
CONTAINS 

  SUBROUTINE calc_refractivity_tl()

    REAL(kind=jprb) :: dry_pp_h((profiles(1)%nlevels))
    REAL(kind=jprb) :: t_norm((profiles(1)%nlevels))
    REAL(kind=jprb) :: dry_pp_h_tl((profiles(1)%nlevels))
    REAL(kind=jprb) :: U((profiles(1)%nlevels)), V_R((profiles(1)%nlevels))
    REAL(kind=jprb) :: DU((profiles(1)%nlevels)), DV((profiles(1)%nlevels))

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

    IF (opts%rt_ir%co2_data) THEN
      DO iprof = 1, nprofiles
        raytracing_tl%dispco2(:, iprof) = disp * &
          dco2 * profiles_dry_tl(iprof)%co2(:) * 1e-6_jprb
      ENDDO
    ELSE
      DISPCO2 = DISP * (1._jprb + dco2 * (CO2_CONC * 1e-6_jprb - 0.0003_jprb))
    ENDIF

    DO iprof = 1, nprofiles
! calculate useful temporary array to break refractivity calculation up 
! in to a better size
      dry_pp_h(:) = htop * (profiles(iprof)%p(:) - raytracing%ppw(:, iprof))
      t_norm(:) = profiles(iprof)%t(:) - t0
      U = 1._jprb + 1.0e-8_jprb * (ed2 - ed3 * t_norm(:)) * dry_pp_h(:)
      CALL reciprocal(1._jprb + ed4 * t_norm(:), V_R)

      IF(opts%interpolation%LGRADP) THEN
        dry_pp_h_tl(:) = htop * &
          (pres_tl(:,iprof) - raytracing_tl%ppw(:, iprof))
      ELSE
        dry_pp_h_tl(:) = -htop * raytracing_tl%ppw(:, iprof)
      ENDIF

      DU = 1.0e-8_jprb * &
           (dry_pp_h(:) * (-ed3) * temp_tl(:,iprof) + &
           (ed2 - ed3 * t_norm(:)) * dry_pp_h_tl(:))

      DV = ed4 * temp_tl(:,iprof)
       
      raytracing_tl%refractivity(:, iprof) = dry_pp_h_tl(:) * (U * V_R)/ed1 + &
        (dry_pp_h(:) * DU/ed1 - raytracing%refractivity(:, iprof) * DV) * V_R

! apply CO2 correction then
! compute moist contribution to refractivity
! and store final refractive index for atmospheric profile
      
      IF (opts%rt_ir%co2_data) THEN
        raytracing_tl%r(:, iprof) = & !refractivity back to refractance
          raytracing_tl%refractivity(:, iprof) * raytracing%dispco2(:, iprof) +& 
          (raytracing%refractivity(:, iprof)) * raytracing_tl%dispco2(:, iprof) - &
          raytracing_tl%ppw(:, iprof) * DISP_MOIST
      ELSE
        raytracing_tl%r(:, iprof) = raytracing_tl%refractivity(:, iprof) * DISPCO2 - &
          raytracing_tl%ppw(:, iprof) * DISP_MOIST
      ENDIF
    ENDDO
! store reciprocal of refractive index for performance - also reverse array index
!    CALL reciprocal_tl(raytracing%r_r, raytracing_tl%r, raytracing_tl%r_r)
  END SUBROUTINE calc_refractivity_tl

  SUBROUTINE calc_hgpl_tl()
    REAL(kind=jprb) :: dp((profiles(1)%nlevels),size(profiles))
    REAL(kind=jprb) :: dp_tl((profiles(1)%nlevels),size(profiles))
    REAL(kind=jprb) :: ztemp((profiles(1)%nlevels),size(profiles))

!-------Compute height of pressure levels:levels above surface. ------------
!       The height of pressure levels H is obtained by integrating             |
!       the hydrostatic equation dPRES=-GRAVL(H)*DMAIR*dH between              |
!       two adjacent pressure levels                                           |
!                                                                              |
!           -P2         - H2                                                   |
!          |           |                                                       |
!          | dP/DMAIR= | GRAVL(H)*dH                                           |
!          |           |                                                       |
!         -  P1       -   H1                                                   |
!                                                                              |
!       The integration id DMAIR is carried out assuming this                  |
!       quantity varies linearly with pressure. The integral is                |
!       then computed applying the trapezoidal formula.                        |
!-----------------------------------------------------------------------------

! unpack input pressure profiles
    DO iprof = 1, nprofiles
!      dp(1:nlevels-1, iprof) = 0.5_jprb * (profiles(iprof)%p(nlevels-1:1:-1)  - &
!                             profiles(iprof)%p(nlevels:2:-1))

      DO lev = 1, nlevels - 1
        dp(lev, iprof) = 0.5_jprb * (profiles(iprof)%p(lev) - &
                             profiles(iprof)%p(lev + 1))
      ENDDO
    ENDDO

    IF(opts%interpolation%LGRADP) THEN
!      dp_tl(1:nlevels-1, :) = 0.5_jprb * (pres_tl(nlevels-1:1:-1, :) - pres_tl(nlevels:2:-1, :))
      DO lev = 1, nlevels - 1
        dp_tl(lev, :) = 0.5_jprb * (pres_tl(lev, :) - pres_tl(lev + 1, :))
      ENDDO
    ENDIF

! Store ztemp from direct temporarily.
    CALL INVSQRT(raytracing%ztemp(:, :), ztemp(:, :))

    DO iprof = 1, nprofiles

      DO lev = aux%s(iprof)%nearestlev_surf - 1, 1, -1
          IF(opts%interpolation%LGRADP) THEN
            raytracing_tl%INT(lev, iprof) = &
              -dp(lev,iprof) * &
              (raytracing_tl%dmair_r(lev, iprof) + &
              raytracing_tl%dmair_r(lev + 1, iprof)) - &
              dp_tl(lev,iprof) * &
              (raytracing%dmair_r(lev, iprof) + &
              raytracing%dmair_r(lev + 1, iprof))
          ELSE
            raytracing_tl%INT(lev, iprof)  = &
              -dp(lev,iprof) *      &
              (raytracing_tl%dmair_r(lev, iprof) + &
              raytracing_tl%dmair_r(lev + 1, iprof))
          ENDIF

          raytracing_tl%ztemp(lev, iprof) = &
            2._jprb * 1.0e3_jprb * raytracing_tl%hgpl(lev + 1, iprof) * &
            (1.0e3_jprb * raytracing%hgpl(lev + 1, iprof) - rlh(iprof)) - &
            2._jprb * 100._jprb * gravh_r(iprof) * raytracing_tl%INT(lev, iprof)

          RAYTRACING_TL%HGPL(lev, iprof) = &
            -0.5_jprb * ztemp(lev,iprof) * 1e-3_jprb * raytracing_tl%ztemp(lev, iprof)
      ENDDO

      DO lev = aux%s(iprof)%nearestlev_surf + 1, nlevels
        IF(opts%interpolation%LGRADP) THEN
          raytracing_tl%INT(lev, iprof) = &
            dp_tl(lev-1,iprof) * &
            (raytracing%dmair_r(lev - 1, iprof) + &
            raytracing%dmair_r(lev, iprof)) + &
            dp(lev-1, iprof) *      &
            (raytracing_tl%dmair_r(lev - 1, iprof) + &
            raytracing_tl%dmair_r(lev, iprof))
        ELSE
          raytracing_tl%INT(lev, iprof)  = &
            dp(lev-1,iprof) *      &
            (raytracing_tl%dmair_r(lev - 1, iprof) + &
            raytracing_tl%dmair_r(lev, iprof))
        ENDIF

        raytracing_tl%ztemp(lev, iprof) = &
          2._jprb * 1.0e3_jprb * raytracing_tl%hgpl(lev - 1, iprof) * &
          (1.0e3_jprb * raytracing%hgpl(lev - 1, iprof) - rlh(iprof)) - &
          2._jprb * 100._jprb * gravh_r(iprof) * raytracing_tl%INT(lev, iprof)

        RAYTRACING_TL%HGPL(lev, iprof) = &
          -0.5_jprb * ztemp(lev,iprof) * 1e-3_jprb * raytracing_tl%ztemp(lev, iprof)
      ENDDO
    ENDDO

  END SUBROUTINE calc_hgpl_tl

!-------Compute secant of the zenith angle at the lower boundary of each -
!       layer for the sat/sun-to-surface line of view.                       |
  SUBROUTINE calc_seczen_tl(ra,za,path,ra_tl,za_tl,path_tl,do_sat)
!       The atmospheric layers are considered as concentric rings. If we trace |
!       a ray across these rings at any angle other than nadir, the local      |
!       angle relative to the outward radial direction at the point of         |
!       intersection will be different at each ring because due to the         |
!       curvature of the Earth and to atmospheric refraction. The secant of    |
!       the local/solar zenith angle PATHSAT/SUN is thus computed taking into  |
!       account the geometry of the situation and the bending of rays as they  |
!       traverse the atmosphere (by application of Snell's law).               |
!------------------------------------------------------------------------------

    LOGICAL(jplm), INTENT(IN) :: do_sat ! do satellite secant or solar secant
    REAL(jprb), INTENT(IN) :: ra(:,:), za(:,:), path(:,:)
    REAL(jprb), INTENT(OUT) :: ra_tl(:,:), za_tl(:,:), path_tl(:,:)
    
    REAL(jprb) :: z(SIZE(ra(:,1)),SIZE(ra(1,:)))
    INTEGER(jpim) :: n
    
    CALL reciprocal_tl(raytracing%z_r,raytracing_tl%hgpl(2:nlayers+1, :),&
         raytracing_tl%z_r)

    IF (do_sat) THEN
      n = nlayers - 1
      DO iprof = 1, nprofiles
        ! z is temporary array because we overwrote ra
        z(:, iprof) = (rearth(iprof) + coef%fc_sat_height) * &
                      raytracing%z_r(:,iprof)
        ra_tl(:,iprof) = (rearth(iprof) + coef%fc_sat_height) * &
                         raytracing_tl%z_r(:,iprof)
     ENDDO
    ELSE ! do sun
      n = 0
      DO iprof = 1, nprofiles
        ! z is temporary array because we overwrote ra
        z(:,iprof) = rearth(iprof) * raytracing%z_r(:,iprof)
        ra_tl(:,iprof) = rearth(iprof) * raytracing_tl%z_r(:,iprof)
      ENDDO
    ENDIF
  
    IF(opts%rt_all%ADDREFRAC .OR. opts%rt_ir%pc%ADDPC) THEN
      DO iprof = 1, nprofiles
! the (nlevels - n) takes in to account the fact that 10.2 and earlier
! had the refraction calculation on unpacked, reversed levels. (nlevels - n) [=2] maps to jplay for sat because a(jplay) -> a(-2)
        ra_tl(:,iprof) = &
          (z(:,iprof) * raytracing_tl%r(nlevels - n, iprof)  - &
          z(:,iprof) * raytracing%r(nlevels - n, iprof) * raytracing_tl%r(2:nlevels,iprof) * raytracing%r_r(2:nlevels, iprof) + &
          ra_tl(:,iprof) * raytracing%r(nlevels - n, iprof)) * raytracing%r_r(2:nlevels,iprof)
      ENDDO
    ENDIF
   
    DO iprof = 1, nprofiles
      IF(do_sat) THEN
        za_tl(1:nlayers,iprof) = ra_tl(:,iprof) * angles(iprof)%sinview
      ELSE ! do_sun
        IF(profiles(iprof)%sunzenangle >= 0.0 .AND. &
          profiles(iprof)%sunzenangle < max_sol_zen) THEN
          za_tl(1:nlayers,iprof) = ra_tl(:,iprof) * angles(iprof)%sinzen_sun
        ELSE
          za_tl(1:nlayers,iprof) = 0.0_jprb
        ENDIF
      ENDIF
    ENDDO

! calculate path from zenith angle using path = 1/sqrt(1-za**2)
    CALL SINTOSEC_TL(za(1:nlayers,:), za_tl(1:nlayers,:), path, path_tl)    

    IF(do_sat .AND. aux%on_coef_levels) THEN
      CALL INVSQRT_TL(raytracing%pathsat_rsqrt, raytracing_tl%pathsat, raytracing_tl%pathsat_rsqrt)

      raytracing_tl%pathsat_sqrt = raytracing_tl%pathsat * raytracing%pathsat_rsqrt + &
                                   raytracing%pathsat * raytracing_tl%pathsat_rsqrt
    ENDIF

  
  END SUBROUTINE CALC_SECZEN_TL
  
  SUBROUTINE calc_ssu_co2_thickness_tl
! pmr
    REAL(kind=jprb) :: pstd
    REAL(kind=jprb) :: tstd
    REAL(kind=jprb) :: pa(profiles(1)%nlevels-1),   pa_tl(profiles(1)%nlevels-1)
    REAL(kind=jprb) :: dpa(profiles(1)%nlevels-1),  dpa_tl(profiles(1)%nlevels-1)
    REAL(kind=jprb) :: rho(profiles(1)%nlevels-1),  rho_tl(profiles(1)%nlevels-1)
    REAL(kind=jprb) :: temp(profiles(1)%nlevels-1), temp_tl(profiles(1)%nlevels-1)
    REAL(kind=jprb) :: co2(profiles(1)%nlevels-1),  co2_tl(profiles(1)%nlevels - 1)
    REAL(kind=jprb) :: dmair(profiles(1)%nlevels),  dmair_tl(profiles(1)%nlevels)
    REAL(kind=jprb) :: dz, dz_tl
    REAL(kind=jprb) :: dzstp_tl
    
    pstd=p0*100.    ! pa
    tstd=t0         ! k

    gravh(:) = 1._jprb / gravh_r(:)

    DO iprof = 1, nprofiles
      CALL reciprocal(raytracing%dmair_r(:,iprof), dmair(:))
      CALL reciprocal_tl(dmair(:), raytracing_tl%dmair_r(:,iprof), dmair_tl(:))
! layer variables
      DO lay = 1, nlayers
! on layers top to bottom
        pa(lay)  = (profiles(iprof)%p(lay) + profiles(iprof)%p(lay+1)) * 100._jprb / 2._jprb  ! pa
        dpa(lay) = (profiles(iprof)%p(lay) - profiles(iprof)%p(lay+1)) * 100._jprb            ! pa
        rho(lay) = (dmair(lay) + dmair(lay+1)) / 2._jprb                                      ! kg/m3
        temp(lay) = (profiles(iprof)%t(lay) + profiles(iprof)%t(lay+1)) / 2._jprb             ! K
        co2(lay)  = (profiles_dry(iprof)%co2(lay) + profiles_dry(iprof)%co2(lay+1)) / 2._jprb ! ppmv dry
! TL
        IF (opts%interpolation%lgradp) THEN
          pa_tl(lay) = (profiles_tl(iprof)%p(lay) + profiles_tl(iprof)%p(lay+1)) * 100._jprb / 2._jprb
          dpa_tl(lay) = (profiles_tl(iprof)%p(lay) - profiles_tl(iprof)%p(lay+1)) * 100._jprb
        ELSE
          pa_tl(lay)   = 0._jprb
          dpa_tl(lay)  = 0._jprb
        ENDIF
        rho_tl(lay)  = (dmair_tl(lay) + &
                        dmair_tl(lay+1)) / 2._jprb
        temp_tl(lay) = (profiles_tl(iprof)%t(lay) + &
                        profiles_tl(iprof)%t(lay+1)) / 2._jprb
        co2_tl(lay)  = (profiles_dry_tl(iprof)%co2(lay) + &
                        profiles_dry_tl(iprof)%co2(lay+1)) / 2._jprb
      ENDDO
     
! calculation of co2 thickness at stp
      raytracing_tl%co2_cm(iprof) = 0._jprb

      DO lay = 1, nlayers
! on layers top to bottom
! fwd
         dz = -dpa(lay)/( (gravl(iprof)-gravh(iprof) * &
              raytracing%hgpl(lay+1, iprof)) * rho(lay) )    ! m
! tl
        IF (opts%interpolation%lgradp) THEN
          dz_tl = -dpa_tl(lay) / &
               ((gravl(iprof) - &
               gravh(iprof) * raytracing%hgpl(lay+1, iprof)) * &
               rho(lay)) + &
               dpa(lay) / &
               ((gravl(iprof)-gravh(iprof) * &
                 raytracing%hgpl(lay+1, iprof)) * rho(lay) )**2_jpim  * &
               ((gravl(iprof) - gravh(iprof) * &
                 raytracing%hgpl(lay+1, iprof)) * rho_tl(lay) - &
                gravh(iprof)*raytracing_tl%hgpl(lay+1, iprof) * rho(lay))

          dzstp_tl = dz_tl * pa(lay)/pstd * tstd/temp(lay)* co2(lay)*1e-6_jprb*100._jprb                    &
            +dz * pa_tl(lay)/pstd * tstd/temp(lay)* co2(lay)*1e-6_jprb*100._jprb                      &
            -dz * pa(lay)/pstd * tstd/temp(lay)**2_jpim * temp_tl(lay) * co2(lay)*1e-6_jprb*100._jprb  &
            +dz * pa(lay)/pstd * tstd/temp(lay)* co2_tl(lay)*1e-6_jprb*100._jprb

        ELSE
           dz_tl = dpa(lay)/&
                ((gravl(iprof)-gravh(iprof)* &
                  raytracing%hgpl(lay+1, iprof)) * rho(lay) )**2_jpim * &
                  ((gravl(iprof) - gravh(iprof)* &
                    raytracing%hgpl(lay+1, iprof))* rho_tl(lay) - &
                    gravh(iprof)*raytracing_tl%hgpl(lay+1, iprof) * &
                    rho(lay))

           dzstp_tl = dz_tl * pa(lay)/pstd * &
                      tstd/temp(lay)* co2(lay)*1e-6_jprb*100._jprb - &
                      dz * pa(lay)/pstd * tstd/temp(lay)**2_jpim * &
                      temp_tl(lay) * co2(lay)*1e-6_jprb*100._jprb +&
                      dz * pa(lay)/pstd * tstd/temp(lay)* &
                      co2_tl(lay)*1e-6_jprb*100._jprb
        ENDIF

        raytracing_tl%co2_cm(iprof) = raytracing_tl%co2_cm(iprof) + & 
                                      dzstp_tl
        
      ENDDO
    ENDDO
  END SUBROUTINE calc_ssu_co2_thickness_tl

END SUBROUTINE rttov_locpat_tl
