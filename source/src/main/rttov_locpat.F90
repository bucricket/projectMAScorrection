!     Calculate the secant of the zenith angle at the lower boundary
!     of each atmospheric layer.
SUBROUTINE rttov_locpat(opts, dosolar, profiles, profiles_dry, &
                        aux, coef, angles, raytracing, do_pmc_calc)
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
  TYPE(profile_type   ), INTENT(IN)    :: profiles_dry(:)
  TYPE(profile_aux    ), INTENT(IN)    :: aux
  TYPE(geometry_type  ), INTENT(IN)    :: angles(:)
  TYPE(rttov_coef     ), INTENT(IN)    :: coef
  TYPE(raytracing_type), INTENT(INOUT) :: raytracing
  LOGICAL(jplm)        , INTENT(IN), OPTIONAL :: do_pmc_calc
!INTF_END

  INTEGER(KIND=jpim) :: iprof, lev, lay
  INTEGER(KIND=jpim) :: nlevels, nlayers
  REAL   (KIND=jprb) :: rearth(SIZE(profiles))
  REAL   (KIND=jprb) :: gravl(SIZE(profiles))   ! gravity at latitude lat [m/s^2]
  REAL   (KIND=jprb) :: gravh(SIZE(profiles)), gravh_r(SIZE(profiles))
  REAL   (KIND=jprb) :: eta, beta ! coefficients of the international gravity formula
  REAL   (KIND=jprb) :: rlh(SIZE(profiles))
  REAL   (KIND=jprb) :: dflat
  REAL   (KIND=jprb) :: fac, c

  INTEGER(KIND=jpim) :: nprofiles
  LOGICAL(KIND=jplm) :: do_pmc_calc1

  REAL   (KIND=jprb) :: scale_eps
  REAL   (KIND=jprb) :: small_val

  REAL   (KIND=jprb) :: ZHOOK_HANDLE
!-----end of header-------------------------------------------------------------
  IF (LHOOK) CALL DR_HOOK('RTTOV_LOCPAT', 0_jpim, ZHOOK_HANDLE)

  nprofiles = SIZE(profiles)
  nlevels     = profiles(1)%nlevels
  nlayers     = nlevels - 1

  small_val = (TINY(1._jprb)) ** (0.333333_jprb)

  do_pmc_calc1 = .FALSE.
  IF (PRESENT(do_pmc_calc)) do_pmc_calc1 = do_pmc_calc

! DAR: Move constant calculations to top of locpat subroutine

! scale water profile in air density calculation
  scale_eps = (1._jprb - mh2o/mair)

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

  c = 1000._jprb * rgc / (100._jprb * mair)

  DO iprof = 1, nprofiles
! first calculate density for dry air (ideal gas)

!     raytracing%dair(:, iprof) = &
!           1.e+02_jprb * profiles(iprof)%p(:) * mair / &
!           (1000._jprb * rgc * profiles(iprof)%t(:)) ! THIS IS DENSITY OF DRY AIR TO HERE
!      CALL reciprocal(&
!        raytracing%dair(:, iprof) * (1._jprb - profiles(iprof)%q(:) * scale_eps), &
!        raytracing%dmair_r(:,iprof))

    ! calculate partial pressure of water vapour
    IF (profiles(iprof)%gas_units == gas_unit_compatibility) THEN
      ! Do nothing in compatibility mode: using input profiles with no conversion
      raytracing%ppw(:, iprof) = profiles(iprof)%p(:) * 1e-6_jprb * profiles_dry(iprof)%q(:)
    ELSE
      ! Convert ppmv dry to ppmv wet
      raytracing%ppw(:, iprof) = profiles(iprof)%p(:) * 1e-6_jprb * &
                                 profiles_dry(iprof)%q(:) / (1._jprb + profiles_dry(iprof)%q(:) * 1e-6_jprb)
    ENDIF

    CALL DIVIDE(c * profiles(iprof)%t(:), &
      profiles(iprof)%p(:) - raytracing%ppw(:,iprof) * scale_eps, &
      raytracing%dmair_r(:, iprof), acc = .FALSE._jplm)
!       CALL reciprocal_1d(raytracing%dmair_r(:, iprof), raytracing%dmair(:, iprof))
  ENDDO

!-------4.  set up the height of the surface pressure level---------------------
! height of nearest level below surf counted from bottom of profile grid
  DO iprof = 1, nprofiles
    raytracing%hgpl(aux%s(iprof)%nearestlev_surf, iprof) = &
      profiles(iprof)%elevation
  ENDDO

!-------5. Compute height of pressure levels
  CALL calc_hgpl()

!-------5. compute atmospheric refractive index (if reqd.)---------------------
! compute the refractivity index for dry air at temperature tempp and pressure
! as a function of wave-number using an updated version of edlen equation
! not yet corrected for CO2 contribution - used in calc_seczen
  IF(opts%rt_all%ADDREFRAC .OR. opts%rt_ir%pc%ADDPC) CALL calc_refractivity()

!-------6.  Compute secant of the zenith angle at the lower boundary of each -
!       layer for the satellite-to-surface line of view and of each             |
!       layer for the sun-to-surface line of view (if possible and reqd.)       |

  CALL calc_seczen(raytracing%ratoesat, raytracing%zasat, raytracing%pathsat,&
                   do_sat=.TRUE._jplm) ! satellite call

  IF (dosolar) THEN
    CALL calc_seczen(raytracing%ratoesun, raytracing%zasun, raytracing%pathsun,&
                     do_sat=.FALSE._jplm) ! solar call
    raytracing%patheff = raytracing%pathsat + raytracing%pathsun
  ENDIF

! Extra bits
  DO iprof = 1, nprofiles
    raytracing%hgpl(1, iprof) = raytracing%hgpl(2, iprof)
    IF (opts%rt_ir%addaerosl .OR. opts%rt_ir%addclouds) THEN
      DO lay = 1, nlayers ! layer thickness
        raytracing%ltick(lay, iprof) = raytracing%hgpl(lay, iprof) - &
          raytracing%hgpl(lay + 1, iprof)
      ENDDO
    ENDIF
  ENDDO

!-------7.  Compute the STP co2 thickness in cm
  IF (coef%pmc_shift .AND. do_pmc_calc1) CALL calc_ssu_co2_thickness()

  IF (LHOOK) CALL DR_HOOK('RTTOV_LOCPAT', 1_jpim, ZHOOK_HANDLE)
CONTAINS

  SUBROUTINE calc_refractivity()

    REAL(KIND=jprb) :: dry_pp_h(profiles(1)%nlevels)
    REAL(KIND=jprb) :: t_norm(profiles(1)%nlevels)

    REAL(KIND=jprb) :: DISP ! the value of the refractive index given by the dispersion equation.
    REAL(KIND=jprb) :: DISP_MOIST
    REAL(KIND=jprb) :: DISPCO2

! calculate constants

! dry air dispersion equation constant in updated Edlen eqn
    DISP = 1.0e-8_jprb * ( &
      d1 + d2 / (d3 - (waver * ctom) ** 2_jpim) + &
           d4 / (d5 - (waver * ctom) ** 2_jpim))

! moist air dispersion constant
    DISP_MOIST = htop * (ew1 - ew2 * (waver * ctom) ** 2_jpim) * 1.0e-10_jprb

! compute and store (multiplicative) correction factor to account for presence of CO2 in the air
! assume constant CO2 using figure from 2005 (376ppmv) unless co2 data supplied
! calculation assumes ppmv over dry air
    IF (opts%rt_ir%co2_data) THEN
      DO iprof = 1, nprofiles
        raytracing%dispco2(:, iprof) = disp * &
          (1._jprb + dco2 * (profiles_dry(iprof)%co2(:) * 1e-6_jprb - 0.0003_jprb))
      ENDDO
    ELSE
      DISPCO2 = DISP * (1._jprb + dco2 * (CO2_CONC * 1e-6_jprb - 0.0003_jprb))
    ENDIF

    DO iprof = 1, nprofiles
! calculate useful temporary array to break refractivity calculation up
! in to a better size

      dry_pp_h(:) = htop * (profiles(iprof)%p(:) - raytracing%ppw(:, iprof))
      t_norm(:) = profiles(iprof)%t(:) - t0

      raytracing%refractivity(:, iprof) = (dry_pp_h(:) / ed1) *         &
        (1._jprb + 1.0e-8_jprb * &
        (ed2 - ed3 * t_norm(:)) * dry_pp_h(:)) /  &
        (1._jprb + ed4 * t_norm(:))

! apply CO2 correction then
! compute moist contribution to refractivity
! and store final refractive index for atmospheric profile

      IF (opts%rt_ir%co2_data) THEN
        raytracing%r(:, iprof) = (raytracing%refractivity(:, iprof) * &
          raytracing%dispco2(:, iprof)) - &
          (raytracing%ppw(:, iprof) * DISP_MOIST) + 1._jprb
      ELSE
        raytracing%r(:, iprof) = raytracing%refractivity(:, iprof) * DISPCO2 - &
          (raytracing%ppw(:, iprof) * DISP_MOIST) + 1._jprb
      ENDIF
    ENDDO
! store reciprocal of refractive index for performance
    CALL reciprocal(raytracing%r, raytracing%r_r)
  END SUBROUTINE calc_refractivity

  SUBROUTINE calc_hgpl()
    REAL(KIND=jprb) :: dp((profiles(1)%nlevels),SIZE(profiles))
    DO iprof = 1, nprofiles

      dp(1:nlevels-1, iprof) = 0.5_jprb * (profiles(iprof)%p(1:nlevels-1)  - &
                                         profiles(iprof)%p(2:nlevels))
    ENDDO

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
!------------------------------------------------------------------------------

!-------compute height of pressure levels above surface. ------------
! lower boundary of first layer is nearest level below surface
! upper boundary of last layer is toa
    DO iprof = 1, nprofiles
      DO lev = aux%s(iprof)%nearestlev_surf - 1, 1, -1
        raytracing%int(lev, iprof)  = & ! integrated layer values for dmair
          -dp(lev,iprof) *      &
          (raytracing%dmair_r(lev, iprof) + &
          raytracing%dmair_r(lev + 1, iprof))

        raytracing%ztemp(lev, iprof) = &
          rlh(iprof)**2_jpim - 1.0e3_jprb * raytracing%hgpl(lev + 1, iprof) * &
          (2.0_jprb * rlh(iprof) - 1.0e3_jprb * raytracing%hgpl(lev + 1, iprof)) - &
           2.0_jprb * raytracing%int(lev, iprof) * 1e2_jprb * gravh_r(iprof)

        IF (raytracing%ztemp(lev, iprof) > 0._jprb) THEN
          raytracing%hgpl(lev, iprof) = (rlh(iprof) - SQRT(raytracing%ztemp(lev, iprof))) * 1.0e-3_jprb
        ELSE
          raytracing%ztemp(lev, iprof) = small_val
          raytracing%hgpl(lev, iprof) = rlh(iprof) * 1.0e-3_jprb
        ENDIF
     ENDDO

!-------compute height of pressure levels below surface-------------
! upper boundary of first layer is nearest level below surface
! lower boundary of last layer is bottom of profile grid
      DO lev = aux%s(iprof)%nearestlev_surf + 1, nlevels
        raytracing%int(lev, iprof)  = &! integrated layer values for dmair
          dp(lev - 1,iprof) *      &
          (raytracing%dmair_r(lev - 1, iprof) + &
          raytracing%dmair_r(lev, iprof))

        raytracing%ztemp(lev, iprof) = &
          rlh(iprof) ** 2_jpim - 1.0e3_jprb * raytracing%hgpl(lev - 1, iprof) * &
          (2.0_jprb * rlh(iprof) - 1.0e3_jprb * raytracing%hgpl(lev - 1, iprof)) - &
          2.0_jprb * raytracing%int(lev, iprof) * 1.0e2_jprb * gravh_r(iprof)

        raytracing%hgpl(lev, iprof) = &
          (rlh(iprof) - SQRT(raytracing%ztemp(lev, iprof))) * 1.0e-3_jprb
      ENDDO
      ! FOR TL/AD/K
      raytracing%ztemp(aux%s(iprof)%nearestlev_surf, iprof) = &
        1._jprb ! DAR: NOT USED
    ENDDO

  END SUBROUTINE calc_hgpl

!-------Compute secant of the zenith angle at the lower boundary of each -
!       layer for the sat/sun-to-surface line of view.                       |
  SUBROUTINE calc_seczen(ra,za,path,do_sat)
!       The atmospheric layers are considered as concentric rings. If we trace |
!       a ray across these rings at any angle other than nadir, the local      |
!       angle relative to the outward radial direction at the point of         |
!       intersection will be different at each ring because due to the         |
!       curvature of the Earth and to atmospheric refraction. The secant of    |
!       the local/solar zenith angle PATHSAT/SUN is thus computed taking into  |
!       account the geometry of the situation and the bending of rays as they  |
!       traverse the atmosphere (by application of Snell's law).               |
!------------------------------------------------------------------------------

    REAL(jprb), INTENT(OUT) :: ra(:,:), za(:,:), path(:,:)
    LOGICAL(jplm), INTENT(IN) :: do_sat ! do satellite secant or solar secant

    INTEGER(jpim) :: n

! calculate and store pressure level reciprocal distance from earth centre
    DO iprof = 1, nprofiles
      CALL reciprocal(rearth(iprof) + raytracing%hgpl(2:nlayers+1, iprof), &
        raytracing%z_r(:,iprof))
    ENDDO

    IF (do_sat) THEN
      n = nlayers - 1
      DO iprof = 1, nprofiles
        ra(:,iprof) = (rearth(iprof) + coef%fc_sat_height) * &
                      raytracing%z_r(:,iprof)
      ENDDO
    ELSE ! do sun
      n = 0
      DO iprof = 1, nprofiles
        ra(:,iprof) = rearth(iprof) * raytracing%z_r(:,iprof)
      ENDDO
    ENDIF

    IF(opts%rt_all%ADDREFRAC .OR. opts%rt_ir%pc%ADDPC) THEN
      DO iprof = 1, nprofiles
! the (nlevels - n) takes in to account the fact that 10.2 and earlier
! had the refraction calculation on unpacked, reversed levels.
        ra(:,iprof) = ra(:,iprof) * raytracing%r(nlevels - n, iprof) * &
                      raytracing%r_r(2:nlevels,iprof)
      ENDDO
    ENDIF

    DO iprof = 1, nprofiles
      IF (do_sat) THEN
        za(1:nlayers,iprof) = ra(:,iprof) * angles(iprof)%sinview
      ELSE ! do_sun
        IF(profiles(iprof)%sunzenangle >= 0.0 .AND. &
          profiles(iprof)%sunzenangle < max_sol_zen) THEN
          za(1:nlayers,iprof) = ra(:,iprof) * angles(iprof)%sinzen_sun
        ELSE
          za(1:nlayers,iprof) = 1._jprb + 5e-39_jprb ! make pathsun ~1e38
        ENDIF
      ENDIF
    ENDDO

! calculate path from zenith angle using path = 1/sqrt(1-za**2)
    CALL SINTOSEC(za(1:nlayers,:), path)

    IF (do_sat .AND. aux%on_coef_levels) THEN
      CALL INVSQRT(raytracing%pathsat, raytracing%pathsat_rsqrt)
      raytracing%pathsat_sqrt = raytracing%pathsat * raytracing%pathsat_rsqrt
    ENDIF

  END SUBROUTINE CALC_SECZEN

  SUBROUTINE calc_ssu_co2_thickness
! pmr
    REAL   (KIND=jprb) :: pstd
    REAL   (KIND=jprb) :: tstd
    REAL   (KIND=jprb) :: pa(profiles(1)%nlevels - 1)
    REAL   (KIND=jprb) :: dpa(profiles(1)%nlevels - 1)
    REAL   (KIND=jprb) :: rho(profiles(1)%nlevels - 1)
    REAL   (KIND=jprb) :: temp(profiles(1)%nlevels - 1)
    REAL   (KIND=jprb) :: co2(profiles(1)%nlevels - 1)
    REAL   (KIND=jprb) :: dmair(profiles(1)%nlevels)
    REAL   (KIND=jprb) :: dz
    REAL   (KIND=jprb) :: dzstp

    pstd=p0*100.    ! pa
    tstd=t0         ! k

! Get density moist air from its reciprocal because it's only used in this calc.

    gravh(:) = 1._jprb / gravh_r(:)

    DO iprof = 1, nprofiles
      CALL reciprocal(raytracing%dmair_r(:,iprof), dmair(:))
! layer variables
      DO lay = 1, nlayers
! on layers top to bottom
        pa(lay)   = (profiles(iprof)%p(lay) + profiles(iprof)%p(lay+1)) * 100._jprb / 2._jprb ! pa
        dpa(lay)  = (profiles(iprof)%p(lay) - profiles(iprof)%p(lay+1)) * 100._jprb           ! pa
        rho(lay)  = (dmair(lay) + dmair(lay+1)) / 2._jprb                                     ! kg/m3
        temp(lay) = (profiles(iprof)%t(lay) + profiles(iprof)%t(lay+1)) / 2._jprb             ! K
        co2(lay)  = (profiles_dry(iprof)%co2(lay) + profiles_dry(iprof)%co2(lay+1)) / 2._jprb ! ppmv dry
      ENDDO

! calculation of co2 thickness at stp
      raytracing%co2_cm(iprof) = 0._jprb

      DO lay = 1, nlayers
! on layers top to bottom
        dz = -dpa(lay)/( (gravl(iprof)-gravh(iprof) * &
          raytracing%hgpl(lay+1, iprof)) * rho(lay) )    ! m
! pv=rt so v=rt/p so h=rt/p when v is height h of unit column
        dzstp = dz * pa(lay)/pstd * tstd/temp(lay) * co2(lay) * 1e-6_jprb * 100._jprb   ! cm
! nb equal v contains equal numbers of molecules (e.g. co2) - result in cm

        raytracing%co2_cm(iprof) = raytracing%co2_cm(iprof) + dzstp
      ENDDO
    ENDDO
  END SUBROUTINE calc_ssu_co2_thickness

END SUBROUTINE rttov_locpat
