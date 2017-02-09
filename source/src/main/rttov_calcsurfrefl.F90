!
SUBROUTINE rttov_calcsurfrefl( &
            & coef,            &
            & profiles,        &
            & sunglint,        &
            & fresnrefl,       &
            & solar,           &
            & chanprof,        &
            & refl_norm,       &
            & calcrefl,        &
            & emissivity,      &
            & reflectance)
! Description:
! To compute short-wave surface reflectances for all channels and all
! profiles where required
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
!    Copyright 2011, EUMETSAT, All Rights Reserved.
!
! Method:
!
! Current Code Owner: SAF NWP
!
! History:
! Version   Date     Comment
! -------   ----     -------
!  1.0       05/10/2011  Created (J Hocking)
!
! Code Description:
!   Language:           Fortran 90.
!   Software Standards: "European Standards for Writing and
!     Documenting Exchangeable Fortran 90 Code".
!
! Declarations:
! Modules used:
! Imported Parameters:
! Imported Type Definitions:
  USE rttov_types, ONLY :  &
       & rttov_chanprof,   &
       & rttov_coef,       &
       & rttov_emissivity, &
       & profile_Type,     &
       & sunglint_type
  USE parkind1, ONLY : jprb, jplm
!INTF_OFF
  USE rttov_const, ONLY :     &
       & pi_r,                &
       & deg2rad,             &
       & surftype_sea,        &
       & surftype_seaice
  USE parkind1, ONLY : jpim
  USE yomhook, ONLY : LHOOK, DR_HOOK
!INTF_ON
  IMPLICIT NONE
!subroutine arguments:
  TYPE(rttov_chanprof   ), INTENT(IN)             :: chanprof(:)
  TYPE(profile_type     ), INTENT(IN)             :: profiles(:)
  TYPE(rttov_coef       ), INTENT(IN)             :: coef
  TYPE(sunglint_type    ), INTENT(IN)             :: sunglint
  REAL(KIND=jprb)        , INTENT(IN)             :: fresnrefl(size(chanprof))
  LOGICAL(KIND=jplm)     , INTENT(IN)             :: solar(size(chanprof))
  LOGICAL(KIND=jplm)     , INTENT(IN)             :: calcrefl(size(chanprof))
  TYPE(rttov_emissivity ), INTENT(IN),   OPTIONAL :: emissivity(size(chanprof))
  REAL(KIND=jprb)        , INTENT(INOUT)          :: reflectance(size(chanprof))
  REAL(KIND=jprb)        , INTENT(OUT)            :: refl_norm(size(chanprof))
!INTF_END
!local variables:
  INTEGER(KIND=jpim) :: j, prof, chan
  INTEGER(KIND=jpim) :: nchanprof
  REAL   (KIND=jprb) :: emis
  REAL   (KIND=JPRB) :: ZHOOK_HANDLE
!- End of header --------------------------------------------------------
  IF (LHOOK) CALL DR_HOOK('RTTOV_CALCSURFREFL', 0_jpim, ZHOOK_HANDLE)

  ! The output reflectances in the radiance structure are BRFs.

  ! rttov_integrate expects the reflectance(:) input array to contain BRDFs. These are
  ! multiplied by refl_norm(:) which is equal to cos(sol_zen_ang).

  ! The reflectance(:)%refl_out should obviously contain the same "type" of reflectance as
  ! reflectance(:)%refl_in.

  ! Currently we assume refl_in and refl_out contain BRDFs (not BRFs).

  ! In the case where calcrefl = FALSE, we copy refl_in to refl_out.

  ! In the case where calcrefl = TRUE, we use a BRF or a directional-hemispherical albedo
  ! and we need to multiply by 1/pi to get the BRDF, or in the case of sea surfaces we
  ! calculate the "BRDF" directly from the sunglint factor and fresnel coefficients, but
  ! this also takes the solar zenith angle into account so we scale it by refl_norm.

  nchanprof = size(chanprof)
  DO j = 1, nchanprof

    IF (.NOT. solar(j)) CYCLE

    prof = chanprof(j)%prof
    chan = chanprof(j)%chan

    refl_norm(j) = COS(profiles(prof)%sunzenangle * deg2rad)

    IF (.NOT. calcrefl(j)) CYCLE

    ! For land and sea-ice, the values chosen represent directional-hemispherical
    ! albedos. The reflectance(:) passed into rttov_integrate must be a BRDF, hence 
    ! the normalisation by 1/pi.

    IF (profiles(prof)%skin%surftype == surftype_sea) THEN

      reflectance(j) = sunglint%s(prof)%glint * fresnrefl(j)

      ! Scale the reflectance by cos(sunzenangle) because the sea surface model
      ! takes the solar zenith angle in account. This ensures the output reflectance
      ! is BRDF-like.
      IF (profiles(prof)%s2m%u /= 0._jprb .OR. profiles(prof)%s2m%v /= 0._jprb) THEN
        reflectance(j) = reflectance(j) / refl_norm(j)
      ENDIF

    ELSE IF (profiles(prof)%skin%surftype == surftype_seaice) THEN

      IF (PRESENT(emissivity)) THEN
        emis = emissivity(j)%emis_out
      ELSE
        emis = 0._jprb
      ENDIF

      IF (coef%ss_val_chn(chan) == 1 .AND. emis > 0.0_jprb) THEN
        reflectance(j) = (1.0_jprb - emis) * pi_r
      ELSE
        reflectance(j) = 0.8_jprb * pi_r  ! JAH - Placeholder rough value for sea-ice
      ENDIF

    ELSE

      IF (PRESENT(emissivity)) THEN
        emis = emissivity(j)%emis_out
      ELSE
        emis = 0._jprb
      ENDIF

      IF (coef%ss_val_chn(chan) == 1 .AND. emis > 0.0_jprb) THEN
        reflectance(j) = (1.0_jprb - emis) * pi_r
      ELSE
        reflectance(j) = 0.3_jprb * pi_r   ! JAH - Placeholder rough value for land
      ENDIF

    ENDIF

  ENDDO
  IF (LHOOK) CALL DR_HOOK('RTTOV_CALCSURFREFL', 1_jpim, ZHOOK_HANDLE)
END SUBROUTINE rttov_calcsurfrefl
