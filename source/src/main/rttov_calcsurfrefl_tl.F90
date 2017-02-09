!
SUBROUTINE rttov_calcsurfrefl_tl( &
            & coef,               &
            & profiles,           &
            & sunglint,           &
            & sunglint_tl,        &
            & fresnrefl,          &
            & fresnrefl_tl,       &
            & solar,              &
            & chanprof,           &
            & refl_norm,          &
            & calcrefl,           &
            & emissivity,         &
            & emissivity_tl,      &
            & reflectance_tl)
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
       & surftype_sea,        &
       & surftype_seaice
  USE parkind1, ONLY : jpim
  USE yomhook, ONLY : LHOOK, DR_HOOK
!INTF_ON
  IMPLICIT NONE
!subroutine arguments:
  TYPE(rttov_chanprof),   INTENT(IN)             :: chanprof(:)
  TYPE(profile_type),     INTENT(IN)             :: profiles(:)
  TYPE(rttov_coef),       INTENT(IN)             :: coef
  TYPE(sunglint_type),    INTENT(IN)             :: sunglint
  TYPE(sunglint_type),    INTENT(IN)             :: sunglint_tl
  REAL(KIND=jprb),        INTENT(IN)             :: fresnrefl(size(chanprof))
  REAL(KIND=jprb),        INTENT(IN)             :: fresnrefl_tl(size(chanprof))
  LOGICAL(KIND=jplm),     INTENT(IN)             :: solar(size(chanprof))
  REAL(KIND=jprb),        INTENT(IN)             :: refl_norm(size(chanprof))
  LOGICAL(KIND=jplm),     INTENT(IN)             :: calcrefl(size(chanprof))
  TYPE(rttov_emissivity), INTENT(IN),   OPTIONAL :: emissivity(size(chanprof))
  TYPE(rttov_emissivity), INTENT(IN),   OPTIONAL :: emissivity_tl(size(chanprof))
  REAL(KIND=jprb),        INTENT(INOUT)          :: reflectance_tl(size(chanprof))
!INTF_END
!local variables:
  INTEGER(KIND=jpim) :: j, prof, chan
  INTEGER(KIND=jpim) :: nchanprof
  REAL   (KIND=jprb) :: emis, emis_tl
  REAL   (KIND=JPRB) :: ZHOOK_HANDLE
!- End of header --------------------------------------------------------
  IF (LHOOK) CALL DR_HOOK('RTTOV_CALCSURFREFL_TL', 0_jpim, ZHOOK_HANDLE)

  nchanprof = size(chanprof)
  DO j = 1, nchanprof

    IF (.NOT. solar(j)) CYCLE
    IF (.NOT. calcrefl(j)) CYCLE

    prof = chanprof(j)%prof
    chan = chanprof(j)%chan

    IF (profiles(prof)%skin%surftype == surftype_sea) THEN

      reflectance_tl(j) = sunglint_tl%s(prof)%glint * fresnrefl(j) + &
                        & sunglint%s(prof)%glint * fresnrefl_tl(j)

      IF (profiles(prof)%s2m%u /= 0._jprb .OR. profiles(prof)%s2m%v /= 0._jprb) THEN
        reflectance_tl(j) = reflectance_tl(j) / refl_norm(j)
      ENDIF

    ELSE IF (profiles(prof)%skin%surftype == surftype_seaice) THEN

      IF (PRESENT(emissivity) .AND. PRESENT(emissivity_tl)) THEN
        emis = emissivity(j)%emis_out
        emis_tl = emissivity_tl(j)%emis_out
      ELSE
        emis = 0._jprb
        emis_tl = 0._jprb
      ENDIF

      IF (coef%ss_val_chn(chan) == 1 .AND. emis > 0.0_jprb) THEN
        reflectance_tl(j) = - emis_tl * pi_r
      ELSE
! User should supply TL, same as for emissivity
!         reflectance_tl(j) = 0.0_jprb
      ENDIF

    ELSE

      IF (PRESENT(emissivity) .AND. PRESENT(emissivity_tl)) THEN
        emis = emissivity(j)%emis_out
        emis_tl = emissivity_tl(j)%emis_out
      ELSE
        emis = 0._jprb
        emis_tl = 0._jprb
      ENDIF

      IF (coef%ss_val_chn(chan) == 1 .AND. emis > 0.0_jprb) THEN
        reflectance_tl(j) = - emis_tl * pi_r
      ELSE
! User should supply TL, same as for emissivity
!         reflectance_tl(j) = 0.0_jprb
      ENDIF

    ENDIF

  ENDDO
  IF (LHOOK) CALL DR_HOOK('RTTOV_CALCSURFREFL_TL', 1_jpim, ZHOOK_HANDLE)
END SUBROUTINE rttov_calcsurfrefl_tl
