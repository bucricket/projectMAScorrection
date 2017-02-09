!
SUBROUTINE rttov_calcsurfrefl_k( &
            & coef,              &
            & profiles,          &
            & sunglint,          &
            & sunglint_k,        &
            & fresnrefl,         &
            & fresnrefl_k,       &
            & solar,             &
            & chanprof,          &
            & refl_norm,         &
            & calcrefl,          &
            & emissivity,        &
            & emissivity_k,      &
            & reflectance_k)
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
  TYPE(rttov_chanprof),   INTENT(IN)              :: chanprof(:)
  TYPE(profile_type),     INTENT(IN)              :: profiles(:)
  TYPE(rttov_coef),       INTENT(IN)              :: coef
  TYPE(sunglint_type),    INTENT(IN)              :: sunglint
  TYPE(sunglint_type),    INTENT(INOUT)           :: sunglint_k
  REAL(KIND=jprb),        INTENT(IN)              :: fresnrefl(size(chanprof))
  REAL(KIND=jprb),        INTENT(INOUT)           :: fresnrefl_k(size(chanprof))
  LOGICAL(KIND=jplm),     INTENT(IN)              :: solar(size(chanprof))
  REAL(KIND=jprb),        INTENT(IN)              :: refl_norm(size(chanprof))
  LOGICAL(KIND=jplm),     INTENT(IN)              :: calcrefl(size(chanprof))
  TYPE(rttov_emissivity), INTENT(IN),    OPTIONAL :: emissivity(size(chanprof))
  TYPE(rttov_emissivity), INTENT(INOUT), OPTIONAL :: emissivity_k(size(chanprof))
  REAL(KIND=jprb),        INTENT(IN)              :: reflectance_k(size(chanprof))
!INTF_END
!local variables:
  INTEGER(KIND=jpim) :: j, prof, chan
  INTEGER(KIND=jpim) :: nchanprof
  REAL   (KIND=jprb) :: emis, tmp_refl_k
  REAL   (KIND=JPRB) :: ZHOOK_HANDLE
!- End of header --------------------------------------------------------
  IF (LHOOK) CALL DR_HOOK('RTTOV_CALCSURFREFL_K', 0_jpim, ZHOOK_HANDLE)

  nchanprof = size(chanprof)
  DO j = 1, nchanprof

    IF (.NOT. solar(j)) CYCLE
    IF (.NOT. calcrefl(j)) CYCLE

    prof = chanprof(j)%prof
    chan = chanprof(j)%chan

    IF (profiles(prof)%skin%surftype == surftype_sea) THEN

      ! The output reflectance_k should not be scaled by cos(sunzenangle)
      ! Only the refl_k used in the Jacobian sea surface model calculation should be scaled
      IF (profiles(prof)%s2m%u /= 0._jprb .OR. profiles(prof)%s2m%v /= 0._jprb) THEN
        tmp_refl_k = reflectance_k(j) / refl_norm(j)
      ELSE
        tmp_refl_k = reflectance_k(j)
      ENDIF

      sunglint_k%s(j)%glint = sunglint_k%s(j)%glint + tmp_refl_k * fresnrefl(j)
      fresnrefl_k(j)        = fresnrefl_k(j)        + tmp_refl_k * sunglint%s(prof)%glint

    ELSE IF (profiles(prof)%skin%surftype == surftype_seaice) THEN

      IF (PRESENT(emissivity_k)) THEN

        IF (PRESENT(emissivity)) THEN
          emis = emissivity(j)%emis_out
        ELSE
          emis = 0._jprb
        ENDIF

        IF (coef%ss_val_chn(chan) == 1 .AND. emis > 0.0_jprb) THEN
          emissivity_k(j)%emis_out = emissivity_k(j)%emis_out - reflectance_k(j) * pi_r
        ENDIF

      ENDIF

    ELSE

      IF (PRESENT(emissivity_k)) THEN

        IF (PRESENT(emissivity)) THEN
          emis = emissivity(j)%emis_out
        ELSE
          emis = 0._jprb
        ENDIF

        IF (coef%ss_val_chn(chan) == 1 .AND. emis > 0.0_jprb) THEN
          emissivity_k(j)%emis_out = emissivity_k(j)%emis_out - reflectance_k(j) * pi_r
        ENDIF

      ENDIF

    ENDIF

  ENDDO
  IF (LHOOK) CALL DR_HOOK('RTTOV_CALCSURFREFL_K', 1_jpim, ZHOOK_HANDLE)
END SUBROUTINE rttov_calcsurfrefl_k
