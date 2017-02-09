SUBROUTINE rttov_fresnel_k( &
            & chanprof,    &
            & profiles,    &
            & solar,       &
            & coef,        &
            & sunglint,    &
            & sunglint_k,  &
            & fresnrefl,   &
            & fresnrefl_k)
!     Description:
!     Compute sea surface emissivity for radiative transfer calculations
!     Set up surface emissivities for all channels and all
!     profiles.
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
!    Copyright 2010, EUMETSAT, All Rights Reserved.
!
!     Method:
!     See K.MASUDA,T.TAKASHIMA,Y.TAKAYAMA.
!     Emissivity of pure water and sea waters for the
!     sea surface in the infrared window regions.
!     Remote Sensing of Environment 24:313-329 (1988)
!
!     Owner:
!     EUMETSAT
!
!     History:
!     Version      Date        Comment
!     1            21/03/99.   Original Code. M. Matricardi. ECMWF.
!     2            14/03/2002  Original Code. M. Matricardi. ECMWF.
!     3            15/07/2003  RTIASI-4. M. Matricardi. ECMWF.
!                              Rewritten to include the computation of the
!                              Fresnel reflectivity for a flat water surface
!     4            06/12/2005  Modified intent of sunglint_k (R Saunders)
!     5            05/02/2007  Remove polarisation (R Saunders)
!     6            05/07/2010  Remove addsolar flag from profiles structure (J Hocking)

! 2010/03 Code cleaning Pascal Brunel, Philippe Marguinaud
!
!     Code description:
!       Language:              Fortran 90.
!       Software Standards:    "European Standards for Writing and Documenting
!                              Exchangeable Fortran 90 code".
!
!     Module used:
!     Imported Type Definitions:
  USE rttov_types, ONLY :  &
       & rttov_chanprof, &
       & profile_Type,   &
       & rttov_coef,     &
       & sunglint_type
!     Imported parameters:
  USE parkind1, ONLY : jprb, jplm
!INTF_OFF
  USE rttov_const, ONLY :       &
       & surftype_sea,          &
       & watertype_fresh_water, &
       & watertype_ocean_water
  USE yomhook, ONLY : LHOOK, DR_HOOK
  USE parkind1, ONLY : jpim
!INTF_ON
  IMPLICIT NONE
!     Subroutine arguments:
  TYPE(rttov_chanprof), INTENT(IN)    :: chanprof(:)
  TYPE(profile_type  ), INTENT(IN)    :: profiles(:)
  LOGICAL(KIND=jplm)  , INTENT(IN)    :: solar(size(chanprof))
  TYPE(rttov_coef    ), INTENT(IN)    :: coef
  TYPE(sunglint_type ), INTENT(IN)    :: sunglint
  TYPE(sunglint_type ), INTENT(INOUT) :: sunglint_k
  REAL(KIND=jprb)     , INTENT(IN)    :: fresnrefl  (size(chanprof))
  REAL(KIND=jprb)     , INTENT(IN)    :: fresnrefl_k(size(chanprof))
!INTF_END
!     End of subroutine arguments
!       Local scalars:
  INTEGER(KIND=jpim) :: i, prof, chan
  REAL   (KIND=jprb) :: sincsi
  REAL   (KIND=jprb) :: coscsi
  COMPLEX(KIND=jprb) :: crefpar
  COMPLEX(KIND=jprb) :: crefparj
  COMPLEX(KIND=jprb) :: crefperp
  COMPLEX(KIND=jprb) :: crefperpj
  COMPLEX(KIND=jprb) :: sincsi1
  COMPLEX(KIND=jprb) :: coscsi1
  REAL   (KIND=jprb) :: sincsi_k
  REAL   (KIND=jprb) :: coscsi_k
  COMPLEX(KIND=jprb) :: crefpar_k
  COMPLEX(KIND=jprb) :: crefparj_k
  COMPLEX(KIND=jprb) :: crefperp_k
  COMPLEX(KIND=jprb) :: crefperpj_k
  COMPLEX(KIND=jprb) :: sincsi1_k
  COMPLEX(KIND=jprb) :: coscsi1_k
  INTEGER(KIND=jpim) :: nchannels            ! Number of radiances computed (channels used * profiles)
!       Local arrays:
  COMPLEX(KIND=jprb) :: waopc(size(chanprof))
  REAL   (KIND=JPRB) :: ZHOOK_HANDLE
!-----End of header-------------------------------------------------------------
  IF (LHOOK) CALL DR_HOOK('RTTOV_FRESNEL_K', 0_jpim, ZHOOK_HANDLE)
  nchannels = size(chanprof)
  DO i = 1, nchannels
    chan = chanprof(i)%chan
    prof = chanprof(i)%prof
    IF (solar(i)) THEN
      IF (profiles(prof)%skin%surftype == surftype_sea) THEN!Water surface type
        IF (profiles(prof)%skin%watertype == watertype_fresh_water) THEN
          waopc(i) = coef%woc_waopc_fw(chan)
        ELSE IF (profiles(prof)%skin%watertype == watertype_ocean_water) THEN
          waopc(i) = coef%woc_waopc_ow(chan)
        ENDIF
!-------Compute infrared sunglint for a flat water surface-------------------
!-------Compute direct quantities-----------------------------------------------
        coscsi                = cos(sunglint%s(prof)%omega)
        sincsi                = sin(sunglint%s(prof)%omega)
        sincsi1               = sincsi / waopc(i)
        coscsi1               = sqrt(1._jprb - sincsi1 ** 2_jpim)
        crefpar               =  - (waopc(i) * coscsi - coscsi1) / (waopc(i) * coscsi + coscsi1)
        crefparj              = conjg(crefpar)
        crefperp              = (coscsi - waopc(i) * coscsi1) / (coscsi + waopc(i) * coscsi1)
        crefperpj             = conjg(crefperp)
!-------------------------------------------------------------------------------
        crefpar_k             = fresnrefl_k(i) * 0.5_jprb * crefparj
        crefparj_k            = fresnrefl_k(i) * 0.5_jprb * crefpar
        crefperp_k            = fresnrefl_k(i) * 0.5_jprb * crefperpj
        crefperpj_k           = fresnrefl_k(i) * 0.5_jprb * crefperp
        crefperp_k            = crefperp_k + conjg(crefperpj_k)
        crefperpj_k           = 0
        coscsi_k              = REAL (crefperp_k / (coscsi + waopc(i) * coscsi1), jprb)
        coscsi1_k             =  - crefperp_k * waopc(i) / (coscsi + waopc(i) * coscsi1)
        coscsi_k              = REAL (coscsi_k - crefperp_k * (coscsi - waopc(i) * coscsi1) / &
                                      (coscsi + waopc(i) * coscsi1) ** 2_jpim, jprb)
        coscsi1_k             = coscsi1_k - crefperp_k * waopc(i) * (coscsi - waopc(i) * coscsi1) / &
                                (coscsi + waopc(i) * coscsi1) ** 2_jpim
        crefperp_k            = 0
        crefpar_k             = crefpar_k + conjg(crefparj_k)
        crefparj_k            = 0
        coscsi_k              = REAL (coscsi_k - crefpar_k * waopc(i) / (waopc(i) * coscsi + coscsi1), jprb)
        coscsi1_k             = coscsi1_k + crefpar_k / (waopc(i) * coscsi + coscsi1)
        coscsi_k              = REAL (coscsi_k + crefpar_k * waopc(i) * ( - coscsi1 + waopc(i) * coscsi) / &
                                      (waopc(i) * coscsi + coscsi1) ** 2_jpim, jprb)
        coscsi1_k             = coscsi1_k + crefpar_k * ( - coscsi1 + waopc(i) * coscsi) / &
                                (waopc(i) * coscsi + coscsi1) ** 2_jpim
        crefpar_k             = 0
        sincsi1_k             =  - coscsi1_k * sincsi1 / coscsi1
        coscsi1_k             = 0
        sincsi_k              = REAL (sincsi1_k / waopc(i), jprb)
        sincsi1_k             = 0
        sunglint_k%s(i)%omega = sunglint_k%s(i)%omega + sincsi_k * coscsi
        sunglint_k%s(i)%omega = sunglint_k%s(i)%omega - coscsi_k * sincsi
        sincsi_k              = 0
        coscsi_k              = 0
      ENDIF
    ENDIF
  ENDDO
!End of channel loop
  IF (LHOOK) CALL DR_HOOK('RTTOV_FRESNEL_K', 1_jpim, ZHOOK_HANDLE)
END SUBROUTINE rttov_fresnel_k
