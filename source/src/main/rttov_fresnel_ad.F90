!     Compute sea surface emissivity for radiative transfer calculations
SUBROUTINE rttov_fresnel_ad( &
            & chanprof,     &
            & profiles,     &
            & solar,        &
            & coef,         &
            & sunglint,     &
            & sunglint_ad,  &
            & fresnrefl,    &
            & fresnrefl_ad)
!     Description:
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
!     4            06/12/2005  Modified intent to inout for sunglint_ad (R Saunders)
!     5            05/02/2007  Removed polarisation indices (R Saunders)
!     6            05/07/2010  Remove addsolar flag from profiles structure (J Hocking)
!
! 2010/03 Code cleaning Pascal Brunel, Philippe Marguinaud
!
!     Code description:
!       Language:              Fortran 90.
!       Software Standards:    "European Standards for Writing and Documenting
!                              Exchangeable Fortran 90 code".
!
!     Module used:
  USE rttov_types, ONLY :  &
       & rttov_chanprof, &
       & profile_Type,   &
       & rttov_coef,     &
       & sunglint_type
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
  TYPE(sunglint_type ), INTENT(INOUT) :: sunglint_ad
  REAL(KIND=jprb)     , INTENT(IN)    :: fresnrefl   (size(chanprof))
  REAL(KIND=jprb)     , INTENT(IN)    :: fresnrefl_ad(size(chanprof))
!INTF_END
!     End of subroutine arguments
!       Local scalars:
  INTEGER(KIND=jpim) :: i, chan, prof
  REAL   (KIND=jprb) :: sincsi
  REAL   (KIND=jprb) :: coscsi
  COMPLEX(KIND=jprb) :: crefpar
  COMPLEX(KIND=jprb) :: crefparj
  COMPLEX(KIND=jprb) :: crefperp
  COMPLEX(KIND=jprb) :: crefperpj
  COMPLEX(KIND=jprb) :: sincsi1
  COMPLEX(KIND=jprb) :: coscsi1
  REAL   (KIND=jprb) :: sincsi_ad
  REAL   (KIND=jprb) :: coscsi_ad
  COMPLEX(KIND=jprb) :: crefpar_ad
  COMPLEX(KIND=jprb) :: crefparj_ad
  COMPLEX(KIND=jprb) :: crefperp_ad
  COMPLEX(KIND=jprb) :: crefperpj_ad
  COMPLEX(KIND=jprb) :: sincsi1_ad
  COMPLEX(KIND=jprb) :: coscsi1_ad
  INTEGER(KIND=jpim) :: nchannels            ! Number of radiances computed (channels used * profiles)
!       Local arrays:
  COMPLEX(KIND=jprb) :: waopc(size(chanprof))
  REAL   (KIND=JPRB) :: ZHOOK_HANDLE
!-----End of header-------------------------------------------------------------
  IF (LHOOK) CALL DR_HOOK('RTTOV_FRESNEL_AD', 0_jpim, ZHOOK_HANDLE)
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
        coscsi = cos(sunglint%s(prof)%omega)
        sincsi = sin(sunglint%s(prof)%omega)
        sincsi1                   = sincsi / waopc(i)
        coscsi1                   = sqrt(1._jprb - sincsi1 ** 2_jpim)
        crefpar                   =  - (waopc(i) * coscsi - coscsi1) / (waopc(i) * coscsi + coscsi1)
        crefparj                  = conjg(crefpar)
        crefperp                  = (coscsi - waopc(i) * coscsi1) / (coscsi + waopc(i) * coscsi1)
        crefperpj                 = conjg(crefperp)
!-------------------------------------------------------------------------------
        crefpar_ad                = fresnrefl_ad(i) * 0.5_jprb * crefparj
        crefparj_ad               = fresnrefl_ad(i) * 0.5_jprb * crefpar
        crefperp_ad               = fresnrefl_ad(i) * 0.5_jprb * crefperpj
        crefperpj_ad              = fresnrefl_ad(i) * 0.5_jprb * crefperp
        crefperp_ad               = crefperp_ad + conjg(crefperpj_ad)
        crefperpj_ad              = 0
        coscsi_ad                 = REAL (crefperp_ad / (coscsi + waopc(i) * coscsi1), jprb)
        coscsi1_ad                =  - crefperp_ad * waopc(i) / (coscsi + waopc(i) * coscsi1)
        coscsi_ad                 = REAL (coscsi_ad - crefperp_ad * (coscsi - waopc(i) * coscsi1) / &
                                          (coscsi + waopc(i) * coscsi1) ** 2_jpim, jprb)
        coscsi1_ad                = coscsi1_ad - crefperp_ad * waopc(i) * (coscsi - waopc(i) * coscsi1) / &
                                    (coscsi + waopc(i) * coscsi1) ** 2_jpim
        crefperp_ad               = 0
        crefpar_ad                = crefpar_ad + conjg(crefparj_ad)
        crefparj_ad               = 0
        coscsi_ad                 = REAL (coscsi_ad - crefpar_ad * waopc(i) / (waopc(i) * coscsi + coscsi1), jprb)
        coscsi1_ad                = coscsi1_ad + crefpar_ad / (waopc(i) * coscsi + coscsi1)
        coscsi_ad                 = REAL (coscsi_ad + crefpar_ad * waopc(i) * ( - coscsi1 + waopc(i) * coscsi) / &
                                          (waopc(i) * coscsi + coscsi1) ** 2_jpim, jprb)
        coscsi1_ad                = coscsi1_ad + crefpar_ad * ( - coscsi1 + waopc(i) * coscsi) / &
                                    (waopc(i) * coscsi + coscsi1) ** 2_jpim
        crefpar_ad                = 0
        sincsi1_ad                =  - coscsi1_ad * sincsi1 / coscsi1
        coscsi1_ad                = 0
        sincsi_ad                 = REAL (sincsi1_ad / waopc(i), jprb)
        sincsi1_ad                = 0
        sunglint_ad%s(prof)%omega = sunglint_ad%s(prof)%omega + sincsi_ad * coscsi
        sunglint_ad%s(prof)%omega = sunglint_ad%s(prof)%omega - coscsi_ad * sincsi
        sincsi_ad                 = 0
        coscsi_ad                 = 0
      ENDIF
    ENDIF
  ENDDO
!End of channel loop
  IF (LHOOK) CALL DR_HOOK('RTTOV_FRESNEL_AD', 1_jpim, ZHOOK_HANDLE)
END SUBROUTINE rttov_fresnel_ad
