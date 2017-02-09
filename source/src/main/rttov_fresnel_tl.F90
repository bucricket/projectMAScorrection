!     Compute sea surface emissivity for radiative transfer calculations
SUBROUTINE rttov_fresnel_tl( &
            & chanprof,     &
            & profiles,     &
            & solar,        &
            & coef,         &
            & sunglint,     &
            & sunglint_tl,  &
            & fresnrefl,    &
            & fresnrefl_tl)
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
!     2            14/03/2002  M. Matricardi. ECMWF.
!                              Rewritten in free format FORTRAN90.
!     3            15/07/2003  RTIASI-4. M. Matricardi. ECMWF.
!                              Rewritten to include the computation of the
!                              Fresnel reflectivity for a flat watersurface
!     4            05/02/2007  Removed polarisation index (R Saunders)
!     5            05/07/2010  Remove addsolar flag from profiles structure (J Hocking)
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
  TYPE(rttov_chanprof), INTENT(IN)  :: chanprof(:)
  TYPE(profile_type  ), INTENT(IN)  :: profiles(:)
  LOGICAL(KIND=jplm)  , INTENT(IN)  :: solar(size(chanprof))
  TYPE(rttov_coef    ), INTENT(IN)  :: coef
  TYPE(sunglint_type ), INTENT(IN)  :: sunglint
  TYPE(sunglint_type ), INTENT(IN)  :: sunglint_tl
  REAL(KIND=jprb)     , INTENT(IN)  :: fresnrefl   (size(chanprof))
  REAL(KIND=jprb)     , INTENT(OUT) :: fresnrefl_tl(size(chanprof))
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
  REAL   (KIND=jprb) :: sincsi_tl
  REAL   (KIND=jprb) :: coscsi_tl
  COMPLEX(KIND=jprb) :: crefpar_tl
  COMPLEX(KIND=jprb) :: crefparj_tl
  COMPLEX(KIND=jprb) :: crefperp_tl
  COMPLEX(KIND=jprb) :: crefperpj_tl
  COMPLEX(KIND=jprb) :: sincsi1_tl
  COMPLEX(KIND=jprb) :: coscsi1_tl
  INTEGER(KIND=jpim) :: nchannels            ! Number of radiances computed (channels used * profiles)
!       Local arrays
  COMPLEX(KIND=jprb) :: waopc(size(chanprof))
  REAL   (KIND=JPRB) :: ZHOOK_HANDLE
!-----End of header-------------------------------------------------------------
!-------------------------------------------------------------------------------
!         1.   SETS UP SURFACE EMISSIVITIES.
!-------------------------------------------------------------------------------
  IF (LHOOK) CALL DR_HOOK('RTTOV_FRESNEL_TL', 0_jpim, ZHOOK_HANDLE)
  nchannels = size(chanprof)
  DO i = 1, nchannels!Loop over channels.
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
        coscsi          = cos(sunglint%s(prof)%omega)
        sincsi          = sin(sunglint%s(prof)%omega)
        sincsi1         = sincsi / waopc(i)
        coscsi1         = sqrt(1._jprb - sincsi1 ** 2_jpim)
        crefpar         =  - (waopc(i) * coscsi - coscsi1) / (waopc(i) * coscsi + coscsi1)
        crefparj        = conjg(crefpar)
        crefperp        = (coscsi - waopc(i) * coscsi1) / (coscsi + waopc(i) * coscsi1)
        crefperpj       = conjg(crefperp)
        coscsi_tl       =  - sincsi * sunglint_tl%s(prof)%omega
        sincsi_tl       = coscsi * sunglint_tl%s(prof)%omega
        sincsi1_tl      = sincsi_tl / waopc(i)
        coscsi1_tl      =  - sincsi1_tl * sincsi1 / coscsi1
        crefpar_tl      =                                                             &
          &  - (waopc(i) * coscsi_tl - coscsi1_tl) / (waopc(i) * coscsi + coscsi1) +  &
          & (waopc(i) * coscsi_tl + coscsi1_tl) *                                     &
          & (waopc(i) * coscsi - coscsi1) / (waopc(i) * coscsi + coscsi1) ** 2_jpim
        crefparj_tl     = conjg(crefpar_tl)
        crefperp_tl     =                                                          &
          & (coscsi_tl - waopc(i) * coscsi1_tl) / (coscsi + waopc(i) * coscsi1) -  &
          & (coscsi_tl + waopc(i) * coscsi1_tl) *                                  &
          & (coscsi - waopc(i) * coscsi1) / (coscsi + waopc(i) * coscsi1) ** 2_jpim
        crefperpj_tl    = conjg(crefperp_tl)
        fresnrefl_tl(i) = REAL (crefpar_tl * 0.5_jprb * crefparj + crefparj_tl * 0.5_jprb * crefpar + &
                                crefperp_tl * 0.5_jprb * crefperpj + crefperpj_tl * 0.5_jprb * crefperp, jprb)
      ENDIF
    ENDIF
  ENDDO
!End of channel loop
  IF (LHOOK) CALL DR_HOOK('RTTOV_FRESNEL_TL', 1_jpim, ZHOOK_HANDLE)
END SUBROUTINE rttov_fresnel_tl
