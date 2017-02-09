SUBROUTINE rttov_calcrad_k( &
            & chanprof,       &
            & profiles,       &
            & profiles_k,    &
            & coeffs,         &
            & auxrad,         &
            & rad_skin_k,    &
            & rad_surfair_k, &
            & rad_air_k)
!
! Description:
! AD code to convert an array of atmospheric temperatures
!   to planck radiances in many channels
! No AD on Rad Cosmic, ad = 0.
!
! derivative of Planck function with respect to temperature is
!
!                                     C2 * Nu
!              C1 * C2 * Nu**4 * Exp( ------- )
!                                        T
! B'(T,Nu) = ------------------------------------- dT
!                     (      C2 * Nu       )**2
!               T**2 *( Exp( ------- ) - 1 )
!                     (         T          )
!
!
! which can be reduced to the following, with
!  C1 = C1 * Nu**3
!  C2 = C2 * Nu
!
!              C2 * B(T,Nu) * (C1 + B(T,Nu))
!  B'(T,Nu) =  ----------------------------- dT
!                        C1 * T**2
!
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
!    Copyright 2002, EUMETSAT, All Rights Reserved.
!
! Method: Uses band correction factors to convert T to radiance
!         which have been precomputed for each channel and are read from
!         the RT coefficient file.
!
! Current Code Owner: SAF NWP
!
! History:
! Version   Date     Comment
! -------   ----     -------
! 1.0    01/12/2002  New F90 code with structures (P Brunel A Smith)
!                      based on PLNCX of previous RTTOV versions
! 1.1    02/01/2003  Comments added (R Saunders)
! 1.2    26/09/2003  Multiple polarisations (S English)
! 1.3    29/03/2005  Add end of header comment (J. Cameron)
! 1.4    05/02/2007  Removed polarisation index (R Saunders)
! 1.5    15/08/2009  User defined ToA. Layers distinct from levels (P.Rayer)
! 1.6    03/11/2009  Transmittances / optical depths on levels (A Geer)
!
! 2010/03 Code cleaning Pascal Brunel, Philippe Marguinaud
!
!
! Code Description:
!   Language:           Fortran 90.
!   Software Standards: "European Standards for Writing and
!     Documenting Exchangeable Fortran 90 Code".
!
! Declarations:
! Modules used:
  USE rttov_types, ONLY : rttov_chanprof, rttov_coef, profile_type, radiance_aux
  USE parkind1, ONLY : jprb
!INTF_OFF
  USE parkind1, ONLY : jpim, jplm
  USE rttov_math_mod, ONLY : PLANCK_AD
!INTF_ON
  IMPLICIT NONE
!subroutine arguments:
  TYPE(rttov_chanprof), INTENT(IN)    :: chanprof      (:)
  TYPE(profile_Type  ), INTENT(IN)    :: profiles  (:)
  TYPE(profile_Type  ), INTENT(INOUT) :: profiles_k(:)
  TYPE(rttov_coef    ), INTENT(IN)    :: coeffs
  TYPE(radiance_aux ),  INTENT(IN)    :: auxrad    ! auxilary profile info
  REAL(KIND=jprb)     , INTENT(IN)    :: rad_skin_k   (size(chanprof))
  REAL(KIND=jprb)     , INTENT(IN)    :: rad_surfair_k(size(chanprof))
  REAL(KIND=jprb),      INTENT(IN)    :: rad_air_k (size(profiles(1)%t), size(chanprof))
!INTF_END
!local variables:
  REAL   (KIND=jprb) :: t_effective_skin_k, t_effective_s2m_k
  REAL   (KIND=jprb) :: t_effective_air_k(size(profiles(1)%t))
  INTEGER(KIND=jpim) :: ichan, chan ! loop indices
  INTEGER(KIND=jpim) :: nchannels  ! Number of radiances computed
!- End of header --------------------------------------------------------
  nchannels = SIZE(chanprof)
  
  DO ichan = 1, nchannels
     chan = chanprof(ichan)%chan

     CALL PLANCK_AD(coeffs%planck1(chan), coeffs%planck2(chan),&
                    auxrad%skin_t_eff(ichan), t_effective_skin_k, &
                    auxrad%skin(ichan), rad_skin_k(ichan), acc= .FALSE._jplm)

     CALL PLANCK_AD(coeffs%planck1(chan), coeffs%planck2(chan), &
                    auxrad%surf_t_eff(ichan), t_effective_s2m_k, &
                    auxrad%surfair(ichan), rad_surfair_k(ichan), acc= .FALSE._jplm)

     CALL PLANCK_AD(coeffs%planck1(chan), coeffs%planck2(chan),&
                    auxrad%air_t_eff(:,ichan), t_effective_air_k(:), &
                    auxrad%air(:,ichan), rad_air_k(:,ichan), acc= .FALSE._jplm)      

     IF(coeffs%ff_val_bc) THEN
       profiles_k(ichan)%skin%t = profiles_k(ichan)%skin%t + & 
         coeffs%ff_bcs(chan) * t_effective_skin_k
       profiles_k(ichan)%s2m%t = profiles_k(ichan)%s2m%t + &
         coeffs%ff_bcs(chan) * t_effective_s2m_k
       profiles_k(ichan)%t(:) = profiles_k(ichan)%t(:) + &
         coeffs%ff_bcs(chan) * t_effective_air_k(:)
     ELSE
       profiles_k(ichan)%skin%t = profiles_k(ichan)%skin%t + & 
         t_effective_skin_k
       profiles_k(ichan)%s2m%t = profiles_k(ichan)%s2m%t + &
         t_effective_s2m_k
       profiles_k(ichan)%t(:) = profiles_k(ichan)%t(:) + &
         t_effective_air_k(:)
     ENDIF

  ENDDO

END SUBROUTINE rttov_calcrad_k
