!
SUBROUTINE rttov_calcrad( &
            & addcosmic,   &
            & chanprof,    &
            & profiles,    &
            & coeffs,      &
            & auxrad)
! Description:
! To convert an array of atmospheric temperatures
!   to planck radiances in many channels
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
! 1.3    11/02/2005  Code vectorisation improved (D Dent)
! 1.4    26/01/2007  Removed polarisation (R Saunders)
! 1.5    15/04/2009  User defined ToA. Layers distinct from levels (P.Rayer)
! 1.6    03/11/2009  Transmittances on levels (A Geer)
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
! Imported Parameters:
! Imported Type Definitions:
  USE rttov_types, ONLY : rttov_chanprof, rttov_coef, profile_Type, radiance_aux
  USE parkind1, ONLY : jplm
  
!INTF_OFF
  USE parkind1, ONLY : jprb, jpim
  USE rttov_const, ONLY : tcosmic
  USE rttov_math_mod, ONLY : PLANCK
!INTF_ON
  IMPLICIT NONE
!subroutine arguments:
  TYPE(profile_Type  ), INTENT(IN)    :: profiles(:) ! profiles
  TYPE(rttov_coef    ), INTENT(IN)    :: coeffs      ! coefficients (Planck)
  TYPE(rttov_chanprof), INTENT(IN)    :: chanprof(:) ! Array of channel indices.
  LOGICAL(KIND=jplm)  , INTENT(IN)    :: addcosmic   ! switch for adding CMB
  TYPE(radiance_aux ),  INTENT(INOUT) :: auxrad      ! auxilary profile info

!INTF_END
! radiances are expressed in mw/cm-1/ster/sq.m
! and temperatures in Kelvin
!local variables:
  INTEGER(KIND=jpim) :: chan, prof, ichan ! loop indices
  INTEGER(KIND=jpim) :: nchannels ! Number of radiances computed
!- End of header ------------------------------------------------------

  nchannels = SIZE(chanprof)
  
  DO ichan = 1, nchannels
     chan = chanprof(ichan)%chan
     prof = chanprof(ichan)%prof

! Calculate effective temperatures, store in auxrad
     IF(coeffs%ff_val_bc) THEN
       auxrad%skin_t_eff(ichan)  = coeffs%ff_bco(chan) + coeffs%ff_bcs(chan) * &
                                   profiles(prof)%skin%t
       auxrad%surf_t_eff(ichan)  = coeffs%ff_bco(chan) + coeffs%ff_bcs(chan) * &
                                   profiles(prof)%s2m%t
       auxrad%air_t_eff(:,ichan) = coeffs%ff_bco(chan) + coeffs%ff_bcs(chan) * &
                                   profiles(prof)%t(:)
     ELSE
       auxrad%skin_t_eff(ichan)  = profiles(prof)%skin%t
       auxrad%surf_t_eff(ichan)  = profiles(prof)%s2m%t
       auxrad%air_t_eff(:,ichan) = profiles(prof)%t(:)
     ENDIF
     
     IF (addcosmic) THEN
       auxrad%cosmic_t_eff(ichan) = tcosmic
       IF(coeffs%ff_val_bc) auxrad%cosmic_t_eff(ichan) = coeffs%ff_bco(chan) + &
         coeffs%ff_bcs(chan) * tcosmic

       CALL PLANCK(coeffs%planck1(chan), coeffs%planck2(chan), &
                   auxrad%cosmic_t_eff(ichan), auxrad%cosmic(ichan))
     ELSE
        auxrad%cosmic(ichan) = 0.0_JPRB
     ENDIF
        
     CALL PLANCK(coeffs%planck1(chan), coeffs%planck2(chan), &
                 auxrad%skin_t_eff(ichan), auxrad%skin(ichan))

     CALL PLANCK(coeffs%planck1(chan), coeffs%planck2(chan), &
                 auxrad%surf_t_eff(ichan), auxrad%surfair(ichan))

     CALL PLANCK(coeffs%planck1(chan), coeffs%planck2(chan), &
                 auxrad%air_t_eff(:,ichan), auxrad%air(:,ichan))
   ENDDO
END SUBROUTINE rttov_calcrad
