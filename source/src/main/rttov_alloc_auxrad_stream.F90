SUBROUTINE rttov_alloc_auxrad_stream( &
            & ERR,           &
            & auxrad_stream, &
            & opts,          &
            & nstreams,      &
            & nlayers,       &
            & nchannels,     &
            & asw,           &
            & init)
! Description:
!   Allocates/deallocates the auxiliary radiance structure 
!   for streams
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
! Method:
!
! Current Code Owner: SAF NWP
!
! History:
! Version   Date     Comment
! -------   ----     -------
! 2010/03 Code cleaning Pascal Brunel, Philippe Marguinaud
!
!
! Code Description:
!   Language:           Fortran 90.
!   Software Standards: "European Standards for Writing and
!     Documenting Exchangeable Fortran 90 Code".

  USE parkind1, ONLY : jpim, jplm
  USE rttov_types, ONLY : rttov_options, radiance_aux
!INTF_OFF
#include "throw.h"
!INTF_ON
  IMPLICIT NONE
  INTEGER(KIND=jpim) , INTENT(OUT)             :: ERR
  TYPE(radiance_aux ), INTENT(INOUT)           :: auxrad_stream
  TYPE(rttov_options), INTENT(IN)              :: opts
  INTEGER(KIND=jpim) , INTENT(IN)              :: nstreams
  INTEGER(KIND=jpim) , INTENT(IN)              :: nlayers
  INTEGER(KIND=jpim) , INTENT(IN)              :: nchannels
  INTEGER(KIND=jpim) , INTENT(IN)              :: asw
  LOGICAL(KIND=jplm) , INTENT(IN)   , OPTIONAL :: init
!INTF_END
#include "rttov_errorreport.interface"
#include "rttov_init_auxrad_stream.interface"
  LOGICAL(KIND=jplm) :: init1
  TRY
  init1 = .FALSE.
  IF (Present(init)) init1 = init
  IF (asw .EQ. 1) THEN
    ALLOCATE (auxrad_stream%up(nlayers, 0:nstreams, nchannels),       &
            & auxrad_stream%down(nlayers, 0:nstreams, nchannels),     &
            & auxrad_stream%down_ref(nlayers, 0:nstreams, nchannels), &
            & auxrad_stream%cloudy(0:nstreams, nchannels),            &
            & auxrad_stream%meanrad_up(0:nstreams, nchannels),        &
            & auxrad_stream%meanrad_down(0:nstreams, nchannels), STAT = ERR)
    THROWM( ERR .NE. 0 , "allocation of auxrad_stream")

    NULLIFY (auxrad_stream%up_solar,           &
           & auxrad_stream%down_solar,         &
           & auxrad_stream%meanrad_up_solar,   &
           & auxrad_stream%meanrad_down_solar, &
           & auxrad_stream%down_ref_solar)
    NULLIFY (auxrad_stream%FAC1_2, auxrad_stream%FAC2_2, &
           & auxrad_stream%FAC3_2, auxrad_stream%FAC4_2, auxrad_stream%FAC5_2, &
           & auxrad_stream%FAC6_2, auxrad_stream%FAC7_2, auxrad_stream%FAC1_3, &
           & auxrad_stream%FAC2_3, auxrad_stream%FAC3_3, auxrad_stream%FAC4_3, &
           & auxrad_stream%FAC5_3, auxrad_stream%FAC6_3, auxrad_stream%FAC7_3)

    IF (opts%rt_ir%addsolar) THEN
      ALLOCATE (auxrad_stream%up_solar(nlayers, 0:nstreams, nchannels),   &
              & auxrad_stream%meanrad_up_solar(0:nstreams, nchannels), STAT = ERR)
      THROWM( ERR .NE. 0 , "allocation of auxrad_stream up solar radiance arrays")
      ALLOCATE (auxrad_stream%down_solar(nlayers, 0:nstreams, nchannels), &
              & auxrad_stream%meanrad_down_solar(0:nstreams, nchannels),  &
              & auxrad_stream%down_ref_solar(nlayers, 0:nstreams, nchannels), STAT = ERR)
      THROWM( ERR .NE. 0 , "allocation of auxrad_stream down solar radiance arrays")
      IF (opts%rt_ir%addaerosl .OR. opts%rt_ir%addclouds) THEN
        ALLOCATE (auxrad_stream%FAC1_2(nlayers, 0:nstreams, nchannels), &
                & auxrad_stream%FAC2_2(nlayers, 0:nstreams, nchannels), &
                & auxrad_stream%FAC3_2(nlayers, 0:nstreams, nchannels), &
                & auxrad_stream%FAC4_2(nlayers, 0:nstreams, nchannels), &
                & auxrad_stream%FAC5_2(nlayers, 0:nstreams, nchannels), &
                & auxrad_stream%FAC6_2(nlayers, 0:nstreams, nchannels), &
                & auxrad_stream%FAC7_2(nlayers, 0:nstreams, nchannels), &
                & auxrad_stream%FAC1_3(0:nstreams, nchannels),          &
                & auxrad_stream%FAC2_3(0:nstreams, nchannels),          &
                & auxrad_stream%FAC3_3(0:nstreams, nchannels),          &
                & auxrad_stream%FAC4_3(0:nstreams, nchannels),          &
                & auxrad_stream%FAC5_3(0:nstreams, nchannels),          &
                & auxrad_stream%FAC6_3(0:nstreams, nchannels),          &
                & auxrad_stream%FAC7_3(0:nstreams, nchannels), STAT = ERR)
        THROWM( ERR .NE. 0 , "deallocation of auxrad_stream%FAC")
      ENDIF
    ENDIF
    IF (init1) THEN
      CALL rttov_init_auxrad_stream(auxrad_stream)
    ENDIF
  ENDIF

  IF (asw .EQ. 0) THEN
    DEALLOCATE (auxrad_stream%up, auxrad_stream%down, auxrad_stream%down_ref, auxrad_stream%cloudy, &
                auxrad_stream%meanrad_up, auxrad_stream%meanrad_down, STAT = ERR)
    THROWM( ERR .NE. 0 , "deallocation of auxrad_stream")

    IF (opts%rt_ir%addsolar) THEN
      DEALLOCATE (auxrad_stream%up_solar,           &
                  auxrad_stream%meanrad_up_solar, STAT = ERR)
      THROWM( ERR .NE. 0 , "deallocation of auxrad_stream up solar radiance arrays")

      DEALLOCATE (auxrad_stream%down_solar,         &
                & auxrad_stream%meanrad_down_solar, &
                & auxrad_stream%down_ref_solar, STAT = ERR)
      THROWM( ERR .NE. 0 , "deallocation of auxrad_stream down solar radiance arrays")
      
      IF (opts%rt_ir%addaerosl .OR. opts%rt_ir%addclouds) THEN
        DEALLOCATE (auxrad_stream%FAC1_2, auxrad_stream%FAC2_2, &
                  & auxrad_stream%FAC3_2, auxrad_stream%FAC4_2, auxrad_stream%FAC5_2, &
                  & auxrad_stream%FAC6_2, auxrad_stream%FAC7_2, auxrad_stream%FAC1_3, &
                  & auxrad_stream%FAC2_3, auxrad_stream%FAC3_3, auxrad_stream%FAC4_3, &
                  & auxrad_stream%FAC5_3, auxrad_stream%FAC6_3, auxrad_stream%FAC7_3, STAT = ERR)
        THROWM( ERR .NE. 0 , "deallocation of auxrad_stream%FAC")
      ENDIF
    ENDIF

    NULLIFY (auxrad_stream%up, auxrad_stream%down, auxrad_stream%down_ref, auxrad_stream%cloudy)
    NULLIFY (auxrad_stream%meanrad_up, auxrad_stream%meanrad_down)
    NULLIFY (auxrad_stream%up_solar,           &
           & auxrad_stream%down_solar,         &
           & auxrad_stream%meanrad_up_solar,   &
           & auxrad_stream%meanrad_down_solar, &
           & auxrad_stream%down_ref_solar)
    NULLIFY (auxrad_stream%FAC1_2, auxrad_stream%FAC2_2, &
           & auxrad_stream%FAC3_2, auxrad_stream%FAC4_2, auxrad_stream%FAC5_2, &
           & auxrad_stream%FAC6_2, auxrad_stream%FAC7_2, auxrad_stream%FAC1_3, &
           & auxrad_stream%FAC2_3, auxrad_stream%FAC3_3, auxrad_stream%FAC4_3, &
           & auxrad_stream%FAC5_3, auxrad_stream%FAC6_3, auxrad_stream%FAC7_3)
  ENDIF
  CATCH
END SUBROUTINE 
