SUBROUTINE rttov_alloc_traj_dyn (err, traj_dyn, opts, nchannels, nlayers, nstreams, ncldtyp, asw, direct)
! Description:
!   Allocates/deallocates dynamic trajectory data
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
USE rttov_types, ONLY : &
    & rttov_traj_dyn,   &
    & rttov_options

USE parkind1, ONLY : &
    & jpim, jplm
#include "throw.h"
IMPLICIT NONE

INTEGER(KIND=jpim),   INTENT(OUT)          :: err
TYPE(rttov_traj_dyn), INTENT(INOUT)        :: traj_dyn
TYPE(rttov_options),  INTENT(IN)           :: opts
INTEGER(KIND=jpim),   INTENT(IN)           :: nchannels
INTEGER(KIND=jpim),   INTENT(IN)           :: nlayers
INTEGER(KIND=jpim),   INTENT(IN)           :: nstreams
INTEGER(KIND=jpim),   INTENT(IN)           :: ncldtyp
INTEGER(KIND=jpim),   INTENT(IN)           :: asw
LOGICAL(KIND=jplm),   INTENT(IN), OPTIONAL :: direct

!INTF_END

#include "rttov_alloc_transmission_aux.interface"
#include "rttov_alloc_auxrad_stream.interface"
#include "rttov_alloc_trans_scatt_ir.interface"

TRY

  CALL rttov_alloc_transmission_aux( &
        & err,                           &
        & opts,                          &
        & traj_dyn%transmission_aux,     &
        & nlayers,                       &
        & nchannels,                     &
        & asw,                           &
        & nstreams,                      &
        & direct = direct)

  THROWM(err .NE. 0, "allocation of transmission_aux")

  CALL rttov_alloc_auxrad_stream( &
        & err,                      &
        & traj_dyn%auxrad_stream,   &
        & opts,                     &
        & nstreams,                 &
        & nlayers,                  &
        & nchannels,                &
        & asw)

  THROWM(err .NE. 0, "allocation of auxrad_stream")

  IF (opts%rt_ir%addaerosl .OR. opts%rt_ir%addclouds) THEN
    CALL rttov_alloc_trans_scatt_ir( &
          & err,                                    &
          & opts,                                   &
          & traj_dyn%transmission_scatt_ir_stream,  &
          & nchannels,                              &
          & ncldtyp,                                &
          & nlayers,                                &
          & asw,                                    &
          & nstreams,                               &
          & stream = .TRUE._jplm)

    THROWM(err .NE. 0, "allocation of trans_scatt_ir")
  ENDIF

  traj_dyn%nstreams = nstreams
CATCH

END SUBROUTINE
