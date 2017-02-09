SUBROUTINE rttov_alloc_aux_prof( &
            & err,       &
            & nprofiles, &
            & nlevels,   &
            & id_sensor, &
            & aux_prof,  &
            & opts,      &
            & asw,       &
            & init,      &
            & alloc_layer_vars)
! Description:
!   Allocates/deallocates the auxiliary profile structure
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
!INTF_OFF
#include "throw.h"
!INTF_ON
  USE rttov_types, ONLY : rttov_options, profile_aux
  USE parkind1, ONLY : jpim, jplm
!INTF_OFF
  USE rttov_const, ONLY : sensor_id_mw, sensor_id_po, naer_max
!INTF_ON
  INTEGER(KIND=jpim) , INTENT(OUT)               :: err
  INTEGER(KIND=jpim) , INTENT(IN)                :: nprofiles
  INTEGER(KIND=jpim) , INTENT(IN)                :: nlevels
  INTEGER(KIND=jpim) , INTENT(IN)                :: id_sensor
  TYPE(profile_aux  ), INTENT(INOUT)             :: aux_prof
  TYPE(rttov_options), INTENT(IN)                :: opts
  INTEGER(KIND=jpim) , INTENT(IN)                :: asw
  LOGICAL(KIND=jpim) , OPTIONAL     , INTENT(IN) :: init
  LOGICAL(KIND=jplm) , OPTIONAL     , INTENT(IN) :: alloc_layer_vars
!INTF_END

#include "rttov_errorreport.interface"
#include "rttov_init_aux_prof.interface"

  LOGICAL(KIND=jplm) :: init1
  LOGICAL(KIND=jplm) :: alloc_layer_vars1
  INTEGER(KIND=jpim) :: nlayers

  TRY
  nlayers = nlevels - 1
  init1 = .FALSE.
  alloc_layer_vars1 = .FALSE.

  IF (PRESENT(alloc_layer_vars)) alloc_layer_vars1 = alloc_layer_vars
  IF (PRESENT(init)) init1 = init

  IF (asw .EQ. 1) THEN
    CALL nullify_struct()

    ALLOCATE (aux_prof%s(nprofiles), STAT = err)
    THROWM(err .NE. 0, "allocation of aux_prof%s")
    IF ((id_sensor == sensor_id_mw .OR. id_sensor == sensor_id_po) .AND. opts%rt_mw%clw_data) THEN
      ALLOCATE (aux_prof%debye_prof(5, nlevels, nprofiles), STAT = err)
      THROWM(err .NE. 0, "allocation of aux_prof%debye_prof")
    ENDIF
    IF (opts%rt_ir%addclouds .AND. .NOT. opts%rt_ir%user_cld_opt_param) THEN
      ALLOCATE (aux_prof%dg(nlevels, nprofiles), STAT = err)
      THROWM(err .NE. 0, "allocation of aux_prof%dg")
      ALLOCATE (aux_prof%fac1_dg(nlevels, nprofiles), STAT = err)
      THROWM(err .NE. 0, "allocation of aux_prof%fac1_dg")
      ALLOCATE (aux_prof%fac2_dg(nlevels, nprofiles), STAT = err)
      THROWM(err .NE. 0, "allocation of aux_prof%fac2_dg")
      ALLOCATE (aux_prof%fac3_dg(nlevels, nprofiles), STAT = err)
      THROWM(err .NE. 0, "allocation of aux_prof%fac3_dg")
    ENDIF
    IF (opts%rt_ir%addaerosl .AND. .NOT. opts%rt_ir%user_aer_opt_param) THEN
      ALLOCATE (aux_prof%iaernum(nlevels, nprofiles), STAT = err)
      THROWM(err .NE. 0, "allocation of aux_prof%iaernum")
      ALLOCATE (aux_prof%iaertyp(naer_max, nlevels, nprofiles), STAT = err)
      THROWM(err .NE. 0, "allocation of aux_prof%iaertyp")
      ALLOCATE (aux_prof%relhum(nlevels, nprofiles), STAT = err)
      THROWM(err .NE. 0, "allocation of aux_prof%relhum")
      ALLOCATE (aux_prof%relhumref(nlevels, nprofiles), STAT = err)
      THROWM(err .NE. 0, "allocation of aux_prof%relhumref")
    ENDIF

    IF (alloc_layer_vars1) THEN
      ALLOCATE(&
        aux_prof%t_layer(nlayers, nprofiles), aux_prof%w_layer(nlayers, nprofiles), &
        aux_prof%o3_layer(nlayers, nprofiles), aux_prof%dt(nlayers, nprofiles), &
        aux_prof%dto(nlayers, nprofiles), aux_prof%tr(nlayers, nprofiles), aux_prof%tr_r(nlayers, nprofiles), &
        aux_prof%tw_sqrt(nlayers, nprofiles), aux_prof%wr_sqrt(nlayers, nprofiles), &
        aux_prof%wr_rsqrt(nlayers, nprofiles), &
        aux_prof%tw_4rt(nlayers, nprofiles), aux_prof%ww_r(nlayers, nprofiles), &
        aux_prof%ow_rsqrt(nlayers, nprofiles), aux_prof%ow_r(nlayers, nprofiles),& 
        aux_prof%or_sqrt(nlayers, nprofiles), aux_prof%ow_sqrt(nlayers, nprofiles), &
        aux_prof%wr(nlayers, nprofiles), aux_prof%or(nlayers, nprofiles), &
        aux_prof%tw(nlayers, nprofiles), aux_prof%ww(nlayers, nprofiles), &
        aux_prof%ow(nlayers, nprofiles), aux_prof%sum(nlayers, 2), STAT = err)
      THROWM(err .NE. 0, "allocation of aux_prof layer quantities")
    ENDIF

    IF (init1) THEN
      CALL rttov_init_aux_prof(aux_prof)
    ENDIF
  ENDIF
  IF (asw .EQ. 0) THEN
    DEALLOCATE (aux_prof%s, STAT = err)
    THROWM(err .NE. 0, "deallocation of aux_prof%s")
    IF ((id_sensor == sensor_id_mw .OR. id_sensor == sensor_id_po) .AND. opts%rt_mw%clw_data) THEN
      DEALLOCATE (aux_prof%debye_prof, STAT = err)
      THROWM(err .NE. 0, "deallocation of aux_prof%debye_prof")
    ENDIF
    IF (opts%rt_ir%addclouds .AND. .NOT. opts%rt_ir%user_cld_opt_param) THEN
      DEALLOCATE (aux_prof%dg, STAT = err)
      THROWM(err .NE. 0, "deallocation of aux_prof%dg")
      DEALLOCATE (aux_prof%fac1_dg, STAT = err)
      THROWM(err .NE. 0, "deallocation of aux_prof%fac1_dg")
      DEALLOCATE (aux_prof%fac2_dg, STAT = err)
      THROWM(err .NE. 0, "deallocation of aux_prof%fac2_dg")
      DEALLOCATE (aux_prof%fac3_dg, STAT = err)
      THROWM(err .NE. 0, "deallocation of aux_prof%fac3_dg")
    ENDIF
    IF (opts%rt_ir%addaerosl .AND. .NOT. opts%rt_ir%user_aer_opt_param) THEN
      DEALLOCATE (aux_prof%iaernum, STAT = err)
      THROWM(err .NE. 0, "deallocation of aux_prof%iaernum")
      DEALLOCATE (aux_prof%iaertyp, STAT = err)
      THROWM(err .NE. 0, "deallocation of aux_prof%iaertyp")
      DEALLOCATE (aux_prof%relhum, STAT = err)
      THROWM(err .NE. 0, "deallocation of aux_prof%relhum")
      DEALLOCATE (aux_prof%relhumref, STAT = err)
      THROWM(err .NE. 0, "deallocation of aux_prof%relhumref")
    ENDIF

    IF (alloc_layer_vars1) THEN
      DEALLOCATE(&
        aux_prof%t_layer, aux_prof%w_layer, aux_prof%o3_layer, aux_prof%dt, &
        aux_prof%dto, aux_prof%tr, aux_prof%tr_r, aux_prof%tw_sqrt, &
        aux_prof%wr_sqrt, aux_prof%wr_rsqrt, aux_prof%tw_4rt, aux_prof%ww_r, aux_prof%ow_rsqrt,&
        aux_prof%ow_r, aux_prof%or_sqrt, aux_prof%ow_sqrt, aux_prof%wr, aux_prof%or, aux_prof%tw, &
        aux_prof%ww, aux_prof%ow, aux_prof%sum, STAT = err)
      THROWM(err .NE. 0, "deallocation of aux_prof layer quantities")
    ENDIF

    CALL nullify_struct()
  ENDIF
  CATCH
CONTAINS
  SUBROUTINE nullify_struct()
    NULLIFY (aux_prof%s)
    NULLIFY (aux_prof%debye_prof)
    NULLIFY (aux_prof%dg)
    NULLIFY (aux_prof%fac1_dg)
    NULLIFY (aux_prof%fac2_dg)
    NULLIFY (aux_prof%fac3_dg)
    NULLIFY (aux_prof%iaernum)
    NULLIFY (aux_prof%iaertyp)
    NULLIFY (aux_prof%relhum)
    NULLIFY (aux_prof%relhumref)

    NULLIFY(&
      aux_prof%t_layer, aux_prof%w_layer, aux_prof%o3_layer, aux_prof%dt, &
      aux_prof%dto, aux_prof%tr, aux_prof%tr_r, aux_prof%tw_sqrt, &
      aux_prof%wr_sqrt, aux_prof%wr_rsqrt, aux_prof%tw_4rt, aux_prof%ww_r, aux_prof%ow_rsqrt,&
      aux_prof%ow_r, aux_prof%or_sqrt, aux_prof%ow_sqrt, aux_prof%wr, aux_prof%or, aux_prof%tw, &
      aux_prof%ww, aux_prof%ow, aux_prof%sum)
  END SUBROUTINE nullify_struct
END SUBROUTINE 
