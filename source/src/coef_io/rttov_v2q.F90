SUBROUTINE rttov_v2q (&
       & h2o_unit,  &! in
       & h2o,       &! in
       & gas_id,    &! in
       & v_gas,     &! in
       & q_gas     ) ! inout
  !
  ! Description:
  !   Conversion of volume mixing ratio (wet air) to specific concentration (wet air).
  !   Gases are defined by the "gas_id_xxx" codes in the rttov_const module.
  !   Method uses an equivalent molecular weight of wet air.
  !
  !   Water vapour value (h2o) must be in units with respect to moist air.
  !   To convert from dry air value: wv_wet = wv_dry / (1 + wv_dry)
  !   for water vapour wv in units of kg/kg or ppmv wrt dry air.
  !
  ! Copyright:
  !
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

  USE parkind1, ONLY : jpim, jprb

!INTF_OFF
  USE rttov_const, ONLY : &
      ngases_max,&
      mair,&
      mh2o,&
      gas_mass,&
      gas_id_mixed,&
      gas_unit_specconc,&
      gas_unit_ppmv
!INTF_ON

  IMPLICIT NONE

  INTEGER(KIND=jpim), INTENT(IN)    :: h2o_unit ! Water vapour input unit
                                                ! 1 = specific concent. (kg/kg over wet air)
                                                ! 2 = volume mixing ratio (ppmv over wet air)
                                                ! (see gas unit id codes in module rttov_const)
  REAL(KIND=jprb),    INTENT(IN)    :: h2o      ! Water Vapour content in unit h2o_unit

  INTEGER(KIND=jpim), INTENT(IN)    :: gas_id   ! Gas identification number
                                                ! (see gas id codes in module rttov_const)
  REAL(KIND=jprb),    INTENT(IN)    :: v_gas    ! volume mixing ratio for gas (ppmv)
  REAL(KIND=jprb),    INTENT(INOUT) :: q_gas    ! specific concentration for gas (kg/kg)
!INTF_END

  REAL(KIND=jprb), PARAMETER :: eps = mh2o / mair

  REAL(KIND=jprb) :: Mwet  ! equivalent molecular weight of wet air (g)
  REAL(KIND=jprb) :: v_h2o ! volume mixing ratio for Water Vapour (v/v)

  !- End of header --------------------------------------------------------

  ! Calculate volume mixing ratio (no unit: v/v) for Water Vapour
  IF (h2o_unit == gas_unit_specconc) THEN
     v_h2o = h2o / (eps * (1 - h2o) + h2o)
  ELSE IF (h2o_unit == gas_unit_ppmv) THEN
     v_h2o = h2o * 1.E-06_jprb
  ELSE
     v_h2o = 0._jprb
  ENDIF

  ! Humid air molar mass
  Mwet = (1 - v_h2o) * Mair + v_h2o * Mh2o

  ! Calculate specific concentration for gas (kg/kg)
  IF (gas_id == gas_id_mixed) THEN
    ! keep same value for Mixed gases
    q_gas = v_gas

  ELSE IF (gas_id > 0 .AND. gas_id < ngases_max) THEN
    q_gas = v_gas * 1.E-06_jprb * gas_mass(gas_id) / Mwet

  ELSE
    q_gas = 0._jprb

  ENDIF

END SUBROUTINE rttov_v2q
