Subroutine rttov789_write_ascii_coef  (&
       & errorstatus,   &! out
       & coef,          &! in
       & file_id        )! in
  ! Description:
  ! Write ASCII coef file in v7/9 format.
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
  !    Copyright 2013, EUMETSAT, All Rights Reserved.
  !
  ! Method:
  !
  ! Current Code Owner: SAF NWP
  !
  ! Code Description:
  !   Language:           Fortran 90.
  !   Software Standards: "European Standards for Writing and
  !     Documenting Exchangeable Fortran 90 Code".
  !
  ! Declarations:
  ! Modules used:
#include "throw.h"

  Use rttov_types, Only :  &
       & rttov_coef

  Use parkind1, Only : jpim
!INTF_OFF
  Use rttov_const, Only :   &
       & errorstatus_success ,&
       & errorstatus_fatal   ,&
       & gas_id_mixed        ,&
       & gas_id_watervapour  ,&
       & gas_id_ozone        ,&
       & gas_id_wvcont       ,&
       & gas_id_co2          ,&
       & gas_id_n2o          ,&
       & gas_id_co           ,&
       & gas_id_ch4          ,&
       & gas_name            ,&
       & gas_unit_name       ,&
       & speedl

  Use mod_rttov_fastem3_coef, Only : &
       & fastem3_coef
!INTF_ON
  Implicit None

  Integer(Kind=jpim), Intent (in)  :: file_id      ! file logical unit number
  Type( rttov_coef ), Intent (in)  :: coef         ! coefficients
  Integer(Kind=jpim), Intent (out) :: errorstatus  ! return code

!INTF_END

#include "rttov_errorreport.interface"

  Integer(Kind=jpim) :: i, j, l, k
  Integer            :: io_status
  Integer(Kind=jpim) :: id_comp_lvl
  Character (len=2) :: sensor
  Character (len=32) :: section

  Character (len=80) :: errMessage
  Character (len=16) :: NameOfRoutine = 'rttov_writecoef '
  Integer(Kind=jpim) :: err
  !- End of header --------------------------------------------------------
TRY
  errorstatus     = errorstatus_success

   Write( errMessage, '( "write coefficient to file_id ", i2, " in ASCII format")' ) &
           & file_id
   INFO(errMessage)

   Write(file_id,'(a)') ' ! RTTOV coefficient file '//Trim(coef % id_Common_name)
   Write(file_id,'(a)') ' ! automatic creation by subroutine Rttov_writecoef '
   Write(file_id,'(a)') ' ! ------------------------------------------------------'

   ! COEF structure (V7, V8, V9)

     ! IDENTIFICATION
     section = 'IDENTIFICATION'
     Select Case (coef % id_sensor)
     Case (1_jpim)
        sensor = 'ir'
     Case (2_jpim)
        sensor = 'mw'
     Case (3_jpim)
        sensor = 'hi'
     Case (4_jpim)
        sensor = 'po'
     End Select

     If( coef % id_comp_lvl == 7 ) then
       Write( errMessage, '( "output will be version 8 format compatible")' )
       INFO(errMessage)
       ! even if input coef file is version 7 format compatible
       ! the output will be version 8 format compatible
       ! The variable coef % fmv_model_ver contains the fast model
       ! version for the predictors and should never be modified
       id_comp_lvl = 8
     Else
       id_comp_lvl = coef % id_comp_lvl
     End If

     Write(file_id,'(a)') Trim(section)
     Write(file_id,'(a)') ' ! '
     Write(file_id,'(3i3,T20,a)')&
           & coef % id_platform, coef % id_sat, coef % id_inst,'! platform sat_id instrument'
     Write(file_id,'(1x,a)')  TRIM(coef % id_Common_name)
     Write(file_id,'(1x,a,T20,a)')  sensor,'! sensor type [ir,mw,hi]'
     Write(file_id,'(1x,i2,T20,a)') id_comp_lvl,'! RTTOV coefficient file version number'
     Write(file_id,'(1x,a)') TRIM(coef % id_creation)
     Write(file_id,'(1x,i4,1x,i2.2,1x,i2.2,t20,a)') coef % id_creation_date,'! creation date'

     ! No LINE-BY-LINE section


     If( coef%line_by_line(1) .ne. 'xxxx' ) Then
       section = 'LINE-BY-LINE'
       Write(file_id,'(a)') ' ! ------------------------------------------------------'
       Write(file_id,'(a)') Trim(section)
       Write(file_id,'(a)') ' ! '
       Do i = 1, 10
         If( coef%line_by_line(i) .eq. 'xxxx' ) Exit
         Write( file_id, '(a)' ) Trim(coef%line_by_line(i))
       EndDo
     EndIf

     ! FAST_MODEL_VARIABLES
     section = 'FAST_MODEL_VARIABLES'
     Write(file_id,'(a)') ' ! ------------------------------------------------------'
     Write(file_id,'(a)') Trim(section)
     Write(file_id,'(a)') ' ! '
     Write(file_id,'(a)') ' !'
     Write(file_id,'(1x,a,t20,a)')   coef % fmv_model_def,    '! fast model name'
     Write(file_id,'(1x,i4,t20,a)')  coef % fmv_model_ver,    '! fast model version compatibility level'
     Write(file_id,'(1x,i6,t20,a)')  coef % fmv_chn ,   '! Number of channels described in the coef file'
     Write(file_id,'(1x,i4,t20,a)')  coef % fmv_gas ,   '! Number of gases described in the coef file'
     If( coef % id_comp_lvl > 8 ) then
!        If( coef % IncTop ) Then
!          Write( file_id, * ) 1
!        Else
         Write( file_id, * ) 0
!        EndIf
     EndIf
     If( coef % id_comp_lvl == 9 ) then         !pb
       Do i = 1, coef % fmv_gas
         Write(file_id,'(1x,a,t20,a)')  Trim(gas_name( coef % fmv_gas_id( i ) ) ),'! gas identification'
         Write(file_id,'(1x,3i4,t20,a)')  coef % fmv_var(i), coef % fmv_coe(i), &
          &      coef % fmv_lvl(i), '! variables/predictors  levels (pressure/absorber)'
       End Do
     Else
       Do i = 1, coef % fmv_gas
         Write(file_id,'(1x,a,t20,a)')  Trim(gas_name( coef % fmv_gas_id( i ) ) ),'! gas identification'
         Write(file_id,'(1x,2i4,t20,a)')  coef % fmv_var(i), coef % fmv_lvl(i), '! variables/predictors  levels (pressure/absorber)'
       End Do
     End If

     ! GAZ_UNITS
     section = 'GAZ_UNITS'
     Write(file_id,'(a)') ' ! ------------------------------------------------------'
     Write(file_id,'(a)') Trim(section)
     Write(file_id,'(a)') ' ! Gaz concentrations can be expressed in '
     Write(file_id,'(a)') ' ! volume mixing ratio (ppmv)'
     Write(file_id,'(a)') ' ! specific concentration (kg/kg)'
     Write(file_id,'(a)') ' ! '
     Do i = 1, coef % fmv_gas
        Write(file_id,'(a)') ' !     '//gas_name( coef % fmv_gas_id( i ) )
        Write(file_id,'(1x,i4,t20,"! ",a)') &
                 & coef % gaz_units( i ), gas_unit_name( coef % gaz_units( i ) )
     End Do

     If( &
       ( coef % id_comp_lvl == 9 ) .and. &
       ( coef %nintmixed > 0 )     .and. &
       ( coef %nintwater > 0 )     .and. &
       ( coef %nintozone > 0 )     .and. &
       ( coef %nintwvcont > 0 )    .and. &
       ( coef %nintco2 > 0 )       .and. &
       ( coef %nintn2o > 0 )       .and. &
       ( coef %nintco > 0 )        .and. &
       ( coef %nintch4 > 0 )             &
     ) then         !pb
       section = 'GAS_SPECTRAL_INTERVAL'
       Write(file_id,'(a)') ' ! ------------------------------------------------------'
       Write(file_id,'(a)') Trim(section)
       Write(file_id,'(a)') ' ! '
       If ( coef %nintmixed > 0 ) Then
         Write(file_id,'(1x,a,t20,a)')  Trim(gas_name( gas_id_mixed ) ),'! gas identification'
         Write(file_id,'(1x,i4)')  coef %nintmixed
         Do i = 1, coef %nintmixed
           Write(file_id,'(1x,2f10.3)')  coef%mixedgasint(1,i),coef%mixedgasint(2,i)
         End Do
       End If
       If ( coef %nintwater > 0 ) Then
         Write(file_id,'(1x,a,t20,a)')  Trim(gas_name( gas_id_watervapour ) ),'! gas identification'
         Write(file_id,'(1x,i4)')  coef %nintwater
         Do i = 1, coef %nintwater
           Write(file_id,'(1x,2f10.3)')  coef%watervapourint(1,i),coef%watervapourint(2,i)
         End Do
       End If
       If ( coef %nintozone > 0 ) Then
         Write(file_id,'(1x,a,t20,a)')  Trim(gas_name( gas_id_ozone ) ),'! gas identification'
         Write(file_id,'(1x,i4)')  coef %nintozone
         Do i = 1, coef %nintozone
           Write(file_id,'(1x,2f10.3)')  coef%ozoneint(1,i),coef%ozoneint(2,i)
         End Do
       End If
       If ( coef %nintwvcont > 0 ) Then
         Write(file_id,'(1x,a,t20,a)')  Trim(gas_name( gas_id_wvcont ) ),'! gas identification'
         Write(file_id,'(1x,i4)')  coef %nintwvcont
         Do i = 1, coef %nintwvcont
           Write(file_id,'(1x,2f10.3)')  coef%wvcontint(1,i),coef%wvcontint(2,i)
         End Do
       End If
       If ( coef %nintco2 > 0 ) Then
         Write(file_id,'(1x,a,t20,a)')  Trim(gas_name( gas_id_co2 ) ),'! gas identification'
         Write(file_id,'(1x,i4)')  coef %nintco2
         Do i = 1, coef %nintco2
           Write(file_id,'(1x,2f10.3)')  coef%co2int(1,i),coef%co2int(2,i)
         End Do
       End If
       If ( coef %nintn2o > 0 ) Then
         Write(file_id,'(1x,a,t20,a)')  Trim(gas_name( gas_id_n2o ) ),'! gas identification'
         Write(file_id,'(1x,i4)')  coef %nintn2o
         Do i = 1, coef %nintn2o
           Write(file_id,'(1x,2f10.3)')  coef%n2oint(1,i),coef%n2oint(2,i)
         End Do
       End If
       If ( coef %nintco > 0 ) Then
         Write(file_id,'(1x,a,t20,a)')  Trim(gas_name( gas_id_co ) ),'! gas identification'
         Write(file_id,'(1x,i4)')  coef %nintco
         Do i = 1, coef %nintco
           Write(file_id,'(1x,2f10.3)')  coef%coint(1,i),coef%coint(2,i)
         End Do
       End If
       If ( coef %nintch4 > 0 ) Then
         Write(file_id,'(1x,a,t20,a)')  Trim(gas_name( gas_id_ch4 ) ),'! gas identification'
         Write(file_id,'(1x,i4)')  coef %nintch4
         Do i = 1, coef %nintch4
           Write(file_id,'(1x,2f10.3)')  coef%ch4int(1,i),coef%ch4int(2,i)
         End Do
       End If
     End If


     If( coef % id_comp_lvl == 9 .and. Associated( coef % tt_chn ) ) then         !pb
       section = 'TRANSMITTANCE_TRESHOLD'
       Write(file_id,'(a)') ' ! ------------------------------------------------------'
       Write(file_id,'(a)') Trim(section)
       Write(file_id,'(a)') ' ! '
       Write(file_id,'(a)') ' ! chan number'
       Write(file_id,'(a)') ' ! validity of channel '
       Write(file_id,'(a)') ' ! central wave number'
       Write(file_id,'(a)') ' ! transmittance treshold'
       Write(file_id,'(a)') ' ! transmittance value   '
       Do i = 1, coef % fmv_chn
         Write(file_id,'(2i5,3e19.10)',iostat=io_status)&
               & coef % tt_chn(i)    ,  &
               & coef % tt_val_chn(i),  &
               & coef % tt_cwn(i)    ,  &
               & coef % tt_a0(i)     ,  &
               & coef % tt_a1(i)
         If(io_status /= 0) Then
           errMessage="write TRANSMITTANCE_TRESHOLD section"
           Call Rttov_ErrorReport (Int(io_status,jpim), errMessage, NameOfRoutine)
           errorstatus = errorstatus_fatal
           Return
         Endif
       End Do
     End If                                     !pb

     If( coef % id_comp_lvl == 9 .and. Associated( coef % ss_chn ) ) then         !pb
       section = 'SOLAR_SPECTRUM'
       Write(file_id,'(a)') ' ! ------------------------------------------------------'
       Write(file_id,'(a)') Trim(section)
       Write(file_id,'(a)') ' ! '
       Write(file_id,'(a)') ' ! chan number'
       Write(file_id,'(a)') ' ! validity of channel'
       Write(file_id,'(a)') ' ! central wave number'
       Write(file_id,'(a)') ' ! solar spectrum'
       Do i = 1, coef % fmv_chn
         Write(file_id,'(2i5,2e19.10)',iostat=io_status)&
                 & coef % ss_chn(i)    ,  &
                 & coef % ss_val_chn(i),  &
                 & coef % ss_cwn(i)    ,  &
                 & coef % ss_solar_spectrum(i)
         If(io_status /= 0) Then
           errMessage="write SOLAR_SPECTRUM section"
           Call Rttov_ErrorReport (Int(io_status,jpim), errMessage, NameOfRoutine)
           errorstatus = errorstatus_fatal
           Return
         Endif
       End Do
     End If                                     !pb

     If( coef % id_comp_lvl == 9 .and. associated( coef % woc_chn ) ) then         !pb
       section = 'WATER_OPTICAL_CONSTANT'
       Write(file_id,'(a)') ' ! ------------------------------------------------------'
       Write(file_id,'(a)') Trim(section)
       Write(file_id,'(a)') ' ! '
       Write(file_id,'(a)') ' ! chan number'
       Write(file_id,'(a)') ' ! central wave number'
       Write(file_id,'(a)') ' ! ocean water optical constants(real and imaginary part)'
       Write(file_id,'(a)') ' ! fresh water optical constants(real and imaginary part)'
       Do i = 1, coef % fmv_chn
         Write(file_id,'(i5,e19.10,2(" (",e17.10,",",e17.10,")"))',iostat=io_status)&
               & coef % woc_chn(i)    ,   &
               & coef % woc_cwn(i)    ,   &
               & coef % woc_waopc_ow(i)  ,&
               & coef % woc_waopc_fw(i)
         If(io_status /= 0) Then
           errMessage="write WATER_OPTICAL_CONSTANT section"
           Call Rttov_ErrorReport (Int(io_status,jpim), errMessage, NameOfRoutine)
           errorstatus = errorstatus_fatal
           Return
         Endif
       End Do
     End If                                     !pb

     If( coef % id_comp_lvl == 9 .and. Associated( coef%ws_npoint ) ) then         !pb
       section = 'WAVE_SPECTRUM'
       Write(file_id,'(a)') ' ! ------------------------------------------------------'
       Write(file_id,'(a)') Trim(section)
       Write(file_id,'(a)') ' ! '
       Write(file_id,'(a)') ' ! Number of points'
       Write(file_id,'(a)') ' ! Point number'
       Write(file_id,'(a)') ' ! WAve spectrum'
       Write(file_id,*,iostat=io_status) coef%ws_nomega
       Do i = 1, coef % ws_nomega
         Write(file_id,'(f10.3,f12.5)',iostat=io_status)coef%ws_npoint(i),coef%ws_k_omega(i)
         If(io_status /= 0) Then
           errMessage="write WAVE_SPECTRUM section"
           Call Rttov_ErrorReport (Int(io_status,jpim), errMessage, NameOfRoutine)
           errorstatus = errorstatus_fatal
           Return
         Endif
       End Do
     End If                                     !pb





     section = 'FILTER_FUNCTIONS'
     Write(file_id,'(a)') ' ! ------------------------------------------------------'
     Write(file_id,'(a)') Trim(section)
     Write(file_id,'(a)') ' ! '
     Write(file_id,'(a)') ' ! Channel Number (from instrument original description)'
     Write(file_id,'(a)') ' ! Channel status '
     Write(file_id,'(a)') ' ! Central Wavenumber'
     Write(file_id,'(a)') ' ! Band Correction coefficients(Offset,Slope)'
     Write(file_id,'(a)') ' ! Gamma correction factor'

     Do i = 1, coef % fmv_chn
        Write(file_id,'(1x,i5,1x,i4,4(1x,e18.10))') &
              & coef % ff_ori_chn(i), coef % ff_val_chn(i), coef % ff_cwn(i),&
              & coef % ff_bco(i), coef % ff_bcs(i), coef % ff_gam(i)
     End Do


     section = 'FUNDAMENTAL_CONSTANTS'
     Write(file_id,'(a)') ' ! ------------------------------------------------------'
     Write(file_id,'(a)') Trim(section)
     Write(file_id,'(a)') ' ! '
     Write(file_id,'(a)') ' ! units of constants for spectral radiance'
     Write(file_id,'(a)') ' ! first radiation constant(mW/(m2.sr.cm-4))'
     Write(file_id,'(a)') ' ! second radiation constant (cm.K)'
     Write(file_id,'(1x,f14.1,t30,a)') speedl,'! speed of light (cm/s)'
     Write(file_id,'(1x,1p,e15.8,0p,f10.6,t30,a)') coef % fc_planck_c1, coef % fc_planck_c2,'! Planck constants'
     Write(file_id,'(1x,f8.1,t30,a)') coef % fc_sat_height,'! nominal satellite height (km)'

     If( coef % fastem_ver >= 1 ) Then
        section = 'FASTEM'
        Write(file_id,'(a)') ' ! ------------------------------------------------------'
        Write(file_id,'(a)') Trim(section)
        Write(file_id,'(a)') ' ! '
        Write(file_id,'(a)') ' ! S. English fast generic millimetre wave ocean emissivity model'
        Write(file_id,'(a)') ' ! Polarisation of each channel', &
              & ' !       MPOL=0:  0.5_JPRB*(V+H)', &
              & ' !       MPOL=1: polarisation angle=90-incidence angle', &
              & ' !       MPOL=2: polarisation angle=incidence angle', &
              & ' !       MPOL=3: vertical polarisation', &
              & ' !       MPOL=4: horizontal polarisation'
        Write(file_id,'(1x,i2,a)') coef % fastem_ver,'   ! version number'
        Write(file_id,'(1x,i3,a)') size(fastem3_coef),'  ! number of coefficients'
        Write(file_id,'(5e14.6)') fastem3_coef
        Write(file_id,'(20i3)') (coef % fastem_polar(i), i= 1, coef % fmv_chn)
     Endif

     If( coef % ssirem_ver >= 1 ) Then
        section = 'SSIREM'
        Write(file_id,'(a)') ' ! ------------------------------------------------------'
        Write(file_id,'(a)') Trim(section)
        Write(file_id,'(a)') ' ! '
        Write(file_id,'(a)') ' ! Channel Number (from instrument original description)'
        Write(file_id,'(a)') ' ! 5 coefficients for emissivity model ssirem'
        Write(file_id,'(1x,i2,a)') coef % ssirem_ver,'   ! version number'

        Do i = 1, coef % fmv_chn
           Write(file_id,'(1x,i5,3f12.7,2f4.1)') &
                 & coef % ssirem_chn(i) , coef % ssirem_a0(i),&
                 & coef % ssirem_a1(i)  , coef % ssirem_a2(i),&
                 & coef % ssirem_xzn1(i), coef % ssirem_xzn2(i)
        End Do
     Endif

     section = 'REFERENCE_PROFILE'
     Write(file_id,'(a)') ' ! ------------------------------------------------------'
     Write(file_id,'(a)') Trim(section)
     Write(file_id,'(a)') ' ! '
     Write(file_id,'(a)') ' ! Ref.pressure (hPa)'
     Write(file_id,'(a)') ' ! Ref.Temp (K) Ref.Volume Mixing Ratio [ppmv] for each gas'
     Write(file_id,'(a)') ' ! Note for MxG that mixing ratio is "missing"'

     Do i = 1, coef % fmv_gas
        Write(file_id,'(a)') ' !     '//gas_name( coef % fmv_gas_id( i ) )
        Do l = 1, coef % fmv_lvl(i)
           Write(file_id,'(1x,f9.4,2x,f7.3,1x,e13.6)')&
                 & coef % ref_prfl_p(l), coef % ref_prfl_t(l,i), coef % ref_prfl_mr(l,i)
  !!$                & coef % ref_prfl_p(l), coef % ref_prfl_t(l,i), ref_mr(l,i)
        End Do
     End Do

     section = 'PROFILE_LIMITS'
     Write(file_id,'(a)') ' ! ------------------------------------------------------'
     Write(file_id,'(a)') Trim(section)
     Write(file_id,'(a)') ' ! '
     Write(file_id,'(a)') ' ! Ref.pressure (hPa)'
     Write(file_id,'(a)') ' ! Temp Max (K) Temp Min (K)'
     Write(file_id,'(a)') ' ! Volume Mixing Ratio for  Max and Min [ppmv] for each gas'
     Write(file_id,'(a)') ' !      Temperature'
     Do l = 1, coef % fmv_lvl(1)
        Write(file_id,'(1x,f9.4,2(1x,f7.2))',iostat=io_status)&
              & coef % lim_prfl_p(l), coef % lim_prfl_tmax(l), coef % lim_prfl_tmin(l)
     End Do

     Do i = 1, coef % fmv_gas
        Write(file_id,'(a)') ' !     '//gas_name( coef % fmv_gas_id( i ) )
        Do l = 1, coef % fmv_lvl(i)
           Write(file_id,'(1x,f9.4,2x,e12.4,e12.4)',iostat=io_status)&
                 & coef % lim_prfl_p(l), coef % lim_prfl_gmax(l,i), coef % lim_prfl_gmin(l,i)
        End Do
     End Do


     section = 'FAST_COEFFICIENTS'
     Write(file_id,'(a)') ' ! ------------------------------------------------------'
     Write(file_id,'(a)') Trim(section)
     Write(file_id,'(a)') ' ! '
     Write(file_id,'(a)') ' ! transmission coefficients'
     Write(file_id,'(a)') ' ! Order of the gases:'
     Do i = 1, coef % fmv_gas
        Write(file_id,'(a)') ' !     '//gas_name( coef % fmv_gas_id ( i ) )
     End Do


     Do l = 1, coef % fmv_gas
        Write(file_id,'(a)') gas_name( coef % fmv_gas_id( l ) )

        Select Case( coef % fmv_gas_id(l) )

        Case(gas_id_mixed)
           Write(file_id,'(5(1x,e15.8))',iostat=io_status)  &
                 & (((coef % thermal % mixedgas(i,j,k)   &
                  & ,i = 1, coef % fmv_lvl(l) ) &
                  & ,j = 1, coef % fmv_chn    ) &
                  & ,k = 1, coef % ncmixed )
        Case(gas_id_watervapour)
           Write(file_id,'(5(1x,e15.8))',iostat=io_status)  &
                 & (((coef % thermal % watervapour(i,j,k)   &
                  & ,i = 1, coef % fmv_lvl(l) ) &
                  & ,j = 1, coef % fmv_chn    ) &
                  & ,k = 1, coef % ncwater )
        Case(gas_id_ozone)
           Write(file_id,'(5(1x,e15.8))',iostat=io_status)  &
                 & (((coef % thermal % ozone(i,j,k)   &
                  & ,i = 1, coef % fmv_lvl(l) ) &
                  & ,j = 1, coef % fmv_chn    ) &
                  & ,k = 1, coef % ncozone )
        Case(gas_id_wvcont)
           Write(file_id,'(5(1x,e15.8))',iostat=io_status)  &
                 & (((coef % thermal % wvcont(i,j,k)   &
                  & ,i = 1, coef % fmv_lvl(l) ) &
                  & ,j = 1, coef % fmv_chn    ) &
                  & ,k = 1, coef % ncwvcont )
        Case(gas_id_co2)
           Write(file_id,'(5(1x,e15.8))',iostat=io_status)  &
                 & (((coef % thermal % co2(i,j,k)   &
                  & ,i = 1, coef % fmv_lvl(l) ) &
                  & ,j = 1, coef % fmv_chn    ) &
                  & ,k = 1, coef % ncco2 )
        Case(gas_id_n2o)
           Write(file_id,'(5(1x,e15.8))',iostat=io_status)  &
                 & (((coef % thermal % n2o(i,j,k)   &
                  & ,i = 1, coef % fmv_lvl(l) ) &
                  & ,j = 1, coef % fmv_chn    ) &
                  & ,k = 1, coef % ncn2o )
        Case(gas_id_co)
           Write(file_id,'(5(1x,e15.8))',iostat=io_status)  &
                 & (((coef % thermal % co(i,j,k)   &
                  & ,i = 1, coef % fmv_lvl(l) ) &
                  & ,j = 1, coef % fmv_chn    ) &
                  & ,k = 1, coef % ncco )
        Case(gas_id_ch4)
           Write(file_id,'(5(1x,e15.8))',iostat=io_status)  &
                 & (((coef % thermal % ch4(i,j,k)   &
                  & ,i = 1, coef % fmv_lvl(l) ) &
                  & ,j = 1, coef % fmv_chn    ) &
                  & ,k = 1, coef % ncch4 )
        End Select
     End Do

     section = 'END'
     Write(file_id,'(a)') ' ! ------------------------------------------------------'
     Write(file_id,'(a)') Trim(section)


  If( io_status /= 0 ) Then
     Write( errMessage, '( "write IO error")' )
     Call Rttov_ErrorReport (errorstatus, errMessage, NameOfRoutine)
     Return
  End If

  Write( errMessage, '( "end of write coefficient")' ) 
  INFO(errMessage)
  CATCH
End Subroutine rttov789_write_ascii_coef
