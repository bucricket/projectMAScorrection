!
SUBROUTINE rttov_read_binary_coef( &
            & ERR,           &
            & coef,          &
            & file_lu,       &
            & channels)
! Description:
!
! Read an binary coefficient file and fills coeff structure
!   arrays according to the optional list of channels.
!
! The user can provide an optional list of channels in "channels" argument
!  array to reduce the output coefficient structure to this list. This
! can be important for reducing the memory allocation required when running
! with advanced IR sounders (e.g. AIRS or IASI). If the user
!  wants all channels the "channels" argument shall not be present.
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
! Method:
!
! Current Code Owner: SAF NWP
!
! History:
! Version         Date         Author    Comment
! -------         ----         ------    -------
!  1.0            01/12/2002             New F90 code with structures (P Brunel A Smith)
!  1.1            02/01/2003             A few comments added (R Saunders)
!  1.2            24/01/2003             add tests on all read statements (P Brunel)
!                                        one record per channel for coefficients in binary format
!                                        New header to allow checking R4<->R8
!  1.3            06/05/2003             Remove "optional" attribute of argument file_lu (P Brunel)
!  1.4            02/06/2004             New format for FMV section with RTTOV8  (P. Brunel)
!  1.5            15/06/2004             Corrected array dimension for coef % fmv_gas_pos (R Saunders)
!  1.6            14/05/2007             Updated for RTTOV-9 ( P Brunel)
!  1.7            12/12/2007             Option of additional layer at top flag read (R Saunders)
!  1.8            27/06/2008             Introduced the case when no channels are present for the
!                                        phase function in the solar range (M. Matricardi)
!  1.9            27/02/2009             Profile levels to include ToA. Allocate coef arrays
!                                        according to number of layers, not levels (P. Rayer)
!  1.10           06/03/2009             Separation of flags for IncZeeman and IncTop.
!                                        Conditionals depending on coef % id_comp_lvl == 9
!                                        extended to >= 9 (P. Rayer)
!  1.11           08/06/2009             Made interim fix to allocate cloud/aerosol arrays with right shape (R Saunders)
!                                        Ideally the channel order in all IR files needs to be in frequency not wavelength
!  1.12           02/12/2009             Introduced principal component capability. Marco Matricardi. ECMWF
!  r1241 (pre-11) 12/2012      DAR       Added support for NLTE coefs
!  1.14           10/01/2013             Add PMC shifts (P Rayer)
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
#include "throw.h"
! Imported Type Definitions:
  USE rttov_types, ONLY :  &
       & rttov_coef
  USE parkind1, ONLY : jpim
!INTF_OFF
  USE rttov_const, ONLY :  &
        rttov_magic_string,     &
        rttov_magic_number,     &
        version_compatible_min, &
        version_compatible_max, &
        gas_id_mixed,           &
        gas_id_watervapour,     &
        gas_id_ozone,           &
        gas_id_wvcont,          &
        gas_id_co2,             &
        gas_id_n2o,             &
        gas_id_co,              &
        gas_id_ch4,             &
        ngases_max

  USE parkind1, ONLY : jprb, jplm
!INTF_ON
  IMPLICIT NONE
! subroutine arguments
! scalar arguments with intent(in):
  INTEGER(KIND=jpim)       , INTENT(IN)                :: file_lu       ! file logical unit number
  INTEGER(KIND=jpim)       , OPTIONAL     , INTENT(IN) :: channels   (:)! list of channels to extract
! scalar arguments with intent(inout):
  TYPE(rttov_coef         ), INTENT(INOUT)             :: coef          ! coefficients
! scalar arguments with intent(out):
  INTEGER(KIND=jpim)       , INTENT(OUT)               :: ERR           ! return code
!INTF_END
#include "rttov_errorreport.interface"
! Local Scalars:
  LOGICAL(KIND=jplm) :: all_channels
  INTEGER(KIND=jpim) :: n
  INTEGER(KIND=jpim) :: chn
  INTEGER(KIND=jpim) :: i, IncZeeman
  INTEGER(KIND=jpim) :: nlte_count, nlte_start, nlte_file_nchan
  INTEGER(KIND=jpim), ALLOCATABLE :: nlte_chans(:)
  REAL   (KIND=jprb) :: fc_speedl
! pointers for generic inputs
  COMPLEX(KIND=jprb), POINTER :: values_c_0(:)
  COMPLEX(KIND=jprb), POINTER :: values_c_1(:)
  REAL   (KIND=jprb), POINTER :: values0(:)
  REAL   (KIND=jprb), POINTER :: values1(:)
  REAL   (KIND=jprb), POINTER :: values2(:)
  REAL   (KIND=jprb), POINTER :: values3(:)
  REAL   (KIND=jprb), POINTER :: values4(:)
  REAL   (KIND=jprb), POINTER :: nlte_values(:,:,:,:)
  INTEGER(KIND=jpim), POINTER :: ivalues0(:)
  INTEGER(KIND=jpim), POINTER :: ivalues1(:)
  CHARACTER(LEN = 16) :: bin_check_string
  REAL(KIND=jprb)     :: bin_check_number
  REAL(KIND=jprb)     :: bin_check_value
  CHARACTER(LEN = 80) :: errMessage
  LOGICAL(KIND=jplm)  :: section_present
!- End of header --------------------------------------------------------
  TRY
! 0 Initialise variables
!---------------------------------------------
! test presence of channels argument

  IF (Present(channels)) THEN
    all_channels = .FALSE.
  ELSE
    all_channels = .TRUE.
  ENDIF

! 3 Read binary file
!-------------------
! Binary file
  READ (file_lu, iostat=ERR)bin_check_string, bin_check_number

  THROWM( ERR .NE. 0, 'io status while reading header')

! Verification of header string
  IF (bin_check_string /= rttov_magic_string) err = errorstatus_fatal

  THROWM( ERR .NE. 0,'Wrong header string in file')

! Verification of single/double precision using a 5 digit number
! with exponent 12, which is always Ok for single precision
  bin_check_value = 1._JPRB - abs(bin_check_number - rttov_magic_number)
  IF (bin_check_value > 1.01_JPRB .OR. bin_check_value < 0.99_JPRB) err = errorstatus_fatal

  THROWM( ERR .NE. 0,'File created with a different real precision (R4<->R8)')

! COEF structure (V10)

  errMessage = 'io status while reading IDENTIFICATION'
  READ (file_lu, iostat=ERR)coef%id_platform, coef%id_sat, coef%id_inst, coef%id_sensor

  THROWM(err.NE.0,errMessage)

  READ (file_lu, iostat=ERR)coef%id_comp_lvl, coef%id_creation_date, coef%id_creation, coef%id_Common_name

  THROWM(err.NE.0,errMessage)

  errMessage = 'io status while reading FAST_MODEL_VARIABLES'

  READ (file_lu, iostat=ERR)coef%fmv_model_def, coef%fmv_model_ver, coef%fmv_ori_nchn, coef%fmv_gas
  
  THROWM(err.NE.0,errMessage)

  coef%id_comp_pc = 0_jpim

  READ (file_lu, iostat=ERR)coef%id_comp_pc

  THROWM(err.NE.0,errMessage)
    
  coef%IncZeeman = .FALSE.

  READ (file_lu, iostat=ERR)IncZeeman

  THROWM(err.NE.0,errMessage)

  IF (IncZeeman == 0) coef%IncZeeman = .FALSE.
  IF (IncZeeman == 1) coef%IncZeeman = .TRUE.

! Error if the compatibility version of the coefficient file
! is not in the range defined by the constant module

  IF (coef%id_comp_lvl < version_compatible_min .OR. coef%id_comp_lvl > version_compatible_max) THEN
    err = errorstatus_fatal

    THROWM(err.NE.0,"Version of coefficient file is incompatible with RTTOV library")
    
  ENDIF

! Take care of the user list of channels
! coef%fmv_ori_nchn store the number of channels in the file
! coef % fmv_chn is the number of channels that the user requests

  IF (all_channels) THEN
    coef%fmv_chn = coef%fmv_ori_nchn
  ELSE
    IF (MAXVAL(channels) > coef%fmv_ori_nchn) THEN
      err = errorstatus_fatal
      THROWM(err.ne.0,"Channel index out of range for coefficient file")
    ENDIF
    coef%fmv_chn = SIZE(channels)
  ENDIF
  
  ALLOCATE (coef%fmv_gas_id(coef%fmv_gas), STAT = ERR)

  THROWM( ERR .NE. 0, "allocation of fmv_gas_id")
  
  ALLOCATE (coef%fmv_gas_pos(ngases_max), STAT = ERR)
  
  THROWM( ERR .NE. 0, "allocation of fmv_gas_pos")
  
  ALLOCATE (coef%fmv_var(coef%fmv_gas), STAT = ERR)
  
  THROWM( ERR .NE. 0, "allocation of fmv_var")
    
  ALLOCATE (coef%fmv_coe(coef%fmv_gas), STAT = ERR)

  THROWM( ERR .NE. 0, "allocation of fmv_coe")
!pb
  ALLOCATE (coef%fmv_lvl(coef%fmv_gas), STAT = ERR)

  THROWM( ERR .NE. 0, "allocation of fmv_lvl")

  ALLOCATE (coef%ff_ori_chn(coef%fmv_chn), STAT = ERR)
  
  THROWM( ERR .NE. 0, "allocation of ff_ori_chn")
  
  ALLOCATE (coef%ff_val_chn(coef%fmv_chn), STAT = ERR)
  
  THROWM( ERR .NE. 0, "allocation of ff_val_chn")
  
  ALLOCATE (coef%ff_cwn(coef%fmv_chn), STAT = ERR)
  
  THROWM( ERR .NE. 0, "allocation of ff_cwn")
  
  ALLOCATE (coef%ff_bco(coef%fmv_chn), STAT = ERR)
  
  THROWM( ERR .NE. 0, "allocation of ff_bco")
  
  ALLOCATE (coef%ff_bcs(coef%fmv_chn), STAT = ERR)
  
  THROWM( ERR .NE. 0, "allocation of ff_bcs")
  
  ALLOCATE (coef%ff_gam(coef%fmv_chn), STAT = ERR)
  
  THROWM( ERR .NE. 0, "allocation of ff_gam")
  
  ALLOCATE (coef%gaz_units(coef%fmv_gas), STAT = ERR)
  
  THROWM( ERR .NE. 0, "allocation of gaz_units")
  
  coef%fmv_gas_id(:)  = 0_jpim
  coef%fmv_gas_pos(:) = 0_jpim
  coef%fmv_var(:)     = 0_jpim
  coef%fmv_lvl(:)     = 0_jpim
  coef%fmv_coe(:)     = 0_jpim
  READ (file_lu, iostat=ERR)coef%fmv_gas_id, coef%fmv_gas_pos, coef%fmv_var, coef%fmv_lvl
  
  THROWM(err.NE.0,errMessage)
    
  READ (file_lu, iostat=ERR)coef%fmv_coe!pb
  
  THROWM(err.NE.0,errMessage)
  
  DO n = 1, coef%fmv_gas
    
    SELECT CASE (coef%fmv_gas_id(n))
    CASE (gas_id_mixed)
        coef%nmixed  = coef%fmv_var(n)
        coef%nlevels = coef%fmv_lvl(n)
        coef%ncmixed = coef%fmv_coe(n)!pb
      CASE (gas_id_watervapour)
        coef%nwater  = coef%fmv_var(n)
        coef%ncwater = coef%fmv_coe(n)!pb
      CASE (gas_id_ozone)
        coef%nozone  = coef%fmv_var(n)
        coef%ncozone = coef%fmv_coe(n)!pb
      CASE (gas_id_wvcont)
        coef%nwvcont  = coef%fmv_var(n)
        coef%ncwvcont = coef%fmv_coe(n)!pb
      CASE (gas_id_co2)
        coef%nco2  = coef%fmv_var(n)
        coef%ncco2 = coef%fmv_coe(n)!pb
      CASE (gas_id_n2o)
        coef%nn2o  = coef%fmv_var(n)
        coef%ncn2o = coef%fmv_coe(n)!pb
      CASE (gas_id_co)
        coef%nco  = coef%fmv_var(n)
        coef%ncco = coef%fmv_coe(n)!pb
      CASE (gas_id_ch4)
        coef%nch4  = coef%fmv_var(n)
        coef%ncch4 = coef%fmv_coe(n)!pb
      END SELECT

    ENDDO

    coef%nlayers = coef%nlevels - 1
    errMessage   = 'io status while reading GAZ_UNITS'
    READ (file_lu, iostat=ERR)coef%gaz_units

    THROWM(err.ne.0,errMessage)


! GAS_SPECTRAL_INTERVAL                  !pb
    READ (file_lu, iostat=ERR)coef%nintmixed, coef%nintwater, coef%nintozone, coef%nintwvcont, coef%nintco2,      &
      & coef%nintn2o, coef%nintco, coef%nintch4!pb

    THROWM(err.ne.0,errMessage)


    IF (coef%nintmixed > 0) THEN!pb
      ALLOCATE (coef%mixedgasint(2, coef%nintmixed), STAT = ERR)

      THROWM( ERR .NE. 0, "allocation of mixedgasint arra")

      READ (file_lu, iostat=ERR)coef%mixedgasint!pb

      THROWM(err.ne.0,errMessage)

    ENDIF


    IF (coef%nintwater > 0) THEN!pb
      ALLOCATE (coef%watervapourint(2, coef%nintwater), STAT = ERR)

      THROWM( ERR .NE. 0, "allocation of watervapourint array")

      READ (file_lu, iostat=ERR)coef%watervapourint!pb

      THROWM(err.ne.0,errMessage)

    ENDIF


    IF (coef%nintozone > 0) THEN!pb
      ALLOCATE (coef%ozoneint(2, coef%nintozone), STAT = ERR)

      THROWM( ERR .NE. 0, "allocation of ozoneint array")

      READ (file_lu, iostat=ERR)coef%ozoneint!pb

      THROWM(err.ne.0,errMessage)

    ENDIF


    IF (coef%nintwvcont > 0) THEN!pb
      ALLOCATE (coef%wvcontint(2, coef%nintwvcont), STAT = ERR)

      THROWM( ERR .NE. 0, "allocation of wvcontint array")

      READ (file_lu, iostat=ERR)coef%wvcontint!pb

      THROWM(err.ne.0,errMessage)

    ENDIF


    IF (coef%nintco2 > 0) THEN!pb
      ALLOCATE (coef%co2int(2, coef%nintco2), STAT = ERR)

      THROWM( ERR .NE. 0, "allocation of co2int array")

      READ (file_lu, iostat=ERR)coef%co2int!pb

      THROWM(err.ne.0,errMessage)

    ENDIF


    IF (coef%nintn2o > 0) THEN!pb
      ALLOCATE (coef%n2oint(2, coef%nintn2o), STAT = ERR)

      THROWM( ERR .NE. 0, "allocation of n2oint array")

      READ (file_lu, iostat=ERR)coef%n2oint!pb

      THROWM(err.ne.0,errMessage)

    ENDIF


    IF (coef%nintco > 0) THEN!pb
      ALLOCATE (coef%coint(2, coef%nintco), STAT = ERR)

      THROWM( ERR .NE. 0, "allocation of coint array")

      READ (file_lu, iostat=ERR)coef%coint!pb

      THROWM(err.ne.0,errMessage)

    ENDIF


    IF (coef%nintch4 > 0) THEN!pb
      ALLOCATE (coef%ch4int(2, coef%nintch4), STAT = ERR)

      THROWM( ERR .NE. 0, "allocation of ch4int array")

      READ (file_lu, iostat=ERR)coef%ch4int!pb

      THROWM(err.ne.0,errMessage)

    ENDIF

      
    READ (file_lu)section_present
    errMessage = 'io status while reading TRANSMITTANCE_TRESHOLD'

    IF (section_present) THEN
!TRANSMITTANCE_TRESHOLD                  !pb
      ALLOCATE (coef%tt_chn(coef%fmv_chn), STAT = ERR)

      THROWM( ERR .NE. 0, "allocation of tt_chn")

      ALLOCATE (coef%tt_val_chn(coef%fmv_chn), STAT = ERR)

      THROWM( ERR .NE. 0, "allocation of tt_val_chn")

      ALLOCATE (coef%tt_cwn(coef%fmv_chn), STAT = ERR)

      THROWM( ERR .NE. 0, "allocation of tt_cwn")

      ALLOCATE (coef%tt_a0(coef%fmv_chn), STAT = ERR)

      THROWM( ERR .NE. 0, "allocation of tt_a0")

      ALLOCATE (coef%tt_a1(coef%fmv_chn), STAT = ERR)

      THROWM( ERR .NE. 0, "allocation of tt_a1")


      IF (all_channels) THEN
        READ (file_lu, iostat=ERR)coef%tt_chn, coef%tt_val_chn, coef%tt_cwn, coef%tt_a0, coef%tt_a1!pb

        THROWM(err.ne.0,errMessage)

      ELSE
        ALLOCATE (ivalues0(coef%fmv_ori_nchn), STAT = ERR)

        THROWM( ERR .NE. 0, "allocation of ivalues0")

        ALLOCATE (ivalues1(coef%fmv_ori_nchn), STAT = ERR)

        THROWM( ERR .NE. 0, "allocation of ivalues1")

        ALLOCATE (values0(coef%fmv_ori_nchn), STAT = ERR)

        THROWM( ERR .NE. 0, "allocation of values0")

        ALLOCATE (values1(coef%fmv_ori_nchn), STAT = ERR)

        THROWM( ERR .NE. 0, "allocation of values1")

        ALLOCATE (values2(coef%fmv_ori_nchn), STAT = ERR)

        THROWM( ERR .NE. 0, "allocation of values2 ")

        READ (file_lu, iostat=ERR)ivalues0, ivalues1, values0, values1, values2!pb

        THROWM(err.ne.0,errMessage)

        coef%tt_chn(:)     = ivalues0(channels(:))
        coef%tt_val_chn(:) = ivalues1(channels(:))
        coef%tt_cwn(:)     = values0(channels(:))
        coef%tt_a0(:)      = values1(channels(:))
        coef%tt_a1(:)      = values2(channels(:))
        DEALLOCATE (ivalues0, STAT = ERR)

        THROWM( ERR .NE. 0, "deallocation of ivalues0")

        DEALLOCATE (ivalues1, STAT = ERR)

        THROWM( ERR .NE. 0, "deallocation of ivalues1")

        DEALLOCATE (values0, STAT = ERR)

        THROWM( ERR .NE. 0, "deallocation of values0")

        DEALLOCATE (values1, STAT = ERR)

        THROWM( ERR .NE. 0, "deallocation of values1")

        DEALLOCATE (values2, STAT = ERR)

        THROWM( ERR .NE. 0, "deallocation of values2")

      ENDIF

    ENDIF

    
    READ (file_lu)section_present
    errMessage = 'io status while reading SOLAR_SPECTRUM'

    IF (section_present) THEN
!SOLAR_SPECTRUM                          !pb
      ALLOCATE (coef%ss_chn(coef%fmv_chn), STAT = ERR)

      THROWM( ERR .NE. 0, "allocation of ss_chn")

      ALLOCATE (coef%ss_val_chn(coef%fmv_chn), STAT = ERR)

      THROWM( ERR .NE. 0, "allocation of ss_val_chn")

      ALLOCATE (coef%ss_cwn(coef%fmv_chn), STAT = ERR)

      THROWM( ERR .NE. 0, "allocation of ss_cwn")

      ALLOCATE (coef%ss_solar_spectrum(coef%fmv_chn), STAT = ERR)

      THROWM( ERR .NE. 0, "allocation of ss_solar_spectrum")


      IF (all_channels) THEN
        READ (file_lu, iostat=ERR)coef%ss_chn, coef%ss_val_chn, coef%ss_cwn, coef%ss_solar_spectrum!pb

        THROWM(err.ne.0,errMessage)

      ELSE
        ALLOCATE (ivalues0(coef%fmv_ori_nchn), STAT = ERR)

        THROWM( ERR .NE. 0, "allocation of ivalues0")

        ALLOCATE (ivalues1(coef%fmv_ori_nchn), STAT = ERR)

        THROWM( ERR .NE. 0, "allocation of ivalues1")

        ALLOCATE (values0(coef%fmv_ori_nchn), STAT = ERR)

        THROWM( ERR .NE. 0, "allocation of values0")

        ALLOCATE (values1(coef%fmv_ori_nchn), STAT = ERR)

        THROWM( ERR .NE. 0, "allocation of values1")

        READ (file_lu, iostat=ERR)ivalues0, ivalues1, values0, values1

        THROWM(err.ne.0,errMessage)

        coef%ss_chn(:)            = ivalues0(channels(:))
        coef%ss_val_chn(:)        = ivalues1(channels(:))
        coef%ss_cwn(:)            = values0(channels(:))
        coef%ss_solar_spectrum(:) = values1(channels(:))
        DEALLOCATE (ivalues0, STAT = ERR)

        THROWM( ERR .NE. 0, "deallocation of ivalues0")

        DEALLOCATE (ivalues1, STAT = ERR)

        THROWM( ERR .NE. 0, "deallocation of ivalues1")

        DEALLOCATE (values0, STAT = ERR)

        THROWM( ERR .NE. 0, "deallocation of values0")

        DEALLOCATE (values1, STAT = ERR)

        THROWM( ERR .NE. 0, "deallocation of values1")

      ENDIF

    ENDIF

    
    READ (file_lu)section_present
    errMessage = 'io status while reading WATER_OPTICAL_CONSTANT'

    IF (section_present) THEN
!WATER_OPTICAL_CONSTANT                 !pb
      ALLOCATE (coef%woc_chn(coef%fmv_chn), STAT = ERR)

      THROWM( ERR .NE. 0, "allocation of woc_chn")

      ALLOCATE (coef%woc_cwn(coef%fmv_chn), STAT = ERR)

      THROWM( ERR .NE. 0, "allocation of woc_cwn")

      ALLOCATE (coef%woc_waopc_ow(coef%fmv_chn), STAT = ERR)

      THROWM( ERR .NE. 0, "allocation of woc_waopc_ow")

      ALLOCATE (coef%woc_waopc_fw(coef%fmv_chn), STAT = ERR)

      THROWM( ERR .NE. 0, "allocation of woc_waopc_fw")


      IF (all_channels) THEN
        READ (file_lu, iostat=ERR)coef%woc_chn, coef%woc_cwn, coef%woc_waopc_ow, coef%woc_waopc_fw!pb

        THROWM(err.ne.0,errMessage)

      ELSE
        ALLOCATE (ivalues0(coef%fmv_ori_nchn), STAT = ERR)

        THROWM( ERR .NE. 0, "allocation of ivalues0")

        ALLOCATE (values1(coef%fmv_ori_nchn), STAT = ERR)

        THROWM( ERR .NE. 0, "allocation of values1")

        ALLOCATE (values_c_0(coef%fmv_ori_nchn), STAT = ERR)

        THROWM( ERR .NE. 0, "allocation of values_c_0")

        ALLOCATE (values_c_1(coef%fmv_ori_nchn), STAT = ERR)

        THROWM( ERR .NE. 0, "allocation of values_c_1")

        READ (file_lu, iostat=ERR)ivalues0, values1, values_c_0, values_c_1

        THROWM(err.ne.0,errMessage)

        coef%woc_chn(:)      = ivalues0(channels(:))
        coef%woc_cwn(:)      = values1(channels(:))
        coef%woc_waopc_ow(:) = values_c_0(channels(:))
        coef%woc_waopc_fw(:) = values_c_1(channels(:))
        DEALLOCATE (ivalues0, STAT = ERR)

        THROWM( ERR .NE. 0, "deallocation of ivalues0")

        DEALLOCATE (values1, STAT = ERR)

        THROWM( ERR .NE. 0, "deallocation of values1")

        DEALLOCATE (values_c_0, STAT = ERR)

        THROWM( ERR .NE. 0, "deallocation of values_c_0")

        DEALLOCATE (values_c_1, STAT = ERR)

        THROWM( ERR .NE. 0, "deallocation of values_c_1")

      ENDIF

    ENDIF


    READ (file_lu)section_present
    errMessage = 'io status while reading WAVE_SPECTRUM'

    IF (section_present) THEN
!WAVE_SPECTRUM                          !pb
      READ (file_lu, iostat=ERR)coef%ws_nomega!pb

      THROWM(err.ne.0,errMessage)

      ALLOCATE (coef%ws_k_omega(coef%ws_nomega), STAT = ERR)

      THROWM( ERR .NE. 0, "allocation of ws_k_omega")

      ALLOCATE (coef%ws_npoint(coef%ws_nomega), STAT = ERR)

      THROWM( ERR .NE. 0, "allocation of ws_npoint")

      READ (file_lu, iostat=ERR)coef%ws_k_omega, coef%ws_npoint!pb

      THROWM(err.ne.0,errMessage)

    ENDIF
    
    errMessage = 'io status while reading FILTER_FUNCTIONS'

    IF (all_channels) THEN
      READ (file_lu, iostat=ERR)coef%ff_ori_chn, coef%ff_val_chn, coef%ff_cwn, coef%ff_bco, coef%ff_bcs, coef%ff_gam

      THROWM(err.ne.0,errMessage)

    ELSE
      ALLOCATE (ivalues0(coef%fmv_ori_nchn), STAT = ERR)

      THROWM( ERR .NE. 0, "allocation of ivalues0")

      ALLOCATE (ivalues1(coef%fmv_ori_nchn), STAT = ERR)

      THROWM( ERR .NE. 0, "allocation of ivalues1")

      ALLOCATE (values0(coef%fmv_ori_nchn), STAT = ERR)

      THROWM( ERR .NE. 0, "allocation of values0")

      ALLOCATE (values1(coef%fmv_ori_nchn), STAT = ERR)

      THROWM( ERR .NE. 0, "allocation of values1")

      ALLOCATE (values2(coef%fmv_ori_nchn), STAT = ERR)

      THROWM( ERR .NE. 0, "allocation of values2")

      ALLOCATE (values3(coef%fmv_ori_nchn), STAT = ERR)

      THROWM( ERR .NE. 0, "allocation of values3")

      READ (file_lu, iostat=ERR)ivalues0, ivalues1, values0, values1, values2, values3

      THROWM(err.ne.0,errMessage)

      coef%ff_ori_chn(:) = ivalues0(channels(:))
      coef%ff_val_chn(:) = ivalues1(channels(:))
      coef%ff_cwn(:)     = values0(channels(:))
      coef%ff_bco(:)     = values1(channels(:))
      coef%ff_bcs(:)     = values2(channels(:))
      coef%ff_gam(:)     = values3(channels(:))
      DEALLOCATE (ivalues0, STAT = ERR)

      THROWM( ERR .NE. 0, "deallocation of ivalues0")

      DEALLOCATE (ivalues1, STAT = ERR)

      THROWM( ERR .NE. 0, "deallocation of ivalues1")

      DEALLOCATE (values0, STAT = ERR)

      THROWM( ERR .NE. 0, "deallocation of values0")

      DEALLOCATE (values1, STAT = ERR)

      THROWM( ERR .NE. 0, "deallocation of values1")

      DEALLOCATE (values2, STAT = ERR)

      THROWM( ERR .NE. 0, "deallocation of values2")

      DEALLOCATE (values3, STAT = ERR)

      THROWM( ERR .NE. 0, "deallocation of values3")

    ENDIF

    errMessage = 'io status while reading FUNDAMENTAL_CONSTANTS'
    READ (file_lu, iostat=ERR)fc_speedl, coef%fc_planck_c1, coef%fc_planck_c2, coef%fc_sat_height

    THROWM(err.ne.0,errMessage)

    errMessage = 'io status while reading EMISSIVITY model versions'
    READ (file_lu, iostat=ERR)coef%fastem_ver, coef%ssirem_ver

    THROWM(err.ne.0,errMessage)

    errMessage = 'io status while reading FASTEM'

    IF (coef%fastem_ver >= 1) THEN
      ALLOCATE (coef%fastem_polar(coef%fmv_chn), STAT = ERR)

      THROWM( ERR .NE. 0, "allocation of fastem_polar")


      IF (all_channels) THEN
        READ (file_lu, iostat=ERR)coef%fastem_polar

        THROWM(err.ne.0,errMessage)

      ELSE
        ALLOCATE (ivalues0(coef%fmv_ori_nchn), STAT = ERR)

        THROWM( ERR .NE. 0, "allocation of ivalues0")

        READ (file_lu, iostat=ERR)ivalues0

        THROWM(err.ne.0,errMessage)

        coef%fastem_polar(:) = ivalues0(channels(:))
        DEALLOCATE (ivalues0, STAT = ERR)

        THROWM( ERR .NE. 0, "deallocation of ivalues0")

      ENDIF

    ENDIF

    errMessage = 'io status while reading SSIREM'

    IF (coef%ssirem_ver >= 1) THEN
      ALLOCATE (coef%ssirem_chn(coef%fmv_chn), STAT = ERR)

      THROWM( ERR .NE. 0, "allocation of ssirem_chn")

      ALLOCATE (coef%ssirem_a0(coef%fmv_chn), STAT = ERR)

      THROWM( ERR .NE. 0, "allocation of ssirem_a0")

      ALLOCATE (coef%ssirem_a1(coef%fmv_chn), STAT = ERR)

      THROWM( ERR .NE. 0, "allocation of ssirem_a1")

      ALLOCATE (coef%ssirem_a2(coef%fmv_chn), STAT = ERR)

      THROWM( ERR .NE. 0, "allocation of ssirem_a2")

      ALLOCATE (coef%ssirem_xzn1(coef%fmv_chn), STAT = ERR)

      THROWM( ERR .NE. 0, "allocation of ssirem_xzn1")

      ALLOCATE (coef%ssirem_xzn2(coef%fmv_chn), STAT = ERR)

      THROWM( ERR .NE. 0, "allocation of ssirem_xzn2")


      IF (all_channels) THEN
        READ (file_lu, iostat=ERR)     &
          & coef%ssirem_chn, coef%ssirem_a0, coef%ssirem_a1, coef%ssirem_a2, coef%ssirem_xzn1, coef%ssirem_xzn2

        THROWM(err.ne.0,errMessage)

      ELSE
        ALLOCATE (ivalues0(coef%fmv_ori_nchn), STAT = ERR)

        THROWM( ERR .NE. 0, "allocation of ivalues0")

        ALLOCATE (values0(coef%fmv_ori_nchn), STAT = ERR)

        THROWM( ERR .NE. 0, "allocation of values0")

        ALLOCATE (values1(coef%fmv_ori_nchn), STAT = ERR)

        THROWM( ERR .NE. 0, "allocation of values1")

        ALLOCATE (values2(coef%fmv_ori_nchn), STAT = ERR)

        THROWM( ERR .NE. 0, "allocation of values2")

        ALLOCATE (values3(coef%fmv_ori_nchn), STAT = ERR)

        THROWM( ERR .NE. 0, "allocation of values3")

        ALLOCATE (values4(coef%fmv_ori_nchn), STAT = ERR)

        THROWM( ERR .NE. 0, "allocation of values4")

        READ (file_lu, iostat=ERR)ivalues0, values0, values1, values2, values3, values4

        THROWM(err.ne.0,errMessage)

        coef%ssirem_chn(:)  = ivalues0(channels(:))
        coef%ssirem_a0(:)   = values0(channels(:))
        coef%ssirem_a1(:)   = values1(channels(:))
        coef%ssirem_a2(:)   = values2(channels(:))
        coef%ssirem_xzn1(:) = values3(channels(:))
        coef%ssirem_xzn2(:) = values4(channels(:))
        DEALLOCATE (ivalues0, STAT = ERR)

        THROWM( ERR .NE. 0, "deallocation of ivalues0")

        DEALLOCATE (values0, STAT = ERR)

        THROWM( ERR .NE. 0, "deallocation of values0")

        DEALLOCATE (values1, STAT = ERR)

        THROWM( ERR .NE. 0, "deallocation of values1")

        DEALLOCATE (values2, STAT = ERR)

        THROWM( ERR .NE. 0, "deallocation of values2")

        DEALLOCATE (values3, STAT = ERR)

        THROWM( ERR .NE. 0, "deallocation of values3")

        DEALLOCATE (values4, STAT = ERR)

        THROWM( ERR .NE. 0, "deallocation of values4")

      ENDIF

    ENDIF

    ALLOCATE (coef%ref_prfl_p(coef%fmv_lvl(gas_id_mixed)), STAT = ERR)

    THROWM( ERR .NE. 0, "allocation of ref_prfl_p")

    ALLOCATE (coef%ref_prfl_t(coef%fmv_lvl(gas_id_mixed), coef%fmv_gas), STAT = ERR)

    THROWM( ERR .NE. 0, "allocation of ref_prfl_t")

    ALLOCATE (coef%ref_prfl_mr(coef%fmv_lvl(gas_id_mixed), coef%fmv_gas), STAT = ERR)

    THROWM( ERR .NE. 0, "allocation of ref_prfl_mr")

    ALLOCATE (coef%lim_prfl_p(coef%fmv_lvl(gas_id_mixed)), STAT = ERR)

    THROWM( ERR .NE. 0, "allocation of lim_prfl_p")

    ALLOCATE (coef%lim_prfl_tmax(coef%fmv_lvl(gas_id_mixed)), STAT = ERR)

    THROWM( ERR .NE. 0, "allocation of lim_prfl_tmax")

    ALLOCATE (coef%lim_prfl_tmin(coef%fmv_lvl(gas_id_mixed)), STAT = ERR)

    THROWM( ERR .NE. 0, "allocation of lim_prfl_tmin")

    ALLOCATE (coef%lim_prfl_gmin(coef%fmv_lvl(gas_id_mixed), coef%fmv_gas), STAT = ERR)

    THROWM( ERR .NE. 0, "allocation of lim_prfl_gmin")

    ALLOCATE (coef%lim_prfl_gmax(coef%fmv_lvl(gas_id_mixed), coef%fmv_gas), STAT = ERR)

    THROWM( ERR .NE. 0, "allocation of lim_prfl_gmax")

    errMessage = 'io status while reading REFERENCE PROFILE'
    READ (file_lu, iostat=ERR)coef%ref_prfl_p, coef%ref_prfl_t, coef%ref_prfl_mr

    THROWM(err.ne.0,errMessage)

    errMessage = 'io status while reading PROFILE LIMITS'
    READ (file_lu, iostat=ERR)     &
      & coef%lim_prfl_p, coef%lim_prfl_tmax, coef%lim_prfl_tmin, coef%lim_prfl_gmax, coef%lim_prfl_gmin

    THROWM(err.ne.0,errMessage)

! FAST COEFFICIENT section

    ALLOCATE (coef%thermal, STAT = ERR)
    THROWM( ERR .NE. 0, "allocation of thermal fast coefs")
    NULLIFY (coef%thermal%mixedgas)
    NULLIFY (coef%thermal%watervapour)
    NULLIFY (coef%thermal%ozone)
    NULLIFY (coef%thermal%wvcont)
    NULLIFY (coef%thermal%co2)
    NULLIFY (coef%thermal%n2o)
    NULLIFY (coef%thermal%co)
    NULLIFY (coef%thermal%ch4)
    
    errMessage = 'io status while reading Mixed gases coefs'

    IF (coef%ncmixed > 0) THEN
      ALLOCATE (coef%thermal%mixedgas(coef%nlayers, coef%fmv_chn, coef%ncmixed), STAT = ERR)

      THROWM( ERR .NE. 0, "allocation of mixedgas")

      DO chn = 1, coef%fmv_ori_nchn
        IF (all_channels) THEN
          READ (file_lu, iostat=ERR)coef%thermal%mixedgas(:, chn, :)
        ELSE
          DO i = 1, coef%fmv_chn
            IF (chn == channels(i)) EXIT
          END DO
          IF (i > coef%fmv_chn) THEN
            READ (file_lu, iostat=ERR)
          ELSE
            READ (file_lu, iostat=ERR)coef%thermal%mixedgas(:, i, :)
          END IF
        END IF
        THROWM(err.ne.0,errMessage)
      ENDDO

    ENDIF

    errMessage = 'io status while reading Water vapour coefs'

    IF (coef%ncwater > 0) THEN
      ALLOCATE (coef%thermal%watervapour(coef%nlayers, coef%fmv_chn, coef%ncwater), STAT = ERR)

      THROWM( ERR .NE. 0, "allocation of watervapour")

      DO chn = 1, coef%fmv_ori_nchn
        IF (all_channels) THEN
          READ (file_lu, iostat=ERR)coef%thermal%watervapour(:, chn, :)
        ELSE
          DO i = 1, coef%fmv_chn
            IF (chn == channels(i)) EXIT
          END DO
          IF (i > coef%fmv_chn) THEN
            READ (file_lu, iostat=ERR)
          ELSE
            READ (file_lu, iostat=ERR)coef%thermal%watervapour(:, i, :)
          END IF
        END IF
        THROWM(err.ne.0,errMessage)
      ENDDO

    ENDIF

    errMessage = 'io status while reading Ozone coefs'

    IF (coef%ncozone > 0) THEN
      ALLOCATE (coef%thermal%ozone(coef%nlayers, coef%fmv_chn, coef%ncozone), STAT = ERR)

      THROWM( ERR .NE. 0, "allocation of ozone")

      DO chn = 1, coef%fmv_ori_nchn
        IF (all_channels) THEN
          READ (file_lu, iostat=ERR)coef%thermal%ozone(:, chn, :)
        ELSE
          DO i = 1, coef%fmv_chn
            IF (chn == channels(i)) EXIT
          END DO
          IF (i > coef%fmv_chn) THEN
            READ (file_lu, iostat=ERR)
          ELSE
            READ (file_lu, iostat=ERR)coef%thermal%ozone(:, i, :)
          END IF
        END IF
        THROWM(err.ne.0,errMessage)
      ENDDO

    ENDIF

    errMessage = 'io status while reading WV continuum coefs'

    IF (coef%ncwvcont > 0) THEN
      ALLOCATE (coef%thermal%wvcont(coef%nlayers, coef%fmv_chn, coef%ncwvcont), STAT = ERR)

      THROWM( ERR .NE. 0, "allocation of wvcont")

      DO chn = 1, coef%fmv_ori_nchn
        IF (all_channels) THEN
          READ (file_lu, iostat=ERR)coef%thermal%wvcont(:, chn, :)
        ELSE
          DO i = 1, coef%fmv_chn
            IF (chn == channels(i)) EXIT
          END DO
          IF (i > coef%fmv_chn) THEN
            READ (file_lu, iostat=ERR)
          ELSE
            READ (file_lu, iostat=ERR)coef%thermal%wvcont(:, i, :)
          END IF
        END IF
        THROWM(err.ne.0,errMessage)
      ENDDO

    ENDIF

    errMessage = 'io status while reading CO2 coefs'

    IF (coef%ncco2 > 0) THEN
      ALLOCATE (coef%thermal%co2(coef%nlayers, coef%fmv_chn, coef%ncco2), STAT = ERR)

      THROWM( ERR .NE. 0, "allocation of co2 ")

      DO chn = 1, coef%fmv_ori_nchn
        IF (all_channels) THEN
          READ (file_lu, iostat=ERR)coef%thermal%co2(:, chn, :)
        ELSE
          DO i = 1, coef%fmv_chn
            IF (chn == channels(i)) EXIT
          END DO
          IF (i > coef%fmv_chn) THEN
            READ (file_lu, iostat=ERR)
          ELSE
            READ (file_lu, iostat=ERR)coef%thermal%co2(:, i, :)
          END IF
        END IF
        THROWM(err.ne.0,errMessage)
      ENDDO

    ENDIF

    errMessage = 'io status while reading N2O coefs'

    IF (coef%ncn2o > 0) THEN
      ALLOCATE (coef%thermal%n2o(coef%nlayers, coef%fmv_chn, coef%ncn2o), STAT = ERR)

      THROWM( ERR .NE. 0, "allocation of n2o")

      DO chn = 1, coef%fmv_ori_nchn
        IF (all_channels) THEN
          READ (file_lu, iostat=ERR)coef%thermal%n2o(:, chn, :)
        ELSE
          DO i = 1, coef%fmv_chn
            IF (chn == channels(i)) EXIT
          END DO
          IF (i > coef%fmv_chn) THEN
            READ (file_lu, iostat=ERR)
          ELSE
            READ (file_lu, iostat=ERR)coef%thermal%n2o(:, i, :)
          END IF
        END IF
        THROWM(err.ne.0,errMessage)
      ENDDO

    ENDIF

    errMessage = 'io status while reading CO coefs'

    IF (coef%ncco > 0) THEN
      ALLOCATE (coef%thermal%co(coef%nlayers, coef%fmv_chn, coef%ncco), STAT = ERR)

      THROWM( ERR .NE. 0, "allocation of co")

      DO chn = 1, coef%fmv_ori_nchn
        IF (all_channels) THEN
          READ (file_lu, iostat=ERR)coef%thermal%co(:, chn, :)
        ELSE
          DO i = 1, coef%fmv_chn
            IF (chn == channels(i)) EXIT
          END DO
          IF (i > coef%fmv_chn) THEN
            READ (file_lu, iostat=ERR)
          ELSE
            READ (file_lu, iostat=ERR)coef%thermal%co(:, i, :)
          END IF
        END IF
        THROWM(err.ne.0,errMessage)
      ENDDO

    ENDIF

    errMessage = 'io status while reading CH4 coefs'

    IF (coef%ncch4 > 0) THEN
      ALLOCATE (coef%thermal%ch4(coef%nlayers, coef%fmv_chn, coef%ncch4), STAT = ERR)

      THROWM( ERR .NE. 0, "allocation of ch4")

      DO chn = 1, coef%fmv_ori_nchn
        IF (all_channels) THEN
          READ (file_lu, iostat=ERR)coef%thermal%ch4(:, chn, :)
        ELSE
          DO i = 1, coef%fmv_chn
            IF (chn == channels(i)) EXIT
          END DO
          IF (i > coef%fmv_chn) THEN
            READ (file_lu, iostat=ERR)
          ELSE
            READ (file_lu, iostat=ERR)coef%thermal%ch4(:, i, :)
          END IF
        END IF
        THROWM(err.ne.0,errMessage)
      ENDDO

    ENDIF

  
! SOLAR_FAST_COEFFICIENTS SECTION
! We need to know if they are present, and if so read them just as for FAST_COEFS
    READ (file_lu, iostat=ERR) coef%solarcoef

    IF (coef%solarcoef) THEN

      ALLOCATE (coef%solar, STAT = ERR)
      THROWM( ERR .NE. 0, "allocation of solar fast coefs")
      NULLIFY (coef%solar%mixedgas)
      NULLIFY (coef%solar%watervapour)
      NULLIFY (coef%solar%ozone)
      NULLIFY (coef%solar%wvcont)
      NULLIFY (coef%solar%co2)
      NULLIFY (coef%solar%n2o)
      NULLIFY (coef%solar%co)
      NULLIFY (coef%solar%ch4)
      
      errMessage = 'io status while reading solar Mixed gases coefs'
  
      IF (coef%ncmixed > 0) THEN
        ALLOCATE (coef%solar%mixedgas(coef%nlayers, coef%fmv_chn, coef%ncmixed), STAT = ERR)
  
        THROWM( ERR .NE. 0, "allocation of solar mixedgas")

        DO chn = 1, coef%fmv_ori_nchn
          IF (all_channels) THEN
            READ (file_lu, iostat=ERR)coef%solar%mixedgas(:, chn, :)
          ELSE
            DO i = 1, coef%fmv_chn
              IF (chn == channels(i)) EXIT
            END DO
            IF (i > coef%fmv_chn) THEN
              READ (file_lu, iostat=ERR)
            ELSE
              READ (file_lu, iostat=ERR)coef%solar%mixedgas(:, i, :)
            END IF
          END IF
          THROWM(err.ne.0,errMessage)
        ENDDO

      ENDIF
  
      errMessage = 'io status while reading solar Water vapour coefs'
  
      IF (coef%ncwater > 0) THEN
        ALLOCATE (coef%solar%watervapour(coef%nlayers, coef%fmv_chn, coef%ncwater), STAT = ERR)
  
        THROWM( ERR .NE. 0, "allocation of solar watervapour")
  
        DO chn = 1, coef%fmv_ori_nchn
          IF (all_channels) THEN
            READ (file_lu, iostat=ERR)coef%solar%watervapour(:, chn, :)
          ELSE
            DO i = 1, coef%fmv_chn
              IF (chn == channels(i)) EXIT
            END DO
            IF (i > coef%fmv_chn) THEN
              READ (file_lu, iostat=ERR)
            ELSE
              READ (file_lu, iostat=ERR)coef%solar%watervapour(:, i, :)
            END IF
          END IF
          THROWM(err.ne.0,errMessage)
        ENDDO
  
      ENDIF
  
      errMessage = 'io status while reading solar Ozone coefs'
  
      IF (coef%ncozone > 0) THEN
        ALLOCATE (coef%solar%ozone(coef%nlayers, coef%fmv_chn, coef%ncozone), STAT = ERR)
  
        THROWM( ERR .NE. 0, "allocation of solar ozone")

        DO chn = 1, coef%fmv_ori_nchn
          IF (all_channels) THEN
            READ (file_lu, iostat=ERR)coef%solar%ozone(:, chn, :)
          ELSE
            DO i = 1, coef%fmv_chn
              IF (chn == channels(i)) EXIT
            END DO
            IF (i > coef%fmv_chn) THEN
              READ (file_lu, iostat=ERR)
            ELSE
              READ (file_lu, iostat=ERR)coef%solar%ozone(:, i, :)
            END IF
          END IF
          THROWM(err.ne.0,errMessage)
        ENDDO

      ENDIF
  
      errMessage = 'io status while reading solar WV continuum coefs'
  
      IF (coef%ncwvcont > 0) THEN
        ALLOCATE (coef%solar%wvcont(coef%nlayers, coef%fmv_chn, coef%ncwvcont), STAT = ERR)
  
        THROWM( ERR .NE. 0, "allocation of solar wvcont")

        DO chn = 1, coef%fmv_ori_nchn
          IF (all_channels) THEN
            READ (file_lu, iostat=ERR)coef%solar%wvcont(:, chn, :)
          ELSE
            DO i = 1, coef%fmv_chn
              IF (chn == channels(i)) EXIT
            END DO
            IF (i > coef%fmv_chn) THEN
              READ (file_lu, iostat=ERR)
            ELSE
              READ (file_lu, iostat=ERR)coef%solar%wvcont(:, i, :)
            END IF
          END IF
          THROWM(err.ne.0,errMessage)
        ENDDO

      ENDIF
  
      errMessage = 'io status while reading solar CO2 coefs'
  
      IF (coef%ncco2 > 0) THEN
        ALLOCATE (coef%solar%co2(coef%nlayers, coef%fmv_chn, coef%ncco2), STAT = ERR)
  
        THROWM( ERR .NE. 0, "allocation of solar co2 ")

        DO chn = 1, coef%fmv_ori_nchn
          IF (all_channels) THEN
            READ (file_lu, iostat=ERR)coef%solar%co2(:, chn, :)
          ELSE
            DO i = 1, coef%fmv_chn
              IF (chn == channels(i)) EXIT
            END DO
            IF (i > coef%fmv_chn) THEN
              READ (file_lu, iostat=ERR)
            ELSE
              READ (file_lu, iostat=ERR)coef%solar%co2(:, i, :)
            END IF
          END IF
          THROWM(err.ne.0,errMessage)
        ENDDO

      ENDIF
  
      errMessage = 'io status while reading solar N2O coefs'
  
      IF (coef%ncn2o > 0) THEN
        ALLOCATE (coef%solar%n2o(coef%nlayers, coef%fmv_chn, coef%ncn2o), STAT = ERR)
  
        THROWM( ERR .NE. 0, "allocation of solar n2o")

        DO chn = 1, coef%fmv_ori_nchn
          IF (all_channels) THEN
            READ (file_lu, iostat=ERR)coef%solar%n2o(:, chn, :)
          ELSE
            DO i = 1, coef%fmv_chn
              IF (chn == channels(i)) EXIT
            END DO
            IF (i > coef%fmv_chn) THEN
              READ (file_lu, iostat=ERR)
            ELSE
              READ (file_lu, iostat=ERR)coef%solar%n2o(:, i, :)
            END IF
          END IF
          THROWM(err.ne.0,errMessage)
        ENDDO

      ENDIF
  
      errMessage = 'io status while reading solar CO coefs'
  
      IF (coef%ncco > 0) THEN
        ALLOCATE (coef%solar%co(coef%nlayers, coef%fmv_chn, coef%ncco), STAT = ERR)
  
        THROWM( ERR .NE. 0, "allocation of solar co")

        DO chn = 1, coef%fmv_ori_nchn
          IF (all_channels) THEN
            READ (file_lu, iostat=ERR)coef%solar%co(:, chn, :)
          ELSE
            DO i = 1, coef%fmv_chn
              IF (chn == channels(i)) EXIT
            END DO
            IF (i > coef%fmv_chn) THEN
              READ (file_lu, iostat=ERR)
            ELSE
              READ (file_lu, iostat=ERR)coef%solar%co(:, i, :)
            END IF
          END IF
          THROWM(err.ne.0,errMessage)
        ENDDO

      ENDIF
  
      errMessage = 'io status while reading solar CH4 coefs'
  
      IF (coef%ncch4 > 0) THEN
        ALLOCATE (coef%solar%ch4(coef%nlayers, coef%fmv_chn, coef%ncch4), STAT = ERR)
  
        THROWM( ERR .NE. 0, "allocation of solar ch4")

        DO chn = 1, coef%fmv_ori_nchn
          IF (all_channels) THEN
            READ (file_lu, iostat=ERR)coef%solar%ch4(:, chn, :)
          ELSE
            DO i = 1, coef%fmv_chn
              IF (chn == channels(i)) EXIT
            END DO
            IF (i > coef%fmv_chn) THEN
              READ (file_lu, iostat=ERR)
            ELSE
              READ (file_lu, iostat=ERR)coef%solar%ch4(:, i, :)
            END IF
          END IF
          THROWM(err.ne.0,errMessage)
        ENDDO

      ENDIF
    
    ENDIF ! ERR == 0 i.e. SOLAR_FAST_COEFS present
    

    READ (file_lu)section_present
    errMessage = 'io status while reading PLANCK_WEIGHTED'

    IF (section_present) THEN
!PLANCK_WEIGHTED
      ALLOCATE (coef%pw_chn(coef%fmv_chn), STAT = ERR)

      THROWM( ERR .NE. 0, "allocation of pw_chn")

      ALLOCATE (coef%pw_val_chn(coef%fmv_chn), STAT = ERR)

      THROWM( ERR .NE. 0, "allocation of pw_val_chn")

      IF (all_channels) THEN
        READ (file_lu, iostat=ERR)coef%pw_chn,coef%pw_val_chn

        THROWM(err.ne.0,errMessage)
      ELSE
        ALLOCATE (ivalues0(coef%fmv_ori_nchn), STAT = ERR)

        THROWM( ERR .NE. 0, "allocation of ivalues0")

        ALLOCATE (ivalues1(coef%fmv_ori_nchn), STAT = ERR)

        THROWM( ERR .NE. 0, "allocation of ivalues1")

        READ (file_lu, iostat=ERR)ivalues0, ivalues1

        THROWM(err.ne.0,errMessage)

        coef%pw_chn(:)            = ivalues0(channels(:))
        coef%pw_val_chn(:)        = ivalues1(channels(:))

        DEALLOCATE (ivalues0, STAT = ERR)

        THROWM( ERR .NE. 0, "deallocation of ivalues0")

        DEALLOCATE (ivalues1, STAT = ERR)

        THROWM( ERR .NE. 0, "deallocation of ivalues1")

      ENDIF

    ENDIF
    
    READ (file_lu)section_present
    errMessage = 'io status while reading README_SPECTRAL_RESPONSE_FUNCTION'

    IF (section_present) THEN
!README_SPECTRAL_RESPONSE_FUNCTION
      READ (file_lu, iostat=ERR)coef%readme_srf
      THROWM(err.ne.0,errMessage)

    ENDIF

    READ (file_lu) section_present
    errMessage = 'io status while reading NLTE_RADIANCE_COEFS'

    IF (section_present) THEN
!NLTE_RADIANCE_COEF

      coef%nltecoef = .TRUE. ! nlte coef present
      ALLOCATE(coef%nlte_coef)

      READ(file_lu, iostat=err) coef%nlte_coef%ncoef, &
      coef%nlte_coef%nsol, coef%nlte_coef%nsat, coef%nlte_coef%nchan, &
      coef%nlte_coef%start_chan

      NULLIFY (coef%nlte_coef%coef)
      NULLIFY (coef%nlte_coef%sol_zen_angle, coef%nlte_coef%cos_sol)
      NULLIFY (coef%nlte_coef%sat_zen_angle, coef%nlte_coef%sec_sat)

      ! For any monotonic channel selection we must find those selected channels
      ! which lie within the range of NLTE channels in the coef file. This
      ! constitutes another contiguous block of channels in the coef structure.
      IF (.NOT. all_channels) THEN
        ALLOCATE(nlte_chans(SIZE(channels)))
        nlte_count = 0
        nlte_start = 0
        DO i = 1, SIZE(channels)
          IF (i > 1) THEN
            IF (channels(i) < channels(i-1)) THEN
              err = errorstatus_fatal
              THROWM( ERR .NE. 0, "non-monotonic channel selection incompatible with NLTE coefficients")
            ENDIF
          ENDIF
          IF (channels(i) >= coef%nlte_coef%start_chan .AND. &
              channels(i) < coef%nlte_coef%start_chan + coef%nlte_coef%nchan) THEN
            nlte_count = nlte_count + 1
            nlte_chans(nlte_count) = channels(i) - coef%nlte_coef%start_chan + 1
            IF (nlte_count == 1) nlte_start = i
          ENDIF
        ENDDO

        IF (nlte_count > 0) THEN
          nlte_file_nchan           = coef%nlte_coef%nchan
          coef%nlte_coef%start_chan = nlte_start
          coef%nlte_coef%nchan      = nlte_count
        ELSE
          coef%nlte_coef%start_chan = 0_jpim
          coef%nlte_coef%nchan = 0_jpim
          DEALLOCATE(coef%nlte_coef, nlte_chans)
          coef%nltecoef = .FALSE.
        ENDIF
      ENDIF

      IF (coef%nltecoef) THEN

        ALLOCATE(coef%nlte_coef%sol_zen_angle(coef%nlte_coef%nsol), STAT = ERR)
        ALLOCATE(coef%nlte_coef%sec_sat(coef%nlte_coef%nsat), STAT = ERR)
        ALLOCATE(coef%nlte_coef%coef(coef%nlte_coef%ncoef, coef%nlte_coef%nsat, &
                 coef%nlte_coef%nsol, coef%nlte_coef%nchan), STAT = ERR)

        THROWM( ERR .NE. 0, "allocation of nlte_coef")

        READ(file_lu, iostat=ERR) coef%nlte_coef%sol_zen_angle
        READ(file_lu, iostat=ERR) coef%nlte_coef%sec_sat

        IF (all_channels) THEN
          READ(file_lu, iostat=ERR) coef%nlte_coef%coef
        ELSE
          ALLOCATE(nlte_values(coef%nlte_coef%ncoef, coef%nlte_coef%nsat, &
                   coef%nlte_coef%nsol, nlte_file_nchan), STAT = ERR)
          THROWM( ERR .NE. 0, "allocation of nlte_values array")
          READ(file_lu, iostat=ERR) nlte_values
          coef%nlte_coef%coef(:,:,:,:) = nlte_values(:,:,:,nlte_chans(1:nlte_count))
          DEALLOCATE(nlte_chans, nlte_values)
        ENDIF

      ENDIF
    ENDIF

    READ (file_lu) section_present
    errMessage = 'io status while reading PRESSURE_MODULATED_CELL'

    IF (section_present) THEN
!PRESSURE_MODULATED_CELL

      coef%pmc_shift = .TRUE. ! pmc_shift coef present

      READ (file_lu, iostat=err) coef%pmc_lengthcell

      ALLOCATE (coef%pmc_pnominal(coef%fmv_chn), STAT = ERR)
      THROWM( ERR .NE. 0, "allocation of coef%pmc_pnominal array")
      IF (all_channels) THEN
        READ (file_lu, iostat=ERR) coef%pmc_pnominal
      ELSE
        ALLOCATE (values1(coef%fmv_ori_nchn), STAT = ERR)
        THROWM( ERR .NE. 0, "allocation of values1")
        READ (file_lu, iostat=ERR) values1
        coef%pmc_pnominal(:) = values1(channels(:))
        DEALLOCATE (values1, STAT = ERR)
        THROWM( ERR .NE. 0, "deallocation of values1")
      ENDIF

      READ (file_lu, iostat=err) coef%pmc_tempcell, coef%pmc_betaplus1, &
          coef%pmc_nlay, coef%pmc_nvar

      ALLOCATE (coef%pmc_coef(coef%pmc_nlay, coef%fmv_chn, coef%pmc_nvar), STAT = ERR)
      THROWM( ERR .NE. 0, "allocation of coef%pmc_coef array")
 
      DO chn = 1, coef%fmv_ori_nchn
        IF (all_channels) THEN
          READ (file_lu, iostat=ERR) coef%pmc_coef(:, chn, :)
        ELSE
          DO i = 1, coef%fmv_chn
            IF (chn == channels(i)) EXIT
          END DO
          IF (i > coef%fmv_chn) THEN
            READ (file_lu, iostat=ERR)
          ELSE
            READ (file_lu, iostat=ERR)coef%pmc_coef(:, i, :)
          END IF
        END IF
        THROWM(err.ne.0,errMessage)
      ENDDO

    ENDIF

!
! Here add reading of new sections for binary format in order to keep compatibility with
! previous versions
!
  CATCH
END SUBROUTINE rttov_read_binary_coef
