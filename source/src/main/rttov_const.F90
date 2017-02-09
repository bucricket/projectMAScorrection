! Description:
!> @file
!!   Defines various constants used by RTTOV
!
!> @brief
!!   Defines various constants used by RTTOV
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
!    Copyright 2015, EUMETSAT, All Rights Reserved.
!
MODULE rttov_const

  USE parkind1, ONLY : jpim, jprb
  IMPLICIT NONE

  ! General
  ! -------
  ! Version number of the current code
  INTEGER(KIND=jpim), PARAMETER :: version = 11
  INTEGER(KIND=jpim), PARAMETER :: release = 3
  INTEGER(KIND=jpim), PARAMETER :: minor_version = 0

  ! Min/max version numbers of compatible coefficient files:
  !   coef files with "id_comp_lvl" outside range will be rejected
  INTEGER(KIND=jpim), PARAMETER :: version_compatible_min = 10
  INTEGER(KIND=jpim), PARAMETER :: version_compatible_max = 11

  CHARACTER(LEN=16),  PARAMETER :: rttov_magic_string = '%RTTOV_COEFF    '
  REAL(KIND=jprb),    PARAMETER :: rttov_magic_number = 1.2345E+12_jprb

  INTEGER(KIND=jpim), PARAMETER :: default_err_unit = 0  ! standard error unit number
                              ! standard error unit number is 7 for HPUX

  ! Precision and numerical constants
  ! ---------------------------------
  ! Try to ensure this is large enough to avoid overflows in reciprocals but small enough to not affect the calculations.
  ! these parameters are defined at bottom of module (because they use
  REAL(jprb), PARAMETER :: max_exp_exponent = 50._jprb ! approx 1e22 which should be sufficiently big for most purposes
  REAL(jprb), PARAMETER :: min_exponent     = 1e-16_jprb ! approx log_10(1+2^-52) - anything raised to this power or smaller
                                                         ! should be approx equal to 1
! small_val is defined in rttov_transmit to avoid compiler incompatibility
!   ! small_val is used in rttov_transmit to ensure small values do not result in underflows. In subsequent calculations
!   ! these values are multiplied together hence the exponent of 1/3.
!   REAL(jprb)            :: small_val = (tiny(min_exponent)) ** (0.333333_jprb) ! XLF doesn't like 1/3

  ! Physical constants
  ! ------------------
  ! Molecular weights  (g/mole) are calculated by adding NIST Standard Atomic Weights
  ! Molecular weight of dry air refers to US standard atmosphere 1976
  ! NIST  Standard Atomic Weight are:
  ! H    1.00794   (7)
  ! C   12.0107    (8)
  ! N   14.0067    (2)
  ! O   15.9994    (3)
  ! http://webbook.nist.gov/chemistry/form-ser.html
  REAL(KIND=jprb), PARAMETER :: mair   = 28.9644_jprb
  REAL(KIND=jprb), PARAMETER :: mh2o   = 18.01528_jprb
  REAL(KIND=jprb), PARAMETER :: mo3    = 47.9982_jprb
  REAL(KIND=jprb), PARAMETER :: mco2   = 44.0095_jprb
  REAL(KIND=jprb), PARAMETER :: mch4   = 16.04246_jprb
  REAL(KIND=jprb), PARAMETER :: mn2o   = 44.0128_jprb
  REAL(KIND=jprb), PARAMETER :: mco    = 28.0101_jprb
  REAL(KIND=jprb), PARAMETER :: mso2   = 64.064_jprb
  REAL(KIND=jprb), PARAMETER :: mo2    = 31.9988_jprb
  REAL(KIND=jprb), PARAMETER :: mno    = 30.0061_jprb
  REAL(KIND=jprb), PARAMETER :: mno2   = 46.0055_jprb
  REAL(KIND=jprb), PARAMETER :: mhno3  = 63.0128_jprb
  REAL(KIND=jprb), PARAMETER :: mocs   = 60.075_jprb
  REAL(KIND=jprb), PARAMETER :: mn2    = 28.0134_jprb
  REAL(KIND=jprb), PARAMETER :: mccl4  = 153.823_jprb
  REAL(KIND=jprb), PARAMETER :: mcfc11 = 137.368_jprb ! CFCl3
  REAL(KIND=jprb), PARAMETER :: mcfc12 = 120.914_jprb ! CF2Cl2
  REAL(KIND=jprb), PARAMETER :: mcfc14 = 88.0043_jprb ! CF4


  ! NB Included in v11.3 for compatibility where other software uses these: to be removed in v12.
  !    Simple units conversion from  mixing ratio to ppmv
  REAL(KIND=jprb), PARAMETER :: q_mixratio_to_ppmv  = 1.60771704e+6_jprb
  REAL(KIND=jprb), PARAMETER :: o3_mixratio_to_ppmv = 6.03504e+5_jprb
  REAL(KIND=jprb), PARAMETER :: co2_mixratio_to_ppmv= 6.58114e+5_jprb
  REAL(KIND=jprb), PARAMETER :: co_mixratio_to_ppmv = 1.0340699e+6_jprb
  REAL(KIND=jprb), PARAMETER :: n2o_mixratio_to_ppmv= 6.58090e+5_jprb
  REAL(KIND=jprb), PARAMETER :: ch4_mixratio_to_ppmv= 1.80548e+6_jprb


  ! Avogadro constant from NIST (mol-1)
  REAL(KIND=jprb), PARAMETER :: na = 6.02214129E23_jprb
  
  ! Gravity from NIST 9.80665 ms-1 (exact)
  REAL(KIND=jprb), PARAMETER :: gravity = 9.80665_jprb

  ! Fundamental constants taken from http://physics.nist.gov/cuu/index.html
  ! Barry N. Taylor (Fundamental Constants Data Center of NIST) and Peter J.
  ! Mohr (Atomic Physics Division of NIST). Values as of 01/12/2010.
  ! c1 = 2hc**2; c2 = hc/k; units are consistent with those used within RTTOV.
  ! NB Planck fn constant values (c1, c2) in core code are still taken from coef files.
  !    Newer values are available from NIST: we will update in a future version.
  REAL(KIND=jprb), PARAMETER :: c1 = .00001191042722_jprb   ! * mW/(m2 sr cm-4)
  REAL(KIND=jprb), PARAMETER :: c2 = 1.4387752_jprb         ! cm K
  REAL(KIND=jprb), PARAMETER :: speedl = 29979245800.0_jprb ! Speed of light cm s-1

  !
  ! Kaye & Laby latest library edition is 16e 1995, and gives
  ! * standard value  g = 9.80665 ms-1 exactly (p.191)
  ! * earth mean radius r= 6371.00 km (p191)
  !    [defined as [(r_equator)^2 (r_pole)]^1/3]
  REAL(KIND=jprb), PARAMETER :: pi      = 3.1415926535_jprb
  REAL(KIND=jprb), PARAMETER :: deg2rad = pi/180.0_jprb
  REAL(KIND=jprb), PARAMETER :: earthradius = 6371.00_jprb
  REAL(KIND=jprb), PARAMETER :: flatt       = 3.3528107E-3_jprb
  REAL(KIND=jprb), PARAMETER :: omega       = 7292115E-11_jprb
  REAL(KIND=jprb), PARAMETER :: eqrad       = 6378.137_jprb
  REAL(KIND=jprb), PARAMETER :: grave       = 9.7803267715_jprb
  REAL(KIND=jprb), PARAMETER :: z4pi_r      = 0.0795774715_jprb
  REAL(KIND=jprb), PARAMETER :: pi_r        = 0.3183098862_jprb
  REAL(KIND=jprb), PARAMETER :: sec_theta_eff = 1.743446796_jprb   ! theta_eff = 55 degrees

  ! The Cosmic Microwave Background Spectrum from the Full COBE FIRAS Data Set
  ! Fixsen D.J. et al Astrophysical Journal v.473, p.576 December 1996
  ! CMBR = 2.728 +- 0.004K
  REAL(KIND=jprb), PARAMETER :: tcosmic     = 2.728_jprb
  !  REAL(KIND=jprb), PARAMETER :: tcosmic     = 0.1_jprb !used for ECMWF tests

  ! Universal gas constant R = 8.314510 J/mol/K
  REAL(KIND=jprb), PARAMETER :: rgp = 8.314510_jprb
  REAL(KIND=jprb), PARAMETER :: rgc = 8.314472_jprb

  ! mean molar mass of dry air rm = 0.0289644 kg.mol^-1
  REAL(KIND=jprb), PARAMETER :: rm = 0.0289644_jprb

  ! zero temperature(K)
  REAL(KIND=jprb), PARAMETER :: t0 = 273.15_jprb
  ! standard pressure
  REAL(KIND=jprb), PARAMETER :: p0 = 1013.25_jprb


  ! Satellite and instrument information
  ! ------------------------------------

  ! Platform ID codes
  INTEGER(KIND=jpim), PARAMETER :: nplatforms = 45
  ! NB ID 22 was changed from "insat_3d" - use ID 40 for INSAT3
  !    IDs 41 and 45 are for ground-based and airborne sensors (experimental)
  INTEGER(KIND=jpim), PARAMETER :: &
       & platform_id_noaa      = 1,  platform_id_dmsp      = 2,  platform_id_meteosat  = 3,  &
       & platform_id_goes      = 4,  platform_id_gms       = 5,  platform_id_fy2       = 6,  &
       & platform_id_trmm      = 7,  platform_id_ers       = 8,  platform_id_eos       = 9,  &
       & platform_id_metop     = 10, platform_id_envisat   = 11, platform_id_msg       = 12, &
       & platform_id_fy1       = 13, platform_id_adeos     = 14, platform_id_mtsat     = 15, &
       & platform_id_coriolis  = 16, platform_id_jpss      = 17, platform_id_gifts     = 18, &
       & platform_id_sentinel3 = 19, platform_id_meghatr   = 20, platform_id_kalpana   = 21, &
       & platform_id_meteor    = 22, platform_id_fy3       = 23, platform_id_coms      = 24, &
       & platform_id_meteorm   = 25, platform_id_gosat     = 26, platform_id_calipso   = 27, &
       & platform_id_dummy     = 28, platform_id_gcomw     = 29, platform_id_nimbus    = 30, &
       & platform_id_himawari  = 31, platform_id_mtg       = 32, platform_id_saral     = 33, &
       & platform_id_metopsg   = 34, platform_id_landsat   = 35, platform_id_jason     = 36, &
       & platform_id_gpm       = 37, platform_id_insat1    = 38, platform_id_insat2    = 39, &
       & platform_id_insat3    = 40, platform_id_ground    = 41, platform_id_dscovr    = 42, &
       & platform_id_clarreo   = 43, platform_id_ticfire   = 44, platform_id_aircraft  = 45

  ! Platform names
  CHARACTER(LEN=9), PARAMETER :: platform_name(nplatforms) = &
       & (/ 'noaa     ', 'dmsp     ', 'meteosat ', 'goes     ', 'gms      ', &
          & 'fy2      ', 'trmm     ', 'ers      ', 'eos      ', 'metop    ', &
          & 'envisat  ', 'msg      ', 'fy1      ', 'adeos    ', 'mtsat    ', &
          & 'coriolis ', 'jpss     ', 'gifts    ', 'sentinel3', 'meghatr  ', &
          & 'kalpana  ', 'meteor   ', 'fy3      ', 'coms     ', 'meteor-m ', &
          & 'gosat    ', 'calipso  ', 'dummy    ', 'gcom-w   ', 'nimbus   ', &
          & 'himawari ', 'mtg      ', 'saral    ', 'metopsg  ', 'landsat  ', &
          & 'jason    ', 'gpm      ', 'insat1   ', 'insat2   ', 'insat3   ', &
          & 'ground   ', 'dscovr   ', 'clarreo  ', 'ticfire  ', 'aircraft ' /)

  ! Instrument ID codes
  INTEGER(KIND=jpim), PARAMETER :: &
       & inst_id_hirs   =  0, inst_id_msu      =  1, inst_id_ssu    =  2, inst_id_amsua   =  3, &
       & inst_id_amsub  =  4, inst_id_avhrr    =  5, inst_id_ssmi   =  6, inst_id_vtpr1   =  7, &
       & inst_id_vtpr2  =  8, inst_id_tmi      =  9, inst_id_ssmis  = 10, inst_id_airs    = 11, &
       & inst_id_hsb    = 12, inst_id_modis    = 13, inst_id_atsr   = 14, inst_id_mhs     = 15, &
       & inst_id_iasi   = 16, inst_id_amsre    = 17, inst_id_gmsim  = 18, inst_id_atms    = 19, &
       & inst_id_mviri  = 20, inst_id_seviri   = 21, inst_id_goesim = 22, inst_id_goessd  = 23, &
       & inst_id_mtsatim= 24, inst_id_vissr    = 25, inst_id_mvisr  = 26, inst_id_cris    = 27, &
       & inst_id_cmis   = 28, inst_id_viirs    = 29, inst_id_windsat= 30, inst_id_gifts   = 31, &
       & inst_id_ssmt1  = 32, inst_id_ssmt2    = 33, inst_id_saphir = 34, inst_id_madras  = 35, &
       & inst_id_ssmisz = 36, inst_id_vhrr     = 37, inst_id_insatim= 38, inst_id_insatsd = 39, &
       & inst_id_mwts   = 40, inst_id_mwhs     = 41, inst_id_iras   = 42, inst_id_mwri    = 43, &
       & inst_id_abi    = 44, inst_id_mi       = 45, inst_id_msumr  = 46, inst_id_tansofts= 47, &
       & inst_id_iir    = 48, inst_id_mwr      = 49, inst_id_dummyir= 50, inst_id_dummymw = 51, &
       & inst_id_dummyhi= 52, inst_id_dummypo  = 53, inst_id_scams  = 54, inst_id_smmr    = 55, &
       & inst_id_ahi    = 56, inst_id_irs      = 57, inst_id_altika = 58, inst_id_iasing  = 59, &
       & inst_id_tm     = 60, inst_id_fci      = 61, inst_id_amsr1  = 62, inst_id_amsr2   = 63, &
       & inst_id_vissr2 = 64, inst_id_slstr    = 65, inst_id_tirs   = 66, inst_id_amr     = 67, &
       & inst_id_oli    = 68, inst_id_iris     = 69, inst_id_ici    = 70, inst_id_gmi     = 71, &
       & inst_id_mwts2  = 72, inst_id_mwhs2    = 73, inst_id_aster  = 74, inst_id_hatpro  = 75, &
       & inst_id_mtvzagy= 76, inst_id_metimage = 77, inst_id_mws    = 78, inst_id_mwi     = 79, &
       & inst_id_epic   = 80, inst_id_mrir     = 81, inst_id_si     = 82, inst_id_mrfirs  = 83, &
       & inst_id_mbfiri = 84, inst_id_lhr      = 85, inst_id_ismar  = 86, inst_id_mersi1  = 87, &
       & inst_id_mersi2 = 88

  INTEGER(KIND=jpim), PARAMETER :: ninst = 89
  ! List of instruments  !!!! HIRS is number 0
  CHARACTER(LEN=8), PARAMETER :: inst_name(0:ninst-1) =       &
        & (/ 'hirs    ', 'msu     ', 'ssu     ', 'amsua   ', 'amsub   ',  &
           & 'avhrr   ', 'ssmi    ', 'vtpr1   ', 'vtpr2   ', 'tmi     ',  &
           & 'ssmis   ', 'airs    ', 'hsb     ', 'modis   ', 'atsr    ',  &
           & 'mhs     ', 'iasi    ', 'amsre   ', 'imager  ', 'atms    ',  &
           & 'mviri   ', 'seviri  ', 'imager  ', 'sounder ', 'imager  ',  &
           & 'vissr   ', 'mvisr   ', 'cris    ', 'cmis    ', 'viirs   ',  &
           & 'windsat ', 'gifts   ', 'ssmt1   ', 'ssmt2   ', 'saphir  ',  &
           & 'madras  ', 'ssmisz  ', 'vhrr    ', 'imager  ', 'sounder ',  &
           & 'mwts    ', 'mwhs    ', 'iras    ', 'mwri    ', 'abi     ',  &
           & 'mi      ', 'msumr   ', 'tansofts', 'iir     ', 'mwr     ',  &
           & 'dummyir ', 'dummymw ', 'dummyhi ', 'dummypo ', 'scams   ',  &
           & 'smmr    ', 'ahi     ', 'irs     ', 'altika  ', 'iasing  ',  &
           & 'tm      ', 'fci     ', 'amsr    ', 'amsr2   ', 'vissr   ',  &
           & 'slstr   ', 'tirs    ', 'amr     ', 'oli     ', 'iris    ',  &
           & 'ici     ', 'gmi     ', 'mwts2   ', 'mwhs2   ', 'aster   ',  &
           & 'hatpro  ', 'mtvzagy ', 'metimage', 'mws     ', 'mwi     ',  &
           & 'epic    ', 'mrir    ', 'si      ', 'mrfirs  ', 'mbfiri  ',  &
           & 'lhr     ', 'ismar   ', 'mersi1  ', 'mersi2  ' /)

  ! Sensor type ID codes
  INTEGER(KIND=jpim), PARAMETER :: nsensors = 4
  INTEGER(KIND=jpim), PARAMETER :: &
       & sensor_id_ir     = 1, &
       & sensor_id_mw     = 2, &
       & sensor_id_hi     = 3, &
       & sensor_id_po     = 4

  ! Sensor type names
  CHARACTER(LEN=2), PARAMETER :: sensor_name(nsensors) = &
       & (/ 'ir', 'mw', 'hi', 'po' /)

  ! Sensor types corresponding to entries in the inst_name array
  INTEGER(KIND=jpim), PARAMETER :: sensor_id(0:ninst-1) = (/ &
    sensor_id_ir, sensor_id_mw, sensor_id_ir, sensor_id_mw, sensor_id_mw,  &
    sensor_id_ir, sensor_id_mw, sensor_id_ir, sensor_id_ir, sensor_id_mw,  &
    sensor_id_mw, sensor_id_hi, sensor_id_mw, sensor_id_ir, sensor_id_ir,  &
    sensor_id_mw, sensor_id_hi, sensor_id_mw, sensor_id_ir, sensor_id_mw,  &
    sensor_id_ir, sensor_id_ir, sensor_id_ir, sensor_id_ir, sensor_id_ir,  &
    sensor_id_ir, sensor_id_ir, sensor_id_hi, sensor_id_mw, sensor_id_ir,  &
    sensor_id_po, sensor_id_hi, sensor_id_mw, sensor_id_mw, sensor_id_mw,  &
    sensor_id_mw, sensor_id_mw, sensor_id_ir, sensor_id_ir, sensor_id_ir,  &
    sensor_id_mw, sensor_id_mw, sensor_id_ir, sensor_id_mw, sensor_id_ir,  &
    sensor_id_ir, sensor_id_ir, sensor_id_hi, sensor_id_ir, sensor_id_mw,  &
    sensor_id_ir, sensor_id_mw, sensor_id_hi, sensor_id_po, sensor_id_mw,  &
    sensor_id_mw, sensor_id_ir, sensor_id_hi, sensor_id_mw, sensor_id_hi,  &
    sensor_id_ir, sensor_id_ir, sensor_id_mw, sensor_id_mw, sensor_id_ir,  &
    sensor_id_ir, sensor_id_ir, sensor_id_mw, sensor_id_ir, sensor_id_hi,  &
    sensor_id_mw, sensor_id_mw, sensor_id_mw, sensor_id_mw, sensor_id_ir,  &
    sensor_id_mw, sensor_id_mw, sensor_id_ir, sensor_id_mw, sensor_id_mw,  &
    sensor_id_ir, sensor_id_ir, sensor_id_hi, sensor_id_hi, sensor_id_ir,  &
    sensor_id_hi, sensor_id_mw, sensor_id_ir, sensor_id_ir /)


  ! Coefficient file Section names
  ! ------------------------------
  INTEGER(KIND=jpim), PARAMETER :: nsections = 43
  INTEGER(KIND=jpim), PARAMETER :: lensection = 34
  CHARACTER(LEN=lensection), PARAMETER :: section_types(nsections) = &
    & (/ 'IDENTIFICATION                    ', 'LINE-BY-LINE                      ', &
       & 'FAST_MODEL_VARIABLES              ', 'FILTER_FUNCTIONS                  ', &
       & 'FUNDAMENTAL_CONSTANTS             ', 'SSIREM                            ', &
       & 'FASTEM                            ', 'REFERENCE_PROFILE                 ', &
       & 'PROFILE_LIMITS                    ', 'FAST_COEFFICIENTS                 ', &
       & 'COEF_SUB_FILES                    ', 'GAZ_UNITS                         ', &
       & 'DIMENSIONS                        ', 'FREQUENCIES                       ', &
       & 'HYDROMETEOR                       ', 'CONVERSIONS                       ', &
       & 'EXTINCTION                        ', 'ALBEDO                            ', &
       & 'ASYMMETRY                         ', 'GAS_SPECTRAL_INTERVAL             ', &
       & 'TRANSMITTANCE_TRESHOLD            ', 'SOLAR_SPECTRUM                    ', &
       & 'WATER_OPTICAL_CONSTANT            ', 'WAVE_SPECTRUM                     ', &
       & 'AEROSOLS_PARAMETERS               ', 'AEROSOLS_COMPONENTS               ', &
       & 'WATERCLOUD_TYPES                  ', 'WATERCLOUD_PARAMETERS             ', &
       & 'ICECLOUD_TYPES                    ', 'HEXAGONAL_PARAMETERS              ', &
       & 'AGGREGATE_PARAMETERS              ', 'PRINCOMP_PREDICTORS               ', &
       & 'PRINCOMP_EIGENVECTORS             ', 'PRINCOMP_COEFFICIENTS             ', &
       & 'EMISSIVITY_COEFFICIENTS           ', 'PC_REFERENCE_PROFILE              ', &
       & 'PC_PROFILE_LIMITS                 ', 'INSTRUMENT_NOISE                  ', &
       & 'PLANCK_WEIGHTED                   ', 'SOLAR_FAST_COEFFICIENTS           ', &
       & 'README_SPECTRAL_RESPONSE_FUNCTION ', 'NLTE_RADIANCE_COEFS               ', &
       & 'PRESSURE_MODULATED_CELL           '/)

  ! Gas ID codes
  INTEGER(KIND=jpim), PARAMETER :: ngases_max = 8
  INTEGER(KIND=jpim), PARAMETER :: &
        & gas_id_mixed       = 1,  &
        & gas_id_watervapour = 2,  &
        & gas_id_ozone       = 3,  &
        & gas_id_wvcont      = 4,  &
        & gas_id_co2         = 5,  &
        & gas_id_n2o         = 6,  &
        & gas_id_co          = 7,  &
        & gas_id_ch4         = 8!,  &
!         & gas_id_so2         = 9,  &
!         & gas_id_o2          = 10, &
!         & gas_id_no          = 11, &
!         & gas_id_no2         = 12, &
!         & gas_id_hno3        = 13, &
!         & gas_id_ocs         = 14, &
!         & gas_id_n2          = 15, &
!         & gas_id_ccl4        = 16, &
!         & gas_id_cfc11       = 17, &
!         & gas_id_cfc12       = 18, &
!         & gas_id_cfc14       = 19

  ! Gas names
  CHARACTER(LEN=12), PARAMETER :: gas_name(ngases_max) = &
        & (/ 'Mixed_gases ', &
           & 'Water_vapour', &
           & 'Ozone       ', &
           & 'WV_Continuum', &
           & 'CO2         ', &
           & 'N2O         ', &
           & 'CO          ', &
           & 'CH4         ' /) !, &
!            & 'SO2         ', &
!            & 'O2          ', &
!            & 'NO          ', &
!            & 'NO2         ', &
!            & 'HNO3        ', &
!            & 'OCS         ', &
!            & 'N2          ', &
!            & 'CCl4        ', &
!            & 'CFC-11      ', &
!            & 'CFC-12      ', &
!            & 'CFC-14      ' /)

  ! Convenient array of gas molecular masses
  REAL(KIND=jprb), PARAMETER :: gas_mass(ngases_max) = &
    & (/ mair, mh2o, mo3, mh2o, mco2, mn2o, mco, mch4 /) !, mso2, mo2, mno, &
!          mno2, mhno3, mocs, mn2, mccl4, mcfc11, mcfc12, mcfc14 /)

  ! Gas units
  INTEGER(KIND=jpim), PARAMETER :: ngases_unit = 2
  INTEGER(KIND=jpim), PARAMETER :: &
        & gas_unit_compatibility = 0, &   ! Used for input gas units to ensure behaviour like v11.2 and earlier
        & gas_unit_specconc  = 1, &     ! specific concentration (kg/kg over wet air)
        & gas_unit_ppmv      = 2        ! volume mixing ratio (ppmv over wet air)
  CHARACTER(LEN=12), PARAMETER :: gas_unit_name(ngases_unit) = &
        & (/ 'spec. concen', &
           & 'ppmv        '  /)

  ! Error reporting
  ! ---------------
  INTEGER(KIND=jpim), PARAMETER :: errorstatus_success = 0
  INTEGER(KIND=jpim), PARAMETER :: errorstatus_fatal   = 1

  ! Surface types
  ! -------------
  INTEGER(KIND=jpim), PARAMETER :: nsurftype = 2
  INTEGER(KIND=jpim), PARAMETER :: surftype_land = 0
  INTEGER(KIND=jpim), PARAMETER :: surftype_sea = 1
  INTEGER(KIND=jpim), PARAMETER :: surftype_seaice = 2

  ! Water types
  ! -----------
  INTEGER(KIND=jpim), PARAMETER :: nwatertype = 1
  INTEGER(KIND=jpim), PARAMETER :: watertype_fresh_water = 0
  INTEGER(KIND=jpim), PARAMETER :: watertype_ocean_water = 1

  ! Hard limits for control of input profile
  ! ----------------------------------------
  ! Temperature
  REAL(KIND=jprb), PARAMETER :: tmax   = 400.0_jprb       ! degK
  REAL(KIND=jprb), PARAMETER :: tmin   = 90.0_jprb        ! degK
  ! Water Vapour
  REAL(KIND=jprb), PARAMETER :: qmax   = 0.60E+06_jprb    ! ppmv 0.373_jprb kg/kg
  REAL(KIND=jprb), PARAMETER :: qmin   = 0.1E-10_jprb     ! ppmv
  ! Ozone
  REAL(KIND=jprb), PARAMETER :: o3max  = 1000.0_jprb      ! ppmv  1.657E-3_jprb kg/kg
  REAL(KIND=jprb), PARAMETER :: o3min  = 0.1E-10_jprb     ! ppmv
  ! CO2
  REAL(KIND=jprb), PARAMETER :: co2max = 1000.0_jprb      ! ppmv
  REAL(KIND=jprb), PARAMETER :: co2min = 0.1E-10_jprb     ! ppmv
  ! CO
  REAL(KIND=jprb), PARAMETER :: comax  = 10.0_jprb        ! ppmv
  REAL(KIND=jprb), PARAMETER :: comin  = 0.1E-10_jprb     ! ppmv
  ! N2O
  REAL(KIND=jprb), PARAMETER :: n2omax = 10.0_jprb        ! ppmv
  REAL(KIND=jprb), PARAMETER :: n2omin = 0.1E-10_jprb     ! ppmv
  ! CH4
  REAL(KIND=jprb), PARAMETER :: ch4max = 50.0_jprb        ! ppmv
  REAL(KIND=jprb), PARAMETER :: ch4min = 0.1E-10_jprb     ! ppmv
  ! Cloud Liquid Water
  REAL(KIND=jprb), PARAMETER :: clwmax = 1.0_jprb         ! kg/kg
  REAL(KIND=jprb), PARAMETER :: clwmin = 0.0_jprb         ! kg/kg
  ! Surface Pressure
  REAL(KIND=jprb), PARAMETER :: pmax   = 1100.0_jprb      ! surface pressure hPa
  REAL(KIND=jprb), PARAMETER :: pmin   = 400.0_jprb       ! hPa
  ! Surface Wind
  REAL(KIND=jprb), PARAMETER :: wmax   =  100.0_jprb      ! surface wind speed (m/s)
  ! Zenith Angle
  REAL(KIND=jprb), PARAMETER :: zenmax = 75.0_jprb        ! zenith angle (Deg) = secant 3.86_jprb
  REAL(KIND=jprb), PARAMETER :: zenmaxv9 = 85.3_jprb      ! larger zenmax for v9 predictors = secant 12
  ! Cloud Top Pressure
  REAL(KIND=jprb), PARAMETER :: ctpmax = 1100.0_jprb      ! (hPa)
  REAL(KIND=jprb), PARAMETER :: ctpmin =   50.0_jprb      ! (hPa)
  ! Magnetic field strength
  REAL(KIND=jprb), PARAMETER :: bemax = 0.7_jprb          ! (Gauss)
  REAL(KIND=jprb), PARAMETER :: bemin = 0.2_jprb          ! (Guass)
  ! Ice Crystal Diameter
  REAL(KIND=jprb), PARAMETER :: dgmin_hex =  12.2_jprb    ! (micron)
  REAL(KIND=jprb), PARAMETER :: dgmax_hex =  118.29_jprb  ! (micron)
  REAL(KIND=jprb), PARAMETER :: dgmin_agg =  5.61_jprb    ! (micron)
  REAL(KIND=jprb), PARAMETER :: dgmax_agg =  166.46_jprb  ! (micron)
  ! Ice Water Content
  REAL(KIND=jprb), PARAMETER :: iwcmin_hex =  0.000608_jprb ! (g.m-3)
  REAL(KIND=jprb), PARAMETER :: iwcmax_hex =  0.254639_jprb ! (g.m-3)
  REAL(KIND=jprb), PARAMETER :: iwcmin_agg =  0.000235_jprb ! (g.m-3)
  REAL(KIND=jprb), PARAMETER :: iwcmax_agg =  0.489046_jprb ! (g.m-3)


  ! Min/max optical depth and transmittance values
  ! ----------------------------------------------
  ! maximum value of optical depth for transmittance calculation
  ! e(-30) -> 10**-14
  ! e(-50) -> 10**-22
  REAL(KIND=jprb), PARAMETER :: max_optical_depth = 50._jprb
  REAL(KIND=jprb), PARAMETER :: min_tau = 1.0e-8_jprb
  REAL(KIND=jprb), PARAMETER :: min_od  = 1.0e-5_jprb


  ! Maximum solar zenith angle for which to apply solar calculation
  ! ---------------------------------------------------------------
  REAL(KIND=jprb), PARAMETER  :: max_sol_zen = 85.3_jprb ! = secant 12


  ! Auxiliary parameters
  ! --------------------
  INTEGER(KIND=jpim), PARAMETER :: max_fastem_version = 6 ! Highest FASTEM version number available
  INTEGER(KIND=jpim), PARAMETER :: fastem_sp = 5          ! Number of fastem surface parameters
  REAL(KIND=jprb), PARAMETER    :: mwcldtp = 322.0_jprb   ! Upper pressure level (HPa) for lwp calcs
  REAL(KIND=jprb), PARAMETER    :: dcoeff(8) = &          ! Debye coefs
         (/ 17.1252_jprb, 134.2450_jprb, 310.2125_jprb,  5.667_jprb, &
           188.7979_jprb,  80.5419_jprb,   0.1157_jprb,  4.8417_jprb/)

  ! Pressure [hPa] top of radiatively insignificant "low cloud" in simple IR cloud streams approach
  REAL(KIND=jprb), PARAMETER :: cldstr_low_cloud_top = 750.0_jprb


  ! Polarisation definitions
  ! ------------------------
  ! == pol_id +1
  !   1 average of vertical and horizontal
  !   2 nominal vertical at nadir, rotating
  !      with view angle
  !   3 nominal horizontal at nadir, rotating
  !      with view angle
  !   4 vertical
  !   5 horizontal
  !   6 + 45 minus -45 (3rd stokes vector)
  !   7 left circular - right circular (4th stokes vector)
  INTEGER(KIND=jpim), PARAMETER :: npolar_compute(7) = &
   & (/ 2, 2, 2, 1, 1, 2, 4/)
  INTEGER(KIND=jpim), PARAMETER :: npolar_return(7) = &
   & (/ 1, 1, 1, 1, 1, 2, 4/)

  ! pol_v and pol_h give proportion of v and h pol to use in emissivity calculation
  ! pol_s3 adds the 3rd/4th stokes vectors
  REAL(KIND=jprb), PARAMETER :: pol_v(3,7) = RESHAPE( &
    & (/ 0.5_jprb, 0.0_jprb, 0.0_jprb, &
       & 0.0_jprb, 0.0_jprb, 1.0_jprb, &
       & 0.0_jprb, 1.0_jprb, 0.0_jprb, &
       & 1.0_jprb, 0.0_jprb, 0.0_jprb, &
       & 0.0_jprb, 0.0_jprb, 0.0_jprb, &
       & 0.0_jprb, 0.0_jprb, 0.0_jprb, &
       & 0.0_jprb, 0.0_jprb, 0.0_jprb  /), (/3,7/) )
  REAL(KIND=jprb), PARAMETER :: pol_h(3,7) = RESHAPE( &
    & (/ 0.5_jprb, 0.0_jprb, 0.0_jprb, &
       & 0.0_jprb, 1.0_jprb, 0.0_jprb, &
       & 0.0_jprb, 0.0_jprb, 1.0_jprb, &
       & 0.0_jprb, 0.0_jprb, 0.0_jprb, &
       & 1.0_jprb, 0.0_jprb, 0.0_jprb, &
       & 0.0_jprb, 0.0_jprb, 0.0_jprb, &
       & 0.0_jprb, 0.0_jprb, 0.0_jprb  /), (/3,7/) )
  REAL(KIND=jprb), PARAMETER :: pol_s3(0:1,7) = RESHAPE( &
    & (/ 0.0_jprb, 0.0_jprb, &
       & 0.0_jprb, 0.0_jprb, &
       & 0.0_jprb, 0.0_jprb, &
       & 0.0_jprb, 0.0_jprb, &
       & 0.0_jprb, 0.0_jprb, &
       & 1.0_jprb, 0.0_jprb, &
       & 0.0_jprb, 1.0_jprb  /), (/2,7/) )

  ! RTTOVSCATT aux parameters
  ! -------------------------
  ! Pressure of top level for line by line calculations (hPa)
  REAL(KIND=jprb), PARAMETER :: pressure_top = 0.004985_jprb
  ! Minimum cloud cover processed by rttov_scatt
  REAL(KIND=jprb), PARAMETER :: ccthres = 0.05_jprb
  ! Minimum single scattering albedo processed by rttov_scatt
  REAL(KIND=jprb), PARAMETER :: min_ssa = 1.0E-03_jprb
  ! Rain density (g.cm-3)
  REAL(KIND=jprb), PARAMETER :: rho_rain = 1.0_jprb
  ! Snow density (g.cm-3)
  REAL(KIND=jprb), PARAMETER :: rho_snow = 0.1_jprb

  ! Flags to identify function in shared K/Adjoint routines
  INTEGER(KIND=jpim), PARAMETER :: adk_adjoint = 0
  INTEGER(KIND=jpim), PARAMETER :: adk_k       = 1

  ! Parameters to compute refractive index of air
  ! ---------------------------------------------
  REAL(KIND=jprb), PARAMETER :: D1   =8341.87_jprb
  REAL(KIND=jprb), PARAMETER :: D2   =2405955.0_jprb
  REAL(KIND=jprb), PARAMETER :: D3   =130.0_jprb
  REAL(KIND=jprb), PARAMETER :: D4   =15996.0_jprb
  REAL(KIND=jprb), PARAMETER :: D5   =38.9_jprb
  REAL(KIND=jprb), PARAMETER :: DCO2 =0.540_jprb
  REAL(KIND=jprb), PARAMETER :: ED1  =96095.43_jprb
  REAL(KIND=jprb), PARAMETER :: ED2  =0.601_jprb
  REAL(KIND=jprb), PARAMETER :: ED3  =0.00972_jprb
  REAL(KIND=jprb), PARAMETER :: ED4  =0.003661_jprb
  REAL(KIND=jprb), PARAMETER :: EW1  =3.7345_jprb
  REAL(KIND=jprb), PARAMETER :: EW2  =0.0401_jprb
  REAL(KIND=jprb), PARAMETER :: HTOP =100.0_jprb
  REAL(KIND=jprb), PARAMETER :: CTOM =1.0E-4_jprb
  REAL(KIND=jprb), PARAMETER :: WAVER=1700.0_jprb

  ! CO2 concentration assumed for atmospheric refractivity calculation
  REAL(KIND=jprb), PARAMETER :: co2_conc = 376._jprb


  ! IR scattering parameters
  ! ------------------------

  ! Aerosols
  INTEGER(KIND=jpim), PARAMETER :: naer_max = 13

  INTEGER(KIND=jpim), PARAMETER :: &
        & aer_id_inso       = 1, &
        & aer_id_waso       = 2, &
        & aer_id_soot       = 3, &
        & aer_id_ssam       = 4, &
        & aer_id_sscm       = 5, &
        & aer_id_minm       = 6, &
        & aer_id_miam       = 7, &
        & aer_id_micm       = 8, &
        & aer_id_mitr       = 9, &
        & aer_id_suso       =10, &
        & aer_id_vola       =11, &
        & aer_id_vapo       =12, &
        & aer_id_asdu       =13

  CHARACTER(LEN=4), PARAMETER :: aer_name(naer_max) = &
        & (/ 'inso', &
           & 'waso', &
           & 'soot', &
           & 'ssam', &
           & 'sscm', &
           & 'minm', &
           & 'miam', &
           & 'micm', &
           & 'mitr', &
           & 'suso', &
           & 'vola', &
           & 'vapo', &
           & 'asdu'  /)

  ! Constants for relative humidity calculation
  REAL(KIND=jprb), PARAMETER :: E00 = 611.21_jprb
  REAL(KIND=jprb), PARAMETER :: T00 = 273.16_jprb
  REAL(KIND=jprb), PARAMETER :: TI  = T00 - 23.0_jprb

  ! Phase functions
  INTEGER(KIND=jpim), PARAMETER :: nphangle = 208

  REAL(KIND=jprb), PARAMETER :: phangle(nphangle) = &
             (/   0.0_jprb,   0.1_jprb,   0.2_jprb,   0.3_jprb,   0.4_jprb,   0.5_jprb,   0.6_jprb, &
              &   0.7_jprb,   0.8_jprb,   0.9_jprb,   1.0_jprb,   1.1_jprb,   1.2_jprb,   1.3_jprb, &
              &   1.4_jprb,   1.5_jprb,   1.6_jprb,   1.7_jprb,   1.8_jprb,   1.9_jprb,   2.0_jprb, &
              &   2.1_jprb,   2.2_jprb,   2.3_jprb,   2.4_jprb,   2.5_jprb,   2.6_jprb,   2.7_jprb, &
              &   2.8_jprb,   2.9_jprb,   3.0_jprb,   4.0_jprb,   5.0_jprb,   6.0_jprb,   7.0_jprb, &
              &   8.0_jprb,   9.0_jprb,  10.0_jprb,  11.0_jprb,  12.0_jprb,  13.0_jprb,  14.0_jprb, &
              &  15.0_jprb,  16.0_jprb,  17.0_jprb,  18.0_jprb,  19.0_jprb,  20.0_jprb,  21.0_jprb, &
              &  22.0_jprb,  23.0_jprb,  24.0_jprb,  25.0_jprb,  26.0_jprb,  27.0_jprb,  28.0_jprb, &
              &  29.0_jprb,  30.0_jprb,  31.0_jprb,  32.0_jprb,  33.0_jprb,  34.0_jprb,  35.0_jprb, &
              &  36.0_jprb,  37.0_jprb,  38.0_jprb,  39.0_jprb,  40.0_jprb,  41.0_jprb,  42.0_jprb, &
              &  43.0_jprb,  44.0_jprb,  45.0_jprb,  46.0_jprb,  47.0_jprb,  48.0_jprb,  49.0_jprb, &
              &  50.0_jprb,  51.0_jprb,  52.0_jprb,  53.0_jprb,  54.0_jprb,  55.0_jprb,  56.0_jprb, &
              &  57.0_jprb,  58.0_jprb,  59.0_jprb,  60.0_jprb,  61.0_jprb,  62.0_jprb,  63.0_jprb, &
              &  64.0_jprb,  65.0_jprb,  66.0_jprb,  67.0_jprb,  68.0_jprb,  69.0_jprb,  70.0_jprb, &
              &  71.0_jprb,  72.0_jprb,  73.0_jprb,  74.0_jprb,  75.0_jprb,  76.0_jprb,  77.0_jprb, &
              &  78.0_jprb,  79.0_jprb,  80.0_jprb,  81.0_jprb,  82.0_jprb,  83.0_jprb,  84.0_jprb, &
              &  85.0_jprb,  86.0_jprb,  87.0_jprb,  88.0_jprb,  89.0_jprb,  90.0_jprb,  91.0_jprb, &
              &  92.0_jprb,  93.0_jprb,  94.0_jprb,  95.0_jprb,  96.0_jprb,  97.0_jprb,  98.0_jprb, &
              &  99.0_jprb, 100.0_jprb, 101.0_jprb, 102.0_jprb, 103.0_jprb, 104.0_jprb, 105.0_jprb, &
              & 106.0_jprb, 107.0_jprb, 108.0_jprb, 109.0_jprb, 110.0_jprb, 111.0_jprb, 112.0_jprb, &
              & 113.0_jprb, 114.0_jprb, 115.0_jprb, 116.0_jprb, 117.0_jprb, 118.0_jprb, 119.0_jprb, &
              & 120.0_jprb, 121.0_jprb, 122.0_jprb, 123.0_jprb, 124.0_jprb, 125.0_jprb, 126.0_jprb, &
              & 127.0_jprb, 128.0_jprb, 129.0_jprb, 130.0_jprb, 131.0_jprb, 132.0_jprb, 133.0_jprb, &
              & 134.0_jprb, 135.0_jprb, 136.0_jprb, 137.0_jprb, 138.0_jprb, 139.0_jprb, 140.0_jprb, &
              & 141.0_jprb, 142.0_jprb, 143.0_jprb, 144.0_jprb, 145.0_jprb, 146.0_jprb, 147.0_jprb, &
              & 148.0_jprb, 149.0_jprb, 150.0_jprb, 151.0_jprb, 152.0_jprb, 153.0_jprb, 154.0_jprb, &
              & 155.0_jprb, 156.0_jprb, 157.0_jprb, 158.0_jprb, 159.0_jprb, 160.0_jprb, 161.0_jprb, &
              & 162.0_jprb, 163.0_jprb, 164.0_jprb, 165.0_jprb, 166.0_jprb, 167.0_jprb, 168.0_jprb, &
              & 169.0_jprb, 170.0_jprb, 171.0_jprb, 172.0_jprb, 173.0_jprb, 174.0_jprb, 175.0_jprb, &
              & 176.0_jprb, 177.0_jprb, 178.0_jprb, 179.0_jprb, 180.0_jprb /)

  ! Water clouds
  INTEGER(KIND=jpim), PARAMETER :: nwcl_max = 5

  INTEGER(KIND=jpim), PARAMETER :: &
        & wcl_id_stco       = 1, &
        & wcl_id_stma       = 2, &
        & wcl_id_cucc       = 3, &
        & wcl_id_cucp       = 4, &
        & wcl_id_cuma       = 5

  CHARACTER(LEN=4), PARAMETER :: wcl_name(nwcl_max) = &
        & (/ 'stco', &
           & 'stma', &
           & 'cucc', &
           & 'cucp', &
           & 'cuma' /)

  INTEGER(KIND=jpim), PARAMETER:: ncldtyp = 6

  ! Ice clouds
  INTEGER(KIND=jpim), PARAMETER :: nish = 4 ! Max valid value for ish
  INTEGER(KIND=jpim), PARAMETER :: nidg = 4 ! Max valid value for idg


  ! Wavenumber intervals for RTTOV9 predictors
  ! -------------------------------------------
  REAL(KIND=jprb), PARAMETER :: rttov9_wv0690_50 =  690.50_jprb, &
                                rttov9_wv1050_00 = 1050.00_jprb, &
                                rttov9_wv1095_25 = 1095.25_jprb, &
                                rttov9_wv1100_25 = 1100.25_jprb, &
                                rttov9_wv1350_25 = 1350.25_jprb, &
                                rttov9_wv1750_25 = 1750.25_jprb, &
                                rttov9_wv1900_25 = 1900.25_jprb, &
                                rttov9_wv1995_00 = 1995.00_jprb, &
                                rttov9_wv2000_00 = 2000.00_jprb, &
                                rttov9_wv2250_00 = 2250.00_jprb, &
                                rttov9_wv2295_25 = 2295.25_jprb, &
                                rttov9_wv2360_00 = 2360.00_jprb, &
                                rttov9_wv2380_25 = 2380.25_jprb, &
                                rttov9_wv2660_25 = 2660.25_jprb, &
                                rttov9_wv2760_25 = 2760.25_jprb


  ! Parameters for solar overcast radiance calculation
  ! --------------------------------------------------
  REAL(KIND=jprb), PARAMETER :: overcast_albedo_wvn = 10000._jprb ! Wavenumber (cm-1) at which albedo changes
  REAL(KIND=jprb), PARAMETER :: overcast_albedo1    = 0.7_jprb    ! Overcast albedo for wvn > limit
  REAL(KIND=jprb), PARAMETER :: overcast_albedo2    = 0.6_jprb    ! Overcast albedo for wvn <= limit

  ! Parameters for Rayleigh cross-section parameterization taken from Bucholzt 1995
  ! -------------------------------------------------------------------------------
  REAL(KIND=jprb), PARAMETER :: ray_min_wvn = 5000.0_jprb,     & ! Min wavenumber (cm-1) for which Rayleigh is calculated
                                ray_scs_wlm = 0.5_jprb,        & ! Wavelength limit: below 0.5um the
                                ray_scs_a1 = 3.01577E-28_jprb, & !   first set of parameters a1-d1 are used
                                ray_scs_b1 = -3.55212_jprb,    & !   while above this a2-d2 are used.
                                ray_scs_c1 = -1.35579_jprb,    &
                                ray_scs_d1 = -0.11563_jprb,    &
                                ray_scs_a2 = 4.01061E-28_jprb, &
                                ray_scs_b2 = -3.99668_jprb,    &
                                ray_scs_c2 = -1.10298E-3_jprb, &
                                ray_scs_d2 = -2.71393E-2_jprb

  ! Interpolation modes
  ! -------------------
  !       MODE                     USER->COEF LEVEL          COEF->USER LEVEL
  !                            (profile interpolation)   (optical depth/weighting fn interpolation)
  ! interp_rochon:                    Rochon                 Rochon,     op dep
  ! interp_loglinear:                 Log-linear             Log-linear, op dep
  ! interp_rochon_loglinear:          Rochon                 Log-linear, op dep
  ! interp_rochon_wfn:                Rochon                 Rochon,     weighting fns
  ! interp_rochon_loglinear_wfn:      Rochon                 Log-linear, weighting fns
  INTEGER(KIND=jpim), PARAMETER :: ninterp_modes               = 5_jpim ! Number of valid interpolation options
  INTEGER(KIND=jpim), PARAMETER :: interp_rochon               = 1_jpim
  INTEGER(KIND=jpim), PARAMETER :: interp_loglinear            = 2_jpim
  INTEGER(KIND=jpim), PARAMETER :: interp_rochon_loglinear     = 3_jpim
  INTEGER(KIND=jpim), PARAMETER :: interp_rochon_wfn           = 4_jpim
  INTEGER(KIND=jpim), PARAMETER :: interp_rochon_loglinear_wfn = 5_jpim

END MODULE rttov_const
