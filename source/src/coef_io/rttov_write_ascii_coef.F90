SUBROUTINE rttov_write_ascii_coef(err, coef, file_id, verbose)
! Description:
! write on unit file_id the coef structure.
! If lbinary is false or not present the file is assumed as
! an ASCII sequential formatted, in other case it is sequential unformatted.
! I/O write status are only tested at the end of the code
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
! Version   Date     Comment
! -------   ----     -------
!  1.0       01/12/2002  New F90 code with structures (P Brunel A Smith)
!  1.1       24/01/2003  insert I/O status (P Brunel)
!                        one record per channel for coefficients in binary format
!                        New header to allow checking R4<->R8
!  1.2       02/06/2004  Update for RTTOV8 coefs (P. Brunel)
!  1.3       02/08/2006  Change format for pressure levels f8.3 -> f9.4 (P. Brunel)
!  1.4       14/05/2007  Updated for RTTOV-9 (P Brunel)
!  1.5       19/02/2008  Another update for RTTOV-9 for inc_top (R Saunders)
!  1.6       27/06/2008  Introduced the case where no channels are available for
!                        the phase function in the solar range (M. Matricardi)
!  1.7       06/03/2009  Conditionals depending on coef % id_comp_lvl == 9
!                        extended to >= 9 (P. Rayer)
!  1.8       02/12/2009  Introduced principal component capability. Marco Matricardi. ECMWF
!  1.9       10/01/2013  Add PMC shifts (P Rayer)
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
  USE rttov_types, ONLY : rttov_coef
  USE parkind1, ONLY : jpim, jplm
!INTF_OFF
  USE rttov_const, ONLY :    &
       & version,            &
       & release,            &
       & minor_version,      &
       & gas_id_mixed,       &
       & gas_id_watervapour, &
       & gas_id_ozone,       &
       & gas_id_wvcont,      &
       & gas_id_co2,         &
       & gas_id_n2o,         &
       & gas_id_co,          &
       & gas_id_ch4,         &
       & gas_name,           &
       & gas_unit_name,      &
       & lensection,         &
       & speedl
!INTF_ON
  IMPLICIT NONE
! subroutine arguments
! scalar arguments with intent(in):
  INTEGER(KIND=jpim), INTENT(IN)            :: file_id! file logical unit number
  TYPE(rttov_coef)  , INTENT(IN)            :: coef   ! coefficients
  LOGICAL(KIND=jplm), INTENT(IN), OPTIONAL  :: verbose! print out information messages or not
! scalar arguments with intent(in):
  INTEGER(KIND=jpim), INTENT(OUT)           :: err    ! return code
!INTF_END
#include "rttov_errorreport.interface"
! local scalars
  INTEGER(KIND=jpim)  :: i, j, l, k
  INTEGER(KIND=jpim)  :: isat, isol, ichan
  LOGICAL(KIND=jplm)  :: lverbose
  CHARACTER(LEN = 2 ) :: sensor
  CHARACTER(LEN = lensection) :: section
  CHARACTER(LEN = 80) :: errMessage
  CHARACTER(LEN = 80) :: version_name
  CHARACTER(LEN = 20) :: FMT_sol, FMT_sat, FMT_coef, FMT_pnom
  INTEGER(KIND=jpim)  :: IncZeeman
!- End of header --------------------------------------------------------
  TRY
  IF (PRESENT(verbose)) THEN
    lverbose = verbose
  ELSE
    lverbose = .TRUE._jplm
  END IF
  
! Consider lbinary option to create the option
!ASCII file
  IF (lverbose) THEN
    WRITE (errMessage, '( "write coefficient to file_id ", i2, " in ASCII format")')file_id
    INFO(errMessage)
  END IF
  WRITE (file_id, '(a)', iostat=err)' ! RTTOV coefficient file '//Trim(coef%id_Common_name)

  THROW(err.ne.0)
 
  WRITE (file_id, '(a)', iostat=err)' ! Automatic creation by subroutine rttov_writecoef'
  
  THROW(err.ne.0)
  
  IF (release < 10 .and. minor_version < 10 ) THEN
    WRITE (version_name, '(I2.2,".",i1,".",i1)', iostat=err) version, release, minor_version
  ELSE
    WRITE (version_name, '(I2.2,".",i2.2,".",i2.2)', iostat=err) version, release, minor_version
  ENDIF

  WRITE (file_id, '(a)', iostat=err)' ! RTTOV library version '//TRIM(version_name)

  THROW(err.ne.0)

  WRITE (file_id, '(a)', iostat=err)' ! ------------------------------------------------------'

  THROW(err.ne.0)

! COEF structure (V10)
! IDENTIFICATION
  section = 'IDENTIFICATION'

  SELECT CASE (coef%id_sensor)
  CASE (1_jpim)
    sensor = 'ir'
  CASE (2_jpim)
    sensor = 'mw'
  CASE (3_jpim)
    sensor = 'hi'
  CASE (4_jpim)
    sensor = 'po'
  END SELECT


  WRITE (file_id, '(a)', iostat=err)Trim(section)

  THROW(err.ne.0)

  WRITE (file_id, '(a)', iostat=err)' ! '

  THROW(err.ne.0)

  WRITE (file_id, '(3i3,T20,a)', iostat=err)coef%id_platform, coef%id_sat, coef%id_inst, '! Platform  sat_id  instrument'

  THROW(err.ne.0)

  WRITE (file_id, '(1x,a)', iostat=err)TRIM(coef%id_Common_name)

  THROW(err.ne.0)

  WRITE (file_id, '(1x,a,T20,a)', iostat=err)sensor, '! Sensor type [ir,mw,hi,po]'

  THROW(err.ne.0)

  WRITE (file_id, '(1x,i2,T20,a)', iostat=err)coef%id_comp_lvl, '! RTTOV coefficient file version number'

  THROW(err.ne.0)

  WRITE (file_id, '(1x,a)', iostat=err)TRIM(coef%id_creation)

  THROW(err.ne.0)

  WRITE (file_id, '(1x,i4,1x,i2.2,1x,i2.2,t20,a)', iostat=err)coef%id_creation_date, '! Creation date'

  THROW(err.ne.0)

! No LINE-BY-LINE section

  IF (coef%line_by_line(1) .NE. 'xxxx') THEN
    section = 'LINE-BY-LINE'
    WRITE (file_id, '(a)', iostat=err)' ! ------------------------------------------------------'

    THROW(err.ne.0)

    WRITE (file_id, '(a)', iostat=err)Trim(section)

    THROW(err.ne.0)

    WRITE (file_id, '(a)', iostat=err)' ! Line-by-line and other information'

    THROW(err.ne.0)

    WRITE (file_id, '(a)', iostat=err)' ! '

    THROW(err.ne.0)

    DO i = 1, Size(coef%line_by_line)
      IF (coef%line_by_line(i) .EQ. 'xxxx') EXIT
      WRITE (file_id, '(a)', iostat=err)Trim(coef%line_by_line(i))

      THROW(err.ne.0)

    ENDDO

  ENDIF

! FAST_MODEL_VARIABLES
  section = 'FAST_MODEL_VARIABLES'
  WRITE (file_id, '(a)', iostat=err)' ! ------------------------------------------------------'

  THROW(err.ne.0)

  WRITE (file_id, '(a)', iostat=err)Trim(section)

  THROW(err.ne.0)

  WRITE (file_id, '(a)', iostat=err)' ! '

  THROW(err.ne.0)

  WRITE (file_id, '(a)', iostat=err)' !'

  THROW(err.ne.0)

  WRITE (file_id, '(1x,a,t20,a)', iostat=err)coef%fmv_model_def, '! Fast model name'

  THROW(err.ne.0)

  WRITE (file_id, '(1x,i4,t20,a)', iostat=err)coef%fmv_model_ver, '! Fast model version compatibility level'

  THROW(err.ne.0)

  WRITE (file_id, '(1x,i6,t20,a)', iostat=err)coef%fmv_chn, '! Number of channels described in the coef file'

  THROW(err.ne.0)

  WRITE (file_id, '(1x,i4,t20,a)', iostat=err)coef%fmv_gas, '! Number of gases described in the coef file'

  THROW(err.ne.0)

  WRITE (file_id, '(1x,i4,t20,a)' , iostat=ERR)coef%id_comp_pc, '! PC compatibility level'

  THROW(err.ne.0)

  IncZeeman = 0_jpim
  IF (coef%IncZeeman) IncZeeman = 1_jpim
  WRITE (file_id, '(1x,i4,t20,a)', iostat=err)IncZeeman, '! Zeeman flag'

  THROW(err.ne.0)


  DO i = 1, coef%fmv_gas
    WRITE (file_id, '(1x,a,t20,a)', iostat=err)Trim(gas_name(coef%fmv_gas_id(i))), '! Gas identification'

    THROW(err.ne.0)

    WRITE (file_id, '(1x,3i4,t20,a)', iostat=err)coef%fmv_var(i), coef%fmv_coe(i), coef%fmv_lvl(i),&
            & '! Variables/predictors  levels (pressure/absorber)'

    THROW(err.ne.0)

  ENDDO


  IF (coef%readme_srf(1) .NE. 'xxxx') THEN
    section = 'README_SPECTRAL_RESPONSE_FUNCTION'
    WRITE (file_id, '(a)', iostat=err)' ! ------------------------------------------------------'

    THROW(err.ne.0)

    WRITE (file_id, '(a)', iostat=err)Trim(section)

    THROW(err.ne.0)

    WRITE (file_id, '(a)', iostat=err)' ! '

    THROW(err.ne.0)

    DO i = 1, Size(coef%readme_srf)
      IF (coef%readme_srf(i) .EQ. 'xxxx') EXIT
      WRITE (file_id, '(a)', iostat=err)Trim(coef%readme_srf(i))

      THROW(err.ne.0)

    ENDDO

  ENDIF


  section = 'FILTER_FUNCTIONS'
  WRITE (file_id, '(a)', iostat=err)' ! ------------------------------------------------------'

  THROW(err.ne.0)

  WRITE (file_id, '(a)', iostat=err)Trim(section)

  THROW(err.ne.0)

  WRITE (file_id, '(a)', iostat=err)' ! '

  THROW(err.ne.0)

  WRITE (file_id, '(a)', iostat=err)' ! Channel number (from instrument original description)'

  THROW(err.ne.0)

  WRITE (file_id, '(a)', iostat=err)' ! Channel status'

  THROW(err.ne.0)

  WRITE (file_id, '(a)', iostat=err)' ! Central wavenumber'

  THROW(err.ne.0)

  WRITE (file_id, '(a)', iostat=err)' ! Band correction coefficients (offset, slope)'

  THROW(err.ne.0)

  WRITE (file_id, '(a)', iostat=err)' ! Gamma correction factor'

  THROW(err.ne.0)


  DO i = 1, coef%fmv_chn
    WRITE (file_id, '(1x,i5,1x,i4,4(1x,e18.10))', iostat=err)     &
      & coef%ff_ori_chn(i), coef%ff_val_chn(i), coef%ff_cwn(i), coef%ff_bco(i), coef%ff_bcs(i), coef%ff_gam(i)

    THROW(err.ne.0)

  ENDDO

  section = 'FUNDAMENTAL_CONSTANTS'
  WRITE (file_id, '(a)', iostat=err)' ! ------------------------------------------------------'

  THROW(err.ne.0)

  WRITE (file_id, '(a)', iostat=err)Trim(section)

  THROW(err.ne.0)

  WRITE (file_id, '(a)', iostat=err)' ! '

  THROW(err.ne.0)

  WRITE (file_id, '(a)', iostat=err)' ! Units of constants for spectral radiance'

  THROW(err.ne.0)

  WRITE (file_id, '(a)', iostat=err)' ! First radiation constant (mW/(m2.sr.cm-4))'

  THROW(err.ne.0)

  WRITE (file_id, '(a)', iostat=err)' ! Second radiation constant (cm.K)'

  THROW(err.ne.0)

  WRITE (file_id, '(1x,f14.1,t30,a)', iostat=err)speedl, '! Speed of light (cm/s)'

  THROW(err.ne.0)

  WRITE (file_id, '(1x,1p,e15.8,0p,f10.6,t30,a)', iostat=err)coef%fc_planck_c1, coef%fc_planck_c2, '! Planck constants'

  THROW(err.ne.0)

  WRITE (file_id, '(1x,f10.1,t30,a)', iostat=err)coef%fc_sat_height, '! Nominal satellite height (km)'

  THROW(err.ne.0)

  IF (Associated(coef%ss_chn)) THEN!pb
    IF (ANY(coef%ss_val_chn > 0)) THEN
      section = 'SOLAR_SPECTRUM'
      WRITE (file_id, '(a)', iostat=err)' ! ------------------------------------------------------'

      THROW(err.ne.0)

      WRITE (file_id, '(a)', iostat=err)Trim(section)

      THROW(err.ne.0)

      WRITE (file_id, '(a)', iostat=err)' ! '

      THROW(err.ne.0)

      WRITE (file_id, '(a)', iostat=err)' ! Channel number'

      THROW(err.ne.0)

      WRITE (file_id, '(a)', iostat=err)' ! Type of channel (0 => thermal; 1 => thermal+solar; 2 => solar)'

      THROW(err.ne.0)

      WRITE (file_id, '(a)', iostat=err)' ! Channel central wavenumber'

      THROW(err.ne.0)

      WRITE (file_id, '(a)', iostat=err)' ! Solar spectrum'

      THROW(err.ne.0)


      DO i = 1, coef%fmv_chn
        WRITE (file_id, '(2i5,2e19.10)', iostat=err)     &
          & coef%ss_chn(i), coef%ss_val_chn(i), coef%ss_cwn(i), coef%ss_solar_spectrum(i)

        THROW(err.ne.0)

      ENDDO
    ENDIF
  ENDIF
!pb

  IF (Any(coef%pw_val_chn > 0)) THEN
    section = 'PLANCK_WEIGHTED'
    WRITE (file_id, '(a)', iostat=err)' ! ------------------------------------------------------'

    THROW(err.ne.0)

    WRITE (file_id, '(a)', iostat=err)Trim(section)

    THROW(err.ne.0)

    WRITE (file_id, '(a)', iostat=err)' ! '

    THROW(err.ne.0)

    WRITE (file_id, '(a)', iostat=err)' ! Channel number'

    THROW(err.ne.0)

    WRITE (file_id, '(a)', iostat=err)' ! Planck-weighted flag (1 => yes; 0 => no)'

    THROW(err.ne.0)

    DO i = 1, coef%fmv_chn
      WRITE (file_id, '(2i5)', iostat=err)     &
        & coef%pw_chn(i), coef%pw_val_chn(i)

      THROW(err.ne.0)

    ENDDO

  ENDIF


! GAZ_UNITS
  section = 'GAZ_UNITS'
  WRITE (file_id, '(a)', iostat=err)' ! ------------------------------------------------------'

  THROW(err.ne.0)

  WRITE (file_id, '(a)', iostat=err)Trim(section)

  THROW(err.ne.0)

  WRITE (file_id, '(a)', iostat=err)' ! Gas concentrations can be expressed in '

  THROW(err.ne.0)

  WRITE (file_id, '(a)', iostat=err)' ! volume mixing ratio (ppmv)'

  THROW(err.ne.0)

  WRITE (file_id, '(a)', iostat=err)' ! specific concentration (kg/kg)'

  THROW(err.ne.0)

  WRITE (file_id, '(a)', iostat=err)' ! '

  THROW(err.ne.0)


  DO i = 1, coef%fmv_gas
    WRITE (file_id, '(a)', iostat=err)' !     '//gas_name(coef%fmv_gas_id(i))

    THROW(err.ne.0)

    WRITE (file_id, '(1x,i4,t20,"! ",a)', iostat=err)coef%gaz_units(i), gas_unit_name(coef%gaz_units(i))

    THROW(err.ne.0)

  ENDDO
  
  
  IF (coef%fastem_ver >= 1) THEN
    section = 'FASTEM'
    WRITE (file_id, '(a)', iostat=err)' ! ------------------------------------------------------'

    THROW(err.ne.0)

    WRITE (file_id, '(a)', iostat=err)Trim(section)

    THROW(err.ne.0)

    WRITE (file_id, '(a)', iostat=err)' ! '

    THROW(err.ne.0)

    WRITE (file_id, '(a)', iostat=err)' ! S. English fast generic millimetre wave ocean emissivity model'

    THROW(err.ne.0)

    WRITE (file_id, '(a)', iostat=err)' ! Polarisation of each channel', &
              & ' !       MPOL=0: Average of vertical and horizontal polarisation ie 0.5(H+V)', &
              & ' !       MPOL=1: Nominal vertical at nadir rotating with view angle QV', &
              & ' !       MPOL=2: Nominal horizontal at nadir rotating with view angle QH', &
              & ' !       MPOL=3: Vertical V', &
              & ' !       MPOL=4: Horizontal H', &
              & ' !       MPOL=5: +45 minus -45 (3rd stokes vector) S3', &
              & ' !       MPOL=6: Left circular - right circular (4th stokes vector) S4'

    THROW(err.ne.0)

    WRITE (file_id, '(1x,i2,a)', iostat=err)coef%fastem_ver, '   ! Version number'

    THROW(err.ne.0)

    WRITE (file_id, '(20i3)', iostat=err)(coef%fastem_polar(i), i = 1, coef%fmv_chn)

    THROW(err.ne.0)

  ENDIF


  IF (coef%ssirem_ver >= 1) THEN
    section = 'SSIREM'
    WRITE (file_id, '(a)', iostat=err)' ! ------------------------------------------------------'

    THROW(err.ne.0)

    WRITE (file_id, '(a)', iostat=err)Trim(section)

    THROW(err.ne.0)

    WRITE (file_id, '(a)', iostat=err)' ! '

    THROW(err.ne.0)

    WRITE (file_id, '(a)', iostat=err)' ! Channel number'

    THROW(err.ne.0)

    WRITE (file_id, '(a)', iostat=err)' ! 5 coefficients for emissivity model SSIREM'

    THROW(err.ne.0)

    WRITE (file_id, '(1x,i2,a)', iostat=err)coef%ssirem_ver, '   ! Version number'

    THROW(err.ne.0)


    DO i = 1, coef%fmv_chn
      WRITE (file_id, '(1x,i5,3f12.7,2f4.1)', iostat=err)coef%ssirem_chn(i), coef%ssirem_a0(i), coef%ssirem_a1(i),      &
        & coef%ssirem_a2(i), coef%ssirem_xzn1(i), coef%ssirem_xzn2(i)

      THROW(err.ne.0)

    ENDDO

  ENDIF

  IF ((coef%nintmixed > 0) .AND. (coef%nintwater > 0) .AND. (coef%nintozone > 0) .AND. (coef%nintwvcont > 0) .AND. &
      (coef%nintco2 > 0) .AND. (coef%nintn2o > 0) .AND. (coef%nintco > 0) .AND. (coef%nintch4 > 0)) THEN!pb
    section = 'GAS_SPECTRAL_INTERVAL'
    WRITE (file_id, '(a)', iostat=err)' ! ------------------------------------------------------'

    THROW(err.ne.0)

    WRITE (file_id, '(a)', iostat=err)Trim(section)

    THROW(err.ne.0)

    WRITE (file_id, '(a)', iostat=err)' ! '

    THROW(err.ne.0)


    IF (coef%nintmixed > 0) THEN
      WRITE (file_id, '(1x,a,t20,a)', iostat=err)Trim(gas_name(gas_id_mixed)), '! Gas identification'

      THROW(err.ne.0)

      WRITE (file_id, '(1x,i4)', iostat=err)coef%nintmixed

      THROW(err.ne.0)


      DO i = 1, coef%nintmixed
        WRITE (file_id, '(1x,2f10.3)', iostat=err)coef%mixedgasint(1, i), coef%mixedgasint(2, i)

        THROW(err.ne.0)

      ENDDO

    ENDIF


    IF (coef%nintwater > 0) THEN
      WRITE (file_id, '(1x,a,t20,a)', iostat=err)Trim(gas_name(gas_id_watervapour)), '! Gas identification'

      THROW(err.ne.0)

      WRITE (file_id, '(1x,i4)', iostat=err)coef%nintwater

      THROW(err.ne.0)


      DO i = 1, coef%nintwater
        WRITE (file_id, '(1x,2f10.3)', iostat=err)coef%watervapourint(1, i), coef%watervapourint(2, i)

        THROW(err.ne.0)

      ENDDO

    ENDIF


    IF (coef%nintozone > 0) THEN
      WRITE (file_id, '(1x,a,t20,a)', iostat=err)Trim(gas_name(gas_id_ozone)), '! Gas identification'

      THROW(err.ne.0)

      WRITE (file_id, '(1x,i4)', iostat=err)coef%nintozone

      THROW(err.ne.0)


      DO i = 1, coef%nintozone
        WRITE (file_id, '(1x,2f10.3)', iostat=err)coef%ozoneint(1, i), coef%ozoneint(2, i)

        THROW(err.ne.0)

      ENDDO

    ENDIF


    IF (coef%nintwvcont > 0) THEN
      WRITE (file_id, '(1x,a,t20,a)', iostat=err)Trim(gas_name(gas_id_wvcont)), '! Gas identification'

      THROW(err.ne.0)

      WRITE (file_id, '(1x,i4)', iostat=err)coef%nintwvcont

      THROW(err.ne.0)


      DO i = 1, coef%nintwvcont
        WRITE (file_id, '(1x,2f10.3)', iostat=err)coef%wvcontint(1, i), coef%wvcontint(2, i)

        THROW(err.ne.0)

      ENDDO

    ENDIF


    IF (coef%nintco2 > 0) THEN
      WRITE (file_id, '(1x,a,t20,a)', iostat=err)Trim(gas_name(gas_id_co2)), '! Gas identification'

      THROW(err.ne.0)

      WRITE (file_id, '(1x,i4)', iostat=err)coef%nintco2

      THROW(err.ne.0)


      DO i = 1, coef%nintco2
        WRITE (file_id, '(1x,2f10.3)', iostat=err)coef%co2int(1, i), coef%co2int(2, i)

        THROW(err.ne.0)

      ENDDO

    ENDIF


    IF (coef%nintn2o > 0) THEN
      WRITE (file_id, '(1x,a,t20,a)', iostat=err)Trim(gas_name(gas_id_n2o)), '! Gas identification'

      THROW(err.ne.0)

      WRITE (file_id, '(1x,i4)', iostat=err)coef%nintn2o

      THROW(err.ne.0)


      DO i = 1, coef%nintn2o
        WRITE (file_id, '(1x,2f10.3)', iostat=err)coef%n2oint(1, i), coef%n2oint(2, i)

        THROW(err.ne.0)

      ENDDO

    ENDIF


    IF (coef%nintco > 0) THEN
      WRITE (file_id, '(1x,a,t20,a)', iostat=err)Trim(gas_name(gas_id_co)), '! Gas identification'

      THROW(err.ne.0)

      WRITE (file_id, '(1x,i4)', iostat=err)coef%nintco

      THROW(err.ne.0)


      DO i = 1, coef%nintco
        WRITE (file_id, '(1x,2f10.3)', iostat=err)coef%coint(1, i), coef%coint(2, i)

        THROW(err.ne.0)

      ENDDO

    ENDIF


    IF (coef%nintch4 > 0) THEN
      WRITE (file_id, '(1x,a,t20,a)', iostat=err)Trim(gas_name(gas_id_ch4)), '! Gas identification'

      THROW(err.ne.0)

      WRITE (file_id, '(1x,i4)', iostat=err)coef%nintch4

      THROW(err.ne.0)


      DO i = 1, coef%nintch4
        WRITE (file_id, '(1x,2f10.3)', iostat=err)coef%ch4int(1, i), coef%ch4int(2, i)

        THROW(err.ne.0)

      ENDDO

    ENDIF

  ENDIF


  IF (Associated(coef%tt_chn)) THEN!pb
    section = 'TRANSMITTANCE_TRESHOLD'
    WRITE (file_id, '(a)', iostat=err)' ! ------------------------------------------------------'

    THROW(err.ne.0)

    WRITE (file_id, '(a)', iostat=err)Trim(section)

    THROW(err.ne.0)

    WRITE (file_id, '(a)', iostat=err)' ! '

    THROW(err.ne.0)

    WRITE (file_id, '(a)', iostat=err)' ! Channel number'

    THROW(err.ne.0)

    WRITE (file_id, '(a)', iostat=err)' ! Apply transmittance threshold (1 => yes; 0 => no)'

    THROW(err.ne.0)

    WRITE (file_id, '(a)', iostat=err)' ! Channel central wavenumber'

    THROW(err.ne.0)

    WRITE (file_id, '(a)', iostat=err)' ! Transmittance threshold'

    THROW(err.ne.0)

    WRITE (file_id, '(a)', iostat=err)' ! Transmittance value'

    THROW(err.ne.0)


    DO i = 1, coef%fmv_chn
      WRITE (file_id, '(2i5,3e19.10)', iostat=err)     &
        & coef%tt_chn(i), coef%tt_val_chn(i), coef%tt_cwn(i), coef%tt_a0(i), coef%tt_a1(i)

      THROW(err.ne.0)

    ENDDO

  ENDIF
!pb




  IF (associated(coef%woc_chn)) THEN!pb
    section = 'WATER_OPTICAL_CONSTANT'
    WRITE (file_id, '(a)', iostat=err)' ! ------------------------------------------------------'

    THROW(err.ne.0)

    WRITE (file_id, '(a)', iostat=err)Trim(section)

    THROW(err.ne.0)

    WRITE (file_id, '(a)', iostat=err)' ! '

    THROW(err.ne.0)

    WRITE (file_id, '(a)', iostat=err)' ! Channel number'

    THROW(err.ne.0)

    WRITE (file_id, '(a)', iostat=err)' ! Channel central wavenumber'

    THROW(err.ne.0)

    WRITE (file_id, '(a)', iostat=err)' ! Ocean water (real and imaginary part)'

    THROW(err.ne.0)

    WRITE (file_id, '(a)', iostat=err)' ! Fresh water (real and imaginary part)'

    THROW(err.ne.0)


    DO i = 1, coef%fmv_chn
      WRITE (file_id, '(i5,e19.10,2(" (",e17.10,",",e17.10,")"))', iostat=err)     &
        & coef%woc_chn(i), coef%woc_cwn(i), coef%woc_waopc_ow(i), coef%woc_waopc_fw(i)

      THROW(err.ne.0)

    ENDDO

  ENDIF
!pb


  section = 'REFERENCE_PROFILE'
  WRITE (file_id, '(a)', iostat=err)' ! ------------------------------------------------------'

  THROW(err.ne.0)

  WRITE (file_id, '(a)', iostat=err)Trim(section)

  THROW(err.ne.0)

  WRITE (file_id, '(a)', iostat=err)' ! '

  THROW(err.ne.0)

  WRITE (file_id, '(a)', iostat=err)' ! Reference pressure (hPa), reference temperature (K) and'

  THROW(err.ne.0)

  WRITE (file_id, '(a)', iostat=err)' ! reference volume mixing ratio (ppmv) for each gas'

  THROW(err.ne.0)

  WRITE (file_id, '(a)', iostat=err)' ! Note that mixing ratio is "missing" for mixed gases'

  THROW(err.ne.0)


  DO i = 1, coef%fmv_gas
    WRITE (file_id, '(a)', iostat=err)' !     '//gas_name(coef%fmv_gas_id(i))

    THROW(err.ne.0)


    DO l = 1, coef%fmv_lvl(i)
      WRITE (file_id, '(1x,f9.4,2x,f7.3,1x,e13.6)')coef%ref_prfl_p(l), coef%ref_prfl_t(l, i), coef%ref_prfl_mr(l, i)
!!$                & coef % ref_prfl_p(l), coef % ref_prfl_t(l,i), ref_mr(l,i)

      THROW(err.ne.0)

    ENDDO

  ENDDO

  section = 'PROFILE_LIMITS'
  WRITE (file_id, '(a)', iostat=err)' ! ------------------------------------------------------'

  THROW(err.ne.0)

  WRITE (file_id, '(a)', iostat=err)Trim(section)

  THROW(err.ne.0)

  WRITE (file_id, '(a)', iostat=err)' ! '

  THROW(err.ne.0)

  WRITE (file_id, '(a)', iostat=err)' ! Reference pressure (hPa), temperature max and min (K) and'

  THROW(err.ne.0)

  WRITE (file_id, '(a)', iostat=err)' ! volume mixing ratio max and min (ppmv) for each gas'

  THROW(err.ne.0)

  WRITE (file_id, '(a)', iostat=err)' !      Temperature'

  THROW(err.ne.0)


  DO l = 1, coef%fmv_lvl(1)
    WRITE (file_id, '(1x,f9.4,2(1x,f7.2))', iostat=err)coef%lim_prfl_p(l), coef%lim_prfl_tmax(l), coef%lim_prfl_tmin(l)

    THROW(err.ne.0)

  ENDDO


  DO i = 1, coef%fmv_gas
    WRITE (file_id, '(a)', iostat=err)' !     '//gas_name(coef%fmv_gas_id(i))

    THROW(err.ne.0)


    DO l = 1, coef%fmv_lvl(i)
      WRITE (file_id, '(1x,f9.4,2x,e12.4,e12.4)', iostat=err)     &
        & coef%lim_prfl_p(l), coef%lim_prfl_gmax(l, i), coef%lim_prfl_gmin(l, i)

      THROW(err.ne.0)

    ENDDO

  ENDDO

  section = 'FAST_COEFFICIENTS'
  WRITE (file_id, '(a)', iostat=err)' ! ------------------------------------------------------'

  THROW(err.ne.0)

  WRITE (file_id, '(a)', iostat=err)Trim(section)

  THROW(err.ne.0)

  WRITE (file_id, '(a)', iostat=err)' ! '

  THROW(err.ne.0)

  WRITE (file_id, '(a)', iostat=err)' ! Transmission coefficients'

  THROW(err.ne.0)

  WRITE (file_id, '(a)', iostat=err)' ! Order of the gases:'

  THROW(err.ne.0)


  DO i = 1, coef%fmv_gas
    WRITE (file_id, '(a)', iostat=err)' !     '//gas_name(coef%fmv_gas_id(i))

    THROW(err.ne.0)

  ENDDO


  DO l = 1, coef%fmv_gas
    WRITE (file_id, '(a)', iostat=err)gas_name(coef%fmv_gas_id(l))

    THROW(err.ne.0)


    SELECT CASE (coef%fmv_gas_id(l))
    CASE (gas_id_mixed)
      WRITE (file_id, '(5(1x,e15.8))', iostat=err)     &
        & (((coef%thermal%mixedgas(i, j, k), i = 1, coef%fmv_lvl(l) - 1), j = 1, coef%fmv_chn), k = 1, coef%ncmixed)

      THROW(err.ne.0)

    CASE (gas_id_watervapour)
      WRITE (file_id, '(5(1x,e15.8))', iostat=err)     &
        & (((coef%thermal%watervapour(i, j, k), i = 1, coef%fmv_lvl(l) - 1), j = 1, coef%fmv_chn), k = 1, coef%ncwater)

      THROW(err.ne.0)

    CASE (gas_id_ozone)
      WRITE (file_id, '(5(1x,e15.8))', iostat=err)     &
        & (((coef%thermal%ozone(i, j, k), i = 1, coef%fmv_lvl(l) - 1), j = 1, coef%fmv_chn), k = 1, coef%ncozone)

      THROW(err.ne.0)

    CASE (gas_id_wvcont)
      WRITE (file_id, '(5(1x,e15.8))', iostat=err)     &
        & (((coef%thermal%wvcont(i, j, k), i = 1, coef%fmv_lvl(l) - 1), j = 1, coef%fmv_chn), k = 1, coef%ncwvcont)

      THROW(err.ne.0)

    CASE (gas_id_co2)
      WRITE (file_id, '(5(1x,e15.8))', iostat=err)     &
        & (((coef%thermal%co2(i, j, k), i = 1, coef%fmv_lvl(l) - 1), j = 1, coef%fmv_chn), k = 1, coef%ncco2)

      THROW(err.ne.0)

    CASE (gas_id_n2o)
      WRITE (file_id, '(5(1x,e15.8))', iostat=err)     &
        & (((coef%thermal%n2o(i, j, k), i = 1, coef%fmv_lvl(l) - 1), j = 1, coef%fmv_chn), k = 1, coef%ncn2o)

      THROW(err.ne.0)

    CASE (gas_id_co)
      WRITE (file_id, '(5(1x,e15.8))', iostat=err)     &
        & (((coef%thermal%co(i, j, k), i = 1, coef%fmv_lvl(l) - 1), j = 1, coef%fmv_chn), k = 1, coef%ncco)

      THROW(err.ne.0)

    CASE (gas_id_ch4)
      WRITE (file_id, '(5(1x,e15.8))', iostat=err)     &
        & (((coef%thermal%ch4(i, j, k), i = 1, coef%fmv_lvl(l) - 1), j = 1, coef%fmv_chn), k = 1, coef%ncch4)

      THROW(err.ne.0)

    END SELECT

  ENDDO

  IF (coef%solarcoef) THEN
    section = 'SOLAR_FAST_COEFFICIENTS'
    WRITE (file_id, '(a)', iostat=err)' ! ------------------------------------------------------'
  
    THROW(err.ne.0)
  
    WRITE (file_id, '(a)', iostat=err)Trim(section)
  
    THROW(err.ne.0)
  
    WRITE (file_id, '(a)', iostat=err)' ! '
  
    THROW(err.ne.0)
  
    WRITE (file_id, '(a)', iostat=err)' ! Transmission coefficients'
  
    THROW(err.ne.0)
  
    WRITE (file_id, '(a)', iostat=err)' ! Order of the gases:'
  
    THROW(err.ne.0)
  
  
    DO i = 1, coef%fmv_gas
      WRITE (file_id, '(a)', iostat=err)' !     '//gas_name(coef%fmv_gas_id(i))
  
      THROW(err.ne.0)
  
    ENDDO
  
  
    DO l = 1, coef%fmv_gas
      WRITE (file_id, '(a)', iostat=err)gas_name(coef%fmv_gas_id(l))
  
      THROW(err.ne.0)
  
  
      SELECT CASE (coef%fmv_gas_id(l))
      CASE (gas_id_mixed)
        WRITE (file_id, '(5(1x,e15.8))', iostat=err)     &
          & (((coef%solar%mixedgas(i, j, k), i = 1, coef%fmv_lvl(l) - 1), j = 1, coef%fmv_chn), k = 1, coef%ncmixed)
  
        THROW(err.ne.0)
  
      CASE (gas_id_watervapour)
        WRITE (file_id, '(5(1x,e15.8))', iostat=err)     &
          & (((coef%solar%watervapour(i, j, k), i = 1, coef%fmv_lvl(l) - 1), j = 1, coef%fmv_chn), k = 1, coef%ncwater)
  
        THROW(err.ne.0)
  
      CASE (gas_id_ozone)
        WRITE (file_id, '(5(1x,e15.8))', iostat=err)     &
          & (((coef%solar%ozone(i, j, k), i = 1, coef%fmv_lvl(l) - 1), j = 1, coef%fmv_chn), k = 1, coef%ncozone)
  
        THROW(err.ne.0)
  
      CASE (gas_id_wvcont)
        WRITE (file_id, '(5(1x,e15.8))', iostat=err)     &
          & (((coef%solar%wvcont(i, j, k), i = 1, coef%fmv_lvl(l) - 1), j = 1, coef%fmv_chn), k = 1, coef%ncwvcont)
  
        THROW(err.ne.0)
  
      CASE (gas_id_co2)
        WRITE (file_id, '(5(1x,e15.8))', iostat=err)     &
          & (((coef%solar%co2(i, j, k), i = 1, coef%fmv_lvl(l) - 1), j = 1, coef%fmv_chn), k = 1, coef%ncco2)
  
        THROW(err.ne.0)
  
      CASE (gas_id_n2o)
        WRITE (file_id, '(5(1x,e15.8))', iostat=err)     &
          & (((coef%solar%n2o(i, j, k), i = 1, coef%fmv_lvl(l) - 1), j = 1, coef%fmv_chn), k = 1, coef%ncn2o)
  
        THROW(err.ne.0)
  
      CASE (gas_id_co)
        WRITE (file_id, '(5(1x,e15.8))', iostat=err)     &
          & (((coef%solar%co(i, j, k), i = 1, coef%fmv_lvl(l) - 1), j = 1, coef%fmv_chn), k = 1, coef%ncco)
  
        THROW(err.ne.0)
  
      CASE (gas_id_ch4)
        WRITE (file_id, '(5(1x,e15.8))', iostat=err)     &
          & (((coef%solar%ch4(i, j, k), i = 1, coef%fmv_lvl(l) - 1), j = 1, coef%fmv_chn), k = 1, coef%ncch4)
  
        THROW(err.ne.0)
  
      END SELECT
  
    ENDDO
  ENDIF
  
  IF (ASSOCIATED(coef%ws_npoint)) THEN!pb
    section = 'WAVE_SPECTRUM'
    WRITE (file_id, '(a)', iostat=err)' ! ------------------------------------------------------'

    THROW(err.ne.0)

    WRITE (file_id, '(a)', iostat=err)Trim(section)

    THROW(err.ne.0)

    WRITE (file_id, '(a)', iostat=err)' ! '

    THROW(err.ne.0)

    WRITE (file_id, '(a)', iostat=err)' ! Number of points'

    THROW(err.ne.0)

    WRITE (file_id, '(a)', iostat=err)' ! Point number'

    THROW(err.ne.0)

    WRITE (file_id, '(a)', iostat=err)' ! Wave spectrum'

    THROW(err.ne.0)

    WRITE (file_id,  * , iostat=err)coef%ws_nomega

    THROW(err.ne.0)


    DO i = 1, coef%ws_nomega
      WRITE (file_id, '(f10.3,f12.5)', iostat=err)coef%ws_npoint(i), coef%ws_k_omega(i)

      THROW(err.ne.0)

    ENDDO

  ENDIF

  IF (coef%nltecoef) THEN
    section = 'NLTE_RADIANCE_COEFS'
    WRITE (file_id, '(a)', iostat=err) ' ! ------------------------------------------------------'  
    WRITE (file_id, '(a)', iostat=err) TRIM(section)  
    WRITE (file_id, '(a)', iostat=err) ' ! '
    WRITE (file_id, '(a)', iostat=err) ' ! Radiance coefficients'
    WRITE (file_id, '(a)', iostat=err) ' ! ncoef nsol nsat nchan start_chan'
    WRITE (file_id, '(a)', iostat=err) ' ! sol_zen_angle'
    WRITE (file_id, '(a)', iostat=err) ' ! sat_zen_angle'
    WRITE (file_id, '(a)', iostat=err) &
      ' ! nsol * nsat * nchan * [co c1 c2] coefficients'

    THROW(err .NE. 0)

    WRITE (file_id, '(5(i6))', iostat=err) coef%nlte_coef%ncoef, &
      coef%nlte_coef%nsol, coef%nlte_coef%nsat, coef%nlte_coef%nchan, &
      coef%nlte_coef%start_chan

    WRITE(FMT_sol,'(a,i1,a)') '(',coef%nlte_coef%nsol,'(f6.0))'
    IF (coef%nlte_coef%nsat < 10) THEN
      WRITE(FMT_sat,'(a,i1,a)') '(',coef%nlte_coef%nsat,'(f6.3))'
    ELSE
      WRITE(FMT_sat,'(a,i2,a)') '(',coef%nlte_coef%nsat,'(f6.3))'
    ENDIF
    WRITE(FMT_coef,'(a,i1,a)') '(',coef%nlte_coef%ncoef,'(e15.7))'

    WRITE (file_id, TRIM(FMT_sol), iostat=err) coef%nlte_coef%sol_zen_angle
    WRITE (file_id, TRIM(FMT_sat), iostat=err) coef%nlte_coef%sec_sat

    DO ichan = 1, coef%nlte_coef%nchan
      DO isol = 1, coef%nlte_coef%nsol
        DO isat = 1, coef%nlte_coef%nsat
          WRITE(file_id, TRIM(FMT_coef), iostat=ERR) &
            coef%nlte_coef%coef(:, isat, isol, ichan)
        ENDDO
      ENDDO
    ENDDO    
    
  ENDIF

  IF (coef%pmc_shift) THEN
    section = 'PRESSURE_MODULATED_CELL'
    WRITE (file_id, '(a)', iostat=err) ' ! ------------------------------------------------------'  
    WRITE (file_id, '(a)', iostat=err) TRIM(section)  
    WRITE (file_id, '(a)', iostat=err) ' ! '

    WRITE (file_id, '(a)', iostat=err) ' ! Cell length (cm)' 
    WRITE (file_id, '(a)', iostat=err) ' ! Nominal cell pressures (hPa) - as used for the above fast coefficients' 
    WRITE (file_id, '(a)', iostat=err) ' ! Temperature of cell (K) - fixed'
    WRITE (file_id, '(a)', iostat=err) ' ! Ratio gamma(co2)/gamma(air) - band-averaged co2 halfwidth ratio'
    WRITE (file_id, '(a)', iostat=err) ' ! Number of layers used  - may be a subset of coef%nlayers'
    WRITE (file_id, '(a)', iostat=err) ' ! Number of variables used'
    WRITE (file_id, '(a)', iostat=err) ' ! Coefficients - lev, ((coef(lev,chan,var), var=1,nvar), chan=1,nchan)'

    THROW(err .NE. 0)

    WRITE (file_id, *, iostat=err)  coef%pmc_lengthcell

    IF(coef%fmv_chn < 10) THEN  
      WRITE(FMT_pnom,'(a,i1,a)') '(',coef%fmv_chn,'(f6.2,1x))'
    ELSE
      WRITE(FMT_pnom,'(a,i2,a)') '(',coef%fmv_chn,'(f6.2,1x))'
    ENDIF
    WRITE (file_id, TRIM(FMT_pnom), iostat=err) (coef%pmc_pnominal(i),i=1,coef%fmv_chn)
    WRITE (file_id, *, iostat=err)  coef%pmc_tempcell
    WRITE (file_id, *, iostat=err)  coef%pmc_betaplus1
    WRITE (file_id, *, iostat=err)  coef%pmc_nlay
    WRITE (file_id, *, iostat=err)  coef%pmc_nvar
    IF(coef%fmv_chn < 10) THEN  
      WRITE(FMT_coef,'(a,i1,a)') '(',coef%pmc_nvar*coef%fmv_chn,'(e17.8))'
    ELSE
      WRITE(FMT_coef,'(a,i2,a)') '(',coef%pmc_nvar*coef%fmv_chn,'(e17.8))'
    ENDIF
    DO i=1,coef%pmc_nlay
      WRITE (file_id, TRIM(FMT_coef), iostat=ERR)     &
        & ((coef%pmc_coef(i, j, k), k = 1, coef%pmc_nvar), j = 1, coef%fmv_chn)
      ENDDO
  ENDIF

  section = 'END'
  WRITE (file_id, '(a)', iostat=err)' ! ------------------------------------------------------'

  THROW(err.ne.0)

  WRITE (file_id, '(a)', iostat=err)Trim(section)

  THROW(err.ne.0)

  IF (lverbose) INFO("end of write coefficient")
  CATCH
END SUBROUTINE 
