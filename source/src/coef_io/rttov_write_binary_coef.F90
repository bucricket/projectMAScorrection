SUBROUTINE rttov_write_binary_coef( &
            & err,           &
            & coef,          &
            & file_id,       &
            & verbose)
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
  USE rttov_types, ONLY :  &
       & rttov_coef
  USE parkind1, ONLY : jpim, jplm
!INTF_OFF
  USE rttov_const, ONLY :  &
       & rttov_magic_string, &
       & rttov_magic_number, &
       & speedl
!INTF_ON
  IMPLICIT NONE
! subroutine arguments
! scalar arguments with intent(in):
  INTEGER(KIND=jpim)       , INTENT(IN)              :: file_id      ! file logical unit number
  TYPE(rttov_coef         ), INTENT(IN)              :: coef         ! coefficients
  LOGICAL(KIND=jplm)       , INTENT(IN), OPTIONAL    :: verbose      ! print out information messages or not
! scalar arguments with intent(in):
  INTEGER(KIND=jpim)       , INTENT(OUT)             :: err          ! return code
!INTF_END
#include "rttov_errorreport.interface"
! local scalars
  LOGICAL(KIND=jplm)  :: section_present
  LOGICAL(KIND=jplm)  :: lverbose
  INTEGER(KIND=jpim)  :: i
  CHARACTER(LEN = 80) :: errMessage
  INTEGER(KIND=jpim)  :: IncZeeman
!- End of header --------------------------------------------------------
  TRY
  IF (PRESENT(verbose)) THEN
    lverbose = verbose
  ELSE
    lverbose = .TRUE._jplm
  END IF
  
! Consider lbinary option to create the option

! Binary file
  IF (lverbose) THEN
    WRITE (errMessage, '( "write coefficient to file_id ", i2, " in binary format")')file_id
    INFO(errMessage)
  END IF
! Write a string that could be displayed
! Write a real number to be able to check single/double precision
  WRITE (file_id, iostat=err)rttov_magic_string, rttov_magic_number

  THROW(err.ne.0)

! COEF structure (V10)

    WRITE (file_id, iostat=err)coef%id_platform, coef%id_sat, coef%id_inst, coef%id_sensor

    THROW(err.ne.0)


    WRITE (file_id, iostat=err)coef%id_comp_lvl, coef%id_creation_date, coef%id_creation, coef%id_Common_name

    THROW(err.ne.0)

    WRITE (file_id, iostat=err)coef%fmv_model_def, coef%fmv_model_ver, coef%fmv_chn, coef%fmv_gas

    THROW(err.ne.0)

    WRITE (file_id, iostat=ERR)coef%id_comp_pc

    THROW(err.ne.0)

    IncZeeman = 0_jpim
    IF (coef%IncZeeman) IncZeeman = 1_jpim
    WRITE (file_id, iostat=err) IncZeeman

    THROW(err.ne.0)

    WRITE (file_id, iostat=err)coef%fmv_gas_id, coef%fmv_gas_pos, coef%fmv_var, coef%fmv_lvl

    THROW(err.ne.0)


    WRITE (file_id, iostat=err)coef%fmv_coe!pb

    THROW(err.ne.0)
    
    WRITE (file_id, iostat=err)coef%gaz_units

    THROW(err.ne.0)


! GAS_SPECTRAL_INTERVAL                  !pb
    WRITE (file_id, iostat=err)coef%nintmixed, coef%nintwater, coef%nintozone, coef%nintwvcont, coef%nintco2,      &
      & coef%nintn2o, coef%nintco, coef%nintch4!pb

    THROW(err.ne.0)


    IF (coef%nintmixed > 0) THEN
      WRITE (file_id, iostat=err)coef%mixedgasint!pb

      THROW(err.ne.0)

    ENDIF


    IF (coef%nintwater > 0) THEN
      WRITE (file_id, iostat=err)coef%watervapourint!pb

      THROW(err.ne.0)

    ENDIF


    IF (coef%nintozone > 0) THEN
      WRITE (file_id, iostat=err)coef%ozoneint!pb

      THROW(err.ne.0)

    ENDIF


    IF (coef%nintwvcont > 0) THEN
      WRITE (file_id, iostat=err)coef%wvcontint!pb

      THROW(err.ne.0)

    ENDIF


    IF (coef%nintco2 > 0) THEN
      WRITE (file_id, iostat=err)coef%co2int!pb

      THROW(err.ne.0)

    ENDIF


    IF (coef%nintn2o > 0) THEN
      WRITE (file_id, iostat=err)coef%n2oint!pb

      THROW(err.ne.0)

    ENDIF


    IF (coef%nintco > 0) THEN
      WRITE (file_id, iostat=err)coef%coint!pb

      THROW(err.ne.0)

    ENDIF


    IF (coef%nintch4 > 0) THEN
      WRITE (file_id, iostat=err)coef%ch4int!pb

      THROW(err.ne.0)

    ENDIF


    section_present = Associated(coef%tt_chn)
!TRANSMITTANCE_TRESHOLD                  !pb
    WRITE (file_id, iostat=err)section_present

    THROW(err.ne.0)


    IF (section_present) THEN
      WRITE (file_id, iostat=err)coef%tt_chn, coef%tt_val_chn, coef%tt_cwn, coef%tt_a0, coef%tt_a1!pb

      THROW(err.ne.0)

    ENDIF

    

!SOLAR_SPECTRUM                          !pb
    section_present = .FALSE.
    IF (ASSOCIATED(coef%ss_chn)) section_present = ANY(coef%ss_val_chn > 0)

    WRITE (file_id, iostat=err)section_present
    THROW(err.ne.0)

    IF (section_present) THEN
      WRITE (file_id, iostat=err)coef%ss_chn, coef%ss_val_chn, coef%ss_cwn, coef%ss_solar_spectrum!pb
    ENDIF


!WATER_OPTICAL_CONSTANT                 !pb
    section_present = Associated(coef%woc_chn)
    WRITE (file_id, iostat=err)section_present

    THROW(err.ne.0)


    IF (section_present) THEN
      WRITE (file_id, iostat=err)coef%woc_chn, coef%woc_cwn, coef%woc_waopc_ow, coef%woc_waopc_fw!pb

      THROW(err.ne.0)

    ENDIF

    
!WAVE_SPECTRUM                          !pb
    section_present = Associated(coef%ws_k_omega)
    WRITE (file_id, iostat=err)section_present

    THROW(err.ne.0)


    IF (section_present) THEN
      WRITE (file_id, iostat=err)coef%ws_nomega!pb

      THROW(err.ne.0)

      WRITE (file_id, iostat=err)coef%ws_k_omega, coef%ws_npoint!pb

      THROW(err.ne.0)

    ENDIF

    
    WRITE (file_id, iostat=err)coef%ff_ori_chn, coef%ff_val_chn, coef%ff_cwn, coef%ff_bco, coef%ff_bcs, coef%ff_gam

    THROW(err.ne.0)

    WRITE (file_id, iostat=err)speedl, coef%fc_planck_c1, coef%fc_planck_c2, coef%fc_sat_height

    THROW(err.ne.0)

    WRITE (file_id, iostat=err)coef%fastem_ver, coef%ssirem_ver

    THROW(err.ne.0)


    IF (coef%fastem_ver >= 1) THEN
      WRITE (file_id, iostat=err)coef%fastem_polar

      THROW(err.ne.0)

    ENDIF


    IF (coef%ssirem_ver >= 1) THEN
      WRITE (file_id, iostat=err)     &
        & coef%ssirem_chn, coef%ssirem_a0, coef%ssirem_a1, coef%ssirem_a2, coef%ssirem_xzn1, coef%ssirem_xzn2

      THROW(err.ne.0)

    ENDIF

    WRITE (file_id, iostat=err)coef%ref_prfl_p, coef%ref_prfl_t, coef%ref_prfl_mr

    THROW(err.ne.0)

    WRITE (file_id, iostat=err)     &
      & coef%lim_prfl_p, coef%lim_prfl_tmax, coef%lim_prfl_tmin, coef%lim_prfl_gmax, coef%lim_prfl_gmin

    THROW(err.ne.0)

! Write coefficients with ONE record per Channel

    IF (coef%nmixed > 0) THEN

      DO i = 1, coef%fmv_chn
        WRITE (file_id, iostat=err)coef%thermal%mixedgas(:, i, :)

        THROW(err.ne.0)

      ENDDO

    ENDIF


    IF (coef%nwater > 0) THEN

      DO i = 1, coef%fmv_chn
        WRITE (file_id, iostat=err)coef%thermal%watervapour(:, i, :)

        THROW(err.ne.0)

      ENDDO

    ENDIF


    IF (coef%nozone > 0) THEN

      DO i = 1, coef%fmv_chn
        WRITE (file_id, iostat=err)coef%thermal%ozone(:, i, :)

        THROW(err.ne.0)

      ENDDO

    ENDIF


    IF (coef%nwvcont > 0) THEN

      DO i = 1, coef%fmv_chn
        WRITE (file_id, iostat=err)coef%thermal%wvcont(:, i, :)

        THROW(err.ne.0)

      ENDDO

    ENDIF


    IF (coef%nco2 > 0) THEN

      DO i = 1, coef%fmv_chn
        WRITE (file_id, iostat=err)coef%thermal%co2(:, i, :)

        THROW(err.ne.0)

      ENDDO

    ENDIF


    IF (coef%nn2o > 0) THEN

      DO i = 1, coef%fmv_chn
        WRITE (file_id, iostat=err)coef%thermal%n2o(:, i, :)

        THROW(err.ne.0)

      ENDDO

    ENDIF


    IF (coef%nco > 0) THEN

      DO i = 1, coef%fmv_chn
        WRITE (file_id, iostat=err)coef%thermal%co(:, i, :)

        THROW(err.ne.0)

      ENDDO

    ENDIF


    IF (coef%nch4 > 0) THEN

      DO i = 1, coef%fmv_chn
        WRITE (file_id, iostat=err)coef%thermal%ch4(:, i, :)

        THROW(err.ne.0)

      ENDDO

    ENDIF
    
  
  ! SOLAR_FAST_COEFFICIENTS
  
  WRITE (file_id, iostat=err)coef%solarcoef
  
  IF (coef%solarcoef) THEN
  
    IF (coef%nmixed > 0) THEN

      DO i = 1, coef%fmv_chn
        WRITE (file_id, iostat=err)coef%solar%mixedgas(:, i, :)

        THROW(err.ne.0)

      ENDDO

    ENDIF


    IF (coef%nwater > 0) THEN

      DO i = 1, coef%fmv_chn
        WRITE (file_id, iostat=err)coef%solar%watervapour(:, i, :)

        THROW(err.ne.0)

      ENDDO

    ENDIF


    IF (coef%nozone > 0) THEN

      DO i = 1, coef%fmv_chn
        WRITE (file_id, iostat=err)coef%solar%ozone(:, i, :)

        THROW(err.ne.0)

      ENDDO

    ENDIF


    IF (coef%nwvcont > 0) THEN

      DO i = 1, coef%fmv_chn
        WRITE (file_id, iostat=err)coef%solar%wvcont(:, i, :)

        THROW(err.ne.0)

      ENDDO

    ENDIF


    IF (coef%nco2 > 0) THEN

      DO i = 1, coef%fmv_chn
        WRITE (file_id, iostat=err)coef%solar%co2(:, i, :)

        THROW(err.ne.0)

      ENDDO

    ENDIF


    IF (coef%nn2o > 0) THEN

      DO i = 1, coef%fmv_chn
        WRITE (file_id, iostat=err)coef%solar%n2o(:, i, :)

        THROW(err.ne.0)

      ENDDO

    ENDIF


    IF (coef%nco > 0) THEN

      DO i = 1, coef%fmv_chn
        WRITE (file_id, iostat=err)coef%solar%co(:, i, :)

        THROW(err.ne.0)

      ENDDO

    ENDIF


    IF (coef%nch4 > 0) THEN

      DO i = 1, coef%fmv_chn
        WRITE (file_id, iostat=err)coef%solar%ch4(:, i, :)

        THROW(err.ne.0)

      ENDDO

    ENDIF
  ENDIF
  
!PLANCK_WEIGHTED
    section_present = Any(coef%pw_val_chn > 0)
    WRITE (file_id, iostat=err)section_present

    THROW(err.ne.0)

    IF (section_present) THEN
      WRITE (file_id, iostat=err)coef%pw_chn, coef%pw_val_chn
    ENDIF
  
!README_SRF
    section_present = coef%readme_srf(1)(1:4) .NE. "xxxx"
    WRITE (file_id, iostat=err)section_present

    THROW(err.ne.0)

    IF (section_present) THEN
      WRITE (file_id, iostat=err)coef%readme_srf
    ENDIF

!NLTE_RADIANCE_COEFS
    section_present = coef%nltecoef
    WRITE (file_id, iostat=err) section_present

    THROW(err.ne.0)

    IF (section_present) THEN
      WRITE (file_id, iostat=err) coef%nlte_coef%ncoef, &
      coef%nlte_coef%nsol, coef%nlte_coef%nsat, coef%nlte_coef%nchan, &
      coef%nlte_coef%start_chan

      WRITE (file_id, iostat=err) coef%nlte_coef%sol_zen_angle
      WRITE (file_id, iostat=err) coef%nlte_coef%sec_sat
      WRITE(file_id, iostat=ERR) coef%nlte_coef%coef
      
    ENDIF

!PRESSURE_MODULATED_CELL
    section_present = coef%pmc_shift
    WRITE (file_id, iostat=err) section_present

    THROW(err.ne.0)

    IF (section_present) THEN
      WRITE (file_id, iostat=err) coef%pmc_lengthcell
      WRITE (file_id, iostat=err) coef%pmc_pnominal 
      WRITE (file_id, iostat=err) coef%pmc_tempcell, coef%pmc_betaplus1, &
          coef%pmc_nlay, coef%pmc_nvar
      DO i = 1, coef%fmv_chn
        WRITE (file_id, iostat=err)coef%pmc_coef(:, i, :)
        THROW(err.ne.0)
      ENDDO
    ENDIF


    IF (lverbose) INFO("end of write coefficient")
    CATCH
  END SUBROUTINE rttov_write_binary_coef
