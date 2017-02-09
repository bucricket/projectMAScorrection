SUBROUTINE rttov_write_binary_scaercoef( &
            & err,           &
            & coef,          &
            & coef_scatt_ir, &
            & optp,          &
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
!  1.7       06/03/2009  Conditionals depending on coef%id_comp_lvl == 9
!                        extended to >= 9 (P. Rayer)
!  1.8       02/12/2009  Introduced principal component capability. Marco Matricardi. ECMWF
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
       & rttov_coef,          &
       & rttov_optpar_ir,     &
       & rttov_coef_scatt_ir
  USE parkind1, ONLY : jpim, jplm
!INTF_OFF
  USE rttov_const, ONLY :  &
       & rttov_magic_string, &
       & rttov_magic_number
!INTF_ON
  IMPLICIT NONE
! subroutine arguments
! scalar arguments with intent(in):
  INTEGER(KIND=jpim)       , INTENT(IN)              :: file_id      ! file logical unit number
  TYPE(rttov_coef         ), INTENT(IN)              :: coef         ! coefficients
  TYPE(rttov_optpar_ir    ), INTENT(IN)              :: optp
  TYPE(rttov_coef_scatt_ir), INTENT(IN)              :: coef_scatt_ir
  LOGICAL(KIND=jplm)       , INTENT(IN), OPTIONAL    :: verbose      ! print out information messages or not
! scalar arguments with intent(in):
  INTEGER(KIND=jpim)       , INTENT(OUT)             :: err          ! return code
!INTF_END
#include "rttov_errorreport.interface"
! local scalars
  INTEGER(KIND=jpim)  :: n
  INTEGER(KIND=jpim)  :: aer_pha_ioff, aer_pha_chn
  LOGICAL(KIND=jplm)  :: lverbose
  CHARACTER(LEN = 80) :: errMessage
!- End of header --------------------------------------------------------
  TRY
  IF (PRESENT(verbose)) THEN
    lverbose = verbose
  ELSE
    lverbose = .TRUE._jplm
  END IF

  IF (ALL(coef%ss_val_chn(:) == 0)) THEN
    aer_pha_ioff = 0
    aer_pha_chn  = 0
  ELSE
    aer_pha_ioff = coef_scatt_ir%fmv_aer_pha_ioff
    aer_pha_chn  = coef_scatt_ir%fmv_aer_pha_chn
  ENDIF

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

    WRITE (file_id, iostat=err)coef_scatt_ir%fmv_aer_chn, aer_pha_chn, aer_pha_ioff     &
      & , coef_scatt_ir%fmv_aer_comp, coef_scatt_ir%fmv_aer_ph

    THROW(err.ne.0)

    IF (aer_pha_ioff == 0 .AND. aer_pha_chn > 0) THEN
      WRITE (file_id, iostat=err)coef_scatt_ir%aer_pha_chanlist(:)

      THROW(err.ne.0)
    ENDIF

    WRITE (file_id, iostat=err)coef_scatt_ir%fmv_aer_ph_val

    THROW(err.ne.0)

    WRITE (file_id, iostat=err)coef_scatt_ir%fmv_aer_rh

    THROW(err.ne.0)


    DO n = 1, coef_scatt_ir%fmv_aer_comp
      WRITE (file_id, iostat=err)optp%optpaer(n)%fmv_aer_rh_val

      THROW(err.ne.0)

    ENDDO

!AEROSOLS_PARAMETERS

    DO n = 1, coef_scatt_ir%fmv_aer_comp
      WRITE (file_id, iostat=err)optp%optpaer(n)%abs

      THROW(err.ne.0)

    ENDDO


    DO n = 1, coef_scatt_ir%fmv_aer_comp
      WRITE (file_id, iostat=err)optp%optpaer(n)%sca

      THROW(err.ne.0)

    ENDDO


    DO n = 1, coef_scatt_ir%fmv_aer_comp
      WRITE (file_id, iostat=err)optp%optpaer(n)%bpr

      THROW(err.ne.0)

    ENDDO


    IF (aer_pha_chn > 0) THEN

      DO n = 1, coef_scatt_ir%fmv_aer_comp
        WRITE (file_id, iostat=err)optp%optpaer(n)%pha

        THROW(err.ne.0)

      ENDDO

    ENDIF

  IF (lverbose) INFO("end of write coefficient")
  CATCH
END SUBROUTINE 
