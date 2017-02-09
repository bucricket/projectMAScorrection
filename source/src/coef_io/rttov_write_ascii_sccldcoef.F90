SUBROUTINE rttov_write_ascii_sccldcoef ( &
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
  INTEGER(KIND=jpim)  :: i, n, nr, nrh
  INTEGER(KIND=jpim)  :: wcl_pha_ioff, wcl_pha_chn
  INTEGER(KIND=jpim)  :: icl_pha_ioff, icl_pha_chn
  LOGICAL(KIND=jplm)  :: lverbose
  CHARACTER(LEN = 32) :: section
  CHARACTER(LEN = 80) :: errMessage
!- End of header --------------------------------------------------------
  TRY
  IF (PRESENT(verbose)) THEN
    lverbose = verbose
  ELSE
    lverbose = .TRUE._jplm
  END IF

  ! Ensure we don't write out phase functions unnecessarily
  IF (ALL(coef%ss_val_chn(:) == 0)) THEN
    wcl_pha_ioff = 0
    wcl_pha_chn  = 0
    icl_pha_ioff = 0
    icl_pha_chn  = 0
  ELSE
    wcl_pha_ioff = coef_scatt_ir%fmv_wcl_pha_ioff
    wcl_pha_chn  = coef_scatt_ir%fmv_wcl_pha_chn
    icl_pha_ioff = coef_scatt_ir%fmv_icl_pha_ioff
    icl_pha_chn  = coef_scatt_ir%fmv_icl_pha_chn
  ENDIF

!ASCII file
  IF (lverbose) THEN
    WRITE (errMessage, '( "write coefficient to file_id ", i2, " in ASCII format")')file_id
    INFO(errMessage)
  END IF
  WRITE (file_id, '(a)', iostat=err)' ! RTTOV coefficient file '//Trim(coef%id_Common_name)

  THROW(err.ne.0)

  WRITE (file_id, '(a)', iostat=err)' ! automatic creation by subroutine Rttov_writecoef '

  THROW(err.ne.0)

  WRITE (file_id, '(a)', iostat=err)' ! ------------------------------------------------------'

  THROW(err.ne.0)

! COEF structure (V10)

  section = 'WATERCLOUD_TYPES'
  WRITE (file_id, '(a)', iostat=err)' ! ------------------------------------------------------'

  THROW(err.ne.0)

  WRITE (file_id, '(a)', iostat=err)Trim(section)

  THROW(err.ne.0)

  WRITE (file_id, '(a)', iostat=err)' !'

  THROW(err.ne.0)

  WRITE (file_id, '(1x,i8,t20,a)', iostat=err)coef_scatt_ir%fmv_wcl_chn, &
      &'! number of channels for which optical parameters are stored'

  THROW(err.ne.0)

  WRITE (file_id, '(1x,i8,t20,a)', iostat=err)wcl_pha_chn,&
      &'! Number of channels for which phase function values are stored'

  THROW(err.ne.0)

  WRITE (file_id, '(1x,i8,t20,a)', iostat=err)wcl_pha_ioff,&
      &'! index of first channel for which phase function values are available'

  THROW(err.ne.0)

  IF (wcl_pha_ioff == 0 .AND. wcl_pha_chn > 0) THEN
    WRITE(file_id,'(1x,a)') '! Channel list for which phase function values are available'
    WRITE(file_id,'(10i6)') coef_scatt_ir%wcl_pha_chanlist(1:wcl_pha_chn)
  ENDIF

  WRITE (file_id, '(1x,i8,t20,a)', iostat=err)coef_scatt_ir%fmv_wcl_comp, &
      &'! number of water cloud types'

  THROW(err.ne.0)

  WRITE (file_id, '(1x,i8,t20,a)', iostat=err)coef_scatt_ir%fmv_wcl_ph, &
      &'! number of angles for phase function for water cloud types'

  THROW(err.ne.0)

  WRITE (file_id, '(10f7.2)', iostat=err)coef_scatt_ir%fmv_wcl_ph_val

  THROW(err.ne.0)


  DO n = 1, coef_scatt_ir%fmv_wcl_comp
    IF (LEN_TRIM(coef_scatt_ir%fmv_wcl_comp_name(n)) > 0) THEN
      WRITE (file_id, '(a5,i2.2)', iostat=err)TRIM(coef_scatt_ir%fmv_wcl_comp_name(n))
    ELSE
      WRITE (file_id, '(a)', iostat=err)' cloud   ! default name for rttov_writecoef'
    ENDIF

    THROW(err.ne.0)

    WRITE (file_id, '(1x,i8,t20,a)', iostat=err)coef_scatt_ir%fmv_wcl_rh(n),&
        &'!RH values for which parameters are available'

    THROW(err.ne.0)

    WRITE (file_id, '(10f7.2)', iostat=err)optp%optpwcl(n)%fmv_wcl_rh_val

    THROW(err.ne.0)

    WRITE (file_id, '(1x,f12.6,t20,a)', iostat=err)coef_scatt_ir%confac(n), &
        & '!Conversion from LWC to particle density'

    THROW(err.ne.0)

  ENDDO

  section = 'WATERCLOUD_PARAMETERS'
  WRITE (file_id, '(a)', iostat=err)' ! ------------------------------------------------------'

  THROW(err.ne.0)

  WRITE (file_id, '(a)', iostat=err)Trim(section)

  THROW(err.ne.0)

  WRITE (file_id, '(a)', iostat=err)' !'

  THROW(err.ne.0)


  DO n = 1, coef_scatt_ir%fmv_wcl_comp

    DO nrh = 1, coef_scatt_ir%fmv_wcl_rh(n)
      IF (LEN_TRIM(coef_scatt_ir%fmv_wcl_comp_name(n)) > 0) THEN
        WRITE (file_id, '(a5,i2.2)', iostat=err)TRIM(coef_scatt_ir%fmv_wcl_comp_name(n)), &
                                                INT(optp%optpwcl(n)%fmv_wcl_rh_val(nrh))
      ELSE
        WRITE (file_id, '(a)', iostat=err)' cloud   ! default name for rttov_writecoef'
      ENDIF

      THROW(err.ne.0)

      WRITE (file_id, '(5e16.8)', iostat=err)optp%optpwcl(n)%abs(:, nrh)

      THROW(err.ne.0)

      WRITE (file_id, '(5e16.8)', iostat=err)optp%optpwcl(n)%sca(:, nrh)

      THROW(err.ne.0)

      WRITE (file_id, '(5e16.8)', iostat=err)optp%optpwcl(n)%bpr(:, nrh)

      THROW(err.ne.0)

!Write(file_id, '(5e16.8)') &
!& ((optp%optpwcl(n)%pha(i,nrh,j),&
!& j = 1, coef_scatt_ir%fmv_wcl_ph),i=1,wcl_pha_chn)

      IF (wcl_pha_chn > 0) THEN

        DO i = 1, wcl_pha_chn
          WRITE (file_id, '(5e16.8)', iostat=err)optp%optpwcl(n)%pha(i, nrh, :)

          THROW(err.ne.0)

        ENDDO

      ENDIF

    ENDDO
  
  ENDDO

  section = 'ICECLOUD_TYPES'
  WRITE (file_id, '(a)', iostat=err)' ! ------------------------------------------------------'

  THROW(err.ne.0)

  WRITE (file_id, '(a)', iostat=err)Trim(section)

  THROW(err.ne.0)

  WRITE (file_id, '(a)', iostat=err)' !'

  THROW(err.ne.0)

  WRITE (file_id, '(1x,i8,t20,a)', iostat=err)coef_scatt_ir%fmv_icl_chn, &
      & '! number of channels for which regression coefficients are stored'

  THROW(err.ne.0)

  WRITE (file_id, '(1x,i8,t20,a)', iostat=err)icl_pha_chn,&
      & '! Number of channels for which phase function values are stored'

  THROW(err.ne.0)

  WRITE (file_id, '(1x,i8,t20,a)', iostat=err)icl_pha_ioff,&
      & '! index of first channel for which phase function values are available'

  THROW(err.ne.0)

  IF (icl_pha_ioff == 0 .AND. icl_pha_chn > 0) THEN
    WRITE(file_id,'(1x,a)') '! Channel list for which phase function values are available'
    WRITE(file_id,'(10i6)') coef_scatt_ir%icl_pha_chanlist(1:icl_pha_chn)
  ENDIF

  WRITE (file_id, '(1x,i8,t20,a)', iostat=err)coef_scatt_ir%icl_nabs,&
      & '! number of coefficients used in the regression for absorption optical depth'

  THROW(err.ne.0)

  WRITE (file_id, '(1x,i8,t20,a)', iostat=err)coef_scatt_ir%icl_nsca, &
      & '! number of coefficients used in the regression for scattering optical depth'

  THROW(err.ne.0)

  WRITE (file_id, '(1x,i8,t20,a)', iostat=err)coef_scatt_ir%icl_nbpr,&
      & '! number of coefficients used in the regression for backscattering parameter' 

  THROW(err.ne.0)

  WRITE (file_id, '(1x,i8,t20,a)', iostat=err)coef_scatt_ir%fmv_icl_comp, &
      & '! number of size distributions used in the regression' 

  THROW(err.ne.0)

  WRITE (file_id, '(1x,i8,t20,a)', iostat=err)coef_scatt_ir%fmv_icl_ishp,&
      & '! number of ice crystal shapes for which parameters are available'

  THROW(err.ne.0)

  WRITE (file_id, '(1x,i8,t20,a)', iostat=err)coef_scatt_ir%fmv_icl_ph,&
      & '! number of angles for phase function for ice clouds' 

  THROW(err.ne.0)

  WRITE (file_id, '(10f7.2)', iostat=err)coef_scatt_ir%fmv_icl_ph_val

  THROW(err.ne.0)

  section = 'HEXAGONAL_PARAMETERS'
  WRITE (file_id, '(a)', iostat=err)' ! ------------------------------------------------------'

  THROW(err.ne.0)

  WRITE (file_id, '(a)', iostat=err)Trim(section)

  THROW(err.ne.0)

  WRITE (file_id, '(a)', iostat=err)' !'

  THROW(err.ne.0)

  WRITE (file_id, '(a)', iostat=err)' ! Effective diameter for each size distribution'

  THROW(err.ne.0)

  WRITE (file_id, '(a)', iostat=err)' !'

  THROW(err.ne.0)

  WRITE (file_id, '(f10.4)', iostat=err)coef_scatt_ir%fmv_icl_dg(:, 1)

  THROW(err.ne.0)

  WRITE (file_id, '(a)', iostat=err)' !'

  THROW(err.ne.0)

  WRITE (file_id, '(a)', iostat=err)' ! Regression coefficients for ice clouds'

  THROW(err.ne.0)

  WRITE (file_id, '(a)', iostat=err)' !'

  THROW(err.ne.0)


  DO i = 1, coef_scatt_ir%fmv_icl_chn
    WRITE (file_id, '(4E16.8)', iostat=err)optp%optpicl(1)%abs(i, :)

    THROW(err.ne.0)

  ENDDO


  DO i = 1, coef_scatt_ir%fmv_icl_chn
    WRITE (file_id, '(4E16.8)', iostat=err)optp%optpicl(1)%sca(i, :)

    THROW(err.ne.0)

  ENDDO


  DO i = 1, coef_scatt_ir%fmv_icl_chn
    WRITE (file_id, '(4E16.8)', iostat=err)optp%optpicl(1)%bpr(i, :)

    THROW(err.ne.0)

  ENDDO


  DO nr = 1, coef_scatt_ir%fmv_icl_comp
    WRITE (file_id, '(a)', iostat=err)' !'

    THROW(err.ne.0)

    WRITE (file_id, '(a)', iostat=err)' ! Phase function values'

    THROW(err.ne.0)

    WRITE (file_id, '(a)', iostat=err)' !'

    THROW(err.ne.0)

    IF (LEN_TRIM(coef_scatt_ir%fmv_icl_comp_name(nr, 1)) > 0) THEN
      WRITE (file_id, '(a17)', iostat=err)TRIM(coef_scatt_ir%fmv_icl_comp_name(nr, 1))
    ELSE
      WRITE (file_id, '(a)', iostat=err)' Hexagonal   ! default name for rttov_writecoef'
    ENDIF

    THROW(err.ne.0)

!Write(file_id, '(5e16.8)') &
!  & ((optp%optpicl(1)%pha(i,nr,j),&
!  & j = 1, coef_scatt_ir%fmv_icl_ph),i=1,icl_pha_chn)

    IF (icl_pha_chn > 0) THEN

      DO i = 1, icl_pha_chn
        WRITE (file_id, '(5e16.8)', iostat=err)optp%optpicl(1)%pha(i, nr, :)

        THROW(err.ne.0)

      ENDDO

    ENDIF

  ENDDO

  section = 'AGGREGATE_PARAMETERS'
  WRITE (file_id, '(a)', iostat=err)' ! ------------------------------------------------------'

  THROW(err.ne.0)

  WRITE (file_id, '(a)', iostat=err)Trim(section)

  THROW(err.ne.0)

  WRITE (file_id, '(a)', iostat=err)' !'

  THROW(err.ne.0)

  WRITE (file_id, '(a)', iostat=err)' !'

  THROW(err.ne.0)

  WRITE (file_id, '(a)', iostat=err)' ! Effective diameter for each size distribution'

  THROW(err.ne.0)

  WRITE (file_id, '(a)', iostat=err)' !'

  THROW(err.ne.0)

  WRITE (file_id, '(f10.4)', iostat=err)coef_scatt_ir%fmv_icl_dg(:, 2)

  THROW(err.ne.0)

  WRITE (file_id, '(a)', iostat=err)' !'

  THROW(err.ne.0)

  WRITE (file_id, '(a)', iostat=err)' ! Regression coefficients for ice clouds'

  THROW(err.ne.0)

  WRITE (file_id, '(a)', iostat=err)' !'

  THROW(err.ne.0)


  DO i = 1, coef_scatt_ir%fmv_icl_chn
    WRITE (file_id, '(4e16.8)', iostat=err)optp%optpicl(2)%abs(i, :)

    THROW(err.ne.0)

  ENDDO


  DO i = 1, coef_scatt_ir%fmv_icl_chn
    WRITE (file_id, '(4e16.8)', iostat=err)optp%optpicl(2)%sca(i, :)

    THROW(err.ne.0)

  ENDDO


  DO i = 1, coef_scatt_ir%fmv_icl_chn
    WRITE (file_id, '(4e16.8)', iostat=err)optp%optpicl(2)%bpr(i, :)

    THROW(err.ne.0)

  ENDDO


  DO nr = 1, coef_scatt_ir%fmv_icl_comp
    WRITE (file_id, '(a)', iostat=err)' !'

    THROW(err.ne.0)

    WRITE (file_id, '(a)', iostat=err)' ! Phase function values'

    THROW(err.ne.0)

    WRITE (file_id, '(a)', iostat=err)' !'

    THROW(err.ne.0)

    IF (LEN_TRIM(coef_scatt_ir%fmv_icl_comp_name(nr, 2)) > 0) THEN
      WRITE (file_id, '(a17)', iostat=err)TRIM(coef_scatt_ir%fmv_icl_comp_name(nr, 2))
    ELSE
      WRITE (file_id, '(a)', iostat=err)' aggregate   ! default name for rttov_writecoef'
    ENDIF

    THROW(err.ne.0)

!Write(file_id, '(5e16.8)') &
!  & ((optp%optpicl(2)%pha(i,nr,j),&
!  & j = 1, coef_scatt_ir%fmv_icl_ph),i=1,icl_pha_chn)

    IF (icl_pha_chn > 0) THEN

      DO i = 1, icl_pha_chn
        WRITE (file_id, '(5e16.8)', iostat=err)optp%optpicl(2)%pha(i, nr, :)

        THROW(err.ne.0)

      ENDDO

    ENDIF

  ENDDO

  IF (lverbose) INFO("end of write coefficient")
  CATCH
END SUBROUTINE 
