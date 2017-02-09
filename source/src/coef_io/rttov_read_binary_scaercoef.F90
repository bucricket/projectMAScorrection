!
SUBROUTINE rttov_read_binary_scaercoef( &
            & ERR,           &
            & coef,          &
            & coef_scatt_ir, &
            & optp,          &
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
! Version   Date     Comment
! -------   ----     -------
!  1.0       01/12/2002  New F90 code with structures (P Brunel A Smith)
!  1.1       02/01/2003  A few comments added (R Saunders)
!  1.2       24/01/2003  add tests on all read statements (P Brunel)
!                        one record per channel for coefficients in binary format
!                        New header to allow checking R4<->R8
!  1.3       06/05/2003  Remove "optional" attribute of argument file_lu (P Brunel)
!  1.4       02/06/2004  New format for FMV section with RTTOV8  (P. Brunel)
!  1.5       15/06/2004  Corrected array dimension for coef % fmv_gas_pos (R Saunders)
!  1.6       14/05/2007  Updated for RTTOV-9 ( P Brunel)
!  1.7       12/12/2007  Option of additional layer at top flag read (R Saunders)
!  1.8       27/06/2008  Introduced the case when no channels are present for the
!                        phase function in the solar range (M. Matricardi)
!  1.9       27/02/2009  Profile levels to include ToA. Allocate coef arrays
!                        according to number of layers, not levels (P. Rayer)
!  1.10      06/03/2009  Separation of flags for IncZeeman and IncTop.
!                        Conditionals depending on coef % id_comp_lvl == 9
!                        extended to >= 9 (P. Rayer)
!  1.11      08/06/2009  Made interim fix to allocate cloud/aerosol arrays with right shape (R Saunders)
!                        Ideally the channel order in all IR files needs to be in frequency not wavelength
!  1.12      02/12/2009  Introduced principal component capability. Marco Matricardi. ECMWF
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
  USE parkind1, ONLY : jpim
!INTF_OFF
  USE rttov_const, ONLY :  &
       & rttov_magic_string,     &
       & rttov_magic_number,     &
       & deg2rad,                &
       & errorstatus_fatal
  USE parkind1, ONLY : jprb, jplm
!INTF_ON
  IMPLICIT NONE
! subroutine arguments
! scalar arguments with intent(in):
  INTEGER(KIND=jpim)       , INTENT(IN)                :: file_lu       ! file logical unit number
  INTEGER(KIND=jpim)       , OPTIONAL     , INTENT(IN) :: channels   (:)! list of channels to extract
! scalar arguments with intent(inout):
  TYPE(rttov_coef         ), INTENT(INOUT)             :: coef          ! coefficients
  TYPE(rttov_optpar_ir    ), INTENT(INOUT)             :: optp
  TYPE(rttov_coef_scatt_ir), INTENT(INOUT)             :: coef_scatt_ir
! scalar arguments with intent(out):
  INTEGER(KIND=jpim)       , INTENT(OUT)               :: ERR           ! return code
!INTF_END
#include "rttov_errorreport.interface"
#include "rttov_nullify_coef_scatt_ir.interface"
! Local Scalars:
  INTEGER(KIND=jpim) :: file_channels
  LOGICAL(KIND=jplm) :: all_channels
  INTEGER(KIND=jpim) :: n
  INTEGER(KIND=jpim) :: i, j, k, icount
  INTEGER(KIND=jpim) :: n_phase_channels
  INTEGER(KIND=jpim), ALLOCATABLE :: list_of_channels(:)
  INTEGER(KIND=jpim), ALLOCATABLE :: phase_channels(:)     ! Solar channel numbers in original file
  INTEGER(KIND=jpim), ALLOCATABLE :: phase_channels_ext(:) ! Solar channel numbers to extract
  INTEGER(KIND=jpim), ALLOCATABLE :: phase_ext_index(:)    ! Indexes of phase data to extract
! pointers for generic inputs
  REAL   (KIND=jprb), POINTER     :: dvalues0        (:, :   )
  REAL   (KIND=jprb), POINTER     :: tvalues0        (:, :, :)
  CHARACTER(LEN = 16) :: bin_check_string
  REAL(KIND=jprb)     :: bin_check_number
  REAL(KIND=jprb)     :: bin_check_value
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


! COEF_SCATT_IR aerosols structure (V9)

!AEROSOLS_COMPONENTS
    READ (file_lu, iostat=ERR)coef_scatt_ir%fmv_aer_chn, n_phase_channels, coef_scatt_ir%fmv_aer_pha_ioff, &
      & coef_scatt_ir%fmv_aer_comp, coef_scatt_ir%fmv_aer_ph

    THROWM( ERR .NE. 0, "reading aerosols components")

    IF (coef%fmv_ori_nchn /= coef_scatt_ir%fmv_aer_chn) THEN
      err = errorstatus_fatal
      THROWM(err.ne.0,"Incompatible channels between rtcoef and scaercoef files")
    ENDIF

    IF (.NOT. all_channels) THEN
      ALLOCATE (list_of_channels(size(channels)))
      list_of_channels = channels
    ELSE
      ALLOCATE (list_of_channels(coef%fmv_chn))
      list_of_channels = (/(i, i = 1, coef%fmv_chn)/)
    ENDIF

! Take care of the user list of channels
! file_channels store the number of channels in the file
! coef % fmv_aer_chn is the number of channels that the user requests
    file_channels = coef_scatt_ir%fmv_aer_chn

    IF (.NOT. all_channels) THEN
      coef_scatt_ir%fmv_aer_chn = Size(channels)
    ENDIF

! Sort out the solar channels/phase functions
    IF (n_phase_channels > 0) THEN
      ALLOCATE(phase_channels(n_phase_channels))

      ! Determine the solar channel numbers
      IF (coef_scatt_ir%fmv_aer_pha_ioff > 0) THEN
        ! This maintains backward compatibility
        phase_channels(:) = (/ (i, i = coef_scatt_ir%fmv_aer_pha_ioff, &
                                       coef_scatt_ir%fmv_aer_pha_ioff + n_phase_channels - 1) /)
        coef_scatt_ir%fmv_aer_pha_ioff = 0   ! Don't use this any more
      ELSE
        ! The solar channel indexes are written in the file
        READ (file_lu, iostat=ERR)phase_channels(:)
        THROWM(err.ne.0,"reading solar channel numbers")
      ENDIF

      ! Determine solar channels/phase functions to be extracted
      ALLOCATE(phase_channels_ext(n_phase_channels))
      ALLOCATE(phase_ext_index(n_phase_channels))
      coef_scatt_ir%fmv_aer_pha_chn = 0 ! count up number of solar chans to extract as we go
      DO i = 1, SIZE(list_of_channels)
        j = 1  ! here j is the index into the phase_channels array (list of solar channels in file)
        DO WHILE (j <= n_phase_channels)
          IF (list_of_channels(i) == phase_channels(j)) THEN
            coef_scatt_ir%fmv_aer_pha_chn = coef_scatt_ir%fmv_aer_pha_chn + 1
            phase_channels_ext(coef_scatt_ir%fmv_aer_pha_chn) = i
            phase_ext_index(coef_scatt_ir%fmv_aer_pha_chn) = j
            EXIT
          ENDIF
          j = j + 1
        ENDDO
      ENDDO

      IF (coef_scatt_ir%fmv_aer_pha_chn > 0) THEN
        ! Copy solar channels to extract into correctly-sized array
        ALLOCATE(coef_scatt_ir%aer_pha_chanlist(coef_scatt_ir%fmv_aer_pha_chn))
        coef_scatt_ir%aer_pha_chanlist(:) = phase_channels_ext(1:coef_scatt_ir%fmv_aer_pha_chn)

        ! Create map from extracted channel list into pha array
        ALLOCATE(coef_scatt_ir%aer_pha_index(SIZE(list_of_channels)))
        coef_scatt_ir%aer_pha_index(:) = 0
        coef_scatt_ir%aer_pha_index(coef_scatt_ir%aer_pha_chanlist(1:coef_scatt_ir%fmv_aer_pha_chn)) = &
              (/ (i, i = 1, coef_scatt_ir%fmv_aer_pha_chn) /)
      ENDIF
      DEALLOCATE(phase_channels_ext)

      ! See rttov_read_ascii_scaercoef.F90 for a description of what the various arrays contain.

      IF (ALLOCATED(phase_channels)) DEALLOCATE(phase_channels)
    ELSE
      ! What do we need to do if there are no phase fns at all?
      coef_scatt_ir%fmv_aer_pha_chn = 0
    ENDIF

    ALLOCATE (coef_scatt_ir%fmv_aer_ph_val(coef_scatt_ir%fmv_aer_ph), STAT = ERR)

    THROWM( ERR .NE. 0, "allocation of fmv_aer_ph_val")

    ALLOCATE (coef_scatt_ir%fmv_aer_ph_val_cos(coef_scatt_ir%fmv_aer_ph), STAT = ERR)

    THROWM( ERR .NE. 0, "allocation of fmv_aer_ph_val_cos")

    READ (file_lu, iostat=ERR)coef_scatt_ir%fmv_aer_ph_val

    THROWM( ERR .NE. 0, "reading aerosols fmv_aer_ph_val")


    DO i = 1, coef_scatt_ir%fmv_aer_ph
      coef_scatt_ir%fmv_aer_ph_val_cos(i) = cos(coef_scatt_ir%fmv_aer_ph_val(i) * deg2rad)
    ENDDO

    coef_scatt_ir%fmv_aer_ph_val_min = 99999.0_JPRB

    DO i = 1, coef_scatt_ir%fmv_aer_ph - 1

      IF (coef_scatt_ir%fmv_aer_ph_val_min > coef_scatt_ir%fmv_aer_ph_val(i + 1) - &
                                             coef_scatt_ir%fmv_aer_ph_val(i)) THEN
        coef_scatt_ir%fmv_aer_ph_val_min = coef_scatt_ir%fmv_aer_ph_val(i + 1) - coef_scatt_ir%fmv_aer_ph_val(i)
      ENDIF

    ENDDO

    coef_scatt_ir%fmv_aer_ph_val_min = coef_scatt_ir%fmv_aer_ph_val_min / 2.0_JPRB
    icount = coef_scatt_ir%fmv_aer_ph_val(coef_scatt_ir%fmv_aer_ph) / coef_scatt_ir%fmv_aer_ph_val_min
    ALLOCATE (coef_scatt_ir%ifmv_aer_ph_val(icount), STAT = ERR)

    THROWM( ERR .NE. 0, "allocation of fmv_aer_ph_val")

    k = 1

    DO i = 1, icount

      IF (coef_scatt_ir%fmv_aer_ph_val(k) >= (i + 1) * coef_scatt_ir%fmv_aer_ph_val_min) THEN
        coef_scatt_ir%ifmv_aer_ph_val(i) = k
      ELSE
        k = k + 1
        coef_scatt_ir%ifmv_aer_ph_val(i) = k
      ENDIF

    ENDDO

! allocate arrays of FAST_MODEL_VARIABLES section
    ! aerosol type names are not stored in binary files, but allocate the array anyway
    ALLOCATE (coef_scatt_ir%fmv_aer_comp_name(coef_scatt_ir%fmv_aer_comp), STAT = ERR)

    THROWM( ERR .NE. 0, "allocation of fmv_aer_comp_name")

    coef_scatt_ir%fmv_aer_comp_name(:) = ""

    ALLOCATE (coef_scatt_ir%fmv_aer_rh(coef_scatt_ir%fmv_aer_comp), STAT = ERR)

    THROWM( ERR .NE. 0, "allocation of fmv_aer_rh")

    ALLOCATE (optp%optpaer(coef_scatt_ir%fmv_aer_comp), STAT = ERR)

    THROWM( ERR .NE. 0, "allocation of optpaer")


    DO n = 1, coef_scatt_ir%fmv_aer_comp
      CALL rttov_nullify_coef_scatt_ir (optp%optpaer(n))
    ENDDO

    READ (file_lu, iostat=ERR)coef_scatt_ir%fmv_aer_rh

    THROWM( ERR .NE. 0, "reading coef_scatt_ir % fmv_aer_rh")


    DO n = 1, coef_scatt_ir%fmv_aer_comp
      ALLOCATE (optp%optpaer(n)%fmv_aer_rh_val(coef_scatt_ir%fmv_aer_rh(n)), STAT = ERR)

      THROWM( ERR .NE. 0, "allocation of optp%optpaer(n)% fmv_aer_rh_val")

      READ (file_lu, iostat=ERR)optp%optpaer(n)%fmv_aer_rh_val

      THROWM( ERR .NE. 0, "reading optp%optpaer(n)% fmv_aer_rh_val")

    ENDDO

!'AEROSOLS_PARAMETERS'

    DO n = 1, coef_scatt_ir%fmv_aer_comp
      ALLOCATE (optp%optpaer(n)%abs(coef_scatt_ir%fmv_aer_chn, coef_scatt_ir%fmv_aer_rh(n)), STAT = ERR)

      THROWM( ERR .NE. 0, "allocation of optp%optpaer(n) % abs")


      IF (all_channels) THEN
        READ (file_lu, iostat=ERR)optp%optpaer(n)%abs

        THROWM( ERR .NE. 0, "reading optp%optpaer(n) % abs")

      ELSE
        ALLOCATE (dvalues0(file_channels, coef_scatt_ir%fmv_aer_rh(n)), STAT = ERR)

        THROWM( ERR .NE. 0, "allocation of dvalues0")

        READ (file_lu, iostat=ERR)dvalues0

        THROWM( ERR .NE. 0, "reading dvalues0")

        optp%optpaer(n)%abs(:,:) = dvalues0(channels(:), :)
        DEALLOCATE (dvalues0, STAT = ERR)

        THROWM( ERR .NE. 0, "deallocation of dvalues0")

      ENDIF

    ENDDO


    DO n = 1, coef_scatt_ir%fmv_aer_comp
      ALLOCATE (optp%optpaer(n)%sca(coef_scatt_ir%fmv_aer_chn, coef_scatt_ir%fmv_aer_rh(n)), STAT = ERR)

      THROWM( ERR .NE. 0, "allocation of optp%optpaer(n) % sca")


      IF (all_channels) THEN
        READ (file_lu, iostat=ERR)optp%optpaer(n)%sca

        THROWM( ERR .NE. 0, "reading optp%optpaer(n) % sca")

      ELSE
        ALLOCATE (dvalues0(file_channels, coef_scatt_ir%fmv_aer_rh(n)), STAT = ERR)

        THROWM( ERR .NE. 0, "allocation of dvalues0")

        READ (file_lu, iostat=ERR)dvalues0

        THROWM( ERR .NE. 0, "reading dvalues0")

        optp%optpaer(n)%sca(:,:) = dvalues0(channels(:), :)
        DEALLOCATE (dvalues0, STAT = ERR)

        THROWM( ERR .NE. 0, "deallocation of dvalues0")

      ENDIF

    ENDDO


    DO n = 1, coef_scatt_ir%fmv_aer_comp
      ALLOCATE (optp%optpaer(n)%bpr(coef_scatt_ir%fmv_aer_chn, coef_scatt_ir%fmv_aer_rh(n)), STAT = ERR)

      THROWM( ERR .NE. 0, "allocation of optp%optpaer(n) % bpr")


      IF (all_channels) THEN
        READ (file_lu, iostat=ERR)optp%optpaer(n)%bpr

        THROWM( ERR .NE. 0, "reading optp%optpaer(n) % bpr")

      ELSE
        ALLOCATE (dvalues0(file_channels, coef_scatt_ir%fmv_aer_rh(n)), STAT = ERR)

        THROWM( ERR .NE. 0, "allocation of dvalues0")

        READ (file_lu, iostat=ERR)dvalues0

        THROWM( ERR .NE. 0, "reading dvalues0")

        optp%optpaer(n)%bpr(:,:) = dvalues0(channels(:), :)
        DEALLOCATE (dvalues0, STAT = ERR)

        THROWM( ERR .NE. 0, "deallocation of dvalues0")

      ENDIF

    ENDDO


    DO n = 1, coef_scatt_ir%fmv_aer_comp
      IF (coef_scatt_ir%fmv_aer_pha_chn > 0) THEN
        ALLOCATE (optp%optpaer(n)%pha(coef_scatt_ir%fmv_aer_pha_chn, coef_scatt_ir%fmv_aer_rh(n), &
                                      coef_scatt_ir%fmv_aer_ph), STAT = ERR)

        THROWM( ERR .NE. 0, "allocation of optp%optpaer(n) % pha")
      ENDIF

      IF (n_phase_channels > 0) THEN

        IF (all_channels) THEN
          READ (file_lu, iostat=ERR)optp%optpaer(n)%pha

          THROWM( ERR .NE. 0, "reading optp%optpaer(n) % pha")

        ELSE
          ALLOCATE (tvalues0(n_phase_channels, coef_scatt_ir%fmv_aer_rh(n), coef_scatt_ir%fmv_aer_ph), STAT = ERR)

          THROWM( ERR .NE. 0, "allocation of tvalues0")

          READ (file_lu, iostat=ERR)tvalues0

          THROWM( ERR .NE. 0, "reading tvalues0")

          IF ( coef_scatt_ir%fmv_aer_pha_chn > 0 ) THEN
            optp%optpaer(n)%pha(:,:,:) = tvalues0(phase_ext_index(1:coef_scatt_ir%fmv_aer_pha_chn), :, :)
          END IF
          DEALLOCATE (tvalues0, STAT = ERR)

          THROWM( ERR .NE. 0, "deallocation of tvalues0")

        ENDIF

      ENDIF

    ENDDO

    IF (ALLOCATED(phase_ext_index)) DEALLOCATE (phase_ext_index)
    DEALLOCATE (list_of_channels)

!
! Here add reading of new sections for binary format in order to keep compatibility with
! previous versions
!
  CATCH
END SUBROUTINE rttov_read_binary_scaercoef
