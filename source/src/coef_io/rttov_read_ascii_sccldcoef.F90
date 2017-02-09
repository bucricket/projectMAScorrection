!
SUBROUTINE rttov_read_ascii_sccldcoef( &
            & err,           &
            & coef,          &
            & coef_scatt_ir, &
            & optp,          &
            & file_lu,       &
            & channels)
! Description:
!
! Read an ASCII coefficient file and fills coeff structure
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
!  1.2       24/01/2003  Add return when section END encountered (P Brunel)
!                        any I/O error is coded as fatal
!                        Add GAZ_UNITS section
!  1.3       02/06/2004  New format for FMV section with RTTOV8 (P. Brunel)
!  1.4       15/06/2004  Corrected array dimension for coef % fmv_gas_pos (R Saunders)
!  1.5        June 2005  Added Additional arrays for RTTOV8M Marco Matricardi
!  1.6       05/12/2005  Corrected some array types for optical  const (R Saunders)
!  1.7       19/03/2007  Reduced no of continuation lines to below 39 (R Saunders)
!  1.8       12/12/2007  Add option of using top level (R Saunders)
!  1.9       01/11/2007  Removed hardcoded section length (A Geer)
!  1.10      26/06/2008  Introduced the case where no channels are available for the
!                        phase function in the solar range (M. Matricardi)
!  1.11      27/02/2009  Profile levels to include ToA. Allocate coef arrays
!                        according to number of layers, not levels (P. Rayer)
!  1.12      06/03/2009  Separation of flags for IncZeeman and IncTop.
!                        Conditionals depending on coef % id_comp_lvl == 9
!                        extended to >= 9 (P. Rayer)
!  1.13      08/06/2009  Made interim fix to allocate cloud/aerosol arrays with right shape (R Saunders)
!                        Ideally the channel order in all IR files needs to be in frequency not wavelength
!  1.14      02/12/2009  Introduced principal component capability. Marco Matricardi. ECMWF
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
       & deg2rad,          &
       & lensection,       &
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
  INTEGER(KIND=jpim)       , INTENT(OUT)               :: err           ! return code
!INTF_END
#include "rttov_errorreport.interface"
#include "rttov_skipcommentline.interface"
#include "rttov_findnextsection.interface"
#include "rttov_nullify_coef_scatt_ir.interface"
! Local Scalars:
  INTEGER(KIND=jpim) :: file_channels
  INTEGER(KIND=jpim) :: n_phase_channels
  LOGICAL(KIND=jplm) :: all_channels
  INTEGER(KIND=jpim) :: io_status
  INTEGER(KIND=jpim) :: i, j, k, n, nrh, nr, icount
  INTEGER(KIND=jpim), ALLOCATABLE :: list_of_channels(:)
  INTEGER(KIND=jpim), ALLOCATABLE :: phase_channels(:)     ! Solar channel numbers in original file
  INTEGER(KIND=jpim), ALLOCATABLE :: phase_channels_ext(:) ! Solar channel numbers to extract
  INTEGER(KIND=jpim), ALLOCATABLE :: phase_ext_index(:)    ! Indexes of phase data to extract
! pointers for generic inputs
  REAL   (KIND=jprb), POINTER :: abs_wcl_array(:, :   )
  REAL   (KIND=jprb), POINTER :: sca_wcl_array(:, :   )
  REAL   (KIND=jprb), POINTER :: bpr_wcl_array(:, :   )
  REAL   (KIND=jprb), POINTER :: pha_wcl_array(:, :, :)
  REAL   (KIND=jprb), POINTER :: abs_icl_array(:, :   )
  REAL   (KIND=jprb), POINTER :: sca_icl_array(:, :   )
  REAL   (KIND=jprb), POINTER :: bpr_icl_array(:, :   )
  REAL   (KIND=jprb), POINTER :: pha_icl_array(:, :, :)
  CHARACTER(LEN = 32)         :: wcl_comp_name
  CHARACTER(LEN = 32)         :: icl_comp_name
  CHARACTER(LEN = lensection) :: section
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

!read the file

  readfile : DO
    CALL rttov_findnextsection(file_lu, io_status, section)
    IF (io_status < 0) EXIT!end-of-file

    SELECT CASE (Trim(section))

    CASE ('WATERCLOUD_TYPES')

      CALL rttov_skipcommentline(file_lu, ERR)

      THROWM(err.ne.0,"io status while reading section "//section)

! number of channels for which optical parameters are stored
      READ (file_lu,  * , iostat=ERR)coef_scatt_ir%fmv_wcl_chn

      THROWM(err.ne.0,"io status while reading section "//section)

      IF (coef%fmv_ori_nchn /= coef_scatt_ir%fmv_wcl_chn) THEN
        err = errorstatus_fatal
        THROWM(err.ne.0,"Incompatible channels between rtcoef and sccldcoef files")
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
! coef % fmv_chn is the number of channels that the user requests
      file_channels = coef_scatt_ir%fmv_wcl_chn

      IF (.NOT. all_channels) THEN
        coef_scatt_ir%fmv_wcl_chn = Size(channels)
      ENDIF

! Number of channels for which phase function values are stored
      READ (file_lu,  * , iostat=ERR)n_phase_channels

      THROWM(err.ne.0,"io status while reading section "//section)

! index of first channel for which phase function values are available
      READ (file_lu,  * , iostat=ERR)coef_scatt_ir%fmv_wcl_pha_ioff

      THROWM(err.ne.0,"io status while reading section "//section)

! Sort out the solar channels/phase functions
      ! If coef_scatt_ir%fmv_wcl_pha_ioff > 0 then the phase functions are stored for a
      ! contiguous set of channels in the coef file.
      ! If coef_scatt_ir%fmv_wcl_pha_ioff == 0 then the channel numbers for which phase
      ! functions are stored are explicitly listed in the coef file.
      IF (n_phase_channels > 0) THEN
        ALLOCATE(phase_channels(n_phase_channels))

        ! Determine the solar channel numbers
        IF (coef_scatt_ir%fmv_wcl_pha_ioff > 0) THEN
          ! This maintains backward compatibility
          phase_channels(:) = (/ (i, i = coef_scatt_ir%fmv_wcl_pha_ioff, &
                                         coef_scatt_ir%fmv_wcl_pha_ioff + n_phase_channels - 1) /)
          coef_scatt_ir%fmv_wcl_pha_ioff = 0   ! Don't use this any more
        ELSE
          ! The solar channel indexes are written in the file
          CALL rttov_skipcommentline(file_lu, ERR)
          THROWM(err.ne.0,"io status while reading section "//section)
          READ (file_lu,  *, iostat=ERR)phase_channels(:)
          THROWM(err.ne.0,"io status while reading section "//section)
        ENDIF

        ! Determine solar channels/phase functions to be extracted
        ALLOCATE (phase_channels_ext(n_phase_channels))
        ALLOCATE (phase_ext_index(n_phase_channels))
        coef_scatt_ir%fmv_wcl_pha_chn = 0 ! count up number of solar chans to extract as we go
        DO i = 1, SIZE(list_of_channels)
          j = 1  ! here j is the index into the phase_channels array (list of solar channels in file)
          DO WHILE (j <= n_phase_channels)
            IF (list_of_channels(i) == phase_channels(j)) THEN
              coef_scatt_ir%fmv_wcl_pha_chn = coef_scatt_ir%fmv_wcl_pha_chn + 1
              phase_channels_ext(coef_scatt_ir%fmv_wcl_pha_chn) = i
              phase_ext_index(coef_scatt_ir%fmv_wcl_pha_chn) = j
              EXIT
            ENDIF
            j = j + 1
          ENDDO
        ENDDO

        IF (coef_scatt_ir%fmv_wcl_pha_chn > 0) THEN
          ! Copy solar channels to extract into correctly-sized array
          ALLOCATE (coef_scatt_ir%wcl_pha_chanlist(coef_scatt_ir%fmv_wcl_pha_chn))
          coef_scatt_ir%wcl_pha_chanlist(:) = phase_channels_ext(1:coef_scatt_ir%fmv_wcl_pha_chn)

          ! Create map from extracted channel list into pha array
          ALLOCATE (coef_scatt_ir%wcl_pha_index(SIZE(list_of_channels)))
          coef_scatt_ir%wcl_pha_index(:) = 0
          coef_scatt_ir%wcl_pha_index(coef_scatt_ir%wcl_pha_chanlist(1:coef_scatt_ir%fmv_wcl_pha_chn)) = &
                (/ (i, i = 1, coef_scatt_ir%fmv_wcl_pha_chn) /)
        ENDIF
        DEALLOCATE(phase_channels_ext)

        ! At this point:
        !   n_phase_channels                  = total number of solar channels in file
        !   phase_channels(:)                 = list of solar channel numbers
        !   coef_scatt_ir%fmv_wcl_pha_chn     = number of solar chans being extracted
        !   coef_scatt_ir%wcl_pha_chanlist(:) = list of solar channels to extract from file (indexes into extracted
        !                                         chan list, NOT original channel numbers)
        !   phase_ext_index(:)                = list of indexes into the phase fns which are to be extracted to pha
        !   coef_scatt_ir%wcl_pha_index(:)    = indexes for each extracted channel into the pha array

        IF (ALLOCATED(phase_channels)) DEALLOCATE(phase_channels)
      ELSE
        coef_scatt_ir%fmv_wcl_pha_chn = 0
      ENDIF

! number of water cloud types
      READ (file_lu,  * , iostat=ERR)coef_scatt_ir%fmv_wcl_comp

      THROWM(err.ne.0,"io status while reading section "//section)

! number of angles for phase function for water cloud types
      READ (file_lu,  * , iostat=ERR)coef_scatt_ir%fmv_wcl_ph

      THROWM(err.ne.0,"io status while reading section "//section)

      ALLOCATE (coef_scatt_ir%fmv_wcl_ph_val(coef_scatt_ir%fmv_wcl_ph), STAT = ERR)

      THROWM( ERR .NE. 0, "allocation of fmv_wcl_ph_val")

      ALLOCATE (coef_scatt_ir%fmv_wcl_ph_val_cos(coef_scatt_ir%fmv_wcl_ph), STAT = ERR)

      THROWM( ERR .NE. 0, "allocation of fmv_wcl_ph_val_cos")

      coef_scatt_ir%fmv_wcl_ph_val = 0.
      READ (file_lu,  * , iostat=ERR)(coef_scatt_ir%fmv_wcl_ph_val(i), i = 1, coef_scatt_ir%fmv_wcl_ph)

      THROWM(err.ne.0,"io status while reading section "//section)


      DO i = 1, coef_scatt_ir%fmv_wcl_ph
        coef_scatt_ir%fmv_wcl_ph_val_cos(i) = cos(coef_scatt_ir%fmv_wcl_ph_val(i) * deg2rad)
      ENDDO

      coef_scatt_ir%fmv_wcl_ph_val_min = 99999.0_JPRB

      DO i = 1, coef_scatt_ir%fmv_wcl_ph - 1

        IF (coef_scatt_ir%fmv_wcl_ph_val_min > coef_scatt_ir%fmv_wcl_ph_val(i + 1) - &
                                               coef_scatt_ir%fmv_wcl_ph_val(i)) THEN
          coef_scatt_ir%fmv_wcl_ph_val_min = coef_scatt_ir%fmv_wcl_ph_val(i + 1) - coef_scatt_ir%fmv_wcl_ph_val(i)
        ENDIF

      ENDDO

      coef_scatt_ir%fmv_wcl_ph_val_min = coef_scatt_ir%fmv_wcl_ph_val_min / 2.0_JPRB
      icount = int (coef_scatt_ir%fmv_wcl_ph_val(coef_scatt_ir%fmv_wcl_ph) / coef_scatt_ir%fmv_wcl_ph_val_min, jpim)
      ALLOCATE (coef_scatt_ir%ifmv_wcl_ph_val(icount), STAT = ERR)

      THROWM( ERR .NE. 0, "allocation of ifmv_wcl_ph_val")

      k = 1

      DO i = 1, icount

        IF (coef_scatt_ir%fmv_wcl_ph_val(k) >= (i + 1) * coef_scatt_ir%fmv_wcl_ph_val_min) THEN
          coef_scatt_ir%ifmv_wcl_ph_val(i) = k
        ELSE
          k = k + 1
          coef_scatt_ir%ifmv_wcl_ph_val(i) = k
        ENDIF

      ENDDO

! allocate arrays of FAST_MODEL_VARIABLES section
      ALLOCATE (coef_scatt_ir%fmv_wcl_comp_name(coef_scatt_ir%fmv_wcl_comp), STAT = ERR)

      THROWM( ERR .NE. 0, "allocation of fmv_wcl_comp_name")

      ALLOCATE (coef_scatt_ir%fmv_wcl_rh(coef_scatt_ir%fmv_wcl_comp), STAT = ERR)

      THROWM( ERR .NE. 0, "allocation of fmv_wcl_rh")

      ALLOCATE (coef_scatt_ir%confac(coef_scatt_ir%fmv_wcl_comp), STAT = ERR)

      THROWM( ERR .NE. 0, "allocation of confac")

      ALLOCATE (optp%optpwcl(coef_scatt_ir%fmv_wcl_comp), STAT = ERR)

      THROWM( ERR .NE. 0, "allocation of fmv_wcl_comp")

      DO n = 1, coef_scatt_ir%fmv_wcl_comp
        CALL rttov_nullify_coef_scatt_ir (optp%optpwcl(n))

        READ (file_lu, '(a5)', iostat=ERR)wcl_comp_name
        coef_scatt_ir%fmv_wcl_comp_name(n) = TRIM(ADJUSTL(wcl_comp_name))
        THROWM(err.ne.0,"io status while reading section "//section)

        READ (file_lu,  * , iostat=ERR)coef_scatt_ir%fmv_wcl_rh(n)

        THROWM(err.ne.0,"io status while reading section "//section)

        ALLOCATE (optp%optpwcl(n)%fmv_wcl_rh_val(coef_scatt_ir%fmv_wcl_rh(n)), STAT = ERR)

        THROWM( ERR .NE. 0, "allocation of fmv_wcl_rh_val")

        READ (file_lu,  * , iostat=ERR)(optp%optpwcl(n)%fmv_wcl_rh_val(i), i = 1, coef_scatt_ir%fmv_wcl_rh(n))

        THROWM(err.ne.0,"io status while reading section "//section)

        READ (file_lu,  * , iostat=ERR)coef_scatt_ir%confac(n)

        THROWM(err.ne.0,"io status while reading section "//section)

! Transfer information to some "classical" variables
! with more common names
! Note that the number of levels is taken from the Mixed Gases line
      ENDDO

      DEALLOCATE (list_of_channels)
    CASE ('WATERCLOUD_PARAMETERS')
      CALL rttov_skipcommentline(file_lu, ERR)

      THROWM(err.ne.0,"io status while reading section "//section)

! loop on cloud types

      DO n = 1, coef_scatt_ir%fmv_wcl_comp
        ALLOCATE (optp%optpwcl(n)%abs(coef_scatt_ir%fmv_wcl_chn, coef_scatt_ir%fmv_wcl_rh(n)), STAT = ERR)

        THROWM( ERR .NE. 0, "allocation of optpwcl(n) % abs")

        ALLOCATE (optp%optpwcl(n)%sca(coef_scatt_ir%fmv_wcl_chn, coef_scatt_ir%fmv_wcl_rh(n)), STAT = ERR)

        THROWM( ERR .NE. 0, "allocation of optpwcl(n) % sca")

        ALLOCATE (optp%optpwcl(n)%bpr(coef_scatt_ir%fmv_wcl_chn, coef_scatt_ir%fmv_wcl_rh(n)), STAT = ERR)

        THROWM( ERR .NE. 0, "allocation of optpwcl(n) % bpr")

        IF (coef_scatt_ir%fmv_wcl_pha_chn > 0) THEN
          ALLOCATE (optp%optpwcl(n)%pha(coef_scatt_ir%fmv_wcl_pha_chn, coef_scatt_ir%fmv_wcl_rh(n), &
                                        coef_scatt_ir%fmv_wcl_ph), STAT = ERR)

          THROWM( ERR .NE. 0, "allocation of optpwcl(n) % pha")
        ENDIF

        IF (all_channels) THEN
          abs_wcl_array => optp%optpwcl(n)%abs
          sca_wcl_array => optp%optpwcl(n)%sca
          bpr_wcl_array => optp%optpwcl(n)%bpr
          pha_wcl_array => optp%optpwcl(n)%pha
        ELSE
          ALLOCATE (abs_wcl_array(file_channels, coef_scatt_ir%fmv_wcl_rh(n)), STAT = ERR)

          THROWM( ERR .NE. 0, "allocation of abs_wcl_array")

          ALLOCATE (sca_wcl_array(file_channels, coef_scatt_ir%fmv_wcl_rh(n)), STAT = ERR)

          THROWM( ERR .NE. 0, "allocation of sca_wcl_array")

          ALLOCATE (bpr_wcl_array(file_channels, coef_scatt_ir%fmv_wcl_rh(n)), STAT = ERR)

          THROWM( ERR .NE. 0, "allocation of bpr_wcl_array")

          IF (n_phase_channels > 0) THEN
            ALLOCATE (pha_wcl_array(n_phase_channels, coef_scatt_ir%fmv_wcl_rh(n), coef_scatt_ir%fmv_wcl_ph), &
                      STAT = ERR)

            THROWM( ERR .NE. 0, "allocation of pha_wcl_array")
          ENDIF
        ENDIF


        DO nrh = 1, coef_scatt_ir%fmv_wcl_rh(n)
          CALL rttov_skipcommentline(file_lu, ERR)

          THROWM(err.ne.0,"io status while reading section "//section)

! read dummy string of gas name or filename of the sub_coefficient file
          READ (file_lu,  * , iostat=ERR)wcl_comp_name

          THROWM(err.ne.0,"io status while reading section "//section)

          READ (file_lu,  * , iostat=ERR)(abs_wcl_array(i, nrh), i = 1, file_channels)

          THROWM(err.ne.0,"io status while reading section "//section)

          READ (file_lu,  * , iostat=ERR)(sca_wcl_array(i, nrh), i = 1, file_channels)

          THROWM(err.ne.0,"io status while reading section "//section)

          READ (file_lu,  * , iostat=ERR)(bpr_wcl_array(i, nrh), i = 1, file_channels)

          THROWM(err.ne.0,"io status while reading section "//section)

          IF (n_phase_channels > 0) THEN
            READ (file_lu,  * , iostat=ERR)     &
              & ((pha_wcl_array(i, nrh, j), j = 1, coef_scatt_ir%fmv_wcl_ph), i = 1, n_phase_channels)

            THROWM(err.ne.0,"io status while reading section "//section)
          ENDIF

        ENDDO


        IF (.NOT. all_channels) THEN
          optp%optpwcl(n)%abs(:,:)   = abs_wcl_array(channels(:), :)
          optp%optpwcl(n)%sca(:,:)   = sca_wcl_array(channels(:), :)
          optp%optpwcl(n)%bpr(:,:)   = bpr_wcl_array(channels(:), :)
          IF ( coef_scatt_ir%fmv_wcl_pha_chn > 0 ) THEN
            optp%optpwcl(n)%pha(:,:,:) = pha_wcl_array(phase_ext_index(1:coef_scatt_ir%fmv_wcl_pha_chn), :, :)
          ENDIF
          DEALLOCATE (abs_wcl_array, STAT = ERR)

          THROWM( ERR .NE. 0, "deallocation of abs_wcl_array")

          DEALLOCATE (sca_wcl_array, STAT = ERR)

          THROWM( ERR .NE. 0, "deallocation of sca_wcl_array")

          DEALLOCATE (bpr_wcl_array, STAT = ERR)

          THROWM( ERR .NE. 0, "deallocation of bpr_wcl_array")

          IF (coef_scatt_ir%fmv_wcl_pha_chn > 0) DEALLOCATE (pha_wcl_array, STAT = ERR)

          THROWM( ERR .NE. 0, "deallocation of pha_wcl_array")

        ENDIF

      ENDDO

      IF (ALLOCATED(phase_ext_index)) DEALLOCATE(phase_ext_index)

!        print*,optpaer%optpar(11) % abs (1,1)
!        print*,optpaer%optpar(11) % sca (1,1)
!        print*,optpaer%optpar(11) % pha (1,1,1:208)
!        print*,optpaer%optpar(11) % pha (3041,1,1:208)
!-------------------------------------------------------
    CASE ('ICECLOUD_TYPES')

      CALL rttov_skipcommentline(file_lu, ERR)

      THROWM(err.ne.0,"io status while reading section "//section)

! number of channels for which regression coefficients are stored
      READ (file_lu,  * , iostat=ERR)coef_scatt_ir%fmv_icl_chn

      THROWM(err.ne.0,"io status while reading section "//section)

      IF (coef%fmv_ori_nchn /= coef_scatt_ir%fmv_icl_chn) THEN
        err = errorstatus_fatal
        THROWM(err.ne.0,"Incompatible channels between rtcoef and sccldcoef files")
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
! coef % fmv_chn is the number of channels that the user requests
      file_channels = coef_scatt_ir%fmv_icl_chn

      IF (.NOT. all_channels) THEN
        coef_scatt_ir%fmv_icl_chn = Size(channels)
      ENDIF

! Number of channels for which phase function values are stored
      READ (file_lu,  * , iostat=ERR)n_phase_channels

      THROWM(err.ne.0,"io status while reading section "//section)

! index of first channel for which phase function values are available
      READ (file_lu,  * , iostat=ERR)coef_scatt_ir%fmv_icl_pha_ioff

      THROWM(err.ne.0,"io status while reading section "//section)

! Sort out the solar channels/phase functions
      ! If coef_scatt_ir%fmv_icl_pha_ioff > 0 then the phase functions are stored for a
      ! contiguous set of channels in the coef file.
      ! If coef_scatt_ir%fmv_icl_pha_ioff == 0 then the channel numbers for which phase
      ! functions are stored are explicitly listed in the coef file.
      IF (n_phase_channels > 0) THEN
        ALLOCATE(phase_channels(n_phase_channels))

        ! Determine the solar channel numbers
        IF (coef_scatt_ir%fmv_icl_pha_ioff > 0) THEN
          ! This maintains backward compatibility
          phase_channels(:) = (/ (i, i = coef_scatt_ir%fmv_icl_pha_ioff, &
                                         coef_scatt_ir%fmv_icl_pha_ioff + n_phase_channels - 1) /)
          coef_scatt_ir%fmv_icl_pha_ioff = 0   ! Don't use this any more
        ELSE
          ! The solar channel indexes are written in the file
          CALL rttov_skipcommentline(file_lu, ERR)
          THROWM(err.ne.0,"io status while reading section "//section)
          READ (file_lu,  * , iostat=ERR)phase_channels(:)
          THROWM(err.ne.0,"io status while reading section "//section)
        ENDIF

        ! Determine solar channels/phase functions to be extracted
        ALLOCATE (phase_channels_ext(n_phase_channels))
        ALLOCATE (phase_ext_index(n_phase_channels))
        coef_scatt_ir%fmv_icl_pha_chn = 0 ! count up number of solar chans to extract as we go
        DO i = 1, SIZE(list_of_channels)
          j = 1  ! here j is the index into the phase_channels array (list of solar channels in file)
          DO WHILE (j <= n_phase_channels)
            IF (list_of_channels(i) == phase_channels(j)) THEN
              coef_scatt_ir%fmv_icl_pha_chn = coef_scatt_ir%fmv_icl_pha_chn + 1
              phase_channels_ext(coef_scatt_ir%fmv_icl_pha_chn) = i
              phase_ext_index(coef_scatt_ir%fmv_icl_pha_chn) = j
              EXIT
            ENDIF
            j = j + 1
          ENDDO
        ENDDO

        IF (coef_scatt_ir%fmv_icl_pha_chn > 0) THEN
          ! Copy solar channels to extract into correctly-sized array
          ALLOCATE (coef_scatt_ir%icl_pha_chanlist(coef_scatt_ir%fmv_icl_pha_chn))
          coef_scatt_ir%icl_pha_chanlist(:) = phase_channels_ext(1:coef_scatt_ir%fmv_icl_pha_chn)

          ! Create map from extracted channel list into pha array
          ALLOCATE (coef_scatt_ir%icl_pha_index(SIZE(list_of_channels)))
          coef_scatt_ir%icl_pha_index(:) = 0
          coef_scatt_ir%icl_pha_index(coef_scatt_ir%icl_pha_chanlist(1:coef_scatt_ir%fmv_icl_pha_chn)) = &
                (/ (i, i = 1, coef_scatt_ir%fmv_icl_pha_chn) /)
        ENDIF
        DEALLOCATE(phase_channels_ext)

        ! At this point:
        !   n_phase_channels                  = total number of solar channels in file
        !   phase_channels(:)                 = list of solar channel numbers
        !   coef_scatt_ir%fmv_icl_pha_chn     = number of solar chans being extracted
        !   coef_scatt_ir%icl_pha_chanlist(:) = list of solar channels to extract from file (indexes into extracted
        !                                         chan list, NOT original channel numbers)
        !   phase_ext_index(:)                = list of indexes into the phase fns which are to be extracted to pha
        !   coef_scatt_ir%icl_pha_index(:)    = indexes for each extracted channel into the pha array

        IF (ALLOCATED(phase_channels)) DEALLOCATE(phase_channels)
      ELSE
        coef_scatt_ir%fmv_icl_pha_chn = 0
      ENDIF

! number of coefficients used in the regression for absorption optical depth
      READ (file_lu,  * , iostat=ERR)coef_scatt_ir%icl_nabs

      THROWM(err.ne.0,"io status while reading section "//section)

! number of coefficients used in the regression for scattering optical depth
      READ (file_lu,  * , iostat=ERR)coef_scatt_ir%icl_nsca

      THROWM(err.ne.0,"io status while reading section "//section)

! number of coefficients used in the regression for backscattering parameter
      READ (file_lu,  * , iostat=ERR)coef_scatt_ir%icl_nbpr

      THROWM(err.ne.0,"io status while reading section "//section)

! number of size distributions used in the regression
      READ (file_lu,  * , iostat=ERR)coef_scatt_ir%fmv_icl_comp

      THROWM(err.ne.0,"io status while reading section "//section)

! number of ice crystal shapes for which parameters are available
      READ (file_lu,  * , iostat=ERR)coef_scatt_ir%fmv_icl_ishp

      THROWM(err.ne.0,"io status while reading section "//section)

! allocation of arrays
      ALLOCATE (coef_scatt_ir%fmv_icl_dg(coef_scatt_ir%fmv_icl_comp, coef_scatt_ir%fmv_icl_ishp), STAT = ERR)

      THROWM( ERR .NE. 0, "allocation of fmv_icl_dg")

      ALLOCATE (coef_scatt_ir%fmv_icl_comp_name(coef_scatt_ir%fmv_icl_comp, 2), STAT = ERR)

      THROWM( ERR .NE. 0, "allocation of fmv_icl_comp_name")

      ALLOCATE (optp%optpicl(coef_scatt_ir%fmv_icl_ishp), STAT = ERR)

      THROWM( ERR .NE. 0, "allocation of optpicl")


      DO n = 1, coef_scatt_ir%fmv_icl_ishp
        CALL rttov_nullify_coef_scatt_ir (optp%optpicl(n))
      ENDDO

! number of angles for phase function for ice clouds
      READ (file_lu,  * , iostat=ERR)coef_scatt_ir%fmv_icl_ph

      THROWM(err.ne.0,"io status while reading section "//section)

      ALLOCATE (coef_scatt_ir%fmv_icl_ph_val(coef_scatt_ir%fmv_icl_ph), STAT = ERR)

      THROWM( ERR .NE. 0, "allocation of fmv_icl_ph_val")

      ALLOCATE (coef_scatt_ir%fmv_icl_ph_val_cos(coef_scatt_ir%fmv_icl_ph), STAT = ERR)

      THROWM( ERR .NE. 0, "allocation of fmv_icl_ph_val_cos")

      coef_scatt_ir%fmv_icl_ph_val = 0.
      READ (file_lu,  * , iostat=ERR)(coef_scatt_ir%fmv_icl_ph_val(i), i = 1, coef_scatt_ir%fmv_icl_ph)

      THROWM(err.ne.0,"io status while reading section "//section)


      DO i = 1, coef_scatt_ir%fmv_icl_ph
        coef_scatt_ir%fmv_icl_ph_val_cos(i) = cos(coef_scatt_ir%fmv_icl_ph_val(i) * deg2rad)
      ENDDO

      coef_scatt_ir%fmv_icl_ph_val_min = 99999.0_JPRB

      DO i = 1, coef_scatt_ir%fmv_icl_ph - 1

        IF (coef_scatt_ir%fmv_icl_ph_val_min > coef_scatt_ir%fmv_icl_ph_val(i + 1) - &
                                               coef_scatt_ir%fmv_icl_ph_val(i)) THEN
          coef_scatt_ir%fmv_icl_ph_val_min = coef_scatt_ir%fmv_icl_ph_val(i + 1) - coef_scatt_ir%fmv_icl_ph_val(i)
        ENDIF

      ENDDO

      coef_scatt_ir%fmv_icl_ph_val_min = coef_scatt_ir%fmv_icl_ph_val_min / 2.0_JPRB
      icount = int (coef_scatt_ir%fmv_icl_ph_val(coef_scatt_ir%fmv_icl_ph) / coef_scatt_ir%fmv_icl_ph_val_min, jpim)
      ALLOCATE (coef_scatt_ir%ifmv_icl_ph_val(icount), STAT = ERR)

      THROWM( ERR .NE. 0, "allocation of ifmv_icl_ph_val")

      k = 1

      DO i = 1, icount

        IF (coef_scatt_ir%fmv_icl_ph_val(k) >= (i + 1) * coef_scatt_ir%fmv_icl_ph_val_min) THEN
          coef_scatt_ir%ifmv_icl_ph_val(i) = k
        ELSE
          k = k + 1
          coef_scatt_ir%ifmv_icl_ph_val(i) = k
        ENDIF

      ENDDO

      DEALLOCATE (list_of_channels)

    CASE ('HEXAGONAL_PARAMETERS')
      CALL rttov_skipcommentline(file_lu, ERR)

      THROWM(err.ne.0,"io status while reading section "//section)

!Read the effective generalized diameter

      DO i = 1, coef_scatt_ir%fmv_icl_comp
        READ (file_lu,  * , iostat=ERR)coef_scatt_ir%fmv_icl_dg(i, 1)

        THROWM(err.ne.0,"io status while reading section "//section)

      ENDDO

      ALLOCATE (optp%optpicl(1)%abs(coef_scatt_ir%fmv_icl_chn, coef_scatt_ir%icl_nabs), STAT = ERR)

      THROWM( ERR .NE. 0, "allocation of optpicl(1) % abs")

      ALLOCATE (optp%optpicl(1)%sca(coef_scatt_ir%fmv_icl_chn, coef_scatt_ir%icl_nsca), STAT = ERR)

      THROWM( ERR .NE. 0, "allocation of optpicl(1) % sca")

      ALLOCATE (optp%optpicl(1)%bpr(coef_scatt_ir%fmv_icl_chn, coef_scatt_ir%icl_nbpr), STAT = ERR)

      THROWM( ERR .NE. 0, "allocation of optpicl(1) % bpr")

      IF (coef_scatt_ir%fmv_icl_pha_chn > 0) THEN
        ALLOCATE (optp%optpicl(1)%pha(coef_scatt_ir%fmv_icl_pha_chn, coef_scatt_ir%fmv_icl_comp, &
                                      coef_scatt_ir%fmv_icl_ph), STAT = ERR)

        THROWM( ERR .NE. 0, "allocation of optpicl(1) % pha")
      ENDIF

      IF (all_channels) THEN
        abs_icl_array => optp%optpicl(1)%abs
        sca_icl_array => optp%optpicl(1)%sca
        bpr_icl_array => optp%optpicl(1)%bpr
        pha_icl_array => optp%optpicl(1)%pha
      ELSE
        ALLOCATE (abs_icl_array(file_channels, coef_scatt_ir%icl_nabs), STAT = ERR)

        THROWM( ERR .NE. 0, "allocation of abs_icl_array")

        ALLOCATE (sca_icl_array(file_channels, coef_scatt_ir%icl_nsca), STAT = ERR)

        THROWM( ERR .NE. 0, "allocation of sca_icl_array")

        ALLOCATE (bpr_icl_array(file_channels, coef_scatt_ir%icl_nbpr), STAT = ERR)

        THROWM( ERR .NE. 0, "allocation of bpr_icl_array")

        IF (n_phase_channels > 0) THEN
          ALLOCATE (pha_icl_array(n_phase_channels, coef_scatt_ir%fmv_icl_comp, coef_scatt_ir%fmv_icl_ph), STAT = ERR)

          THROWM( ERR .NE. 0, "allocation of pha_icl_array")
        ENDIF
      ENDIF

      CALL rttov_skipcommentline(file_lu, ERR)

      THROWM(err.ne.0,"io status while reading section "//section)


      DO i = 1, file_channels
        READ (file_lu,  * , iostat=ERR)(abs_icl_array(i, n), n = 1, coef_scatt_ir%icl_nabs)

        THROWM(err.ne.0,"io status while reading section "//section)

      ENDDO


      DO i = 1, file_channels
        READ (file_lu,  * , iostat=ERR)(sca_icl_array(i, n), n = 1, coef_scatt_ir%icl_nsca)

        THROWM(err.ne.0,"io status while reading section "//section)

      ENDDO


      DO i = 1, file_channels
        READ (file_lu,  * , iostat=ERR)(bpr_icl_array(i, n), n = 1, coef_scatt_ir%icl_nbpr)

        THROWM(err.ne.0,"io status while reading section "//section)

      ENDDO


      DO nr = 1, coef_scatt_ir%fmv_icl_comp
        CALL rttov_skipcommentline(file_lu, ERR)

        THROWM(err.ne.0,"io status while reading section "//section)

        READ (file_lu,  '(a17)', iostat=ERR)icl_comp_name
        coef_scatt_ir%fmv_icl_comp_name(nr, 1) = TRIM(ADJUSTL(icl_comp_name))
        THROWM(err.ne.0,"io status while reading section "//section)


        IF (n_phase_channels > 0) THEN
          READ (file_lu,  * , iostat=ERR)     &
            & ((pha_icl_array(i, nr, j), j = 1, coef_scatt_ir%fmv_icl_ph), i = 1, n_phase_channels)

          THROWM(err.ne.0,"io status while reading section "//section)
        ENDIF

      ENDDO


      IF (.NOT. all_channels) THEN
        optp%optpicl(1)%abs(:,:)   = abs_icl_array(channels(:), :)
        optp%optpicl(1)%sca(:,:)   = sca_icl_array(channels(:), :)
        optp%optpicl(1)%bpr(:,:)   = bpr_icl_array(channels(:), :)
        IF ( coef_scatt_ir%fmv_icl_pha_chn > 0 ) THEN
          optp%optpicl(1)%pha(:,:,:) = pha_icl_array(phase_ext_index(1:coef_scatt_ir%fmv_icl_pha_chn), :, :)
        ENDIF
        DEALLOCATE (abs_icl_array, STAT = ERR)

        THROWM( ERR .NE. 0, "deallocation of abs_icl_array")

        DEALLOCATE (sca_icl_array, STAT = ERR)

        THROWM( ERR .NE. 0, "deallocation of sca_icl_array")

        DEALLOCATE (bpr_icl_array, STAT = ERR)

        THROWM( ERR .NE. 0, "deallocation of bpr_icl_array")

        IF (coef_scatt_ir%fmv_icl_pha_chn > 0) DEALLOCATE (pha_icl_array, STAT = ERR)

        THROWM( ERR .NE. 0, "deallocation of pha_icl_array")

      ENDIF

    CASE ('AGGREGATE_PARAMETERS')
      CALL rttov_skipcommentline(file_lu, ERR)

      THROWM(err.ne.0,"io status while reading section "//section)

!Read the effective generalized diameter

      DO i = 1, coef_scatt_ir%fmv_icl_comp
        READ (file_lu,  * , iostat=ERR)coef_scatt_ir%fmv_icl_dg(i, 2)

        THROWM(err.ne.0,"io status while reading section "//section)

      ENDDO

      ALLOCATE (optp%optpicl(2)%abs(coef_scatt_ir%fmv_icl_chn, coef_scatt_ir%icl_nabs), STAT = ERR)

      THROWM( ERR .NE. 0, "allocation of optpicl(2) % abs")

      ALLOCATE (optp%optpicl(2)%sca(coef_scatt_ir%fmv_icl_chn, coef_scatt_ir%icl_nsca), STAT = ERR)

      THROWM( ERR .NE. 0, "allocation of optpicl(2) % sca")

      ALLOCATE (optp%optpicl(2)%bpr(coef_scatt_ir%fmv_icl_chn, coef_scatt_ir%icl_nbpr), STAT = ERR)

      THROWM( ERR .NE. 0, "allocation of optpicl(2) % bpr")

      IF (coef_scatt_ir%fmv_icl_pha_chn > 0) THEN
        ALLOCATE (optp%optpicl(2)%pha(coef_scatt_ir%fmv_icl_pha_chn, coef_scatt_ir%fmv_icl_comp, &
                                      coef_scatt_ir%fmv_icl_ph), STAT = ERR)

        THROWM( ERR .NE. 0, "allocation of optpicl(2) % pha")
      ENDIF

      IF (all_channels) THEN
        abs_icl_array => optp%optpicl(2)%abs
        sca_icl_array => optp%optpicl(2)%sca
        bpr_icl_array => optp%optpicl(2)%bpr
        pha_icl_array => optp%optpicl(2)%pha
      ELSE
        ALLOCATE (abs_icl_array(file_channels, coef_scatt_ir%icl_nabs), STAT = ERR)

        THROWM( ERR .NE. 0, "allocation of abs_icl_array")

        ALLOCATE (sca_icl_array(file_channels, coef_scatt_ir%icl_nsca), STAT = ERR)

        THROWM( ERR .NE. 0, "allocation of sca_icl_array")

        ALLOCATE (bpr_icl_array(file_channels, coef_scatt_ir%icl_nbpr), STAT = ERR)

        THROWM( ERR .NE. 0, "allocation of bpr_icl_array")

        IF (n_phase_channels > 0) THEN
          ALLOCATE (pha_icl_array(n_phase_channels, coef_scatt_ir%fmv_icl_comp, coef_scatt_ir%fmv_icl_ph), STAT = ERR)

          THROWM( ERR .NE. 0, "allocation of pha_icl_array")
        ENDIF
      ENDIF

      CALL rttov_skipcommentline(file_lu, ERR)

      THROWM(err.ne.0,"io status while reading section "//section)


      DO i = 1, file_channels
        READ (file_lu,  * , iostat=ERR)(abs_icl_array(i, n), n = 1, coef_scatt_ir%icl_nabs)

        THROWM(err.ne.0,"io status while reading section "//section)

      ENDDO


      DO i = 1, file_channels
        READ (file_lu,  * , iostat=ERR)(sca_icl_array(i, n), n = 1, coef_scatt_ir%icl_nabs)

        THROWM(err.ne.0,"io status while reading section "//section)

      ENDDO


      DO i = 1, file_channels
        READ (file_lu,  * , iostat=ERR)(bpr_icl_array(i, n), n = 1, coef_scatt_ir%icl_nsca)

        THROWM(err.ne.0,"io status while reading section "//section)

      ENDDO


      DO nr = 1, coef_scatt_ir%fmv_icl_comp
        CALL rttov_skipcommentline(file_lu, ERR)

        THROWM(err.ne.0,"io status while reading section "//section)

        READ (file_lu,  '(a17)', iostat=ERR)icl_comp_name
        coef_scatt_ir%fmv_icl_comp_name(nr, 2) = TRIM(ADJUSTL(icl_comp_name))
        THROWM(err.ne.0,"io status while reading section "//section)


        IF (n_phase_channels > 0) THEN
          READ (file_lu,  * , iostat=ERR)     &
            & ((pha_icl_array(i, nr, j), j = 1, coef_scatt_ir%fmv_icl_ph), i = 1, n_phase_channels)
        ENDIF

        THROWM(err.ne.0,"io status while reading section "//section)

      ENDDO


      IF (.NOT. all_channels) THEN
        optp%optpicl(2)%abs(:,:)   = abs_icl_array(channels(:), :)
        optp%optpicl(2)%sca(:,:)   = sca_icl_array(channels(:), :)
        optp%optpicl(2)%bpr(:,:)   = bpr_icl_array(channels(:), :)
        IF ( coef_scatt_ir%fmv_icl_pha_chn > 0 ) THEN
          optp%optpicl(2)%pha(:,:,:) = pha_icl_array(phase_ext_index(1:coef_scatt_ir%fmv_icl_pha_chn), :, :)
        ENDIF
        DEALLOCATE (abs_icl_array, STAT = ERR)

        THROWM( ERR .NE. 0, "deallocation of abs_icl_array")

        DEALLOCATE (sca_icl_array, STAT = ERR)

        THROWM( ERR .NE. 0, "deallocation of sca_icl_array")

        DEALLOCATE (bpr_icl_array, STAT = ERR)

        THROWM( ERR .NE. 0, "deallocation of bpr_icl_array")

        IF (coef_scatt_ir%fmv_icl_pha_chn > 0) DEALLOCATE (pha_icl_array, STAT = ERR)

        THROWM( ERR .NE. 0, "deallocation of pha_icl_array")

      ENDIF

      IF (ALLOCATED(phase_ext_index)) DEALLOCATE(phase_ext_index)

    CASE ('END')
      RETURN
    CASE DEFAULT
      CYCLE readfile
    END SELECT

  ENDDO readfile

  CATCH
END SUBROUTINE rttov_read_ascii_sccldcoef
