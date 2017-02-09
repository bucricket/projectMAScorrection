! Description:
!> @file
!!   Read coefficient file(s) into RTTOV coefficients structure.
!
!> @brief
!!   Read coefficient file(s) into RTTOV coefficients structure.
!!
!! @details
!!   You can read files either by specifying filenames explicitly
!!   (the recommended method), by opening the file(s) and passing the
!!   relevant logical unit number(s), or by specifying the instrument
!!   ID triplet so that RTTOV can construct the filenames (see
!!   rttov_const.F90 for the platform and instrument IDs).
!!
!!   RTTOV can normally detect the format of coefficient files so
!!   the file format arguments are not usually required.
!!
!!   If you read a subset of n channels from the coefficient file
!!   these will be identified by indexes 1...n in chanprof(:)\%chan
!!   when calling RTTOV, not by the original channel numbers.
!!
!!   A similar consideration applies to the channels_rec channel list
!!   for PC reconstructed radiances.
!!
!!   Cloud coefficient (sccldcoef) files will only be read if the
!!   addclouds option is true and the user_cld_opt_param option is false.
!!
!!   Aerosol coefficient (scaercoef) files will only be read if the
!!   addaerosl option is true and the user_aer_opt_param option is false.
!!
!!   PC coefficient (pccoef) files will only be read if the addpc option
!!   is true.
!!
!!   All coefficient files being read must be in the same format.
!!
!! @param[out]    err             status on exit
!! @param[in,out] coefs           RTTOV coefs structure
!! @param[in]     opts            options to configure the simulations
!! @param[in]     channels        list of channels to read from file(s), optional
!! @param[in]     channels_rec    list of channels for which PC reconstructed radiances are required, optional
!! @param[in]     form_coef       format of rtcoef file, optional
!! @param[in]     form_scaer      format of IR aerosol scaer file, optional
!! @param[in]     form_sccld      format of IR cloud scaer file, optional
!! @param[in]     form_pccoef     format of PC-RTTOV pccoef file, optional
!! @param[in]     file_coef       file name of rtcoef file, optional
!! @param[in]     file_scaer      file name of IR aerosol scaer file, optional
!! @param[in]     file_sccld      file name of IR cloud scaer file, optional
!! @param[in]     file_pccoef     file name of PC-RTTOV pccoef file, optional
!! @param[in]     file_id_coef    logical unit for rtcoef file, optional
!! @param[in]     file_id_scaer   logical unit for IR aerosol scaer file, optional
!! @param[in]     file_id_sccld   logical unit for IR cloud scaer file, optional
!! @param[in]     file_id_pccoef  logical unit for PC-RTTOV pccoef file, optional
!! @param[in]     instrument      (platform,satellite,instrument) ID triplet, optional
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
SUBROUTINE rttov_read_coefs(err, coefs, opts,  &
                            channels,          &
                            channels_rec,      &
                            form_coef,         &
                            form_scaer,        &
                            form_sccld,        &
                            form_pccoef,       &
                            file_coef,         &
                            file_scaer,        &
                            file_sccld,        &
                            file_pccoef,       &
                            file_id_coef,      &
                            file_id_scaer,     &
                            file_id_sccld,     &
                            file_id_pccoef,    &
                            instrument)
#include "throw.h"

  USE rttov_types, ONLY : rttov_coefs, rttov_options
  USE parkind1, ONLY : jpim
!INTF_OFF
  USE parkind1, ONLY : jprb, jplm
  USE rttov_coef_io_mod, ONLY : getlun, closelun
  USE rttov_const, ONLY : inst_name
  USE yomhook, ONLY : LHOOK, DR_HOOK
#ifdef _RTTOV_HDF
  USE rttov_hdf_mod, ONLY : &
      open_hdf,             &
      close_hdf,            &
      is_hdf_open,          &
      is_hdf_64bit_reals
#endif
!INTF_ON

  IMPLICIT NONE

  INTEGER(KIND=jpim), INTENT(OUT) :: err
  TYPE(rttov_coefs),  INTENT(OUT) :: coefs
  TYPE(rttov_options),INTENT(IN)  :: opts
  INTEGER(KIND=jpim), INTENT(IN), OPTIONAL :: channels(:)
  INTEGER(KIND=jpim), INTENT(IN), OPTIONAL :: channels_rec(:)
  CHARACTER(len=*),   INTENT(IN), OPTIONAL :: form_coef
  CHARACTER(len=*),   INTENT(IN), OPTIONAL :: form_scaer
  CHARACTER(len=*),   INTENT(IN), OPTIONAL :: form_sccld
  CHARACTER(len=*),   INTENT(IN), OPTIONAL :: form_pccoef
  CHARACTER(len=*),   INTENT(IN), OPTIONAL :: file_coef
  CHARACTER(len=*),   INTENT(IN), OPTIONAL :: file_scaer
  CHARACTER(len=*),   INTENT(IN), OPTIONAL :: file_sccld
  CHARACTER(len=*),   INTENT(IN), OPTIONAL :: file_pccoef
  INTEGER(KIND=jpim), INTENT(IN), OPTIONAL :: file_id_coef
  INTEGER(KIND=jpim), INTENT(IN), OPTIONAL :: file_id_scaer
  INTEGER(KIND=jpim), INTENT(IN), OPTIONAL :: file_id_sccld
  INTEGER(KIND=jpim), INTENT(IN), OPTIONAL :: file_id_pccoef
  INTEGER(KIND=jpim), INTENT(IN), OPTIONAL :: instrument(3)
!INTF_END

#include "rttov_errorreport.interface"
#include "rttov_cmpuc.interface"

#include "rttov_read_ascii_coef.interface"
#include "rttov_read_ascii_pccoef.interface"
#include "rttov_read_ascii_scaercoef.interface"
#include "rttov_read_ascii_sccldcoef.interface"

#include "rttov_read_binary_coef.interface"
#include "rttov_read_binary_pccoef.interface"
#include "rttov_read_binary_scaercoef.interface"
#include "rttov_read_binary_sccldcoef.interface"

#include "rttov_channel_extract_coef.interface"
#include "rttov_channel_extract_scaercoef.interface"
#include "rttov_channel_extract_sccldcoef.interface"
#include "rttov_channel_extract_pccoef.interface"

#include "rttov_check_channels_pc.interface"

#include "rttov_init_coefs.interface"
#include "rttov_nullify_coefs.interface"

#include "rttov_dealloc_coef.interface"
#include "rttov_dealloc_coef_scatt_ir.interface"
#include "rttov_dealloc_optpar_ir.interface"
#include "rttov_dealloc_coef_pccomp.interface"

#ifdef _RTTOV_HDF
#include "rttov_hdf_load.interface"
#endif

  INTEGER(KIND=jpim) :: file_id_coef1=0_jpim
  INTEGER(KIND=jpim) :: file_id_scaer1=0_jpim
  INTEGER(KIND=jpim) :: file_id_sccld1=0_jpim
  INTEGER(KIND=jpim) :: file_id_pccoef1=0_jpim

  CHARACTER(LEN=256) :: file_coef1
  CHARACTER(LEN=256) :: file_scaer1
  CHARACTER(LEN=256) :: file_sccld1
  CHARACTER(LEN=256) :: file_pccoef1

  CHARACTER(LEN=32) :: form_coef1=""
  CHARACTER(LEN=32) :: form_scaer1=""
  CHARACTER(LEN=32) :: form_sccld1=""
  CHARACTER(LEN=32) :: form_pccoef1=""

  TYPE(rttov_coefs) :: coefs_all_chan

  CHARACTER(LEN=200) :: zeemanmsg
  REAL(KIND=jprb)    :: ZHOOK_HANDLE

#ifdef _RTTOV_HDF
  LOGICAL(KIND=jplm) :: hdf_was_open, hdf_was_64bit_reals
#endif
!- End of header --------------------------------------------------------

  TRY


  IF (LHOOK) CALL DR_HOOK('RTTOV_READ_COEFS', 0_jpim, ZHOOK_HANDLE)

  CALL rttov_nullify_coefs(coefs)
  IF (PRESENT(channels)) CALL rttov_nullify_coefs(coefs_all_chan)

  ! Get current status of HDF5 library
#ifdef _RTTOV_HDF
  hdf_was_open = is_hdf_open
  hdf_was_64bit_reals = is_hdf_64bit_reals
  CALL open_hdf(.TRUE._jplm, err)
  THROW(err.NE.0)
#endif

  !--------------------------
  ! Read rtcoef file
  !--------------------------
  CALL getlun(err, file_id_coef1, file_coef1, form_coef1, .FALSE._jplm, "rtcoef", &
                   file_id_coef, file_coef, form_coef, instrument )
  THROW(err.NE.0)

  IF (rttov_cmpuc(form_coef1,"unformatted")) THEN

    CALL rttov_read_binary_coef(err, coefs%coef, file_id_coef1, channels)
    THROW(err.NE.0)
    CALL closelun(err, file_id_coef1, file_id_coef)
    THROW(err.NE.0)

  ELSE IF (rttov_cmpuc(form_coef1,"formatted")) THEN

    CALL rttov_read_ascii_coef(err, coefs%coef, file_id_coef1, channels)
    THROW(err.NE.0)
    CALL closelun(err, file_id_coef1, file_id_coef)
    THROW(err.NE.0)

  ELSE IF (rttov_cmpuc(form_coef1,"hdf5")) THEN
#ifndef _RTTOV_HDF
    err =errorstatus_fatal
    THROWM(err.NE.0,"This program is not compiled with HDF5 capability; use RTTOV_HDF=1 with Makefile.PL")
#else
    IF (PRESENT(channels)) THEN

      CALL rttov_hdf_load(err, file_coef1, "/COEF", coef=coefs_all_chan%coef)
      THROW(err.NE.0)

      CALL rttov_channel_extract_coef(err, coefs_all_chan%coef, coefs%coef, channels)
      THROW(err.NE.0)

    ELSE

      CALL rttov_hdf_load(err, file_coef1, "/COEF", coef=coefs%coef)
      THROW(err.NE.0)

    ENDIF
#endif
  ELSE
    err = errorstatus_fatal
    THROWM(err.NE.0,"Unknown format "//TRIM(form_coef1))
  ENDIF
  file_id_coef1 = 0_jpim

  IF (coefs%coef%fmv_chn <= 0_jpim ) err = errorstatus_fatal
  THROWM(err.NE.0,"File does not contain RTTOV coefficients")

  IF (coefs%coef%inczeeman .AND. opts%config%verbose) THEN
    WRITE (zeemanmsg,"(a,i3,a)") "Zeeman coefficient file for "// &
      TRIM(inst_name(coefs%coef%id_inst))//" with ", coefs%coef%fmv_var(1), " predictors"
    INFO(zeemanmsg)
  ENDIF


  !--------------------------
  ! Read scaer file
  !--------------------------
  IF (opts%rt_ir%addaerosl .AND. .NOT. opts%rt_ir%user_aer_opt_param) THEN

    CALL getlun(err, file_id_scaer1, file_scaer1, form_scaer1, .FALSE._jplm, "scaercoef", &
                     file_id_scaer,  file_scaer, form_scaer, instrument )
    THROW(err.NE.0)

    IF (rttov_cmpuc(form_scaer1,"unformatted")) THEN

      CALL rttov_read_binary_scaercoef(err, coefs%coef, coefs%coef_scatt_ir, &
                                       coefs%optp, file_id_scaer1, channels)
      THROW(err.NE.0)
      CALL closelun(err, file_id_scaer1, file_id_scaer)
      THROW(err.NE.0)

    ELSE IF (rttov_cmpuc(form_scaer1,"formatted")) THEN

      CALL rttov_read_ascii_scaercoef(err, coefs%coef, coefs%coef_scatt_ir, &
                                      coefs%optp, file_id_scaer1, channels)
      THROW(err.NE.0)
      CALL closelun(err, file_id_scaer1, file_id_scaer)
      THROW(err.NE.0)

    ELSE IF (rttov_cmpuc(form_scaer1,"hdf5")) THEN
#ifndef _RTTOV_HDF
    err =errorstatus_fatal
    THROWM(err.NE.0,"This program is not compiled with HDF5 capability; use RTTOV_HDF=1 with Makefile.PL")
#else

      IF (PRESENT(channels)) THEN

        CALL rttov_hdf_load(err, file_scaer1, "/SCAER", &
                    scaercoef=coefs_all_chan%coef_scatt_ir, &
                    optp=coefs_all_chan%optp, coef=coefs_all_chan%coef)
        THROW(err.NE.0)

        CALL rttov_channel_extract_scaercoef(err, coefs_all_chan%coef_scatt_ir, coefs_all_chan%optp, &
                                             coefs%coef_scatt_ir, coefs%optp, channels)
        THROW(err.NE.0)

        CALL rttov_dealloc_coef_scatt_ir(err, coefs_all_chan%coef_scatt_ir)
        THROW(err.NE.0)

        CALL rttov_dealloc_optpar_ir(err, coefs_all_chan%optp)
        THROW(err.NE.0)

      ELSE

        CALL rttov_hdf_load(err, file_scaer1, "/SCAER", &
                    scaercoef=coefs%coef_scatt_ir, &
                    optp=coefs%optp, coef=coefs%coef)
        THROW(err.NE.0)

      ENDIF

#endif
    ELSE
      err = errorstatus_fatal
      THROWM(err.NE.0,"Unknown format "//TRIM(form_scaer1))
    ENDIF

  ENDIF
  file_id_scaer1 = 0_jpim


  !--------------------------
  ! Read sccld file
  !--------------------------
  IF (opts%rt_ir%addclouds .AND. .NOT. opts%rt_ir%user_cld_opt_param) THEN

    CALL getlun(err, file_id_sccld1, file_sccld1, form_sccld1, .FALSE._jplm, "sccldcoef", &
                     file_id_sccld, file_sccld, form_sccld, instrument )
    THROW(err.NE.0)

    IF (rttov_cmpuc(form_sccld1,"unformatted")) THEN

      CALL rttov_read_binary_sccldcoef(err, coefs%coef, coefs%coef_scatt_ir, &
                                       coefs%optp, file_id_sccld1, channels)
      THROW(err.NE.0)
      CALL closelun(err, file_id_sccld1, file_id_sccld)
      THROW(err.NE.0)

    ELSE IF (rttov_cmpuc(form_sccld1,"formatted")) THEN

      CALL rttov_read_ascii_sccldcoef(err, coefs%coef, coefs%coef_scatt_ir, &
                                      coefs%optp, file_id_sccld1, channels)
      THROW(err.NE.0)
      CALL closelun(err, file_id_sccld1, file_id_sccld)
      THROW(err.NE.0)

    ELSE IF (rttov_cmpuc(form_sccld1,"hdf5")) THEN
#ifndef _RTTOV_HDF
    err = errorstatus_fatal
    THROWM(err.NE.0,"This program is not compiled with HDF5 capability; use RTTOV_HDF=1 with Makefile.PL")
#else

      IF (PRESENT(channels)) THEN

        CALL rttov_hdf_load(err, file_sccld1, "/SCCLD", &
                    sccldcoef=coefs_all_chan%coef_scatt_ir, &
                    optp=coefs_all_chan%optp, coef=coefs_all_chan%coef)
        THROW(err.NE.0)

        CALL rttov_channel_extract_sccldcoef(err, coefs_all_chan%coef_scatt_ir, coefs_all_chan%optp, &
                                             coefs%coef_scatt_ir, coefs%optp, channels)
        THROW(err.NE.0)

        CALL rttov_dealloc_coef_scatt_ir(err, coefs_all_chan%coef_scatt_ir)
        THROW(err.NE.0)

        CALL rttov_dealloc_optpar_ir(err, coefs_all_chan%optp)
        THROW(err.NE.0)

      ELSE

        CALL rttov_hdf_load(err, file_sccld1, "/SCCLD", &
                    sccldcoef=coefs%coef_scatt_ir, &
                    optp=coefs%optp, coef=coefs%coef)
        THROW(err.NE.0)

      ENDIF

#endif
    ELSE
      err = errorstatus_fatal
      THROWM(err.NE.0,"Unknown format "//TRIM(form_sccld1))
    ENDIF

  ENDIF
  file_id_sccld1 = 0_jpim


  !--------------------------
  ! Read pccoef file
  !--------------------------
  IF (opts%rt_ir%pc%addpc) THEN

    IF (opts%rt_ir%addclouds) THEN
      CALL getlun(err, file_id_pccoef1, file_pccoef1, form_pccoef1, .FALSE._jplm, "pccoefcld", &
                       file_id_pccoef, file_pccoef, form_pccoef, instrument )
    ELSE
      CALL getlun(err, file_id_pccoef1, file_pccoef1, form_pccoef1, .FALSE._jplm, "pccoef", &
                       file_id_pccoef, file_pccoef, form_pccoef, instrument )
    ENDIF
    THROW(err.NE.0)

    IF (rttov_cmpuc(form_pccoef1,"unformatted")) THEN

      IF (opts%rt_ir%pc%addradrec) THEN
        CALL rttov_read_binary_pccoef(err, opts, coefs%coef, coefs%coef_pccomp, &
                                      file_id_pccoef1, channels, channels_rec)
        THROW(err.NE.0)
      ELSE
        CALL rttov_read_binary_pccoef(err, opts, coefs%coef, coefs%coef_pccomp, &
                                      file_id_pccoef1, channels)
        THROW(err.NE.0)
      ENDIF
      CALL closelun(err, file_id_pccoef1, file_id_pccoef)
      THROW(err.NE.0)

    ELSE IF (rttov_cmpuc(form_pccoef1,"formatted")) THEN

      IF (opts%rt_ir%pc%addradrec) THEN
        CALL rttov_read_ascii_pccoef(err, opts, coefs%coef, coefs%coef_pccomp, &
                                     file_id_pccoef1, channels, channels_rec)
        THROW(err.NE.0)
      ELSE
        CALL rttov_read_ascii_pccoef(err, opts, coefs%coef, coefs%coef_pccomp, &
                                     file_id_pccoef1, channels)
        THROW(err.NE.0)
      ENDIF
      CALL closelun(err, file_id_pccoef1, file_id_pccoef)
      THROW(err.NE.0)

    ELSE IF (rttov_cmpuc(form_pccoef1,"hdf5")) THEN
#ifndef _RTTOV_HDF
      err =errorstatus_fatal
      THROWM(err.NE.0,"This program is not compiled with HDF5 capability; use RTTOV_HDF=1 with Makefile.PL")
#else

      IF (PRESENT(channels)) THEN

        CALL rttov_hdf_load(err, file_pccoef1, "/PC", pccoef=coefs_all_chan%coef_pccomp, coef=coefs_all_chan%coef)
        THROW(err.NE.0)

        IF (opts%rt_ir%pc%addradrec) THEN
          CALL rttov_channel_extract_pccoef(err, coefs_all_chan%coef_pccomp, &
                                            coefs%coef_pccomp, channels, channels_rec)
        ELSE
          CALL rttov_channel_extract_pccoef(err, coefs_all_chan%coef_pccomp, &
                                            coefs%coef_pccomp, channels)
        ENDIF
        THROW(err.NE.0)

        CALL rttov_dealloc_coef_pccomp(err, coefs_all_chan%coef_pccomp)
        THROW(err.NE.0)

      ELSE

        CALL rttov_hdf_load(err, file_pccoef1, "/PC", pccoef=coefs%coef_pccomp, coef=coefs%coef)
        THROW(err.NE.0)

      ENDIF

#endif
    ELSE
      err = errorstatus_fatal
      THROWM(err.NE.0,"Unknown format "//TRIM(form_pccoef1))
    ENDIF

    IF (PRESENT(channels)) THEN
      CALL rttov_check_channels_pc(err, opts, coefs, channels)
      THROW(err.NE.0)
    ENDIF

  ENDIF
  file_id_pccoef1 = 0_jpim


#ifdef _RTTOV_HDF
  ! If HDF5 lib was open before rttov_read_coefs was called, make sure the real
  ! kind is the same as it was previously. Otherwise close the library.
  IF (hdf_was_open) THEN
    CALL open_hdf(hdf_was_64bit_reals, err)
    THROW(err.NE.0)
  ELSE
    CALL close_hdf(err)
    THROW(err.NE.0)
  ENDIF
#endif

  IF (PRESENT(channels)) THEN
    ! We could call dealloc_coefs here for all sub-structures in coefs_all_chan
    ! but it's more efficient in memory to deallocate them as soon as possible above.
    ! coefs_all_chan%coef is required for reading other coef types so must be
    ! deallocated last.
    CALL rttov_dealloc_coef(err, coefs_all_chan%coef)
    THROW(err.NE.0)
  ENDIF

  !--------------------------
  ! Finalise the coefficients
  !--------------------------
  CALL rttov_init_coefs(err, opts, coefs)
  THROWM(err.NE.0,"Failed to initialise coefficients")

  IF (LHOOK) CALL DR_HOOK('RTTOV_READ_COEFS', 1_jpim, ZHOOK_HANDLE)

  RETURN
  CATCH_L(999)
  CALL cleanup(coefs,             &
               form_coef1,        &
               form_scaer1,       &
               form_sccld1,       &
               form_pccoef1,      &
               file_id_coef1,     &
               file_id_scaer1,    &
               file_id_sccld1,    &
               file_id_pccoef1,   &
               file_id_coef,      &
               file_id_scaer,     &
               file_id_sccld,     &
               file_id_pccoef)
    err = errorstatus_fatal  ! err is reset by the try statement in cleanup routine

  IF (LHOOK) CALL DR_HOOK('RTTOV_READ_COEFS', 1_jpim, ZHOOK_HANDLE)

CONTAINS
  SUBROUTINE cleanup(coefs,             &
                     form_coef1,        &
                     form_scaer1,       &
                     form_sccld1,       &
                     form_pccoef1,      &
                     file_id_coef1,     &
                     file_id_scaer1,    &
                     file_id_sccld1,    &
                     file_id_pccoef1,   &
                     file_id_coef,      &
                     file_id_scaer,     &
                     file_id_sccld,     &
                     file_id_pccoef)
    IMPLICIT NONE

    TYPE(rttov_coefs),  INTENT(INOUT)        :: coefs

    CHARACTER(LEN=*),   INTENT(IN)           :: form_coef1
    CHARACTER(LEN=*),   INTENT(IN)           :: form_scaer1
    CHARACTER(LEN=*),   INTENT(IN)           :: form_sccld1
    CHARACTER(LEN=*),   INTENT(IN)           :: form_pccoef1
    INTEGER(KIND=jpim), INTENT(IN)           :: file_id_coef1
    INTEGER(KIND=jpim), INTENT(IN)           :: file_id_scaer1
    INTEGER(KIND=jpim), INTENT(IN)           :: file_id_sccld1
    INTEGER(KIND=jpim), INTENT(IN)           :: file_id_pccoef1
    INTEGER(KIND=jpim), INTENT(IN), OPTIONAL :: file_id_coef
    INTEGER(KIND=jpim), INTENT(IN), OPTIONAL :: file_id_scaer
    INTEGER(KIND=jpim), INTENT(IN), OPTIONAL :: file_id_sccld
    INTEGER(KIND=jpim), INTENT(IN), OPTIONAL :: file_id_pccoef

    TRY

    IF (file_id_coef1 > 0) THEN
      IF (.NOT. rttov_cmpuc(form_coef1, "hdf5")) THEN
        CALL closelun(err, file_id_coef1,   file_id_coef)
        THROW(err.NE.0)
      ENDIF
    ENDIF

    IF (file_id_scaer1 > 0) THEN
      IF (.NOT. rttov_cmpuc(form_scaer1, "hdf5")) THEN
        CALL closelun(err, file_id_scaer1,  file_id_scaer)
        THROW(err.NE.0)
      ENDIF
    ENDIF

    IF (file_id_sccld1 > 0) THEN
      IF (.NOT. rttov_cmpuc(form_sccld1, "hdf5")) THEN
        CALL closelun(err, file_id_sccld1,  file_id_sccld)
        THROW(err.NE.0)
      ENDIF
    ENDIF

    IF (file_id_pccoef1 > 0) THEN
      IF (.NOT. rttov_cmpuc(form_pccoef1, "hdf5")) THEN
        CALL closelun(err, file_id_pccoef1, file_id_pccoef)
        THROW(err.NE.0)
      ENDIF
    ENDIF

    CALL rttov_dealloc_coefs(err, coefs)
    THROWM(err.NE.0, "error deallocating coefficients")

    CATCH
  END SUBROUTINE
END SUBROUTINE
