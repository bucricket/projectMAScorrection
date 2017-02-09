subroutine rttov_write_coefs( err, coefs, opts, &
                             form_coef,     &
                             form_scaer,    &
                             form_sccld,    &
                             form_pccoef,   &
                             file_coef,     &
                             file_scaer,    &
                             file_sccld,    &
                             file_pccoef,   &
                             file_id_coef,  &
                             file_id_scaer, &
                             file_id_sccld, &
                             file_id_pccoef,&
                             instrument, &
                             compress)
! Description:
!   Write the coefficients structure
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
!    Copyright 2010, EUMETSAT, All Rights Reserved.
!
! Method:
!
! Current Code Owner: SAF NWP
!
! History:
! Version   Date     Comment
! -------   ----     -------
! 2010/03 Code cleaning Pascal Brunel, Philippe Marguinaud
!
!
! Code Description:
!   Language:           Fortran 90.
!   Software Standards: "European Standards for Writing and
!     Documenting Exchangeable Fortran 90 Code".
#include "throw.h"

Use rttov_types, Only : rttov_coefs, rttov_options
Use parkind1, Only : jpim
!INTF_OFF
Use parkind1, Only : jprb, jplm
Use rttov_coef_io_mod, only : getlun, closelun
Use yomhook, only : LHOOK, DR_HOOK
!INTF_ON

Implicit None

Integer(Kind=jpim), Intent(out) :: err
Type(rttov_coefs),  Intent(in)  :: coefs
Type(rttov_options),Intent(in)  :: opts
Character(len=*),   Intent(in), Optional :: form_coef
Character(len=*),   Intent(in), Optional :: form_scaer
Character(len=*),   Intent(in), Optional :: form_sccld
Character(len=*),   Intent(in), Optional :: form_pccoef
Character(len=*),   Intent(in), Optional :: file_coef
Character(len=*),   Intent(in), Optional :: file_scaer
Character(len=*),   Intent(in), Optional :: file_sccld
Character(len=*),   Intent(in), Optional :: file_pccoef
Integer(Kind=jpim), Intent(in), Optional :: file_id_coef
Integer(Kind=jpim), Intent(in), Optional :: file_id_scaer
Integer(Kind=jpim), Intent(in), Optional :: file_id_sccld
Integer(Kind=jpim), Intent(in), Optional :: file_id_pccoef
Integer(Kind=jpim), Intent(in), Optional :: instrument(3)
Logical,            Intent(in), Optional :: compress

!INTF_END

#include "rttov_errorreport.interface"
#include "rttov_cmpuc.interface"


#include "rttov_write_ascii_coef.interface"
#include "rttov_write_ascii_pccoef.interface"
#include "rttov_write_ascii_scaercoef.interface"
#include "rttov_write_ascii_sccldcoef.interface"
#include "rttov_write_binary_coef.interface"
#include "rttov_write_binary_pccoef.interface"
#include "rttov_write_binary_scaercoef.interface"
#include "rttov_write_binary_sccldcoef.interface"

#ifdef _RTTOV_HDF
#include "rttov_hdf_save.interface"
#endif

Integer(Kind=jpim) :: file_id_coef1
Integer(Kind=jpim) :: file_id_scaer1
Integer(Kind=jpim) :: file_id_sccld1
Integer(Kind=jpim) :: file_id_pccoef1

Character(len=256) :: file_coef1
Character(len=256) :: file_scaer1
Character(len=256) :: file_sccld1
Character(len=256) :: file_pccoef1

character(len=32) :: form_coef1
character(len=32) :: form_scaer1
character(len=32) :: form_sccld1
character(len=32) :: form_pccoef1

Logical :: CREATE
Real(Kind=jprb)   :: ZHOOK_HANDLE

TRY

!
IF (LHOOK) CALL DR_HOOK('RTTOV_WRITE_COEFS', 0_jpim, ZHOOK_HANDLE)

call getlun(err, file_id_coef1, file_coef1, form_coef1, .true._jplm, "rtcoef", &
               & file_id_coef, file_coef, form_coef, instrument )
THROW(err.ne.0)

if (rttov_cmpuc(form_coef1,"unformatted")) then
  call rttov_write_binary_coef (err, coefs%coef, file_id_coef1, opts%config%verbose)
  THROW(err.ne.0)
  call closelun(err, file_id_coef1, file_id_coef)
  THROW(err.ne.0)
else if (rttov_cmpuc(form_coef1,"formatted")) then
  call rttov_write_ascii_coef (err, coefs%coef, file_id_coef1, opts%config%verbose)
  THROW(err.ne.0)
  call closelun(err, file_id_coef1, file_id_coef)
  THROW(err.ne.0)
else if (rttov_cmpuc(form_coef1,"hdf5")) then
#ifndef _RTTOV_HDF
  err =errorstatus_fatal
  THROWM(err.ne.0,"This program is not compiled with HDF5 capability; use RTTOV_HDF=1 with Makefile.PL")
#else
  CALL RTTOV_HDF_SAVE( ERR, file_coef1, "/COEF", CREATE=.TRUE., COEF=coefs%coef, COMPRESS=compress)
  THROW(err.ne.0)
#endif
else
  err = errorstatus_fatal
  THROWM(err.ne.0,"Unknown format "//Trim(form_coef1))
endif



if(opts%rt_ir%addaerosl .and. .not. opts%rt_ir%user_aer_opt_param) then

  call getlun(err, file_id_scaer1, file_scaer1, form_scaer1, .true._jplm, "scaer", &
                 & file_id_scaer, file_scaer, form_scaer, instrument )
  THROW(err.ne.0)
  if (rttov_cmpuc(form_scaer1,"unformatted")) then
    call rttov_write_binary_scaercoef (err, coefs%coef, coefs%coef_scatt_ir, &
                                      coefs%optp, file_id_scaer1, opts%config%verbose)
    THROW(err.ne.0)
    call closelun(err, file_id_scaer1, file_id_scaer)
    THROW(err.ne.0)
  else if (rttov_cmpuc(form_scaer1,"formatted")) then
    call rttov_write_ascii_scaercoef (err, coefs%coef, coefs%coef_scatt_ir, &
                                     coefs%optp, file_id_scaer1, opts%config%verbose)
    THROW(err.ne.0)
    call closelun(err, file_id_scaer1, file_id_scaer)
    THROW(err.ne.0)
  else if (rttov_cmpuc(form_coef1,"hdf5")) then
#ifndef _RTTOV_HDF
    err =errorstatus_fatal
    THROWM(err.ne.0,"This program is not compiled with HDF5 capability; use RTTOV_HDF=1 with Makefile.PL")
#else
    if( TRIM(file_coef) .EQ. TRIM(file_scaer) .AND. TRIM(file_coef) .NE. "") THEN
      CREATE=.FALSE.
    else
      CREATE=.TRUE.
    endif
    CALL RTTOV_HDF_SAVE( ERR, file_scaer, "/SCAER", CREATE=CREATE, &
              & SCAERCOEF = COEFS%COEF_SCATT_IR, &                            
              & OPTP=COEFS%OPTP, COMPRESS=compress)
    THROW(err.ne.0)
#endif
  else
    err = errorstatus_fatal
    THROWM(err.ne.0,"Unknown format "//Trim(form_scaer1))
  endif
endif

!

if(opts%rt_ir%addclouds .and. .not. opts%rt_ir%user_cld_opt_param) then
  call getlun(err, file_id_sccld1, file_sccld1, form_sccld1, .true._jplm, "sccld", &
                 & file_id_sccld, file_sccld, form_sccld, instrument)
  THROW(err.ne.0)
  if (rttov_cmpuc(form_sccld1,"unformatted")) then
    call rttov_write_binary_sccldcoef (err, coefs%coef, coefs%coef_scatt_ir, &                            
                                      coefs%optp, file_id_sccld1, opts%config%verbose)
    THROW(err.ne.0)
    call closelun(err, file_id_sccld1, file_id_sccld)
    THROW(err.ne.0)
  else if (rttov_cmpuc(form_sccld1,"formatted")) then
    call rttov_write_ascii_sccldcoef (err, coefs%coef, coefs%coef_scatt_ir, &                             
                                     coefs%optp, file_id_sccld1, opts%config%verbose)
    THROW(err.ne.0)   
    call closelun(err, file_id_sccld1, file_id_sccld)
    THROW(err.ne.0)
else if (rttov_cmpuc(form_coef1,"hdf5")) then
#ifndef _RTTOV_HDF
    err =errorstatus_fatal
    THROWM(err.ne.0,"This program is not compiled with HDF5 capability; use RTTOV_HDF=1 with Makefile.PL")
#else 
    if( TRIM(file_coef) .EQ. TRIM(file_sccld) .AND. TRIM(file_coef) .NE. "") THEN
      CREATE=.FALSE.
    else
      CREATE=.TRUE.
    endif
    CALL RTTOV_HDF_SAVE( ERR, file_sccld, "/SCCLD", CREATE=CREATE, &
              & SCCLDCOEF = COEFS%COEF_SCATT_IR, &                            
              & OPTP=COEFS%OPTP, COMPRESS=compress)
    THROW(err.ne.0)
#endif                                                                                   
  else
    err = errorstatus_fatal
    THROWM(err.ne.0,"Unknown format "//Trim(form_sccld1))
  endif
endif
!

if(opts%rt_ir%pc%addpc) then
  call getlun(err, file_id_pccoef1, file_pccoef1, form_pccoef1, .true._jplm, "pccoef", &
                 & file_id_pccoef, file_pccoef, form_pccoef, instrument)
  THROW(err.ne.0)
  if (rttov_cmpuc(form_pccoef1,"unformatted")) then
    call rttov_write_binary_pccoef (err, coefs%coef_pccomp, file_id_pccoef1, opts%config%verbose)
    THROW(err.ne.0)
    call closelun(err, file_id_pccoef1, file_id_pccoef)
    THROW(err.ne.0)
  else if (rttov_cmpuc(form_pccoef1,"formatted")) then
    call rttov_write_ascii_pccoef (err,  coefs%coef_pccomp, file_id_pccoef1, opts%config%verbose)
    THROW(err.ne.0) 
    call closelun(err, file_id_pccoef1, file_id_pccoef)
    THROW(err.ne.0)
  else if (rttov_cmpuc(form_coef1,"hdf5")) then
#ifndef _RTTOV_HDF
    err =errorstatus_fatal
    THROWM(err.ne.0,"This program is not compiled with HDF5 capability; use RTTOV_HDF=1 with Makefile.PL")
#else
    if( TRIM(file_coef) .EQ. TRIM(file_pccoef) .AND. TRIM(file_coef) .NE. "") THEN
      CREATE=.FALSE.
    else
      CREATE=.TRUE.
    endif
    CALL RTTOV_HDF_SAVE( ERR, file_pccoef, "/PC", CREATE=CREATE, &
              & PCCOEF = COEFS%COEF_PCCOMP, COMPRESS=compress)
    THROW(err.ne.0)
#endif   
  else
    err = errorstatus_fatal
    THROWM(err.ne.0,"Unknown format "//Trim(form_pccoef1))
  endif
endif

IF (LHOOK) CALL DR_HOOK('RTTOV_WRITE_COEFS', 1_jpim, ZHOOK_HANDLE)
CATCH
IF (LHOOK) CALL DR_HOOK('RTTOV_WRITE_COEFS', 1_jpim, ZHOOK_HANDLE)

end subroutine
