Program rttov_test

  !
  !
  ! Copyright:
  !
  !    This software was developed within the context of
  !    the EUMETSAT Satellite Application Facility on
  !    Numerical Weather Prediction (NWP SAF), under the
  !    Cooperation Agreement dated 25 November 1998, between
  !    EUMETSAT and the Met Office, UK, by one or more partners
  !    within the NWP SAF. The partners in the NWP SAF are
  !    the Met Office, ECMWF, KNMI and MeteoFrance.
  !
  !    Copyright 2007, EUMETSAT, All Rights Reserved.
  !
  !     *************************************************************
  !
  !
  ! Current Code Owner: SAF NWP
  !
  ! History:
  ! Version   Date     Comment
  ! -------   ----     -------
  !          06/2008    Creation P. Marguinaud & P. Brunel
  !
  ! 2010/03 Code cleaning Pascal Brunel, Philippe Marguinaud
  !

  Use rttov_const, Only :      &
       & errorstatus_success,  &
       & errorstatus_fatal,    &
       & sensor_id_mw,         &
       & sensor_id_po,         &
       & sensor_id_ir,         &
       & sensor_id_hi,         &
       & surftype_sea,         &
       & pi_r,                 &
       & gas_id_watervapour,   &
       & gas_id_ozone,         &
       & gas_id_co2,           &
       & gas_id_n2o,           &
       & gas_id_co,            &
       & gas_id_ch4,           &
       & mair,                 &
       & gas_mass,             &
       & gas_unit_specconc,    &
       & gas_unit_ppmv,        &
       & gas_unit_compatibility

  Use rttov_types, Only :     &
       & rttov_coefs,         &
       & rttov_options,       &
       & rttov_pccomp,        &
       & profile_Type,        &
       & transmission_Type,   &
       & radiance_Type,       &
       & radiance2_Type,      &
       & rttov_chanprof,      &
       & rttov_emissivity,    &
       & rttov_reflectance,   &
       & rttov_opt_param,     &
       & sskin_type,          &
       & s2m_type,            &
       & rttov_traj

#include "throw.h"

  Use parkind1, Only : jpim, jprb, jplm

  Use rttov_getoptions, Only : getoption
  Use yomhook, Only : lhook, dr_hook
  Use rttov_lun
  Use rttov_unix_env, Only: rttov_exit

#ifdef _RTTOV_HDF
  Use hdf5
  Use rttov_hdf_mod
#endif

  Implicit None



#include "rttov_direct.interface"
#include "rttov_tl.interface"
#include "rttov_ad.interface"
#include "rttov_k.interface"
#include "rttov_k_tl.interface"
#include "rttov_k_ad.interface"
#include "rttov_k_bf.interface"


#include "rttov_parallel_direct.interface"
#include "rttov_parallel_tl.interface"
#include "rttov_parallel_ad.interface"
#include "rttov_parallel_k.interface"

#include "rttov_make_profile_inc.interface"
#include "rttov_scale_profile_inc.interface"
#include "rttov_make_radiance_inc.interface"
#include "rttov_scale_radiance_inc.interface"
#include "rttov_make_pccomp_inc.interface"
#include "rttov_scale_pccomp_inc.interface"

#include "rttov_read_coefs.interface"
#include "rttov_dealloc_coefs.interface"
#include "rttov_alloc_rad.interface"
#include "rttov_alloc_transmission.interface"
#include "rttov_alloc_prof.interface"
#include "rttov_alloc_pccomp.interface"
#include "rttov_init_rad.interface"
#include "rttov_init_prof.interface"
#include "rttov_add_prof.interface"
#include "rttov_copy_prof.interface"
#include "rttov_copy_rad.interface"
#include "rttov_copy_pccomp.interface"
#include "rttov_alloc_traj.interface"
#include "rttov_alloc_opt_param.interface"
#include "rttov_init_opt_param.interface"

#include "rttov_user_options_checkinput.interface"
#include "rttov_user_profile_checkinput.interface"

#include "rttov_errorreport.interface"

#ifdef _RTTOV_HDF
#include "rttov_hdf_save.interface"
#endif

  Character(len=*), Parameter :: Gformat = 'G15.6E3'


! Test definition; we have here all dimensions, flags, etc...
  Type rttov_test_defn
! dimensions
    Integer(Kind=jpim) :: nlevels = -1_jpim, nchannels = -1_jpim, nprofiles = -1_jpim, mult = 1_jpim
    Integer(Kind=jpim) :: nchannels_rec = -1_jpim, npcscores = -1_jpim, prof_ndigits = -1_jpim

! coefficients files
    Character(len=256) :: f_coef = "", f_scat = "", f_optp = "", f_pccomp = ""
    Character(len=256) :: coef_prefix

! flags

    Type(rttov_options) :: opts

! input files
    Character(len=256) :: f_p    = "", f_t   = "", f_q      = "", f_o3      = "", f_co2    = "", f_co    = "",  &
                          f_ch4  = "", f_clw = "", f_aerosl = "", f_cloud   = "", f_cfrac  = "", f_n2o   = "",  &
                          f_skin = "", f_s2m = "", f_angles = "", f_aerosli = "", f_cloud0 = "", f_flags = "",  &
                          f_icede = "", f_be   = "", f_elevation = "", f_datetime = "", f_gas_units = ""
    Character(len=256) :: f_channels = "", f_pmc   = "", f_emissivity = "", f_calcemis = "", &
                          f_reflectance = "", f_calcrefl = "", f_reflectance_cloud_top = "", &
                          f_lprofiles = "", f_pcscores = "", f_channels_rec = ""
    Character(len=256) :: f_aer_opt_param = "", f_cld_opt_param = ""

! what we will do
   Logical(Kind=jplm) :: do_direct = .FALSE., do_tl   = .FALSE., do_ad   = .FALSE., do_k = .FALSE., &
                         do_k_bf   = .FALSE., do_k_tl = .FALSE., do_k_ad = .FALSE.,                 &
                         do_taylor = .FALSE., taylor_by_chan = .FALSE.
   Logical(Kind=jplm) :: do_print  = .FALSE., calc_rad2 = .FALSE., ltemp = .FALSE., lallchans = .FALSE.
   Logical(Kind=jplm) :: savehdf5 = .FALSE.
   Integer(Kind=jpim) :: ntimes = 1, nthreads = 0
   Logical(Kind=jplm) :: prof_by_prof = .FALSE. ! run RTTOV profile by profile
   Logical(Kind=jplm) :: chan_by_chan = .FALSE. ! run RTTOV channel by channel
   Logical(Kind=jplm) :: user_check_opts = .FALSE., user_check_prof = .FALSE.
   Integer(Kind=jpim) :: input_gas_units = -1           ! Gas units in input file
   Integer(Kind=jpim) :: run_gas_units = gas_unit_ppmv  ! Gas units RTTOV is run with
   Real(Kind=jprb)    :: scale_inc = 1.0_jprb
   Real(Kind=jprb)    :: scale_out = 1.0_jprb
  End Type

! The basic data structures passed to rttov_direct/rttov_tl/rttov_ad/rttov_k
  Type rttov_input_data
   Type(profile_type),      Pointer :: profiles(:)         => NULL()
   Type(profile_type),      Pointer :: profiles_pert(:)    => NULL()
   Type(rttov_emissivity),  Pointer :: emissivity(:)       => NULL()
   Type(rttov_emissivity),  Pointer :: emissivity_pert(:)  => NULL()
   Type(rttov_reflectance), Pointer :: reflectance(:)      => NULL()
   Type(rttov_reflectance), Pointer :: reflectance_pert(:) => NULL()
   Type(transmission_Type)          :: transmission
   Type(radiance_Type)              :: radiance
   Type(radiance_Type)              :: radiance_saved
   Type(radiance2_Type)             :: radiance2
   Type(rttov_pccomp)               :: pccomp
   Type(rttov_pccomp)               :: pccomp_saved
  End Type
 
  Type rttov_k_input_data
   Type(profile_type),   Pointer :: profiles_k_pc(:) => NULL()
   Type(profile_type),   Pointer :: profiles_k_rec(:) => NULL()
  End Type

! Test data. All dynamically allocated data goes here.
  Type rttov_test_data
   Type(rttov_coefs)             :: coefs
   Type(rttov_chanprof), Pointer :: chanprof(:)          => NULL()
   Type(rttov_opt_param)         :: aer_opt_param
   Type(rttov_opt_param)         :: cld_opt_param
   Integer(Kind=jpim)            :: rttov_errorstatus
   Type(rttov_input_data)        :: direct, tl, ad, k, k_bf, k_tl, k_ad
   Type(rttov_k_input_data)      :: data_k, data_k_bf, data_k_tl, data_k_ad
   Logical(Kind=jplm),   Pointer :: calcemis(:)          => NULL()
   Logical(Kind=jplm),   Pointer :: calcrefl(:)          => NULL()
   Type(rttov_traj)              :: traj, traj_tl, traj_ad, traj_k
   Integer(Kind=jpim),   Pointer :: channels_rec(:)      => NULL()
  End Type

  Type rttov_test_struct
    Integer(Kind=jpim)    :: rank
    Character(len=256)    :: path
    Type(rttov_test_defn) :: defn
    Type(rttov_test_data) :: data
    Integer(Kind=jpim)    :: err = 0_jpim
  End Type


  Type(rttov_test_struct), Pointer :: ts(:)     => NULL()
  Character(len=256),      Pointer :: paths(:)  => NULL()
  Integer(Kind=jpim) :: ip, np
  Real(Kind=jprb)    :: zhook_handle
  Integer(Kind=jpim) :: err, rttov_run_err

TRY

  If(lhook) Call dr_hook('rttov_test',0_jpim,zhook_handle)

#ifdef _RTTOV_HDF
  CALL OPEN_HDF( .FALSE., ERR ) ! 32 bits
  THROW(ERR.NE.0)
#endif


  Call GetOption( "--test-path-list", paths )

  If( .not. associated( paths ) ) Then
    Allocate( paths(1), stat = err )
    THROWM(err.ne.0,"Allocation of paths failed")
    paths(1) = "."
  EndIf

  np = Size(paths)

  Allocate( ts( np ), stat = err )
  THROWM(err.ne.0,"Allocation of ts failed")
  Do ip = 1, np
    ts%rank = ip
  EndDo


  Do ip = 1, np
    call setup( ts(ip), paths(ip), err )
    THROW(err.ne.0)
  EndDo
  Deallocate( paths, stat = err )
  THROWM(err.ne.0,"DeAllocation of paths failed")

  rttov_run_err = errorstatus_success

  Do ip = 1, np

    if(err.ne.0)cycle

    If( ts(ip)%defn%ltemp ) Then
      call run(                                   &
        ts(ip),                                   &
        ts(ip)%err,                               &
        ts(ip)%data%traj,    ts(ip)%data%traj_tl, &
        ts(ip)%data%traj_ad, ts(ip)%data%traj_k )
    Else
      call run( ts(ip), ts(ip)%err )
    EndIf

    THROW_L(ts(ip)%err.ne.0,800)

    cycle

    CATCH_L(800)

    if( ts(ip)%err .ne. 0 ) then
!$OMP FLUSH(err)
      rttov_run_err = errorstatus_fatal
    endif

  EndDo

  Do ip = 1, np
    call cleanup( ts(ip), err )
    THROW(err.ne.0)
  EndDo

  Deallocate( ts, stat = err )
  THROWM(err.ne.0,"Deallocation of ts")

  THROW(rttov_run_err.ne.0)

#ifdef _RTTOV_HDF
  CALL CLOSE_HDF( ERR )
  THROW(ERR.NE.0)
#endif

  If(lhook) Call dr_hook('rttov_test',1_jpim,zhook_handle)

PCATCH
  If(lhook) Call dr_hook('rttov_test',1_jpim,zhook_handle)

  Call rttov_exit(1_jpim)

  
  Contains

  Subroutine run( ts, err, traj, traj_tl, traj_ad, traj_k )
!
  Implicit None
  
  Type(rttov_test_struct),    Target, Intent(inout) :: ts
  Integer(Kind=jpim),                 Intent(out)   :: err
  Type(rttov_traj), Optional, Target, Intent(inout) :: traj, traj_tl, traj_ad, traj_k
!
  Integer(Kind=jpim), Parameter :: niter = 8
  Integer(Kind=jpim)            :: itime, iter, i, j, chan
  Integer(Kind=jpim)            :: nchanrec, npcscores
  Integer(Kind=jpim)            :: file_id
  Real(Kind=jprb)               :: ratio
  Character(Len=256)            :: taylor_filename
  Real(Kind=jprb), Allocatable  :: taylor_direct(:,:), taylor_tl(:,:)
  
!
TRY

  ! If we do the user checks here and throw any resulting errors then the cleanup routine
  ! will ensure all memory is deallocated
  If( ts%defn%user_check_opts ) Then
    Call rttov_user_options_checkinput(err, ts%defn%opts, ts%data%coefs)
    THROW(err.ne.0)
  EndIf
  If( ts%defn%user_check_prof ) Then
    Do i = 1, ts%defn%nprofiles
      Call rttov_user_profile_checkinput(err, ts%defn%opts, ts%data%coefs, &
          & ts%data%direct%profiles(i), ts%data%aer_opt_param, ts%data%cld_opt_param)
      THROW(err.ne.0)
    EndDo
  EndIf
  
  Do itime = 1, ts%defn%ntimes

!    verbose = ts%verbose
  
    If( ts%defn%do_direct ) Then

      If (ts%defn%calc_rad2) Then

        If( (ts%defn%nthreads .le. 0) .and. ts%defn%prof_by_prof ) Then

          call rttov_parallel_direct(                                         &
            & ts%data%rttov_errorstatus,       ts%data%chanprof,              &
            & ts%defn%opts,                    ts%data%direct%profiles,       &
            & ts%data%coefs,                   ts%data%direct%transmission,   &
            & ts%data%direct%radiance,         ts%data%direct%radiance2,      &
            & ts%data%calcemis,                ts%data%direct%emissivity,     &
            & ts%data%calcrefl,                ts%data%direct%reflectance,    &
            & ts%data%aer_opt_param,           ts%data%cld_opt_param,         &
            & nthreads = 1_jpim,               strategy = 1_jpim )

        Else If( (ts%defn%nthreads .le. 0) .and. ts%defn%chan_by_chan ) Then

          call rttov_parallel_direct(                                         &
            & ts%data%rttov_errorstatus,       ts%data%chanprof,              &
            & ts%defn%opts,                    ts%data%direct%profiles,       &
            & ts%data%coefs,                   ts%data%direct%transmission,   &
            & ts%data%direct%radiance,         ts%data%direct%radiance2,      &
            & ts%data%calcemis,                ts%data%direct%emissivity,     &
            & ts%data%calcrefl,                ts%data%direct%reflectance,    &
            & ts%data%aer_opt_param,           ts%data%cld_opt_param,         &
            & nthreads = 1_jpim,               strategy = 2_jpim )

        Else If( ts%defn%nthreads .le. 0 ) Then

          call rttov_direct(                       &
            & ts%data%rttov_errorstatus,           &
            & ts%data%chanprof,                    &
            & ts%defn%opts,                        &
            & ts%data%direct%profiles,             &
            & ts%data%coefs,                       &
            & ts%data%direct%transmission,         &
            & ts%data%direct%radiance,             &
            & ts%data%direct%radiance2,            &
            & ts%data%calcemis,                    &
            & ts%data%direct%emissivity,           &
            & ts%data%calcrefl,                    &
            & ts%data%direct%reflectance,          &
            & ts%data%aer_opt_param,               &
            & ts%data%cld_opt_param,               &
            & pccomp       = ts%data%direct%pccomp,&
            & channels_rec = ts%data%channels_rec, &
            & traj         = traj )

        Else

          call rttov_parallel_direct(                                         &
            & ts%data%rttov_errorstatus,       ts%data%chanprof,              &
            & ts%defn%opts,                    ts%data%direct%profiles,       &
            & ts%data%coefs,                   ts%data%direct%transmission,   &
            & ts%data%direct%radiance,         ts%data%direct%radiance2,      &
            & ts%data%calcemis,                ts%data%direct%emissivity,     &
            & ts%data%calcrefl,                ts%data%direct%reflectance,    &
            & ts%data%aer_opt_param,           ts%data%cld_opt_param,         &
            & pccomp       = ts%data%direct%pccomp,                           &
            & channels_rec = ts%data%channels_rec,                            &
            & nthreads     = ts%defn%nthreads)

        EndIf

      Else

        If( (ts%defn%nthreads .le. 0) .and. ts%defn%prof_by_prof ) Then

          call rttov_parallel_direct(                                         &
            & ts%data%rttov_errorstatus,       ts%data%chanprof,              &
            & ts%defn%opts,                    ts%data%direct%profiles,       &
            & ts%data%coefs,                   ts%data%direct%transmission,   &
            & ts%data%direct%radiance,                                        &
            & calcemis      = ts%data%calcemis,                               &
            & emissivity    = ts%data%direct%emissivity,                      &
            & calcrefl      = ts%data%calcrefl,                               &
            & reflectance   = ts%data%direct%reflectance,                     &
            & aer_opt_param = ts%data%aer_opt_param,                          &
            & cld_opt_param = ts%data%cld_opt_param,                          &
            & nthreads = 1_jpim,               strategy = 1_jpim )

        Else If( (ts%defn%nthreads .le. 0) .and. ts%defn%chan_by_chan ) Then

          call rttov_parallel_direct(                                         &
            & ts%data%rttov_errorstatus,       ts%data%chanprof,              &
            & ts%defn%opts,                    ts%data%direct%profiles,       &
            & ts%data%coefs,                   ts%data%direct%transmission,   &
            & ts%data%direct%radiance,                                        &
            & calcemis      = ts%data%calcemis,                               &
            & emissivity    = ts%data%direct%emissivity,                      &
            & calcrefl      = ts%data%calcrefl,                               &
            & reflectance   = ts%data%direct%reflectance,                     &
            & aer_opt_param = ts%data%aer_opt_param,                          &
            & cld_opt_param = ts%data%cld_opt_param,                          &
            & nthreads = 1_jpim,               strategy = 2_jpim )

        Else If( ts%defn%nthreads .le. 0 ) Then

          call rttov_direct(                                                  &
            & ts%data%rttov_errorstatus,                                      &
            & ts%data%chanprof,                                               &
            & ts%defn%opts,                                                   &
            & ts%data%direct%profiles,                                        &
            & ts%data%coefs,                                                  &
            & ts%data%direct%transmission,                                    &
            & ts%data%direct%radiance,                                        &
            & calcemis      = ts%data%calcemis,                               &
            & emissivity    = ts%data%direct%emissivity,                      &
            & calcrefl      = ts%data%calcrefl,                               &
            & reflectance   = ts%data%direct%reflectance,                     &
            & aer_opt_param = ts%data%aer_opt_param,                          &
            & cld_opt_param = ts%data%cld_opt_param,                          &
            & pccomp       = ts%data%direct%pccomp,                           &
            & channels_rec = ts%data%channels_rec,                            &
            & traj         = traj )

        Else

          call rttov_parallel_direct(                                         &
            & ts%data%rttov_errorstatus,       ts%data%chanprof,              &
            & ts%defn%opts,                    ts%data%direct%profiles,       &
            & ts%data%coefs,                   ts%data%direct%transmission,   &
            & ts%data%direct%radiance,                                        &
            & calcemis      = ts%data%calcemis,                               &
            & emissivity    = ts%data%direct%emissivity,                      &
            & calcrefl      = ts%data%calcrefl,                               &
            & reflectance   = ts%data%direct%reflectance,                     &
            & aer_opt_param = ts%data%aer_opt_param,                          &
            & cld_opt_param = ts%data%cld_opt_param,                          &
            & pccomp        = ts%data%direct%pccomp,                          &
            & channels_rec  = ts%data%channels_rec,                           &
            & nthreads      = ts%defn%nthreads)

        EndIf

      EndIf

      If( ts%data%rttov_errorstatus /= 0 ) Then
        err = errorstatus_fatal
        THROW(err.ne.0)
      EndIf
      
      If( ts%defn%do_taylor ) Then
        If( ts%defn%opts%rt_ir%pc%addpc ) Then
          Call rttov_copy_pccomp(ts%data%direct%pccomp_saved, ts%data%direct%pccomp)
        Else
          Call rttov_copy_rad(ts%data%direct%radiance_saved, ts%data%direct%radiance)
        EndIf
      EndIf

      If( ts%defn%do_print ) Then
        If( ts%defn%calc_rad2 ) Then
          Call PrintRadiance (err, ts%data%direct%radiance, "RADIANCE", Trim(ts%path)//"/direct/radiance.txt", &
                              ts%data%direct%radiance2 )
        Else
          Call PrintRadiance (err, ts%data%direct%radiance, "RADIANCE", Trim(ts%path)//"/direct/radiance.txt" )
        EndIf
        THROW(err.ne.0)
        Call PrintTransmission (err, ts%data%direct%transmission, "TRANSMISSION", &
          & Trim(ts%path)//"/direct/transmission.txt" )
        THROW(err.ne.0)
        Call PrintEmissivity (err, ts%data%direct%emissivity%emis_out, "EMISSIVITY_OUT", &
          & Trim(ts%path)//"/direct/emissivity_out.txt" )
        THROW(err.ne.0)
        If (ts%defn%opts%rt_ir%addsolar) Then
          Call PrintReflectance (err, ts%data%direct%reflectance%refl_out, "REFLECTANCE_OUT", &
            & Trim(ts%path)//"/direct/reflectance_out.txt" )
          THROW(err.ne.0)
        EndIf
        If( ts%defn%opts%rt_ir%pc%addpc ) Then
          Call PrintPCScores (err, ts%data%direct%pccomp, "PCSCORES", Trim(ts%path)//"/direct/pcscores.txt" )
          THROW(err.ne.0)
        EndIf
      EndIf


#ifdef _RTTOV_HDF
      If( ts%defn%savehdf5 ) Then
       CALL RTTOV_HDF_SAVE( ERR, TRIM(ts%path)//"/direct/options.H5",      '/OPTIONS', CREATE=.true., &
       &  OPTIONS = ts%defn%opts )
                 THROW(err.ne.0)

       CALL RTTOV_HDF_SAVE( ERR, TRIM(ts%path)//"/direct/chanprof.H5",      '/CHANPROF', CREATE=.true., &
       &  CHANPROF = ts%data%chanprof )
                 THROW(err.ne.0)

       CALL RTTOV_HDF_SAVE( ERR, TRIM(ts%path)//"/direct/radiance.H5",      '/RADIANCE', CREATE=.true., &
       &  RADIANCE = ts%data%direct%radiance )
                 THROW(err.ne.0)

       If( ts%defn%calc_rad2 ) Then
         CALL RTTOV_HDF_SAVE( ERR, TRIM(ts%path)//"/direct/radiance2.H5",     '/RADIANCE2', CREATE=.true., &
         &  RADIANCE2 = ts%data%direct%radiance2 )
                   THROW(err.ne.0)
       EndIf

       CALL RTTOV_HDF_SAVE( ERR, TRIM(ts%path)//"/direct/transmission.H5",  '/TRANSMISSION', CREATE=.true., &
       &  TRANSMISSION = ts%data%direct%transmission )
                 THROW(err.ne.0)

       CALL RTTOV_HDF_SAVE( ERR, TRIM(ts%path)//"/direct/emissivity.H5",     '/EMISSIVITY', CREATE=.true., &
       &  EMISSIVITY = ts%data%direct%emissivity)
                 THROW(err.ne.0)
       CALL RTTOV_HDF_SAVE( ERR, TRIM(ts%path)//"/direct/emissivity.H5",     '/EMISSIVITY', CREATE=.false., &
       &  l1 = ts%data%calcemis, SNAME='CALCEMIS')
       THROW(err.ne.0)
         
       If (ts%defn%opts%rt_ir%addsolar) Then
         CALL RTTOV_HDF_SAVE( ERR, TRIM(ts%path)//"/direct/reflectance.H5",  '/REFLECTANCE', CREATE=.true., &
         &  REFLECTANCE = ts%data%direct%reflectance)
                   THROW(err.ne.0)
         CALL RTTOV_HDF_SAVE( ERR, TRIM(ts%path)//"/direct/reflectance.H5",     '/REFLECTANCE', CREATE=.false., &
         &  l1 = ts%data%calcrefl, SNAME='CALCREFL')
         THROW(err.ne.0)
       EndIf
      EndIf 
#endif

      ts%defn%opts%config%verbose = .false.
    
    EndIf 
    
    If( ts%defn%do_tl ) Then
  
      If( (ts%defn%nthreads .le. 0) .and. ts%defn%prof_by_prof ) Then

        call rttov_parallel_tl(                                             &
          & ts%data%rttov_errorstatus,      ts%data%chanprof,               &
          & ts%defn%opts,                   ts%data%direct%profiles,        &
          & ts%data%tl%profiles,            ts%data%coefs,                  &
          & ts%data%direct%transmission,    ts%data%tl%transmission,        &
          & ts%data%direct%radiance,        ts%data%tl%radiance,            &
          & ts%data%calcemis,               ts%data%direct%emissivity,      &
          & ts%data%tl%emissivity,          ts%data%calcrefl,               &
          & ts%data%direct%reflectance,     ts%data%tl%reflectance,         &
          & ts%data%aer_opt_param,          ts%data%cld_opt_param,          &
          & nthreads = 1_jpim,              strategy = 1_jpim )

      Else If( (ts%defn%nthreads .le. 0) .and. ts%defn%chan_by_chan ) Then

        call rttov_parallel_tl(                                             &
          & ts%data%rttov_errorstatus,      ts%data%chanprof,               &
          & ts%defn%opts,                   ts%data%direct%profiles,        &
          & ts%data%tl%profiles,            ts%data%coefs,                  &
          & ts%data%direct%transmission,    ts%data%tl%transmission,        &
          & ts%data%direct%radiance,        ts%data%tl%radiance,            &
          & ts%data%calcemis,               ts%data%direct%emissivity,      &
          & ts%data%tl%emissivity,          ts%data%calcrefl,               &
          & ts%data%direct%reflectance,     ts%data%tl%reflectance,         &
          & ts%data%aer_opt_param,          ts%data%cld_opt_param,          &
          & nthreads = 1_jpim,              strategy = 2_jpim )
          
      Else If( ts%defn%nthreads .le. 0 ) Then
    
        call rttov_tl(                             &
          & ts%data%rttov_errorstatus,             &
          & ts%data%chanprof,                      &
          & ts%defn%opts,                          &
          & ts%data%direct%profiles,               &
          & ts%data%tl%profiles,                   &
          & ts%data%coefs,                         &
          & ts%data%direct%transmission,           &
          & ts%data%tl%transmission,               &
          & ts%data%direct%radiance,               &
          & ts%data%tl%radiance,                   &
          & ts%data%calcemis,                      &
          & ts%data%direct%emissivity,             &
          & ts%data%tl%emissivity,                 &
          & ts%data%calcrefl,                      &
          & ts%data%direct%reflectance,            &
          & ts%data%tl%reflectance,                &
          & ts%data%aer_opt_param,                 &
          & ts%data%cld_opt_param,                 &
          & channels_rec = ts%data%channels_rec,   &
          & pccomp       = ts%data%direct%pccomp,  &
          & pccomp_tl    = ts%data%tl%pccomp,      &
          & traj         = traj,                   &
          & traj_tl      = traj_tl )

      Else

        call rttov_parallel_tl(                                             &
          & ts%data%rttov_errorstatus,      ts%data%chanprof,               &
          & ts%defn%opts,                   ts%data%direct%profiles,        &
          & ts%data%tl%profiles,            ts%data%coefs,                  &
          & ts%data%direct%transmission,    ts%data%tl%transmission,        &
          & ts%data%direct%radiance,        ts%data%tl%radiance,            &
          & ts%data%calcemis,               ts%data%direct%emissivity,      &
          & ts%data%tl%emissivity,          ts%data%calcrefl,               &
          & ts%data%direct%reflectance,     ts%data%tl%reflectance,         &
          & ts%data%aer_opt_param,          ts%data%cld_opt_param,          &
          & channels_rec = ts%data%channels_rec,                            &
          & pccomp       = ts%data%direct%pccomp,                           &
          & pccomp_tl    = ts%data%tl%pccomp,                               &
          & nthreads     = ts%defn%nthreads )

      EndIf
        
      If( ts%data%rttov_errorstatus /= 0 ) Then
        err = errorstatus_fatal
        THROW(err.ne.0)
      EndIf
        
      If( ts%defn%do_print ) Then
        If( ts%defn%scale_out /= 1.0_jprb ) Then
          ts%data%tl%radiance%total = ts%data%tl%radiance%total * ts%defn%scale_out
          ts%data%tl%radiance%bt    = ts%data%tl%radiance%bt * ts%defn%scale_out
          ts%data%tl%radiance%refl  = ts%data%tl%radiance%refl * ts%defn%scale_out
          ts%data%tl%radiance%clear      = ts%data%tl%radiance%clear * ts%defn%scale_out
          ts%data%tl%radiance%bt_clear   = ts%data%tl%radiance%bt_clear * ts%defn%scale_out
          ts%data%tl%radiance%refl_clear = ts%data%tl%radiance%refl_clear * ts%defn%scale_out
          ts%data%tl%radiance%overcast = ts%data%tl%radiance%overcast * ts%defn%scale_out
          ts%data%tl%radiance%cloudy = ts%data%tl%radiance%cloudy * ts%defn%scale_out

          ts%data%tl%transmission%tau_total = ts%data%tl%transmission%tau_total * ts%defn%scale_out
          ts%data%tl%transmission%tau_levels = ts%data%tl%transmission%tau_levels * ts%defn%scale_out
          ts%data%tl%transmission%tausun_total_path1 = ts%data%tl%transmission%tausun_total_path1 * ts%defn%scale_out
          ts%data%tl%transmission%tausun_levels_path1 = ts%data%tl%transmission%tausun_levels_path1 * ts%defn%scale_out
          ts%data%tl%transmission%tausun_total_path2 = ts%data%tl%transmission%tausun_total_path2 * ts%defn%scale_out
          ts%data%tl%transmission%tausun_levels_path2 = ts%data%tl%transmission%tausun_levels_path2 * ts%defn%scale_out
          
          ts%data%tl%emissivity%emis_out = ts%data%tl%emissivity%emis_out * ts%defn%scale_out
          ts%data%tl%reflectance%refl_out = ts%data%tl%reflectance%refl_out * ts%defn%scale_out
          
          If( ts%defn%opts%rt_ir%pc%addpc ) Then
            If(Associated(ts%data%tl%pccomp%pcscores)) &
             & ts%data%tl%pccomp%pcscores = ts%data%tl%pccomp%pcscores * ts%defn%scale_out
            If(Associated(ts%data%tl%pccomp%total_pccomp)) &
             & ts%data%tl%pccomp%total_pccomp = ts%data%tl%pccomp%total_pccomp * ts%defn%scale_out
            If(Associated(ts%data%tl%pccomp%bt_pccomp)) &
             & ts%data%tl%pccomp%bt_pccomp = ts%data%tl%pccomp%bt_pccomp * ts%defn%scale_out
          End If
        End If

        Call PrintRadiance (err, ts%data%direct%radiance, "RADIANCE", Trim(ts%path)//"/tl/radiance.txt")
        THROW(err.ne.0)
        Call PrintTransmission (err, ts%data%direct%transmission, "TRANSMISSION", Trim(ts%path)//"/tl/transmission.txt" )
        THROW(err.ne.0)
        Call PrintEmissivity (err, ts%data%direct%emissivity%emis_out, "EMISSIVITY_OUT", Trim(ts%path)//"/tl/emissivity_out.txt" )
        THROW(err.ne.0)
        If (ts%defn%opts%rt_ir%addsolar) Then
          Call PrintReflectance (err, ts%data%direct%reflectance%refl_out, "REFLECTANCE_OUT", &
            & Trim(ts%path)//"/tl/reflectance_out.txt" )
          THROW(err.ne.0)
        EndIf
        Call PrintRadiance (err, ts%data%tl%radiance, "RADIANCE_TL", Trim(ts%path)//"/tl/radiance_tl.txt")
        THROW(err.ne.0)
        Call PrintTransmission (err, ts%data%tl%transmission, "TRANSMISSION_TL", Trim(ts%path)//"/tl/transmission_tl.txt" )
        THROW(err.ne.0)
        Call PrintEmissivity (err, ts%data%tl%emissivity%emis_out, "EMISSIVITY_OUT_TL", Trim(ts%path)//"/tl/emissivity_out_tl.txt" )
        THROW(err.ne.0)
        If (ts%defn%opts%rt_ir%addsolar) Then
          Call PrintReflectance (err, ts%data%tl%reflectance%refl_out, "REFLECTANCE_OUT_TL", &
            & Trim(ts%path)//"/tl/reflectance_out_tl.txt" )
          THROW(err.ne.0)
        EndIf
        If( ts%defn%opts%rt_ir%pc%addpc ) Then
          Call PrintPCScores (err, ts%data%direct%pccomp, "PCSCORES", Trim(ts%path)//"/tl/pcscores.txt" )
          THROW(err.ne.0)
          Call PrintPCScores (err, ts%data%tl%pccomp, "PCSCORES_TL", Trim(ts%path)//"/tl/pcscores_tl.txt" )
          THROW(err.ne.0)
        EndIf
      EndIf
        
      ts%defn%opts%config%verbose = .false.
    
    EndIf 
    
    If( ts%defn%do_ad ) Then
      
      ! Initialise the adjoint variables
      call rttov_init_prof(ts%data%ad%profiles)
      ts%data%ad%emissivity(:)%emis_in = 0._jprb
      ts%data%ad%emissivity(:)%emis_out = 0._jprb
      ts%data%ad%reflectance(:)%refl_in = 0._jprb
      ts%data%ad%reflectance(:)%refl_out = 0._jprb
      If( ts%defn%opts%rt_ir%pc%addpc )Then
        call rttov_copy_pccomp(ts%data%ad%pccomp, ts%data%ad%pccomp_saved)
        call rttov_init_rad(ts%data%ad%radiance)
      Else
        ts%data%ad%radiance%total(:) = ts%data%ad%radiance_saved%total(:)
        ts%data%ad%radiance%bt(:)    = ts%data%ad%radiance_saved%bt(:)
      EndIf
      
      If( (ts%defn%nthreads .le. 0) .and. ts%defn%prof_by_prof ) Then

        call rttov_parallel_ad(                                               &
          & ts%data%rttov_errorstatus,       ts%data%chanprof,                &
          & ts%defn%opts,                    ts%data%direct%profiles,         &
          & ts%data%ad%profiles,             ts%data%coefs,                   &
          & ts%data%direct%transmission,     ts%data%ad%transmission,         &
          & ts%data%direct%radiance,         ts%data%ad%radiance,             &
          & ts%data%calcemis,                ts%data%direct%emissivity,       &
          & ts%data%ad%emissivity,           ts%data%calcrefl,                &
          & ts%data%direct%reflectance,      ts%data%ad%reflectance,          &
          & ts%data%aer_opt_param,           ts%data%cld_opt_param,           &
          & nthreads = 1_jpim,               strategy = 1_jpim )

      Else If( (ts%defn%nthreads .le. 0) .and. ts%defn%chan_by_chan ) Then

        call rttov_parallel_ad(                                               &
          & ts%data%rttov_errorstatus,       ts%data%chanprof,                &
          & ts%defn%opts,                    ts%data%direct%profiles,         &
          & ts%data%ad%profiles,             ts%data%coefs,                   &
          & ts%data%direct%transmission,     ts%data%ad%transmission,         &
          & ts%data%direct%radiance,         ts%data%ad%radiance,             &
          & ts%data%calcemis,                ts%data%direct%emissivity,       &
          & ts%data%ad%emissivity,           ts%data%calcrefl,                &
          & ts%data%direct%reflectance,      ts%data%ad%reflectance,          &
          & ts%data%aer_opt_param,           ts%data%cld_opt_param,           &
          & nthreads = 1_jpim,               strategy = 2_jpim )
          
      Else If( ts%defn%nthreads .le. 0 ) Then
    
        call rttov_ad(                             &
          & ts%data%rttov_errorstatus,             &
          & ts%data%chanprof,                      &
          & ts%defn%opts,                          &
          & ts%data%direct%profiles,               &
          & ts%data%ad%profiles,                   &
          & ts%data%coefs,                         &
          & ts%data%direct%transmission,           &
          & ts%data%ad%transmission,               &
          & ts%data%direct%radiance,               &
          & ts%data%ad%radiance,                   &
          & ts%data%calcemis,                      &
          & ts%data%direct%emissivity,             &
          & ts%data%ad%emissivity,                 &
          & ts%data%calcrefl,                      &
          & ts%data%direct%reflectance,            &
          & ts%data%ad%reflectance,                &
          & ts%data%aer_opt_param,                 &
          & ts%data%cld_opt_param,                 &
          & channels_rec = ts%data%channels_rec,   &
          & pccomp       = ts%data%direct%pccomp,  &
          & pccomp_ad    = ts%data%ad%pccomp,      &
          & traj         = traj,                   &
          & traj_ad      = traj_ad )

      Else

        call rttov_parallel_ad(                                               &
          & ts%data%rttov_errorstatus,       ts%data%chanprof,                &
          & ts%defn%opts,                    ts%data%direct%profiles,         &
          & ts%data%ad%profiles,             ts%data%coefs,                   &
          & ts%data%direct%transmission,     ts%data%ad%transmission,         &
          & ts%data%direct%radiance,         ts%data%ad%radiance,             &
          & ts%data%calcemis,                ts%data%direct%emissivity,       &
          & ts%data%ad%emissivity,           ts%data%calcrefl,                &
          & ts%data%direct%reflectance,      ts%data%ad%reflectance,          &
          & ts%data%aer_opt_param,           ts%data%cld_opt_param,           &
          & channels_rec = ts%data%channels_rec,                              &
          & pccomp       = ts%data%direct%pccomp,                             &
          & pccomp_ad    = ts%data%ad%pccomp,                                 &
          & nthreads     = ts%defn%nthreads )
        
      EndIf
      
      If( ts%data%rttov_errorstatus /= 0 ) Then
        err = errorstatus_fatal
        THROW(err.ne.0)
      EndIf
        
      If( ts%defn%do_print ) Then

        If( ts%defn%scale_out /= 1.0_jprb )Then
          Call rttov_scale_profile_inc(ts%data%ad%profiles, ts%defn%scale_out)
        End If

        Call PrintRadiance (err, ts%data%direct%radiance, "RADIANCE", Trim(ts%path)//"/ad/radiance.txt")
        THROW(err.ne.0)
        Call PrintTransmission (err, ts%data%direct%transmission, "TRANSMISSION", Trim(ts%path)//"/ad/transmission.txt" )
        THROW(err.ne.0)
        Call PrintEmissivity (err, ts%data%direct%emissivity%emis_out, "EMISSIVITY_OUT", Trim(ts%path)//"/ad/emissivity_out.txt" )
        THROW(err.ne.0)
        If (ts%defn%opts%rt_ir%addsolar) Then
          Call PrintReflectance (err, ts%data%direct%reflectance%refl_out, "REFLECTANCE_OUT", &
            & Trim(ts%path)//"/ad/reflectance_out.txt" )
          THROW(err.ne.0)
        EndIf
        Call PrintProfiles ( err, ts%data%ad%profiles, "PROFILES_AD", Trim(ts%path)//"/ad/profiles_ad.txt" )
        THROW(err.ne.0)
        Call PrintEmissivity (err, ts%data%ad%emissivity%emis_in * ts%defn%scale_out, "EMISSIVITY_AD", &
          Trim(ts%path)//"/ad/emissivity_ad.txt")
        THROW(err.ne.0)
        If (ts%defn%opts%rt_ir%addsolar) Then
          Call PrintReflectance (err, ts%data%ad%reflectance%refl_in * ts%defn%scale_out, "REFLECTANCE_AD", &
            & Trim(ts%path)//"/ad/reflectance_ad.txt" )
          THROW(err.ne.0)
        EndIf
        If( ts%defn%opts%rt_ir%pc%addpc ) Then
          Call PrintPCScores (err, ts%data%direct%pccomp, "PCSCORES", Trim(ts%path)//"/ad/pcscores.txt" )
          THROW(err.ne.0)
          Call PrintPCScores (err, ts%data%ad%pccomp, "PCSCORES_AD", Trim(ts%path)//"/ad/pcscores_ad.txt" )
          THROW(err.ne.0)
        EndIf
      EndIf

      ts%defn%opts%config%verbose = .false.

    EndIf

    If( ts%defn%do_k ) Then

      call rttov_init_prof(ts%data%k%profiles)
      ts%data%k%emissivity(:)%emis_in = 0._jprb
      ts%data%k%emissivity(:)%emis_out = 0._jprb
      ts%data%k%reflectance(:)%refl_in = 0._jprb
      ts%data%k%reflectance(:)%refl_out = 0._jprb
      If( ts%defn%opts%rt_ir%pc%addpc )Then
        If( ts%defn%opts%rt_ir%pc%addradrec )Then
          If ( ts%defn%opts%rt_all%switchrad )Then
            ts%data%k%pccomp%bt_pccomp = 1._jprb
          Else
            ts%data%k%pccomp%total_pccomp = 1._jprb
          EndIf
        Else
          ts%data%k%pccomp%pcscores = 1._jprb
        EndIf
      Else
        ts%data%k%radiance%total = 1._jprb ! For any pure-solar channels, increment is always in radiance
        If( ts%defn%opts%rt_all%switchrad ) ts%data%k%radiance%bt = 1._jprb
      EndIf

      If( (ts%defn%nthreads .le. 0) .and. ts%defn%prof_by_prof ) Then

        Call rttov_parallel_k(                                           &
          & ts%data%rttov_errorstatus,     ts%data%chanprof,             &
          & ts%defn%opts,                  ts%data%direct%profiles,      & 
          & ts%data%k%profiles,            ts%data%coefs,                & 
          & ts%data%direct%transmission,   ts%data%k%transmission,       &
          & ts%data%direct%radiance,       ts%data%k%radiance,           &
          & ts%data%calcemis,              ts%data%direct%emissivity,    &
          & ts%data%k%emissivity,          ts%data%calcrefl,             &
          & ts%data%direct%reflectance,    ts%data%k%reflectance,        &
          & ts%data%aer_opt_param,         ts%data%cld_opt_param,        &
          & nthreads = 1_jpim,             strategy = 1_jpim )


      Else If( (ts%defn%nthreads .le. 0) .and. ts%defn%chan_by_chan ) Then

        Call rttov_parallel_k(                                           &
          & ts%data%rttov_errorstatus,     ts%data%chanprof,             &
          & ts%defn%opts,                  ts%data%direct%profiles,      &
          & ts%data%k%profiles,            ts%data%coefs,                &
          & ts%data%direct%transmission,   ts%data%k%transmission,       &
          & ts%data%direct%radiance,       ts%data%k%radiance,           &
          & ts%data%calcemis,              ts%data%direct%emissivity,    &
          & ts%data%k%emissivity,          ts%data%calcrefl,             &
          & ts%data%direct%reflectance,    ts%data%k%reflectance,        &
          & ts%data%aer_opt_param,         ts%data%cld_opt_param,        &
          & nthreads = 1_jpim,             strategy = 2_jpim )
          
      Else If( ts%defn%nthreads .le. 0 ) Then
    
        Call rttov_k(                                       &
          & ts%data%rttov_errorstatus,                      &
          & ts%data%chanprof,                               &
          & ts%defn%opts,                                   &
          & ts%data%direct%profiles,                        &
          & ts%data%k%profiles,                             &
          & ts%data%coefs,                                  &
          & ts%data%direct%transmission,                    &
          & ts%data%k%transmission,                         &
          & ts%data%direct%radiance,                        &
          & ts%data%k%radiance,                             &
          & ts%data%calcemis,                               &
          & ts%data%direct%emissivity,                      &
          & ts%data%k%emissivity,                           &
          & ts%data%calcrefl,                               &
          & ts%data%direct%reflectance,                     &
          & ts%data%k%reflectance,                          &
          & ts%data%aer_opt_param,                          &
          & ts%data%cld_opt_param,                          &
          & traj           = traj,                          &
          & traj_k         = traj_k,                        &
          & pccomp         = ts%data%direct%pccomp,         &
          & pccomp_k       = ts%data%k%pccomp,              &
          & profiles_k_rec = ts%data%data_k%profiles_k_rec, &
          & profiles_k_pc  = ts%data%data_k%profiles_k_pc,  &
          & channels_rec   = ts%data%channels_rec )

      Else

        Call rttov_parallel_k(                                           &
          & ts%data%rttov_errorstatus,     ts%data%chanprof,             &
          & ts%defn%opts,                  ts%data%direct%profiles,      &
          & ts%data%k%profiles,            ts%data%coefs,                &
          & ts%data%direct%transmission,   ts%data%k%transmission,       &
          & ts%data%direct%radiance,       ts%data%k%radiance,           &
          & ts%data%calcemis,              ts%data%direct%emissivity,    &
          & ts%data%k%emissivity,          ts%data%calcrefl,             &
          & ts%data%direct%reflectance,    ts%data%k%reflectance,        &
          & ts%data%aer_opt_param,         ts%data%cld_opt_param,        &
          & pccomp         = ts%data%direct%pccomp,                      &
          & pccomp_k       = ts%data%k%pccomp,                           &
          & profiles_k_rec = ts%data%data_k%profiles_k_rec,              &
          & profiles_k_pc  = ts%data%data_k%profiles_k_pc,               &
          & channels_rec   = ts%data%channels_rec,                       &
          & nthreads       = ts%defn%nthreads)
      
      EndIf
    
      If( ts%data%rttov_errorstatus /= 0 ) Then
        err = errorstatus_fatal
        THROW(err.ne.0)
      EndIf
        
      If( ts%defn%do_print ) Then
        Call PrintRadiance (err, ts%data%direct%radiance, "RADIANCE",  Trim(ts%path)//"/k/radiance.txt")
        THROW(err.ne.0)
        Call PrintTransmission (err, ts%data%direct%transmission, "TRANSMISSION",  Trim(ts%path)//"/k/transmission.txt" )
        THROW(err.ne.0)
        Call PrintEmissivity (err, ts%data%direct%emissivity%emis_out, "EMISSIVITY_OUT", Trim(ts%path)//"/k/emissivity_out.txt" )
        THROW(err.ne.0)
        If (ts%defn%opts%rt_ir%addsolar) Then
          Call PrintReflectance (err, ts%data%direct%reflectance%refl_out, "REFLECTANCE_OUT", &
            & Trim(ts%path)//"/k/reflectance_out.txt" )
          THROW(err.ne.0)
        EndIf
        Call PrintProfiles ( err, ts%data%k%profiles, "PROFILES_K",  Trim(ts%path)//"/k/profiles_k.txt" )
        THROW(err.ne.0)
        Call PrintEmissivity (err, ts%data%k%emissivity%emis_in, "EMISSIVITY_K",  Trim(ts%path)//"/k/emissivity_k.txt" )
        THROW(err.ne.0)
        If (ts%defn%opts%rt_ir%addsolar) Then
          Call PrintReflectance (err, ts%data%k%reflectance%refl_in, "REFLECTANCE_K", &
            & Trim(ts%path)//"/k/reflectance_k.txt" )
          THROW(err.ne.0)
        EndIf
        If( ts%defn%opts%rt_ir%pc%addpc ) Then
          Call PrintPCScores (err, ts%data%direct%pccomp, "PCSCORES", Trim(ts%path)//"/k/pcscores.txt" )
          THROW(err.ne.0)
          Call PrintPCScores (err, ts%data%k%pccomp, "PCSCORES_K", Trim(ts%path)//"/k/pcscores_k.txt" )
          THROW(err.ne.0)
          If( ts%defn%opts%rt_ir%pc%addradrec ) then
            Call PrintProfiles ( err, ts%data%data_k%profiles_k_rec, "PROFILES_K_REC",  Trim(ts%path)//"/k/profiles_k_rec.txt" )
            THROW(err.ne.0)
          else
            Call PrintProfiles ( err, ts%data%data_k%profiles_k_pc, "PROFILES_K_PC",  Trim(ts%path)//"/k/profiles_k_pc.txt" )
            THROW(err.ne.0)
          endif
        EndIf
      EndIf


#ifdef _RTTOV_HDF
      If( ts%defn%savehdf5 ) Then
       INFO("save kmatrix")
       CALL RTTOV_HDF_SAVE( ERR, TRIM(ts%path)//"/k/kmatrix.H5",     '/KMATRIX', CREATE=.true., &
       &  KMATRIX = ts%data%k%profiles, OPTS= ts%defn%opts)
       THROW(err.ne.0)
                 
       CALL RTTOV_HDF_SAVE( ERR, TRIM(ts%path)//"/k/kmatrix.H5",     '/MISC', CREATE=.false., &
       &  R1 = ts%data%coefs%coef%ff_cwn(ts%data%chanprof(:)%chan), SNAME='WAVENUMBERS', UNITS='cm-1' )
       THROW(err.ne.0)
       
       CALL RTTOV_HDF_SAVE( ERR, TRIM(ts%path)//"/k/kmatrix.H5",     '/MISC', CREATE=.false., &
       &  c0 = ts%data%coefs%coef%ID_COMMON_NAME, SNAME='ID_COMMON_NAME')
       THROW(err.ne.0)
       
       CALL RTTOV_HDF_SAVE( ERR, TRIM(ts%path)//"/k/kmatrix.H5",     '/MISC', CREATE=.false., &
       &  c0 = ts%defn%f_coef, SNAME='COEF_FILENAME')
       THROW(err.ne.0)
       
       CALL RTTOV_HDF_SAVE( ERR, TRIM(ts%path)//"/k/kmatrix.H5",     '/PROFILES', CREATE=.false., &
       &  PROFILES = ts%data%direct%profiles)
       THROW(err.ne.0)
                 
       CALL RTTOV_HDF_SAVE( ERR, TRIM(ts%path)//"/k/kmatrix.H5",     '/CHANPROF', CREATE=.false., &
       &  CHANPROF = ts%data%chanprof)
       THROW(err.ne.0)
                 
       CALL RTTOV_HDF_SAVE( ERR, TRIM(ts%path)//"/k/kmatrix.H5",     '/EMISSIVITY_K', CREATE=.false., &
       &  EMISSIVITY = ts%data%k%emissivity)
       THROW(err.ne.0)
       CALL RTTOV_HDF_SAVE( ERR, TRIM(ts%path)//"/k/kmatrix.H5",     '/EMISSIVITY_K', CREATE=.false., &
       &  l1 = ts%data%calcemis, SNAME='CALCEMIS')
       THROW(err.ne.0)
         
       CALL RTTOV_HDF_SAVE( ERR, TRIM(ts%path)//"/k/kmatrix.H5",     '/OPTIONS', CREATE=.false., &
       &  OPTIONS = ts%defn%opts)
       THROW(err.ne.0)
       
       If (ts%defn%opts%rt_ir%addsolar) Then
         CALL RTTOV_HDF_SAVE( ERR, TRIM(ts%path)//"/k/kmatrix.H5",     '/REFLECTANCE_K', CREATE=.false., &
       &    REFLECTANCE = ts%data%k%reflectance)
         THROW(err.ne.0)
         CALL RTTOV_HDF_SAVE( ERR, TRIM(ts%path)//"/k/kmatrix.H5",     '/REFLECTANCE_K', CREATE=.false., &
       &  l1 = ts%data%calcrefl, SNAME='CALCREFL')
       THROW(err.ne.0)
       EndIf
       INFO("fin save")
      EndIf
#endif

      ts%defn%opts%config%verbose = .false.
  
    EndIf

  EndDo
    
  If( ts%defn%do_k_bf ) Then
  
    Call rttov_k_bf(                                       &
      & ts%data%rttov_errorstatus,                         &
      & ts%data%chanprof,                                  &
      & ts%defn%opts,                                      &
      & ts%data%direct%profiles,                           &
      & ts%data%k_bf%profiles,                             &
      & ts%data%coefs,                                     &
      & ts%data%direct%transmission,                       &
      & ts%data%direct%radiance,                           &
      & ts%data%calcemis,                                  &
      & ts%data%direct%emissivity,                         &
      & ts%data%k_bf%emissivity,                           &
      & ts%data%calcrefl,                                  &
      & ts%data%direct%reflectance,                        &
      & ts%data%k_bf%reflectance,                          &
      & ts%data%aer_opt_param,                             &
      & ts%data%cld_opt_param,                             &
      & pccomp         = ts%data%direct%pccomp,            &
      & profiles_k_rec = ts%data%data_k_bf%profiles_k_rec, &
      & profiles_k_pc  = ts%data%data_k_bf%profiles_k_pc,  &
      & channels_rec   = ts%data%channels_rec,             &
      & nthreads       = ts%defn%nthreads )

    If( ts%data%rttov_errorstatus /= 0 ) Then
      err = errorstatus_fatal
      THROW(err.ne.0)
    EndIf
      
    If( ts%defn%do_print ) Then
      Call PrintRadiance (err, ts%data%direct%radiance, "RADIANCE",  Trim(ts%path)//"/k_bf/radiance.txt")
      THROW(err.ne.0)
      Call PrintTransmission (err, ts%data%direct%transmission, "TRANSMISSION",  Trim(ts%path)//"/k_bf/transmission.txt" )
      THROW(err.ne.0)
      Call PrintEmissivity (err, ts%data%direct%emissivity%emis_out, "EMISSIVITY_OUT", Trim(ts%path)//"/k_bf/emissivity_out.txt" )
      THROW(err.ne.0)
      If (ts%defn%opts%rt_ir%addsolar) Then
        Call PrintReflectance (err, ts%data%direct%reflectance%refl_out, "REFLECTANCE_OUT", &
          & Trim(ts%path)//"/k_bf/reflectance_out.txt" )
        THROW(err.ne.0)
      EndIf
      Call PrintProfiles ( err, ts%data%k_bf%profiles, "PROFILES_K",  Trim(ts%path)//"/k_bf/profiles_k.txt" )
      THROW(err.ne.0)
      Call PrintEmissivity (err, ts%data%k_bf%emissivity%emis_in, "EMISSIVITY_K",  Trim(ts%path)//"/k_bf/emissivity_k.txt" )
      THROW(err.ne.0)
      If (ts%defn%opts%rt_ir%addsolar) Then
        Call PrintReflectance (err, ts%data%k_bf%reflectance%refl_in, "REFLECTANCE_K", &
          & Trim(ts%path)//"/k_bf/reflectance_k.txt" )
        THROW(err.ne.0)
      EndIf
      If( ts%defn%opts%rt_ir%pc%addpc ) Then
        Call PrintPCScores (err, ts%data%direct%pccomp, "PCSCORES", Trim(ts%path)//"/k_bf/pcscores.txt" )
        THROW(err.ne.0)
        Call PrintPCScores (err, ts%data%k_bf%pccomp, "PCSCORES_K", Trim(ts%path)//"/k_bf/pcscores_k.txt" )
        THROW(err.ne.0)
        If( ts%defn%opts%rt_ir%pc%addradrec ) then
          Call PrintProfiles ( err, ts%data%data_k_bf%profiles_k_rec, "PROFILES_K_REC",  Trim(ts%path)//"/k_bf/profiles_k_rec.txt" )
          THROW(err.ne.0)
        else
          Call PrintProfiles ( err, ts%data%data_k_bf%profiles_k_pc, "PROFILES_K_PC",  Trim(ts%path)//"/k_bf/profiles_k_pc.txt" )
          THROW(err.ne.0)
        endif
      EndIf
    EndIf
  
    ts%defn%opts%config%verbose = .false.
  
  EndIf
  
  
  If( ts%defn%do_k_tl ) Then
  
    Call rttov_k_tl(                                       &
      & ts%data%rttov_errorstatus,                         &
      & ts%data%chanprof,                                  &
      & ts%defn%opts,                                      &
      & ts%data%direct%profiles,                           &
      & ts%data%k_tl%profiles,                             &
      & ts%data%coefs,                                     &
      & ts%data%direct%transmission,                       &
      & ts%data%direct%radiance,                           &
      & ts%data%calcemis,                                  &
      & ts%data%direct%emissivity,                         &
      & ts%data%k_tl%emissivity,                           &
      & ts%data%calcrefl,                                  &
      & ts%data%direct%reflectance,                        &
      & ts%data%k_tl%reflectance,                          &
      & ts%data%aer_opt_param,                             &
      & ts%data%cld_opt_param,                             &
      & pccomp         = ts%data%direct%pccomp,            &
      & profiles_k_rec = ts%data%data_k_tl%profiles_k_rec, &
      & profiles_k_pc  = ts%data%data_k_tl%profiles_k_pc,  &
      & channels_rec   = ts%data%channels_rec,             &
      & nthreads       = ts%defn%nthreads )
  
    If( ts%data%rttov_errorstatus /= 0 ) Then
      err = errorstatus_fatal
      THROW(err.ne.0)
    EndIf
      
    If( ts%defn%do_print ) Then
      Call PrintRadiance (err, ts%data%direct%radiance, "RADIANCE",  Trim(ts%path)//"/k_tl/radiance.txt")
      THROW(err.ne.0)
      Call PrintTransmission (err, ts%data%direct%transmission, "TRANSMISSION",  Trim(ts%path)//"/k_tl/transmission.txt" )
      THROW(err.ne.0)
      Call PrintEmissivity (err, ts%data%direct%emissivity%emis_out, "EMISSIVITY_OUT", Trim(ts%path)//"/k_tl/emissivity_out.txt" )
      THROW(err.ne.0)
      If (ts%defn%opts%rt_ir%addsolar) Then
        Call PrintReflectance (err, ts%data%direct%reflectance%refl_out, "REFLECTANCE_OUT", &
          & Trim(ts%path)//"/k_tl/reflectance_out.txt" )
        THROW(err.ne.0)
      EndIf
      Call PrintProfiles ( err, ts%data%k_tl%profiles, "PROFILES_K",  Trim(ts%path)//"/k_tl/profiles_k.txt" )
      THROW(err.ne.0)
      Call PrintEmissivity (err, ts%data%k_tl%emissivity%emis_in, "EMISSIVITY_K",  Trim(ts%path)//"/k_tl/emissivity_k.txt" )
      THROW(err.ne.0)
      If (ts%defn%opts%rt_ir%addsolar) Then
        Call PrintReflectance (err, ts%data%k_tl%reflectance%refl_in, "REFLECTANCE_K", &
          & Trim(ts%path)//"/k_tl/reflectance_k.txt" )
        THROW(err.ne.0)
      EndIf
      If( ts%defn%opts%rt_ir%pc%addpc ) Then
        Call PrintPCScores (err, ts%data%direct%pccomp, "PCSCORES", Trim(ts%path)//"/k_tl/pcscores.txt" )
        THROW(err.ne.0)
        Call PrintPCScores (err, ts%data%k_tl%pccomp, "PCSCORES_K", Trim(ts%path)//"/k_tl/pcscores_k.txt" )
        THROW(err.ne.0)
        If( ts%defn%opts%rt_ir%pc%addradrec ) then
          Call PrintProfiles ( err, ts%data%data_k_tl%profiles_k_rec, "PROFILES_K_REC",  Trim(ts%path)//"/k_tl/profiles_k_rec.txt" )
          THROW(err.ne.0)
        else
          Call PrintProfiles ( err, ts%data%data_k_tl%profiles_k_pc, "PROFILES_K_PC",  Trim(ts%path)//"/k_tl/profiles_k_pc.txt" )
          THROW(err.ne.0)
        endif
      EndIf
    EndIf
  
    ts%defn%opts%config%verbose = .false.
  
  EndIf
  
  
  If( ts%defn%do_k_ad ) Then
  
    Call rttov_k_ad(                                       &
      & ts%data%rttov_errorstatus,                         &
      & ts%data%chanprof,                                  &
      & ts%defn%opts,                                      &
      & ts%data%direct%profiles,                           &
      & ts%data%k_ad%profiles,                             &
      & ts%data%coefs,                                     &
      & ts%data%direct%transmission,                       &
      & ts%data%direct%radiance,                           &
      & ts%data%calcemis,                                  &
      & ts%data%direct%emissivity,                         &
      & ts%data%k_ad%emissivity,                           &
      & ts%data%calcrefl,                                  &
      & ts%data%direct%reflectance,                        &
      & ts%data%k_ad%reflectance,                          &
      & ts%data%aer_opt_param,                             &
      & ts%data%cld_opt_param,                             &
      & pccomp         = ts%data%direct%pccomp,            &
      & profiles_k_rec = ts%data%data_k_ad%profiles_k_rec, &
      & profiles_k_pc  = ts%data%data_k_ad%profiles_k_pc,  &
      & channels_rec   = ts%data%channels_rec,             &
      & nthreads       = ts%defn%nthreads )
  
    If( ts%data%rttov_errorstatus /= 0 ) Then
      err = errorstatus_fatal
      THROW(err.ne.0)
    EndIf
      
    If( ts%defn%do_print ) Then
      Call PrintRadiance (err, ts%data%direct%radiance, "RADIANCE",  Trim(ts%path)//"/k_ad/radiance.txt")
      THROW(err.ne.0)
      Call PrintTransmission (err, ts%data%direct%transmission, "TRANSMISSION", Trim(ts%path)//"/k_ad/transmission.txt" )
      THROW(err.ne.0)
      Call PrintEmissivity (err, ts%data%direct%emissivity%emis_out, "EMISSIVITY_OUT", Trim(ts%path)//"/k_ad/emissivity_out.txt" )
      THROW(err.ne.0)
      If (ts%defn%opts%rt_ir%addsolar) Then
        Call PrintReflectance (err, ts%data%direct%reflectance%refl_out, "REFLECTANCE_OUT", &
          & Trim(ts%path)//"/k_ad/reflectance_out.txt" )
        THROW(err.ne.0)
      EndIf
      Call PrintProfiles ( err, ts%data%k_ad%profiles, "PROFILES_K",  Trim(ts%path)//"/k_ad/profiles_k.txt" )
      THROW(err.ne.0)
      Call PrintEmissivity (err, ts%data%k_ad%emissivity%emis_in, "EMISSIVITY_K",  Trim(ts%path)//"/k_ad/emissivity_k.txt" )
      THROW(err.ne.0)
      If (ts%defn%opts%rt_ir%addsolar) Then
        Call PrintReflectance (err, ts%data%k_ad%reflectance%refl_in, "REFLECTANCE_K", &
          & Trim(ts%path)//"/k_ad/reflectance_k.txt" )
        THROW(err.ne.0)
      EndIf
      If( ts%defn%opts%rt_ir%pc%addpc ) Then
        Call PrintPCScores (err, ts%data%direct%pccomp, "PCSCORES", Trim(ts%path)//"/k_ad/pcscores.txt" )
        THROW(err.ne.0)
        Call PrintPCScores (err, ts%data%k_ad%pccomp, "PCSCORES_K", Trim(ts%path)//"/k_ad/pcscores_k.txt" )
        THROW(err.ne.0)
        If( ts%defn%opts%rt_ir%pc%addradrec ) then
          Call PrintProfiles ( err, ts%data%data_k_ad%profiles_k_rec, "PROFILES_K_REC",  Trim(ts%path)//"/k_ad/profiles_k_rec.txt" )
          THROW(err.ne.0)
        else
          Call PrintProfiles ( err, ts%data%data_k_ad%profiles_k_pc, "PROFILES_K_PC",  Trim(ts%path)//"/k_ad/profiles_k_pc.txt" )
          THROW(err.ne.0)
        endif
      EndIf
    EndIf
  
    ts%defn%opts%config%verbose = .false.
  
  EndIf
  
  If( ts%defn%do_taylor ) Then
    
    ! For IR sensors when addpc is false or the surface type is not sea,
    ! the calculated emissivity is independent of the profile and so
    ! will not reflect the perturbed profile. So we need to manually
    ! perturb the emissivity in this case.

    ! For IR sensors, rttov_calcsurfrefl_tl only calculates the TL for sea surfaces 
    ! when calcrefl is true, so need to manually perturb reflectance otherwise
    ! (over sea, the refl TL comes from the perturbed profile wind variables).
    ! However, over land and sea-ice, for thermal+solar channels, if there is a valid
    ! emissivity, this will be used to calculate a reflectance (and reflectance_tl)
    ! when calcrefl is true.
    
    If (ts%data%coefs%coef%id_sensor == sensor_id_ir .or. &
        ts%data%coefs%coef%id_sensor == sensor_id_hi) Then
      Do i = 1, ts%defn%nchannels
        j = ts%data%chanprof(i)%prof
        chan = ts%data%chanprof(i)%chan
        If (.not. (ts%defn%opts%rt_ir%pc%addpc .and. &
            ts%data%direct%profiles(j)%skin%surftype == surftype_sea)) Then
          
          ! Save the calculated emissivity for this channel and turn calcemis
          ! off so we can supply the perturbed emissivity to rttov_direct
          ts%data%direct%emissivity(i)%emis_in = ts%data%direct%emissivity(i)%emis_out
          ts%data%calcemis(i) = .false.
          
        EndIf
      
        If (Associated(ts%data%coefs%coef%ss_val_chn)) Then
          If (.not. (ts%data%calcrefl(i) .and. &
                      (ts%data%direct%profiles(j)%skin%surftype == surftype_sea .or. &
                      (ts%data%direct%profiles(j)%skin%surftype /= surftype_sea .and. &
                      ts%data%coefs%coef%ss_val_chn(chan) == 1_jpim .and. &
                      ts%data%direct%emissivity(i)%emis_in > 0._jprb)))) Then
            
            ! Save the calculated reflectance for this channel and turn calcrefl
            ! off so we can supply the perturbed reflectance to rttov_direct
            ts%data%direct%reflectance(i)%refl_in = ts%data%direct%reflectance(i)%refl_out
            ts%data%calcrefl(i) = .false.
            
          EndIf
        EndIf
        
      EndDo
    EndIf

    ! If TAYLOR_BY_CHAN is set Taylor ratios are calculated for every channel individually.
    ! Otherwise the ratios are calculated based on radiances/BTs/refls summed for each profile.
    ! For PC-RTTOV the Taylor test is always calculated per profile.
    If( ts%defn%taylor_by_chan .and. .not. ts%defn%opts%rt_ir%pc%addpc ) Then
      Allocate(taylor_direct(niter,ts%defn%nchannels), taylor_tl(niter,ts%defn%nchannels))
    Else
      Allocate(taylor_direct(niter,ts%defn%nprofiles), taylor_tl(niter,ts%defn%nprofiles))
    EndIf
    Do iter = 1, niter

      ! Scale the profile increments in place by a factor of 0.1 at each iteration
      If (iter > 1) Then 
        Call rttov_scale_profile_inc( ts%data%tl%profiles, 0.1_jprb )
        ts%data%tl%emissivity(:)%emis_in = ts%data%tl%emissivity(:)%emis_in * 0.1_jprb
        ts%data%tl%reflectance(:)%refl_in = ts%data%tl%reflectance(:)%refl_in * 0.1_jprb
      EndIf

      Call rttov_copy_prof(ts%data%direct%profiles_pert, ts%data%direct%profiles)
      Call rttov_add_prof(ts%data%direct%profiles_pert, ts%data%direct%profiles, ts%data%tl%profiles)

      ts%data%direct%emissivity_pert(:)%emis_in  = ts%data%direct%emissivity(:)%emis_in +  &
                                                 & ts%data%tl%emissivity(:)%emis_in
      ts%data%direct%reflectance_pert(:)%refl_in = ts%data%direct%reflectance(:)%refl_in + &
                                                 & ts%data%tl%reflectance(:)%refl_in
      ts%data%direct%reflectance_pert(:)%refl_cloud_top = ts%data%direct%reflectance(:)%refl_cloud_top

      Call rttov_direct(                                 &
        & ts%data%rttov_errorstatus,                     &
        & ts%data%chanprof,                              &
        & ts%defn%opts,                                  &
        & ts%data%direct%profiles_pert,                  &
        & ts%data%coefs,                                 &
        & ts%data%direct%transmission,                   &
        & ts%data%direct%radiance,                       &
        & calcemis = ts%data%calcemis,                   &
        & emissivity = ts%data%direct%emissivity_pert,   &
        & calcrefl = ts%data%calcrefl,                   &
        & reflectance = ts%data%direct%reflectance_pert, &
        & aer_opt_param = ts%data%aer_opt_param,         &
        & cld_opt_param = ts%data%cld_opt_param,         &
        & pccomp       = ts%data%direct%pccomp,          &
        & channels_rec = ts%data%channels_rec,           &
        & traj         = traj )

      If( ts%data%rttov_errorstatus /= 0 ) Then
        err = errorstatus_fatal
        THROW(err.ne.0)
      EndIf

      If( ts%defn%opts%rt_ir%pc%addpc ) Then
        Do i = 1, ts%defn%nprofiles
          If( ts%defn%opts%rt_ir%pc%addradrec ) Then
            nchanrec = ts%defn%nchannels_rec
            j = (i-1)*nchanrec+1
            taylor_direct(iter,i) = sum(ts%data%direct%pccomp%bt_pccomp(j:j+nchanrec-1) - &
                                        ts%data%direct%pccomp_saved%bt_pccomp(j:j+nchanrec-1))
            taylor_tl(iter,i) = sum(ts%data%tl%pccomp%bt_pccomp(j:j+nchanrec-1))
          Else
            npcscores = ts%defn%npcscores
            j = (i-1)*npcscores+1
            taylor_direct(iter,i) = sum(ts%data%direct%pccomp%pcscores(j:j+npcscores-1) - &
                                        ts%data%direct%pccomp_saved%pcscores(j:j+npcscores-1))
            taylor_tl(iter,i) = sum(ts%data%tl%pccomp%pcscores(j:j+npcscores-1))
          EndIf
        EndDo
      Else
        If( ts%defn%taylor_by_chan ) Then
          Do i = 1, ts%defn%nchannels
!             If (ts%data%direct%radiance%bt(i) > 0._jprb) Then
!               taylor_direct(iter,i) = ts%data%direct%radiance%bt(i) - ts%data%direct%radiance_saved%bt(i)
!               taylor_tl(iter,i) = ts%data%tl%radiance%bt(i)
!             ElseIf (ts%data%direct%radiance%refl(i) > 0._jprb) Then
!               taylor_direct(iter,i) = ts%data%direct%radiance%refl(i) - ts%data%direct%radiance_saved%refl(i)
!               taylor_tl(iter,i) = ts%data%tl%radiance%refl(i)
!             Else
              taylor_direct(iter,i) = ts%data%direct%radiance%total(i) - ts%data%direct%radiance_saved%total(i)
              taylor_tl(iter,i) = ts%data%tl%radiance%total(i)
!             EndIf
          EndDo
        Else
          j = 1
          Do i = 1, ts%defn%nprofiles
            taylor_direct(iter,i) = 0._jprb
            taylor_tl(iter,i) = 0._jprb
            Do While (ts%data%chanprof(j)%prof == i)
!               If (ts%data%direct%radiance%bt(i) > 0._jprb) Then
!                 taylor_direct(iter,i) = taylor_direct(iter,i) + ts%data%direct%radiance%bt(j) - &
!                                                                 ts%data%direct%radiance_saved%bt(j)
!                 taylor_tl(iter,i) = taylor_tl(iter,i) + ts%data%tl%radiance%bt(j)
!               Else If (ts%data%direct%radiance%refl(j) > 0._jprb) Then
!                 taylor_direct(iter,i) = taylor_direct(iter,i) + ts%data%direct%radiance%refl(j) - &
!                                                                 ts%data%direct%radiance_saved%refl(j)
!                 taylor_tl(iter,i) = taylor_tl(iter,i) + ts%data%tl%radiance%refl(j)
!               Else
                taylor_direct(iter,i) = taylor_direct(iter,i) + ts%data%direct%radiance%total(j) - &
                                                                ts%data%direct%radiance_saved%total(j)
                taylor_tl(iter,i) = taylor_tl(iter,i) + ts%data%tl%radiance%total(j)
!               EndIf
              j = j + 1
              If( j > size(ts%data%chanprof) ) Exit
            EndDo
            If( j > size(ts%data%chanprof) ) Exit
          EndDo
        EndIf
      EndIf
    EndDo

    Call rttov_get_lun (file_id)
    taylor_filename = Trim(ts%path)//'/taylor_test.log'
    Open( file_id, file = taylor_filename, form = 'formatted', iostat = err )
    THROWM(err.ne.0,"Cannot open "//Trim(taylor_filename))
    
    Write( file_id, '(a,i6)' ) 'nprof: ', size(taylor_direct(1,:))
    Write( file_id, '(a,i6)' ) 'niter: ', niter
    Write( file_id, '(a)' ) 'Taylor test: profile, channel, iteration, (H(x*)-H(x))/H(x*-x)'
    Do i = 1, size(taylor_direct(1,:))
      Do iter = 1, niter
        If (taylor_tl(iter,i) /= 0._jprb) Then
          ratio = taylor_direct(iter,i) * (10._jprb**(iter-1)) / taylor_tl(iter,i)
        Else
          ratio = 0._jprb
        EndIf
        If( ts%defn%taylor_by_chan .and. .not. ts%defn%opts%rt_ir%pc%addpc ) Then
          Write( file_id, '(14x,i5,4x,i5,4x,i6,5x,7f20.16)' ) &
              ts%data%chanprof(i)%prof, ts%data%chanprof(i)%chan, iter, ratio
        Else
          Write( file_id, '(14x,i5,4x,a9,i6,5x,7f20.16)' ) &
              i, '   all   ', iter, ratio
        EndIf
      EndDo
      Write( file_id, '(a)' ) ''
    EndDo
    
    Close( file_id, iostat = err )
    THROWM(err.ne.0,"Cannot close "//Trim(taylor_filename))
    Call rttov_put_lun (file_id)

    Deallocate(taylor_direct, taylor_tl)
  EndIf
  
  ts%defn%do_print = .false.
  

CATCH
  End Subroutine

  
  Subroutine alloc_traj( err, ts, asw )
!
    Implicit None
    
    Integer(Kind=jpim),      Intent(out)           :: err
    Type(rttov_test_struct), Intent(inout), Target :: ts
    Integer(Kind=jpim),      Intent(in)            :: asw
!
TRY

!     If( ts%defn%do_direct ) Then                   
      call rttov_alloc_traj(                               &
          & err,                  ts%defn%nprofiles,       &
          & ts%defn%nchannels,    ts%defn%opts,            &
          & ts%defn%nlevels,      ts%data%coefs,           &
          & asw,                  traj = ts%data%traj )
      THROW(err.ne.0)
!     EndIf
    
    If( ts%defn%do_tl ) Then                       
      call rttov_alloc_traj(                               &
          & err,                  ts%defn%nprofiles,       &
          & ts%defn%nchannels,    ts%defn%opts,            &
          & ts%defn%nlevels,      ts%data%coefs,           &
          & asw,                  traj_tl = ts%data%traj_tl )
      THROW(err.ne.0)
    EndIf
    
    If( ts%defn%do_ad ) Then                       
      call rttov_alloc_traj(                               &
          & err,                  ts%defn%nprofiles,       &
          & ts%defn%nchannels,    ts%defn%opts,            &
          & ts%defn%nlevels,      ts%data%coefs,           &
          & asw,                  traj_ad = ts%data%traj_ad )
      THROW(err.ne.0)
    EndIf
    
    If( ts%defn%do_k ) Then                        
      call rttov_alloc_traj(                               &
          & err,                  ts%defn%nprofiles,       &
          & ts%defn%nchannels,    ts%defn%opts,            &
          & ts%defn%nlevels,      ts%data%coefs,           &
          & asw,                  traj_k = ts%data%traj_k )
      THROW(err.ne.0)
    EndIf
    

CATCH

 
  End Subroutine

  Subroutine Setup( ts, path, err )
 !
    Implicit None
    
    Type(rttov_test_struct), Intent(inout), Target :: ts
    Character(len=*),        Intent(in)            :: path
    Integer(Kind=jpim),      Intent(out)           :: err
 !
    Integer(Kind=jpim), Pointer :: kchannels(:), & ! absolute channel number in coef file, full list of channels
                                   ichannels(:), & ! restricted list; each channel number appears once
                                   lchannels(:)    ! absolute channel number -> relative channel number
    Integer(Kind=jpim) :: i, j, lchannelsmax, ichannelsmax, chan, ss_val
    Integer(Kind=jpim) :: nprofiles1, nchannels1
    Integer(Kind=jpim) :: imult, k1, k2
    Integer(Kind=jpim) :: file_id
    Integer(Kind=jpim) :: file_id_nml
    Logical(Kind=jplm) :: file_exists
    Type(rttov_test_defn) :: defn
    Character(Len=256) :: interp_filename, interp_msg
    Namelist / coef_nml / defn
    Namelist / rttov_test_nml / defn
TRY

    ! pgf90 requires these optional variables to be initialised for multiple-instrument calls
    ! otherwise it will "remember" the values from the previous instrument if they are not
    ! overwritten with new values when the namelists are read in.
    defn%f_o3            = ""
    defn%f_co2           = ""
    defn%f_n2o           = ""
    defn%f_co            = ""
    defn%f_ch4           = ""
    defn%f_clw           = ""
    defn%f_aer_opt_param = ""
    defn%f_cld_opt_param = ""
    defn%f_aerosl        = ""
    defn%f_cloud         = ""
    defn%f_cfrac         = ""
    defn%f_pcscores      = ""
    defn%f_channels_rec  = ""
    defn%f_scat          = ""
    defn%f_optp          = ""
    defn%f_pccomp        = ""

    call rttov_get_lun(file_id_nml )

    Open( file_id_nml, file = Trim(path)//'/rttov_test.txt', iostat = err)
    THROWM(err.ne.0,"Cannot open "//Trim(path)//'/rttov_test.txt')
    Read( file_id_nml, nml = rttov_test_nml, iostat = err)
    THROWM(err.ne.0,"Cannot read "//Trim(path)//'/rttov_test.txt')
    Close( file_id_nml, iostat = err )
    THROWM(err.ne.0,"Cannot close "//Trim(path)//'/rttov_test.txt')


    Open( file_id_nml, file = Trim(path)//'/../in/coef.txt', iostat = err)
    THROWM(err.ne.0,"Cannot open "//Trim(path)//'/../in/coef.txt')
    Read( file_id_nml, nml = coef_nml, iostat = err)
    THROWM(err.ne.0,"Cannot read "//Trim(path)//'/../in/coef.txt')
    Close( file_id_nml, iostat = err )
    THROWM(err.ne.0,"Cannot close "//Trim(path)//'/../in/coef.txt')

    call rttov_put_lun (file_id_nml)

    defn%opts%rt_ir%ozone_data = defn%f_o3     .ne. ""
    defn%opts%rt_ir%co2_data   = defn%f_co2    .ne. ""
    defn%opts%rt_ir%n2o_data   = defn%f_n2o    .ne. ""
    defn%opts%rt_ir%co_data    = defn%f_co     .ne. ""
    defn%opts%rt_ir%ch4_data   = defn%f_ch4    .ne. ""
    defn%opts%rt_mw%clw_data   = defn%f_clw    .ne. ""
    defn%opts%rt_ir%user_aer_opt_param = defn%f_aer_opt_param .ne. ""
    defn%opts%rt_ir%addaerosl  = defn%f_aerosl .ne. "" .or. defn%opts%rt_ir%user_aer_opt_param
    defn%opts%rt_ir%user_cld_opt_param = defn%f_cld_opt_param .ne. ""
    defn%opts%rt_ir%addclouds  = (defn%f_cloud .ne. "" .or. defn%opts%rt_ir%user_cld_opt_param) &
                                  .and. defn%f_cfrac  .ne. ""
    defn%opts%rt_ir%pc%addpc     = defn%f_pcscores     .ne. ""
    defn%opts%rt_ir%pc%addradrec = defn%f_channels_rec .ne. "" .and. defn%opts%rt_ir%pc%addpc

    ! Switchrad doesn't apply for PC-RTTOV unless addradrec is true: this avoids check_opts reporting an error
    if (defn%opts%rt_ir%pc%addpc .and. .not. defn%opts%rt_ir%pc%addradrec) defn%opts%rt_all%switchrad = .false.

    nprofiles1 = defn%nprofiles
    nchannels1 = defn%nchannels

    defn%nprofiles = defn%mult * defn%nprofiles
    defn%nchannels = defn%mult * defn%nchannels

    ts%defn = defn
    ts%path = path
    
    if(ts%defn%do_taylor) then
      if( ts%defn%do_ad   .or. &
          ts%defn%do_k    .or. &
          ts%defn%do_k_tl .or. &
          ts%defn%do_k_ad .or. &
          ts%defn%do_k_bf ) then
        err = errorstatus_fatal
        THROWM(err.ne.0,"Taylor test must be run alone")
      endif
      ts%defn%do_direct = .true. 
      ts%defn%do_tl     = .true.
    endif
    
! read input data

    if (ts%defn%opts%rt_ir%pc%addpc) then
      if (ts%defn%opts%rt_ir%pc%addradrec) then
        allocate(ts%data%channels_rec(ts%defn%nchannels_rec), stat = err)
        THROW(err.ne.0)
        call reada( err, Trim(ts%path)//'/../in/channels_rec.txt',  n1 = ts%data%channels_rec, foe = .true._jplm )
        THROW(err.ne.0)
      endif
    endif

    allocate(                                         &
      ts%data%chanprof( defn%nchannels ),             &
      ts%data%calcemis( defn%nchannels ),             &
      ts%data%calcrefl( defn%nchannels ),             &
      kchannels( nchannels1 ),                        &
      stat = err )
    THROWM(err.ne.0,"Cannot allocate ts%data")

    call reada( err, Trim(ts%path)//'/../in/channels.txt',  n1 = kchannels,                       foe = .true._jplm )
    THROW(err.ne.0)
    call reada( err, Trim(ts%path)//'/../in/lprofiles.txt', n1 = ts%data%chanprof(1:nchannels1)%prof, foe = .true._jplm )
    THROW(err.ne.0)


    if (ts%defn%opts%rt_ir%pc%addpc) then
      call reada( err, Trim(ts%path)//'/../in/pcscores.txt', npcscores = ts%defn%npcscores, &
                  ipcreg = ts%defn%opts%rt_ir%pc%ipcreg, ipcbnd = ts%defn%opts%rt_ir%pc%ipcbnd, foe = .true._jplm )
      THROW(err.ne.0)
    endif


! Channel indices calculations

    ! kchannels(:) holds the full list of channels from channels.txt
    ! lchannels(:) is of size MAXVAL(kchannels): for every channel appearing
    !   in channels.txt, the corresponding entry in lchannels is set non-zero.
    !   This provides a lookup mapping actual channel numbers to the extracted
    !   channel indices.
    ! ichannelsmax is the number of channels to extract
    ! ichannels(:) contains the list of channel numbers to extract

    lchannelsmax = MAXVAL(kchannels)
    ALLOCATE(lchannels(lchannelsmax), stat = err)
    THROWM(err.ne.0,"Cannot allocate lchannels")

    ! First set elements of lchannels to 1 for each channel present
    lchannels = 0
    DO i = 1, nchannels1
      lchannels(kchannels(i)) = 1
    ENDDO

    ! This gives the number of channels to extract
    ichannelsmax = SUM(lchannels)

    ! Now renumber the requested channels from 1 to ichannelsmax
    ! By doing it this way we allow for channels.txt to contain channels
    ! in non-consecutive order (e.g. chan 3 for prof 1, chan 2 for prof 2)
    j = 0
    DO i = 1, lchannelsmax
      IF (lchannels(i) > 0) THEN
        j = j + 1
        lchannels(i) = j
      ENDIF
    ENDDO

    ! Now populate the list of channels to extract by looping over
    ! lchannels and setting ichannel for every non-zero value
    ALLOCATE(ichannels(ichannelsmax), stat = err)
    THROWM(err.ne.0,"Cannot allocate ichannels")
    j = 0
    DO i = 1, lchannelsmax
      IF (lchannels(i) > 0) THEN
        j = j + 1
        ichannels(j) = i
      ENDIF
    ENDDO

!
!
! read coefficients


    If( ts%defn%lallchans ) Then
      Call rttov_read_coefs( &
        & err, ts%data%coefs, ts%defn%opts, &
        & file_coef   = Trim(ts%defn%coef_prefix)//'/'//Trim(ts%defn%f_coef),      &
        & file_scaer  = Trim(ts%defn%coef_prefix)//'/'//Trim(ts%defn%f_scat),      &
        & file_sccld  = Trim(ts%defn%coef_prefix)//'/'//Trim(ts%defn%f_optp),      &
        & file_pccoef = Trim(ts%defn%coef_prefix)//'/'//Trim(ts%defn%f_pccomp))

      THROW(err.ne.0)
    Else
      Call rttov_read_coefs( &
        & err, ts%data%coefs, ts%defn%opts, &
        & channels     = ichannels(1:ichannelsmax), &
        & channels_rec = ts%data%channels_rec,                                      &
        & file_coef    = Trim(ts%defn%coef_prefix)//'/'//Trim(ts%defn%f_coef),      &
        & file_scaer   = Trim(ts%defn%coef_prefix)//'/'//Trim(ts%defn%f_scat),      &
        & file_sccld   = Trim(ts%defn%coef_prefix)//'/'//Trim(ts%defn%f_optp),      &
        & file_pccoef  = Trim(ts%defn%coef_prefix) //'/'//Trim(ts%defn%f_pccomp))

      THROW(err.ne.0)
    EndIf
  
    If( ts%defn%lallchans ) Then
      ts%data%chanprof(1:nchannels1)%chan = kchannels
    Else
!
! Renumber all channels
!
      ts%data%chanprof(1:nchannels1)%chan = lchannels(kchannels)
      If (ts%defn%opts%rt_ir%pc%addpc) Then
        If (ts%defn%opts%rt_ir%pc%addradrec) Then
          Do i = 1, ts%defn%nchannels_rec
            ts%data%channels_rec( i ) = i
          EndDo
        EndIf
      EndIf
    EndIf

    Deallocate( kchannels, lchannels, ichannels, stat = err )
    THROWM(err.ne.0,"Deallocation of kchannels, lchannels, ... failed")

! alloc data
  
    call alloc_data( err, ts, ts%defn%nprofiles, ts%data%direct, do_taylor=ts%defn%do_taylor )
      
    If( ts%defn%do_tl    ) Then
      call alloc_data( err, ts, ts%defn%nprofiles, ts%data%tl )
      THROW(err.ne.0)
    endif
    If( ts%defn%do_ad    ) Then
      call alloc_data( err, ts, ts%defn%nprofiles, ts%data%ad,   do_ad=.true._jplm )
      THROW(err.ne.0)
    endif
    If( ts%defn%do_k     ) Then
      call alloc_data( err, ts, ts%defn%nchannels, ts%data%k,    ts%data%data_k )
      THROW(err.ne.0)
    endif
    If( ts%defn%do_k_bf  ) Then
      call alloc_data( err, ts, ts%defn%nchannels, ts%data%k_bf, ts%data%data_k_bf )
      THROW(err.ne.0)
    endif
    If( ts%defn%do_k_tl  ) Then
      call alloc_data( err, ts, ts%defn%nchannels, ts%data%k_tl, ts%data%data_k_tl )
      THROW(err.ne.0)
    endif
    If( ts%defn%do_k_ad  ) Then
      call alloc_data( err, ts, ts%defn%nchannels, ts%data%k_ad, ts%data%data_k_ad )
      THROW(err.ne.0)
    endif
   

!
! emissivity
!
    ts%data%calcemis = .FALSE.
    ts%data%direct%emissivity%emis_in = 1.0_jprb
    inquire( file=Trim(ts%path)//'/../in/calcemis.txt', exist=file_exists )
    if ( file_exists ) then
      call reada( err, Trim(ts%path)//'/../in/calcemis.txt',    l1 = ts%data%calcemis(1:nchannels1) )
      THROW(err.ne.0)
    endif
    inquire( file=Trim(ts%path)//'/../in/emissivity.txt', exist=file_exists )
    if ( file_exists ) then
      Call reada( err, Trim(ts%path)//'/../in/emissivity.txt',  r1 = ts%data%direct%emissivity(1:nchannels1)%emis_in )
      THROW(err.ne.0)
    endif

!
! reflectance
!
    ts%data%calcrefl = .TRUE.
    ts%data%direct%reflectance%refl_in = 0.3_jprb * pi_r
    ts%data%direct%reflectance%refl_cloud_top = 0._jprb   ! If zero, RTTOV uses default values
    inquire( file=Trim(ts%path)//'/../in/calcrefl.txt', exist=file_exists )
    if ( file_exists ) then
      call reada( err, Trim(ts%path)//'/../in/calcrefl.txt', l1 = ts%data%calcrefl(1:nchannels1) )
      THROW(err.ne.0)
    endif
    inquire( file=Trim(ts%path)//'/../in/reflectance.txt', exist=file_exists )
    if ( file_exists ) then
      Call reada( err, Trim(ts%path)//'/../in/reflectance.txt', r1 = ts%data%direct%reflectance(1:nchannels1)%refl_in )
      THROW(err.ne.0)
    endif
    inquire( file=Trim(ts%path)//'/../in/reflectance_cloud_top.txt', exist=file_exists )
    if ( file_exists ) then
      Call reada( err, Trim(ts%path)//'/../in/reflectance_cloud_top.txt', &
                  r1 = ts%data%direct%reflectance(1:nchannels1)%refl_cloud_top )
      THROW(err.ne.0)
    endif

!
! replicate emissivity and reflectance data
!
    Do imult = 2, ts%defn%mult
      k1 = 1 + ( imult - 1 ) * nchannels1
      k2 = k1 + nchannels1 - 1
      ts%data%direct%emissivity(k1:k2)%emis_in         = ts%data%direct%emissivity(1:nchannels1)%emis_in
      ts%data%direct%reflectance(k1:k2)%refl_in        = ts%data%direct%reflectance(1:nchannels1)%refl_in
      ts%data%direct%reflectance(k1:k2)%refl_cloud_top = ts%data%direct%reflectance(1:nchannels1)%refl_cloud_top
      ts%data%chanprof(k1:k2)%chan                     = ts%data%chanprof(1:nchannels1)%chan
      ts%data%chanprof(k1:k2)%prof                     = ts%data%chanprof(1:nchannels1)%prof + (imult - 1) * nprofiles1
      ts%data%calcemis(k1:k2)                          = ts%data%calcemis(1:nchannels1)
      ts%data%calcrefl(k1:k2)                          = ts%data%calcrefl(1:nchannels1)
    EndDo

    Call read_profile_data( ts, err )
    THROW(err.ne.0)

!
! Check pressure levels to determine whether interpolation is required
!

    ! User can set addinterp on commandline (false by default) and this is read
    ! into ts%defn%opts%interpolation%addinterp, but is overridden if profile
    ! pressure levels differ sufficiently from coef levels (tolerance is more
    ! strict than in rttov_checkinput to ensure we turn interpolator on for
    ! certain tests which require it).

    Call rttov_get_lun (file_id)
    interp_filename = Trim(ts%path)//'/interpolation.log'
    Open( file_id, file = interp_filename, form = 'formatted', iostat = err )
    THROWM(err.ne.0,"Cannot open "//Trim(interp_filename))

    interp_msg = 'Interpolation is OFF'
    If (ts%defn%opts%interpolation%addinterp) Then
      interp_msg = 'User has set interpolation ON'
    Else
      Do i = 1, nprofiles1
        If( ts%data%direct%profiles(i)%nlevels /= size(ts%data%coefs%coef%ref_prfl_p(:)) ) Then
          If (.Not. ts%defn%opts%interpolation%addinterp) &
            interp_msg = 'Number of input levels differs to coef file: switching interpolation ON'
          ts%defn%opts%interpolation%addinterp = .TRUE.
          exit
        Else If( Any(Abs(ts%data%direct%profiles(i)%p(:) - ts%data%coefs%coef%ref_prfl_p(:)) / &
                                   ts%data%coefs%coef%ref_prfl_p(:) > 0.005_jprb) ) Then
          If (.Not. ts%defn%opts%interpolation%addinterp) &
            interp_msg = 'Input pressures differ to coef levels: switching interpolation ON'
          ts%defn%opts%interpolation%addinterp = .TRUE.
          exit
        EndIf
      EndDo
    EndIf

    Write( file_id, '(a)' ) Trim(interp_msg)

    ! lgradp is only valid/useful if interpolation is switched on
    If(.Not. ts%defn%opts%interpolation%addinterp .and. ts%defn%opts%interpolation%lgradp)Then
      Write( file_id, '(a)' ) 'Interpolation is off: switching lgradp off'
      ts%defn%opts%interpolation%lgradp = .FALSE.
    End If

    Close( file_id, iostat = err )
    THROWM(err.ne.0,"Cannot close "//Trim(interp_filename))
    Call rttov_put_lun (file_id)


    If( ts%defn%ltemp ) Then
      Call alloc_traj( err, ts, 1_jpim )
      THROW(err.ne.0)
    endif

    
    If( ts%defn%do_tl ) Then
      call rttov_make_profile_inc( ts%data%tl%profiles, ts%data%direct%profiles, ts%defn%opts )

      ts%data%tl%emissivity(:)%emis_in = 0._jprb
      ts%data%tl%reflectance(:)%refl_in = 0._jprb
      Do i = 1, ts%defn%nchannels
        j = ts%data%chanprof(i)%prof
        chan = ts%data%chanprof(i)%chan
        
        If (ts%data%coefs%coef%id_sensor == sensor_id_mw .or. &
            ts%data%coefs%coef%id_sensor == sensor_id_po) Then
          
          ! For MW sensors rttov_calcemis_mw_tl calculates the TL if calcemis(:) is true
          If (.not. ts%data%calcemis(i)) Then
            ts%data%tl%emissivity(i)%emis_in = -0.01_jprb * (Modulo(j, 3_jpim) + 1_jpim)
          End If
          
        Else
          If (Associated(ts%data%coefs%coef%ss_val_chn)) Then
            ss_val = ts%data%coefs%coef%ss_val_chn(chan)
          Else
            ss_val = 0_jpim
          EndIf

          ! For IR sensors rttov_calcemis_ir_tl only calculates the TL if addpc is true 
          ! and the surface type is sea. Exclude pure-solar channels.
          If (.not. (ss_val == 2_jpim .or. (ts%defn%opts%rt_ir%pc%addpc .and. &
                       ts%data%direct%profiles(j)%skin%surftype == surftype_sea))) Then
            ts%data%tl%emissivity(i)%emis_in = -0.01_jprb * (Modulo(j, 3_jpim) + 1_jpim)
          End If

          ! For IR sensors, rttov_calcsurfrefl_tl only calculates the TL for sea surfaces 
          ! when calcrefl is true. However, over land and sea-ice, for thermal+solar channels, 
          ! if there is a valid emissivity it will use this to calculate a reflectance_tl.
          ! If the input reflectance is set to zero, the reflectance_tl is set to zero.

          If (.not. (ts%data%calcrefl(i) .and. &
                      (ts%data%direct%profiles(j)%skin%surftype == surftype_sea .or. &
                      (ts%data%direct%profiles(j)%skin%surftype /= surftype_sea .and. &
                      ss_val < 2_jpim .and. &
                      (ts%data%calcemis(i) .or. ts%data%direct%emissivity(i)%emis_in > 0._jprb))))) Then
            If (ts%data%calcrefl(i) .or. ts%data%direct%reflectance(i)%refl_in > 0._jprb) Then
              ts%data%tl%reflectance(i)%refl_in = 0.01_jprb * (Modulo(j, 3_jpim) + 1_jpim)
            End If
          End If

        End If

      End Do

      If( ts%defn%scale_inc /= 1.0_jprb ) Then
        Call rttov_scale_profile_inc( ts%data%tl%profiles, ts%defn%scale_inc )
        ts%data%tl%emissivity(:)%emis_in = ts%data%tl%emissivity(:)%emis_in * ts%defn%scale_inc
        ts%data%tl%reflectance(:)%refl_in = ts%data%tl%reflectance(:)%refl_in * ts%defn%scale_inc
      End If
    End If

    If( ts%defn%do_ad ) Then
      If( ts%defn%opts%rt_ir%pc%addpc )Then
        Call rttov_make_pccomp_inc( ts%data%ad%pccomp_saved, ts%defn%opts ) 
        If( ts%defn%scale_inc /= 1.0_jprb ) &
          Call rttov_scale_pccomp_inc( ts%data%ad%pccomp_saved, ts%defn%scale_inc, & 
         &                   ts%defn%opts ) 
      Else
        call rttov_make_radiance_inc( ts%data%coefs%coef, ts%data%ad%radiance_saved, &
           ts%data%chanprof(:)%chan, ts%defn%nchannels, ts%defn%opts )
    
        If( ts%defn%scale_inc /= 1.0_jprb ) &
          Call rttov_scale_radiance_inc( ts%data%ad%radiance_saved, &
            ts%defn%scale_inc, ts%defn%opts )
      End If
    End If
    
!     If (ts%data%coefs%coef%inczeeman) Then
!       print *,"Zeeman coefficient file with    Be   cosbk"
!       Do i = 1, nprofiles1
!         print "(24x,i4,2f8.4)",i,ts%data%direct%profiles(i)%be,ts%data%direct%profiles(i)%cosbk
!       EndDo
!     EndIf

    If( ts%data%coefs%coef%pmc_shift )Then
      inquire( file=Trim(ts%path)//'/../in/pmc.txt', exist=file_exists )
      If( file_exists )Then
        Call rttov_get_lun (file_id )
        Open( file_id, file = Trim(path)//'/../in/pmc.txt', iostat = err)
        THROWM(err.ne.0,"Cannot open "//Trim(path)//'/../in/pmc.txt')
        Read( file_id, *, iostat = err) ts%data%coefs%coef%pmc_ppmc
        THROWM(err.ne.0,"Cannot read "//Trim(path)//'/../in/pmc.txt')
        Close( file_id, iostat = err )
        THROWM(err.ne.0,"Cannot close "//Trim(path)//'/../in/pmc.txt')
      Else
        err = errorstatus_fatal
        THROWM(err.ne.0,'PMC coef file requires input pmc.txt file')
      End If
    End If

CATCH

  End Subroutine  

  Subroutine read_profile_data( ts, err )
  
    Implicit None
    
    Type(rttov_test_struct), Intent(inout), Target :: ts
    Integer(Kind=jpim),      Intent(out)           :: err
!
    Real(Kind=jprb)  :: becosbk(2)
    Integer(Kind=jpim) :: datetime(6)
    Type(rttov_test_defn), Pointer :: defn
    Type(rttov_input_data), Pointer :: direct
    Integer(Kind=jpim) :: iprof
    Integer(Kind=jpim) :: kprof
    Integer(Kind=jpim) :: imult
    Character(Len=6)   :: sproffmt
    Character(Len=8)   :: sprof
    Character(Len=256) :: in_dir
    Integer(Kind=jpim) :: nprofiles1
    Logical(Kind=jplm) :: foe
    Logical(Kind=jplm) :: file_exists
    Integer(Kind=jpim) :: file_id
    Logical(Kind=jplm) :: gas_conv_flag
    Character(Len=256) :: gas_units_filename
!
    Type(rttov_options) :: opts
#ifdef _RTTOV_HDF
#include "rttov_hdf_save.interface"
#endif

TRY


    defn => ts%defn
    direct => ts%data%direct

    nprofiles1 = ts%defn%nprofiles / ts%defn%mult



    Do iprof = 1, ts%defn%nprofiles / ts%defn%mult


      foe = iprof .eq. 1

      write( sproffmt, "('(i',i1,'.',i1,')')" ) ts%defn%prof_ndigits, ts%defn%prof_ndigits
      write( sprof, sproffmt ) iprof
      in_dir =  Trim(ts%path)//'/../in/profiles/'//Trim(sprof)


      call reada( err, Trim(in_dir)//'/ground/skin.txt', k0 = direct%profiles(iprof)%skin, &
                  foe = foe, k0_a = direct%profiles(1)%skin )
      THROW(err.ne.0)
      call reada( err, Trim(in_dir)//'/ground/s2m.txt',  s0 = direct%profiles(iprof)%s2m,  &
                  foe = foe, s0_a = direct%profiles(1)%s2m )
      THROW(err.ne.0)
      call reada( err, Trim(in_dir)//'/ground/elevation.txt', &
                  elevation = direct%profiles(iprof)%elevation, &
                  foe = foe, elevation_a = direct%profiles(1)%elevation )
      THROW(err.ne.0)
      call reada( err, Trim(in_dir)//'/angles.txt', &
        zenangle    = direct%profiles(iprof)%zenangle,    azangle    = direct%profiles(iprof)%azangle,    &
        sunzenangle = direct%profiles(iprof)%sunzenangle, sunazangle = direct%profiles(iprof)%sunazangle, &
        latitude    = direct%profiles(iprof)%latitude,                                                    &
        foe = foe,                          &
        zenangle_a    = direct%profiles(1)%zenangle,    azangle_a    = direct%profiles(1)%azangle,    &
        sunzenangle_a = direct%profiles(1)%sunzenangle, sunazangle_a = direct%profiles(1)%sunazangle, &
        latitude_a    = direct%profiles(1)%latitude                                                   &
         )
      THROW(err.ne.0)
      call reada( err, Trim(in_dir)//'/atm/cloud0.txt', &
        ctp = direct%profiles(iprof)%ctp, cfraction = direct%profiles(iprof)%cfraction, &
        ctp_a = direct%profiles(1)%ctp, cfraction_a = direct%profiles(1)%cfraction, &
        foe = foe )
      THROW(err.ne.0)
      call reada( err, Trim(in_dir)//'/atm/aerosli.txt', &
                  idg = direct%profiles(iprof)%idg, ish = direct%profiles(iprof)%ish, &
                  idg_a = direct%profiles(1)%idg, ish_a = direct%profiles(1)%ish, &
                  foe = foe )
      THROW(err.ne.0)
      Call reada( err, Trim(in_dir)//'/atm/p.txt',  r1 = direct%profiles(iprof)%p, foe = foe )
      THROW(err.ne.0)
      Call reada( err, Trim(in_dir)//'/atm/t.txt',  r1 = direct%profiles(iprof)%t, foe = foe )
      THROW(err.ne.0)
      Call reada( err, Trim(in_dir)//'/atm/q.txt',  r1 = direct%profiles(iprof)%q, foe = foe )
      THROW(err.ne.0)

      direct%profiles(iprof)%nlevels     = ts%defn%nlevels

      becosbk(1) = direct%profiles(1)%be  
      becosbk(2) = direct%profiles(1)%cosbk

      Call reada( err, Trim(in_dir)//'/be.txt', r1 = becosbk, foe = foe )
      THROW(err.ne.0)

      direct%profiles(iprof)%be    = becosbk(1)
      direct%profiles(iprof)%cosbk = becosbk(2)


      If( defn%opts%rt_ir%ozone_data ) then
        Call reada( err, Trim(in_dir)//'/atm/o3.txt',     &
          r1 = direct%profiles(iprof)%o3,  foe = foe )
        THROW(err.ne.0)
      endif
      If( defn%opts%rt_ir%co2_data   ) then
        Call reada( err, Trim(in_dir)//'/atm/co2.txt',    &
          r1 = direct%profiles(iprof)%co2, foe = foe ) 
        THROW(err.ne.0)
      endif
      If( defn%opts%rt_ir%n2o_data   ) then
        Call reada( err, Trim(in_dir)//'/atm/n2o.txt',    &
          r1 = direct%profiles(iprof)%n2o, foe = foe ) 
        THROW(err.ne.0)
      endif
      If( defn%opts%rt_ir%co_data    ) then
        Call reada( err, Trim(in_dir)//'/atm/co.txt',     &
          r1 = direct%profiles(iprof)%co,  foe = foe ) 
        THROW(err.ne.0)
      endif
      If( defn%opts%rt_ir%ch4_data   ) then
        Call reada( err, Trim(in_dir)//'/atm/ch4.txt',    &
          r1 = direct%profiles(iprof)%ch4, foe = foe ) 
        THROW(err.ne.0)
      endif
      If( defn%opts%rt_mw%clw_data   ) then
        Call reada( err, Trim(in_dir)//'/atm/clw.txt',    &
          r1 = direct%profiles(iprof)%clw, foe = foe ) 
        THROW(err.ne.0)
      endif
      If( defn%opts%rt_ir%addaerosl .and. .not. defn%opts%rt_ir%user_aer_opt_param ) then
        Call reada( err, Trim(in_dir)//'/atm/aerosl.txt', &
          r2 = direct%profiles(iprof)%aerosols, foe = foe )  
        THROW(err.ne.0)
      endif
      If( defn%opts%rt_ir%addclouds ) Then
        Call reada( err, Trim(in_dir)//'/atm/cfrac.txt', &
          r1 = direct%profiles(iprof)%cfrac, foe = foe )  
        THROW(err.ne.0)
        If ( .not. defn%opts%rt_ir%user_cld_opt_param ) Then
          Call reada( err, Trim(in_dir)//'/atm/cloud.txt', &
            r2 = direct%profiles(iprof)%cloud, foe = foe )  
          THROW(err.ne.0)
          Inquire(file=Trim(in_dir)//'/atm/icede.txt', exist=file_exists)
          If (file_exists) Then
            Call reada( err, Trim(in_dir)//'/atm/icede.txt', &
              r1 = direct%profiles(iprof)%icede, foe = foe )  
            THROW(err.ne.0)
          End If
        End If
      EndIf

      Inquire(file=Trim(in_dir)//'/datetime.txt', exist=file_exists)
      If (file_exists) Then
        Call reada( err, Trim(in_dir)//'/datetime.txt', &
          n1 = datetime(:), foe = foe )
        THROW(err.ne.0)
        direct%profiles(iprof)%date(1:3) = datetime(1:3)
        direct%profiles(iprof)%time(1:3) = datetime(4:6)
      End If

      ! Sort out gas units:
      !   defn%input_gas_units: the units of the input gas profiles in q.txt/etc
      !   defn%gas_units: the gas units with which RTTOV is to be run
      !   Default value for both is ppmv wet

      ! input_gas_units: valid values 0 (ppmv dry), 1 (kg/kg wet), 2 (ppmv wet)
      ! run_gas_units:   valid values -1 (ppmv dry), 0 (compat. mode), 1 (kg/kg wet), 2 (ppmv wet)
      ! NB compatibility mode treats profiles similar to ppmv dry

      ! If defn%input_gas_units was specified on commandline that value takes precedence;
      !   otherwise if gas_units.txt present then this is defn%input_gas_units;
      !   otherwise defn%input_gas_units = ppmv wet

      If (defn%input_gas_units == -1) Then
        Inquire(file=Trim(in_dir)//'/gas_units.txt', exist=file_exists)
        If (file_exists) Then
          call reada( err, Trim(in_dir)//'/gas_units.txt', &
                      gas_units = defn%input_gas_units, &
                      foe = .FALSE._jplm, gas_units_a = gas_unit_ppmv )
          THROW(err.ne.0)
        Else
          defn%input_gas_units = gas_unit_ppmv
        EndIf
      EndIf

      ! Default for gas_units is ppmv wet: all profiles get the same value
      direct%profiles(iprof)%gas_units = defn%run_gas_units

      ! If input_gas_units and run_gas_units differ do conversion
      gas_conv_flag = .TRUE.
      If (defn%input_gas_units == gas_unit_ppmv) Then
        If (direct%profiles(iprof)%gas_units <= gas_unit_compatibility) Then
          ! ppmv wet -> ppmv dry
          Call ppmvwet2ppmvdry(direct%profiles(iprof))
        ElseIf (direct%profiles(iprof)%gas_units == gas_unit_specconc) Then
          ! ppmv wet -> kgkg wet
          Call ppmvwet2kgkgwet(direct%profiles(iprof))
        Else
          gas_conv_flag = .FALSE.
        EndIf
      ElseIf (defn%input_gas_units == gas_unit_specconc) Then
        If (direct%profiles(iprof)%gas_units <= gas_unit_compatibility) Then
          ! kgkg wet -> ppmv dry
          Call kgkgwet2ppmvdry(direct%profiles(iprof))
        ElseIf (direct%profiles(iprof)%gas_units == gas_unit_ppmv) Then
          ! kgkg wet -> ppmv wet
          Call kgkgwet2ppmvwet(direct%profiles(iprof))
        Else
          gas_conv_flag = .FALSE.
        EndIf
      Else
        If (direct%profiles(iprof)%gas_units == gas_unit_ppmv) Then
          ! ppmv dry -> ppmv wet
          Call ppmvdry2ppmvwet(direct%profiles(iprof))
        ElseIf (direct%profiles(iprof)%gas_units == gas_unit_specconc) Then
          ! ppmv dry -> kgkg wet
          Call ppmvdry2kgkgwet(direct%profiles(iprof))
        Else
          gas_conv_flag = .FALSE.
        EndIf
      EndIf

      If (iprof == 1) Then
        Call rttov_get_lun (file_id)
        gas_units_filename = Trim(ts%path)//'/gas_units.log'
        Open( file_id, file = gas_units_filename, form = 'formatted', iostat = err )
        THROWM(err.ne.0,"Cannot open "//Trim(gas_units_filename))

        If (defn%input_gas_units == gas_unit_ppmv) Then
          Write( file_id, '(a)' ) Trim('Input files have gas units  : ppmv over moist air')
        Else If (defn%input_gas_units == gas_unit_specconc) Then
          Write( file_id, '(a)' ) Trim('Input files have gas units  : kg/kg over moist air')
        Else
          Write( file_id, '(a)' ) Trim('Input files have gas units  : ppmv over dry air')
        EndIf
        If (defn%run_gas_units == gas_unit_ppmv) Then
          Write( file_id, '(a)' ) Trim('Running RTTOV with gas units: ppmv over moist air')
        Else If (defn%run_gas_units == gas_unit_specconc) Then
          Write( file_id, '(a)' ) Trim('Running RTTOV with gas units: kg/kg over moist air')
        Else If (defn%run_gas_units == gas_unit_compatibility) Then
          Write( file_id, '(a)' ) Trim('Running RTTOV with gas units: compatibility mode')
        Else
          Write( file_id, '(a)' ) Trim('Running RTTOV with gas units: ppmv over dry air')
        EndIf
        If (gas_conv_flag) Write( file_id, '(a)' ) 'Test suite converted units of input profiles'

        Close( file_id, iostat = err )
        THROWM(err.ne.0,"Cannot close "//Trim(gas_units_filename))
        Call rttov_put_lun (file_id)
      EndIf

!
! if multiplicity > 1, then replicate profile data
!
      Do imult = 2, ts%defn%mult
        kprof = iprof + (imult-1) * nprofiles1
        Call rttov_copy_prof( direct%profiles(kprof:kprof), direct%profiles(iprof:iprof) )
      EndDo

    EndDo

! initialise data structures and read the aer/cld optical parameter files
    If (defn%opts%rt_ir%user_aer_opt_param) Then
      Call read_opt_param( err, ts, Trim(ts%path)//'/../in/aer_opt_param.txt', ts%data%aer_opt_param )
      THROWM(err.ne.0,"Error initialising aerosol optical parameters")
    End If
    If (defn%opts%rt_ir%user_cld_opt_param) Then
      Call read_opt_param( err, ts, Trim(ts%path)//'/../in/cld_opt_param.txt', ts%data%cld_opt_param )
      THROWM(err.ne.0,"Error initialising cloud optical parameters")
    End If


#ifdef _RTTOV_HDF
! aer/cld params are saved for all channels and multipicity (shall we remove multiplicity?)
    If( ts%defn%savehdf5 ) Then
      If (defn%opts%rt_ir%user_aer_opt_param) Then
        CALL RTTOV_HDF_SAVE( ERR, TRIM(ts%path)//"/AER_OPT_PARAM.H5",'/USER_AEROSOL_OPT_PARAM', &
            &   CREATE=.true., OPT_PARAM = ts%data%aer_opt_param)
        THROW(err.ne.0)
        CALL RTTOV_HDF_SAVE( ERR, TRIM(ts%path)//"/AER_OPT_PARAM.H5",'/MISC', CREATE=.false., &
            &   c0 = ts%defn%f_coef, SNAME='COEF_FILENAME')
        THROW(err.ne.0)      
        CALL RTTOV_HDF_SAVE( ERR, TRIM(ts%path)//"/AER_OPT_PARAM.H5",'/CHANPROF', CREATE=.false., &
            &   CHANPROF = ts%data%chanprof)
        THROW(err.ne.0)
      End If

      If (defn%opts%rt_ir%user_cld_opt_param) Then
        CALL RTTOV_HDF_SAVE( ERR, TRIM(ts%path)//"/CLD_OPT_PARAM.H5",'/USER_CLOUD_OPT_PARAM', &
               CREATE=.true., OPT_PARAM = ts%data%cld_opt_param)
        THROW(err.ne.0)
        CALL RTTOV_HDF_SAVE( ERR, TRIM(ts%path)//"/CLD_OPT_PARAM.H5",'/MISC', CREATE=.false., &
            &   c0 = ts%defn%f_coef, SNAME='COEF_FILENAME')
        THROW(err.ne.0)      
        CALL RTTOV_HDF_SAVE( ERR, TRIM(ts%path)//"/CLD_OPT_PARAM.H5",'/CHANPROF', CREATE=.false., &
            &   CHANPROF = ts%data%chanprof)
        THROW(err.ne.0)
      End If

      CALL RTTOV_HDF_SAVE( ERR, TRIM(ts%path)//"/PROFILES.H5",  '/PROFILES', &
                 CREATE=.true., PROFILES = direct%profiles(1:nprofiles1))
      THROW(err.ne.0)

! opts is a local variable, so all fields have the default value
      opts%interpolation%addinterp = .true.
      opts%rt_ir%addaerosl  = defn%opts%rt_ir%addaerosl
      opts%rt_ir%addclouds  = defn%opts%rt_ir%addclouds
      opts%rt_ir%ozone_data = defn%opts%rt_ir%ozone_data
      opts%rt_ir%co2_data   = defn%opts%rt_ir%co2_data
      opts%rt_ir%n2o_data   = defn%opts%rt_ir%n2o_data
      opts%rt_ir%co_data    = defn%opts%rt_ir%co_data
      opts%rt_ir%ch4_data   = defn%opts%rt_ir%ch4_data
      opts%rt_mw%clw_data   = defn%opts%rt_mw%clw_data
      opts%rt_ir%user_aer_opt_param = defn%opts%rt_ir%user_aer_opt_param
      opts%rt_ir%user_cld_opt_param = defn%opts%rt_ir%user_cld_opt_param
      CALL RTTOV_HDF_SAVE( ERR, TRIM(ts%path)//"/PROFILES.H5",  '/OPTIONS', CREATE=.false., OPTIONS = opts)
      THROW(err.ne.0)
      EndIf
#endif

CATCH
  End Subroutine

  Subroutine Cleanup( ts, err )
!
    Implicit None
    
    Type(rttov_test_struct), Intent(inout), Target :: ts
    Integer(Kind=jpim),      Intent(out)           :: err
!
TRY

    deallocate( ts%data%chanprof, ts%data%calcemis, ts%data%calcrefl, stat = err )  
    THROWM(err.ne.0,"Deallocation of ts%data failed")

    if (ts%defn%opts%rt_ir%user_aer_opt_param) then
      call rttov_alloc_opt_param( err, ts%data%aer_opt_param, ts%defn%nchannels, &
                                  ts%defn%nlevels - 1_jpim, size(ts%data%aer_opt_param%phangle), 0_jpim )
      THROW(err.ne.0)
    endif
    if (ts%defn%opts%rt_ir%user_cld_opt_param) then
      call rttov_alloc_opt_param( err, ts%data%cld_opt_param, ts%defn%nchannels, &
                                  ts%defn%nlevels - 1_jpim, size(ts%data%cld_opt_param%phangle), 0_jpim )
      THROW(err.ne.0)
    endif

    if(ts%defn%opts%rt_ir%pc%addpc) then
     if(ts%defn%opts%rt_ir%pc%addradrec) then
       deallocate( ts%data%channels_rec, stat = err)
       THROW(err.ne.0)
     endif
    endif

    If( ts%defn%ltemp ) Then
      Call alloc_traj( err, ts, 0_jpim )
      THROW(err.ne.0)
    endif

    If( ts%defn%do_tl    ) Then
      call dealloc_data( err, ts, ts%defn%nprofiles, ts%data%tl )
      THROW(err.ne.0)
    endif
    If( ts%defn%do_ad    ) Then
      call dealloc_data( err, ts, ts%defn%nprofiles, ts%data%ad,   do_ad=.true._jplm )
      THROW(err.ne.0)
    endif
    If( ts%defn%do_k     ) Then
      call dealloc_data( err, ts, ts%defn%nchannels, ts%data%k,    ts%data%data_k )
      THROW(err.ne.0)
    endif
    If( ts%defn%do_k_bf  ) Then
      call dealloc_data( err, ts, ts%defn%nchannels, ts%data%k_bf, ts%data%data_k_bf )
      THROW(err.ne.0)
    endif
    If( ts%defn%do_k_tl  ) Then
      call dealloc_data( err, ts, ts%defn%nchannels, ts%data%k_tl, ts%data%data_k_tl )
      THROW(err.ne.0)
    endif
    If( ts%defn%do_k_ad  ) Then
      call dealloc_data( err, ts, ts%defn%nchannels, ts%data%k_ad, ts%data%data_k_ad )
      THROW(err.ne.0)
    endif

    call dealloc_data( err, ts, ts%defn%nprofiles, ts%data%direct, do_taylor=ts%defn%do_taylor )
    THROW(err.ne.0)
    
    call rttov_dealloc_coefs( err, ts%data%coefs )
    THROW(err.ne.0)

CATCH
    
  End Subroutine
  

  Subroutine reada( err, f,                                          &
    r1, r1_a, r2, r2_a, n1, n1_a, l1, l1_a, k0, k0_a, s0, s0_a,      &
    gas_units, gas_units_a,                                          &
    zenangle,   azangle,   sunzenangle,   sunazangle,   elevation,   latitude,   &
    zenangle_a, azangle_a, sunzenangle_a, sunazangle_a, elevation_a, latitude_a, &
    idg,   ish,                                                          &
    idg_a, ish_a,                                                        &
    ctp,   cfraction,                                                    &
    ctp_a, cfraction_a,                                                  &
    addinterp,                                                           &
    addinterp_a,                                                         &
    ipcreg, ipcbnd, npcscores,                                           &
    foe                                                                  &
  )
    
    Implicit None
    
    Integer(Kind=jpim), Intent(out)   :: err
    Character(len=*),   Intent(in)    :: f
    Real(Kind=jprb),    Intent(inout), Optional :: r1(:),   r2(:,:)
    Real(Kind=jprb),    Intent(in),    Optional :: r1_a(:), r2_a(:,:)
    Integer(kind=jpim), Intent(inout), Optional :: n1(:)
    Integer(kind=jpim), Intent(in),    Optional :: n1_a(:)
    Logical(Kind=jplm), Intent(inout), Optional :: l1(:)
    Logical(Kind=jplm), Intent(in),    Optional :: l1_a(:)
    Type(sskin_type),   Intent(inout), Optional :: k0
    Type(sskin_type),   Intent(in),    Optional :: k0_a
    Type(s2m_type),     Intent(inout), Optional :: s0
    Type(s2m_type),     Intent(in),    Optional :: s0_a
    Integer(Kind=jpim), Intent(inout), Optional :: gas_units
    Integer(Kind=jpim), Intent(in),    Optional :: gas_units_a
    Real(Kind=jprb),    Intent(inout), Optional :: zenangle,   azangle,   sunzenangle,   sunazangle,   elevation,   latitude
    Real(Kind=jprb),    Intent(in),    Optional :: zenangle_a, azangle_a, sunzenangle_a, sunazangle_a, elevation_a, latitude_a
    Real(Kind=jprb),    Intent(inout), Optional :: ctp,   cfraction
    Real(Kind=jprb),    Intent(in),    Optional :: ctp_a, cfraction_a
    Integer(Kind=jpim), Intent(inout), Optional :: idg,   ish
    Integer(Kind=jpim), Intent(in),    Optional :: idg_a, ish_a
    Logical(Kind=jplm), Intent(inout), Optional :: addinterp
    Logical(Kind=jplm), Intent(in),    Optional :: addinterp_a
    Integer(Kind=jpim), Intent(inout), Optional :: ipcreg
    Integer(Kind=jpim), Intent(inout), Optional :: ipcbnd
    Integer(Kind=jpim), Intent(inout), Optional :: npcscores
!
    Logical(Kind=jplm), Intent(in),    Optional :: foe ! fail on data not present
!
    Integer(Kind=jpim) :: ipcreg1, ipcbnd1, npcscores1
    Namelist / pcscores / ipcreg, ipcbnd, npcscores
    Namelist / flags / addinterp
    Namelist / units / gas_units
    Namelist / angles / zenangle, azangle, sunzenangle, sunazangle, latitude
    Namelist / s2m / s0
    Namelist / aerosli / idg, ish
    Namelist / cloud / ctp, cfraction
    Namelist / elev / elevation
    Namelist / skin / k0
!
    Logical(Kind=jplm) :: foe1
    Integer(Kind=jpim) :: file_id
!
TRY

    foe1 = .false.
    If( Present(foe) ) foe1 = foe
    

    call rttov_get_lun (file_id )
    If( foe1 ) Then
      Open( file_id, file = Trim( f ), form = 'formatted', status = 'old', iostat = err )
      THROWM(err.ne.0,"Cannot open "//Trim( f ))
    Else
      Open( file_id, file = Trim( f ), form = 'formatted', status = 'old', iostat = err )
      If( err .ne. 0 ) Then
        If( Present( r1 ) .and. Present( r1_a ) ) Then
          r1 = r1_a
        Else If( Present( r2 ) .and. Present( r2_a ) ) Then
          r2 = r2_a
        Else If( Present( l1 ) .and. Present( l1_a ) ) Then
          l1 = l1_a
        Else If( Present( n1 ) .and. Present( n1_a ) ) Then
          n1 = n1_a
        Else If( Present( k0 ) .and. Present( k0_a ) ) Then
          k0 = k0_a
        Else If( Present( gas_units ) .and. Present( gas_units_a ) ) Then
          gas_units = gas_units_a
        Else If( Present( elevation ) .and. Present( elevation_a ) ) Then
          elevation = elevation_a
        Else If( Present( s0 ) .and. Present( s0_a ) ) Then
          s0 = s0_a
        Else If( Present( zenangle ) .and. Present( zenangle_a ) ) Then
          zenangle    = zenangle_a
          azangle     = azangle_a
          sunzenangle = sunzenangle_a
          sunazangle  = sunazangle_a
          latitude    = latitude_a
        Else If( Present( ish ) .and. Present( ish_a ) ) Then
          ish = ish_a
          idg = idg_a
        Else If( Present( ctp ) .and. Present( ctp_a ) ) Then
          ctp       = ctp_a
          cfraction = cfraction_a
        Else If( Present( addinterp ) .and. Present( addinterp_a ) ) Then
          addinterp = addinterp_a
        Else
          THROWM(err.ne.0,"Cannot open "//Trim( f ))
        EndIf
        err = 0
      endif
    Endif
    
    
    If( Present( r1 ) ) Then
      Read( file_id, *, iostat = err ) r1
      THROWM(err.ne.0,"Cannot read from "//Trim( f ))
    Else If( Present( r2 ) ) Then
      Read( file_id, *, iostat = err ) r2
      THROWM(err.ne.0,"Cannot read from "//Trim( f ))
    Else If( Present( l1 ) ) Then
      Read( file_id, *, iostat = err ) l1
      THROWM(err.ne.0,"Cannot read from "//Trim( f ))
    Else If( Present( n1 ) ) Then
      Read( file_id, *, iostat = err ) n1
      THROWM(err.ne.0,"Cannot read from "//Trim( f ))
    Else If( Present( k0 ) ) Then
      Read( file_id, nml = skin, iostat = err ) 
      THROWM(err.ne.0,"Cannot read from "//Trim( f ))
    Else If( Present( gas_units ) ) Then
      Read( file_id, nml = units, iostat = err )
      THROWM(err.ne.0,"Cannot read from "//Trim( f ))
    Else If( Present( elevation ) ) Then
      Read( file_id, nml = elev, iostat = err ) 
      THROWM(err.ne.0,"Cannot read from "//Trim( f ))
    Else If( Present( s0 ) ) Then
      Read( file_id, nml = s2m, iostat = err ) 
      THROWM(err.ne.0,"Cannot read from "//Trim( f ))
    Else If( Present( zenangle ) ) Then
      Read( file_id, nml = angles, iostat = err ) 
      THROWM(err.ne.0,"Cannot read from "//Trim( f ))
    Else If( Present( ish ) ) Then
      Read( file_id, nml = aerosli, iostat = err ) 
      THROWM(err.ne.0,"Cannot read from "//Trim( f ))
    Else If( Present( ctp ) ) Then
      Read( file_id, nml = cloud, iostat = err ) 
      THROWM(err.ne.0,"Cannot read from "//Trim( f ))
    Else If( Present( addinterp ) ) Then
      Read( file_id, nml = flags, iostat = err ) 
      THROWM(err.ne.0,"Cannot read from "//Trim( f ))
    Else If( Present( ipcreg ) .and. Present( ipcbnd ) .and. Present( npcscores ) ) Then
      ipcreg1 = ipcreg
      ipcbnd1 = ipcbnd
      npcscores1 = npcscores
      Read( file_id, nml = pcscores, iostat = err )
      THROWM(err.ne.0,"Cannot read from "//Trim( f ))
      If( ipcreg1 .gt. 0 ) ipcreg = ipcreg1
      If( ipcbnd1 .gt. 0 ) ipcbnd = ipcbnd1
      If( npcscores1 .gt. 0 ) npcscores = npcscores1
      If( ipcbnd .lt. 0 ) ipcbnd = 1_jpim
    EndIf

    Close( file_id, iostat = err )
    THROWM(err.ne.0,"Cannot close "//Trim( f ))
    call rttov_put_lun (file_id)

CATCH
  End Subroutine

  Subroutine read_opt_param( err, ts, f, opt_param )
    Integer(Kind=jpim),      Intent(out)   :: err
    Type(rttov_test_struct), Intent(in)    :: ts
    Character(len=*),        Intent(in)    :: f
    Type(rttov_opt_param),   Intent(inout) :: opt_param

    Integer(Kind=jpim) :: file_id, imult, chanlo, chanhi, nchannels1, nphangle

    Call rttov_get_lun (file_id)
    Open( file_id, file = Trim(f), form = 'formatted', status = 'old', iostat = err )
    THROWM(err.ne.0,"Cannot open "//Trim(f))

    Read( file_id, *) nphangle
    Call rttov_alloc_opt_param( err, opt_param, ts%defn%nchannels, &
                                ts%defn%nlevels-1_jpim, nphangle, 1_jpim )
    THROWM(err.ne.0,"Error allocating opt_param")

    nchannels1 = ts%defn%nchannels / ts%defn%mult
    Read( file_id, *) opt_param%abs(1:nchannels1,:)
    Read( file_id, *) opt_param%sca(1:nchannels1,:)
    Read( file_id, *) opt_param%bpr(1:nchannels1,:)
    If (ts%defn%opts%rt_ir%addsolar) Read( file_id, *) opt_param%phangle(:)
    If (ts%defn%opts%rt_ir%addsolar) Read( file_id, *) opt_param%pha(1:nchannels1,:,:)

    Close( file_id, iostat = err )
    THROWM(err.ne.0,"Cannot close "//Trim(f))
    Call rttov_put_lun (file_id)

    Do imult = 2, ts%defn%mult
       chanlo = (imult - 1) * nchannels1 + 1
       chanhi = imult * nchannels1
       opt_param%abs(chanlo:chanhi,:)   = opt_param%abs(1:nchannels1,:)
       opt_param%sca(chanlo:chanhi,:)   = opt_param%sca(1:nchannels1,:)
       opt_param%bpr(chanlo:chanhi,:)   = opt_param%bpr(1:nchannels1,:)
       If (ts%defn%opts%rt_ir%addsolar) opt_param%pha(chanlo:chanhi,:,:) = opt_param%pha(1:nchannels1,:,:)
    End Do

    If (ts%defn%opts%rt_ir%addsolar) Then
      call rttov_init_opt_param( err, ts%defn%opts, opt_param )
      THROWM(err.ne.0,"Error initialising opt_param")
    End If
 
CATCH
  End Subroutine
  
  Subroutine alloc_data( err, ts, nprofiles, dat, dat_k, do_ad, do_taylor )
!
    Implicit None
    
    Integer(Kind=jpim),      Intent(out)   :: err
    Type(rttov_test_struct), Intent(in)    :: ts
    Integer(Kind=jpim),      Intent(in)    :: nprofiles
    Type(rttov_input_data),  Intent(inout) :: dat
    Type(rttov_k_input_data),Intent(inout), Optional :: dat_k
    Logical(Kind=jplm),      Intent(in),    Optional :: do_ad
    Logical(Kind=jplm),      Intent(in),    Optional :: do_taylor
!
    Logical(Kind=jplm) :: ldo_ad, ldo_taylor
!
TRY
    ldo_ad = .false._jplm
    if (present(do_ad)) ldo_ad = do_ad
    ldo_taylor = .false._jplm
    if (present(do_taylor)) ldo_taylor = do_taylor
    
    allocate(                                &
      dat%emissivity( ts%defn%nchannels ),   &
      dat%reflectance( ts%defn%nchannels ),  &
      dat%profiles( nprofiles ),  stat = err )
    THROWM(err.ne.0,"Cannot allocate dat")
    
    if( ldo_taylor ) then
      allocate( dat%profiles_pert( nprofiles ),                       &
                dat%emissivity_pert( ts%defn%nchannels ),             &
                dat%reflectance_pert( ts%defn%nchannels ), stat = err )
      THROWM(err.ne.0,"Cannot allocate dat")
    end if    
   
    dat%emissivity%emis_in      = 0._jprb
    dat%emissivity%emis_out     = 0._jprb
    dat%reflectance%refl_in     = 0._jprb
    dat%reflectance%refl_out    = 0._jprb
   
    call rttov_alloc_prof( err, nprofiles, dat%profiles, ts%defn%nlevels, ts%defn%opts, 1_jpim, & 
                           & init = .true._jplm, coefs = ts%data%coefs )
    THROW(err.ne.0)

    call rttov_alloc_rad( err, ts%defn%nchannels, dat%radiance, ts%defn%nlevels-1_jpim, 1_jpim, dat%radiance2, init = .true._jplm )
    THROW(err.ne.0)

    call rttov_alloc_transmission( err, dat%transmission, ts%defn%nlevels-1_jpim, ts%defn%nchannels, &
                                   & 1_jpim, init = .true._jplm )
    THROW(err.ne.0)
    
    if(ldo_taylor) then
      call rttov_alloc_prof( err, ts%defn%nprofiles, dat%profiles_pert, ts%defn%nlevels, ts%defn%opts, 1_jpim, & 
                           & init = .true._jplm, coefs = ts%data%coefs )
      THROW(err.ne.0)
      dat%emissivity_pert(:)%emis_in = 0._jprb
      dat%reflectance_pert(:)%refl_in = 0._jprb
      
      if(ts%defn%opts%rt_ir%pc%addpc)then
        if(ts%defn%opts%rt_ir%pc%addradrec)then
          call rttov_alloc_pccomp( err, dat%pccomp_saved, ts%defn%npcscores * ts%defn%nprofiles, 1_jpim, &
                                 & init = .true._jplm, nchannels_rec = ts%defn%nchannels_rec * ts%defn%nprofiles )
        else
          call rttov_alloc_pccomp( err, dat%pccomp_saved, ts%defn%npcscores * ts%defn%nprofiles, 1_jpim, &
                                 & init = .true._jplm )
        endif
      else
        call rttov_alloc_rad( err, ts%defn%nchannels, dat%radiance_saved, ts%defn%nlevels-1_jpim, 1_jpim )
      endif
      THROW(err.ne.0)
    end if
    
    if(ts%defn%opts%rt_ir%pc%addpc)then
      if(ts%defn%opts%rt_ir%pc%addradrec)then
        call rttov_alloc_pccomp( err, dat%pccomp, ts%defn%npcscores * ts%defn%nprofiles, 1_jpim, init = .true._jplm, &
                               & nchannels_rec = ts%defn%nchannels_rec * ts%defn%nprofiles )
      else
        call rttov_alloc_pccomp( err, dat%pccomp, ts%defn%npcscores * ts%defn%nprofiles, 1_jpim, init = .true._jplm )
      endif
      THROW(err.ne.0)
    endif

    if( ldo_ad ) then
      if( ts%defn%opts%rt_ir%pc%addpc ) then
        if(ts%defn%opts%rt_ir%pc%addradrec)then
          call rttov_alloc_pccomp( err, dat%pccomp_saved, ts%defn%npcscores * ts%defn%nprofiles, 1_jpim, &
                                 & init = .true._jplm, nchannels_rec = ts%defn%nchannels_rec * ts%defn%nprofiles )
        else
          call rttov_alloc_pccomp( err, dat%pccomp_saved, ts%defn%npcscores * ts%defn%nprofiles, 1_jpim, &
                                 & init = .true._jplm )
        endif
      else
        call rttov_alloc_rad( err, ts%defn%nchannels, dat%radiance_saved, ts%defn%nlevels-1_jpim, 1_jpim, &
                            & init=.true._jplm )
      end if
      THROW(err.ne.0)
    end if
    
    if(present(dat_k)) then
      if(ts%defn%opts%rt_ir%pc%addpc) then
        if(ts%defn%opts%rt_ir%pc%addradrec) then

          allocate(dat_k%profiles_k_rec( ts%defn%nchannels_rec * ts%defn%nprofiles ), stat = err)
          THROW(err.ne.0)
          call rttov_alloc_prof( err, ts%defn%nchannels_rec * ts%defn%nprofiles, dat_k%profiles_k_rec, ts%defn%nlevels, &
           & ts%defn%opts, 1_jpim, init = .true._jplm, coefs = ts%data%coefs )
          THROW(err.ne.0)

        else

          allocate(dat_k%profiles_k_pc( ts%defn%npcscores * ts%defn%nprofiles ), stat = err)
          THROW(err.ne.0)
          call rttov_alloc_prof( err, ts%defn%npcscores * ts%defn%nprofiles, &
           & dat_k%profiles_k_pc, &
           & ts%defn%nlevels, ts%defn%opts, 1_jpim, init = .true._jplm, coefs = ts%data%coefs )
          THROW(err.ne.0)

        endif
      endif

    endif

CATCH
     
  End Subroutine
  
  Subroutine dealloc_data( err, ts, nprofiles, dat, dat_k, do_ad, do_taylor )
!
    Implicit None
    
    Integer(Kind=jpim),      Intent(out)   :: err
    Type(rttov_test_struct), Intent(in)    :: ts
    Integer(Kind=jpim),      Intent(in)    :: nprofiles
    Type(rttov_input_data),  Intent(inout) :: dat
    Type(rttov_k_input_data),Intent(inout), Optional :: dat_k
    Logical(Kind=jplm),      Intent(in),    Optional :: do_ad
    Logical(Kind=jplm),      Intent(in),    Optional :: do_taylor
!
    Logical(Kind=jplm) :: ldo_ad, ldo_taylor
!
TRY
    ldo_ad = .false._jplm
    if (present(do_ad)) ldo_ad = do_ad
    ldo_taylor = .false._jplm
    if (present(do_taylor)) ldo_taylor = do_taylor

    call rttov_alloc_prof( err, nprofiles, dat%profiles, ts%defn%nlevels, ts%defn%opts, 0_jpim )
    THROW(err.ne.0)
    call rttov_alloc_rad( err, ts%defn%nchannels, dat%radiance, ts%defn%nlevels-1_jpim, 0_jpim, dat%radiance2 )
    THROW(err.ne.0)
    call rttov_alloc_transmission( err, dat%transmission, ts%defn%nlevels-1_jpim, ts%defn%nchannels, 0_jpim )
    THROW(err.ne.0)
    
    if (ldo_taylor) then
      call rttov_alloc_prof( err, nprofiles, dat%profiles_pert, ts%defn%nlevels, ts%defn%opts, 0_jpim )
      THROW(err.ne.0)
      if( ts%defn%opts%rt_ir%pc%addpc ) then
        call rttov_alloc_pccomp (err, dat%pccomp_saved, ts%defn%npcscores, 0_jpim)
      else
        call rttov_alloc_rad( err, ts%defn%nchannels, dat%radiance_saved, ts%defn%nlevels-1_jpim, 0_jpim )
      end if
      THROW(err.ne.0)
    end if
    
    deallocate(                                      &
      dat%emissivity,                                &
      dat%reflectance,                               &
      dat%profiles, stat = err )
    THROWM(err.ne.0,"Deallocation of dat failed")

    if( ldo_taylor ) then
      deallocate( dat%profiles_pert,               &
                  dat%emissivity_pert,             &
                  dat%reflectance_pert, stat = err )
      THROWM(err.ne.0,"Cannot deallocate dat")
    end if  
    
    if(ts%defn%opts%rt_ir%pc%addpc)then
      call rttov_alloc_pccomp (err, dat%pccomp, ts%defn%npcscores, 0_jpim)
      THROW(err.ne.0)
    endif

    if( ldo_ad ) then
      if( ts%defn%opts%rt_ir%pc%addpc ) then
        call rttov_alloc_pccomp (err, dat%pccomp_saved, ts%defn%npcscores, 0_jpim)
      else
        call rttov_alloc_rad( err, ts%defn%nchannels, dat%radiance_saved, ts%defn%nlevels-1_jpim, 0_jpim )
      end if
      THROW(err.ne.0)
    end if
    
    if(present(dat_k)) then
      if(ts%defn%opts%rt_ir%pc%addpc) then
        if(ts%defn%opts%rt_ir%pc%addradrec) then

          call rttov_alloc_prof( err, ts%defn%nchannels_rec * ts%defn%nprofiles, dat_k%profiles_k_rec, &
                                 ts%defn%nlevels, ts%defn%opts, 0_jpim)
          THROW(err.ne.0)

          deallocate(dat_k%profiles_k_rec, stat = err)
          THROW(err.ne.0)

        else

          call rttov_alloc_prof( err, ts%defn%npcscores * ts%defn%nprofiles, dat_k%profiles_k_pc, &
                                 ts%defn%nlevels, ts%defn%opts, 0_jpim )
          THROW(err.ne.0)

          deallocate(dat_k%profiles_k_pc, stat = err)
          THROW(err.ne.0)

        endif
      endif

    endif

CATCH
   
  End Subroutine

! Print routines  
  
  Subroutine PrintProfiles( err, profiles, name, f )
    Use rttov_chain, only : chain, chain_profile_type, print_chain, free_chain

    Implicit None

    Integer(Kind=jpim), Intent(out) :: err
    Type(profile_type), Intent(in), Target :: profiles(:) ! target because of chain_profile_type
    Character(len=*),   Intent(in) :: name, f
    
    Type(chain) :: chain_prof
    Integer(Kind=jpim) :: file_id

TRY
    
    call rttov_get_lun (file_id )
    Open( file_id, file = f, form = 'formatted', iostat = err )
    THROWM(err.ne.0,"Cannot open "//Trim(f))
    
    Call chain_profile_type( err, chain_prof, name, a1 = profiles )
    THROWM(err.ne.0,"Cannot chain profiles")
    
    Call print_chain( Int(file_id,jpim), chain_prof, Gformat )
    
    Call free_chain( chain_prof )

    Close( file_id, iostat = err )
    THROWM(err.ne.0,"Cannot close "//Trim(f))
    call rttov_put_lun (file_id)

CATCH

  End Subroutine

  Subroutine PrintRadiance (err, radiance, name, f, radiance2 )
    Use rttov_chain, only : print_array

    Implicit None

    Integer(Kind=jpim),   Intent(out) :: err
    Type(radiance_type),  Intent(in) :: radiance
    Character(len=*),     Intent(in) :: name, f
    Type(radiance2_type), Intent(in), Optional :: radiance2

    Integer(Kind=jpim) :: file_id

TRY

    call rttov_get_lun (file_id )
    Open( file_id, file = f, form = 'formatted', iostat = err )
    THROWM(err.ne.0,"Cannot open "//Trim(f))

    Write( file_id, '(A," = (")' ) name//'%TOTAL'
    Call Print_Array( Int(file_id,jpim), '(10000(' // Trim( Gformat ) // '))', A1 = radiance%total )
    Write( file_id, * ) ')'

    If (Any(radiance%bt /= 0._jprb)) Then
      Write( file_id, '(A," = (")' ) name//'%BT'
      Call Print_Array( Int(file_id,jpim), '(10000(' // Trim( Gformat ) // '))', A1 = radiance%bt )
      Write( file_id, * ) ')'
    End If

    If (Any(radiance%refl /= 0._jprb)) Then
      Write( file_id, '(A," = (")' ) name//'%REFL'
      Call Print_Array( Int(file_id,jpim), '(10000(' // Trim( Gformat ) // '))', A1 = radiance%refl )
      Write( file_id, * ) ')'
    End If

    Write( file_id, '(A," = (")' ) name//'%CLEAR'
    Call Print_Array( Int(file_id,jpim), '(10000(' // Trim( Gformat ) // '))', A1 = radiance%clear )
    Write( file_id, * ) ')'

    If (Any(radiance%bt_clear /= 0._jprb)) Then
      Write( file_id, '(A," = (")' ) name//'%BT_CLEAR'
      Call Print_Array( Int(file_id,jpim), '(10000(' // Trim( Gformat ) // '))', A1 = radiance%bt_clear )
      Write( file_id, * ) ')'
    End If

    If (Any(radiance%refl_clear /= 0._jprb)) Then
      Write( file_id, '(A," = (")' ) name//'%REFL_CLEAR'
      Call Print_Array( Int(file_id,jpim), '(10000(' // Trim( Gformat ) // '))', A1 = radiance%refl_clear )
      Write( file_id, * ) ')'
    EndIf

    Write( file_id, '(A," = (")' ) name//'%OVERCAST'
    Call Print_Array( Int(file_id,jpim), '(10000(' // Trim( Gformat ) // '))', A2 = radiance%overcast )
    Write( file_id, * ) ')'

    Write( file_id, '(A," = (")' ) name//'%CLOUDY100%'
    Call Print_Array( Int(file_id,jpim), '(10000(' // Trim( Gformat ) // '))', A1 = radiance%cloudy )
    Write( file_id, * ) ')'

    If (PRESENT(radiance2)) Then

      Write( file_id, '(A," = (")' ) name//'%UPCLEAR'
      Call Print_Array( Int(file_id,jpim), '(10000(' // Trim( Gformat ) // '))', A1 = radiance2%upclear )
      Write( file_id, * ) ')'

      Write( file_id, '(A," = (")' ) name//'%DNCLEAR'
      Call Print_Array( Int(file_id,jpim), '(10000(' // Trim( Gformat ) // '))', A1 = radiance2%dnclear )
      Write( file_id, * ) ')'

      Write( file_id, '(A," = (")' ) name//'%REFLDNCLEAR'
      Call Print_Array( Int(file_id,jpim), '(10000(' // Trim( Gformat ) // '))', A1 = radiance2%refldnclear )
      Write( file_id, * ) ')'

      Write( file_id, '(A," = (")' ) name//'%UP'
      Call Print_Array( Int(file_id,jpim), '(10000(' // Trim( Gformat ) // '))', A2 = radiance2%up )
      Write( file_id, * ) ')'

      Write( file_id, '(A," = (")' ) name//'%DOWN'
      Call Print_Array( Int(file_id,jpim), '(10000(' // Trim( Gformat ) // '))', A2 = radiance2%down )
      Write( file_id, * ) ')'

      Write( file_id, '(A," = (")' ) name//'%SURF'
      Call Print_Array( Int(file_id,jpim), '(10000(' // Trim( Gformat ) // '))', A2 = radiance2%surf )
      Write( file_id, * ) ')'

    EndIf

    Close( file_id, iostat = err )
    THROWM(err.ne.0,"Cannot close "//Trim(f))
    call rttov_put_lun (file_id)

CATCH
  
  End Subroutine
  
  Subroutine PrintTransmission (err, transmission, name, f )
    Use rttov_chain, only : print_array
!
    Implicit None
    
    Integer(Kind=jpim),      Intent(out) :: err
    Type(transmission_type), Intent(in) :: transmission
    Character(len=*),        Intent(in) :: name, f
!

    Integer(Kind=jpim) :: file_id

TRY

    call rttov_get_lun (file_id )
    Open( file_id, file = f, form = 'formatted', iostat = err )
    THROWM(err.ne.0,"Cannot open "//Trim(f))

    Write( file_id, '(A," = (")' ) name//'%TAU_TOTAL'
    Call Print_Array( Int(file_id,jpim), '(10000(' // Trim( Gformat ) // '))', A1 = transmission%tau_total )
    Write( file_id, * ) ')'

    If (Any(transmission%tausun_total_path1 /= 0._jprb)) Then
      Write( file_id, '(A," = (")' ) name//'%TAUSUN_TOTAL_PATH1'
      Call Print_Array( Int(file_id,jpim), '(10000(' // Trim( Gformat ) // '))', A1 = transmission%tausun_total_path1 )
      Write( file_id, * ) ')'
    EndIf

    If (Any(transmission%tausun_total_path2 /= 0._jprb)) Then
      Write( file_id, '(A," = (")' ) name//'%TAUSUN_TOTAL_PATH2'
      Call Print_Array( Int(file_id,jpim), '(10000(' // Trim( Gformat ) // '))', A1 = transmission%tausun_total_path2 )
      Write( file_id, * ) ')'
    EndIf

    Write( file_id, '(A," = (")' ) name//'%TAU_LAYERS'
    Call Print_Array( Int(file_id,jpim), '(10000(' // Trim( Gformat ) // '))', A2 = transmission%tau_levels )
    Write( file_id, * ) ')'

    If (Any(transmission%tausun_levels_path1 /= 0._jprb)) Then
      Write( file_id, '(A," = (")' ) name//'%TAUSUN_LAYERS_PATH1'
      Call Print_Array( Int(file_id,jpim), '(10000(' // Trim( Gformat ) // '))', A2 = transmission%tausun_levels_path1 )
      Write( file_id, * ) ')'
    EndIf

    If (Any(transmission%tausun_levels_path2 /= 0._jprb)) Then
      Write( file_id, '(A," = (")' ) name//'%TAUSUN_LAYERS_PATH2'
      Call Print_Array( Int(file_id,jpim), '(10000(' // Trim( Gformat ) // '))', A2 = transmission%tausun_levels_path2 )
      Write( file_id, * ) ')'
    EndIf

    Close( file_id, iostat = err )
    THROWM(err.ne.0,"Cannot close "//Trim(f))
    call rttov_put_lun (file_id)

CATCH
  
  End Subroutine
  
  
  Subroutine PrintEmissivity (err, emissivity, name, f )
    Use rttov_chain, only : print_array
!
    Implicit None
    
    Integer(Kind=jpim), Intent(out) :: err
    Real(Kind=jprb),    Intent(in)  :: emissivity(:)
    Character(len=*),   Intent(in)  :: name, f
!

    Integer(Kind=jpim) :: file_id

TRY


    call rttov_get_lun (file_id )
    Open( file_id, file = f, form = 'formatted', iostat = err )
    THROWM(err.ne.0,"Cannot open "//Trim(f))
    
    Write( file_id, '(A," = (")' ) name
    Call Print_Array( Int(file_id,jpim), '(10000(' // Trim( Gformat ) // '))', A1 = emissivity )
    Write( file_id, * ) ')'
    
    Close( file_id, iostat = err )
    THROWM(err.ne.0,"Cannot close "//Trim(f))
    call rttov_put_lun (file_id)

CATCH
  
  End Subroutine

  Subroutine PrintReflectance (err, reflectance, name, f )
    Use rttov_chain, only : print_array
!
    Implicit None
    
    Integer(Kind=jpim), Intent(out) :: err
    Real(Kind=jprb),    Intent(in)  :: reflectance(:)
    Character(len=*),   Intent(in)  :: name, f
!

    Integer(Kind=jpim) :: file_id

TRY


    call rttov_get_lun (file_id )
    Open( file_id, file = f, form = 'formatted', iostat = err )
    THROWM(err.ne.0,"Cannot open "//Trim(f))
    
    Write( file_id, '(A," = (")' ) name
    Call Print_Array( Int(file_id,jpim), '(10000(' // Trim( Gformat ) // '))', A1 = reflectance )
    Write( file_id, * ) ')'
    
    Close( file_id, iostat = err )
    THROWM(err.ne.0,"Cannot close "//Trim(f))
    call rttov_put_lun (file_id)

CATCH
  
  End Subroutine
  
  Subroutine PrintPCScores(err, pccomp, name, f)
    Use rttov_chain, only : print_array
!
    Implicit None
    
    Integer(Kind=jpim),      Intent(out) :: err
    Type(rttov_pccomp),      Intent(in) :: pccomp
    Character(len=*),        Intent(in) :: name, f
!

    Integer(Kind=jpim) :: file_id

TRY

    call rttov_get_lun (file_id )
    Open( file_id, file = f, form = 'formatted', iostat = err )
    THROWM(err.ne.0,"Cannot open "//Trim(f))
    
    If(associated(pccomp%pcscores)) then
      Write( file_id, '(A," = (")' ) name//'%PCSCORES'
      Call Print_Array( Int(file_id,jpim), '(10000(' // Trim( Gformat ) // '))', A1 = pccomp%pcscores )
      Write( file_id, * ) ')'
    endif
    
    If(associated(pccomp%total_pccomp)) then
      Write( file_id, '(A," = (")' ) name//'%TOTAL_PCCOMP'
      Call Print_Array( Int(file_id,jpim), '(10000(' // Trim( Gformat ) // '))', A1 = pccomp%total_pccomp )
      Write( file_id, * ) ')'
    endif
    
    If(associated(pccomp%bt_pccomp)) then
      Write( file_id, '(A," = (")' ) name//'%BT_PCCOMP'
      Call Print_Array( Int(file_id,jpim), '(10000(' // Trim( Gformat ) // '))', A1 = pccomp%bt_pccomp )
      Write( file_id, * ) ')'
    endif
    
    Close( file_id, iostat = err )
    THROWM(err.ne.0,"Cannot close "//Trim(f))
    call rttov_put_lun (file_id)

CATCH
  

  End Subroutine


  SUBROUTINE ppmvdry2ppmvwet(prof)
    TYPE(profile_type), INTENT(INOUT) :: prof

    ! Convert q first
    prof%q =  prof%q / (1._jprb + 1.E-06_jprb * prof%q)
    prof%s2m%q =  prof%s2m%q / (1._jprb + 1.E-06_jprb * prof%s2m%q)

    ! q is now ppmv wet
    IF (ASSOCIATED(prof%o3)) THEN
      prof%o3 = prof%o3 * (1._jprb - 1.E-06_jprb * prof%q)
      prof%s2m%o = prof%s2m%o * (1._jprb - 1.E-06_jprb * prof%s2m%q)
    ENDIF
    IF (ASSOCIATED(prof%co2)) THEN
      prof%co2 = prof%co2 * (1._jprb - 1.E-06_jprb * prof%q)
    ENDIF
    IF (ASSOCIATED(prof%co)) THEN
      prof%co = prof%co * (1._jprb - 1.E-06_jprb * prof%q)
    ENDIF
    IF (ASSOCIATED(prof%ch4)) THEN
      prof%ch4 = prof%ch4 * (1._jprb - 1.E-06_jprb * prof%q)
    ENDIF
    IF (ASSOCIATED(prof%n2o)) THEN
      prof%n2o = prof%n2o * (1._jprb - 1.E-06_jprb * prof%q)
    ENDIF
  END SUBROUTINE

  SUBROUTINE ppmvdry2kgkgwet(prof)
    TYPE(profile_type), INTENT(INOUT) :: prof

    prof%q = prof%q / (Mair * 1.E06_jprb / gas_mass(gas_id_watervapour) + prof%q)
    prof%s2m%q = prof%s2m%q / (Mair * 1.E06_jprb / gas_mass(gas_id_watervapour) + prof%s2m%q)

    IF (ASSOCIATED(prof%o3)) THEN
      prof%o3 = prof%o3 / (Mair * 1.E06_jprb / gas_mass(gas_id_ozone) + prof%o3)
      prof%s2m%o = prof%s2m%o / (Mair * 1.E06_jprb / gas_mass(gas_id_ozone) + prof%s2m%o)
    ENDIF
    IF (ASSOCIATED(prof%co2)) THEN
      prof%co2 = prof%co2 / (Mair * 1.E06_jprb / gas_mass(gas_id_co2) + prof%co2)
    ENDIF
    IF (ASSOCIATED(prof%co)) THEN
      prof%co = prof%co / (Mair * 1.E06_jprb / gas_mass(gas_id_co) + prof%co)
    ENDIF
    IF (ASSOCIATED(prof%ch4)) THEN
      prof%ch4 = prof%ch4 / (Mair * 1.E06_jprb / gas_mass(gas_id_ch4) + prof%ch4)
    ENDIF
    IF (ASSOCIATED(prof%n2o)) THEN
      prof%n2o = prof%n2o / (Mair * 1.E06_jprb / gas_mass(gas_id_n2o) + prof%n2o)
    ENDIF

  END SUBROUTINE

  SUBROUTINE ppmvwet2ppmvdry(prof)
    TYPE(profile_type), INTENT(INOUT) :: prof

    ! q is ppmv wet
    IF (ASSOCIATED(prof%o3)) THEN
      prof%o3 = prof%o3 / (1._jprb - 1.E-06_jprb * prof%q)
      prof%s2m%o = prof%s2m%o / (1._jprb - 1.E-06_jprb * prof%s2m%q)
    ENDIF
    IF (ASSOCIATED(prof%co2)) THEN
      prof%co2 = prof%co2 / (1._jprb - 1.E-06_jprb * prof%q)
    ENDIF
    IF (ASSOCIATED(prof%co)) THEN
      prof%co = prof%co / (1._jprb - 1.E-06_jprb * prof%q)
    ENDIF
    IF (ASSOCIATED(prof%ch4)) THEN
      prof%ch4 = prof%ch4 / (1._jprb - 1.E-06_jprb * prof%q)
    ENDIF
    IF (ASSOCIATED(prof%n2o)) THEN
      prof%n2o = prof%n2o / (1._jprb - 1.E-06_jprb * prof%q)
    ENDIF

    ! Convert q last
    prof%q =  prof%q / (1._jprb - 1.E-06_jprb * prof%q)
    prof%s2m%q =  prof%s2m%q / (1._jprb - 1.E-06_jprb * prof%s2m%q)

  END SUBROUTINE

  SUBROUTINE ppmvwet2kgkgwet(prof)
    TYPE(profile_type), INTENT(INOUT) :: prof

    CALL ppmvwet2ppmvdry(prof)
    CALL ppmvdry2kgkgwet(prof)

  END SUBROUTINE

  SUBROUTINE kgkgwet2ppmvdry(prof)
    TYPE(profile_type), INTENT(INOUT) :: prof

    prof%q = (1.E06_jprb * Mair / gas_mass(gas_id_watervapour)) * prof%q / (1._jprb - prof%q)
    prof%s2m%q = (1.E06_jprb * Mair / gas_mass(gas_id_watervapour)) * prof%s2m%q / (1._jprb - prof%s2m%q)

    IF (ASSOCIATED(prof%o3)) THEN
      prof%o3 = (1.E06_jprb * Mair / gas_mass(gas_id_ozone)) * prof%o3 / (1._jprb - prof%o3)
      prof%s2m%o = (1.E06_jprb * Mair / gas_mass(gas_id_ozone)) * prof%s2m%o / (1._jprb - prof%s2m%o)
    ENDIF
    IF (ASSOCIATED(prof%co2)) THEN
      prof%co2 = (1.E06_jprb * Mair / gas_mass(gas_id_co2)) * prof%co2 / (1._jprb - prof%co2)
    ENDIF
    IF (ASSOCIATED(prof%co)) THEN
      prof%co = (1.E06_jprb * Mair / gas_mass(gas_id_co)) * prof%co / (1._jprb - prof%co)
    ENDIF
    IF (ASSOCIATED(prof%ch4)) THEN
      prof%ch4 = (1.E06_jprb * Mair / gas_mass(gas_id_ch4)) * prof%ch4 / (1._jprb - prof%ch4)
    ENDIF
    IF (ASSOCIATED(prof%n2o)) THEN
      prof%n2o = (1.E06_jprb * Mair / gas_mass(gas_id_n2o)) * prof%n2o / (1._jprb - prof%n2o)
    ENDIF

  END SUBROUTINE

  SUBROUTINE kgkgwet2ppmvwet(prof)
    TYPE(profile_type), INTENT(INOUT) :: prof

    CALL kgkgwet2ppmvdry(prof)
    CALL ppmvdry2ppmvwet(prof)

  END SUBROUTINE

End Program





