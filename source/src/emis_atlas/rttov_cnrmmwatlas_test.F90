!
program rttov_cnrmmwatlas_test
  ! Description:
  !   Test program for CNRM MW atlas.
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
  !  1.0      27/05/2010  Created (J. Hocking P. Brunel)
  !
  ! Code Description:
  !   Language:           Fortran 90.
  !   Software Standards: "European Standards for Writing and
  !     Documenting Exchangeable Fortran 90 Code".
  !
  ! Declarations:
  ! Modules used:
  !
#include "throw.h"
  use parkind1, only : jpim, jprb

  use rttov_unix_env, only : rttov_exit

  use rttov_const, only :    & 
        errorstatus_success, &
        surftype_land

  use rttov_types, only : &
        rttov_coefs,      &
        rttov_options,    &
        profile_type,     &
        rttov_chanprof
        
  use mod_cnrm_mw_atlas, only : &
        cnrm_mw_atlas_version => mw_atlas_version 
  
  implicit none
  
  integer(kind=jpim), parameter :: ioout = 51 ! Output file unit

  character(len=256) :: coef_filename ! Coefficient filename
  character(len=256) :: atlas_path    ! path to atlas data
  integer(kind=jpim) :: imonth        ! Month for which to load data
  
  integer(kind=jpim) :: lo,hi,j,k
  integer(kind=jpim) :: err
  
  integer(kind=jpim) :: nprof  ! Number of profiles
  integer(kind=jpim) :: nchan  ! Number of channels per profile
  real(kind=jprb)    :: latitude_start,  latitude_step
  real(kind=jprb)    :: longitude_start, longitude_step
  real(kind=jprb)    :: zenangle_start,  zenangle_step
  integer(kind=jpim), allocatable :: channels(:)
  
  type(rttov_coefs)                 :: coefs
  type(rttov_options)               :: opts
  type(profile_type),   allocatable :: profiles(:)
  type(rttov_chanprof), allocatable :: chanprof(:)
  
  real(kind=jprb),    allocatable :: emissivity(:)

#include "rttov_read_coefs.interface"
#include "rttov_dealloc_coefs.interface"

#include "rttov_setup_emis_atlas.interface"
#include "rttov_get_emis.interface"
#include "rttov_deallocate_emis_atlas.interface"

!-----------------------------------------------------------------------------------
TRY

! Read input data
  read(*,'(a)') coef_filename
  read(*,'(a)') atlas_path
  read(*,*) imonth
  read(*,*) nprof
  read(*,*) latitude_start,  latitude_step
  read(*,*) longitude_start, longitude_step
  read(*,*) zenangle_start,  zenangle_step
  read(*,*) nchan

  allocate(channels(nchan))
  read(*,*) channels(1:nchan)

  allocate(profiles(nprof))
  do k = 1, nprof
    profiles(k)%latitude      = latitude_start  + (k-1)* latitude_step
    profiles(k)%longitude     = longitude_start + (k-1)* longitude_step
    profiles(k)%zenangle      = zenangle_start  + (k-1)* zenangle_step
    ! do not forget to initialise surface type as the rttov_get_emis is testing it
    profiles(k)%skin%surftype = surftype_land
  end do
 
    
! Generate the chanprof array
  allocate(chanprof(nchan*nprof))
  do k = 1, nprof
    lo = (k-1)*nchan+1
    hi = lo+nchan-1
    chanprof(lo:hi)%prof = k
    chanprof(lo:hi)%chan = (/ (j, j=1,nchan) /)
  end do
  
  allocate(emissivity(nchan*nprof))
    
! Set up RTTOV coefficients
  call rttov_read_coefs(              &
              err,                    &
              coefs,                  &
              opts,                   &
              channels=channels(:),   &
              file_coef = coef_filename)

!-------------------------
! Set up emissivity atlas
!-------------------------
  call rttov_setup_emis_atlas(   &
              err,               &
              opts,              &
              imonth,            &
              coefs,             &
              path = atlas_path, &
              mw_atlas_ver = cnrm_mw_atlas_version)
  THROWM( ERR .NE. 0, "Failure in setting up emissivity atlas" )

!----------------------------
! Retrieve values from atlas
!----------------------------
  call rttov_get_emis(                &
              err,                    &
              opts,                   &
              chanprof,               &
              profiles,               &
              coefs,                  &
              emissivity = emissivity )
  THROWM( ERR .NE. 0, "Failure in retrieving emissivity values" )
  
! Write out emissivity data
  open(ioout, file='output_cnrmmwatlas.ascii', &
     & form='formatted', status='replace', action='write')
  
  do k = 1, nprof
    write(ioout,'(a8,i4,3(1x,a10,f10.3))') 'Profile ',k, &
      & 'longitude',profiles(k)%longitude, &
      & 'latitude ',profiles(k)%latitude,  &
      & 'zenangle ',profiles(k)%zenangle
    write(ioout,'(a)') ' Chan     Emissivity'
    do j = 1, nchan
      lo = (k-1)*nchan
      write(ioout,'(i5,5x,f8.5)') channels(j), emissivity(lo+j)
    end do
  end do
  close(ioout)
  
!-----------------------
! Deallocate atlas data
!-----------------------
  call rttov_deallocate_emis_atlas(coefs)
    
! Tidy up the remaining data
  call rttov_dealloc_coefs(err, coefs)

  deallocate(emissivity)
  deallocate(chanprof)
  deallocate(channels)
  deallocate(profiles)

! End of test
PCATCH
end program rttov_cnrmmwatlas_test


