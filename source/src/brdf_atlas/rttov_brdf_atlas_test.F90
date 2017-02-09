!
program rttov_brdf_atlas_test
  ! Description:
  !   Test program for VIS/NIR BRDF atlas.
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
  !    Copyright 2012, EUMETSAT, All Rights Reserved.
  !
  ! Method:
  !
  ! Current Code Owner: SAF NWP
  !
  ! History:
  ! Version   Date     Comment
  ! -------   ----     -------
  !  1.0      25/12/2012  Created (J. Vidot)
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
  use parkind1, only : jpim, jprb, jplm

  use rttov_unix_env, only : rttov_exit
  
  use rttov_const, only : errorstatus_success
  
  use rttov_types, only : &
        rttov_coefs,      &
        rttov_options,    &
        profile_type,     &
        rttov_chanprof
  
  implicit none

  real(kind=jprb),    parameter :: difftol = 1.E-12_jprb

  integer(kind=jpim), parameter :: ioin  = 50 ! Input file unit
  integer(kind=jpim), parameter :: ioout = 51 ! Output file unit

  character(len=256) :: coef_filename ! Coefficient filename
  character(len=256) :: atlas_path    ! Path to atlas data
  integer(kind=jpim) :: imonth        ! Month for which to load data
  
  integer(kind=jpim) :: lo,hi,i,j,k
  integer(kind=jpim) :: err
  
  integer(kind=jpim) :: nprof  ! Number of profiles
  integer(kind=jpim) :: nchan  ! Number of channels per profile

  logical(kind=jplm) :: single_instrument
  logical(kind=jplm) :: diffflag

  integer(kind=jpim), allocatable :: channels(:)
  
  type(rttov_coefs)                 :: coefs
  type(rttov_options)               :: opts
  type(profile_type),   allocatable :: profiles(:)
  type(rttov_chanprof), allocatable :: chanprof(:)
  
  real(kind=jprb),    allocatable :: reflectance(:,:)
  real(kind=jprb),    allocatable :: bh_albedo(:,:)
  integer(kind=jpim), allocatable :: brdf_flag(:,:)


#include "rttov_read_coefs.interface"
#include "rttov_dealloc_coefs.interface"

#include "rttov_setup_brdf_atlas.interface"
#include "rttov_get_brdf.interface"
#include "rttov_deallocate_brdf_atlas.interface"

!-----------------------------------------------------------------------------------
TRY

! This test uses a simple input file for 'profiles' consisting of lat, lon, surface
! type, water type, viewing zenith angle, viewing azimuth angle, solar zenith angle
! and solar azimuth angle. The controlling shell script passes the month, instrument,
! and channels for which to retrieve brdf.

! BRDF values and the associated atlas flags are returned for all 
! profiles and are written to an output file for comparison with reference output.

! The test is repeated for two cases, where the atlas is initialised for use with
! multiple instruments, and for use with a single instrument.

! Read profiles input file
  open(ioin, file='profiles_visnir', status='old')
  call skipcomment(ioin)
  
  read(ioin,*) nprof
  call skipcomment(ioin)
  
  allocate(profiles(nprof))
  do k = 1, nprof
    read(ioin,*) profiles(k)%latitude,      &
                 profiles(k)%longitude,     &
                 profiles(k)%skin%surftype, &
                 profiles(k)%skin%watertype, &
                 profiles(k)%zenangle,      &
                 profiles(k)%azangle,       &
                 profiles(k)%sunzenangle,   &
                 profiles(k)%sunazangle
    call skipcomment(ioin)
  end do
  
  close(ioin)
  
! Read data from control script
  read(*,*) imonth
  read(*,'(a)') coef_filename
  read(*,*) nchan
  allocate(channels(nchan))
  read(*,*) channels(1:nchan)
  read(*,'(a)') atlas_path
  
! Generate the chanprof array
  allocate(chanprof(nchan*nprof))
  do k = 1, nprof
    lo = (k-1)*nchan+1
    hi = lo+nchan-1
    chanprof(lo:hi)%prof = k
    chanprof(lo:hi)%chan = (/ (j, j=1,nchan) /)
  end do
  
  allocate(reflectance(nchan*nprof,2))
  allocate(bh_albedo(nchan*nprof,2))
  allocate(brdf_flag(nchan*nprof,2))
    
! Set up RTTOV coefficients
  call rttov_read_coefs(                &
              err,                      &
              coefs,                    &
              opts,                     &
              channels = channels(:),   &
              file_coef = coef_filename)

! Open file for output
  open(ioout, file='output_brdf_atlas.ascii', form='formatted', status='replace', action='write')

! Run twice to test both ways of initialising the atlas
  do i = 1, 2

    if (i == 1) then
      print *, 'Initialising atlas for multiple instruments'
      single_instrument = .FALSE.   ! initialise for use with multiple instruments
    else
      print *, 'Initialising atlas for a single instrument'
      single_instrument = .TRUE.    ! initialise for use with a single instrument
    endif

!-------------------------
! Set up emissivity atlas
!-------------------------
    call rttov_setup_brdf_atlas(      &
                err,                  &
                opts,                 &
                imonth,               &
                coefs,                &
                path = atlas_path,    &
                brdf_atlas_single_instrument = single_instrument)
    THROWM( ERR .NE. 0, "Failure in setting up BRDF atlas" )

!----------------------------
! Retrieve values from atlas
!----------------------------
    call rttov_get_brdf(                    &
                err,                        &
                chanprof,                   &
                profiles,                   &
                coefs,                      &
                brdf = reflectance(:,i),    &
                brdf_flag = brdf_flag(:,i), &
                bh_albedo = bh_albedo(:,i))
    THROWM( ERR .NE. 0, "Failure in retrieving BRDF values" )

! Write out emissivity data
    write(ioout,'(a,l1)') 'Init for single inst? ', single_instrument
    do k = 1, nprof
      write(ioout,'(a,i2.2)') 'Profile ',k
      write(ioout,'(a)') ' Chan     brdf     BH albedo    Flag'
      do j = 1, nchan
        lo = (k-1)*nchan
        write(ioout,'(i5,2f11.4,i6)') channels(j), reflectance(lo+j,i), bh_albedo(lo+j,i), brdf_flag(lo+j,i)
      end do
      write(ioout,'(a)') ''
    end do

!-----------------------
! Deallocate atlas data
!-----------------------
    call rttov_deallocate_brdf_atlas(coefs)

  end do

! Ensure single and multi instrument runs gave same results
  diffflag = .false.
  if (any(abs(reflectance(:,1) - reflectance(:,2)) > difftol)) then
    write(ioout,'(a)') 'WARNING: difference in BRDF between single and multi-inst runs'
    diffflag = .true.
  endif
  if (any(abs(brdf_flag(:,1) - brdf_flag(:,2)) > difftol)) then
    write(ioout,'(a)') 'WARNING: difference in BRDF flag between single and multi-instrument results'
    diffflag = .true.
  endif
  if (any(abs(bh_albedo(:,1) - bh_albedo(:,2)) > difftol)) then
    write(ioout,'(a)') 'WARNING: difference in bh_albedo between single and multi-instrument results'
    diffflag = .true.
  endif
  if (diffflag) then
    write(*,*) 'WARNING: difference in results between single- and multi-instrument results'
  endif

! Close the output file
  close(ioout)

! Tidy up the remaining data
  call rttov_dealloc_coefs(err, coefs)

  deallocate(reflectance)
  deallocate(bh_albedo)
  deallocate(brdf_flag)
  deallocate(chanprof)
  deallocate(channels)
  deallocate(profiles)

! End of test
PCATCH
    
contains
  subroutine skipcomment(fileunit)
    use parkind1, only : jpim
  
    integer(kind=jpim), intent(in) :: fileunit
  
    integer(kind=jpim) :: readstatus
    character(len=80)  :: line  ! input line
  
    do
      read( unit=fileunit, fmt='(a)', iostat=readstatus ) line
      if ( readstatus /= 0 ) exit
  
      line = adjustl(line)
      if ( line(1:1) == '!' .or. line(1:1) == '#' .or. line == '' ) then
          cycle !skip blank/comment lines
      else
          !reposition file at the start of the line and exit
          backspace( fileunit )
          exit
      end if
    end do
  end subroutine
  
end program rttov_brdf_atlas_test


