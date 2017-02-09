!
program rttov_mwatlas_test
  ! Description:
  !   Test program for MW atlas.
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
  !  1.0      02/06/2010  Created (J. Hocking)
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
  
  use rttov_const, only : errorstatus_success
  
  use rttov_types, only : &
        rttov_coefs,      &
        rttov_options,    &
        profile_type,     &
        rttov_chanprof
  
  implicit none
  
  integer(kind=jpim), parameter :: ioin  = 50 ! Input file unit
  integer(kind=jpim), parameter :: ioout = 51 ! Output file unit

  character(len=256) :: coef_filename ! Coefficient filename
  character(len=256) :: atlas_path    ! Path to atlas data
  integer(kind=jpim) :: imonth        ! Month for which to load data
  real(kind=jprb)    :: resol = 0.5   ! 'User-defined' resolution for testing
  
  integer(kind=jpim) :: lo,hi,j,k
  integer(kind=jpim) :: err
  character(len=12)  :: formatstr
  
  real(kind=jprb)    :: zenangle
  integer(kind=jpim) :: nprof  ! Number of profiles
  integer(kind=jpim) :: nchan  ! Number of channels per profile
  
  type(rttov_coefs)                 :: coefs
  type(rttov_options)               :: opts
  type(profile_type),   allocatable :: profiles(:)
  type(rttov_chanprof), allocatable :: chanprof(:)
  
  real(kind=jprb), allocatable :: emissivity(:)
  real(kind=jprb), allocatable :: emis_std(:)
  real(kind=jprb), allocatable :: emis_cov(:,:,:)

#include "rttov_read_coefs.interface"
#include "rttov_dealloc_coefs.interface"

#include "rttov_setup_emis_atlas.interface"
#include "rttov_get_emis.interface"
#include "rttov_deallocate_emis_atlas.interface"

!-----------------------------------------------------------------------------------
TRY

! This test uses a simple input file for 'profiles' consisting of lat, lon and 
! surface type. The controlling shell script passes the instrument and zenith angle 
! for which to retrieve emissivities.

! Emissivity values, errors and covariances for all channels are returned for all 
! profiles and are written to an output file for comparison with reference output.

! Read profiles input file
  open(ioin, file='profiles_mw', status='old')
  call skipcomment(ioin)
  
  read(ioin,*) nprof
  call skipcomment(ioin)
  
  allocate(profiles(nprof))
  do k = 1, nprof
    read(ioin,*) profiles(k)%latitude,      &
                 profiles(k)%longitude,     &
                 profiles(k)%skin%surftype
    call skipcomment(ioin)
  end do
  
  close(ioin)
  
! Read data from control script
  read(*,*) imonth
  read(*,'(a)') coef_filename
  read(*,*) zenangle
  read(*,'(a)') atlas_path
  
  profiles(:)%zenangle = zenangle
  
! Set up RTTOV coefficients
  call rttov_read_coefs(                &
              err,                      &
              coefs,                    &
              opts,                     &
              file_coef = coef_filename)

  nchan = coefs%coef%fmv_chn
  
! Generate the chanprof array
  allocate(chanprof(nchan*nprof))
  do k = 1, nprof
    lo = (k-1)*nchan+1
    hi = lo+nchan-1
    chanprof(lo:hi)%prof = k
    chanprof(lo:hi)%chan = (/ (j, j=1,nchan) /)
  end do
  
  allocate(emissivity(nchan*nprof))
  allocate(emis_std(nchan*nprof))
  allocate(emis_cov(nprof,nchan,nchan))
    
!-------------------------
! Set up emissivity atlas
!-------------------------
  call rttov_setup_emis_atlas( &
              err,             &
              opts,            &
              imonth,          &
              coefs,           &
              path = atlas_path)
  THROWM( ERR .NE. 0, "Failure in setting up emissivity atlas" )
  
!----------------------------
! Retrieve values from atlas
!----------------------------

! Open file for output
  open(ioout, file='output_mwatlas.ascii', form='formatted', status='replace', action='write')
  
! Nominal atlas resolution
  call rttov_get_emis(                 &
              err,                     &
              opts,                    &
              chanprof,                &
              profiles,                &
              coefs,                   &
              emissivity = emissivity, &
              emis_std = emis_std,     &
              emis_cov = emis_cov)
  THROWM( ERR .NE. 0, "Failure in retrieving emissivity values" )  

  write(formatstr,'(a,i2.2,a)') '(',nchan,'f11.7)'
  do k = 1, nprof
    write(ioout,'(a,i2.2)') 'Nominal atlas resol. Profile ',k
    write(ioout,'(a)') ' Chan  Emissivity  Standard Dev'
    do j = 1, nchan
      lo = (k-1)*nchan
      write(ioout,'(i5,f11.4,f13.6)') j, emissivity(lo+j), emis_std(lo+j)
    end do
    write(ioout,'(a)') 'Emissivity covariance matrix'  
    do j = 1, nchan
      write(ioout,formatstr) emis_cov(k,j,1:nchan)
    end do
  end do
  
! User-defined resolution
  call rttov_get_emis(                 &
              err,                     &
              opts,                    &
              chanprof,                &
              profiles,                &
              coefs,                   &
              resolution = resol,      &
              emissivity = emissivity, &
              emis_std = emis_std,     &
              emis_cov = emis_cov)
  THROWM( ERR .NE. 0, "Failure in retrieving emissivity values" )
    
  write(formatstr,'(a,i2.2,a)') '(',nchan,'f11.7)'
  do k = 1, nprof
    write(ioout,'(a,i2.2)') 'User-defined resol. Profile ',k
    write(ioout,'(a)') ' Chan  Emissivity  Standard Dev'
    do j = 1, nchan
      lo = (k-1)*nchan
      write(ioout,'(i5,f11.4,f13.6)') j, emissivity(lo+j), emis_std(lo+j)
    end do
    write(ioout,'(a)') 'Emissivity covariance matrix'  
    do j = 1, nchan
      write(ioout,formatstr) emis_cov(k,j,1:nchan)
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
  deallocate(emis_std)
  deallocate(emis_cov)
  deallocate(chanprof)
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

end program rttov_mwatlas_test
