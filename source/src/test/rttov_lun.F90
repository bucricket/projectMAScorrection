module rttov_lun
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

!$ use omp_lib

!
! Fortran logical units management with OpenMP
!

use parkind1,       only : jpim, jplm
use rttov_unix_env, only : rttov_exit

implicit none

private

logical(Kind=jplm) :: lunf(10:25)
logical(Kind=jplm) :: init = .false.

integer(Kind=jpim), parameter :: lun_min = lbound(lunf,1), &
                                 lun_max = ubound(lunf,1)

public :: rttov_get_lun, rttov_put_lun

contains

subroutine rttov_dump_lun()
  integer(Kind=jpim) :: lun
!$write(0,'(I3," ")',advance='no') omp_get_thread_num ()
  do lun = lun_min, lun_max
    write(0,'(L1)',advance='no') lunf(lun)
  enddo
  write(0,*)
end subroutine

subroutine rttov_get_lun(lun_out)
integer(Kind=jpim), intent(out) :: lun_out
integer(Kind=jpim) :: lun


  lun_out = -1

!$OMP CRITICAL (rttov_lun_lock)
!$OMP FLUSH(lunf)
  if(.not.init) then
    lunf = .true.
    init = .true.
  endif
  do lun = lun_min, lun_max
    if(lun_out .le. 0 .and. lunf(lun)) then
      lunf(lun) = .false.
      lun_out = lun
    endif
  enddo
!$OMP FLUSH(lunf)
!$OMP END CRITICAL (rttov_lun_lock)
! write(0,*) " get = ", lun_out, omp_get_thread_num ()
end subroutine

subroutine rttov_put_lun(lun)
integer(Kind=jpim), intent(in) :: lun
!$OMP CRITICAL (rttov_lun_lock)
!$OMP FLUSH(lunf)
  if (lunf(lun)) then
    write(0,*) lun, "was already free"
    call rttov_exit(1_jpim)
  endif
  lunf(lun) = .true.
!$OMP FLUSH(lunf)
!$OMP END CRITICAL (rttov_lun_lock)
! write(0,*) " put = ", lun, omp_get_thread_num ()
end subroutine

end module
