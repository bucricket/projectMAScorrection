Module rttov_unix_env
! Description:
!   Wraps all non portable features and other useful functions
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

#ifdef RTTOV_USE_F90_UNIX_ENV
use f90_unix_env
use f90_unix_proc
use f90_unix_errno
#endif

Use parkind1, only : jpim, jprb, jplm

Contains

Subroutine rttov_getenv( key, val )
  Character(Len=*), Intent(in) :: key
  Character(Len=*), Intent(out) :: val
#ifdef RTTOV_USE_F90_UNIX_ENV
#ifdef RTTOV_NAG51
  Integer :: errno
#else
  Integer(error_kind) :: errno
#endif
  Call getenv( key, val, errno = errno )
  If (errno .ne. 0) Then
    val = ""
  EndIf
#else
  Call getenv( key, val )
#endif
End Subroutine


Integer(Kind=jpim) Function rttov_iargc()
  rttov_iargc = iargc()
End Function

Subroutine rttov_getarg( key, val )
  Integer(Kind=jpim), Intent(in) :: key
  Character(Len=*), Intent(out) :: val
  Call getarg( Int(key,selected_int_kind(9)), val )
End Subroutine

Subroutine rttov_exit( status )
  Integer(Kind=jpim), Intent(in) :: status
  Call exit( Int(status,selected_int_kind(9)) )
End Subroutine

Subroutine rttov_mkdir( path )
  Character(len=*), Intent(in) :: path
  Call system( "mkdir -p " // Trim(path) )
End Subroutine

Character(len=256) Function rttov_dirname( path )
  Character(len=*), Intent(in) :: path
  
  Integer(Kind=jpim) :: i
  rttov_dirname = ""
  i = Len( Trim( path ) ) - 1
  Do
    If( i .le. 0 ) Return
    If( path(i:i) .eq. '/' ) Exit
    i = i - 1
  EndDo
  rttov_dirname = path(1:i)
End Function

Character(len=256) Function rttov_basename( path )
  Character(len=*), Intent(in) :: path
  
  Integer(Kind=jpim) :: i
  rttov_basename = ""
  i = Len( Trim( path ) ) - 1
  Do
    If( i .le. 0 ) Then
      i = 0
      Exit
    EndIf
    If( path(i:i) .eq. '/' ) Exit
    i = i - 1
  EndDo
  rttov_basename = path(i+1:)
End Function

elemental subroutine rttov_lower_case(ous,ins)
! convert a word to lower case
character (len=*) , intent(out) :: ous
character (len=*) , intent(in) :: ins
integer :: i,ic,nlen
nlen = len(ins)
ous = ''
do i=1,nlen
   ic = ichar(ins(i:i))
   if (ic >= 65 .and. ic < 90) then
     ous(i:i) = char(ic+32)
   else
     ous(i:i) = ins(i:i)
   endif
end do
end subroutine rttov_lower_case 

elemental subroutine rttov_upper_case(ous,ins)
! convert a word to upper case
character (len=*) , intent(out) :: ous
character (len=*) , intent(in) :: ins
integer :: i,ic,nlen
nlen = len(ins)
ous = ''
do i=1,nlen
   ic = ichar(ins(i:i))
   if (ic >= 97 .and. ic < 122) then
     ous(i:i) = char(ic-32)
   else
     ous(i:i) = ins(i:i)
   endif
end do
end subroutine rttov_upper_case 

logical(Kind=jplm) function rttov_isalpha(c)
character, intent(in) :: c

rttov_isalpha = ((c.ge.'a').and.(c.le.'z'))&
            .or.((c.ge.'A').and.(c.le.'Z'))

end function

logical(Kind=jplm) function rttov_isdigit(c)
character, intent(in) :: c

rttov_isdigit = (c.ge.'0').and.(c.le.'9')

end function

subroutine rttov_date_and_time( vl )
integer(kind=jpim), intent(out) :: vl(8)
!
integer :: vlx(8)

  call date_and_time( values = vlx )
  
  vl = vlx
end subroutine

subroutine rttov_cpu_time( t )
  real, intent(out) :: t
  call cpu_time( t )
end subroutine

subroutine rttov_countlines( nlines, f, err )
integer(kind=jpim), intent(out) :: nlines
character(len=*),   intent(in) :: f
integer(kind=jpim), intent(out) :: err
character(len=32) :: str

nlines = 0
open( 77, file = f, err = 888 )

do 
  read( 77, *, err = 888, end = 777 ) str
  nlines = nlines + 1
enddo

777 continue

close( 77 )

return
888 continue
  err = 1_jpim
end subroutine

Integer(Kind=jpim) Function Rttov_CountWords( s )
  Implicit None
  Character(len=*), Intent(in) :: s
  Integer(Kind=jpim) :: n, i, l
  Logical(Kind=jplm) :: in
  n = 0_jpim
  in = .false.
  l = Len( Trim( s ) )
  Do i = 1, l
    If( s(i:i) .eq. ' ' ) Then
      in = .false.
    Else If( .not. in ) Then
      n = n + 1
      in = .true.
    EndIf
  EndDo
  Rttov_CountWords = n
End Function


subroutine rttov_banner( lun, text, err )

  integer(kind=jpim), intent(in) :: lun
  character(len=*),   intent(in) :: text(:)
  integer(kind=jpim), intent(out) :: err

  integer(kind=jpim), parameter :: xcol = 8_jpim
  character(len=*),   parameter :: cstar = '**********'

  integer(kind=jpim) :: ln, il, nl, ic, nc1, nc2, kc

  err = 0

  ln = len( text(1) )
  nl = size( text )

  do ic = 1, xcol
    write( lun, '(a)', advance = 'no', err = 888 ) cstar
  enddo
  write( lun, *, err = 888 )

  kc = xcol * 10 - 2 - ln

  if( kc .lt. 0 ) then
    ln = xcol * 10 - 2
    kc = 0
  endif

  nc1 = kc / 2
  nc2 = kc - nc1

  do il = 1, nl + 2
    write( lun, '(a)', advance = 'no', err = 888 ) '*'
    if( ( il .eq. 1 ) .or. ( il .eq. nl + 2 ) ) then
      do ic = 1, xcol * 10 - 2
        write( lun, '(" ")', advance = 'no', err = 888 )
      enddo
    else
      do ic = 1, nc1
        write( lun, '(" ")', advance = 'no', err = 888 )
      enddo
      write( lun, '(a)', advance = 'no', err = 888 ) text(il-1)(1:ln)
      do ic = 1, nc2
        write( lun, '(" ")', advance = 'no', err = 888 )
      enddo
    endif
    write( lun, '(a)', advance = 'no', err = 888 ) '*'
    write( lun, *, err = 888 )
  enddo

  do ic = 1, xcol
    write( lun, '(a)', advance = 'no', err = 888 ) cstar
  enddo
  write( lun, *, err = 888 )


  return
888  continue
  err = 1_jpim
  
end subroutine

subroutine rttov_unlink( f, err )
character(len=*), intent(in) :: f
integer(Kind=jpim), intent(out) :: err
err = 0
open(77, file = f, iostat = err )
if( err .ne. 0 ) return
close(77, status = 'delete', iostat = err)
end subroutine
End Module
