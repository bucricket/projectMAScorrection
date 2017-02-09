Module Rttov_GetOptions
! Description:
!   Functions to handle program arguments
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
Use parkind1, only: jpim, jprb, jplm

Use rttov_unix_env, Only: rttov_iargc, rttov_getarg, &
  rttov_basename, rttov_countwords, rttov_getenv,    &
  rttov_isalpha, rttov_isdigit, rttov_exit

Implicit None

Interface GetOption
  Module Procedure GetOptionS, GetOptionSL, &
                   GetOptionI, GetOptionIL, &
                   GetOptionR, GetOptionRL, &
                   GetOptionB
                   
End Interface

!! @todo : list with fixed size

Public :: GetOption, InitOptions, InitOptionsStr, CheckOptions, AddGroup

Integer, Parameter :: argsizemax = 256

Character(len=argsizemax), Pointer :: myargs(:) => NULL()
Logical(Kind=jplm), Pointer :: check_args(:) => NULL()
Logical(Kind=jplm) :: lhelp  = .false., lshell = .false.

Character(Len=1056) :: message_opt = ""
Integer(kind=jpim) :: sepind(0:1024)

Type rttov_opt
  Character(Len=32) :: key, type
  Character(Len=1024) :: use
  Logical(Kind=jplm) :: group = .false.
End Type

Integer(Kind=jpim) :: nopt_seen
Type(rttov_opt), Pointer :: opt_seen(:) => NULL()

Private

Contains

Subroutine AddGroup( use )
Character(Len=*), Intent(in) :: use

Call init_opt_seen()
nopt_seen = nopt_seen + 1
Call grow_opt_seen()

opt_seen(nopt_seen)%group = .true.
opt_seen(nopt_seen)%use = use


End Subroutine

Character(len=argsizemax) Function get_env_opt( key )
Character(Len=*), Intent(in) :: key
Character(len=argsizemax) :: key_env, val_env
Integer(Kind=jpim) :: i, n
Character :: c

key_env = key(3:)

n = len(Trim(key_env))
do i = 1, n
  c = key_env(i:i)
  if((.not.rttov_isalpha(c)) .and. &
     (.not.rttov_isdigit(c)) .and. &
     (c .ne. '_' )) then
    key_env(i:i) = '_'
  endif
enddo

val_env = ""
Call rttov_getenv( 'rttov_opt_'//Trim(key_env), val_env )

!Print *, " key = ", Trim(key_env), " val = ", Trim(val_env)

get_env_opt = val_env

End Function

Subroutine mygetarg( i, s )
  Integer(Kind=jpim), Intent(in)  :: i
  Character(Len=*),      Intent(out) :: s
!
  If( associated( myargs ) ) Then
    if( i .le. ubound( myargs, 1 ) ) then
      s = myargs(i)
    else
      s = ""
    endif
  Else
    call rttov_getarg( i, s )
  Endif
End Subroutine

Integer function myiargc()
  Integer :: n
  if( associated(myargs) ) then
    n = ubound( myargs, 1 )
  else
    n = rttov_iargc()
  endif
  myiargc = n
End function

Subroutine addopt_shell( key, type, mnd, use )
  character(len=*), intent(in) :: key, type, use
  logical(kind=jplm), intent(in) :: mnd
  optional :: use, mnd
!
  character(len=argsizemax) :: str
  integer :: nn, n, n1, i1, i2, k
  character(len=argsizemax), pointer :: myargs1(:)

  myargs1 => NULL()

  If( Present( use ) ) Write( *, '("> ",A)' ) Trim(use)
  If( Present( mnd ) ) Then
    If( mnd ) Write( *, * ) "[Mandatory]"
  endif 
  Write( *, * ) "* Option: [", type, "]", " ", Trim(key)
  Read( *, '(A)' ) str

! Print *, "str = ",Trim(str)
  if( Trim(str) .ne. "" )  then
    If( type .eq. 'flag' ) Then
      nn = 0
    Else
      nn = rttov_countwords( str )
    EndIf
    n  = ubound( myargs, 1 )
    n1 = n + nn + 1

!
! realloc myargs
!
    allocate( myargs1(0:n1) )
    myargs1(0:n) = myargs(0:n)
    deallocate( myargs )
    myargs => myargs1
    myargs(n+1) = key

!
! parse argument list
!
    if( type .ne. 'flag' ) Then
      k = 1
      i1 = 1
      loop_i1 : do 
        do
          if( i1 .gt. len(str)) exit loop_i1
          if( str(i1:i1) .ne. ' ' ) exit
          i1 = i1+1
        enddo
        i2 = i1+1
        do
          if( i2 .gt. len(str)) exit
          if( str(i2:i2) .eq. ' ' ) exit
          i2 = i2+1
        enddo
!Print *, i1, i2
        myargs(n+1+k) = str(i1:i2-1)
!Print *, k, Trim(myargs(n+1+k))
        k = k+1
        i1 = i2+1
      enddo loop_i1
    endif
  endif

End Subroutine

Subroutine init_opt_seen()

  If( .not. associated( opt_seen ) ) Then
    nopt_seen = 0
    Allocate( opt_seen( 32 ) )
  EndIf

End Subroutine

Subroutine grow_opt_seen()
  Integer(Kind=jpim) :: n
  Type(rttov_opt), Pointer :: opt_seen1(:)

  n = size( opt_seen )
  If( nopt_seen .ge. n ) Then ! realloc data
    opt_seen1 => opt_seen
    allocate( opt_seen( 2 * n ) )
    opt_seen(1:nopt_seen) = opt_seen1(1:nopt_seen)
    deallocate( opt_seen1 )
  EndIf

End Subroutine

Subroutine addopt( key, type, use )
  character(len=*), intent(in) :: key, type, use
  optional :: use

  Call init_opt_seen()

  nopt_seen = nopt_seen + 1

  Call grow_opt_seen()

  opt_seen(nopt_seen)%key  = key
  opt_seen(nopt_seen)%type = type

  If( Present( use ) ) then
    opt_seen(nopt_seen)%use  = use
  else
    opt_seen(nopt_seen)%use  = ''
  endif

End Subroutine

Subroutine InitOptions( Message )
 Character(Len=*), optional, intent(in) :: Message
 integer(kind=jpim) :: n, i
 character(len=32) :: str
 n = rttov_iargc()

 allocate( myargs(0:n) )
 do i = 0, n
   call rttov_getarg( i, myargs(i) )
 enddo

 If( Present( Message ) ) Then
  message_opt = Message
 Else
  message_opt = ""
 EndIf

 if( n .eq. 1 ) then
   call mygetarg( 1_jpim, str )
   if( trim( str ) .eq. '--help' ) then
     lhelp = .true.
     return
   else if( trim( str ) .eq. '--shell' ) then
     lshell = .true.
     return
   endif
 endif

 lhelp = .false.
 allocate( check_args( n ) )
 check_args = .FALSE.

End Subroutine

Subroutine InitOptionsStr( fromstr, Message )
 Character(Len=*), optional, intent(in) :: Message
 Character(Len=*), Intent(in) :: fromstr

 integer(kind=jpim) :: n, i
 character(len=32) :: str
 character(len=len(fromstr)) :: nstr
 
 sepind(0)=0

 nstr = AdjustL(fromstr)
 call nseparators(' ', Trim(nstr), n, sepind(1:100))

 allocate( myargs(0:n) )
 myargs(0) = "fromstring"
 do i = 1, n   
   myargs(i)= AdjustL(nstr(sepind(i-1_jpim)+1:sepind(i)))
 enddo

 If( Present( Message ) ) Then
  message_opt = Message
 Else
  message_opt = ""
 EndIf

 if( n .eq. 1 ) then
   call mygetarg( 1_jpim, str )
   if( trim( str ) .eq. '--help' ) then
     lhelp = .true.
     return
   else if( trim( str ) .eq. '--shell' ) then
     lshell = .true.
     return
   endif
 endif

 lhelp = .false.
 allocate( check_args( n ) )
 check_args = .FALSE.
End Subroutine

Subroutine CheckOptions()
 integer(kind=jpim) :: i, n, is, ns, ks
 character(len=argsizemax) :: opt, prog
 logical(kind=jplm) :: pb
 character(len=10) :: fmt

 integer(kind=jpim), parameter :: lgline = 130
 integer(kind=jpim) :: halflgline = lgline/2
 character(len=lgline+14) :: buf
 character(len=3) :: buffmt

 character(len=lgline*20) line_opt
 integer(kind=jpim) :: sepind(lgline*20)
 integer(kind=jpim) :: nbsep
 integer(kind=jpim) :: l, i1, i2

 call mygetarg( 0_jpim, prog )

 if( lhelp ) then
   Print *, "Program: ", Trim(rttov_basename( prog ))
   if( Trim(message_opt) .ne. "" ) Then
     call nseparators(';', message_opt, nbsep, sepind)
     i1 = 0_jpim
     i2 = 0_jpim
     do l = 1, nbsep
       i1 = i2 + 1_jpim
       i2 = sepind(l)
       If( message_opt(i2:i2) .eq. ";" ) Then
         line_opt=message_opt(i1:i2-1)
       Else
         line_opt=message_opt(i1:i2)
       EndIf
       ns = Len(line_opt)

       do is = 1, ns / lgline
         ks = Len( Trim(line_opt(1+(is-1)*lgline:is*lgline)) )
         If( ks .gt. 0 ) Then
           If( is .eq. 1 ) Then
             Write( *, '("    ")', advance = 'no' )
           Else
             Write( *, '("  > ")', advance = 'no' )
           EndIf
           Write( fmt, '("(A",I3.3,")")' ) ks
           Write( *, fmt ) Trim(line_opt(1+(is-1)*lgline:is*lgline))
         EndIf
       enddo
     enddo
   endif

   write(buffmt, '(i3)') len(buf)

   do i = 1, nopt_seen

     if(opt_seen(i)%group) then
       Write( *, * ) 
       If( Trim(opt_seen(i)%use) .ne. "" ) &
         Write( *, * ) '* '//Trim(opt_seen(i)%use)
       cycle
     endif

     buf = ""

     write( buf, '(A32,"   ",A15)' ) &
         Trim(opt_seen(i)%key), &
         Trim(opt_seen(i)%type)

     If( Trim(opt_seen(i)%use) .ne. '' ) Then
       call nseparators(';',opt_seen(i)%use, nbsep, sepind)
       i1 = 0_jpim
       i2 = 0_jpim
       do l = 1, nbsep
         i1 = i2 + 1_jpim
         i2 = sepind(l)
       
         If( opt_seen(i)%use(i2:i2) .eq. ";" ) Then
           line_opt=opt_seen(i)%use(i1:i2-1)
         Else
           line_opt=opt_seen(i)%use(i1:i2)
         EndIf
         ns = Len(line_opt)

         do is = 1, ns / halflgline
           ks = Len(Trim(line_opt(1+(is-1)*halflgline:is*halflgline)))
           If( ks .gt. 0 ) Then
             If( is .eq. 1 .and. l .eq. 1) Then
               buf = Trim(buf)//" :   "//Trim(line_opt(1+(is-1)*halflgline:is*halflgline))
             Elseif( is.eq.1 ) Then
!                     000000000011111111112222222222333333333344444444445555555555
!                     012345678901234567890123456789012345678901234567890123456789
               buf = "                                                       "&
                     //Trim(line_opt(1+(is-1)*halflgline:is*halflgline))
             Else
               buf = "                                                     > "&
                     //Trim(line_opt(1+(is-1)*halflgline:is*halflgline))
             EndIf
             Write( *, '(A'//Trim(buffmt)//')' ) buf
           EndIf
         enddo
       enddo
     Else
       write( *, '(A'//Trim(buffmt)//')' ) buf
       write( *, * )
     EndIf

   enddo
   stop
 else if( associated( check_args ) ) then
   n = size( check_args )
   pb = .FALSE.
   do i = 1, n
     if( .not. check_args(i) ) then
       call mygetarg( i, opt )
       if( opt(1:2) .eq. '--' ) then
         print *, 'Invalid option: ', Trim(opt)
         pb = .TRUE.
         check_args(i) = .True.
       endif
     endif
   enddo

   do i = 1, n
     if( .not. check_args(i) ) then
       call mygetarg( i, opt )
       print *, 'Garbage in options:`', Trim(opt), "'"
       pb = .TRUE.
       Exit
     endif
   enddo

   if( pb ) call rttov_exit(1_jpim)

   deallocate( check_args )
 else if( lshell ) then
   open( 77, file = Trim(prog)//'.sh', form = 'formatted' )
   write( 77, '("#!/bin/sh")' )
   write( 77, * )
   write( 77, '(A)', advance = 'no' ) Trim(prog)
   n = ubound( myargs, 1 )
   do i = 1, n
     if( myargs(i) .eq. '--shell' ) cycle
     if( myargs(i)(1:2) .eq. '--' ) then
       write( 77, '(" \")' )
       write( 77, '("    ")', advance = 'no' )
     endif
     write( 77, '(" ",A)', advance = 'no' ) Trim(myargs(i))
   enddo
   write( 77, * )
   close(77)
 endif



 if( associated( opt_seen ) ) deallocate( opt_seen )
 if( associated( myargs ) ) deallocate( myargs )
End Subroutine


Subroutine Check_mnd( key, mnd, use )
 Character(Len=*),           Intent(in) :: key
 Character(Len=*), optional, Intent(in) :: use
 Logical(Kind=jplm), optional, Intent(in) :: mnd
!
 Character(len=argsizemax) :: prog

 If( Present( mnd ) ) Then
   If( mnd ) Then
     Call mygetarg( 0_jpim, prog )
     Write( *, '("PROGRAM: ",(a))' ) Trim( prog )
     Write( *, '("ERROR:   Option `",(a),"'' is mandatory")' ) Trim( key )
     If( Present( use ) ) Write( *, '("         ",(a)," : ",(a))' ) Trim( key ), Trim( use )
     Call rttov_exit(1_jpim)
   EndIf
 EndIf

End Subroutine

Subroutine FindArgIndex( key, i, n )
 Character(Len=*),      Intent(in)  :: key
 Integer(Kind=jpim), Intent(out) :: i, n
 Character(len=argsizemax) :: arg

 n = myiargc()
 Do i = 1, n
   Call mygetarg( i, arg )
   If( Trim( Adjustl(arg) ) .eq. Trim( Adjustl(key) ) ) Return
 EndDo
 i = -1_jpim
End Subroutine

Subroutine FindNextArgIndex( i, j )
 Integer(Kind=jpim), Intent(in)  :: i
 Integer(Kind=jpim), Intent(out) :: j
!
 Character(len=argsizemax) :: arg
 Integer(Kind=jpim) :: n

 n = myiargc()
 Do j = i+1, n
  Call mygetarg( j, arg )
  If( arg(1:2) .eq. '--' ) Exit
 EndDo

End Subroutine

Subroutine GetOptionS( key, val, mnd, use, exists )
!
 Character(Len=*), Intent(in)  :: key
 Character(Len=*), Intent(inout) :: val
 Logical(Kind=jplm), Intent(in), optional :: mnd
 Character(Len=*), Intent(in), optional :: use
 Logical(Kind=jplm), optional :: exists
!
 Integer(Kind=jpim) :: i, n
 Character(len=argsizemax) :: arg
 Logical(Kind=jplm) :: lshell1
 Logical(Kind=jplm) :: found

 lshell1 = lshell

 if( lhelp ) then
   Call addopt( key, 'string', use )
   return
 else if( lshell ) then
   lshell = .false.
   call addopt_shell( key, 'string', mnd, use )
 endif

 Call FindArgIndex( key, i, n )

 found = ( 0 .lt. i ) .and. ( i .lt. n )

 If( found ) Then
   if( associated( check_args ) ) Then
     check_args(i)   = .TRUE.
     check_args(i+1) = .TRUE.
   EndIf
   Call mygetarg( i+1_jpim, val )
 Else
   arg = get_env_opt( key )
   found = arg .ne. ""
   If( found ) val = arg
 EndIf

 If( .not. found ) &
   Call check_mnd( key, mnd, use )

 If( present( exists ) ) exists = found

 lshell = lshell1

End Subroutine

Subroutine GetOptionI( key, val, mnd, use, exists )
!
 Character(Len=*),      Intent(in)  :: key
 Integer(Kind=jpim), Intent(inout) :: val
 Logical(Kind=jplm), optional, Intent(in) :: mnd
 Character(Len=*), optional, Intent(in) :: use
 Logical(Kind=jplm), optional :: exists
!
 Character(len=argsizemax) :: sval
 Integer :: err
 Logical(Kind=jplm) :: lshell1

 lshell1 = lshell
 
 if( lhelp ) then
   Call addopt( key, 'integer', use )
   return
 else if( lshell ) then
   lshell = .false.
   call addopt_shell( key, 'integer', mnd, use )
 endif
 
 sval = ""
 Call GetOptionS( key, sval, mnd, use, exists )
 If( Trim( sval ) .ne. "" ) Then
   Read( sval, *, iostat = err ) val
   If( err .ne. 0 ) Then
     Print *, "Error while parsing option "//Trim(key)
     Call rttov_exit(1_jpim)
   EndIf
 EndIf

 lshell = lshell1

End Subroutine

Subroutine GetOptionR( key, val, mnd, use, exists )
!
 Character(Len=*),   Intent(in)  :: key
 Real(Kind=jprb), Intent(inout) :: val
 Logical(Kind=jplm), optional, Intent(in) :: mnd
 Character(Len=*), optional, Intent(in) :: use
 Logical(Kind=jplm), optional :: exists
!
 Character(len=argsizemax) :: sval
 Integer :: err
 Logical(Kind=jplm) :: lshell1

 lshell1 = lshell
 
 if( lhelp ) then
   Call addopt( key, 'real', use )
   return
 else if( lshell ) then
   lshell = .false.
   call addopt_shell( key, 'real', mnd, use )
 endif
 
 sval = ""
 Call GetOptionS( key, sval, mnd, use, exists )
 If( Trim( sval ) .ne. "" ) Then
   Read( sval, *, iostat = err ) val
   If( err .ne. 0 ) Then
     Print *, "Error while parsing option "//Trim(key)
     Call rttov_exit(1_jpim)
   EndIf
 EndIf

 lshell = lshell1

End Subroutine

Subroutine ReadASLFromString( val, sval )
 Character(len=*), Intent(out) :: val(:)
 Character(Len=*), Intent(in) :: sval
!
 Integer(Kind=jpim) :: i, j, k, n

 n = len( sval )

 i = 1
 k = 1
 do1 : do 
   do
     if( i .gt. n ) exit do1
     if( sval(i:i) .ne. ' ' ) exit
     i = i + 1
   enddo
   j = i
   do
     if( j .gt. n ) exit 
     if( sval(j:j) .eq. ' ' ) exit
     j = j + 1
   enddo

   val(k) = sval(i:j-1)
   i = j
   k = k + 1
 enddo do1


End Subroutine

Subroutine ReadSLFromString( val, sval )
 Character(len=*), Pointer :: val(:)
 Character(Len=*), Intent(in) :: sval
!
 Integer(Kind=jpim) :: n

 n = rttov_countwords( sval )
 allocate( val( n ) )

 Call ReadASLFromString( val, sval )

End Subroutine

Subroutine ReadSLFromFile( val, sval )
 Character(len=*), Pointer :: val(:)
 Character(Len=*), Intent(in) :: sval
!
 Integer(Kind=jpim) :: k, n
 Integer(Kind=jpim) :: ioerr
 Character(len=4096) :: buffer


 Open( 77, file = Trim(sval), form = 'formatted', status = 'old', iostat = ioerr )
 If( ioerr .ne. 0 ) Then
   Print '( "Could not open ",A, " for reading")', Trim(sval)
   Call rttov_exit(1_jpim)
 EndIf
 n = 0_jpim
 Do 
   Read( 77, '(A)', end = 500 ) buffer
   n = n + Rttov_CountWords( buffer )
 EndDo

 500 Continue

 Rewind( 77 )

 Allocate( val( n ) )
 
 k = 1
 Do 
   Read( 77, '(A)', end = 600 ) buffer
   n = rttov_countwords( buffer )
   Call ReadASLFromString( val(k:k+n-1), buffer )
   k = k + n
 EndDo

 600 Continue


 Close( 77 )
 
End Subroutine

Subroutine GetOptionSL( key, val, mnd, use, exists )
!
 Character(Len=*), Intent(in) :: key
 Character(len=*), Pointer :: val(:)
 Logical(Kind=jplm), optional, Intent(in) :: mnd
 Character(Len=*), optional, Intent(in) :: use
 Logical(Kind=jplm), optional :: exists
!
 Integer(Kind=jpim) :: i, j, k, n
 Character(len=argsizemax) :: arg
 Character(len=argsizemax) :: sval
 Logical(Kind=jplm) :: lshell1
 Logical(Kind=jplm) :: found

 lshell1 = lshell
 
 if( lhelp ) then
   Call addopt( key, 'string-list', use )
   return
 else if( lshell ) then
   lshell = .false.
   call addopt_shell( key, 'string-list', mnd, use )
 endif
 
 Call FindArgIndex( key, i, n )

 found = i >= 0

 If( found ) Then

   Call FindNextArgIndex( i, j )

   Allocate( val( j - i - 1 ) )
   
   if( associated( check_args ) ) &
     check_args(i) = .TRUE.
   
   Do k = i+1, j-1
     if( associated( check_args ) ) &
       check_args(k)   = .TRUE.
     Call mygetarg( k, arg )
     val(k-i) = arg
   EndDo

   If((size(val) .eq. 1).and.(val(1)(1:7).eq.'file://')) Then
     arg = val(1)(8:)
     deallocate(val)
     Call ReadSLFromFile( val, arg )
   EndIf

 EndIf
 
 If(.not. found) Then
   sval = get_env_opt( key )
   found = sval .ne. ""
   if( found ) &
     call ReadSLFromString( val, sval )
 EndIf

 If( .not. found ) &
   Call check_mnd( key, mnd, use )

 If( present( exists ) ) exists = found

 lshell = lshell1

End Subroutine

Subroutine GetOptionIL( key, val, mnd, use, exists )
!
 Character(Len=*),      Intent(in) :: key
 Integer(Kind=jpim), Pointer    :: val(:)
 Logical(Kind=jplm), optional, Intent(in) :: mnd
 Character(Len=*), optional, Intent(in) :: use
 Logical(Kind=jplm), optional :: exists
!
 Character(len=argsizemax), Pointer :: sval(:) => NULL()
 Integer(kind=jpim) :: i, n
 Integer :: err
 Logical(Kind=jplm) :: lshell1

 lshell1 = lshell

 if( lhelp ) then
   Call addopt( key, 'integer-list', use )
   return
 else if( lshell ) then
   lshell = .false.
   call addopt_shell( key, 'integer-list', mnd, use )
 endif
 
 Call GetOptionSL( key, sval, mnd, use, exists )

 If( .Not. Associated( sval ) ) goto 999

 n = Size( sval )
 Allocate( val( n ) )
 Do i = 1, n
   Read( sval( i ), *, iostat = err ) val( i )
   If( err .ne. 0 ) Then
     Print *, "Error while parsing option "//Trim(key)
     Call rttov_exit(1_jpim)
   EndIf
 EndDo

 Deallocate( sval )

999 continue
 lshell = lshell1

End Subroutine

Subroutine GetOptionRL( key, val, mnd, use, exists )
!
 Character(Len=*),   Intent(in) :: key
 Real(Kind=jprb), Pointer    :: val(:)
 Logical(Kind=jplm), optional, Intent(in) :: mnd
 Character(Len=*), optional, Intent(in) :: use
 Logical(Kind=jplm), optional :: exists
!
 Character(len=argsizemax), Pointer :: sval(:) => NULL()
 Integer(kind=jpim) :: i, n
 Integer :: err
 Logical(Kind=jplm) :: lshell1

 lshell1 = lshell

 if( lhelp ) then
   Call addopt( key, 'real-list', use )
   return
 else if( lshell ) then
   lshell = .false.
   call addopt_shell( key, 'real-list', mnd, use )
 endif
 
 Call GetOptionSL( key, sval, mnd, use, exists )

 If( .Not. Associated( sval ) ) goto 999

 n = Size( sval )
 Allocate( val( n ) )
 Do i = 1, n
   Read( sval( i ), *, iostat = err ) val( i )
   If( err .ne. 0 ) Then
     Print *, "Error while parsing option "//Trim(key)
     Call rttov_exit(1_jpim)
   EndIf
 EndDo

 Deallocate( sval )

999 continue
 lshell = lshell1

End Subroutine

Subroutine GetOptionB( key, val, use, exists )
!
 Character(Len=*), Intent(in)  :: key
 Logical(Kind=jplm), Intent(inout) :: val
 Character(Len=*), optional, Intent(in) :: use
 Logical(Kind=jplm), optional :: exists
!
 Logical(Kind=jplm) :: lshell1
 Logical(Kind=jplm) :: found
 Character(len=argsizemax) :: sval
 Integer(Kind=jpim) :: i, n
 
 lshell1 = lshell

 if( lhelp ) then
   Call addopt( key, '(flag)', use )
   return
 else if( lshell ) then
   lshell = .false.
   call addopt_shell( key, '(flag)', .false._jplm, use )
 endif
 
 Call FindArgIndex( key, i, n )
 found = i > 0
 if( found .and. associated( check_args ) ) Then
   check_args(i)   = .TRUE.
   val = .true.
 Else
   sval = get_env_opt( key )
   If( sval .ne. "" ) &
     Read( sval, * ) val
 EndIf

 If( present( exists ) ) exists = found

 lshell = lshell1

End Subroutine



Subroutine nseparators(sep, line, nbsep, sepind)

 character(len=1), intent(in) :: sep
 character(len=*), intent(in) :: line
 integer, intent(out) :: nbsep
 integer, intent(out) :: sepind(:)

 integer(kind=jpim) :: l, last, ns

 nbsep     = 0_jpim
 sepind(:) = 0_jpim

 if( Trim(line) .ne. "" ) Then
   ns = Len(line)

   l         = 1_jpim
   last      = 0_jpim

   do while ( l .gt. 0_jpim )
     l = index(line(last+1_jpim:ns), sep)
     if (l .gt. 0) then
       nbsep = nbsep + 1_jpim
       sepind(nbsep) = l + last
       last = sepind(nbsep)
       if( nbsep .gt. 1 ) then
         if( last .eq. sepind(nbsep-1)+1 ) then
           sepind(nbsep) = 0
           nbsep = nbsep -1
           !last = sepind(nbsep)
         endif
       endif
     endif
   enddo
   
   if( last .ne. ns ) then
     nbsep = nbsep + 1_jpim
     sepind(nbsep) = ns
   endif
   
 endif
 
End Subroutine


End Module
