! Description:
!> @file
!!   Program to generate climatological aerosol profiles.
!
!> @brief
!!   Program to generate climatological aerosol profiles.
!!
!! @details
!!   Given input pressure, temperature and water vapour
!!   profiles, generates 10 sets of climatological aerosol
!!   profiles in output file 'prof_aerosl_cl.dat'.
!!
!!   The input files plevs.dat and prof.dat must be in the
!!   current directory. Examples of these files can be found
!!   in the RTTOV data/ directory. The program prompts for
!!   user input for the other parameters.
!!
!!   Inputs:
!!   plevs.dat - number of levels followed by pressure profile (hPa)
!!   prof.dat - temperature (K) and water vapour (ppmv) profiles
!!   latitude - latitude of profile (-90 to +90 degrees)
!!   elevation - surface elevation (km)
!!   surface level - level number corresponding to the surface
!!   scale factor - factor by which to scale output aerosol profiles
!!
!!   See the rttov_aer_clim_prof.F90 subroutine for more information.
!
PROGRAM create_aer_clim_prof

  USE parkind1, ONLY : jprb, jpim
  USE rttov_const, ONLY: naer_max

  IMPLICIT NONE

  INTEGER(jpim)           :: i, j
  INTEGER(jpim)           :: levsurf
  INTEGER(jpim)           :: nlev
  REAL(jprb)              :: scalefactor
  REAL(jprb)              :: latitude
  REAL(jprb)              :: elevation
  REAL(jprb), ALLOCATABLE :: p(:), t(:), q(:)
  REAL(jprb), ALLOCATABLE :: aerprof(:,:,:)

#include "rttov_aer_clim_prof.interface"

!-----End of header-------------------------------------------------------------

  WRITE(6,*)'ENTER LATITUDE'
  READ (5,*)latitude
  WRITE(6,*)'ENTER ELEVATION'
  READ (5,*)elevation
  WRITE(6,*)'ENTER SURFACE LEVEL'
  READ (5,*)levsurf
  WRITE(6,*)'ENTER SCALING FACTOR'
  READ (5,*)scalefactor

  OPEN(10,FILE='prof_aerosl_cl.dat')

  OPEN(1,FILE='plevs.dat')

  READ(1,*)nlev

  ALLOCATE(p(nlev), t(nlev), q(nlev), aerprof(nlev-1, 10, naer_max))

  READ(1,*)(p(i),i=1,nlev)

  OPEN(2,FILE='prof.dat')

  READ(2,*)(t(i),i=1,nlev)
  READ(2,*)(q(i),i=1,nlev)

  CALL rttov_aer_clim_prof(p, t, q, levsurf, latitude, elevation, scalefactor, aerprof)

  DO j = 1, 10
    DO i = 1, nlev-1
      WRITE(10,'(13f11.3)') aerprof(i,j,:)
    ENDDO
  ENDDO

  DEALLOCATE(p, t, q, aerprof)

END PROGRAM create_aer_clim_prof

