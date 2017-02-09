! Description:
!> @file
!!   Example which converts obs to PC-space.
!
!> @brief
!!   Example program demonstrating how observations may
!!   be converted to PC-space for assimilation/retrieval
!!   applications using PC-RTTOV.
!!
!! @details
!!   For usage details see user guide or run:
!!   $ rttov_obs_to_pc.exe
!!
!!   See the source code for more details: this is intended
!!   as example code which should be modified for your own
!!   purposes.
!
PROGRAM rttov_obs_to_pc

#include "throw.h"

  USE rttov_types, ONLY : &
    rttov_options,        &
    rttov_coefs
  USE parkind1, ONLY : jplm, jpim, jprb
  USE rttov_getoptions, ONLY : getoption, initoptions
  USE rttov_unix_env, ONLY : rttov_iargc, rttov_exit
  USE rttov_math_mod, ONLY : planck

  IMPLICIT NONE

#include "rttov_read_coefs.interface"
#include "rttov_dealloc_coefs.interface"

  INTEGER(jpim), PARAMETER :: lu = 70

  TYPE(rttov_options)      :: opts
  TYPE(rttov_coefs)        :: coefs

  CHARACTER(256)           :: rtcoef_file
  CHARACTER(256)           :: pccoef_file
  CHARACTER(256)           :: obs_file

  INTEGER(jpim)            :: err
  INTEGER(jpim)            :: i, j
  INTEGER(jpim)            :: nchannels
  INTEGER(jpim)            :: npcscores, ipcbnd, ipcreg
  LOGICAL(jplm)            :: obs_bt
  LOGICAL(jplm)            :: exists

  REAL(jprb), ALLOCATABLE  :: obs(:), pcscores(:), bt_eff(:)

! --End of header--------------------------------------------------------------

TRY

  ! ===========================================================================
  ! This is demonstration code: you would typically want to adapt this for
  ! your own uses.

  ! This code reads in a single observation from an ASCII file: this consists
  ! of radiances or BTs for *all* instrument channels separated by white-space.

  ! You must specify the same PC band (ipcbnd) and regression predictor set
  ! (ipcreg) used in your PC-RTTOV calculations.

  ! The number of PC scores required is specified on the commandline.

  ! The computed PC scores are written to the commandline.
  ! ===========================================================================

  ! ----------------------------------------
  ! Handle program arguments
  ! ----------------------------------------

  IF (rttov_iargc() == 0) THEN
    PRINT *, 'Usage: --rtcoef_file ...  input rtcoef file name'
    PRINT *, '       --pccoef_file ...  input pccoef file name'
    PRINT *, '       --obs_file    ...  input ASCII file containing observations'
    PRINT *, '       --ipcbnd      ...  PC band'
    PRINT *, '       --ipcreg      ...  PC regression set'
    PRINT *, '       --npcscores   ...  number of PC scores to output'
    PRINT *, '       --obs_bt           optional flag to set obs units'
    PRINT *, '  Observations values should be separated by white-space.'
    PRINT *, '  Units are Kelvin if --obs_bt is present, otherwise mW/m-2/sr-1/cm-1'
    STOP
  ENDIF

  CALL initoptions()

  CALL getoption("--rtcoef_file", rtcoef_file, mnd=.TRUE._jplm)
  INQUIRE(FILE=rtcoef_file, EXIST=exists)
  IF (.NOT. exists) THEN
    PRINT *, 'Cannot find rtcoef file: '//TRIM(rtcoef_file)
    STOP
  ENDIF

  CALL getoption("--pccoef_file", pccoef_file, mnd=.TRUE._jplm)
  INQUIRE(FILE=pccoef_file, EXIST=exists)
  IF (.NOT. exists) THEN
    PRINT *, 'Cannot find pccoef file: '//TRIM(pccoef_file)
    STOP
  ENDIF

  CALL getoption("--obs_file", obs_file, mnd=.TRUE._jplm)
  INQUIRE(FILE=obs_file, EXIST=exists)
  IF (.NOT. exists) THEN
    PRINT *, 'Cannot find observations file: '//TRIM(obs_file)
    STOP
  ENDIF

  CALL getoption("--ipcbnd", ipcbnd, mnd=.TRUE._jplm)
  CALL getoption("--ipcreg", ipcreg, mnd=.TRUE._jplm)
  CALL getoption("--npcscores", npcscores, mnd=.TRUE._jplm)

  obs_bt = .FALSE.
  CALL getoption( "--obs_bt", obs_bt)

  ! ----------------------------------------
  ! Set options
  ! ----------------------------------------

  opts%rt_ir%pc%addpc = .TRUE.
  opts%rt_ir%pc%ipcbnd = ipcbnd
  opts%rt_ir%pc%ipcreg = ipcreg

  ! ----------------------------------------
  ! Read coefficients and observations
  ! ----------------------------------------

  CALL rttov_read_coefs(err, coefs, opts, file_coef=rtcoef_file, file_pccoef=pccoef_file)
  THROWM(err .NE. 0, 'Error reading coefficients')

  nchannels = coefs%coef%fmv_chn
  IF (nchannels /= MAXVAL(coefs%coef%ff_ori_chn)) THEN
    err = errorstatus_fatal
    THROWM(err .NE. 0, 'The full coefficient file containing all instrument channels is required')
  ENDIF

  OPEN(lu, FILE=obs_file, IOSTAT=err)
  THROWM(err .NE. 0, 'Error opening observations file')

  ALLOCATE(obs(nchannels))
  READ(lu, *, IOSTAT=err) obs
  THROWM(err .NE. 0, 'Error reading observations')

  CLOSE(lu)

  ! ----------------------------------------
  ! Convert BTs to radiances (if necessary)
  ! ----------------------------------------

  ! Coefficient structure contains correction scaling and offset for BTs
  ! and constants for use in Planck function for each channel

  IF (obs_bt) THEN
    ALLOCATE(bt_eff(nchannels))
    DO j = 1, nchannels
      bt_eff(j) = coefs%coef%ff_bco(j) + coefs%coef%ff_bcs(j) * obs(j)
      CALL planck(coefs%coef%planck1(j), coefs%coef%planck2(j), bt_eff(j), obs(j))
    ENDDO
    DEALLOCATE(bt_eff)
  ENDIF

  ! ----------------------------------------
  ! Convert observations to PCs
  ! ----------------------------------------

  ALLOCATE(pcscores(npcscores))
  pcscores = 0._jprb

  DO i = 1, npcscores
    DO j = 1, nchannels
      pcscores(i)= pcscores(i) + &
        coefs%coef_pccomp%eigen(ipcbnd)%eigenvectors(j,i) * &
        obs(j) / coefs%coef_pccomp%noise_in(j)
    ENDDO
  ENDDO

  ! ----------------------------------------
  ! Write out the PC scores
  ! ----------------------------------------

  WRITE(*,'(a)') 'PC scores:'
  WRITE(*,'(6e13.5)') pcscores

  ! ----------------------------------------
  ! Tidy up
  ! ----------------------------------------

  IF (ALLOCATED(obs))      DEALLOCATE(obs)
  IF (ALLOCATED(pcscores)) DEALLOCATE(pcscores)

  CALL rttov_dealloc_coefs(err, coefs)
  THROWM(err .NE. 0, 'Error deallocating coefficients')

PCATCH
END PROGRAM rttov_obs_to_pc