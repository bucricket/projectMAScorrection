Subroutine rttov_checkpcchan( &
       & nprofiles,         &
       & nchannels,         &
       & opts,              &
       & chanprof,          &
       & coefs,             &
       & ERR                &
       & )
! Description:
!   Checks input channel list against PC regression channel set
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
!INTF_OFF
#include "throw.h"
!INTF_ON

Use parkind1, only : jpim

Use rttov_types, only: rttov_options, rttov_chanprof, rttov_coefs

!INTF_OFF
Use parkind1, Only : &
    jpim

!INTF_ON

Implicit None
INTEGER(KIND=jpim)     , INTENT(IN)    :: nprofiles
INTEGER(KIND=jpim)     , INTENT(IN)    :: nchannels
TYPE(rttov_options)    , INTENT(IN)    :: opts
TYPE(rttov_chanprof)   , INTENT(IN)    :: chanprof(:)
TYPE(rttov_coefs   )   , INTENT(IN)    :: coefs
INTEGER(KIND=jpim)     , INTENT(OUT)   :: ERR
!INTF_END

INTEGER(KIND=jpim) :: prof
INTEGER(KIND=jpim) :: i
INTEGER(KIND=jpim) :: dc(nchannels/nprofiles)
!-----End of header-------------------------------------------------------------

TRY

    ! test if the number of channels per profile is compatible
    ! with the regression set
    IF( nchannels/nprofiles .ne. coefs%coef_pccomp%pcreg(opts%rt_ir%pc%ipcbnd,opts%rt_ir%pc%ipcreg)%fmv_pc_npred ) THEN
      ERR = 1_jpim
    ENDIF
    THROWM( ERR .NE. 0 , "PC calc; invalid number of regression channels")

    IF ( chanprof(1)%chan .eq. 1 .AND. chanprof(nchannels)%chan .eq. nchannels/nprofiles ) THEN

      ! test if channels are [1..N] 
      DO prof = 1, nprofiles
        i=1 + (prof-1)*(nchannels/nprofiles)
        dc = coefs%coef%ff_ori_chn(chanprof(i:i+(nchannels/nprofiles)-1)%chan) - &
           & coefs%coef_pccomp%pcreg(opts%rt_ir%pc%ipcbnd,opts%rt_ir%pc%ipcreg)%predictindex(:)
        IF( ANY (dc .ne. 0_jpim )) THEN
          ERR = 1_jpim
        ENDIF
      ENDDO

    ELSE

      ! test if channels are exactly the same as the predictors
      DO prof = 1, nprofiles
        i=1 + (prof-1)*(nchannels/nprofiles)
        dc = chanprof(i:i+(nchannels/nprofiles)-1)%chan - &
           & coefs%coef_pccomp%pcreg(opts%rt_ir%pc%ipcbnd,opts%rt_ir%pc%ipcreg)%predictindex(:)
        IF( ANY (dc .ne. 0_jpim )) THEN
          ERR = 1_jpim
        ENDIF
      ENDDO

    ENDIF

    THROWM( ERR .NE. 0 , "PC calc; invalid regression channels indices")

CATCH
End Subroutine
