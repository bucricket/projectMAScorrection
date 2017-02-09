!
MODULE mod_rttov_emis_atlas
  ! Description:
  !   module for emissivity atlas.
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
  USE parkind1, Only : jpim, jplm

  ! Disable implicit typing
  IMPLICIT NONE

  ! Atlas version is set up by rttov_setup_emis_atlas
  INTEGER(KIND=jpim) :: ir_atlas_version   ! Version of atlas for IR

  INTEGER(KIND=jpim) :: mw_atlas_version   ! Version of atlas for MW

  ! Flags to indicate whether atlases have been initialised
  LOGICAL(KIND=jplm) :: ir_atlas_init = .FALSE.
  LOGICAL(KIND=jplm) :: mw_atlas_init = .FALSE.

  ! Flag to indicate if error (stdev) data for the IR atlas was read in
  LOGICAL(KIND=jplm) :: ir_atlas_std_init = .FALSE.

  ! Flag to indicate if IR atlas was initialised for a single instrument
  LOGICAL(KIND=jplm) :: ir_atlas_single_inst = .FALSE.

  ! Flag to indicate if angular correction should be applied in IR atlas
  LOGICAL(KIND=jplm) :: ir_atlas_do_ang_corr = .FALSE.

END MODULE  mod_rttov_emis_atlas
