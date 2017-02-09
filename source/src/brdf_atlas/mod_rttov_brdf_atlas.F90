!
MODULE mod_rttov_brdf_atlas
  ! Description:
  !   module for brdf atlas.
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
  !  1.0      11/01/2012  Created (J. Vidot)
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

  ! Atlas version is set up by rttov_setup_brdf_atlas
  INTEGER(KIND=jpim) :: vn_atlas_version   ! Version of atlas for vis/nir

  ! Flags to indicate whether atlases have been initialised
  LOGICAL(KIND=jplm) :: vn_atlas_init = .FALSE.

  ! Flag to indicate if BRDF atlas was initialised for a single instrument
  LOGICAL(KIND=jplm) :: brdf_atlas_single_inst = .FALSE.

END MODULE  mod_rttov_brdf_atlas
