!
Module rttov_global
  ! Description:
  ! Definition of global variables) for RTTOV
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
  !    Copyright 2002, EUMETSAT, All Rights Reserved.
  !
  ! History:
  ! Version   Date     Comment
  ! -------   ----     -------
  !  1.0       10/01/2003  Original (P Brunel)
  !            07/12/2007 add variable for checkinput P.B.
  !

  Use parkind1, Only : jpim, jplm
  Implicit None
  !- End of header --------------------------------------------------------

  !1. error reporting
  !------------------
  Integer(Kind=jpim) :: error_unit   ! logical unit for the error messages

  Logical(Kind=jplm) :: err_init     ! true if module already initialised


  Data err_init /.false./

End Module rttov_global
