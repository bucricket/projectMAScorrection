Interface
Subroutine rttov_check_traj( &
       & ERR,             &
       & nprofiles,       &
       & nchannels,       &
       & opts,            &
       & nlevels,         &
       & coefs,           &
       & asw,             &
       & traj0,    traj1,    traj2,    &
       & traj0_tl, traj1_tl, traj2_tl, &
       & traj0_ad, traj1_ad, traj2_ad, &
       & traj0_k,  traj1_k,  traj2_k   &
       & )
  Use parkind1, only : jpim
  Use rttov_types, Only : &
          & rttov_coefs,    &
          & rttov_options,  &
          & rttov_traj
  Implicit None
  Integer(Kind=jpim),       Intent(in)    :: nprofiles
  Integer(Kind=jpim),       Intent(in)    :: nchannels
  Type(rttov_coefs),       Intent(in), Target :: coefs
  Type(rttov_options),     Intent(in)    :: opts
  Integer(Kind=jpim),       Intent(in)   :: nlevels
  Integer(Kind=jpim),       Intent(out)   :: ERR
  Integer(Kind=jpim),       Intent(in)    :: asw
  Type(rttov_traj),   Pointer, Optional                :: traj0, traj0_tl, traj0_ad, traj0_k
  Type(rttov_traj),   Target,  Optional, Intent(inout) :: traj1, traj1_tl, traj1_ad, traj1_k
  Type(rttov_traj),   Target,  Optional, Intent(inout) :: traj2, traj2_tl, traj2_ad, traj2_k
End Subroutine
End Interface
