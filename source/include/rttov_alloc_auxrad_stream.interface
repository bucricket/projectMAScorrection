Interface
SUBROUTINE rttov_alloc_auxrad_stream( &
            & ERR,           &
            & auxrad_stream, &
            & opts,          &
            & nstreams,      &
            & nlayers,       &
            & nchannels,     &
            & asw,           &
            & init)
  USE parkind1, ONLY : jpim, jplm
  USE rttov_types, ONLY : rttov_options, radiance_aux
  IMPLICIT NONE
  INTEGER(KIND=jpim) , INTENT(OUT)             :: ERR
  TYPE(radiance_aux ), INTENT(INOUT)           :: auxrad_stream
  TYPE(rttov_options), INTENT(IN)              :: opts
  INTEGER(KIND=jpim) , INTENT(IN)              :: nstreams
  INTEGER(KIND=jpim) , INTENT(IN)              :: nlayers
  INTEGER(KIND=jpim) , INTENT(IN)              :: nchannels
  INTEGER(KIND=jpim) , INTENT(IN)              :: asw
  LOGICAL(KIND=jplm) , INTENT(IN)   , OPTIONAL :: init
End Subroutine
End Interface
