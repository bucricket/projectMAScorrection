! compute optical parameters for Ice cloud layer using Baran data
SUBROUTINE rttov_baran2013_calc_optpar_tl (optp, ichn, &
                       & T_in, IWC_in, T_tl_in, IWC_tl_in, &
                       & abso, sca, bpr, asym ,&
                       & abso_tl, sca_tl, bpr_tl, asym_tl)

! Description:
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
!  1.0       Oct 2012    Original (P. Brunel J. Vidot
!  1.1       Mar 2013    New parameterization with 6 parameters 
!			 for absorption and scattering coefficients (J.Vidot)
!  1.2       May 2013    Add iwc and temperature test from baran database (J. Vidot)   
!
! Code Description:
!   Language:           Fortran 90.
!   Software Standards: "European Standards for Writing and
!     Documenting Exchangeable Fortran 90 Code".
!
! Declarations:
! Modules used:
! Imported Parameters:
! Imported Type Definitions:
  USE rttov_types, ONLY :  &
       & rttov_optpar_ir
  USE parkind1, ONLY : jpim, jprb
!INTF_OFF 
  USE mod_rttov_baran2013_icldata
!INTF_ON
  IMPLICIT NONE

  TYPE(rttov_optpar_ir) ,INTENT(IN)  :: optp
  INTEGER(KIND=jpim)    ,INTENT(IN)  :: ichn
  REAL(KIND=JPRB)       ,INTENT(IN)  :: t_in, t_tl_in
  REAL(KIND=JPRB)       ,INTENT(IN)  :: iwc_in, iwc_tl_in
  REAL(KIND=JPRB)       ,INTENT(OUT) :: abso, abso_tl
  REAL(KIND=JPRB)       ,INTENT(OUT) :: sca, sca_tl
  REAL(KIND=JPRB)       ,INTENT(OUT) :: bpr, bpr_tl
  REAL(KIND=JPRB)       ,INTENT(OUT) :: asym, asym_tl

!INTF_END

! Local Scalars:
  REAL(KIND=JPRB)  :: absi, absj, bpri, asymi
  REAL(KIND=JPRB)  :: scai, scaj, bprj, asymj
  REAL(KIND=JPRB)  :: absi_tl, absj_tl, bpri_tl, asymi_tl
  REAL(KIND=JPRB)  :: scai_tl, scaj_tl, bprj_tl, asymj_tl


  REAL(KIND=JPRB)     :: dx_dwn
  REAL(KIND=JPRB)     :: liwc, liwc_tl
  INTEGER(KIND=JPIM)  :: iwn, jwn
  REAL(KIND=JPRB)     :: iwc,t,iwc_tl,t_tl

!- End of header --------------------------------------------------------

    if (iwc_in < baran2013_iwc_min) then
      iwc = baran2013_iwc_min
      iwc_tl = 0._JPRB
    else if (iwc_in > baran2013_iwc_max) then
      iwc = baran2013_iwc_max
      iwc_tl = 0._JPRB
    else
      iwc = iwc_in
      iwc_tl = iwc_tl_in
    endif

    if (T_in < baran2013_temp_min) then 
      T = baran2013_temp_min
      T_tl = 0._JPRB
    else if (T_in > baran2013_temp_max) then 
      T = baran2013_temp_max
      T_tl = 0._JPRB
    else
      T = T_in
      T_tl = T_tl_in
    endif

    LIWC = LOG10(iwc)
    LIWC_TL = iwc_tl / ( iwc * log(10.0_JPRB) )
    
    iwn = optp%optpiclb%iwn2013(ichn)  
    jwn = optp%optpiclb%jwn2013(ichn)  
    dx_dwn = optp%optpiclb%dx_dwn2013(ichn)  

    absi = baran2013_regcoef_abs(1_JPIM,iwn)    + &
         & baran2013_regcoef_abs(2_JPIM,iwn)*T + &
         & baran2013_regcoef_abs(3_JPIM,iwn)*LIWC + &
         & baran2013_regcoef_abs(4_JPIM,iwn)*T*T + &
         & baran2013_regcoef_abs(5_JPIM,iwn)*LIWC*LIWC + &
         & baran2013_regcoef_abs(6_JPIM,iwn)*T*LIWC
    absi_tl =  &
         & baran2013_regcoef_abs(2_JPIM,iwn)*T_tl + &
         & baran2013_regcoef_abs(3_JPIM,iwn)*LIWC_tl + &
         & baran2013_regcoef_abs(4_JPIM,iwn)*2.0_JPRB*T*T_tl + &
         & baran2013_regcoef_abs(5_JPIM,iwn)*2.0_JPRB*LIWC*LIWC_tl + &
         & baran2013_regcoef_abs(6_JPIM,iwn)*(T*LIWC_tl + T_tl*LIWC )

    absj = baran2013_regcoef_abs(1_JPIM,jwn)    + &
         & baran2013_regcoef_abs(2_JPIM,jwn)*T + &
         & baran2013_regcoef_abs(3_JPIM,jwn)*LIWC + &
         & baran2013_regcoef_abs(4_JPIM,jwn)*T*T + &
         & baran2013_regcoef_abs(5_JPIM,jwn)*LIWC*LIWC + &
         & baran2013_regcoef_abs(6_JPIM,jwn)*T*LIWC
    absj_tl =  &
         & baran2013_regcoef_abs(2_JPIM,jwn)*T_tl + &
         & baran2013_regcoef_abs(3_JPIM,jwn)*LIWC_tl + &
         & baran2013_regcoef_abs(4_JPIM,jwn)*2.0_JPRB*T*T_tl + &
         & baran2013_regcoef_abs(5_JPIM,jwn)*2.0_JPRB*LIWC*LIWC_tl + &
         & baran2013_regcoef_abs(6_JPIM,jwn)*(T*LIWC_tl + T_tl*LIWC )
         
    absi = 10**absi
    absi_tl = absi * absi_tl * log(10.0_JPRB)
    absj = 10**absj
    absj_tl = absj * absj_tl * log(10.0_JPRB)

    abso = absi + (absj - absi) * dx_dwn
    abso_tl = absi_tl + (absj_tl - absi_tl) * dx_dwn

    
    scai = baran2013_regcoef_sca(1_JPIM,iwn)    + &
         & baran2013_regcoef_sca(2_JPIM,iwn)*T + &
         & baran2013_regcoef_sca(3_JPIM,iwn)*LIWC + &
         & baran2013_regcoef_sca(4_JPIM,iwn)*T*T + &
         & baran2013_regcoef_sca(5_JPIM,iwn)*LIWC*LIWC + &
         & baran2013_regcoef_sca(6_JPIM,iwn)*T*LIWC
    scai_tl =  &
         & baran2013_regcoef_sca(2_JPIM,iwn)*T_tl + &
         & baran2013_regcoef_sca(3_JPIM,iwn)*LIWC_tl + &
         & baran2013_regcoef_sca(4_JPIM,iwn)*2.0_JPRB*T*T_tl + &
         & baran2013_regcoef_sca(5_JPIM,iwn)*2.0_JPRB*LIWC*LIWC_tl + &
         & baran2013_regcoef_sca(6_JPIM,iwn)*(T*LIWC_tl + T_tl*LIWC )

    scaj = baran2013_regcoef_sca(1_JPIM,jwn)    + &
         & baran2013_regcoef_sca(2_JPIM,jwn)*T + &
         & baran2013_regcoef_sca(3_JPIM,jwn)*LIWC + &
         & baran2013_regcoef_sca(4_JPIM,jwn)*T*T + &
         & baran2013_regcoef_sca(5_JPIM,jwn)*LIWC*LIWC + &
         & baran2013_regcoef_sca(6_JPIM,jwn)*T*LIWC
    scaj_tl =  &
         & baran2013_regcoef_sca(2_JPIM,jwn)*T_tl + &
         & baran2013_regcoef_sca(3_JPIM,jwn)*LIWC_tl + &
         & baran2013_regcoef_sca(4_JPIM,jwn)*2.0_JPRB*T*T_tl + &
         & baran2013_regcoef_sca(5_JPIM,jwn)*2.0_JPRB*LIWC*LIWC_tl + &
         & baran2013_regcoef_sca(6_JPIM,jwn)*(T*LIWC_tl + T_tl*LIWC )

    scai = 10**scai
    scai_tl = scai * scai_tl * log(10.0_JPRB)
    scaj = 10**scaj
    scaj_tl = scaj * scaj_tl * log(10.0_JPRB)

    sca = scai + (scaj - scai) * dx_dwn
    sca_tl = scai_tl + (scaj_tl - scai_tl) * dx_dwn


    bpri = baran2013_regcoef_bpr(1_JPIM,iwn)    + &
         & baran2013_regcoef_bpr(2_JPIM,iwn)*T + &
         & baran2013_regcoef_bpr(3_JPIM,iwn)*LIWC
    bpri_tl =  &
         & baran2013_regcoef_bpr(2_JPIM,iwn)*T_tl + &
         & baran2013_regcoef_bpr(3_JPIM,iwn)*LIWC_tl
    bprj = baran2013_regcoef_bpr(1_JPIM,jwn)    + &
         & baran2013_regcoef_bpr(2_JPIM,jwn)*T + &
         & baran2013_regcoef_bpr(3_JPIM,jwn)*LIWC
    bprj_tl =  &
         & baran2013_regcoef_bpr(2_JPIM,jwn)*T_tl + &
         & baran2013_regcoef_bpr(3_JPIM,jwn)*LIWC_tl
    if( bpri < baran2013_corcoef_bpr(1_JPIM,iwn) ) then
      bpri = ( bpri - baran2013_corcoef_bpr(2_JPIM,iwn) ) / baran2013_corcoef_bpr(3_JPIM,iwn)
      bpri_tl = bpri_tl / baran2013_corcoef_bpr(3_JPIM,iwn)
    elseif( bpri > baran2013_corcoef_bpr(4_JPIM,iwn) ) then
      bpri = ( bpri - baran2013_corcoef_bpr(5_JPIM,iwn) ) / baran2013_corcoef_bpr(6_JPIM,iwn)
      bpri_tl = bpri_tl / baran2013_corcoef_bpr(6_JPIM,iwn)
    endif
    if( bprj < baran2013_corcoef_bpr(1_JPIM,jwn) ) then
      bprj = ( bprj - baran2013_corcoef_bpr(2_JPIM,jwn) ) / baran2013_corcoef_bpr(3_JPIM,jwn)
      bprj_tl = bprj_tl / baran2013_corcoef_bpr(3_JPIM,jwn)
    elseif( bprj > baran2013_corcoef_bpr(4_JPIM,jwn) ) then
      bprj = ( bprj - baran2013_corcoef_bpr(5_JPIM,jwn) ) / baran2013_corcoef_bpr(6_JPIM,jwn)
      bprj_tl = bprj_tl / baran2013_corcoef_bpr(6_JPIM,jwn)
    endif

    bpr = bpri + (bprj - bpri) * dx_dwn
    bpr_tl = bpri_tl + (bprj_tl - bpri_tl) * dx_dwn


    asymi = baran2013_regcoef_asym(1_JPIM,iwn)    + &
         & baran2013_regcoef_asym(2_JPIM,iwn)*T + &
         & baran2013_regcoef_asym(3_JPIM,iwn)*LIWC
    asymi_tl =  &
         & baran2013_regcoef_asym(2_JPIM,iwn)*T_tl + &
         & baran2013_regcoef_asym(3_JPIM,iwn)*LIWC_tl
    asymj = baran2013_regcoef_asym(1_JPIM,jwn)    + &
         & baran2013_regcoef_asym(2_JPIM,jwn)*T + &
         & baran2013_regcoef_asym(3_JPIM,jwn)*LIWC
    asymj_tl =  &
         & baran2013_regcoef_asym(2_JPIM,jwn)*T_tl + &
         & baran2013_regcoef_asym(3_JPIM,jwn)*LIWC_tl
    if( asymi < baran2013_corcoef_asym(1_JPIM,iwn) ) then
      asymi = ( asymi - baran2013_corcoef_asym(2_JPIM,iwn) ) / baran2013_corcoef_asym(3_JPIM,iwn)
      asymi_tl = asymi_tl / baran2013_corcoef_asym(3_JPIM,iwn)
    elseif( asymi > baran2013_corcoef_asym(4_JPIM,iwn) ) then
      asymi = ( asymi - baran2013_corcoef_asym(5_JPIM,iwn) ) / baran2013_corcoef_asym(6_JPIM,iwn)
      asymi_tl = asymi_tl / baran2013_corcoef_asym(6_JPIM,iwn)
    endif
    if( asymj < baran2013_corcoef_asym(1_JPIM,jwn) ) then
      asymj = ( asymj - baran2013_corcoef_asym(2_JPIM,jwn) ) / baran2013_corcoef_asym(3_JPIM,jwn)
      asymj_tl = asymj_tl / baran2013_corcoef_asym(3_JPIM,jwn)
    elseif( asymj > baran2013_corcoef_asym(4_JPIM,jwn) ) then
      asymj = ( asymj - baran2013_corcoef_asym(5_JPIM,jwn) ) / baran2013_corcoef_asym(6_JPIM,jwn)
      asymj_tl = asymj_tl / baran2013_corcoef_asym(6_JPIM,jwn)
    endif

    asym = asymi + (asymj - asymi) * dx_dwn
    asym_tl = asymi_tl + (asymj_tl - asymi_tl) * dx_dwn

    if( asym .GT. 1.0_JPRB) then
      asym = 1.0_JPRB
      asym_tl = 0.0_JPRB
    endif

END SUBROUTINE
