! compute optical parameters for Ice cloud layer using Baran data
SUBROUTINE rttov_baran2013_calc_optpar_ad (optp, ichn, &
                       & T_in, IWC_in, T_in_ad, IWC_in_ad, &
                       & abs_ad, sca_ad, bpr_ad, asym_ad)

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
!                        for absorption and scattering coefficients (J.Vidot)
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
  REAL(KIND=JPRB)       ,INTENT(IN)  :: t_in
  REAL(KIND=JPRB)       ,INTENT(IN)  :: iwc_in
  REAL(KIND=JPRB)       ,INTENT(INOUT):: t_in_ad
  REAL(KIND=JPRB)       ,INTENT(INOUT):: iwc_in_ad

  REAL(KIND=JPRB)       ,INTENT(INOUT) :: abs_ad
  REAL(KIND=JPRB)       ,INTENT(INOUT) :: sca_ad
  REAL(KIND=JPRB)       ,INTENT(INOUT) :: bpr_ad
  REAL(KIND=JPRB)       ,INTENT(INOUT) :: asym_ad

!INTF_END

! Local Scalars:
  REAL(KIND=JPRB)  :: absi, absj, bpri, bprj
  REAL(KIND=JPRB)  :: scai, scaj, asym, asymi, asymj
  REAL(KIND=JPRB)  :: absi_ad, absj_ad, bpri_ad, asymi_ad
  REAL(KIND=JPRB)  :: scai_ad, scaj_ad, bprj_ad, asymj_ad


  REAL(KIND=JPRB)     :: dx_dwn
  REAL(KIND=JPRB)     :: liwc, liwc_ad
  INTEGER(KIND=JPIM)  :: iwn, jwn
  REAL(KIND=JPRB)     :: iwc,t
  REAL(KIND=JPRB)     :: iwc_ad,t_ad

  LOGICAL :: test1, test2, test3, test4
  LOGICAL :: test5, test6, test7, test8
  

!- End of header --------------------------------------------------------
    test5 =.FALSE.
    test6 =.FALSE.
    test7 =.FALSE.
    test8 =.FALSE.
    

    if (iwc_in < baran2013_iwc_min) then
      iwc = baran2013_iwc_min
      test5=.TRUE.
    else if (iwc_in > baran2013_iwc_max) then
      iwc = baran2013_iwc_max 
      test6=.TRUE.
    else
      iwc = iwc_in 
    endif

    if (T_in < baran2013_temp_min) then 
      T = baran2013_temp_min
      test7=.TRUE.      
    else if (T_in > baran2013_temp_max) then 
      T = baran2013_temp_max
      test8=.TRUE.
    else
      T = T_in 
    endif

    LIWC = LOG10(iwc)

    iwn = optp%optpiclb%iwn2013(ichn)  
    jwn = optp%optpiclb%jwn2013(ichn)  
    dx_dwn = optp%optpiclb%dx_dwn2013(ichn)
    
    t_ad     = 0.0_JPRB
    iwc_ad   = 0.0_JPRB
    liwc_ad  = 0.0_JPRB
    absi_ad  = 0.0_JPRB
    absj_ad  = 0.0_JPRB
    bpri_ad  = 0.0_JPRB
    bprj_ad  = 0.0_JPRB
    scai_ad  = 0.0_JPRB
    scaj_ad  = 0.0_JPRB
    asymi_ad = 0.0_JPRB
    asymj_ad = 0.0_JPRB


! ASYM
    test1 =.FALSE.
    test2 =.FALSE.
    test3 =.FALSE.
    test4 =.FALSE.

    asymi = baran2013_regcoef_asym(1_JPIM,iwn)    + &
         &  baran2013_regcoef_asym(2_JPIM,iwn)*T  + &
         &  baran2013_regcoef_asym(3_JPIM,iwn)*LIWC
    asymj = baran2013_regcoef_asym(1_JPIM,jwn)    + &
         &  baran2013_regcoef_asym(2_JPIM,jwn)*T  + &
         &  baran2013_regcoef_asym(3_JPIM,jwn)*LIWC
    if( asymi < baran2013_corcoef_asym(1_JPIM,iwn) ) then
      asymi = ( asymi - baran2013_corcoef_asym(2_JPIM,iwn) ) / baran2013_corcoef_asym(3_JPIM,iwn)
      test1 = .TRUE.
    elseif( asymi > baran2013_corcoef_asym(4_JPIM,iwn) ) then
      asymi = ( asymi - baran2013_corcoef_asym(5_JPIM,iwn) ) / baran2013_corcoef_asym(6_JPIM,iwn)
      test2 = .TRUE.
    endif
    if( asymj < baran2013_corcoef_asym(1_JPIM,jwn) ) then
      asymj = ( asymj - baran2013_corcoef_asym(2_JPIM,jwn) ) / baran2013_corcoef_asym(3_JPIM,jwn)
      test3 = .TRUE.
    elseif( asymj > baran2013_corcoef_asym(4_JPIM,jwn) ) then
      asymj = ( asymj - baran2013_corcoef_asym(5_JPIM,jwn) ) / baran2013_corcoef_asym(6_JPIM,jwn)
      test4 = .TRUE.
    endif

    asym = asymi + (asymj - asymi) * dx_dwn

    if( asym .GT. 1.0_JPRB) then
      asym_ad = 0.0_JPRB
    endif
    
    asymi_ad = asymi_ad + (1.0_JPRB - dx_dwn)*asym_ad
    asymj_ad = asymj_ad + asym_ad * dx_dwn

    if( test1 ) then
      asymi_ad = asymi_ad / baran2013_corcoef_asym(3_JPIM,iwn)
    elseif( test2) then
      asymi_ad = asymi_ad / baran2013_corcoef_asym(6_JPIM,iwn)
    endif
    if( test3 ) then
      asymj_ad = asymj_ad / baran2013_corcoef_asym(3_JPIM,jwn)
    elseif( test4 ) then
      asymj_ad = asymj_ad / baran2013_corcoef_asym(6_JPIM,jwn)
    endif

    T_ad =     T_ad    + baran2013_regcoef_asym(2_JPIM,iwn)*asymi_ad 
    LIWC_ad =  LIWC_ad + baran2013_regcoef_asym(3_JPIM,iwn)*asymi_ad
    asymi_ad = 0.0_JPRB
    T_ad =     T_ad    + baran2013_regcoef_asym(2_JPIM,jwn)*asymj_ad 
    LIWC_ad =  LIWC_ad + baran2013_regcoef_asym(3_JPIM,jwn)*asymj_ad
    asymj_ad = 0.0_JPRB



! BPR
    test1 =.FALSE.
    test2 =.FALSE.
    test3 =.FALSE.
    test4 =.FALSE.

    bpri = baran2013_regcoef_bpr(1_JPIM,iwn)    + &
         & baran2013_regcoef_bpr(2_JPIM,iwn)*T  + &
         & baran2013_regcoef_bpr(3_JPIM,iwn)*LIWC
    bprj = baran2013_regcoef_bpr(1_JPIM,jwn)    + &
         & baran2013_regcoef_bpr(2_JPIM,jwn)*T  + &
         & baran2013_regcoef_bpr(3_JPIM,jwn)*LIWC
    if( bpri < baran2013_corcoef_bpr(1_JPIM,iwn) ) then
      bpri = ( bpri - baran2013_corcoef_bpr(2_JPIM,iwn) ) / baran2013_corcoef_bpr(3_JPIM,iwn)
      test1 = .TRUE.
    elseif( bpri > baran2013_corcoef_bpr(4_JPIM,iwn) ) then
      bpri = ( bpri - baran2013_corcoef_bpr(5_JPIM,iwn) ) / baran2013_corcoef_bpr(6_JPIM,iwn)
      test2 = .TRUE.
    endif
    if( bprj < baran2013_corcoef_bpr(1_JPIM,jwn) ) then
      bprj = ( bprj - baran2013_corcoef_bpr(2_JPIM,jwn) ) / baran2013_corcoef_bpr(3_JPIM,jwn)
      test3 = .TRUE.
    elseif( bprj > baran2013_corcoef_bpr(4_JPIM,jwn) ) then
      bprj = ( bprj - baran2013_corcoef_bpr(5_JPIM,jwn) ) / baran2013_corcoef_bpr(6_JPIM,jwn)
      test4 = .TRUE.
    endif

!     bpr = bpri + (bprj - bpri) * dx_dwn

    bpri_ad = bpri_ad + (1.0_JPRB - dx_dwn)*bpr_ad
    bprj_ad = bprj_ad + bpr_ad * dx_dwn
    
    if( test1 ) then
      bpri_ad = bpri_ad / baran2013_corcoef_bpr(3_JPIM,iwn)
    elseif( test2) then
      bpri_ad = bpri_ad / baran2013_corcoef_bpr(6_JPIM,iwn)
    endif
    if( test3 ) then
      bprj_ad = bprj_ad / baran2013_corcoef_bpr(3_JPIM,jwn)
    elseif( test4 ) then
      bprj_ad = bprj_ad / baran2013_corcoef_bpr(6_JPIM,jwn)
    endif

    T_ad =     T_ad    + baran2013_regcoef_bpr(2_JPIM,iwn)*bpri_ad 
    LIWC_ad =  LIWC_ad + baran2013_regcoef_bpr(3_JPIM,iwn)*bpri_ad
    bpri_ad = 0.0_JPRB
    T_ad =     T_ad    + baran2013_regcoef_bpr(2_JPIM,jwn)*bprj_ad 
    LIWC_ad =  LIWC_ad + baran2013_regcoef_bpr(3_JPIM,jwn)*bprj_ad
    bprj_ad = 0.0_JPRB


! SCA
    scai = baran2013_regcoef_sca(1_JPIM,iwn)    + &
         & baran2013_regcoef_sca(2_JPIM,iwn)*T + &
         & baran2013_regcoef_sca(3_JPIM,iwn)*LIWC + &
         & baran2013_regcoef_sca(4_JPIM,iwn)*T*T + &
         & baran2013_regcoef_sca(5_JPIM,iwn)*LIWC*LIWC + &
         & baran2013_regcoef_sca(6_JPIM,iwn)*T*LIWC         
    scaj = baran2013_regcoef_sca(1_JPIM,jwn)    + &
         & baran2013_regcoef_sca(2_JPIM,jwn)*T + &
         & baran2013_regcoef_sca(3_JPIM,jwn)*LIWC + &
         & baran2013_regcoef_sca(4_JPIM,jwn)*T*T + &
         & baran2013_regcoef_sca(5_JPIM,jwn)*LIWC*LIWC + &
         & baran2013_regcoef_sca(6_JPIM,jwn)*T*LIWC
    scai = 10**scai
    scaj = 10**scaj

    scai_ad = scai_ad + (1.0_JPRB - dx_dwn)*sca_ad
    scaj_ad = scaj_ad + sca_ad * dx_dwn
    scai_ad = scai * scai_ad * log(10.0_JPRB)
    scaj_ad = scaj * scaj_ad * log(10.0_JPRB)

    T_ad =     T_ad    + (baran2013_regcoef_sca(2_JPIM,iwn) + 2.0_JPRB*baran2013_regcoef_sca(4_JPIM,iwn)*T + &
                          baran2013_regcoef_sca(6_JPIM,iwn)*LIWC)*scai_ad
    LIWC_ad =  LIWC_ad + (baran2013_regcoef_sca(3_JPIM,iwn) + 2.0_JPRB*baran2013_regcoef_sca(5_JPIM,iwn)*LIWC + &
                          baran2013_regcoef_sca(6_JPIM,iwn)*T)*scai_ad
    scai_ad = 0.0_JPRB
    T_ad =     T_ad    + (baran2013_regcoef_sca(2_JPIM,jwn) + 2.0_JPRB*baran2013_regcoef_sca(4_JPIM,jwn)*T + &
                          baran2013_regcoef_sca(6_JPIM,jwn)*LIWC)*scaj_ad
    LIWC_ad =  LIWC_ad + (baran2013_regcoef_sca(3_JPIM,jwn) + 2.0_JPRB*baran2013_regcoef_sca(5_JPIM,jwn)*LIWC + &
                          baran2013_regcoef_sca(6_JPIM,jwn)*T)*scaj_ad
    scaj_ad = 0.0_JPRB


! ABS
    absi = baran2013_regcoef_abs(1_JPIM,iwn)    + &
         & baran2013_regcoef_abs(2_JPIM,iwn)*T + &
         & baran2013_regcoef_abs(3_JPIM,iwn)*LIWC + &
         & baran2013_regcoef_abs(4_JPIM,iwn)*T*T + &
         & baran2013_regcoef_abs(5_JPIM,iwn)*LIWC*LIWC + &
         & baran2013_regcoef_abs(6_JPIM,iwn)*T*LIWC         
    absj = baran2013_regcoef_abs(1_JPIM,jwn)    + &
         & baran2013_regcoef_abs(2_JPIM,jwn)*T + &
         & baran2013_regcoef_abs(3_JPIM,jwn)*LIWC + &
         & baran2013_regcoef_abs(4_JPIM,jwn)*T*T + &
         & baran2013_regcoef_abs(5_JPIM,jwn)*LIWC*LIWC + &
         & baran2013_regcoef_abs(6_JPIM,jwn)*T*LIWC
    absi = 10**absi
    absj = 10**absj


    absi_ad = absi_ad + (1.0_JPRB - dx_dwn)*abs_ad
    absj_ad = absj_ad + abs_ad * dx_dwn
    absi_ad = absi * absi_ad * log(10.0_JPRB)
    absj_ad = absj * absj_ad * log(10.0_JPRB)

    T_ad =     T_ad    + (baran2013_regcoef_abs(2_JPIM,iwn) + 2.0_JPRB*baran2013_regcoef_abs(4_JPIM,iwn)*T + &
                          baran2013_regcoef_abs(6_JPIM,iwn)*LIWC)*absi_ad
    LIWC_ad =  LIWC_ad + (baran2013_regcoef_abs(3_JPIM,iwn) + 2.0_JPRB*baran2013_regcoef_abs(5_JPIM,iwn)*LIWC + &
                          baran2013_regcoef_abs(6_JPIM,iwn)*T)*absi_ad
    absi_ad = 0.0_JPRB
    T_ad =     T_ad    + (baran2013_regcoef_abs(2_JPIM,jwn) + 2.0_JPRB*baran2013_regcoef_abs(4_JPIM,jwn)*T + &
                          baran2013_regcoef_abs(6_JPIM,jwn)*LIWC)*absj_ad
    LIWC_ad =  LIWC_ad + (baran2013_regcoef_abs(3_JPIM,jwn) + 2.0_JPRB*baran2013_regcoef_abs(5_JPIM,jwn)*LIWC + &
                          baran2013_regcoef_abs(6_JPIM,jwn)*T)*absj_ad
    absj_ad = 0.0_JPRB


    iwc_ad = iwc_ad + liwc_ad / ( iwc * log(10.0_JPRB) )
    liwc_ad = 0.0_JPRB

    if (.not. (test7 .or. test8)) then
      T_in_ad = T_in_ad + T_ad 
    endif
    if (.not. (test5 .or. test6)) then
      iwc_in_ad = iwc_in_ad + iwc_ad
    endif
    iwc_ad=0.0_JPRB 
    T_ad=0.0_JPRB 

END SUBROUTINE
