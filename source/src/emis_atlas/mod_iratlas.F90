!
MODULE mod_iratlas
  ! Description:
  !   Data and routines for IR emissivity atlas.
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
  ! Current Code Owner: SAF NWP
  !
  ! History:
  ! Version   Date     Comment
  ! -------   ----     -------
  !  1.0      02/06/2010  Based on UW IR atlas code (E. Borbas, B. Ruston, J. Hocking)
  !  1.1      02/06/2014  satellite zenith angle correction added (E. Borbas)
  !
  ! Code Description:
  !   Language:           Fortran 90.
  !   Software Standards: "European Standards for Writing and
  !     Documenting Exchangeable Fortran 90 Code".

#include "throw.h"

  USE parkind1, ONLY : &
    jpim, &   ! 32-bit int
    jpis, &   ! 16-bit
    jpit, &   ! 8-bit
    jprb, &   ! 64-bit real
    jprm, &   ! 32-bit
    jplm      ! logical

  USE mod_rttov_emis_atlas, ONLY : &
    ir_atlas_std_init,    &
    ir_atlas_single_inst, &
    ir_atlas_do_ang_corr

#ifdef _RTTOV_HDF
  USE hdf5
  USE rttov_hdf_mod, ONLY : &
    open_hdf,             &
    close_hdf,            &
    is_hdf_open,          &
    is_hdf_64bit_reals
#endif

  IMPLICIT NONE

#include "rttov_errorreport.interface"

#ifndef _RTTOV_HDF
  INCLUDE 'netcdf.inc'
#endif

  ! User can specify atlas version at set-up. Version 100 is the default.
  INTEGER(KIND=jpim), PARAMETER :: ir_atlas_version = 100

  ! Atlas constants

  INTEGER(KIND=jpim), PARAMETER :: numpcs = 6
  INTEGER(KIND=jpim), PARAMETER :: hngpnts = 10
  INTEGER(KIND=jpim), PARAMETER :: numwave = 416

!   INTEGER(KIND=jpim), PARAMETER :: nb_lats = 1800
!   INTEGER(KIND=jpim), PARAMETER :: nb_lons = 3600
!   ! INTEGER(KIND=jpim), PARAMETER :: nb_pack=2298394  !ver2.1  2006
!   INTEGER(KIND=jpim), PARAMETER :: nb_pack = 2250931  !ver2.1  2007
  INTEGER(KIND=jpim) :: nb_lats
  INTEGER(KIND=jpim) :: nb_lons
  INTEGER(KIND=jpim) :: nb_pack

!   INTEGER(KIND=jpim), PARAMETER :: cv_lats = 360
!   INTEGER(KIND=jpim), PARAMETER :: cv_lons = 720
!   INTEGER(KIND=jpim), PARAMETER :: cv_pack = 98008
  INTEGER(KIND=jpim) :: cv_lats
  INTEGER(KIND=jpim) :: cv_lons
  INTEGER(KIND=jpim) :: cv_pack

!   INTEGER(KIND=jpim), PARAMETER :: nb_05lats = 3600
!   INTEGER(KIND=jpim), PARAMETER :: nb_05lons = 7200
!   INTEGER(KIND=jpim), PARAMETER :: nb_igbp = 18
  INTEGER(KIND=jpim) :: igbp_lats
  INTEGER(KIND=jpim) :: igbp_lons
  INTEGER(KIND=jpim) :: nb_igbp

  INTEGER(KIND=jpim), PARAMETER :: db_ver_year = 2007

  INTEGER(KIND=jpim), PARAMETER :: seaice_flag = 70           ! flag value returned for sea-ice
  REAL(KIND=jprb),    PARAMETER :: default_std = 0.05_jprb    ! default standard deviation

  INTEGER(KIND=jpim), PARAMETER :: bfemis_gridres = 100       ! 0.1 deg
  INTEGER(KIND=jpim), PARAMETER :: bfemis_ygrid1 = 89950      ! 89.95 deg
  INTEGER(KIND=jpim), PARAMETER :: bfemis_xgrid1 = -179950    ! -179.95 deg

  INTEGER(KIND=jpim), PARAMETER :: cov_emis_gridres = 500     ! 0.5 deg
  INTEGER(KIND=jpim), PARAMETER :: cov_emis_ygrid1 = 89750    ! 89.75 deg
  INTEGER(KIND=jpim), PARAMETER :: cov_emis_xgrid1 = -179750  ! -179.75 deg

  INTEGER(KIND=jpim), PARAMETER :: igbp_gridres = 50          ! 0.05 deg
  INTEGER(KIND=jpim), PARAMETER :: igbp_ygrid1 = 89950        ! 89.95 deg
  INTEGER(KIND=jpim), PARAMETER :: igbp_xgrid1 = -179950      ! -179.95 deg

  REAL(KIND=jprb),    PARAMETER :: angcorrminzen = 5._jprb
  REAL(KIND=jprb),    PARAMETER :: angcorrterm = 85._jprb     ! Solar zenith angle of terminator

  ! Atlas data loaded by initialisation routine

  INTEGER(KIND=jpis), POINTER :: bfemis_flag(:,:)  ! dims are (nb_lats,nb_lons)
  INTEGER(KIND=jpim), POINTER :: bfemis_lut(:,:)   ! dims are (nb_lats,nb_lons)
  REAL(KIND=jprm),    POINTER :: pca_coef(:,:)     ! dims are (nb_pack,numpcs)

  INTEGER(KIND=jpim), POINTER :: cov_emis_lut(:,:) ! dims are (cv_lats,cv_lons)
  INTEGER(KIND=jpis), POINTER :: cov_emis(:,:)     ! dims are (cv_pack,numwave)

  INTEGER(KIND=jpit), POINTER :: igbp(:,:)         ! dims are (igbp_lats,igbp_lons)

  REAL(KIND=jprm),    POINTER :: p1d(:,:)          ! dims are (numwave,nb_igbp)
  REAL(KIND=jprm),    POINTER :: p2d(:,:)          ! dims are (numwave,nb_igbp)
  REAL(KIND=jprm),    POINTER :: p3d(:,:)          ! dims are (numwave,nb_igbp)
  REAL(KIND=jprm),    POINTER :: p1n(:,:)          ! dims are (numwave,nb_igbp)
  REAL(KIND=jprm),    POINTER :: p2n(:,:)          ! dims are (numwave,nb_igbp)
  REAL(KIND=jprm),    POINTER :: p3n(:,:)          ! dims are (numwave,nb_igbp)

  REAL(KIND=jprm) :: pcu(numpcs,numwave)           ! pcu need to be read in from flat file

  ! Data to allow more efficient memory usage (convert jprm to jprb only at point of use)
#ifndef _RTTOV_HDF
  REAL(KIND=jprm)             :: pca_sfac(numpcs)  ! PCA coef scale factors
  REAL(KIND=jprm)             :: pca_offs(numpcs)  ! PCA coef offsets
#else
  REAL(KIND=jprm),    POINTER :: pca_sfac(:)       ! PCA coef scale factors
  REAL(KIND=jprm),    POINTER :: pca_offs(:)       ! PCA coef offsets
#endif
  REAL(KIND=jprm) :: cov_sfac                      ! Stdev scale factor

  ! Atlas data contained in this file

  REAL(KIND=jprb) :: sice_em(numwave), snow_em(numwave)
  REAL(KIND=jprb) :: sice_stdv, snow_stdv
  REAL(KIND=jprb) :: pcm(numwave)
  REAL(KIND=jprb) :: hsr_wavenum(numwave)

  ! Arrays to hold hsr data interpolated onto channel wavenumbers

  INTEGER(KIND=jpim) :: ncoefchans                    ! Number of channels in coef file
  REAL(KIND=jprb), POINTER :: cov_emis_int(:,:)       ! dims are (cv_pack,nchannels)
  REAL(KIND=jprb), POINTER :: pcu_int(:,:,:)          ! dims are (numpcs,nchannels,1/2)
  REAL(KIND=jprb), POINTER :: pcm_int(:,:)
  REAL(KIND=jprb), POINTER :: sice_em_int(:), snow_em_int(:)
  REAL(KIND=jprm), POINTER :: p1d_int(:,:,:)
  REAL(KIND=jprm), POINTER :: p2d_int(:,:,:)
  REAL(KIND=jprm), POINTER :: p3d_int(:,:,:)
  REAL(KIND=jprm), POINTER :: p1n_int(:,:,:)
  REAL(KIND=jprm), POINTER :: p2n_int(:,:,:)
  REAL(KIND=jprm), POINTER :: p3n_int(:,:,:)

  DATA hsr_wavenum / &
    699.3_jprb,  704.3_jprb,  709.3_jprb,  714.3_jprb,  719.3_jprb,  724.3_jprb,  729.3_jprb,  734.3_jprb,  &
    739.3_jprb,  744.3_jprb,  749.3_jprb,  754.3_jprb,  759.3_jprb,  764.3_jprb,  769.3_jprb,  774.3_jprb,  &
    779.3_jprb,  784.3_jprb,  789.3_jprb,  794.3_jprb,  799.3_jprb,  804.3_jprb,  809.3_jprb,  814.3_jprb,  &
    819.3_jprb,  824.3_jprb,  829.3_jprb,  834.3_jprb,  839.3_jprb,  844.3_jprb,  849.3_jprb,  854.3_jprb,  &
    859.3_jprb,  864.3_jprb,  869.3_jprb,  874.3_jprb,  879.3_jprb,  884.3_jprb,  889.3_jprb,  894.3_jprb,  &
    899.3_jprb,  904.3_jprb,  909.3_jprb,  914.3_jprb,  919.3_jprb,  924.3_jprb,  929.3_jprb,  934.3_jprb,  &
    939.3_jprb,  944.3_jprb,  949.3_jprb,  954.3_jprb,  959.3_jprb,  964.3_jprb,  969.3_jprb,  974.3_jprb,  &
    979.3_jprb,  984.3_jprb,  989.3_jprb,  994.3_jprb,  999.3_jprb, 1004.3_jprb, 1009.3_jprb, 1014.3_jprb,  &
    1019.3_jprb, 1024.3_jprb, 1029.3_jprb, 1034.3_jprb, 1039.3_jprb, 1044.3_jprb, 1049.3_jprb, 1054.3_jprb,  &
    1059.3_jprb, 1064.3_jprb, 1069.3_jprb, 1074.3_jprb, 1079.3_jprb, 1084.3_jprb, 1089.3_jprb, 1094.3_jprb,  &
    1099.3_jprb, 1104.3_jprb, 1109.3_jprb, 1114.3_jprb, 1119.3_jprb, 1124.3_jprb, 1129.3_jprb, 1134.3_jprb,  &
    1139.3_jprb, 1144.3_jprb, 1149.3_jprb, 1154.3_jprb, 1159.3_jprb, 1164.3_jprb, 1169.3_jprb, 1174.3_jprb,  &
    1179.3_jprb, 1184.3_jprb, 1189.3_jprb, 1194.3_jprb, 1199.3_jprb, 1204.3_jprb, 1209.3_jprb, 1214.3_jprb,  &
    1219.3_jprb, 1224.3_jprb, 1229.3_jprb, 1234.3_jprb, 1239.3_jprb, 1244.3_jprb, 1249.3_jprb, 1254.3_jprb,  &
    1259.3_jprb, 1264.3_jprb, 1269.3_jprb, 1274.3_jprb, 1279.3_jprb, 1284.3_jprb, 1289.3_jprb, 1294.3_jprb,  &
    1299.3_jprb, 1304.3_jprb, 1309.3_jprb, 1314.3_jprb, 1319.3_jprb, 1324.3_jprb, 1329.3_jprb, 1334.3_jprb,  &
    1339.3_jprb, 1344.3_jprb, 1349.3_jprb, 1354.3_jprb, 1359.3_jprb, 1364.3_jprb, 1369.3_jprb, 1374.3_jprb,  &
    1379.3_jprb, 1384.3_jprb, 1389.3_jprb, 1394.3_jprb, 1399.3_jprb, 1404.3_jprb, 1409.3_jprb, 1414.3_jprb,  &
    1419.3_jprb, 1424.3_jprb, 1429.3_jprb, 1434.3_jprb, 1439.3_jprb, 1444.3_jprb, 1449.3_jprb, 1454.3_jprb,  &
    1459.3_jprb, 1464.3_jprb, 1469.3_jprb, 1474.3_jprb, 1479.3_jprb, 1484.3_jprb, 1489.3_jprb, 1494.3_jprb,  &
    1499.3_jprb, 1504.3_jprb, 1509.3_jprb, 1514.3_jprb, 1519.3_jprb, 1524.3_jprb, 1529.3_jprb, 1534.3_jprb,  &
    1539.3_jprb, 1544.3_jprb, 1549.3_jprb, 1554.3_jprb, 1559.3_jprb, 1564.3_jprb, 1569.3_jprb, 1574.3_jprb,  &
    1579.3_jprb, 1584.3_jprb, 1589.3_jprb, 1594.3_jprb, 1599.3_jprb, 1604.3_jprb, 1609.3_jprb, 1614.3_jprb,  &
    1619.3_jprb, 1624.3_jprb, 1629.3_jprb, 1634.3_jprb, 1639.3_jprb, 1644.3_jprb, 1649.3_jprb, 1654.3_jprb,  &
    1659.3_jprb, 1664.3_jprb, 1669.3_jprb, 1674.3_jprb, 1679.3_jprb, 1684.3_jprb, 1689.3_jprb, 1694.3_jprb,  &
    1699.3_jprb, 1704.3_jprb, 1709.3_jprb, 1714.3_jprb, 1719.3_jprb, 1724.3_jprb, 1729.3_jprb, 1734.3_jprb,  &
    1739.3_jprb, 1744.3_jprb, 1749.3_jprb, 1754.3_jprb, 1759.3_jprb, 1764.3_jprb, 1769.3_jprb, 1774.3_jprb,  &
    1779.3_jprb, 1784.3_jprb, 1789.3_jprb, 1794.3_jprb, 1799.3_jprb, 1804.3_jprb, 1809.3_jprb, 1814.3_jprb,  &
    1819.3_jprb, 1824.3_jprb, 1829.3_jprb, 1834.3_jprb, 1839.3_jprb, 1844.3_jprb, 1849.3_jprb, 1854.3_jprb,  &
    1859.3_jprb, 1864.3_jprb, 1869.3_jprb, 1874.3_jprb, 1879.3_jprb, 1884.3_jprb, 1889.3_jprb, 1894.3_jprb,  &
    1899.3_jprb, 1904.3_jprb, 1909.3_jprb, 1914.3_jprb, 1919.3_jprb, 1924.3_jprb, 1929.3_jprb, 1934.3_jprb,  &
    1939.3_jprb, 1944.3_jprb, 1949.3_jprb, 1954.3_jprb, 1959.3_jprb, 1964.3_jprb, 1969.3_jprb, 1974.3_jprb,  &
    1979.3_jprb, 1984.3_jprb, 1989.3_jprb, 1994.3_jprb, 1999.3_jprb, 2004.3_jprb, 2009.3_jprb, 2014.3_jprb,  &
    2019.3_jprb, 2024.3_jprb, 2029.3_jprb, 2034.3_jprb, 2039.3_jprb, 2044.3_jprb, 2049.3_jprb, 2054.3_jprb,  &
    2059.3_jprb, 2064.3_jprb, 2069.3_jprb, 2074.3_jprb, 2079.3_jprb, 2084.3_jprb, 2089.3_jprb, 2094.3_jprb,  &
    2099.3_jprb, 2104.3_jprb, 2109.3_jprb, 2114.3_jprb, 2119.3_jprb, 2124.3_jprb, 2129.3_jprb, 2134.3_jprb,  &
    2139.3_jprb, 2144.3_jprb, 2149.3_jprb, 2154.3_jprb, 2159.3_jprb, 2164.3_jprb, 2169.3_jprb, 2174.3_jprb,  &
    2179.3_jprb, 2184.3_jprb, 2189.3_jprb, 2194.3_jprb, 2199.3_jprb, 2204.3_jprb, 2209.3_jprb, 2214.3_jprb,  &
    2219.3_jprb, 2224.3_jprb, 2229.3_jprb, 2234.3_jprb, 2239.3_jprb, 2244.3_jprb, 2249.3_jprb, 2254.3_jprb,  &
    2259.3_jprb, 2264.3_jprb, 2269.3_jprb, 2274.3_jprb, 2279.3_jprb, 2284.3_jprb, 2289.3_jprb, 2294.3_jprb,  &
    2299.3_jprb, 2304.3_jprb, 2309.3_jprb, 2314.3_jprb, 2319.3_jprb, 2324.3_jprb, 2329.3_jprb, 2334.3_jprb,  &
    2339.3_jprb, 2344.3_jprb, 2349.3_jprb, 2354.3_jprb, 2359.3_jprb, 2364.3_jprb, 2369.3_jprb, 2374.3_jprb,  &
    2379.3_jprb, 2384.3_jprb, 2389.3_jprb, 2394.3_jprb, 2399.3_jprb, 2404.3_jprb, 2409.3_jprb, 2414.3_jprb,  &
    2419.3_jprb, 2424.3_jprb, 2429.3_jprb, 2434.3_jprb, 2439.3_jprb, 2444.3_jprb, 2449.3_jprb, 2454.3_jprb,  &
    2459.3_jprb, 2464.3_jprb, 2469.3_jprb, 2474.3_jprb, 2479.3_jprb, 2484.3_jprb, 2489.3_jprb, 2494.3_jprb,  &
    2499.3_jprb, 2504.3_jprb, 2509.3_jprb, 2514.3_jprb, 2519.3_jprb, 2524.3_jprb, 2529.3_jprb, 2534.3_jprb,  &
    2539.3_jprb, 2544.3_jprb, 2549.3_jprb, 2554.3_jprb, 2559.3_jprb, 2564.3_jprb, 2569.3_jprb, 2574.3_jprb,  &
    2579.3_jprb, 2584.3_jprb, 2589.3_jprb, 2594.3_jprb, 2599.3_jprb, 2604.3_jprb, 2609.3_jprb, 2614.3_jprb,  &
    2619.3_jprb, 2624.3_jprb, 2629.3_jprb, 2634.3_jprb, 2639.3_jprb, 2644.3_jprb, 2649.3_jprb, 2654.3_jprb,  &
    2659.3_jprb, 2664.3_jprb, 2669.3_jprb, 2674.3_jprb, 2679.3_jprb, 2684.3_jprb, 2689.3_jprb, 2694.3_jprb,  &
    2699.3_jprb, 2704.3_jprb, 2709.3_jprb, 2714.3_jprb, 2719.3_jprb, 2724.3_jprb, 2729.3_jprb, 2734.3_jprb,  &
    2739.3_jprb, 2744.3_jprb, 2749.3_jprb, 2754.3_jprb, 2759.3_jprb, 2764.3_jprb, 2769.3_jprb, 2774.3_jprb /

  DATA pcm / &
        0.9782182_jprb,   0.9770744_jprb,   0.9763290_jprb,   0.9763215_jprb,   0.9760258_jprb,  &
        0.9763704_jprb,   0.9767076_jprb,   0.9763077_jprb,   0.9758835_jprb,   0.9753462_jprb,  &
        0.9748067_jprb,   0.9734465_jprb,   0.9721510_jprb,   0.9717180_jprb,   0.9714773_jprb,  &
        0.9706340_jprb,   0.9710826_jprb,   0.9722888_jprb,   0.9731166_jprb,   0.9732918_jprb,  &
        0.9736975_jprb,   0.9751787_jprb,   0.9770049_jprb,   0.9773170_jprb,   0.9765164_jprb,  &
        0.9759824_jprb,   0.9750199_jprb,   0.9746831_jprb,   0.9738413_jprb,   0.9731615_jprb,  &
        0.9720387_jprb,   0.9716908_jprb,   0.9708628_jprb,   0.9705366_jprb,   0.9697853_jprb,  &
        0.9694459_jprb,   0.9688896_jprb,   0.9688236_jprb,   0.9689180_jprb,   0.9692774_jprb,  &
        0.9693237_jprb,   0.9692513_jprb,   0.9689918_jprb,   0.9686664_jprb,   0.9684489_jprb,  &
        0.9681804_jprb,   0.9672847_jprb,   0.9667084_jprb,   0.9661347_jprb,   0.9655386_jprb,  &
        0.9650131_jprb,   0.9641176_jprb,   0.9628995_jprb,   0.9620982_jprb,   0.9605948_jprb,  &
        0.9590283_jprb,   0.9572537_jprb,   0.9552648_jprb,   0.9529146_jprb,   0.9505763_jprb,  &
        0.9486620_jprb,   0.9468448_jprb,   0.9446425_jprb,   0.9428397_jprb,   0.9415421_jprb,  &
        0.9398234_jprb,   0.9378662_jprb,   0.9358756_jprb,   0.9338515_jprb,   0.9317511_jprb,  &
        0.9296144_jprb,   0.9274116_jprb,   0.9248639_jprb,   0.9219664_jprb,   0.9197029_jprb,  &
        0.9187206_jprb,   0.9195539_jprb,   0.9211251_jprb,   0.9227578_jprb,   0.9242273_jprb,  &
        0.9256495_jprb,   0.9265392_jprb,   0.9276078_jprb,   0.9279289_jprb,   0.9282181_jprb,  &
        0.9284544_jprb,   0.9289097_jprb,   0.9299400_jprb,   0.9314128_jprb,   0.9329405_jprb,  &
        0.9349486_jprb,   0.9377099_jprb,   0.9380918_jprb,   0.9354525_jprb,   0.9330018_jprb,  &
        0.9316696_jprb,   0.9308965_jprb,   0.9296793_jprb,   0.9282659_jprb,   0.9273711_jprb,  &
        0.9268156_jprb,   0.9265846_jprb,   0.9264724_jprb,   0.9278417_jprb,   0.9298262_jprb,  &
        0.9342009_jprb,   0.9397170_jprb,   0.9451398_jprb,   0.9501663_jprb,   0.9547508_jprb,  &
        0.9586911_jprb,   0.9618842_jprb,   0.9649577_jprb,   0.9675525_jprb,   0.9696881_jprb,  &
        0.9708689_jprb,   0.9717879_jprb,   0.9722518_jprb,   0.9724457_jprb,   0.9728941_jprb,  &
        0.9731293_jprb,   0.9731925_jprb,   0.9730867_jprb,   0.9733831_jprb,   0.9735166_jprb,  &
        0.9740434_jprb,   0.9742066_jprb,   0.9746855_jprb,   0.9748268_jprb,   0.9749292_jprb,  &
        0.9751188_jprb,   0.9752902_jprb,   0.9751062_jprb,   0.9751985_jprb,   0.9752622_jprb,  &
        0.9750626_jprb,   0.9755121_jprb,   0.9755228_jprb,   0.9760818_jprb,   0.9759580_jprb,  &
        0.9758280_jprb,   0.9755163_jprb,   0.9754220_jprb,   0.9750829_jprb,   0.9743836_jprb,  &
        0.9745844_jprb,   0.9742978_jprb,   0.9740397_jprb,   0.9744191_jprb,   0.9745796_jprb,  &
        0.9749123_jprb,   0.9750853_jprb,   0.9746974_jprb,   0.9747824_jprb,   0.9746920_jprb,  &
        0.9735873_jprb,   0.9733123_jprb,   0.9725510_jprb,   0.9718717_jprb,   0.9713586_jprb,  &
        0.9706160_jprb,   0.9701124_jprb,   0.9698699_jprb,   0.9698430_jprb,   0.9694992_jprb,  &
        0.9691019_jprb,   0.9690002_jprb,   0.9678345_jprb,   0.9668854_jprb,   0.9659764_jprb,  &
        0.9666998_jprb,   0.9669611_jprb,   0.9665817_jprb,   0.9679645_jprb,   0.9695909_jprb,  &
        0.9711555_jprb,   0.9724632_jprb,   0.9737635_jprb,   0.9746142_jprb,   0.9748497_jprb,  &
        0.9752109_jprb,   0.9752749_jprb,   0.9754022_jprb,   0.9753313_jprb,   0.9746057_jprb,  &
        0.9745884_jprb,   0.9747860_jprb,   0.9752877_jprb,   0.9753085_jprb,   0.9759305_jprb,  &
        0.9752344_jprb,   0.9748027_jprb,   0.9757417_jprb,   0.9751943_jprb,   0.9748128_jprb,  &
        0.9743713_jprb,   0.9741939_jprb,   0.9725359_jprb,   0.9723988_jprb,   0.9716700_jprb,  &
        0.9708291_jprb,   0.9705051_jprb,   0.9699901_jprb,   0.9689955_jprb,   0.9683419_jprb,  &
        0.9684200_jprb,   0.9672046_jprb,   0.9660766_jprb,   0.9658424_jprb,   0.9648336_jprb,  &
        0.9640325_jprb,   0.9642861_jprb,   0.9636880_jprb,   0.9638920_jprb,   0.9638573_jprb,  &
        0.9641714_jprb,   0.9648057_jprb,   0.9648220_jprb,   0.9639065_jprb,   0.9635883_jprb,  &
        0.9626419_jprb,   0.9616417_jprb,   0.9600965_jprb,   0.9587714_jprb,   0.9576451_jprb,  &
        0.9557189_jprb,   0.9545730_jprb,   0.9550443_jprb,   0.9551759_jprb,   0.9560625_jprb,  &
        0.9576327_jprb,   0.9587138_jprb,   0.9594474_jprb,   0.9598546_jprb,   0.9601094_jprb,  &
        0.9601356_jprb,   0.9597549_jprb,   0.9590299_jprb,   0.9581512_jprb,   0.9572046_jprb,  &
        0.9557602_jprb,   0.9538486_jprb,   0.9521495_jprb,   0.9503905_jprb,   0.9491790_jprb,  &
        0.9485527_jprb,   0.9479896_jprb,   0.9475234_jprb,   0.9468080_jprb,   0.9469628_jprb,  &
        0.9469683_jprb,   0.9465806_jprb,   0.9468755_jprb,   0.9466828_jprb,   0.9471480_jprb,  &
        0.9470276_jprb,   0.9470209_jprb,   0.9468378_jprb,   0.9464890_jprb,   0.9462101_jprb,  &
        0.9459322_jprb,   0.9449111_jprb,   0.9435923_jprb,   0.9416961_jprb,   0.9401403_jprb,  &
        0.9387150_jprb,   0.9374595_jprb,   0.9347988_jprb,   0.9319339_jprb,   0.9295776_jprb,  &
        0.9268476_jprb,   0.9243815_jprb,   0.9224647_jprb,   0.9208075_jprb,   0.9195780_jprb,  &
        0.9183103_jprb,   0.9171674_jprb,   0.9164810_jprb,   0.9160877_jprb,   0.9151877_jprb,  &
        0.9148492_jprb,   0.9142842_jprb,   0.9142084_jprb,   0.9138089_jprb,   0.9137760_jprb,  &
        0.9137531_jprb,   0.9141592_jprb,   0.9136598_jprb,   0.9125727_jprb,   0.9108481_jprb,  &
        0.9093652_jprb,   0.9080561_jprb,   0.9062355_jprb,   0.9046820_jprb,   0.9028210_jprb,  &
        0.9018152_jprb,   0.9008504_jprb,   0.9000632_jprb,   0.8995758_jprb,   0.8989593_jprb,  &
        0.8987811_jprb,   0.8992507_jprb,   0.8999549_jprb,   0.9013391_jprb,   0.9020863_jprb,  &
        0.9025120_jprb,   0.9023982_jprb,   0.9015658_jprb,   0.9008633_jprb,   0.8996401_jprb,  &
        0.8981582_jprb,   0.8969440_jprb,   0.8946483_jprb,   0.8925536_jprb,   0.8906261_jprb,  &
        0.8889833_jprb,   0.8870751_jprb,   0.8845615_jprb,   0.8825631_jprb,   0.8811586_jprb,  &
        0.8796447_jprb,   0.8779839_jprb,   0.8765292_jprb,   0.8754975_jprb,   0.8739760_jprb,  &
        0.8725729_jprb,   0.8714029_jprb,   0.8706908_jprb,   0.8710466_jprb,   0.8699325_jprb,  &
        0.8697992_jprb,   0.8718969_jprb,   0.8713725_jprb,   0.8701416_jprb,   0.8695096_jprb,  &
        0.8698574_jprb,   0.8700698_jprb,   0.8694080_jprb,   0.8693934_jprb,   0.8693246_jprb,  &
        0.8698239_jprb,   0.8696592_jprb,   0.8681608_jprb,   0.8656288_jprb,   0.8654716_jprb,  &
        0.8640761_jprb,   0.8639477_jprb,   0.8635154_jprb,   0.8630069_jprb,   0.8623275_jprb,  &
        0.8623751_jprb,   0.8627441_jprb,   0.8630516_jprb,   0.8638958_jprb,   0.8644919_jprb,  &
        0.8655882_jprb,   0.8666160_jprb,   0.8676174_jprb,   0.8692035_jprb,   0.8695340_jprb,  &
        0.8703975_jprb,   0.8714244_jprb,   0.8715467_jprb,   0.8713564_jprb,   0.8712272_jprb,  &
        0.8714187_jprb,   0.8701625_jprb,   0.8697796_jprb,   0.8688766_jprb,   0.8682391_jprb,  &
        0.8680181_jprb,   0.8676605_jprb,   0.8672657_jprb,   0.8679592_jprb,   0.8675538_jprb,  &
        0.8686572_jprb,   0.8682060_jprb,   0.8688578_jprb,   0.8693632_jprb,   0.8689557_jprb,  &
        0.8681611_jprb,   0.8684876_jprb,   0.8680010_jprb,   0.8675498_jprb,   0.8675414_jprb,  &
        0.8677824_jprb,   0.8665875_jprb,   0.8668503_jprb,   0.8665696_jprb,   0.8671130_jprb,  &
        0.8669835_jprb,   0.8671956_jprb,   0.8683699_jprb,   0.8685648_jprb,   0.8682314_jprb,  &
        0.8683055_jprb,   0.8694246_jprb,   0.8689486_jprb,   0.8693868_jprb,   0.8694460_jprb,  &
        0.8701811_jprb,   0.8704424_jprb,   0.8709887_jprb,   0.8712862_jprb,   0.8721344_jprb,  &
        0.8724745_jprb,   0.8727338_jprb,   0.8740577_jprb,   0.8748575_jprb,   0.8747587_jprb,  &
        0.8762293_jprb,   0.8772818_jprb,   0.8779803_jprb,   0.8791369_jprb,   0.8807610_jprb,  &
        0.8813813_jprb/

  DATA sice_stdv / 0.015_jprb /
  DATA sice_em / &
      0.9370_jprb, 0.9370_jprb, 0.9370_jprb, 0.9370_jprb, 0.9367_jprb, 0.9367_jprb, 0.9366_jprb, 0.9365_jprb, &
      0.9365_jprb, 0.9365_jprb, 0.9365_jprb, 0.9366_jprb, 0.9367_jprb, 0.9370_jprb, 0.9374_jprb, 0.9381_jprb, &
      0.9386_jprb, 0.9393_jprb, 0.9401_jprb, 0.9408_jprb, 0.9415_jprb, 0.9427_jprb, 0.9440_jprb, 0.9452_jprb, &
      0.9464_jprb, 0.9481_jprb, 0.9496_jprb, 0.9511_jprb, 0.9525_jprb, 0.9544_jprb, 0.9563_jprb, 0.9582_jprb, &
      0.9602_jprb, 0.9620_jprb, 0.9640_jprb, 0.9658_jprb, 0.9678_jprb, 0.9702_jprb, 0.9725_jprb, 0.9748_jprb, &
      0.9770_jprb, 0.9792_jprb, 0.9814_jprb, 0.9836_jprb, 0.9856_jprb, 0.9872_jprb, 0.9885_jprb, 0.9897_jprb, &
      0.9905_jprb, 0.9911_jprb, 0.9913_jprb, 0.9913_jprb, 0.9912_jprb, 0.9910_jprb, 0.9907_jprb, 0.9904_jprb, &
      0.9901_jprb, 0.9897_jprb, 0.9893_jprb, 0.9889_jprb, 0.9885_jprb, 0.9880_jprb, 0.9876_jprb, 0.9871_jprb, &
      0.9867_jprb, 0.9864_jprb, 0.9861_jprb, 0.9858_jprb, 0.9854_jprb, 0.9852_jprb, 0.9849_jprb, 0.9846_jprb, &
      0.9844_jprb, 0.9842_jprb, 0.9840_jprb, 0.9838_jprb, 0.9836_jprb, 0.9834_jprb, 0.9832_jprb, 0.9831_jprb, &
      0.9829_jprb, 0.9828_jprb, 0.9826_jprb, 0.9824_jprb, 0.9822_jprb, 0.9821_jprb, 0.9820_jprb, 0.9819_jprb, &
      0.9817_jprb, 0.9816_jprb, 0.9814_jprb, 0.9813_jprb, 0.9811_jprb, 0.9810_jprb, 0.9808_jprb, 0.9807_jprb, &
      0.9805_jprb, 0.9804_jprb, 0.9803_jprb, 0.9801_jprb, 0.9799_jprb, 0.9797_jprb, 0.9796_jprb, 0.9794_jprb, &
      0.9792_jprb, 0.9791_jprb, 0.9789_jprb, 0.9787_jprb, 0.9786_jprb, 0.9785_jprb, 0.9784_jprb, 0.9784_jprb, &
      0.9783_jprb, 0.9782_jprb, 0.9782_jprb, 0.9782_jprb, 0.9781_jprb, 0.9781_jprb, 0.9781_jprb, 0.9781_jprb, &
      0.9781_jprb, 0.9781_jprb, 0.9781_jprb, 0.9780_jprb, 0.9780_jprb, 0.9780_jprb, 0.9780_jprb, 0.9780_jprb, &
      0.9780_jprb, 0.9780_jprb, 0.9779_jprb, 0.9779_jprb, 0.9778_jprb, 0.9778_jprb, 0.9777_jprb, 0.9777_jprb, &
      0.9777_jprb, 0.9776_jprb, 0.9776_jprb, 0.9776_jprb, 0.9776_jprb, 0.9776_jprb, 0.9776_jprb, 0.9776_jprb, &
      0.9776_jprb, 0.9776_jprb, 0.9776_jprb, 0.9776_jprb, 0.9776_jprb, 0.9775_jprb, 0.9775_jprb, 0.9775_jprb, &
      0.9776_jprb, 0.9776_jprb, 0.9776_jprb, 0.9776_jprb, 0.9777_jprb, 0.9777_jprb, 0.9777_jprb, 0.9777_jprb, &
      0.9777_jprb, 0.9777_jprb, 0.9777_jprb, 0.9777_jprb, 0.9777_jprb, 0.9776_jprb, 0.9776_jprb, 0.9776_jprb, &
      0.9775_jprb, 0.9775_jprb, 0.9774_jprb, 0.9773_jprb, 0.9773_jprb, 0.9773_jprb, 0.9773_jprb, 0.9773_jprb, &
      0.9774_jprb, 0.9774_jprb, 0.9775_jprb, 0.9776_jprb, 0.9777_jprb, 0.9778_jprb, 0.9779_jprb, 0.9780_jprb, &
      0.9781_jprb, 0.9782_jprb, 0.9783_jprb, 0.9785_jprb, 0.9786_jprb, 0.9788_jprb, 0.9790_jprb, 0.9792_jprb, &
      0.9793_jprb, 0.9795_jprb, 0.9797_jprb, 0.9799_jprb, 0.9801_jprb, 0.9802_jprb, 0.9803_jprb, 0.9805_jprb, &
      0.9806_jprb, 0.9807_jprb, 0.9808_jprb, 0.9809_jprb, 0.9810_jprb, 0.9811_jprb, 0.9811_jprb, 0.9811_jprb, &
      0.9810_jprb, 0.9810_jprb, 0.9810_jprb, 0.9809_jprb, 0.9808_jprb, 0.9808_jprb, 0.9807_jprb, 0.9807_jprb, &
      0.9806_jprb, 0.9805_jprb, 0.9805_jprb, 0.9804_jprb, 0.9803_jprb, 0.9802_jprb, 0.9802_jprb, 0.9801_jprb, &
      0.9800_jprb, 0.9799_jprb, 0.9798_jprb, 0.9797_jprb, 0.9797_jprb, 0.9795_jprb, 0.9795_jprb, 0.9794_jprb, &
      0.9793_jprb, 0.9792_jprb, 0.9791_jprb, 0.9791_jprb, 0.9789_jprb, 0.9789_jprb, 0.9788_jprb, 0.9787_jprb, &
      0.9786_jprb, 0.9785_jprb, 0.9785_jprb, 0.9783_jprb, 0.9783_jprb, 0.9782_jprb, 0.9781_jprb, 0.9781_jprb, &
      0.9780_jprb, 0.9779_jprb, 0.9779_jprb, 0.9778_jprb, 0.9777_jprb, 0.9777_jprb, 0.9776_jprb, 0.9775_jprb, &
      0.9774_jprb, 0.9774_jprb, 0.9773_jprb, 0.9772_jprb, 0.9771_jprb, 0.9771_jprb, 0.9770_jprb, 0.9769_jprb, &
      0.9769_jprb, 0.9768_jprb, 0.9767_jprb, 0.9766_jprb, 0.9765_jprb, 0.9765_jprb, 0.9764_jprb, 0.9764_jprb, &
      0.9763_jprb, 0.9762_jprb, 0.9762_jprb, 0.9761_jprb, 0.9761_jprb, 0.9760_jprb, 0.9759_jprb, 0.9758_jprb, &
      0.9757_jprb, 0.9757_jprb, 0.9756_jprb, 0.9756_jprb, 0.9755_jprb, 0.9755_jprb, 0.9754_jprb, 0.9754_jprb, &
      0.9754_jprb, 0.9754_jprb, 0.9754_jprb, 0.9753_jprb, 0.9753_jprb, 0.9753_jprb, 0.9753_jprb, 0.9752_jprb, &
      0.9752_jprb, 0.9752_jprb, 0.9753_jprb, 0.9754_jprb, 0.9755_jprb, 0.9756_jprb, 0.9757_jprb, 0.9757_jprb, &
      0.9758_jprb, 0.9758_jprb, 0.9759_jprb, 0.9759_jprb, 0.9760_jprb, 0.9760_jprb, 0.9760_jprb, 0.9760_jprb, &
      0.9761_jprb, 0.9761_jprb, 0.9761_jprb, 0.9761_jprb, 0.9761_jprb, 0.9761_jprb, 0.9761_jprb, 0.9761_jprb, &
      0.9761_jprb, 0.9760_jprb, 0.9760_jprb, 0.9759_jprb, 0.9759_jprb, 0.9758_jprb, 0.9758_jprb, 0.9757_jprb, &
      0.9757_jprb, 0.9757_jprb, 0.9756_jprb, 0.9756_jprb, 0.9755_jprb, 0.9754_jprb, 0.9753_jprb, 0.9753_jprb, &
      0.9752_jprb, 0.9751_jprb, 0.9751_jprb, 0.9750_jprb, 0.9750_jprb, 0.9749_jprb, 0.9749_jprb, 0.9748_jprb, &
      0.9747_jprb, 0.9746_jprb, 0.9746_jprb, 0.9746_jprb, 0.9745_jprb, 0.9744_jprb, 0.9743_jprb, 0.9742_jprb, &
      0.9742_jprb, 0.9741_jprb, 0.9740_jprb, 0.9739_jprb, 0.9739_jprb, 0.9739_jprb, 0.9738_jprb, 0.9737_jprb, &
      0.9736_jprb, 0.9735_jprb, 0.9735_jprb, 0.9734_jprb, 0.9733_jprb, 0.9732_jprb, 0.9731_jprb, 0.9731_jprb, &
      0.9730_jprb, 0.9729_jprb, 0.9728_jprb, 0.9727_jprb, 0.9726_jprb, 0.9725_jprb, 0.9724_jprb, 0.9723_jprb, &
      0.9723_jprb, 0.9722_jprb, 0.9721_jprb, 0.9720_jprb, 0.9719_jprb, 0.9718_jprb, 0.9717_jprb, 0.9716_jprb, &
      0.9715_jprb, 0.9714_jprb, 0.9713_jprb, 0.9712_jprb, 0.9711_jprb, 0.9709_jprb, 0.9708_jprb, 0.9706_jprb, &
      0.9705_jprb, 0.9704_jprb, 0.9703_jprb, 0.9702_jprb, 0.9700_jprb, 0.9699_jprb, 0.9698_jprb, 0.9696_jprb, &
      0.9695_jprb, 0.9693_jprb, 0.9691_jprb, 0.9690_jprb, 0.9688_jprb, 0.9686_jprb, 0.9685_jprb, 0.9683_jprb, &
      0.9682_jprb, 0.9681_jprb, 0.9679_jprb, 0.9677_jprb, 0.9676_jprb, 0.9674_jprb, 0.9671_jprb, 0.9669_jprb/

  DATA snow_stdv / 0.015_jprb /
  DATA snow_em / &
      0.9716_jprb, 0.9716_jprb, 0.9716_jprb, 0.9716_jprb, 0.9713_jprb, 0.9710_jprb, 0.9708_jprb, 0.9706_jprb, &
      0.9705_jprb, 0.9705_jprb, 0.9705_jprb, 0.9703_jprb, 0.9701_jprb, 0.9700_jprb, 0.9699_jprb, 0.9700_jprb, &
      0.9702_jprb, 0.9703_jprb, 0.9705_jprb, 0.9707_jprb, 0.9710_jprb, 0.9714_jprb, 0.9717_jprb, 0.9722_jprb, &
      0.9728_jprb, 0.9734_jprb, 0.9740_jprb, 0.9746_jprb, 0.9753_jprb, 0.9759_jprb, 0.9765_jprb, 0.9771_jprb, &
      0.9778_jprb, 0.9784_jprb, 0.9792_jprb, 0.9798_jprb, 0.9806_jprb, 0.9814_jprb, 0.9824_jprb, 0.9833_jprb, &
      0.9842_jprb, 0.9852_jprb, 0.9863_jprb, 0.9873_jprb, 0.9882_jprb, 0.9891_jprb, 0.9901_jprb, 0.9908_jprb, &
      0.9914_jprb, 0.9920_jprb, 0.9925_jprb, 0.9926_jprb, 0.9928_jprb, 0.9927_jprb, 0.9926_jprb, 0.9926_jprb, &
      0.9923_jprb, 0.9920_jprb, 0.9918_jprb, 0.9916_jprb, 0.9915_jprb, 0.9913_jprb, 0.9911_jprb, 0.9907_jprb, &
      0.9905_jprb, 0.9903_jprb, 0.9902_jprb, 0.9900_jprb, 0.9897_jprb, 0.9896_jprb, 0.9894_jprb, 0.9892_jprb, &
      0.9890_jprb, 0.9889_jprb, 0.9886_jprb, 0.9884_jprb, 0.9883_jprb, 0.9884_jprb, 0.9885_jprb, 0.9885_jprb, &
      0.9884_jprb, 0.9883_jprb, 0.9881_jprb, 0.9880_jprb, 0.9880_jprb, 0.9880_jprb, 0.9880_jprb, 0.9879_jprb, &
      0.9879_jprb, 0.9879_jprb, 0.9879_jprb, 0.9879_jprb, 0.9879_jprb, 0.9879_jprb, 0.9878_jprb, 0.9877_jprb, &
      0.9876_jprb, 0.9876_jprb, 0.9877_jprb, 0.9876_jprb, 0.9875_jprb, 0.9874_jprb, 0.9873_jprb, 0.9873_jprb, &
      0.9873_jprb, 0.9874_jprb, 0.9875_jprb, 0.9875_jprb, 0.9874_jprb, 0.9874_jprb, 0.9874_jprb, 0.9874_jprb, &
      0.9874_jprb, 0.9874_jprb, 0.9874_jprb, 0.9874_jprb, 0.9874_jprb, 0.9874_jprb, 0.9874_jprb, 0.9873_jprb, &
      0.9873_jprb, 0.9873_jprb, 0.9873_jprb, 0.9873_jprb, 0.9873_jprb, 0.9872_jprb, 0.9871_jprb, 0.9872_jprb, &
      0.9871_jprb, 0.9870_jprb, 0.9870_jprb, 0.9870_jprb, 0.9870_jprb, 0.9869_jprb, 0.9868_jprb, 0.9868_jprb, &
      0.9867_jprb, 0.9866_jprb, 0.9866_jprb, 0.9865_jprb, 0.9865_jprb, 0.9865_jprb, 0.9866_jprb, 0.9866_jprb, &
      0.9865_jprb, 0.9865_jprb, 0.9865_jprb, 0.9866_jprb, 0.9866_jprb, 0.9866_jprb, 0.9866_jprb, 0.9867_jprb, &
      0.9868_jprb, 0.9868_jprb, 0.9868_jprb, 0.9867_jprb, 0.9867_jprb, 0.9867_jprb, 0.9866_jprb, 0.9867_jprb, &
      0.9867_jprb, 0.9867_jprb, 0.9867_jprb, 0.9867_jprb, 0.9867_jprb, 0.9868_jprb, 0.9868_jprb, 0.9868_jprb, &
      0.9869_jprb, 0.9869_jprb, 0.9870_jprb, 0.9872_jprb, 0.9873_jprb, 0.9874_jprb, 0.9874_jprb, 0.9874_jprb, &
      0.9875_jprb, 0.9875_jprb, 0.9875_jprb, 0.9875_jprb, 0.9876_jprb, 0.9876_jprb, 0.9876_jprb, 0.9876_jprb, &
      0.9877_jprb, 0.9877_jprb, 0.9877_jprb, 0.9877_jprb, 0.9878_jprb, 0.9879_jprb, 0.9879_jprb, 0.9879_jprb, &
      0.9878_jprb, 0.9878_jprb, 0.9878_jprb, 0.9879_jprb, 0.9879_jprb, 0.9879_jprb, 0.9879_jprb, 0.9878_jprb, &
      0.9877_jprb, 0.9876_jprb, 0.9876_jprb, 0.9877_jprb, 0.9877_jprb, 0.9876_jprb, 0.9876_jprb, 0.9876_jprb, &
      0.9876_jprb, 0.9876_jprb, 0.9876_jprb, 0.9875_jprb, 0.9874_jprb, 0.9873_jprb, 0.9873_jprb, 0.9872_jprb, &
      0.9870_jprb, 0.9869_jprb, 0.9869_jprb, 0.9869_jprb, 0.9869_jprb, 0.9869_jprb, 0.9868_jprb, 0.9867_jprb, &
      0.9867_jprb, 0.9866_jprb, 0.9866_jprb, 0.9865_jprb, 0.9865_jprb, 0.9865_jprb, 0.9864_jprb, 0.9863_jprb, &
      0.9862_jprb, 0.9862_jprb, 0.9862_jprb, 0.9862_jprb, 0.9862_jprb, 0.9861_jprb, 0.9860_jprb, 0.9860_jprb, &
      0.9860_jprb, 0.9859_jprb, 0.9859_jprb, 0.9859_jprb, 0.9858_jprb, 0.9858_jprb, 0.9857_jprb, 0.9857_jprb, &
      0.9857_jprb, 0.9856_jprb, 0.9855_jprb, 0.9855_jprb, 0.9854_jprb, 0.9854_jprb, 0.9853_jprb, 0.9853_jprb, &
      0.9852_jprb, 0.9852_jprb, 0.9852_jprb, 0.9852_jprb, 0.9852_jprb, 0.9852_jprb, 0.9851_jprb, 0.9850_jprb, &
      0.9850_jprb, 0.9850_jprb, 0.9850_jprb, 0.9851_jprb, 0.9850_jprb, 0.9850_jprb, 0.9849_jprb, 0.9849_jprb, &
      0.9850_jprb, 0.9851_jprb, 0.9851_jprb, 0.9850_jprb, 0.9850_jprb, 0.9849_jprb, 0.9849_jprb, 0.9849_jprb, &
      0.9849_jprb, 0.9848_jprb, 0.9848_jprb, 0.9848_jprb, 0.9849_jprb, 0.9849_jprb, 0.9849_jprb, 0.9849_jprb, &
      0.9849_jprb, 0.9849_jprb, 0.9849_jprb, 0.9848_jprb, 0.9848_jprb, 0.9849_jprb, 0.9849_jprb, 0.9850_jprb, &
      0.9850_jprb, 0.9850_jprb, 0.9851_jprb, 0.9851_jprb, 0.9852_jprb, 0.9852_jprb, 0.9853_jprb, 0.9853_jprb, &
      0.9854_jprb, 0.9854_jprb, 0.9854_jprb, 0.9854_jprb, 0.9854_jprb, 0.9854_jprb, 0.9854_jprb, 0.9854_jprb, &
      0.9855_jprb, 0.9856_jprb, 0.9856_jprb, 0.9856_jprb, 0.9856_jprb, 0.9856_jprb, 0.9856_jprb, 0.9856_jprb, &
      0.9856_jprb, 0.9855_jprb, 0.9855_jprb, 0.9855_jprb, 0.9854_jprb, 0.9853_jprb, 0.9853_jprb, 0.9853_jprb, &
      0.9853_jprb, 0.9853_jprb, 0.9853_jprb, 0.9853_jprb, 0.9853_jprb, 0.9853_jprb, 0.9853_jprb, 0.9853_jprb, &
      0.9852_jprb, 0.9851_jprb, 0.9851_jprb, 0.9850_jprb, 0.9849_jprb, 0.9849_jprb, 0.9849_jprb, 0.9848_jprb, &
      0.9848_jprb, 0.9848_jprb, 0.9848_jprb, 0.9848_jprb, 0.9848_jprb, 0.9848_jprb, 0.9847_jprb, 0.9846_jprb, &
      0.9846_jprb, 0.9846_jprb, 0.9847_jprb, 0.9846_jprb, 0.9845_jprb, 0.9844_jprb, 0.9844_jprb, 0.9843_jprb, &
      0.9842_jprb, 0.9842_jprb, 0.9842_jprb, 0.9842_jprb, 0.9841_jprb, 0.9841_jprb, 0.9840_jprb, 0.9839_jprb, &
      0.9838_jprb, 0.9838_jprb, 0.9837_jprb, 0.9837_jprb, 0.9837_jprb, 0.9836_jprb, 0.9836_jprb, 0.9835_jprb, &
      0.9835_jprb, 0.9834_jprb, 0.9833_jprb, 0.9832_jprb, 0.9832_jprb, 0.9832_jprb, 0.9831_jprb, 0.9831_jprb, &
      0.9830_jprb, 0.9829_jprb, 0.9828_jprb, 0.9828_jprb, 0.9827_jprb, 0.9827_jprb, 0.9826_jprb, 0.9826_jprb, &
      0.9825_jprb, 0.9824_jprb, 0.9824_jprb, 0.9823_jprb, 0.9821_jprb, 0.9821_jprb, 0.9820_jprb, 0.9821_jprb, &
      0.9820_jprb, 0.9820_jprb, 0.9818_jprb, 0.9817_jprb, 0.9817_jprb, 0.9816_jprb, 0.9815_jprb, 0.9815_jprb, &
      0.9815_jprb, 0.9814_jprb, 0.9813_jprb, 0.9813_jprb, 0.9812_jprb, 0.9811_jprb, 0.9811_jprb, 0.9810_jprb/

CONTAINS

!------------------------------------------
! Routines for initialising database
!------------------------------------------

  SUBROUTINE rttov_uwiremis_init(    &
        &             path,          &! in
        &             imonth,        &! in
        &             verbose,       &! in
        &             err,           &! out
        &             instr_wavenum  )! in, optional

    ! Description:
    ! initialize the rttov_uwiremis algorithm by (1) reading in the UW BF IR Global
    ! Emissivty data and (2) the eigenvectors of the laboratory spectra, and (2) make some
    ! precalculations for the PCA regression
    !
    ! History:
    ! Version   Date     Comment
    ! -------   ----     -------
    !  0.9    03/31/2009   origianl code E. Borbas UW-Madison/CIMSS
    !  1.0    03/31/2009  New F90 code with structures (E Borbas B Ruston)

    USE rttov_const, ONLY : &
      errorstatus_success, &
      errorstatus_fatal

    CHARACTER(LEN=*),   INTENT(IN)           :: path
    INTEGER(KIND=jpim), INTENT(IN)           :: imonth
    LOGICAL(KIND=jplm), INTENT(IN)           :: verbose
    INTEGER(KIND=jpim), INTENT(OUT)          :: err
    REAL(KIND=jprb),    INTENT(IN), OPTIONAL :: instr_wavenum(:)

    CHARACTER(LEN=300) :: fn
    CHARACTER(LEN=4)   :: cyear
    CHARACTER(LEN=2)   :: cmonth
    LOGICAL(KIND=jplm) :: file_exists
    CHARACTER(LEN=128) :: errmsg

#ifdef _RTTOV_HDF
    LOGICAL(KIND=jplm) :: hdf_was_open, hdf_was_64bit_reals
#endif

    TRY

    CALL rttov_uwiremis_nullify_pointers()

    WRITE(cyear, '(i4)') db_ver_year
    WRITE(cmonth, '(i2.2)') imonth

#ifdef _RTTOV_HDF
    hdf_was_open = is_hdf_open
    hdf_was_64bit_reals = is_hdf_64bit_reals
    CALL open_hdf(.TRUE._jplm, err)
    THROW(err.NE.0)
    WRITE(errmsg,'(a)') ', RTTOV was compiled with HDF5 so HDF5 atlas files are required'
#else
    WRITE(errmsg,'(a)') ', RTTOV was not compiled with HDF5 so netCDF atlas files are required'
#endif

    !----------------------------------------------------------------------------
    ! reading the 0.1 degree resolution PCA Coefficients of UW BF IR Land Surface Global Emissivity
    !----------------------------------------------------------------------------
    fn = TRIM(path)//'UWirbfemis_COEF_V2.1_0.1deg_'//cyear//cmonth//'_mask'
#ifdef _RTTOV_HDF
    fn = TRIM(fn)//'.H5'
#else
    fn = TRIM(fn)//'.nc'
#endif
    INQUIRE(FILE=fn, EXIST=file_exists)
    IF (file_exists) THEN
      CALL rttov_uwiremis_read_coefs(err, TRIM(fn))
      THROW(err.NE.0)
    ELSE
      err = errorstatus_fatal
      THROWM(err .NE. errorstatus_success, 'UWiremis PCA coefs file not found: '//TRIM(fn)//TRIM(errmsg))
    ENDIF
    IF (verbose) INFO('Using UWiremis coefs: '//TRIM(fn))

    !----------------------------------------------------------------------------
    ! reading the 0.5 degree resolution UW IR Land Surface Global Emissivty STD DEV
    !----------------------------------------------------------------------------
    IF (ir_atlas_std_init) THEN
      fn = TRIM(path)//'UWiremis_hsremis_covmat_V1.0_deg0.5_month'//cmonth//'_mask'
#ifdef _RTTOV_HDF
      fn = TRIM(fn)//'.H5'
#else
      fn = TRIM(fn)//'.nc'
#endif
      INQUIRE(FILE=fn, EXIST=file_exists)
      IF (file_exists) THEN
        CALL rttov_uwiremis_read_cov(err, TRIM(fn))
        THROW(err.NE.0)
      ELSE
        err = errorstatus_fatal
        THROWM(err .NE. errorstatus_success, 'UWiremis covariances file not found: '//TRIM(fn)//TRIM(errmsg))
      ENDIF
      IF (verbose) INFO('Using UWiremis covariances: '//TRIM(fn))
    ENDIF

    !-------------------------------------------------------------------
    !  reading the eigienvectors of the 128 selected laboratory spectra 
    !-------------------------------------------------------------------
    fn = TRIM(path)//'UWiremis_labeigvects'
#ifdef _RTTOV_HDF
    fn = TRIM(fn)//'.H5'
#else
    fn = TRIM(fn)//'.nc'
#endif
    INQUIRE(FILE=fn, EXIST=file_exists)
    IF (file_exists) THEN
      CALL rttov_uwiremis_read_labevecs(err, TRIM(fn))
      THROW(err.NE.0)
    ELSE
      err = errorstatus_fatal
      THROWM(err .NE. errorstatus_success, 'UWiremis eigenvector file not found: '//TRIM(fn)//TRIM(errmsg))
    ENDIF
    IF (verbose) INFO('Using UWiremis eigenvector file: '//TRIM(fn))


    IF (ir_atlas_do_ang_corr) THEN
      !-------------------------------------------------------------------
      !  reading the IGBP ecosystem map file
      !-------------------------------------------------------------------
      fn = TRIM(path)//'uwiremis_igbpmap'
#ifdef _RTTOV_HDF
      fn = TRIM(fn)//'.H5'
#else
      fn = TRIM(fn)//'.nc'
#endif
      INQUIRE(FILE=fn, EXIST=file_exists)
      IF (file_exists) THEN
        CALL rttov_uwiremis_read_igbp(err, TRIM(fn))
        THROW(err.NE.0)
      ELSE
        err = errorstatus_fatal
        THROWM(err .NE. errorstatus_success, 'UWiremis IGBP atlas file not found: '//TRIM(fn)//TRIM(errmsg))
      ENDIF
      IF (verbose) INFO('Using UWiremis IGBP atlas file: '//TRIM(fn))

      !-------------------------------------------------------------------
      !  reading the angular correcction coefficient file
      !-------------------------------------------------------------------
      IF (imonth == 12 .OR. imonth == 1 .OR. imonth == 2) THEN
        fn = TRIM(path)//'uwiremis_angcorr_201301_V1.5'
      ELSE IF (imonth >= 3 .AND. imonth <= 5) THEN
        fn = TRIM(path)//'uwiremis_angcorr_201204_V1.5'
      ELSE IF (imonth >= 6 .AND. imonth <= 8) THEN
        fn = TRIM(path)//'uwiremis_angcorr_201207_V1.5'
      ELSE
        fn = TRIM(path)//'uwiremis_angcorr_201210_V1.5'
      ENDIF

#ifdef _RTTOV_HDF
      fn = TRIM(fn)//'.H5'
#else
      fn = TRIM(fn)//'.nc'
#endif
      INQUIRE(FILE=fn, EXIST=file_exists)
      IF (file_exists) THEN
        CALL rttov_uwiremis_read_angfunc(err, TRIM(fn))
        THROW(err.NE.0)
      ELSE
        err = errorstatus_fatal
        THROWM(err .NE. errorstatus_success, 'UWiremis angle correction file not found: '//TRIM(fn)//TRIM(errmsg))
      ENDIF
      IF (verbose) INFO('Using UWiremis angle correction file: '//TRIM(fn))
    ENDIF

#ifdef _RTTOV_HDF
    ! If HDF5 lib was open before rttov_read_coefs was called, make sure the real
    ! kind is the same as it was previously. Otherwise close the library.
    IF (hdf_was_open) THEN
      CALL open_hdf(hdf_was_64bit_reals, err)
      THROW(err.NE.0)
    ELSE
      CALL close_hdf(err)
      THROW(err.NE.0)
    ENDIF
#endif

    IF (PRESENT(instr_wavenum)) THEN
      ! If a channel wavenumber list is supplied for a particular instrument
      ! we can retrieve emissivities much faster.

      ! Interpolate hsr data onto the channel wavenumbers
      ncoefchans = SIZE(instr_wavenum)
      IF (ir_atlas_do_ang_corr) THEN
        ALLOCATE(pcu_int(numpcs,ncoefchans,2), pcm_int(ncoefchans,2))
        ALLOCATE(p1d_int(ncoefchans,nb_igbp,2), &
                 p2d_int(ncoefchans,nb_igbp,2), &
                 p3d_int(ncoefchans,nb_igbp,2), &
                 p1n_int(ncoefchans,nb_igbp,2), &
                 p2n_int(ncoefchans,nb_igbp,2), &
                 p3n_int(ncoefchans,nb_igbp,2))
      ELSE
        ALLOCATE(pcu_int(numpcs,ncoefchans,1), pcm_int(ncoefchans,1))
      ENDIF
      ALLOCATE(sice_em_int(ncoefchans), snow_em_int(ncoefchans))
      IF (ir_atlas_std_init) ALLOCATE(cov_emis_int(cv_pack,ncoefchans))

      CALL rttov_uwiremis_hsr_interp(instr_wavenum)

      ! HSR data is no longer required
      IF (ir_atlas_std_init) DEALLOCATE(cov_emis)
      IF (ir_atlas_do_ang_corr) DEALLOCATE(p1d, p2d, p3d, p1n, p2n, p3n)
    ENDIF

    CATCH
  END SUBROUTINE rttov_uwiremis_init

#ifndef _RTTOV_HDF
  SUBROUTINE  rttov_uwiremis_read_cov(err, fn)

    ! Description:
    ! read the 0.1 degree resolution UW BF IR Global Emissivty data
    ! from the netCDF file into memory.
    ! 
    ! History:
    ! Version   Date     Comment
    ! -------   ----     -------
    !  1.0    03/31/2009   origianl code B. Ruston 

    INTEGER(KIND=jpim), INTENT(OUT) :: err
    CHARACTER(LEN=*),   INTENT(IN)  :: fn    ! filename including full path

    INTEGER(KIND=jpim) :: nvars     ! number of variables
    INTEGER(KIND=jpim) :: ndims     ! number of dimensions
    INTEGER(KIND=jpim) :: errstat   ! error code
    INTEGER(KIND=jpim) :: recdim    ! record dimension
    INTEGER(KIND=jpim) :: nc_dim(4)

    INTEGER(KIND=jpim) :: i, j
    INTEGER(KIND=jpim) :: ncid, ngatts, nrecs, varid

    INTEGER(KIND=jpis), ALLOCATABLE :: emis_cov(:,:)
    INTEGER(KIND=jpit), ALLOCATABLE :: pack_cov(:,:)

    CHARACTER(LEN=1024) :: strbuf ! string buffer for var
    INTEGER(KIND=jpim)  :: indexlut

    TRY

    ! Open netCDF file.
    errstat = nf_open(TRIM(fn),nf_nowrite,ncid)

    ! Get info on the record dimension for this file.
    errstat = nf_inq(ncid, ndims, nvars, ngatts, recdim)

    DO recdim=1,ndims
      errstat = nf_inq_dim(ncid,recdim,strbuf,nrecs)
      nc_dim(recdim) = nrecs
    ENDDO

    ! Retrieve database of diagonal of covariance of MODIS emissivity values
    cv_lats = nc_dim(2)
    cv_lons = nc_dim(3)
    ALLOCATE(pack_cov(cv_lats,cv_lons), cov_emis_lut(cv_lats,cv_lons))
    errstat = nf_inq_varid (ncid, 'mask', varid)
    errstat = nf_get_var_int1(ncid, varid, pack_cov)

    ! Generate the look-up table into the covariance data
    cov_emis_lut(:,:) = -1
    indexlut = 1
    DO i = 1, nc_dim(3)
      DO j = 1, nc_dim(2)
        IF (pack_cov(j,i) > 0) THEN
          cov_emis_lut(j,i) = indexlut
          indexlut = indexlut + 1
        ENDIF
      ENDDO
    ENDDO
    DEALLOCATE(pack_cov)

    ! Retrieve database of diagonal of covariance of MODIS emissivity values
    IF (nc_dim(1) .NE. numwave) THEN
      err = errorstatus_fatal
      THROWM(err.NE.0,"Inconsistent dimensions in covariance data: numwave")
    ENDIF
    cv_pack = nc_dim(4)
    ALLOCATE(emis_cov(numwave,cv_pack), cov_emis(cv_pack,numwave))
    errstat = nf_inq_varid (ncid, 'emis_diagCov', varid)
    errstat = nf_get_var_int2(ncid, varid, emis_cov)
    errstat = nf_get_att_real(ncid, varid, 'scale_factor', cov_sfac)

    cov_emis = TRANSPOSE(emis_cov)
    DEALLOCATE(emis_cov)

    errstat = nf_close(ncid)

    CATCH
  END SUBROUTINE rttov_uwiremis_read_cov
#else
  SUBROUTINE rttov_uwiremis_read_cov(err, fn)

    ! Description:
    ! read the 0.1 degree resolution MODIS BRDF kernel parameters data
    ! from the HDF5 file into memory.
    !
    ! History:
    ! Version   Date     Comment
    ! -------   ----     -------
    !  1.1      03/06/2014 original code J. Vidot

    INTEGER(KIND=jpim), INTENT(OUT) :: err
    CHARACTER(LEN=*),   INTENT(IN)  :: fn

    INTEGER(KIND=jpim)          :: indexlut
    INTEGER(KIND=jpim)          :: i,j
    INTEGER(KIND=jpis), POINTER :: emis_cov(:,:)
    INTEGER(KIND=jpit), POINTER :: pack_cov(:,:)

#include "rttov_hdf_load.interface"

    TRY

    CALL RTTOV_HDF_LOAD(err, fn, "/EMIS", SNAME="MASK", PIT2=pack_cov)
    THROWM(ERR.NE.0,"Cannot load Emissivity Covariance flag from "//TRIM(FN))

    cv_lats = SIZE(pack_cov,1)
    cv_lons = SIZE(pack_cov,2)
    ALLOCATE(cov_emis_lut(cv_lats,cv_lons))

    ! Generate the look-up table into the covariance data
    cov_emis_lut(:,:) = -1
    indexlut = 1
    DO i = 1, cv_lons
      DO j = 1, cv_lats
        IF (pack_cov(j,i) > 0) THEN
          cov_emis_lut(j,i) = indexlut
          indexlut = indexlut + 1
        ENDIF
      ENDDO
    ENDDO
    DEALLOCATE(pack_cov)
    NULLIFY(pack_cov)

    CALL rttov_hdf_load(err, fn, "/EMIS", SNAME="DIAGCOV", PIS2=emis_cov)
    THROWM(ERR.NE.0,"Cannot get covariance values from "//TRIM(fn))

    IF (SIZE(emis_cov,1) .NE. numwave) THEN
      DEALLOCATE(emis_cov)
      err = errorstatus_fatal
      THROWM(err.NE.0,"Inconsistent dimensions in covariance data: numwave")
    ENDIF
    cv_pack = SIZE(emis_cov,2)
    ALLOCATE(cov_emis(cv_pack,numwave))

    cov_emis = TRANSPOSE(emis_cov)
    DEALLOCATE(emis_cov)
    NULLIFY(emis_cov)

    CALL rttov_hdf_load(err, fn, "/EMIS", SNAME="DIAGCOV_SFAC", RM0=cov_sfac)
    THROWM(ERR.NE.0,"Cannot get covariance scale factor from "//TRIM(fn))

    CATCH
  END SUBROUTINE rttov_uwiremis_read_cov
#endif

#ifndef _RTTOV_HDF
  SUBROUTINE  rttov_uwiremis_read_coefs(err, fn)

    ! Description:
    ! read the 0.1 degree resolution UW BF IR Global Emissivty data
    ! from the netCDF file into memory. 
    ! 
    ! History:
    ! Version   Date     Comment
    ! -------   ----     -------
    !  1.1    11/30/2012   modified for coef files by E. Borbas
    !  1.0    03/31/2009   original code B. Ruston 

    INTEGER(KIND=jpim), INTENT(OUT) :: err
    CHARACTER(LEN=*),   INTENT(IN)  :: fn    ! filename including full path

    INTEGER(KIND=jpim)  :: nvars     ! number of variables
    INTEGER(KIND=jpim)  :: ndims     ! number of dimensions
    INTEGER(KIND=jpim)  :: errstat   ! error code
    INTEGER(KIND=jpim)  :: recdim    ! record dimension
    INTEGER(KIND=jpim)  :: nc_dim(4) ! hng_pnt, lats, lons, pack_len

    INTEGER(KIND=jpim)  :: i, j, k
    INTEGER(KIND=jpim)  :: ncid, ngatts, nrecs, varid

    CHARACTER(LEN=1024) :: strbuf ! string buffer for var
    CHARACTER(LEN=6)    :: cfld
    INTEGER(KIND=jpim)  :: indexlut

    TRY

    ! Open netCDF file.
    errstat = nf_open(TRIM(fn), nf_nowrite, ncid)

    ! Get info on the record dimension for this file.
    errstat = nf_inq(ncid, ndims, nvars, ngatts, recdim)

    DO recdim = 1, ndims
      errstat = nf_inq_dim(ncid, recdim, strbuf, nrecs)
      nc_dim(recdim) = nrecs
    ENDDO

    nb_lats = nc_dim(2)
    nb_lons = nc_dim(3)
    ALLOCATE(bfemis_flag(nb_lats,nb_lons), bfemis_lut(nb_lats,nb_lons))

    ! Retrieve emissivity database flag value
    errstat = nf_inq_varid(ncid, 'emis_flag', varid)
    errstat = nf_get_var_int2(ncid, varid, bfemis_flag)

    ! Generate the look-up table into the emissivity data
    bfemis_lut(:,:) = -1_jpim
    indexlut = 1_jpim
    DO i = 1, nb_lons
      DO j = 1, nb_lats
        IF (bfemis_flag(j,i) > 0) THEN
          bfemis_lut(j,i) = indexlut
          indexlut = indexlut + 1_jpim
        ENDIF
      ENDDO
    ENDDO

    ! Retrieve database of 6 coefs
    nb_pack = nc_dim(4)
    IF (nc_dim(1) .NE. numpcs) THEN
      err = errorstatus_fatal
      THROWM(err.NE.0,"Inconsistent dimensions in emissivity data: numpcs")
    ENDIF
    ALLOCATE(pca_coef(nb_pack,numpcs))
    DO k = 1, numpcs
      IF (k < 10) THEN
        WRITE(cfld, '("coef",i1," ")') k
      ELSE
        WRITE(cfld, '("coef",i2.2)') k
      ENDIF
      errstat = nf_inq_varid (ncid, cfld, varid)
      errstat = nf_get_var_real(ncid, varid, pca_coef(:,k))
      errstat = nf_get_att_real(ncid, varid, 'scale_factor', pca_sfac(k))
      errstat = nf_get_att_real(ncid, varid, 'add_offset', pca_offs(k))
    ENDDO

    errstat = nf_close(ncid)

    CATCH
  END SUBROUTINE rttov_uwiremis_read_coefs
#else
  SUBROUTINE rttov_uwiremis_read_coefs(err, fn)

    ! Description:
    ! read the 0.1 degree resolution MODIS BRDF kernel parameters data
    ! from the HDF5 file into memory.
    !
    ! History:
    ! Version   Date     Comment
    ! -------   ----     -------
    !  1.1      03/06/2014 original code J. Vidot

    INTEGER(KIND=jpim), INTENT(OUT) :: err
    CHARACTER(LEN=*),   INTENT(IN)  :: fn

    INTEGER(KIND=jpim)          :: indexlut
    INTEGER(KIND=jpim)          :: i, j

#include "rttov_hdf_load.interface"
    TRY

    CALL RTTOV_HDF_LOAD(ERR, fn, "/EMIS", SNAME="FLAG", PIS2=bfemis_flag)
    THROWM(ERR.NE.0,"Cannot load emissivity flag from "//TRIM(fn))

    nb_lats = SIZE(bfemis_flag,1)
    nb_lons = SIZE(bfemis_flag,2)
    ALLOCATE(bfemis_lut(nb_lats,nb_lons))

    ! Generate the look-up table into the emissivity data
    bfemis_lut(:,:) = -1_jpim
    indexlut = 1_jpim
    DO i = 1, nb_lons
      DO j = 1, nb_lats
        IF (bfemis_flag(j,i) > 0) THEN
          bfemis_lut(j,i) = indexlut
          indexlut = indexlut + 1_jpim
        ENDIF
      ENDDO
    ENDDO

    CALL rttov_hdf_load(err, fn, "/EMIS", SNAME="COEF", PRM2=pca_coef)
    THROWM(ERR.NE.0,"Cannot get PCA coeffcients from "//TRIM(fn))

    nb_pack = SIZE(pca_coef,1)
    IF (SIZE(pca_coef,2) .NE. numpcs) THEN
      DEALLOCATE(pca_coef)
      err = errorstatus_fatal
      THROWM(err.NE.0,"Inconsistent dimensions in emissivity data: numpcs")
    ENDIF

    CALL rttov_hdf_load(err, fn, "/EMIS", SNAME="COEF_SFAC", PRM1=pca_sfac)
    THROWM(ERR.NE.0,"Cannot get PCA coeffcients scale factor from "//TRIM(fn))

    CALL rttov_hdf_load(err, fn, "/EMIS", SNAME="COEF_OFFS", PRM1=pca_offs)
    THROWM(ERR.NE.0,"Cannot get PCA coeffcients offset from "//TRIM(fn))

    CATCH
  END SUBROUTINE rttov_uwiremis_read_coefs
#endif

#ifndef _RTTOV_HDF
  SUBROUTINE  rttov_uwiremis_read_labevecs(err, fn)

    ! Description:
    ! read the eigienvectors of the 128 selected laboratory spectra
    !(created by E borbas) from the netCDF file into memory.
    !
    ! History:
    ! Version   Date     Comment
    ! -------   ----     -------
    !  1.0    03/31/2009   origianl code B. Ruston

    INTEGER(KIND=jpim), INTENT(OUT) :: err
    CHARACTER(LEN=*),   INTENT(IN)  :: fn

    INTEGER(KIND=jpim) :: errstat
    INTEGER(KIND=jpim) :: ncid, varid
    INTEGER(KIND=jpim) :: ndims, nvars, ngatts, recdim, nrecs
    INTEGER(KIND=jpim) :: nc_dim(2)

    CHARACTER (len=1024) :: strbuf ! string buffer for var

    TRY

    ! Open netCDF file.
    errstat = nf_open(TRIM(fn), nf_nowrite, ncid)

    ! Get info on the record dimension for this file.
    errstat = nf_inq(ncid, ndims, nvars, ngatts, recdim)

    DO recdim = 1, ndims
      errstat = nf_inq_dim(ncid, recdim, strbuf, nrecs)
      nc_dim(recdim) = nrecs
    ENDDO

    IF (nc_dim(1) < numpcs) THEN
      err = errorstatus_fatal
      THROWM(err.NE.0,"Inconsistent dimensions in eigenvector data: numpcs")
    ENDIF
    IF (nc_dim(2) .NE. numwave) THEN
      err = errorstatus_fatal
      THROWM(err.NE.0,"Inconsistent dimensions in eigenvector data: numwave")
    ENDIF

    ! Read the laboratory eigenvalues into array pcu
    errstat = nf_inq_varid(ncid, 'PC_scores', varid)

    ! Specify integer kinds for compatibility
    errstat = nf_get_vara_real(ncid, varid, INT((/1,1/), KIND=jpim), INT((/numpcs,numwave/), KIND=jpim), pcu)

    errstat = nf_close(ncid)

    CATCH
  END SUBROUTINE rttov_uwiremis_read_labevecs
#else
  SUBROUTINE rttov_uwiremis_read_labevecs(err, fn)

    ! Description:
    ! read the eigenvectors of the 126 selected laboratory spectra
    !(created by E borbas) from the HDF5 file into memory.
    !
    ! History:
    ! Version   Date     Comment
    ! -------   ----     -------
    !  1.1      03/06/2014 original code J. Vidot

    INTEGER(KIND=jpim), INTENT(OUT) :: err
    CHARACTER(LEN=*),   INTENT(IN)  :: fn

    REAL(KIND=jprm),  POINTER    :: pcu_x(:,:)

#include "rttov_hdf_load.interface"

    TRY

    CALL rttov_hdf_load(err, fn, "/EMIS", SNAME="PC_SCORES", PRM2=pcu_x)
    THROWM(ERR.NE.0,"Cannot load PC scores from "//TRIM(fn))

    IF (SIZE(pcu_x,1) < numpcs) THEN
      DEALLOCATE(pcu_x)
      err = errorstatus_fatal
      THROWM(err.NE.0,"Inconsistent dimensions in eigenvector data: numpcs")
    ENDIF
    IF (SIZE(pcu_x,2) .NE. numwave) THEN
      DEALLOCATE(pcu_x)
      err = errorstatus_fatal
      THROWM(err.NE.0,"Inconsistent dimensions in eigenvector data: numwave")
    ENDIF
    pcu = pcu_x(1:numpcs,1:numwave)
    DEALLOCATE(pcu_x)

    CATCH
  END SUBROUTINE rttov_uwiremis_read_labevecs
#endif

#ifndef _RTTOV_HDF
  SUBROUTINE  rttov_uwiremis_read_angfunc(err, fn)

    ! Description:
    ! read the three day-night coefficients for angular correcton
    !(created by E borbas) from the netCDF file into memory.
    !
    ! History:
    ! Version   Date     Comment
    ! -------   ----     -------
    !  1.0    02/04/2014   Eva Borbas

    INTEGER(KIND=jpim), INTENT(OUT) :: err
    CHARACTER(len=*),   INTENT(IN)  :: fn

    INTEGER(KIND=jpim) :: nvars     ! number of variables
    INTEGER(KIND=jpim) :: ndims     ! number of dimensions
    INTEGER(KIND=jpim) :: errstat   ! error code
    INTEGER(KIND=jpim) :: recdim    ! record dimension
    INTEGER(KIND=jpim) :: nc_dim(2)

    INTEGER(KIND=jpim) :: ncid, ngatts, nrecs, varid

    CHARACTER (len=1024) :: strbuf ! string buffer for var

    TRY

   ! Open netCDF file
    errstat = nf_open(TRIM(fn), nf_nowrite, ncid)

    ! Get info on the record dimension for this file.
    errstat = nf_inq(ncid, ndims, nvars, ngatts, recdim)

    DO recdim = 1, ndims
      errstat = nf_inq_dim(ncid, recdim, strbuf, nrecs)
      nc_dim(recdim) = nrecs
    ENDDO

    IF (nc_dim(1) .NE. numwave) THEN
      err = errorstatus_fatal
      THROWM(err.NE.0,"Inconsistent dimensions in angular correction data: numwave")
    ENDIF
    nb_igbp = nc_dim(2)
    ALLOCATE(p1d(numwave,nb_igbp))
    ALLOCATE(p2d(numwave,nb_igbp))
    ALLOCATE(p3d(numwave,nb_igbp))
    ALLOCATE(p1n(numwave,nb_igbp))
    ALLOCATE(p2n(numwave,nb_igbp))
    ALLOCATE(p3n(numwave,nb_igbp))

    ! Retrieve p(3) for day and night
    errstat = nf_inq_varid (ncid, 'p1d', varid)
    errstat = nf_get_var_real(ncid, varid, p1d)
    errstat = nf_inq_varid (ncid, 'p2d', varid)
    errstat = nf_get_var_real(ncid, varid, p2d)
    errstat = nf_inq_varid (ncid, 'p3d', varid)
    errstat = nf_get_var_real(ncid, varid, p3d)
    errstat = nf_inq_varid (ncid, 'p1n', varid)
    errstat = nf_get_var_real(ncid, varid, p1n)
    errstat = nf_inq_varid (ncid, 'p2n', varid)
    errstat = nf_get_var_real(ncid, varid, p2n)
    errstat = nf_inq_varid (ncid, 'p3n', varid)
    errstat = nf_get_var_real(ncid, varid, p3n)

    errstat = nf_close(ncid)

    CATCH
  END SUBROUTINE rttov_uwiremis_read_angfunc
#else
  SUBROUTINE  rttov_uwiremis_read_angfunc(err, fn)

    ! Description:
    ! read the three day-night coefficients for angular correcton
    !(created by E borbas) from the HDF5 file into memory.
    !
    ! History:
    ! Version   Date     Comment
    ! -------   ----     -------
    !  1.0    02/04/2014   Eva Borbas

    INTEGER(KIND=jpim), INTENT(OUT) :: err
    CHARACTER(len=*),   INTENT(IN)  :: fn

#include "rttov_hdf_load.interface"

    TRY

    CALL rttov_hdf_load(err, fn, "/EMIS", SNAME="p1d", PRM2=p1d)
    THROWM(ERR.NE.0,"Cannot load angular correction day 1 coefs data from "//TRIM(fn))

    IF (SIZE(p1d,1) .NE. numwave) THEN
      DEALLOCATE(p1d)
      err = errorstatus_fatal
      THROWM(err.NE.0,"Inconsistent dimensions in angular correction data: numwave")
    ENDIF

    nb_igbp = SIZE(p1d,2)

    CALL rttov_hdf_load(err, fn, "/EMIS", SNAME="p2d", PRM2=p2d)
    THROWM(ERR.NE.0,"Cannot load angular correction day 2 coefs data from "//TRIM(fn))

    CALL rttov_hdf_load(err, fn, "/EMIS", SNAME="p3d", PRM2=p3d)
    THROWM(ERR.NE.0,"Cannot load angular correction day 3 coefs data from "//TRIM(fn))

    CALL rttov_hdf_load(err, fn, "/EMIS", SNAME="p1n", PRM2=p1n)
    THROWM(ERR.NE.0,"Cannot load angular correction night 1 coefs data from "//TRIM(fn))

    CALL rttov_hdf_load(err, fn, "/EMIS", SNAME="p2n", PRM2=p2n)
    THROWM(ERR.NE.0,"Cannot load angular correction night 2 coefs data from "//TRIM(fn))

    CALL rttov_hdf_load(err, fn, "/EMIS", SNAME="p3n", PRM2=p3n)
    THROWM(ERR.NE.0,"Cannot load angular correction night 3 coefs data from "//TRIM(fn))

    CATCH
  END SUBROUTINE rttov_uwiremis_read_angfunc
#endif

#ifndef _RTTOV_HDF
  SUBROUTINE  rttov_uwiremis_read_igbp(err, fn)

    ! Description:
    ! read in the IGBP ecosystem map interpolated for the 0.05x0.05 degree grid and
    ! modified for IR emissivity application
    !(created by E borbas) from the netCDF file into memory.
    !
    ! History:
    ! Version   Date     Comment
    ! -------   ----     -------
    !  1.0    04/04/2014   Eva Borbas

    INTEGER(KIND=jpim), INTENT(OUT) :: err
    CHARACTER(LEN=*),   INTENT(IN)  :: fn

    INTEGER(KIND=jpim) :: nvars     ! number of variables
    INTEGER(KIND=jpim) :: ndims     ! number of dimensions
    INTEGER(KIND=jpim) :: errstat   ! error code
    INTEGER(KIND=jpim) :: recdim    ! record dimension
    INTEGER(KIND=jpim) :: nc_dim(2)

    INTEGER(KIND=jpim) :: ncid, varid, ngatts, nrecs

    CHARACTER(LEN=1024) :: strbuf ! string buffer for var

    TRY

    ! Open netCDF file
    errstat = nf_open(TRIM(fn), nf_nowrite, ncid)

    ! Get info on the record dimension for this file.
    errstat = nf_inq(ncid, ndims, nvars, ngatts, recdim)

    DO recdim = 1, ndims
      errstat = nf_inq_dim(ncid, recdim, strbuf, nrecs)
      nc_dim(recdim) = nrecs
    ENDDO

    igbp_lats = nc_dim(1)
    igbp_lons = nc_dim(2)
    ALLOCATE(igbp(igbp_lats,igbp_lons))

    ! Retrieve igbp
    errstat = nf_inq_varid (ncid, 'IGBP', varid)
    errstat = nf_get_var_int1(ncid, varid, igbp)

    errstat = nf_close(ncid)

    CATCH
  END SUBROUTINE rttov_uwiremis_read_igbp
#else
  SUBROUTINE  rttov_uwiremis_read_igbp(err, fn)

    ! Description:
    ! read in the IGBP ecosystem map interpolated for the 0.05x0.05 degree grid and
    ! modified for IR emissivity application from the HDF5 file into memory.

    INTEGER(KIND=jpim), INTENT(OUT) :: err
    CHARACTER(LEN=*),   INTENT(IN)  :: fn

#include "rttov_hdf_load.interface"

    TRY

    CALL rttov_hdf_load(err, fn, "/EMIS", SNAME="IGBP", PIT2=igbp)
    THROWM(ERR.NE.0,"Cannot load IGBP data from "//TRIM(fn))

    igbp_lats = SIZE(igbp,1)
    igbp_lons = SIZE(igbp,2)

    CATCH
  END SUBROUTINE rttov_uwiremis_read_igbp
#endif

!------------------------------------------
! Routines for returning emissivity values
!------------------------------------------

  SUBROUTINE rttov_uwiremis( &
        & verbose,           &! in
        & nchs,              &! in
        & lat,               &! in
        & lon,               &! in
        & satzen,            &! in
        & solzen,            &! in
        & surfacetype,       &! in
        & snowfrac,          &! in
        & instr_wavenum,     &! in
        & channels,          &! in
        & instr_emis,        &! out
        & instr_emis_cov,    &! out
        & instr_emis_flag)    ! out

    ! Description:
    ! To compute IR emissivty for a given location and frequency
    ! from the 0.1 degree resolution UW BF IR Global Emissivty data
    ! (at 10 hinge points) (http://cimss.ssec.wisc.edu/iremis/)
    ! and labratory measurements using principal component analyses
    !
    ! History:
    ! Version   Date     Comment
    ! -------   ----     -------
    !  0.9       03/31/2009  Original code E Borbas UW-Madison/CIMSS
    !  1.0       03/31/2009  New F90 code with structures (E Borbas B Ruston)

    USE rttov_const, ONLY : surftype_land, surftype_seaice

    LOGICAL(KIND=jplm), INTENT(IN)  :: verbose
    INTEGER(KIND=jpim), INTENT(IN)  :: nchs
    INTEGER(KIND=jpim), INTENT(IN)  :: surfacetype
    INTEGER(KIND=jpim), INTENT(OUT) :: instr_emis_flag

    REAL(KIND=jprb),    INTENT(IN)  :: lat, lon
    REAL(KIND=jprb),    INTENT(IN)  :: satzen, solzen
    REAL(KIND=jprb),    INTENT(IN)  :: snowfrac
    REAL(KIND=jprb),    INTENT(IN)  :: instr_wavenum(nchs)
    INTEGER(KIND=jpim), INTENT(IN)  :: channels(nchs)
    REAL(KIND=jprb),    INTENT(OUT) :: instr_emis(nchs)
    REAL(KIND=jprb),    INTENT(OUT) :: instr_emis_cov(nchs)

    REAL(KIND=jprb) :: coeff(numpcs)

    REAL(KIND=jprb) :: angcorr(numwave)
    REAL(KIND=jprb) :: angcorrchn(nchs,2)
    REAL(KIND=jprb) :: instr_emis_angcorr(nchs,2)

    REAL(KIND=jprb) :: hsremis(numwave)
    REAL(KIND=jprb) :: emis_cov(numwave)

    REAL(KIND=jprb) :: long
    INTEGER(KIND=jpim) :: gridy, gridx, rnd_x, rnd_y, i, ilat, ilon, j

    INTEGER(KIND=jpim) :: igbp_type
    INTEGER(KIND=jpim) :: grid05y, grid05x

    REAL(KIND=jprb) :: hkod = -999._jprb

    REAL(KIND=jprb) :: hsr_ir_emis(nchs)
    REAL(KIND=jprb) :: hsr_ir_emis_cov(nchs)
    REAL(KIND=jprb) :: cov_buff(numwave)

    CHARACTER(LEN=128) :: msg

    instr_emis(:) = hkod
    instr_emis_cov(:) = hkod

    IF (surfacetype == surftype_land) THEN

    !------------------------------------------------------------
    ! find the closest grid point from the uwiremis database
    !------------------------------------------------------------

      long = MODULO(lon, 360.0_jprb)
      IF (long >= 180.0) THEN
        long = long - 360.0_jprb
      ENDIF

      ilat = NINT(lat * 1000._jprb, KIND=jpim)
      ilon = NINT(long * 1000._jprb, KIND=jpim)

      gridy = NINT(ABS(bfemis_ygrid1 - ilat) * 1._jprb / bfemis_gridres, KIND=jpim) + 1_jpim
      gridx = NINT(ABS(bfemis_xgrid1 - ilon) * 1._jprb / bfemis_gridres, KIND=jpim) + 1_jpim

      IF (ir_atlas_do_ang_corr) THEN
        grid05y = NINT(ABS(igbp_ygrid1 - ilat) * 1._jprb / igbp_gridres, KIND=jpim) + 1_jpim
        grid05x = NINT(ABS(igbp_xgrid1 - ilon) * 1._jprb / igbp_gridres, KIND=jpim) + 1_jpim
        igbp_type = igbp(grid05y,grid05x)
      ENDIF

      IF (ir_atlas_std_init) THEN
        rnd_y = NINT(ABS(cov_emis_ygrid1 - ilat) * 1._jprb / cov_emis_gridres, KIND=jpim) + 1_jpim
        rnd_x = NINT(ABS(cov_emis_xgrid1 - ilon) * 1._jprb / cov_emis_gridres, KIND=jpim) + 1_jpim
      ENDIF

      instr_emis_flag = bfemis_flag(gridy,gridx)

    !------------------------------
    ! check if it is a land pixel
    !------------------------------

      IF (instr_emis_flag > 0) THEN

        ! Find the emissivity or coefs and covariances

        IF (bfemis_lut(gridy,gridx) > 0) THEN
          coeff(:) = REAL(pca_coef(bfemis_lut(gridy,gridx),:), KIND=jprb) * pca_sfac(:) + pca_offs(:)
        ELSE
          RETURN
        ENDIF

        IF (ir_atlas_std_init) THEN
          IF (ir_atlas_single_inst) THEN
            IF (cov_emis_lut(rnd_y,rnd_x) > 0) THEN
              ! Note cov_emis_int contains the standard deviations (i.e. sqrt is already taken)
              instr_emis_cov(:) = cov_emis_int(cov_emis_lut(rnd_y,rnd_x),channels(:))
            ELSE
              instr_emis_cov(:) = default_std
            ENDIF
          ELSE
            IF (cov_emis_lut(rnd_y,rnd_x) > 0) THEN
              emis_cov(:) = SQRT(REAL(cov_emis(cov_emis_lut(rnd_y,rnd_x),:), KIND=jprb) * cov_sfac)
            ELSE
              emis_cov(:) = default_std
            ENDIF
          ENDIF
        ENDIF

        IF (ir_atlas_single_inst) THEN
          !--------------------------------------------------------------------------------------------------------
          ! compute the emissivity from the PCs
          !-------------------------------------------------------------------------------------------------------

          IF (ir_atlas_do_ang_corr) THEN

            ! The angular correction is not linear so we must apply it to the HSR
            ! emissivities before they are linearly interpolated to the channel wavenumbers

            ! Reconstruct the two nearest HSR emissivities
            DO j = 1, 2
              CALL rttov_uwiremis_recon_emis( &
                & coeff,                      &
                & channels,                   &
                & instr_emis_angcorr(:,j), j)
            ENDDO

            IF (satzen > angcorrminzen) THEN

              ! Calculate the corresponding angular correction factors
              DO j = 1, 2
                CALL rttov_uwiremis_angcorr( &
                    & p1d_int(:,:,j), p2d_int(:,:,j), p3d_int(:,:,j), &
                    & p1n_int(:,:,j), p2n_int(:,:,j), p3n_int(:,:,j), &
                    & solzen,                                         &
                    & satzen,                                         &
                    & igbp_type,                                      &
                    & angcorrchn(:,j))
              ENDDO

              ! Apply angular correction
              instr_emis_angcorr = instr_emis_angcorr * angcorrchn
            ENDIF

            ! Interpolate to channel wavenumbers (note the reconstructed emissivities
            ! already include the relevant interpolation weights)
            instr_emis = instr_emis_angcorr(:,1) + instr_emis_angcorr(:,2)

          ELSE

            ! With no angular correction the emissivity PCs have been interpolated
            ! to the channel central wavenumbers so we can just reconstruct the
            ! emissivities directly

            CALL rttov_uwiremis_recon_emis( &
              & coeff,                      &! in
              & channels,                   &! in
              & instr_emis)                  ! out

          ENDIF

          !---------------------------------------------------
          ! Linearly blend avg snow emis with snow cover frac
          !---------------------------------------------------
          IF (snowfrac > 0.0_jprb) THEN

            IF (snowfrac > 1.0_jprb) THEN
              instr_emis(:) = snow_em_int(channels(:))
              IF (ir_atlas_std_init) instr_emis_cov(:) = snow_stdv
            ELSE
              instr_emis(:) = snowfrac * snow_em_int(channels(:)) + (1.0_jprb - snowfrac) * instr_emis(:)
              IF (ir_atlas_std_init) instr_emis_cov(:) = &
                    (/ (snowfrac * snow_stdv, i = 1, nchs) /) + (1.0_jprb - snowfrac) * instr_emis_cov(:)
            ENDIF

          ENDIF  ! snow chk

        ELSE

          !--------------------------------------------------------------------------------------------------------
          ! compute the hsr emissivity spectra at 416 wavenumber points from the 10 BF emissivity hinge points
          !-------------------------------------------------------------------------------------------------------

          CALL rttov_uwiremis_recon_hsremis( &
            & coeff,                         &! in
            & hsremis)                        ! out

          !--------------------------------------------------------------------------------------------------------
          ! apply angular correction to the hsr emissivity spectra at 416 wavenumber points
          !-------------------------------------------------------------------------------------------------------

          IF (ir_atlas_do_ang_corr .AND. satzen > angcorrminzen) THEN
            CALL rttov_uwiremis_angcorr( &
                & p1d, p2d, p3d,     &
                & p1n, p2n, p3n,     &
                & solzen,            &! in
                & satzen,            &! in
                & igbp_type,         &! in
                & angcorr)            ! out
            hsremis = hsremis * angcorr
          ENDIF

          !--------------------------------------------------------------------------------
          ! create instrument specific emis/stdv by finidng the closest wavenumber value
          !--------------------------------------------------------------------------------

          CALL rttov_uwiremis_select_wavenum( &
            & hsremis,                        &! in
            & emis_cov,                       &! in
            & nchs,                           &! in
            & instr_wavenum(1:nchs),          &! in
            & instr_emis,                     &! out
            & instr_emis_cov)                  ! out

          !---------------------------------------------------
          ! Linearly blend avg snow emis with snow cover frac
          !---------------------------------------------------
          IF (snowfrac > 0.0_jprb) THEN

            ! snow_stdv is a fixed stdv
            cov_buff(:) = snow_stdv
            CALL rttov_uwiremis_select_wavenum( &
                    & snow_em,                  &! in
                    & cov_buff,                 &! in
                    & nchs,                     &! in
                    & instr_wavenum,            &! in
                    & hsr_ir_emis,              &! out
                    & hsr_ir_emis_cov)           ! out

            IF (snowfrac > 1.0_jprb) THEN
              instr_emis(:) = hsr_ir_emis(:)
              ! Note: a stdv was passed into hsr_ir_emis_cov -- no sqrt
              IF (ir_atlas_std_init) instr_emis_cov(:) = hsr_ir_emis_cov(:)
            ELSE
              instr_emis(:) = snowfrac * hsr_ir_emis(:) + (1._jprb - snowfrac) * instr_emis(:)
              ! Note: a stdv was passed into hsr_ir_emis_cov -- no sqrt
              IF (ir_atlas_std_init) instr_emis_cov(:) = &
                    snowfrac * hsr_ir_emis_cov(:) + (1._jprb - snowfrac) * instr_emis_cov(:)
            ENDIF

          ENDIF  ! snow chk

        ENDIF ! single_inst

      ENDIF  ! emis flag > 0

    ELSE IF (surfacetype == surftype_seaice) THEN

      !---------------------------------------
      ! Return emissivity and stdv for seaice
      !---------------------------------------

      IF (ir_atlas_single_inst) THEN

        instr_emis(:) = sice_em_int(channels(:))
        instr_emis_cov(:) = sice_stdv
        instr_emis_flag = seaice_flag

      ELSE

        ! sice_stdv is a fixed stdv
        cov_buff(:) = sice_stdv
        CALL rttov_uwiremis_select_wavenum( &
                  & sice_em,                &! in
                  & cov_buff,               &! in
                  & nchs,                   &! in
                  & instr_wavenum,          &! out
                  & instr_emis,             &! out
                  & instr_emis_cov        )  ! out

        instr_emis_flag = seaice_flag

      ENDIF

    ELSE
      IF (verbose) THEN
        WRITE (msg, '(a)') 'Warning: IR emissivity atlas should only be called for land and seaice surface types'
        CALL rttov_errorreport(errorstatus_success, msg)
      ENDIF
    ENDIF

    ! Cap the final emissivities here for consistency between single-
    ! and multi-instrument initialisation.
    DO i = 1, nchs
      instr_emis(i) = MIN(instr_emis(i), 1._jprb)
    END DO

  END SUBROUTINE rttov_uwiremis

  SUBROUTINE rttov_uwiremis_hsr_interp(instr_wavenum)

    ! Description:
    ! Initialisation for a single instrument.
    ! Interpolate PC data onto a specific set of wavenumbers:
    ! this can be precomputed for a given instrument during
    ! atlas initialisation to enable very rapid emissivity
    ! calculations.

    REAL(KIND=jprb), INTENT(IN) :: instr_wavenum(:)

    INTEGER(KIND=jpim) :: j, k, nchs

    REAL(KIND=jprb) :: dist(numwave)
    REAL(KIND=jprb) :: mindist
    INTEGER(KIND=jpim) :: ind_mindist

    REAL(KIND=jprb) :: dwvnum1, dwvnum2
    REAL(KIND=jprb) :: pcu1(numpcs), pcu2(numpcs), pcm1, pcm2
    REAL(KIND=jprb) :: sice_em1, sice_em2, snow_em1, snow_em2
    REAL(KIND=jprb) :: cov_emis1(cv_pack), cov_emis2(cv_pack)

    !---------------------------------------------------------------
    ! finding the closest frequency from the hsr emissivity spectr
    !--------------------------------------------------------------

    nchs = SIZE(instr_wavenum)

    DO j = 1, nchs

      ! The angular correction is not linear so in this case we need to
      ! store data for the nearest two HSR wavenumbers so that the
      ! corresponding angcorr values can be applied to the reconstructed
      ! HSR emissivities and then the linearinterpolation to channel
      ! wavenumber can be done.

      IF (instr_wavenum(j) <= hsr_wavenum(1)) THEN

        pcu_int(:,j,1)   = REAL(pcu(:,1), KIND=jprb)
        pcm_int(j,1)     = pcm(1)
        IF (ir_atlas_do_ang_corr) THEN
          pcu_int(:,j,2) = 0._jprb
          pcm_int(j,2)   = 0._jprb
          p1d_int(j,:,1) = p1d(1,:)
          p2d_int(j,:,1) = p2d(1,:)
          p3d_int(j,:,1) = p3d(1,:)
          p1n_int(j,:,1) = p1n(1,:)
          p2n_int(j,:,1) = p2n(1,:)
          p3n_int(j,:,1) = p3n(1,:)
          p1d_int(j,:,2) = 0._jprb
          p2d_int(j,:,2) = 0._jprb
          p3d_int(j,:,2) = 0._jprb
          p1n_int(j,:,2) = 0._jprb
          p2n_int(j,:,2) = 0._jprb
          p3n_int(j,:,2) = 0._jprb
        ENDIF
        sice_em_int(j) = sice_em(1)
        snow_em_int(j) = snow_em(1)
        IF (ir_atlas_std_init) cov_emis_int(:,j) = SQRT(cov_emis(:,1) * cov_sfac)

      ELSE IF (instr_wavenum(j) >= hsr_wavenum(numwave)) THEN

        pcu_int(:,j,1)   = REAL(pcu(:,numwave), KIND=jprb)
        pcm_int(j,1)     = pcm(numwave)
        IF (ir_atlas_do_ang_corr) THEN
          pcu_int(:,j,2) = 0._jprb
          pcm_int(j,2)   = 0._jprb
          p1d_int(j,:,1) = p1d(numwave,:)
          p2d_int(j,:,1) = p2d(numwave,:)
          p3d_int(j,:,1) = p3d(numwave,:)
          p1n_int(j,:,1) = p1n(numwave,:)
          p2n_int(j,:,1) = p2n(numwave,:)
          p3n_int(j,:,1) = p3n(numwave,:)
          p1d_int(j,:,2) = 0._jprb
          p2d_int(j,:,2) = 0._jprb
          p3d_int(j,:,2) = 0._jprb
          p1n_int(j,:,2) = 0._jprb
          p2n_int(j,:,2) = 0._jprb
          p3n_int(j,:,2) = 0._jprb
        ENDIF
        sice_em_int(j) = sice_em(numwave)
        snow_em_int(j) = snow_em(numwave)
        IF (ir_atlas_std_init) cov_emis_int(:,j) = SQRT(cov_emis(:,numwave) * cov_sfac)

      ELSE ! within wavenumber compute range

        mindist = 100._jprb
        ind_mindist = 100000_jpim

        DO k = 1, numwave
          ! calculate distances between the instr freq and hsr emissivities
          dist(k) = ABS(instr_wavenum(j) - hsr_wavenum(k))

          ! finding the closest frequency from the hsr emissivity
          IF (dist(k) <= mindist) THEN
            mindist = dist(k)
            ind_mindist = k
          ENDIF
        ENDDO

  !--------------------------------
  ! Interpolate values
  !--------------------------------

  ! Bilinear mean of the two closest spectral points
        k = 1
        IF (ir_atlas_std_init) THEN
          cov_emis1(:) = 0._jprb
          cov_emis2(:) = 0._jprb
        ENDIF

        IF (instr_wavenum(j) <= hsr_wavenum(ind_mindist)) k = -1

        dwvnum1 = dist(ind_mindist) / (dist(ind_mindist) + dist(ind_mindist + k))
        dwvnum2 = 1._jprb - dwvnum1

        pcu1(:)  = dwvnum1 * REAL(pcu(:,ind_mindist + k), KIND=jprb)
        pcu2(:)  = dwvnum2 * REAL(pcu(:,ind_mindist), KIND=jprb)
        pcm1     = dwvnum1 * pcm(ind_mindist + k)
        pcm2     = dwvnum2 * pcm(ind_mindist)
        sice_em1 = dwvnum1 * sice_em(ind_mindist + k)
        sice_em2 = dwvnum2 * sice_em(ind_mindist)
        snow_em1 = dwvnum1 * snow_em(ind_mindist + k)
        snow_em2 = dwvnum2 * snow_em(ind_mindist)

        IF (ir_atlas_do_ang_corr) THEN
          pcu_int(:,j,1) = pcu1(:)
          pcm_int(j,1)   = pcm1
          pcu_int(:,j,2) = pcu2(:)
          pcm_int(j,2)   = pcm2
          p1d_int(j,:,1) = p1d(ind_mindist + k,:)
          p2d_int(j,:,1) = p2d(ind_mindist + k,:)
          p3d_int(j,:,1) = p3d(ind_mindist + k,:)
          p1n_int(j,:,1) = p1n(ind_mindist + k,:)
          p2n_int(j,:,1) = p2n(ind_mindist + k,:)
          p3n_int(j,:,1) = p3n(ind_mindist + k,:)
          p1d_int(j,:,2) = p1d(ind_mindist,:)
          p2d_int(j,:,2) = p2d(ind_mindist,:)
          p3d_int(j,:,2) = p3d(ind_mindist,:)
          p1n_int(j,:,2) = p1n(ind_mindist,:)
          p2n_int(j,:,2) = p2n(ind_mindist,:)
          p3n_int(j,:,2) = p3n(ind_mindist,:)
        ELSE
          pcu_int(:,j,1) = pcu1(:) + pcu2(:)
          pcm_int(j,1)   = pcm1 + pcm2
        ENDIF
        sice_em_int(j) = sice_em1 + sice_em2
        snow_em_int(j) = snow_em1 + snow_em2

        IF (ir_atlas_std_init) THEN
          cov_emis1(:) = dwvnum1 * SQRT(REAL(cov_emis(:,ind_mindist + k), KIND=jprb) * cov_sfac)
          cov_emis2(:) = dwvnum2 * SQRT(REAL(cov_emis(:,ind_mindist), KIND=jprb) * cov_sfac)
          cov_emis_int(:,j) = cov_emis1(:) + cov_emis2(:)
        ENDIF

      ENDIF

    ENDDO
  END SUBROUTINE rttov_uwiremis_hsr_interp

  SUBROUTINE rttov_uwiremis_recon_emis(coef, channels, emis, ind)

    ! Description:
    ! Used with single-instrument initialisation.
    ! Reconstruct the emissivities at the instrument wavenumbers
    ! from interpolated Principal Components.

    REAL(KIND=jprb),    INTENT(IN)    :: coef(numpcs)
    INTEGER(KIND=jpim), INTENT(IN)    :: channels(:)
    REAL(KIND=jprb),    INTENT(OUT)   :: emis(SIZE(channels))
    INTEGER(KIND=jpim), INTENT(IN), OPTIONAL :: ind

    INTEGER(KIND=jpim) :: j, k, nchn

    !-----------------------------------
    ! apply regcoef to get the emissivities
    !-----------------------------------

    IF (coef(1) .NE. -999._jprb) THEN
      j = 1
      IF (PRESENT(ind)) j = ind
      nchn = SIZE(emis)
      DO k = 1, nchn
        emis(k) = SUM(coef(:) * pcu_int(:,channels(k),j)) + pcm_int(channels(k),j)
      ENDDO
    ELSE
      emis = -999._jprb
    ENDIF

  END SUBROUTINE rttov_uwiremis_recon_emis

  SUBROUTINE rttov_uwiremis_recon_hsremis(coef, hsremis)

    ! Description:
    ! Used with multiple-instrument initialisation.
    ! To creates high spectra resolution emissivties at 416 wavenumbers
    ! from the PCA Coefficitents of the UW BF IR Global Emissivty data
    ! and labratory measurements using principal component analyses
    !
    ! History:
    ! Version   Date     Comment
    ! -------   ----     -------
    !  0.9       03/31/2009  Original code E Borbas UW-Madison/CIMSS
    !  1.0       03/31/2009  New F90 code with structures (E Borbas B Ruston)
    !  1.1       11/30/2012  Removed the coef calcualtion part (E Borbas )

    REAL(KIND=jprb), INTENT(IN)  :: coef(numpcs)
    REAL(KIND=jprb), INTENT(OUT) :: hsremis(numwave)

    INTEGER(KIND=jpim) :: k

    !-----------------------------------
    ! apply regcoef to get the hsr dataset
    !-----------------------------------

    IF (coef(1) .NE. -999._jprb) THEN
      DO k = 1, numwave
        hsremis(k) = SUM(coef(:) * REAL(pcu(:,k), KIND=JPRB)) + pcm(k)
      ENDDO
    ELSE
      hsremis = -999._jprb
    ENDIF

  END SUBROUTINE rttov_uwiremis_recon_hsremis

  SUBROUTINE rttov_uwiremis_select_wavenum ( &
        & hsremis,                           &! in
        & emis_cov,                          &! in
        & nchs,                              &! in
        & instr_wavenum,                     &! in
        & instr_emis,                        &! out
        & instr_emis_cov)                     ! out

    ! Description:
    ! Used with multiple-instrument initialisation.
    ! Subroutine to find the closest wavenumber from the UW HSR emissivity spectra
    ! for the instrument frequency and assign the instrument emissivity by choosing the
    ! closest spectral point value or bilinear interpolating  between the two
    ! closest spectral point values
    !
    ! History:
    ! Version   Date     Comment
    ! -------   ----     -------
    !  0.9       03/31/2009  Original code E Borbas UW-Madison/CIMSS
    !  1.0       03/31/2009  New F90 code with structures (E Borbas B Ruston)

    INTEGER(KIND=jpim), INTENT(IN) :: nchs

    REAL(KIND=jprb), INTENT(IN) :: hsremis(numwave)
    REAL(KIND=jprb), INTENT(IN) :: emis_cov(numwave)
    REAL(KIND=jprb), INTENT(IN) :: instr_wavenum(nchs)

    REAL(KIND=jprb), INTENT(OUT) :: instr_emis(nchs)
    REAL(KIND=jprb), INTENT(OUT) :: instr_emis_cov(nchs)


    INTEGER(KIND=jpim) :: j, k

    REAL(KIND=jprb) :: dist(numwave)
    REAL(KIND=jprb) :: mindist
    INTEGER(KIND=jpim) :: ind_mindist

    REAL(KIND=jprb) :: dwvnum1, dwvnum2, dwvsum
    REAL(KIND=jprb) :: hsremis1, hsremis2, emis_cov1, emis_cov2
    LOGICAL(KIND=jplm) :: lcpu_emis, lcpu_cov

    REAL(KIND=jprb) :: hkod = -999._jprb

    ! initialize instrument emissivity

    !---------------------------------------------------------------
    ! finding the closest frequency from the hsr emissivity spectra
    !---------------------------------------------------------------
    lcpu_emis = .TRUE.
    lcpu_cov  = .TRUE.
    IF (ALL(hsremis == hsremis(1))) lcpu_emis = .FALSE.
    IF (ir_atlas_std_init) THEN
      IF (ALL(emis_cov == emis_cov(1))) lcpu_cov  = .FALSE.
    ELSE
      lcpu_cov  = .FALSE.
    ENDIF

    IF (lcpu_emis .or. lcpu_cov) THEN
      instr_emis(:)       = hkod
      instr_emis_cov(:)   = hkod
      DO j = 1, nchs

        IF (instr_wavenum(j) <= hsr_wavenum(1)) THEN

          instr_emis(j) = hsremis(1)
          IF (ir_atlas_std_init) instr_emis_cov(j) = emis_cov(1)

        ELSEif(instr_wavenum(j) >= hsr_wavenum(numwave)) THEN

          instr_emis(j) = hsremis(numwave)
          IF (ir_atlas_std_init) instr_emis_cov(j) = emis_cov(numwave)

        ELSE ! within wavenumber compute range

          mindist = 100._jprb
          ind_mindist = 100000_jpim

          DO k = 1, numwave

            ! calculate distances between the instr freq end hsr emissivities
            dist(k) = ABS(instr_wavenum(j) - hsr_wavenum(k))

            ! finding the closest frequency from the hsr emissivity
            IF(dist(k) <=  mindist) THEN
              mindist = dist(k)
              ind_mindist = k
            ENDIF
          ENDDO

    !--------------------------------
    ! assign instrument emissivity
    !--------------------------------
    !  closest spectral point
    !                       instr_emis(j)=hsremis(ind_mindist)
    !                       instr_emis_cov(j)=cov_emis(ind_mindist)

    ! or bilinear mean of the two closest spectral points

          k = 1
          IF (instr_wavenum(j) <= hsr_wavenum(ind_mindist)) k = -1

          dwvnum1 = dist(ind_mindist)
          dwvnum2 = dist(ind_mindist + k)
          dwvsum = dwvnum1 + dwvnum2

          IF (lcpu_emis) THEN
            hsremis1 = dwvnum1 * hsremis(ind_mindist + k)
            hsremis2 = dwvnum2 * hsremis(ind_mindist)
            instr_emis(j) = (hsremis1 + hsremis2) / dwvsum
          ELSE
            instr_emis(j) = hsremis(1)
          ENDIF

          IF (ir_atlas_std_init) THEN
            IF (lcpu_cov) THEN
              emis_cov1 = dwvnum1 * emis_cov(ind_mindist + k)
              emis_cov2 = dwvnum2 * emis_cov(ind_mindist)
              instr_emis_cov(j) = (emis_cov1 + emis_cov2) / dwvsum
            ELSE
              instr_emis_cov(j) = emis_cov(1)
            ENDIF
          ENDIF

        ENDIF    !==  (instr_wavenum(j) <= hsr_wavenum(1))

      ENDDO

    ELSE  ! all logical computes (lcpu_emis, lcpu_cov) are false
      instr_emis(:) = hsremis(1)
      IF (ir_atlas_std_init) instr_emis_cov(:) = emis_cov(1)
    ENDIF

  END SUBROUTINE rttov_uwiremis_select_wavenum

  SUBROUTINE rttov_uwiremis_angcorr(  &
      & p1d, p2d, p3d, p1n, p2n, p3n, &
      & solzen,                       &! in
      & satzen,                       &! in
      & igbp_type,                    &! in
      & angcorr)                       ! out

    ! Description:
    ! Subroutine to calculate the satellite zenith angle correction
    !
    ! History:
    ! Version   Date     Comment
    ! -------   ----     -------
    !  0.9       01/22/2014  Original code E Borbas UW-Madison/CIMSS

    REAL(KIND=jprm),    INTENT(IN)  :: p1d(:,:), p2d(:,:), p3d(:,:)
    REAL(KIND=jprm),    INTENT(IN)  :: p1n(:,:), p2n(:,:), p3n(:,:)
    REAL(KIND=jprb),    INTENT(IN)  :: satzen, solzen
    INTEGER(KIND=jpim), INTENT(IN)  :: igbp_type
    REAL(KIND=jprb),    INTENT(OUT) :: angcorr(:)

    INTEGER(KIND=jpim) :: indclst
    REAL(KIND=jprb)    :: pv(SIZE(angcorr)), dpv(SIZE(angcorr)), p1(SIZE(angcorr))
!     CHARACTER(LEN=128) :: msg

    IF (igbp_type > 0) THEN

      indclst = NINT(ABS(satzen))

      ! Note from E Borbas:
      ! We have the nonlinear fitting function (with three coefficients) and then we
      ! calculate the differences of the values between the nadir point (p3) and at
      ! the giving viewing angle. In the future we may develop the method further
      ! and p1 may not be equal to p3, but also added some bias.
      ! Therefore current code (where p3 is not actually used) is kept for clarity.

      IF (solzen <= angcorrterm) THEN
        pv(:) = p1d(:,igbp_type) * indclst**2 + p2d(:,igbp_type) * indclst + p3d(:,igbp_type)
        p1(:) = p3d(:,igbp_type)
      ELSE
        pv(:) = p1n(:,igbp_type) * indclst**2 + p2n(:,igbp_type) * indclst + p3n(:,igbp_type)
        p1(:) = p3n(:,igbp_type)
      ENDIF

      dpv(:) = p1(:) - pv(:)
      angcorr(:) = (1._jprb - dpv(:))

    ELSE
!       IF (verbose) &
!         WRITE (msg,'(a)') &
!           'Warning: IGBP type = ocean water or coast line: no IR emissivity angular correction is applied'
!         CALL rttov_errorreport(errorstatus_success, msg)
      angcorr(:) = 1._jprb
    ENDIF

  END SUBROUTINE rttov_uwiremis_angcorr

!------------------------------------
! Routine to deallocate atlas arrays
!------------------------------------
  SUBROUTINE rttov_uwiremis_close_atlas
    IF ( ASSOCIATED(bfemis_flag)  ) DEALLOCATE(bfemis_flag)
    IF ( ASSOCIATED(bfemis_lut)   ) DEALLOCATE(bfemis_lut)
    IF ( ASSOCIATED(pca_coef)     ) DEALLOCATE(pca_coef)
    IF ( ASSOCIATED(cov_emis_lut) ) DEALLOCATE(cov_emis_lut)
    IF ( ASSOCIATED(cov_emis)     ) DEALLOCATE(cov_emis)
    IF ( ASSOCIATED(pcu_int)      ) DEALLOCATE(pcu_int)
    IF ( ASSOCIATED(pcm_int)      ) DEALLOCATE(pcm_int)
    IF ( ASSOCIATED(sice_em_int)  ) DEALLOCATE(sice_em_int)
    IF ( ASSOCIATED(snow_em_int)  ) DEALLOCATE(snow_em_int)
    IF ( ASSOCIATED(cov_emis_int) ) DEALLOCATE(cov_emis_int)

    IF ( ASSOCIATED(igbp)      ) DEALLOCATE(igbp)
    IF ( ASSOCIATED(p1d)       ) DEALLOCATE(p1d)
    IF ( ASSOCIATED(p2d)       ) DEALLOCATE(p2d)
    IF ( ASSOCIATED(p3d)       ) DEALLOCATE(p3d)
    IF ( ASSOCIATED(p1n)       ) DEALLOCATE(p1n)
    IF ( ASSOCIATED(p2n)       ) DEALLOCATE(p2n)
    IF ( ASSOCIATED(p3n)       ) DEALLOCATE(p3n)
    IF ( ASSOCIATED(p1d_int)   ) DEALLOCATE(p1d_int)
    IF ( ASSOCIATED(p2d_int)   ) DEALLOCATE(p2d_int)
    IF ( ASSOCIATED(p3d_int)   ) DEALLOCATE(p3d_int)
    IF ( ASSOCIATED(p1n_int)   ) DEALLOCATE(p1n_int)
    IF ( ASSOCIATED(p2n_int)   ) DEALLOCATE(p2n_int)
    IF ( ASSOCIATED(p3n_int)   ) DEALLOCATE(p3n_int)

#ifdef _RTTOV_HDF
    IF ( ASSOCIATED(pca_offs)   ) DEALLOCATE(pca_offs)
    IF ( ASSOCIATED(pca_sfac)   ) DEALLOCATE(pca_sfac)
#endif

    CALL rttov_uwiremis_nullify_pointers()

  END SUBROUTINE

  SUBROUTINE rttov_uwiremis_nullify_pointers

    NULLIFY(bfemis_flag)
    NULLIFY(bfemis_lut)
    NULLIFY(pca_coef)
    NULLIFY(cov_emis_lut)
    NULLIFY(cov_emis)
    NULLIFY(pcu_int)
    NULLIFY(pcm_int)
    NULLIFY(sice_em_int)
    NULLIFY(snow_em_int)
    NULLIFY(cov_emis_int)

    NULLIFY(igbp)
    NULLIFY(p1d)
    NULLIFY(p2d)
    NULLIFY(p3d)
    NULLIFY(p1n)
    NULLIFY(p2n)
    NULLIFY(p3n)
    NULLIFY(p1d_int)
    NULLIFY(p2d_int)
    NULLIFY(p3d_int)
    NULLIFY(p1n_int)
    NULLIFY(p2n_int)
    NULLIFY(p3n_int)

#ifdef _RTTOV_HDF
    NULLIFY(pca_offs)
    NULLIFY(pca_sfac)
#endif

  END SUBROUTINE

END MODULE mod_iratlas
