/*
! Copyright:
!    This software was developed within the context of
!    the EUMETSAT Satellite Application Facility on
!    Numerical Weather Prediction (NWP SAF), under the
!    Cooperation Agreement dated 25 November 1998, between
!    EUMETSAT and the Met Office, UK, by one or more partners
!    within the NWP SAF. The partners in the NWP SAF are
!    the Met Office, ECMWF, KNMI and MeteoFrance.
!
!    Copyright 2015, EUMETSAT, All Rights Reserved.
*/

#ifndef RTTOV_COMMON_H_
#define RTTOV_COMMON_H_
#include <string>
#include <map>
#include <vector>

namespace rttov {
enum gasUnitType {unknown=-2, ppmv_dry=-1, compatibility_mode=0, kg_per_kg=1, ppmv_wet=2};

enum itemIdType {Q=1, O3, CO2, N2O, CO, CH4, CLW=15, CFRAC=20, STCO, STMA, CUCC, CUCP, CUMA, CIRR=30, ICEDE,
    INSO=41, WASO, SOOT, SSAM, SSCM, MINM, MIAM, MICM, MITR, SUSO, VOLA, VAPO, ASDU};

typedef std::vector <enum itemIdType> itemIdVector;
typedef std::map<itemIdType, std::vector <double>> ItemIdPointerMap;
typedef std::map<itemIdType, int> ItemIdIndexMap;
typedef std::map<std::string, bool> StrBoolMap;

const itemIdVector itemIds {Q, O3, CO2, N2O, CO, CH4, CLW, CFRAC, STCO, STMA, CUCC, CUCP, CUMA, CIRR, ICEDE,
    INSO, WASO, SOOT, SSAM, SSCM, MINM, MIAM, MICM, MITR, SUSO, VOLA, VAPO, ASDU};
}

#endif /* RTTOV_COMMON_H_ */
