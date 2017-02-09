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

///@class Profiles
///  This class represents multiple atmospheric profiles

#ifndef PROFILES_H_
#define PROFILES_H_


#include <rttov_cc_interface.h>
#include <stdio.h>

#include <iostream>

#include <cstdlib>

#include <exception>
#include <stdexcept>
#include <Rttov_common.h>

namespace rttov {

typedef std::map<rttov::itemIdType,double *> GasIdPointerMap;

class Profiles {
public:
    Profiles(int nbprofiles, const int nblevels, int nbgases);
    Profiles(int nbprofiles, const int nblevels);
    Profiles();
    virtual ~Profiles();

//     void setNprofiles(int nprofiles);
//     void setNlevels(int nlevels);
//     void setNgases(int ngases);

    void setGasUnits(int gasUnits);

    void setP(double* p);
    void setT(double* t);
    void setQ(double* q);
    void setO3(double* o3);
    void setCO2(double* co2);
    void setCO(double* co);
    void setN2O(double* n2o);
    void setCH4(double* ch4);
    void setCLW(double* clw);

    void setAngles(double* angles);
    void setS2m(double* s2m);
    void setSkin(double* skin);
    void setSurfType(int* surftype);
    void setSurfGeom(double* surfgeom);
    void setDateTimes(int* datetimes);
    void setSimpleCloud(double* simplecloud);
    void setIceCloud(int* icecloud);
    void setZeeman(double* zeeman);
    void setGasItem(double* gasItem, rttov::itemIdType item_id);

    void setGasId(int* gasId);
    void setGases(double* gases);
    ItemIdIndexMap gas_index;

    int getNprofiles() const;
    int getNlevels() const;
    int getGasUnits() const;

    double* getP() const;
    double* getT() const;

    double* getQ();
    double* getO3();
    double* getCO2();
    double* getCO();
    double* getN2O();
    double* getCH4();
    double* getCLW();

    double* getAngles() const;
    double* getS2m() const;
    double* getSkin() const;
    int* getSurfType() const;
    double* getSurfGeom() const;
    int* getDateTimes() const;
    double* getSimpleCloud() const;
    int* getIceCloud() const;
    double* getZeeman() const;

    int getNgases() const;
    int* getGasId() const;
    double* getGases() const;

    bool isDefaultPressureLevels() const;
    bool check();

private :

    int nprofiles;
    int nlevels;
    int ngases;
    int gas_units;
    double * p;//[nprofiles][nlevels];                  // Input pressure profiles
    double * t;//[nprofiles][nlevels];                  // Input temperature profiles
    double * gases;//[ngases][nprofiles][nlevels];      // Input gas profiles
    bool allocatedGases;

    int * gas_id;// { q=1, o3, co2, n2o, co,  ch4, cfrac=20, lwc1=21, iwc=30}

    // datetimes: yy, mm, dd, hh, mm, ss
    int * datetimes; //[nprofiles][6]

    // angles: satzen, satazi, sunzen, sunazi
    double * angles;//[nprofiles][4]

    // surftype: surftype, watertype
    int * surftype;//[nprofiles][2]

    // surfgeom: lat, lon, elev
    double * surfgeom ;//[nprofiles][3]

    // s2m: 2m p, 2m t, 2m q, 10m wind u, v, wind fetch
    double * s2m;//[nprofiles][6]

    // skin: skin T, salinity, snow_frac, foam_frac, fastem_coefsx5
    double * skin; //[nprofiles][9]

    // simplecloud: ctp, cfraction
    double * simplecloud ;//[nprofiles][2]

    // icecloud: ish, idg
    int * icecloud ;//[nprofiles][2]

    // zeeman: be, cosbk
    double * zeeman;//[nprofiles][2]

    rttov::StrBoolMap allocatedPointers;
    void initialisePointers();
    bool buildGasesArray();
    GasIdPointerMap myGases;
};
} /* namespace rttov */
#endif /* PROFILES_H_ */
