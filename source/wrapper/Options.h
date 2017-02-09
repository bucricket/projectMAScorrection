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

///@class Options
///@brief This class allows the RTTOV and wrapper options to be specified

#ifndef OPTIONS_H_
#define OPTIONS_H_
#include <map>
#include <string>

namespace rttov {

typedef std::map<std::string,bool>StrBoolMap;
enum ipcbndType {band1=1,band2,band3};
enum ipcregType {index1=1,index2,index3,index4};
enum fastemType {fastem1=1,fastem2,fastem3,fastem4,fastem5,fastem6};

class Options {
public:
    Options();
    virtual ~Options();

    // RTTOV opts%config%
    void setApplyRegLimits(bool applyRegLimts);
    void setDoCheckinput(bool doCheckinput);
    void setVerbose(bool verbose);

    //RTTOV opts%interpolation%
    void setAddInterp(bool addinterp);
    void setInterpMode(int interpMode);
    void setRegLimitExtrap(bool regLimitExtrap);
    void setSpacetop(bool spacetop);
    void setLgradp(bool lgradp);

    //RTTOV opts%rt_all%
    void setDoLambertian(bool doLambertian);
    void setUseQ2m(bool useQ2m);
    void setSwitchrad(bool switchrad);
    void setAddRefrac(bool addRefrac);

    //RTTOV opts%rt_mw%
    void setCLWData(bool clwData);
    void setFastemVersion(int fastemVersion);
    void setSupplyFoamFraction(bool supplyFoamFraction);

    //RTTOV opts%rt_ir%
    void setOzoneData(bool ozoneData);
    void setCO2Data(bool co2Data);
    void setCH4Data(bool ch4Data);
    void setCOData(bool coData);
    void setN2OData(bool n2oData);

    void setAddSolar(bool addsolarl);
    void setDoNlteCorrection(bool doNlteCorrection);

    void setAddAerosl(bool addaerosl);
    void setAddClouds(bool addclouds);
    void setCldstrSimple(bool cldstrCimple);
    void setCldstrThreshold(double cldstrThreshold);

//     void setUserAerOptParam(bool userAerOptParam);
//     void setUserCldOptParam(bool userCldOptParam);

//     void setAddPC(bool addpc);
//     void setAddRadrec(bool addradrec);
//     void setIpcreg(int ipcreg);
//     void setIpcbnd(int ipcbnd);

    // Wrapper options
    void setNthreads(int nthreads);
    void setNprofsPerCall(int nprofsPerCall);
    void setVerboseWrapper(bool verboseWrapper);
    void setCheckOpts(bool checkOpts);
    void setStoreRad(bool storeRad);
    void setStoreRad2(bool storeRad2);
    void setStoreTrans(bool storeTrans);


    // RTTOV opts%config%
    bool isApplyRegLimits();
    bool isDoCheckinput();
    bool isVerbose();

    //RTTOV opts%interpolation%
    bool isAddInterp();
    int getInterpMode() const;
    bool isRegLimitExtrap();
    bool isSpacetop();
    bool isLgradp();

    //RTTOV opts%rt_all%
    bool isDoLambertian();
    bool isUseQ2m();
    bool isSwitchrad();
    bool isAddRefrac();

    //RTTOV opts%rt_mw%
    bool isCLWData();
    int getFastemVersion() const;
    bool isSupplyFoamFraction();

    //RTTOV opts%rt_ir%
    bool isOzoneData();
    bool isCO2Data();
    bool isCH4Data();
    bool isCOData();
    bool isN2OData();

    bool isAddSolar();
    bool isDoNlteCorrection();

    bool isAddAerosl();
    bool isAddClouds();
    bool isCldstrSimple();
    double getCldstrThreshold() const;

    bool isUserAerOptParam();
    bool isUserCldOptParam();

    bool isAddPC();
    bool isAddRadrec();
    int getIpcreg() const;
    int getIpcbnd() const;

    // Wrapper options
    int getNthreads() const;
    int getNprofsPerCall() const;
    bool isVerboseWrapper() const;
    bool isCheckOpts() const;
    bool isStoreRad() const;
    bool isStoreRad2() const;
    bool isStoreTrans() const;

    std::string defineStrOptions();

private :

    StrBoolMap config;
    StrBoolMap ir;
    StrBoolMap mw;
    StrBoolMap pc;
    StrBoolMap interp;
    StrBoolMap rt_all;
    bool verbose_wrapper;
    bool check_opts;
    bool store_trans;
    bool store_rad;
    bool store_rad2;

    int nthreads;
    int nprofs_per_call;
    int cldstr_threshold;
    int fastem_version;
    int interp_mode;
    int ipcbnd;
    int ipcreg;
};
} /* namespace rttov */
#endif /* OPTIONS_H_ */
