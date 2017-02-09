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

#include <Options.h>
#include <sstream>
/**
 * @file Options.cpp
 * Options class definitions
 *
 * @brief Options class definitions for managing rttov options.
 *
 * @author Pascale

 * @version
 */

namespace rttov {

///@brief Constructor method
Options::Options(): nprofs_per_call(1),
                    nthreads(1),
                    cldstr_threshold(-1.),
                    fastem_version(fastem5),
                    interp_mode(1),
                    ipcbnd(band1),
                    ipcreg(index1){
    // General configuration options
    this->config["apply_reg_limits"]=false;
    this->config["verbose"]=true;
    this->config["do_checkinput"]=true;
    // Interpolation options
    this->interp["addinterp"]=true ; // to prevent crashing if user forget to set it if profile number of level different from coeffs
    this->interp["spacetop"]=true;
    this->interp["lgradp"]=false;
    this->interp["reg_limit_extrap"]=false;
    // General RT options
    this->rt_all["addrefrac"]=false;
    this->rt_all["switchrad"]=false;
    this->rt_all["use_q2m"]=true;
    this->rt_all["addrefrac"]=false;
    this->rt_all["do_lambertian"]=false;
    // MW-only RT options
    this->mw["clw_data"]=false;
    this->mw["supply_foam_fraction"]=false;
    // VIS/IR-only RT options
    this->ir["addsolar"]=false;
    this->ir["do_nlte_correction"]=false;
    this->ir["addaerosl"]=false;
    this->ir["addclouds"]=false;
    this->ir["cldstr_simple"]=false;
    this->ir["user_aer_opt_param"]=false;
    this->ir["user_cld_opt_param"]=false;
    this->ir["ozone_data"]=false;
    this->ir["co2_data"]=false;
    this->ir["n2o_data"]=false;
    this->ir["co_data"]=false;
    this->ir["ch4_data"]=false;
    // PC-RTTOV options
    this->pc["addpc"]=false;
//     this->pc["addradrec"]=false;
    // Wrapper options
    this->verbose_wrapper=false;
    this->check_opts=true;
    this->store_rad=false;
    this->store_rad2=false;
    this->store_trans=false;
}

///@brief Destructor method
Options::~Options() {

}

// RTTOV opts%config%

///@brief Return the opts\%config\%apply_reg_limits option
bool Options::isApplyRegLimits() {
    return this->config["apply_reg_limits"] ;
}

///@brief Set the opts\%config\%apply_reg_limits option
void Options::setApplyRegLimits(bool applyRegLimits) {
    this->config["apply_reg_limits"] = applyRegLimits;
}

///@brief Return the opts\%config\%do_checkinput option
bool Options::isDoCheckinput() {
    return this->config["do_checkinput"];
}

///@brief Set the opts\%config\%do_checkinput option
void Options::setDoCheckinput(bool doCheckinput) {
    this->config["do_checkinput"] = doCheckinput;
}

///@brief Return the opts\%config\%verbose option
bool Options::isVerbose() {
    return this->config["verbose"];
}

///@brief Set the opts\%config\%verbose option
void Options::setVerbose(bool verbose) {
    this->config["verbose"] = verbose;
}


// RTTOV opts%interpolation%

///@brief Return the opts\%interpolation\%addinterp option
bool Options::isAddInterp() {
    return this->interp["addinterp"];
}

///@brief Set the opts\%interpolation\%addinterp option
void Options::setAddInterp(bool addinterp) {
    this->interp["addinterp"]= addinterp;
}

///@brief Return the opts\%interpolation\%interp_mode option
int Options::getInterpMode() const {
    return interp_mode;
}

///@brief Set the opts\%interpolation\%interp_mode option
void Options::setInterpMode(int interpMode) {
    interp_mode = interpMode;
}

///@brief Return the opts\%interpolation\%reg_limit_extrap option
bool Options::isRegLimitExtrap() {
    return this->interp["reg_limit_extrap"];
}

///@brief Set the opts\%interpolation\%reg_limit_extrap option
void Options::setRegLimitExtrap(bool regLimitExtrap) {
    this->interp["reg_limit_extrap"] = regLimitExtrap;
}

///@brief Return the opts\%interpolation\%spacetop option
bool Options::isSpacetop() {
    return this->interp["spacetop"];
}

///@brief Set the opts\%interpolation\%spacetop option
void Options::setSpacetop(bool spacetop) {
    this->interp["spacetop"] = spacetop;
}

///@brief Return the opts\%interpolation\%lgradp option
bool Options::isLgradp() {
    return this->interp["lgradp"];
}

///@brief Set the opts\%interpolation\%lgradp option
void Options::setLgradp(bool lgradp) {
    this->interp["lgradp"]= lgradp;
}


// RTTOV opts%rt_all%

///@brief Return the opts\%rt_all\%switchrad option
bool Options::isSwitchrad() {
    return this->rt_all["switchrad"];
}

///@brief Set the opts\%rt_all\%switchrad option
void Options::setSwitchrad(bool switchrad) {
    this->rt_all["switchrad"] = switchrad;
}

///@brief Return the opts\%rt_all\%use_q2m option
bool Options::isUseQ2m() {
    return this->rt_all["use_q2m"];
}

///@brief Set the opts\%rt_all\%use_q2m option
void Options::setUseQ2m(bool useQ2m) {
    this->rt_all["use_q2m"] = useQ2m;
}

///@brief Return the opts\%rt_all\%addrefrac option
bool Options::isAddRefrac() {
    return this->rt_all["addrefrac"];
}

///@brief Set the opts\%rt_all\%addrefrac option
void Options::setAddRefrac(bool addRefrac) {
    this->rt_all["addrefrac"] = addRefrac;
}

///@brief Return the opts\%rt_all\%do_lambertian option
bool Options::isDoLambertian() {
    return this->rt_all["do_lambertian"];
}

///@brief Set the opts\%rt_all\%do_lambertian option
void Options::setDoLambertian(bool doLambertian) {
    this->rt_all["do_lambertian"] = doLambertian;
}


// RTTOV opts%rt_mw%

///@brief Return the opts\%rt_mw\%supply_foam_fraction option
bool Options::isSupplyFoamFraction() {
    return this->mw["supply_foam_fraction"];
}

///@brief Set the opts\%rt_mw\%supply_foam_fraction option
void Options::setSupplyFoamFraction(bool supplyFoamFraction) {
    this->mw["supply_foam_fraction"] = supplyFoamFraction;
}

///@brief Return the opts\%rt_mw\%clw_data option
bool Options::isCLWData() {
    return this->mw["clw_data"];
}

///@brief Set the opts\%rt_mw\%clw_data option
void Options::setCLWData(bool clwData) {
    this->mw["clw_data"] = clwData;
}

///@brief Return the opts\%rt_mw\%fastem_version option
int Options::getFastemVersion() const {
    return fastem_version;
}

///@brief Set the opts\%rt_mw\%fastem_version option
void Options::setFastemVersion(int fastemVersion) {
    fastem_version = fastemVersion;
}


// RTTOV opts%rt_ir%

///@brief Return the opts\%rt_ir\%ozone_data option
bool Options::isOzoneData() {
    return this->ir["ozone_data"];
}

///@brief Set  the opts\%rt_ir\%ozone_data option
void Options::setOzoneData(bool ozoneData) {
    this->ir["ozone_data"] = ozoneData;
}

///@brief Return the opts\%rt_ir\%co2_data option
bool Options::isCO2Data() {
    return this->ir["co2_data"];
}

///@brief Set the opts\%rt_ir\%co2_data option
void Options::setCO2Data(bool co2Data) {
    this->ir["co2_data"] = co2Data;
}

///@brief Return the opts\%rt_ir\%ch4_data option
bool Options::isCH4Data() {
    return this->ir["ch4_data"];
}

///@brief Set the opts\%rt_ir\%ch4_data option
void Options::setCH4Data(bool ch4Data) {
    this->ir["ch4_data"] = ch4Data;
}

///@brief Return the opts\%rt_ir\%co_data option
bool Options::isCOData() {
    return this->ir["co_data"];
}

///@brief Set the opts\%rt_ir\%co_data option
void Options::setCOData(bool coData) {
    this->ir["co_data"] = coData;
}

///@brief Return the opts\%rt_ir\%n2o_data option
bool Options::isN2OData() {
    return this->ir["n2o_data"];
}

///@brief Set the opts\%rt_ir\%n2o_data option
void Options::setN2OData(bool n2oData) {
    this->ir["n2o_data"] = n2oData;
}

///@brief Return the opts\%rt_ir\%addsolar option
bool Options::isAddSolar() {
    return this->ir["addsolar"];
}

///@brief Set the opts\%rt_ir\%addsolar option
void Options::setAddSolar(bool addsolar) {
    this->ir["addsolar"] = addsolar;
}

///@brief Return the opts\%rt_ir\%do_nlte_correction option
bool Options::isDoNlteCorrection() {
    return this->ir["do_nlte_correction"];
}

///@brief Set the opts\%rt_ir\%do_nlte_correction option
void Options::setDoNlteCorrection(bool doNlteCorrection) {
    this->ir["do_nlte_correction"] = doNlteCorrection;
}

///@brief Return the opts\%rt_ir\%addaerosl option
bool Options::isAddAerosl() {
    return this->ir["addaerosl"];
}

///@brief Set the opts\%rt_ir\%addaerosl option
void Options::setAddAerosl(bool addaerosl) {
    this->ir["addaerosl"] = addaerosl;
}

///@brief Return the opts\%rt_ir\%addclouds option
bool Options::isAddClouds() {
    return this->ir["addclouds"];
}

///@brief Set the opts\%rt_ir\%addclouds option
void Options::setAddClouds(bool addclouds) {
    this->ir["addclouds"]= addclouds;
}

///@brief Return the opts\%rt_ir\%user_aer_opt_param option
bool Options::isUserAerOptParam() {
    return this->ir["user_aer_opt_param"];
}

/////@brief Set the opts\%rt_ir\%user_aer_opt_param option
// void Options::setUserAerOptParam(bool userAerOptParam) {
//     this->ir["user_aer_opt_param"] = userAerOptParam;
// }

///@brief Return the opts\%rt_ir\%user_cld_opt_param option
bool Options::isUserCldOptParam() {
    return this->ir["user_cld_opt_param"];
}

/////@brief Set the opts\%rt_ir\%user_cld_opt_param option
// void Options::setUserCldOptParam(bool userCldOptParam) {
//     this->ir["user_cld_opt_param"] = userCldOptParam;
// }

///@brief Return the opts\%rt_mw\%cldstr_threshold option
double Options::getCldstrThreshold() const {
    return cldstr_threshold;
}

///@brief Set the opts\%rt_mw\%cldstr_threshold option
void Options::setCldstrThreshold(double cldstrThreshold) {
    cldstr_threshold = cldstrThreshold;
}

///@brief Return the opts\%rt_ir\%cldstr_simple option
bool Options::isCldstrSimple() {
    return this->ir["cldstr_simple"];
}

///@brief Set the opts\%rt_ir\%cldstr_simple option
void Options::setCldstrSimple(bool cldstrSimple) {
    this->ir["cldstr_simple"] = cldstrSimple;
}


// RTTOV opts%rt_ir%pc%

///@brief Return the opts\%rt_ir\%pc\%addpc option
bool Options::isAddPC() {
    return this->pc["addpc"];
}

/////@brief Set the opts\%rt_ir\%pc\%addpc option
// void Options::setAddPC(bool addpc) {
//     this->pc["addpc"]= addpc;
// }

///@brief Return the opts\%rt_ir\%pc\%addradrec option
bool Options::isAddRadrec() {
    return this->pc["addradrec"];
}

/////@brief Set the opts\%rt_ir\%pc\%addradrec option
// void Options::setAddRadrec(bool addradrec) {
//     this->pc["addradrec"]= addradrec;
// }

///@brief Return the opts\%rt_ir\%pc\%ipcreg option
int Options::getIpcreg() const {
    return ipcreg;
}

/////@brief Set the opts\%rt_ir\%pc\%ipcreg option
// void Options::setIpcreg(int ipcregx) {
//     ipcreg = ipcregx;
// }

///@brief Return the opts\%rt_ir\%pc\%ipcbnd option
int Options::getIpcbnd() const {
    return ipcbnd;
}

/////@brief Set the opts\%rt_ir\%pc\%ipcbnd option
// void Options::setIpcbnd(int ipcbndx) {
//     ipcbnd = ipcbndx;
// }


// RTTOV wrapper options

///@brief Return the number of profiles passed into rttov_direct or rttov_k per call
int Options::getNprofsPerCall() const {
    return nprofs_per_call;
}

///@brief Set the number of profiles passed into rttov_direct or rttov_k per call
void Options::setNprofsPerCall(int nprofsPerCall) {
    nprofs_per_call = nprofsPerCall;
}

///@brief Return the number of threads RTTOV will use (compile RTTOV with OpenMP to make use of this)
int Options::getNthreads() const {
    return nthreads;
}

///@brief Set the number of threads RTTOV will use (compile RTTOV with OpenMP to make use of this)
void Options::setNthreads(int nthreads) {
    this->nthreads = nthreads;
}

///@brief Return set the verbose_wrapper option
bool Options::isVerboseWrapper() const {
    return verbose_wrapper;
}

///@brief Set the verbose_wrapper option
void Options::setVerboseWrapper(bool verboseWrapper) {
    verbose_wrapper = verboseWrapper;
}

///@brief Return set the check_opts option
bool Options::isCheckOpts() const {
    return check_opts;
}

///@brief Set the check_opts option
void Options::setCheckOpts(bool checkOpts) {
    check_opts = checkOpts;
}

///@brief Return the store_rad wrapper option
bool Options::isStoreRad() const {
    return store_rad;
}

///@brief Set the store_rad wrapper option
void Options::setStoreRad(bool storeRad) {
    store_rad = storeRad;
}

///@brief Return the store_rad2 wrapper option
bool Options::isStoreRad2() const {
    return store_rad2;
}

///@brief Set the store_rad2 wrapper option
void Options::setStoreRad2(bool storeRad2) {
    store_rad2 = storeRad2;
}

///@brief Return the store_trans wrapper option
bool Options::isStoreTrans() const {
    return store_trans;
}

///@brief Set the store_trans wrapper option
void Options::setStoreTrans(bool storeTrans) {
    store_trans = storeTrans;
}


///@internal generate the options string needed by the wrapper rttov_set_options subroutine
///does not perform any control upon the options
///return a string option example :
///... opts\%interpolation\%addinterp 1 opts\%rt_ir\%co2_data 1 opts\%rt_ir\%addclouds 1 nprofs_per_call 2 nthreads 1 verbose_wrapper 0 ..."
std::string Options::defineStrOptions() {

    std::string strOptions="";
    std::string key;
    bool val;

    for(StrBoolMap::iterator iter = this->config.begin(); iter != this->config.end(); ++iter)
    {
        std::string strval;
        key =  iter->first;
        val=this->config[key];
        if (val) strval.append(" 1"); else strval.append(" 0");
        strOptions.append(" opts%config%");
        strOptions.append(key);
        strOptions.append(strval);
    }
    for(StrBoolMap::iterator iter = this->interp.begin(); iter != this->interp.end(); ++iter)
    {
        std::string strval;
        key =  iter->first;
        val=this->interp[key];
        if (val) strval.append(" 1"); else strval.append(" 0");
        strOptions.append(" opts%interpolation%");
        strOptions.append(key);
        strOptions.append(strval);
    }
    for(StrBoolMap::iterator iter = this->rt_all.begin(); iter != this->rt_all.end(); ++iter)
    {
        std::string strval;
        key =  iter->first;
        val=this->rt_all[key];
        if (val) strval.append(" 1"); else strval.append(" 0");
        strOptions.append(" opts%rt_all%");
        strOptions.append(key);
        strOptions.append(strval);
    }
    for(StrBoolMap::iterator iter = this->mw.begin(); iter != this->mw.end(); ++iter)
    {
        std::string strval;
        key =  iter->first;
        val=this->mw[key];
        if (val) strval.append(" 1"); else strval.append(" 0");
        strOptions.append(" opts%rt_mw%");
        strOptions.append(key);
        strOptions.append(strval);
    }
    for(StrBoolMap::iterator iter = this->ir.begin(); iter != this->ir.end(); ++iter)
    {
        std::string strval;
        key =  iter->first;
        val=this->ir[key];
        if (val) strval.append(" 1"); else strval.append(" 0");
        strOptions.append(" opts%rt_ir%");
        strOptions.append(key);
        strOptions.append(strval);
    }
    for(StrBoolMap::iterator iter = this->pc.begin(); iter != this->pc.end(); ++iter)
    {
        std::string strval;
        key =  iter->first;
        val=this->pc[key];
        if (val) strval.append(" 1"); else strval.append(" 0");
        strOptions.append(" opts%rt_ir%pc%");
        strOptions.append(key);
        strOptions.append(strval);
    }

    std::ostringstream oss;
    if (this->isAddClouds()) {
        oss << " opts%rt_ir%cldstr_threshold "<< this->cldstr_threshold;
    }
    if (this->isAddPC()) {
        oss << " opts%rt_ir%pc%ipcbnd "<< this->ipcbnd;
        oss << " opts%rt_ir%pc%ipcreg "<< this->ipcreg;
    }
    oss << " opts%interpolation%interp_mode "<< this->interp_mode;
    oss << " opts%rt_mw%fastem_version "<< this->fastem_version;
    oss << " nthreads "<< this->nthreads;
    oss << " nprofs_per_call "<< this->nprofs_per_call;
    strOptions.append(oss.str());

    strOptions.append(" verbose_wrapper");
    val=this->verbose_wrapper;
    if (val) strOptions.append(" 1"); else strOptions.append(" 0");

    strOptions.append(" check_opts");
    val=this->check_opts;
    if (val) strOptions.append(" 1"); else strOptions.append(" 0");

    strOptions.append(" store_trans");
    val=this->store_trans;
    if (val) strOptions.append(" 1"); else strOptions.append(" 0");

    strOptions.append(" store_rad");
    val=this->store_rad;
    if (val) strOptions.append(" 1"); else strOptions.append(" 0");

    strOptions.append(" store_rad2");
    val=this->store_rad2;
    if (val) strOptions.append(" 1"); else strOptions.append(" 0");

    return strOptions;
}
} /* namespace rttov */
