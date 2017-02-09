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

///
/// @file Rttov.cpp
///
/// @brief Wrapper class for rttov
///

#include <Rttov.h>

namespace rttov {

bool is_readable( const std::string & file )
{
    std::ifstream fichier( file.c_str() );
    return fichier != 0;
}

/// @brief Rttov class constructor method
Rttov::Rttov() : coeffsLoaded(false),profileSet(false),
        p(NULL),t(NULL),gases(NULL),angles(NULL),
        datetimes(NULL),gas_id(NULL),gas_units(NULL),icecloud(NULL),inst_id(-1),
        nchannels(0),s2m(NULL),simplecloud(NULL),skin(NULL),surfgeom(NULL),surftype(NULL),
        zeeman(NULL),ngases(0),nprofiles(0),nlevels(0),rads(NULL),
        btrefl(NULL),surfemisrefl(NULL),wavenumbers(NULL),refPressures(NULL),
        allocatedSurfemisrefl(NULL),allocatedSurfemisreflk(NULL),rad_k(NULL),
        p_k(NULL),t_k(NULL),gases_k(NULL),bt_k(NULL),skin_k(NULL),s2m_k(NULL),
        surfemisrefl_k(NULL),simplecloud_k(NULL),profiles(NULL),
        allocatedBtk(NULL),allocatedRadk(NULL),nchannelsForLastRun(0),debug(false),
        isIrAtlasLoaded(false),isMwAtlasLoaded(false),isBrdfAtlasLoaded(false),allocatedP(false),
        tautotal(NULL),taulevels(NULL),tausuntotalpath1(NULL),tausunlevelspath1(NULL),
        tausuntotalpath2(NULL),tausunlevelspath2(NULL),radclear(NULL),radtotal(NULL),
        btclear(NULL),bt(NULL),reflclear(NULL),refl(NULL),radcloudy(NULL),overcast(NULL),
        rad2upclear(NULL),rad2dnclear(NULL),rad2refldnclear(NULL),rad2up(NULL),rad2down(NULL),
        rad2surf(NULL)
{

}

/// @brief Deallocate RTTOV structures btrefl and rads and drop the coefficients and atlases
Rttov::~Rttov() {
    // Deallocate RTTOV structures
    int err;

    if ( inst_id > 0 ) {
        rttov_drop_inst_(&err, &inst_id);
        if ( err!=0) {
            throw runtime_error("Error in rttov_drop_inst");
        }
    }

    delete [] btrefl;
    delete [] rads;
    delete [] bt_k;
    delete [] rad_k;
    delete [] surfemisrefl_k;
    if ( wavenumbers != NULL ) delete [] wavenumbers;
    if ( refPressures != NULL ) delete [] refPressures;
    if (this->allocatedSurfemisrefl) delete [] this->surfemisrefl;
    if (this->allocatedP) delete [] this->p;
    delete [] skin_k;
    delete [] s2m_k;
    delete [] simplecloud_k;
    delete [] p_k;
    delete [] t_k;
    delete [] gases_k;

    delete [] tautotal;
    delete [] taulevels;
    delete [] tausuntotalpath1;
    delete [] tausunlevelspath1;
    delete [] tausuntotalpath2;
    delete [] tausunlevelspath2;
    delete [] radclear;
    delete [] radtotal;
    delete [] btclear;
    delete [] bt;
    delete [] reflclear;
    delete [] refl;
    delete [] radcloudy;
    delete [] overcast;
    delete [] rad2upclear;
    delete [] rad2dnclear;
    delete [] rad2refldnclear;
    delete [] rad2up;
    delete [] rad2down;
    delete [] rad2surf;

    if (this->isBrdfAtlasLoaded) this->deallocBrdfAtlas();
    if (this->isIrAtlasLoaded)   this->deallocIrEmisAtlas();
    if (this->isMwAtlasLoaded)   this->deallocMwEmisAtlas();
}


///@brief Return the coefficient filename
const string& Rttov::getFileCoef() const {
    return file_coef;
}

///@brief Set the coefficient filename
///@param fileCoef string containing full path to "rtcoef" coefficient file
void Rttov::setFileCoef(const string& fileCoef) {
    file_coef = fileCoef;
}

///@brief Return the cloud coefficient filename
const string& Rttov::getFileSccld() const {
    return file_sccld;
}

///@brief Set the cloud coefficient filename
///@param fileSccld string containing full path to cloud coefficient file
void Rttov::setFileSccld(const string& fileSccld) {
    file_sccld = fileSccld;
}

///@brief Return the aerosol coefficient filename
const string& Rttov::getFileScaer() const {
    return file_scaer;
}

///@brief Set the aerosol coefficient filename
///@param fileScaer string containing full path to aerosol coefficient file
void Rttov::setFileScaer(const string& fileScaer) {
    file_scaer = fileScaer;
}

///@brief Return the path for the BRDF atlas files
const string& Rttov::getBrdfAtlasPath() const {
    return brdf_atlas_path;
}

///@brief Set the path for the BRDF atlas files
void Rttov::setBrdfAtlasPath(const string& brdfAtlasPath) {
    brdf_atlas_path = brdfAtlasPath;
}

///@brief Return the path for the emissivity atlas files
const string& Rttov::getEmisAtlasPath() const {
    return emis_atlas_path;
}

///@brief Set the path for the emissivity atlas files
void Rttov::setEmisAtlasPath(const string& emisAtlasPath) {
    emis_atlas_path = emisAtlasPath;
}


///@brief Load instrument with all channels
void Rttov::loadInst() {
    std::vector <int> cl(1);
    cl[0]=0;
    this->loadInst(cl);
}

///@brief Load instrument for a list of channels;
///       the method setFileCoef() must have been called previously
///@param channels vector containing the list of channels to be loaded
void Rttov::loadInst(const vector<int>& channels) {

    bool flagClouds=this->options.isAddClouds();
    bool flagAer=this->options.isAddAerosl();
    bool flagPC=this->options.isAddPC();
    string filesPaths;
    if (file_coef.empty()) {
        throw logic_error("Error: " + file_coef + " not set");
    }
    if (is_readable(file_coef)) {
        filesPaths.append(" file_coef ");
        filesPaths.append(file_coef + " ");
    }
    else
    {
        throw logic_error("Error: " + file_coef + " does not exist or is not readable");
    }

    if (! file_sccld.empty() ) {
        if ( is_readable(file_sccld))  {
            filesPaths.append(" file_sccld ");
            filesPaths.append(file_sccld + " ");
            this->options.setAddClouds(true);
        }
        else {
            throw logic_error("Error: " + file_sccld + " does not exist or is not readable");
        }
    }
    if (! file_scaer.empty() ) {
        if ( is_readable(file_scaer))  {
            filesPaths.append(" file_scaer ");
            filesPaths.append(file_scaer + " ");
            this->options.setAddAerosl(true);
        }
        else {
            throw logic_error("Error: " + file_scaer + " does not exist or is not readable");
        }
    }
    if (! file_pccoef.empty() ) {
        if ( is_readable(file_pccoef))  {
            filesPaths.append(" file_pccoef ");
            filesPaths.append(file_pccoef + " ");
        }
        else {
            throw logic_error("Error: " + file_pccoef + " does not exist or is not readable");
        }
    }

    this->coeffsLoaded = false;

    string optionsAndPaths;
    strOptions=this->options.defineStrOptions();
    optionsAndPaths.append(strOptions);
    optionsAndPaths.append(" ");
    optionsAndPaths.append(filesPaths);

    int len_opts_str=optionsAndPaths.length();
    char opts_str[len_opts_str];
    for (int c = 0; c < len_opts_str; c++ ) {
        opts_str[c]=optionsAndPaths[c];
    }
    if (this->debug) {
        std::cerr << optionsAndPaths << std::endl;
    }

    int nchan=channels.size();
    int c_channel_list[nchan];
    for (int i=0; i < nchan; i++){
        c_channel_list[i]=channels[i];
    }

    rttov_load_inst_(&inst_id, opts_str, &nchan, c_channel_list, len_opts_str);
    if (inst_id > 0) {
        int err=0;
        rttov_get_coef_val_i0_(&err,&inst_id,"NCHANNELS",&nchannels,9);
        if ( err != 0 ) throw runtime_error("Error getting number of channels");
        this->coeffsLoaded=true;
        this->options.setAddClouds(flagClouds);
        this->options.setAddAerosl(flagAer);
    }
    else {
        throw runtime_error("Error in rttov_load_inst");
    }
    if (this->options.isVerboseWrapper()) {
        std::cerr << "Load successful >>>>> inst_id : " << inst_id << ", nchannels : " << this->nchannels << std::endl;
    }
}


///@brief Return true if instrument is loaded
bool Rttov::isCoeffsLoaded() const {
    return coeffsLoaded;
}

///@brief Return the number of loaded channels
int Rttov::getNchannels() const {
    return nchannels;
}

///@brief Return the pressure levels of the coefficient file
double* Rttov::getRefPressures() {
    if ( this->isCoeffsLoaded() ){
        int err=0;
        int coeffNlevels = this->getCoeffsNlevels();
        if (refPressures != NULL) delete [] refPressures;
        refPressures=new double[coeffNlevels]();
        rttov_get_coef_val_r1_(&err,&inst_id,"REF_PRESSURE",&coeffNlevels,refPressures,12);
        if ( err != 0 ) throw runtime_error("Error getting pressures from coefficient");
        return refPressures;
    }
    else {
        return NULL;
    }
}

///@brief Return the number of levels of the coefficient file
int Rttov::getCoeffsNlevels() {
    if ( this->isCoeffsLoaded() ){
        int coeffNlevels;
        int err=0;
        rttov_get_coef_val_i0_(&err,&inst_id,"SIZE_REF_PRFL_P",&coeffNlevels,15);
        if ( err != 0 ) throw runtime_error("Error getting number of levels");
        return coeffNlevels;
    }
    else {
        return 0;
    }
}

///@brief Return the channel central wavenumbers of the coefficient file
double* Rttov::getWaveNumbers() {
    if ( this->isCoeffsLoaded() ){
        int err=0;
        if ( wavenumbers != NULL ) delete [] wavenumbers;
        wavenumbers=new double[nchannels]();
        if (this->debug) std::cerr << "Get wave numbers" << nchannels << std::endl;
        rttov_get_coef_val_r1_(&err,&inst_id,"WAVENUMBERS",&nchannels,wavenumbers,11);
        if ( err != 0 ) throw runtime_error("Error getting channel central wave numbers");
        return wavenumbers;
    }
    else {
        return NULL;
    }
}


///@brief Update RTTOV options for the currently loaded instrument
void Rttov::updateOptions() {
    if (! this-> coeffsLoaded  ) {
        throw logic_error( "Error: coefficients not loaded");
    }
    /// set some option in order to prevent crash
    if ( this->file_sccld.size() == 0  ) {
        this->options.setAddClouds(false);
    }
    if ( this->file_scaer .size() == 0  ) {
        this->options.setAddAerosl(false);
    }

    int err=0;
    strOptions=this->options.defineStrOptions();
    int len_opts_str = strOptions.length();
    char opts_str[len_opts_str];
    for (int c = 0; c < len_opts_str; c++ ) {
        opts_str[c]=strOptions[c];
    }
    if (this->debug) std::cerr << "Set RTTOV options " << strOptions << std::endl;
    rttov_set_options_(&err, &inst_id, opts_str, len_opts_str);
    if ( err != 0 ) {
        throw runtime_error("Error in rttov_set_options");
    }
}

///@brief Print RTTOV options for the currently loaded instrument
void Rttov::printOptions() {
    if (! this-> coeffsLoaded  ) {
        throw logic_error( "Error: coefficients not loaded");
    }
    int err=0;
    rttov_print_options_(&err, &inst_id);
    if ( err != 0 ) {
        throw runtime_error("Error in rttov_print_options");
    }
}


///@brief Initialise the BRDF atlas (all options available)
///@param[in] month month (1-12 => Jan-Dec) for which atlas should be initialised
///@param[in] single_inst bool initialise atlas for a single instrument, optional
///@param[in] version version number for IR atlas, optional
bool Rttov::brdfAtlasSetup(int month, bool single_inst, int version) {
    if (this->isBrdfAtlasLoaded) {
        if (this->options.isVerboseWrapper()) std::cerr << "Warning: BRDF atlas already initialised" << std::endl;
        return false;
    }
    if (! this->isCoeffsLoaded()) {
        if (this->options.isVerboseWrapper()) std::cerr << "Warning: instrument must be loaded first, BRDF atlas not initialised" << std::endl;
        return false;
    }
    if (this->debug) std::cerr << "brdf_atlas_path " << this->brdf_atlas_path << std::endl;
    int len_path = this->brdf_atlas_path.length();
    if (len_path == 0) {
        if (this->options.isVerboseWrapper()) std::cerr << "Warning: brdf_atlas_path not set, BRDF atlas not initialised" << std::endl;
        return false;
    }
    char path_str[len_path];
    for (int c = 0; c < len_path; c++ ) {
        path_str[c]=this->brdf_atlas_path[c];
    }
    int err;
    int inst_id_tmp=0;
    if (single_inst) inst_id_tmp=inst_id;
    rttov_brdf_atlas_setup_ (
            &err,
            path_str,
            &month,
            &version,
            &inst_id_tmp,
            len_path);
    if (err == 0) {
        this->isBrdfAtlasLoaded=true;
        if (this->options.isVerboseWrapper()) std::cerr << "BRDF atlas loaded successfully" << std::endl;
        return true;
    }
    else {
        if (this->options.isVerboseWrapper()) std::cerr << "Warning: BRDF atlas not (re-)initialised" << std::endl;
        return false;
    }
}

///@brief Initialise the IR emissivity atlas (all options available)
///@param[in] month month (1-12 => Jan-Dec) for which atlas should be initialised
///@param[in] ang_corr bool apply zenith angle correction to IR emissivities, optional
///@param[in] single_inst bool initialise atlas for a single instrument, optional
///@param[in] version version number for IR atlas, optional
bool Rttov::irEmisAtlasSetup(int month, bool ang_corr, bool single_inst, int version) {
    if (this->isIrAtlasLoaded) {
        if (this->options.isVerboseWrapper()) std::cerr << "Warning: IR emissivity atlas already initialised" << std::endl;
        return false;
    }
    if (! this->isCoeffsLoaded()) {
        if (this->options.isVerboseWrapper()) std::cerr << "Warning: instrument must be loaded first, IR emissivity atlas not initialised" << std::endl;
        return false;
    }
    if (this->debug) std::cerr << "emis_atlas_path " << this->emis_atlas_path << std::endl;
    int len_path = this->emis_atlas_path.length();
    if (len_path == 0) {
        if (this->options.isVerboseWrapper()) std::cerr << "Warning: emis_atlas_path not set, IR emissivity atlas not initialised" << std::endl;
        return false;
    }
    char path_str[len_path];
    for (int c = 0; c < len_path; c++ ) {
        path_str[c]=this->emis_atlas_path[c];
    }
    int err;
    int inst_id_tmp=0;
    if (single_inst) inst_id_tmp=inst_id;
    int int_ang_corr=0;
    if (ang_corr) int_ang_corr=1;
    rttov_ir_emis_atlas_setup_ (
            &err,
            path_str,
            &month,
            &version,
            &inst_id_tmp,
            &int_ang_corr,
            len_path);
    if (err == 0) {
        this->isIrAtlasLoaded=true;
        if (this->options.isVerboseWrapper()) std::cerr << "IR emissivity atlas loaded successfully" << std::endl;
        return true;
    }
    else {
        if (this->options.isVerboseWrapper()) std::cerr<< "Warning: IR emissivity atlas not (re-)initialised" << std::endl;
        return false;
    }
}

///@brief Initialise the MW emissivity atlas
///@param month month (1-12) of the atlas
///@param version optional (default value : 100)
bool Rttov::mwEmisAtlasSetup(int month, int version) {
    if (this->isMwAtlasLoaded) {
        if (this->options.isVerboseWrapper()) std::cerr << "Warning: MW emissivity atlas already initialised" << std::endl;
        return false;
    }
    if (! this->isCoeffsLoaded()) {
        if (this->options.isVerboseWrapper()) std::cerr<< "Warning: instrument must be loaded first, MW emissivity atlas not initialised" << std::endl;
        return false;
    }
    if (this->debug) std::cerr << "emis_atlas_path " << this->emis_atlas_path << std::endl;
    int len_path = this->emis_atlas_path.length();
    if (len_path == 0) {
        if (this->options.isVerboseWrapper()) std::cerr<< "Warning: emis_atlas_path not set, MW emissivity atlas not initialised" << std::endl;
        return false;
    }
    char path_str[len_path];
    for (int c = 0; c < len_path; c++ ) {
        path_str[c]=this->emis_atlas_path[c];
    }
    int err;
    rttov_mw_emis_atlas_setup_ (
            &err,
            path_str,
            &month,
            &version,
            &inst_id,
            len_path);
    if (err == 0) {
        this->isMwAtlasLoaded=true;
        if (this->options.isVerboseWrapper()) std::cerr << "MW emissivity atlas loaded successfully" << std::endl;
        return true;
    }
    else {
        if (this->options.isVerboseWrapper()) std::cerr << "Warning: MW emissivity atlas not (re-)initialised" << std::endl;
        return false;
    }
}

///@brief Initialise the BRDF atlas for the month of the first profile
bool Rttov::brdfAtlasSetup() {
    if (!  this->profileSet) {
        if (this->options.isVerboseWrapper()) std::cerr << "Warning: profiles must be set with valid month, BRDF atlas not initialised" << std::endl;
        return false;
    }
    int month=this->datetimes[1];
    brdfAtlasSetup(month);
}

///@brief Initialise the IR emissivity atlas for the month of the first profile
bool Rttov::irEmisAtlasSetup() {
    if (!  this->profileSet) {
        if (this->options.isVerboseWrapper()) std::cerr << "Warning: profiles must be set with valid month, IR emissivity atlas not initialised" << std::endl;
        return false;
    }
    int month=this->datetimes[1];
    irEmisAtlasSetup(month);
}

///@brief Initialise the MW emissivity atlas for the month of the first profile
bool Rttov::mwEmisAtlasSetup() {
    if (!  this->profileSet) {
        if (this->options.isVerboseWrapper()) std::cerr << "Warning: profiles must be set with valid month, MW emissivity atlas not initialised" << std::endl;
        return false;
    }
    int month=this->datetimes[1];
    mwEmisAtlasSetup(month);
}

///@brief Deallocate memory for the BRDF atlas
void Rttov::deallocBrdfAtlas() {
    if (this->isBrdfAtlasLoaded) rttov_brdf_atlas_dealloc_();
}

///@brief Deallocate memory for the IR emissivity atlas
void Rttov::deallocIrEmisAtlas() {
    if (this->isIrAtlasLoaded) rttov_ir_emis_atlas_dealloc_();
}

///@brief Deallocate memory for the MW emissivity atlas
void Rttov::deallocMwEmisAtlas() {
    if (this->isMwAtlasLoaded) rttov_mw_emis_atlas_dealloc_();
}


///@brief Set pointer to array containing input/output surface emissivity and reflectance values;
///this must be previously allocated a double array of dimensions [2][nprofiles][nchannels]; this
///is used to pass emissivity/reflectance values into RTTOV; if this is not called the Rttov object
///will allocate an array containing the values used by RTTOV which can be accessed by getSurfEmisRefl
///@param surfemisrefl pointer to a previously allocated double array [2][nprofiles][nchannels]
void Rttov::setSurfEmisRefl(double* surfemisrefl) {
    this->surfemisrefl = surfemisrefl;
}


/// @brief Associate a Profiles object with this Rttov object; 
///        this is fast, but does not carry out any checks on profiles before calling RTTOV
void  Rttov::setProfiles(rttov::Profiles * profiles) {
    this->profiles = profiles;
    if (allocatedP) delete[] p;
    if (this->profiles->isDefaultPressureLevels())
        if ( this->coeffsLoaded ) {
            int nlevcoeff=this->getCoeffsNlevels();
            if (nlevcoeff == this->profiles->getNlevels()) {
                if (this->options.isVerboseWrapper()) std::cerr << "Using pressure levels from coefficient file" << std::endl;
                double * pressure=this->getRefPressures();
                delete [] p;
                int nprof=this->profiles->getNprofiles();
                p= new double[nlevcoeff*this->profiles->getNprofiles()];
                for (int prof=0;prof<nprof;prof++)
                    for (int l=0;l<nlevcoeff;l++) p[prof*nlevcoeff+l]=pressure[l];
                allocatedP=true;
            }
            else throw logic_error("Error: number of levels differs to coefficient file, cannot use coefficient file pressure levels");
        }
        else throw logic_error("Error: instrument not loaded, cannot use coefficient file pressure levels");


    bool profileOK=this->profiles->check();
    if (! profileOK )  {
        this->profileSet=false;
        throw logic_error( "Error: some mandatory profile fields are missing");
    }

    // Make some pointers to simplify code
    this->nlevels=this->profiles->getNlevels();
    this->ngases=this->profiles->getNgases();
    this->nprofiles=this->profiles->getNprofiles();

    this->zeeman=this->profiles->getZeeman();
    this->skin=this->profiles->getSkin();
    this->simplecloud=this->profiles->getSimpleCloud();
    this->angles=this->profiles->getAngles();
    this->gas_id=this->profiles->getGasId();
    if (this->debug) for (int i=0; i < this->ngases; i++) std::cerr << "gasid" << " " << i << " : " << this->gas_id[i] << std::endl;
    this->datetimes=this->profiles->getDateTimes();
    this->surfgeom=this->profiles->getSurfGeom();
    this->s2m=this->profiles->getS2m();
    this->icecloud=this->profiles->getIceCloud();
    this->surftype=this->profiles->getSurfType();
    this->surfgeom=this->profiles->getSurfGeom();
    if (! allocatedP ) this->p=this->profiles->getP();
    this->gases=this->profiles->getGases();
    this->t=this->profiles->getT();
    this->gas_units=this->profiles->getGasUnits();

    this->profileSet=true;
}

///@brief Print gases array contents on standard output
void Rttov::printGases() {
    if (this->gas_id != NULL) {
        std::cout << "gas_id : " << std::endl;
        for (int g=0; g < ngases; g++) std::cout << "gases " << g << " gas_id " << gas_id[g] << std::endl;
    }
    if (this->p != NULL && this->t != NULL) {
        std::cout << "profile pressure t" << std::endl;
        for  (int prof=0; prof < nprofiles; prof++)
            for (int lev=0; lev < nlevels; lev++) {
                std::cout << prof << " " << p[prof*nlevels+lev] << " " << t[prof*nlevels+lev] << std::endl;
            }
    }
    if (this->gases != NULL) {
        std::cout << "profile gas level value" << std::endl;
        for  (int prof=0; prof < nprofiles; prof++)
            for (int g=0; g < ngases; g++)
                for (int lev=0; lev < nlevels; lev++) {
                    std::cout << prof << " " << g << " " << lev << " " << gases[g*nprofiles*nlevels+nlevels*prof+lev] << std::endl;
                }
    }
    if (this->surfemisrefl != NULL) {
        for  (int prof=0; prof < nprofiles; prof++){
            std::cout << "surfemisrefl prof " << prof << " " << this->surfemisrefl[0] << " " << this->surfemisrefl[1] << " " << std::endl;
        }
    }
}


///@brief Run the RTTOV direct model for all channels
void Rttov::runDirect() {
    vector<int> all_channels;
    for (int c=0; c < this->getNchannels(); c++) all_channels.push_back(c+1);
    this->runDirect(all_channels);
}

///@brief Run the RTTOV direct model for a list of channels
///@param channels vector containing the list of channels to simulate
void Rttov::runDirect(const vector<int>& channels) {
    if (! this-> coeffsLoaded  ) {
        throw logic_error( "Error: coefficients not loaded");
    }
    if (! this->profileSet) {
        throw logic_error( "Error: profiles not set");
    }

    int err = 0;

    // update the options which are locally stored in the class instance for RTTOV
    this->updateOptions();

    // when doing the run "chanprof" channel list must be now 1...nchan2
    int nchannels = channels.size();
    int channel_list[nchannels];
    for (int i = 0; i < nchannels; i++) {
        channel_list[i] = channels[i];
    }

    if (this->surfemisrefl == NULL) {
        if (this->options.isVerboseWrapper())
            std::cerr << "No surface emissivity/reflectance supplied: setting calcemis and calcrefl to true" << std::endl;
        this->allocatedSurfemisrefl=true;
        this->surfemisrefl = new double [2*nprofiles*nchannels]();
        for (int i=0;i<2*nprofiles*nchannels;i++ ) this->surfemisrefl[i]=-1.;
    }
    else if (allocatedSurfemisrefl) {
        // Reallocate in case nprofiles or nchannels changed since the previous call
        delete [] this->surfemisrefl;
        surfemisrefl = new double [2*nprofiles*nchannels]();
        for (int i=0; i < 2*nprofiles*nchannels; i++) surfemisrefl[i]=-1;
    }

    // alloc memory for btrefl and rads
    delete [] btrefl;
    delete [] rads;
    if (this->debug)
        std::cerr << "Allocating btrefl and rads for " << nprofiles << " profiles and " << nchannels << " channels" << std::endl;
    btrefl = new double[nprofiles*nchannels]();
    rads = new double[nprofiles*nchannels]();
    if (this->debug) {
        std::cerr << "gas_units= " << this->gas_units << std::endl;
        for (int p=0; p< nprofiles;p++) {
            std::cerr << "Profile " << p << std::endl;
            std::cerr << "datetimes ";
            for (int i=0;i<6;i++) std::cerr << datetimes[6*p+i]<<" ";
            std::cerr <<std::endl;
            std::cerr << "angles ";
            for (int i=0;i<4;i++) std::cerr << angles[4*p+i]<<" ";
            std::cerr <<std::endl;
            std::cerr << "surfgeom ";
            for (int i=0;i<3;i++) std::cerr << surfgeom[3*p+i]<<" ";
            std::cerr <<std::endl;
            std::cerr << "surftype ";
            for (int i=0;i<2;i++) std::cerr << surftype[2*p+i]<<" ";
            std::cerr <<std::endl;
            std::cerr << "skin ";
            for (int i=0;i<9;i++) std::cerr << skin[9*p+i]<<" ";
            std::cerr <<std::endl;
            std::cerr << "s2m ";
            for (int i=0;i<6;i++) std::cerr << s2m[6*p+i]<<" ";
            std::cerr <<std::endl;
            std::cerr << "simplecloud ";
            for (int i=0;i<2;i++) std::cerr << simplecloud[2*p+i]<<" ";
            std::cerr <<std::endl;
            std::cerr << "icecloud ";
            for (int i=0;i<2;i++) std::cerr << icecloud[2*p+i]<<" ";
            std::cerr <<std::endl;
            std::cerr << "zeeman ";
            for (int i=0;i<2;i++) std::cerr << zeeman[2*p+i]<<" ";
            std::cerr <<std::endl;
            std::cerr <<std::endl;
        }
    }

    rttov_call_direct_(
            &err,
            &inst_id,
            channel_list,
            datetimes,         // profile dates/times                                     [nprofiles][6]
            angles,            // satzen, satazi, sunzen, sunazi angles                   [nprofiles][4]
            surfgeom,          // lat, lon, elevation                                     [nprofiles][3]
            surftype,          // surftype, watertype                                     [nprofiles][2]
            skin,              // skin T, salinity, snow_frac, foam_frac, fastem_coefsx5  [nprofiles][9]
            s2m,               // 2m p, 2m t, 2m q, 10m wind u, v, wind-fetch             [nprofiles][6]
            simplecloud,       // ctp, cfraction                                          [nprofiles][2]
            icecloud,          // ish, idg                                                [nprofiles][2]
            zeeman,            // Be, cosbk                                               [nprofiles][2]
            p,                 // pressure                                                [nprofiles][nlevels]
            t,                 // temperature                                             [nprofiles][nlevels]
            &gas_units,        // units for gas profiles
            gas_id,            // gas ID list                                             [ngases]
            gases,             // gas profiles                                            [ngases][nprofiles][nlevels]
            surfemisrefl,      // input/output surface emissivities/BRDFs                 [2][nprofiles][nchannels]
            btrefl,            // output BTs/refls (for thermal/solar chans)              [nprofiles][nchannels]
            rads,              // output radiances                                        [nprofiles][nchannels]
            &nchannels, &ngases, &nlevels, &nprofiles);
    if ( err != 0 ) throw runtime_error ("Error in rttov_call_direct");
    this->nchannelsForLastRun=nchannels;

    if (options.isStoreTrans()) doStoreTrans();
    if (options.isStoreRad())   doStoreRad();
    if (options.isStoreRad2())  doStoreRad2();
}


///@brief Run the RTTOV K model for all channels
void Rttov::runK() {
    vector<int> all_channels;
    for (int c=0; c < this->getNchannels(); c++) all_channels.push_back(c+1);
    this->runK(all_channels);
}

///@brief Run the RTTOV K model for a list of channels
///@param channels vector containing the list of channels to simulate
void Rttov::runK(const vector<int>& channels) {
    int err=0;

    if (! this-> coeffsLoaded  ) {
        throw logic_error( "Error: coefficients not loaded");
    }
    if (! this->profileSet) {
        throw logic_error( "Error: profiles not set");
    }

    // update the options which are locally stored in the class instance for RTTOV
    this->updateOptions();

    //   // when doing the run "chanprof" channel list must be now 1...nchan2
    int nchannels = channels.size();
    int channel_list[nchannels];
    for (int i = 0; i < nchannels; i++) {
        channel_list[i] = channels[i];
    }
    // allocate arrays for rttov_call_k output
    delete [] skin_k;
    skin_k= new double [nprofiles*nchannels*9]();
    delete [] s2m_k;
    s2m_k= new double [nprofiles*nchannels*6]();
    delete [] simplecloud_k;
    simplecloud_k = new double [nprofiles*nchannels*2]();
    delete [] p_k;
    p_k= new double [nprofiles*nchannels*nlevels]();
    delete [] t_k;
    t_k = new double [nprofiles*nchannels*nlevels]();
    delete [] gases_k;
    gases_k= new double [ngases*nprofiles*nchannels*nlevels]();

    if (surfemisrefl == NULL) {
        if (this->options.isVerboseWrapper())
            std::cerr << "No surface emissivity/reflectance supplied: setting calcemis and calcrefl to true" << std::endl;
        surfemisrefl = new double [2*nprofiles*nchannels]();
        for (int i=0; i < 2*nprofiles*nchannels; i++) surfemisrefl[i]=-1;
        allocatedSurfemisrefl=true;
    }
    else if (allocatedSurfemisrefl) {
        // Reallocate in case nprofiles or nchannels changed since the previous call
        delete [] this->surfemisrefl;
        surfemisrefl = new double [2*nprofiles*nchannels]();
        for (int i=0; i < 2*nprofiles*nchannels; i++) surfemisrefl[i]=-1;
    }

    delete [] this->surfemisrefl_k;
    surfemisrefl_k = new double [2*nprofiles*nchannels]();
    for (int i=0; i < 2*nprofiles*nchannels; i++) surfemisrefl_k[i]=0;

    delete [] this->bt_k;
    bt_k = new double [nprofiles*nchannels]();
    for (int i=0; i < nprofiles*nchannels; i++) bt_k[i]=1;

    delete [] this->rad_k;
    rad_k = new double [nprofiles*nchannels]();
    for (int i=0; i < nprofiles*nchannels; i++) rad_k[i]=1;

    if (this->debug) {
        std::cerr << "inst_id " << inst_id << std::endl;
        std::cerr << "channels " << nchannels << std::endl;
        std::cerr << "ngases " << ngases << std::endl;
        std::cerr << "nlevels " <<nlevels << std::endl;
        std::cerr << "nprofiles " <<nprofiles << std::endl;
    }

    delete [] btrefl;
    delete [] rads;
    if (this->debug)
        std::cerr << "Allocating btrefl and rads for " << nprofiles << " profiles and " << nchannels << " channels" << std::endl;
    btrefl = new double[nprofiles*nchannels]();
    rads = new double[nprofiles*nchannels]();
    rttov_call_k_(
            &err,
            &inst_id,
            channel_list,
            datetimes,         // profile dates/times                                     [nprofiles][6]
            angles,            // satzen, satazi, sunzen, sunazi angles                   [nprofiles][4]
            surfgeom,          // lat, lon, elevation                                     [nprofiles][3]
            surftype,          // surftype, watertype                                     [nprofiles][2]
            skin,              // skin T, salinity, snow_frac, foam_frac, fastem_coefsx5  [nprofiles][9]
            skin_k,            // output skin K                                           [nprofiles][nchannels][9]
            s2m,               // 2m p, 2m t, 2m q, 10m wind u, v, wind-fetch             [nprofiles][6]
            s2m_k,             // output 2m K                                             [nprofiles][nchannels][6]
            simplecloud,       // ctp, cfraction                                          [nprofiles][2]
            simplecloud_k,     // output ctp, cfraction K                                 [nprofiles][nchannels][2]
            icecloud,          // ish, idg                                                [nprofiles][2]
            zeeman,            // Be, cosbk                                               [nprofiles][2]
            p,                 // pressure                                                [nprofiles][nlevels]
            p_k,               // output pressure K                                       [nprofiles][nchannels][nlevels]
            t,                 // temperature                                             [nprofiles][nlevels]
            t_k,               // output temperature K                                    [nprofiles][nchannels][nlevels]
            &gas_units,        // units for gas profiles
            gas_id,            // gas ID list                                             [ngases]
            gases,             // gas profiles                                            [ngases][nprofiles][nlevels]
            gases_k,           // output gas profiles K                                   [ngases][nprofiles][nchannels][nlevels]
            surfemisrefl,      // input/output surface emissivities/BRDFs                 [2][nprofiles][nchannels]
            surfemisrefl_k,    // input/output surface emissivities/BRDFs K               [2][nprofiles][nchannels]
            btrefl,            // output BTs/refls (for thermal/solar chans)              [nprofiles][nchannels]
            rads,              // output radiances                                        [nprofiles][nchannels]
            bt_k,              // input BT perturbations                                  [nprofiles][nchannels]
            rad_k,             // input radiance perturbations                            [nprofiles][nchannels]
            &nchannels, &ngases, &nlevels, &nprofiles);

    if ( err != 0 ) throw runtime_error ("Error in rttov_call_k");
    this->nchannelsForLastRun=nchannels;

    if (options.isStoreTrans()) doStoreTrans();
    if (options.isStoreRad())   doStoreRad();
}


///@brief Return a pointer to an array of dimensions [2][nprofiles][nchannels]
/// containing output values of surface emissivity and reflectance;
/// this array can be initialised by the user and set by calling the setSurfEmisRefl method;
/// alternatively if the emissivity/reflectance array is allocated by the Rttov object
/// it is deleted at the next run or when the Rttov instance is destroyed.
const double* Rttov::getSurfEmisRefl() const {
    return surfemisrefl;
}


///@brief Return a pointer to an array of dimensions [nprofiles][nchannels]
/// filled with computed brightness temperatures and reflectances by the previous run;
/// this array is allocated by the Rttov object and is destroyed when
/// a new run is performed or if the instance is destroyed
const double* Rttov::getBtRefl() const {
    return btrefl;
}

///@brief Return a pointer to an array of dimensions [nprofiles][nchannels]
/// filled with computed radiances by the previous run;
/// this array is allocated by the Rttov object and is destroyed when
/// a new run is performed or if the instance is destroyed
const double* Rttov::getRads() const {
    return rads;
}


///@brief Return a vector for x y z from a pointer to [dim1][dim2][dim3][dim4]
std::vector<double> Rttov::convertPointer4D2Vector(double ptr[], int x, int y,
        int z, int dim1, int dim2, int dim3, int dim4) {
    std::vector <double>myvector{};
    if (ptr == NULL) throw logic_error("Error: data not present");
    if (x < 0 || x >= dim1) throw logic_error("Error: wrong call to convertPointer4D2Vector");
    if (y < 0 || y >= dim2) throw logic_error("Error: wrong call to convertPointer4D2Vector");
    if (z < 0 || z >= dim3) throw logic_error("Error: wrong call to convertPointer4D2Vector");
    if (ptr == NULL) throw logic_error("Error: wrong call to convertPointer4D2Vector");
    for (int w=0; w < dim4; w++) myvector.push_back(ptr[x*dim4*dim3*dim2 + dim4*dim3*y + dim4*z+w]);
    return myvector;
}

///@brief Return a vector for x and y from a pointer to [dim1][dim2][dim3]
std::vector<double> Rttov::convertPointer3D2Vector(double ptr[], const int x,
        const int y, const int dim1, const int dim2, const int dim3) {
    std::vector <double>myvector{};
    if (ptr == NULL) throw logic_error("Error: data not present");
    if (x < 0 || x >= dim1) throw logic_error("Error: wrong call to convertPointer3D2Vector");
    if (y < 0 || y >= dim2) throw logic_error("Error: wrong call to convertPointer3D2Vector");
    if (ptr == NULL) throw logic_error("Error: wrong call to convertPointer3D2Vector");
    for (int z=0; z < dim3; z++) myvector.push_back(ptr[x*dim3*dim2 + dim3*y + z]);
    return myvector;
}

///@brief Return a vector for x from a pointer to [dim1][dim2]
std::vector<double> Rttov::convertPointer2D2Vector(double ptr[], int x, int dim1, int dim2) {
    if (ptr == NULL) throw logic_error("Error: data not present");
    if (x < 0 || x >= dim1) throw logic_error("Error: wrong call to convertPointer2D2Vector");
    if (ptr == NULL) throw logic_error("Error: wrong call to convertPointer2D2Vector");
    std::vector <double>myvector{};
    for (int y=0; y < dim2; y++) {
        myvector.push_back(ptr[x*dim2 + y]);
    }
    return myvector;
}

///@brief Return vector of brightness temperatures/reflectances
///       computed by the previous run for the given profile number
///@param profile index of the profile (0..nprofiles-1)
std::vector<double> Rttov::getBtRefl(const int profile)  {
    return this->convertPointer2D2Vector(this->btrefl,
            profile,
            this->nprofiles,
            this->nchannelsForLastRun);
}

///@brief Return a vector of radiances computed by the previous run
///       for the given profile number
///@param profile index of the profile (0..nprofiles-1)
std::vector<double> Rttov::getRads(int profile)  {
    return this->convertPointer2D2Vector(this->rads,
            profile,
            this->nprofiles,
            this->nchannelsForLastRun);
}

///@brief Return the computed pressure Jacobians for a given profile and channel
///@param profile index of the profile (0..nprofiles-1)
///@param channel index of the channel (0..nchannels-1)
std::vector<double> Rttov::getPK(int profile, int channel) {
    return this->convertPointer3D2Vector(this->p_k,
            profile,
            channel,
            this->nprofiles,
            this->nchannelsForLastRun,
            this->nlevels);
}

///@brief Return computed temperature Jacobians for a given profile and channel
///@param profile index of the profile (0..nprofiles-1)
///@param channel index of the channel (0..nchannels-1)
std::vector<double> Rttov::getTK(int profile, int channel) {
    return this->convertPointer3D2Vector(this->t_k,
            profile,
            channel,
            this->nprofiles,
            this->nchannelsForLastRun,
            this->nlevels);
}

///@brief Return computed skin variable Jacobians for a given profile and channel
///@param profile index of the profile (0..nprofiles-1)
///@param channel index of the channel (0..nchannels-1)
std::vector<double> Rttov::getSkinK(int profile, int channel) {
    return this->convertPointer3D2Vector(this->skin_k,
            profile,
            channel,
            this->nprofiles,
            this->nchannelsForLastRun,
            9);
}

///@brief Return computed 2m variable Jacobian for a given profile and channel
///@param profile index of the profile (0..nprofiles-1)
///@param channel index of the channel (0..nchannels-1)
std::vector<double> Rttov::getS2mK(int profile, int channel) {
    return this->convertPointer3D2Vector(this->s2m_k,
            profile,
            channel,
            this->nprofiles,
            this->nchannelsForLastRun,
            6);
}

///@brief Return computed simple cloud variable Jacobians for a given profile and channel
///@param profile index of the profile (0..nprofiles-1)
///@param channel index of the channel (0..nchannels-1)
std::vector<double> Rttov::getSimpleCloudK(int profile, int channel) {
    return this->convertPointer3D2Vector(this->simplecloud_k,
            profile,
            channel,
            this->nprofiles,
            this->nchannelsForLastRun,
            2);
}

///@brief Return computed gas, cloud and aerosol Jacobian values for a given profile and channel
///@param itemIdType item type
///@param profile index of the profile (0..nprofiles-1)
///@param channel index of the channel (0..nchannels-1)
std::vector<double> Rttov::getItemK(rttov::itemIdType itemIdType, int profile,
        int channel) {
    return this->convertPointer4D2Vector(this->gases_k,
            this->profiles->gas_index[itemIdType],
            profile,
            channel,
            this->ngases,
            this->nprofiles,
            this->nchannelsForLastRun,
            this->nlevels);
}

///@brief Return computed surface emissivity Jacobians for a given profile
///@param profile index of the profile (0..nprofiles-1)
std::vector<double> Rttov::getSurfEmisK(int profile) {
    return this->convertPointer3D2Vector(this->surfemisrefl_k,
            0,
            profile,
            2,
            this->nprofiles,
            this->nchannelsForLastRun);
}

///@brief Return computed surface reflectance Jacobians for a given profile
///@param profile index of the profile (0..nprofiles-1)
std::vector<double> Rttov::getSurfReflK(int profile) {
    return this->convertPointer3D2Vector(this->surfemisrefl_k,
            1,
            profile,
            2,
            this->nprofiles,
            this->nchannelsForLastRun);
}


void Rttov::doStoreTrans() {
    int err = 0;
    int nchanprof = nprofiles * nchannelsForLastRun;

    delete [] tautotal;
    tautotal = new double [nchanprof]();
    rttov_get_tau_total_(&err, &inst_id, tautotal, &nchanprof);
    if ( err != 0 ) throw runtime_error("Error getting tautotal");

    delete [] taulevels;
    taulevels = new double [nchanprof*nlevels]();
    rttov_get_tau_levels_(&err, &inst_id, taulevels, &nchanprof, &nlevels);
    if ( err != 0 ) throw runtime_error("Error getting taulevels");

    delete [] tausuntotalpath1;
    tausuntotalpath1 = new double [nchanprof]();
    rttov_get_tausun_total_path1_(&err, &inst_id, tausuntotalpath1, &nchanprof);
    if ( err != 0 ) throw runtime_error("Error getting tausuntotalpath1");

    delete [] tausunlevelspath1;
    tausunlevelspath1 = new double [nchanprof*nlevels]();
    rttov_get_tausun_levels_path1_(&err, &inst_id, tausunlevelspath1, &nchanprof, &nlevels);
    if ( err != 0 ) throw runtime_error("Error getting tausunlevelspath1");

    delete [] tausuntotalpath2;
    tausuntotalpath2 = new double [nchanprof]();
    rttov_get_tausun_total_path2_(&err, &inst_id, tausuntotalpath2, &nchanprof);
    if ( err != 0 ) throw runtime_error("Error getting tausuntotalpath2");

    delete [] tausunlevelspath2;
    tausunlevelspath2 = new double [nchanprof*nlevels]();
    rttov_get_tausun_levels_path2_(&err, &inst_id, tausunlevelspath2, &nchanprof, &nlevels);
    if ( err != 0 ) throw runtime_error("Error getting tausunlevelspath2");
}

void Rttov::doStoreRad() {
    int err = 0;
    int nchanprof = nprofiles * nchannelsForLastRun;
    int nlayers = nlevels-1;

    delete [] radclear;
    radclear = new double [nchanprof]();
    rttov_get_rad_clear_(&err, &inst_id, radclear, &nchanprof);
    if ( err != 0 ) throw runtime_error("Error getting radclear");

    delete [] radtotal;
    radtotal = new double [nchanprof]();
    rttov_get_rad_total_(&err, &inst_id, radtotal, &nchanprof);
    if ( err != 0 ) throw runtime_error("Error getting radtotal");

    delete [] btclear;
    btclear = new double [nchanprof]();
    rttov_get_bt_clear_(&err, &inst_id, btclear, &nchanprof);
    if ( err != 0 ) throw runtime_error("Error getting btclear");

    delete [] bt;
    bt = new double [nchanprof]();
    rttov_get_bt_(&err, &inst_id, bt, &nchanprof);
    if ( err != 0 ) throw runtime_error("Error getting bt");

    delete [] reflclear;
    reflclear = new double [nchanprof]();
    rttov_get_refl_clear_(&err, &inst_id, reflclear, &nchanprof);
    if ( err != 0 ) throw runtime_error("Error getting reflclear");

    delete [] refl;
    refl = new double [nchanprof]();
    rttov_get_refl_(&err, &inst_id, refl, &nchanprof);
    if ( err != 0 ) throw runtime_error("Error getting refl");

    delete [] radcloudy;
    radcloudy = new double [nchanprof]();
    rttov_get_rad_cloudy_(&err, &inst_id, radcloudy, &nchanprof);
    if ( err != 0 ) throw runtime_error("Error getting radcloudy");

    delete [] overcast;
    overcast = new double [nchanprof*nlayers]();
    rttov_get_overcast_(&err, &inst_id, overcast, &nchanprof, &nlayers);
    if ( err != 0 ) throw runtime_error("Error getting overcast");
}

void Rttov::doStoreRad2() {
    int err = 0;
    int nchanprof = nprofiles * nchannelsForLastRun;
    int nlayers = nlevels-1;

    delete [] rad2upclear;
    rad2upclear = new double [nchanprof]();
    rttov_get_rad2_upclear_(&err, &inst_id, rad2upclear, &nchanprof);
    if ( err != 0 ) throw runtime_error("Error getting rad2upclear");

    delete [] rad2dnclear;
    rad2dnclear = new double [nchanprof]();
    rttov_get_rad2_dnclear_(&err, &inst_id, rad2dnclear, &nchanprof);
    if ( err != 0 ) throw runtime_error("Error getting rad2dnclear");

    delete [] rad2refldnclear;
    rad2refldnclear = new double [nchanprof]();
    rttov_get_rad2_refldnclear_(&err, &inst_id, rad2refldnclear, &nchanprof);
    if ( err != 0 ) throw runtime_error("Error getting rad2refldnclear");

    delete [] rad2up;
    rad2up = new double [nchanprof*nlayers]();
    rttov_get_rad2_up_(&err, &inst_id, rad2up, &nchanprof, &nlayers);
    if ( err != 0 ) throw runtime_error("Error getting rad2up");

    delete [] rad2down;
    rad2down = new double [nchanprof*nlayers]();
    rttov_get_rad2_down_(&err, &inst_id, rad2down, &nchanprof, &nlayers);
    if ( err != 0 ) throw runtime_error("Error getting rad2down");

    delete [] rad2surf;
    rad2surf = new double [nchanprof*nlayers]();
    rttov_get_rad2_surf_(&err, &inst_id, rad2surf, &nchanprof, &nlayers);
    if ( err != 0 ) throw runtime_error("Error getting rad2surf");
}

///@brief Return RTTOV transmission tau_total output array of size [nchannels] for given profile, requires store_trans true
std::vector <double> Rttov::getTauTotal(int profile) {
    return this->convertPointer2D2Vector(this->tautotal,
            profile,
            this->nprofiles,
            this->nchannelsForLastRun);
}

///@brief Return RTTOV transmission tau_levels output array of size [nlevels] for given profile and channel, requires store_trans true
std::vector <double> Rttov::getTauLevels(int profile, int channel) {
    return this->convertPointer3D2Vector(this->taulevels,
            profile,
            channel,
            this->nprofiles,
            this->nchannelsForLastRun,
            this->nlevels);
}

///@brief Return RTTOV transmission tausun_total_path1 output array of size [nchannels] for given profile, requires store_trans true
std::vector <double> Rttov::getTauSunTotalPath1(int profile) {
    return this->convertPointer2D2Vector(this->tausuntotalpath1,
            profile,
            this->nprofiles,
            this->nchannelsForLastRun);
}

///@brief Return RTTOV transmission tausun_levels_path1 output array of size [nlevels] for given profile and channel, requires store_trans true
std::vector <double> Rttov::getTauSunLevelsPath1(int profile, int channel) {
    return this->convertPointer3D2Vector(this->tausunlevelspath1,
            profile,
            channel,
            this->nprofiles,
            this->nchannelsForLastRun,
            this->nlevels);
}

///@brief Return RTTOV transmission tausun_total_path2 output array of size [nchannels] for given profile, requires store_trans true
std::vector <double> Rttov::getTauSunTotalPath2(int profile) {
    return this->convertPointer2D2Vector(this->tausuntotalpath2,
            profile,
            this->nprofiles,
            this->nchannelsForLastRun);
}

///@brief Return RTTOV transmission tausun_levels_path2 output array of size [nlevels] for given profile and channel, requires store_trans true
std::vector <double> Rttov::getTauSunLevelsPath2(int profile, int channel) {
    return this->convertPointer3D2Vector(this->tausunlevelspath2,
            profile,
            channel,
            this->nprofiles,
            this->nchannelsForLastRun,
            this->nlevels);
}


///@brief Return RTTOV radiance clear output array of size [nchannels] for given profile, requires store_rad true
std::vector <double> Rttov::getRadClear(int profile) {
    return this->convertPointer2D2Vector(this->radclear,
            profile,
            this->nprofiles,
            this->nchannelsForLastRun);
}

///@brief Return RTTOV radiance total output array of size [nchannels] for given profile, requires store_rad true
std::vector <double> Rttov::getRadTotal(int profile) {
    return this->convertPointer2D2Vector(this->radtotal,
            profile,
            this->nprofiles,
            this->nchannelsForLastRun);
}

///@brief Return RTTOV radiance bt_clear output array of size [nchannels] for given profile, requires store_rad true
std::vector <double> Rttov::getBtClear(int profile) {
    return this->convertPointer2D2Vector(this->btclear,
            profile,
            this->nprofiles,
            this->nchannelsForLastRun);
}

///@brief Return RTTOV radiance bt output array of size [nchannels] for given profile, requires store_rad true
std::vector <double> Rttov::getBt(int profile) {
    return this->convertPointer2D2Vector(this->bt,
            profile,
            this->nprofiles,
            this->nchannelsForLastRun);
}

///@brief Return RTTOV radiance refl_clear output array of size [nchannels] for given profile, requires store_rad true
std::vector <double> Rttov::getReflClear(int profile) {
    return this->convertPointer2D2Vector(this->reflclear,
            profile,
            this->nprofiles,
            this->nchannelsForLastRun);
}

///@brief Return RTTOV radiance refl output array of size [nchannels] for given profile, requires store_rad true
std::vector <double> Rttov::getRefl(int profile) {
    return this->convertPointer2D2Vector(this->refl,
            profile,
            this->nprofiles,
            this->nchannelsForLastRun);
}

///@brief Return RTTOV radiance cloudy output array of size [nchannels] for given profile, requires store_rad true
std::vector <double> Rttov::getRadCloudy(int profile) {
    return this->convertPointer2D2Vector(this->radcloudy,
            profile,
            this->nprofiles,
            this->nchannelsForLastRun);
}

///@brief Return RTTOV radiance overcast output array of size [nlayers] for given profile and channel, requires store_rad true
std::vector <double> Rttov::getOvercast(int profile, int channel) {
    return this->convertPointer3D2Vector(this->overcast,
            profile,
            channel,
            this->nprofiles,
            this->nchannelsForLastRun,
            this->nlevels-1);
}


///@brief Return RTTOV radiance2 upclear output array of size [nchannels] for given profile, requires store_rad2 true
std::vector <double> Rttov::getRad2UpClear(int profile) {
    return this->convertPointer2D2Vector(this->rad2upclear,
            profile,
            this->nprofiles,
            this->nchannelsForLastRun);
}

///@brief Return RTTOV radiance2 dnclear output array of size [nchannels] for given profile, requires store_rad2 true
std::vector <double> Rttov::getRad2DnClear(int profile) {
    return this->convertPointer2D2Vector(this->rad2dnclear,
            profile,
            this->nprofiles,
            this->nchannelsForLastRun);
}

///@brief Return RTTOV radiance2 refldnclear output array of size [nchannels] for given profile, requires store_rad2 true
std::vector <double> Rttov::getRad2ReflDnClear(int profile) {
    return this->convertPointer2D2Vector(this->rad2refldnclear,
            profile,
            this->nprofiles,
            this->nchannelsForLastRun);
}

///@brief Return RTTOV radiance2 up output array of size [nlayers] for given profile and channel, requires store_rad2 true
std::vector <double> Rttov::getRad2Up(int profile, int channel) {
    return this->convertPointer3D2Vector(this->rad2up,
            profile,
            channel,
            this->nprofiles,
            this->nchannelsForLastRun,
            this->nlevels-1);
}

///@brief Return RTTOV radiance2 down output array of size [nlayers] for given profile and channel, requires store_rad2 true
std::vector <double> Rttov::getRad2Down(int profile, int channel) {
    return this->convertPointer3D2Vector(this->rad2down,
            profile,
            channel,
            this->nprofiles,
            this->nchannelsForLastRun,
            this->nlevels-1);
}

///@brief Return RTTOV radiance2 surf output array of size [nlayers] for given profile and channel, requires store_rad2 true
std::vector <double> Rttov::getRad2Surf(int profile, int channel) {
    return this->convertPointer3D2Vector(this->rad2surf,
            profile,
            channel,
            this->nprofiles,
            this->nchannelsForLastRun,
            this->nlevels-1);
}
} /* namespace rttov */
