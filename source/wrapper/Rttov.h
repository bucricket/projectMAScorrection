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

/// @class Rttov
///  This class contains a set of functions for running Rttov.
///  An instance of an Rttov Class is defined to handle one instrument.
///  It contains an Options object which allows the Rttov options to
///  be specified and can be associated with a Profiles object.


#ifndef RTTOV_H_
#define RTTOV_H_


using namespace std;
#include <rttov_cc_interface.h>

#include <fstream>
#include <iostream>
#include <cstdlib>
#include <exception>
#include <stdexcept>
#include <Rttov_common.h>
#include <Profiles.h>
#include <Profile.h>
#include <Options.h>

namespace rttov {

class Rttov {

public:
    Rttov();
    virtual ~Rttov();

    Options::Options options;

    const string& getFileCoef() const;
    const string& getFileSccld() const;
    const string& getFileScaer() const;
    const string& getBrdfAtlasPath() const;
    const string& getEmisAtlasPath() const;

    void setFileCoef(const string& fileCoef);
    void setFileSccld(const string& fileSccld);
    void setFileScaer(const string& fileScaer);
    void setBrdfAtlasPath(const string& brdfAtlasPath);
    void setEmisAtlasPath(const string& emisAtlasPath);

    void loadInst(); // load instrument for all channels
    void loadInst(const vector<int>& channels);

    bool isCoeffsLoaded() const;
    int getNchannels() const;
    double * getRefPressures();
    int getCoeffsNlevels();
    double * getWaveNumbers();

    void updateOptions();
    void printOptions();

    bool brdfAtlasSetup(int month, bool single_inst=false, int version=-1);
    bool irEmisAtlasSetup(int month, bool ang_corr=false, bool single_inst=false, int version=-1);
    bool mwEmisAtlasSetup(int month, int version=-1);
    bool brdfAtlasSetup();
    bool irEmisAtlasSetup();
    bool mwEmisAtlasSetup();
    void deallocBrdfAtlas();
    void deallocIrEmisAtlas();
    void deallocMwEmisAtlas();

    void setSurfEmisRefl(double* surfemisrefl);

    void setProfiles(rttov::Profiles * profiles);
    void printGases();

    void runDirect(); // run direct for all loaded channels
    void runDirect(const vector<int>& channels);
    void runK() ;
    void runK(const vector<int>& channels);

    const double* getBtRefl() const;
    const double* getRads() const;
    std::vector<double> getBtRefl(const int profile);
    std::vector<double> getRads(const int profile);
    const double* getSurfEmisRefl() const;

    std::vector <double> getPK(int profile, int channel);
    std::vector <double> getTK(int profile, int channel);
    std::vector <double> getSkinK(int profile, int channel);
    std::vector <double> getS2mK(int profile, int channel);
    std::vector <double> getSimpleCloudK(int profile, int channel);
    std::vector <double> getItemK(rttov::itemIdType, int profile, int channel);
    std::vector <double> getSurfEmisK(int profile);
    std::vector <double> getSurfReflK(int profile);

    std::vector <double> getTauTotal(int profile);
    std::vector <double> getTauLevels(int profile, int channel);
    std::vector <double> getTauSunTotalPath1(int profile);
    std::vector <double> getTauSunLevelsPath1(int profile, int channel);
    std::vector <double> getTauSunTotalPath2(int profile);
    std::vector <double> getTauSunLevelsPath2(int profile, int channel);

    std::vector <double> getRadClear(int profile);
    std::vector <double> getRadTotal(int profile);
    std::vector <double> getBtClear(int profile);
    std::vector <double> getBt(int profile);
    std::vector <double> getReflClear(int profile);
    std::vector <double> getRefl(int profile);
    std::vector <double> getRadCloudy(int profile);
    std::vector <double> getOvercast(int profile, int channel);

    std::vector <double> getRad2UpClear(int profile);
    std::vector <double> getRad2DnClear(int profile);
    std::vector <double> getRad2ReflDnClear(int profile);
    std::vector <double> getRad2Up(int profile, int channel);
    std::vector <double> getRad2Down(int profile, int channel);
    std::vector <double> getRad2Surf(int profile, int channel);

protected :

    Profiles::Profiles * profiles;
    string strOptions;
    int nprofiles;
    int nlevels;
    int nchannels; // number of channels loaded from coefficients
    int nchannelsForLastRun;  // number of channels used for the last run
    int ngases;
    int gas_units;
    int inst_id;
    double * p;//[nprofiles][nlevels];                  // Input pressure profiles
    double * t;//[nprofiles][nlevels];                  // Input temperature profiles
    double * gases;//[ngases][nprofiles][nlevels];      // Input gas profiles
    double * surfemisrefl;//[2][nprofiles][nchannels];  // Input/output surface emissivities and BRDFs
    double * btrefl;//[nprofiles][nchannels];           // Output BTs/refls (for thermal/solar chans)
    double * rads;//[nprofiles][nchannels];             // Output radiances

    double * tautotal;
    double * taulevels;
    double * tausuntotalpath1;
    double * tausunlevelspath1;
    double * tausuntotalpath2;
    double * tausunlevelspath2;

    double * radclear;
    double * radtotal;
    double * btclear;
    double * bt;
    double * reflclear;
    double * refl;
    double * radcloudy;
    double * overcast;

    double * rad2upclear;
    double * rad2dnclear;
    double * rad2refldnclear;
    double * rad2up;
    double * rad2down;
    double * rad2surf;

    int * gas_id;// {gas_id_q, gas_id_co2, gas_id_cfrac, gas_id_lwc1, gas_id_iwc};

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

    string file_coef;
    string file_sccld;
    string file_pccoef;
    string file_scaer;
    string emis_atlas_path;
    string brdf_atlas_path;

    bool coeffsLoaded;
    bool profileSet;
    std::vector<int> vector_channels;

    double * wavenumbers;
    double * refPressures;
    bool allocatedSurfemisrefl;

    double * skin_k; //skin T, salinity, snow_frac, foam_frac, fastem_coefsx5  [nprofiles][nchannels][9]
    double * s2m_k; //2m p, 2m t, 2m q, 10m wind u, v, wind-fetch    [nprofiles][nchannels][6]
    double * simplecloud_k ; // [nprofiles][nchannels][2]
    double * p_k; // [nprofiles][nchannels][nlevels]
    double * t_k; // [nprofiles][nchannels][nlevels]
    double * gases_k; // [ngases][nprofiles][nchannels][nlevels]
    double * surfemisrefl_k; //[2][nprofiles][nchannels]
    double * bt_k; //[nprofiles][nchannels]
    double * rad_k; //[nprofiles][nchannels]
    bool allocatedSurfemisreflk;
    bool allocatedBtk;
    bool allocatedRadk;
    std::vector <double> convertPointer4D2Vector(double  ptr[],int x, int y,int z, int dim1, int dim2, int dim3, int dim4);
    std::vector <double> convertPointer3D2Vector(double  ptr[],int x, int y, int dim1, int dim2, int dim3);
    std::vector <double> convertPointer2D2Vector(double  ptr[],  int x,  int dim2, int dim3);
    bool debug;
    bool isIrAtlasLoaded;
    bool isMwAtlasLoaded;
    bool isBrdfAtlasLoaded;
    bool allocatedP;

    void doStoreTrans();
    void doStoreRad();
    void doStoreRad2();
};
} /* namespace rttov */
#endif /* RTTOV_H_ */
