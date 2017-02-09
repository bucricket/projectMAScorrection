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

/// @class RttovSafe
///  This class contains a set of functions for running Rttov.
///  An instance of a RttovSafe Class is defined to handle one instrument.
///  It contains an Options object which allows the Rttov options to
///  be specified and can be associated with a vector of Profile objects.
///  Profile objects are a safe way to initialise the profile values
///  therefore the RttovSafe presents a safe access to Rttov.

#ifndef RTTOVSAFE_H_
#define RTTOVSAFE_H_

#include <Rttov.h>
#include <Profile.h>
namespace rttov {

class RttovSafe: private Rttov {
  public:
    RttovSafe();
    virtual ~RttovSafe();

    Options::Options &options;

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
    int getCoeffsNlevels();

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

    void setTheProfiles(std::vector <rttov::Profile>& theProfiles);
    ItemIdIndexMap gas_index;
    void printGases();

    void runDirect(); // run direct for all loaded channels
    void runDirect(const vector<int>& channels);
    void runK();
    void runK(const vector<int>& channels);

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
};

} /* namespace rttov */

#endif /* RTTOVSAFE_H_ */
