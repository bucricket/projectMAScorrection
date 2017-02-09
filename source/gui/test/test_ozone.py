'''
Created on Nov 25, 2013

@author: pascale
'''
import rmodel
import rttov
import h5py
import sys
import os
import rttov_gui_f2py
import logging
import numpy
import unittest
import rttovgui_unittest_class


class Test(rttovgui_unittest_class.RttovGuiUnitTest):

    def test_ozone(self):

        p = rmodel.project.Project()
        p.debug = True
        print "Configuration : ", p.config.ENV['RTTOV_GUI_PREFIX']
        print "open profile"
        # profileName=p.config.ENV["RTTOV_GUI_PREFIX"]+"/PROFILES/cld101lev_allgas.H5"
        profileName = p.config.ENV[
            "RTTOV_GUI_PREFIX"] + "/PROFILES/standard54lev_allgas.H5"
        print profileName
        err = p.openProfile(profileName)
        self.assertEqual(err, 0)
        # p.myProfile.display()
        print "save profile"
        print "load coefficient"
        coefFile = p.config.ENV['RTTOV_GUI_COEFF_DIR'] + \
            "/rttov7pred54L/rtcoef_dmsp_18_ssmis.dat"

        print coefFile

        p.myCoeffs.fileName["standard"] = coefFile
        err = p.loadCoefficients()
        self.assertEqual(err, 0)
        self.assertTrue(p.myCoeffs.isMW())
        self.assertFalse(p.myCoeffs.hasPC())
        self.assertFalse(p.myCoeffs.isHighResolution())
        self.assertTrue(p.myCoeffs.isMW())
        self.assertFalse(p.myCoeffs.hasSolar())
        print ">>>>>>>>>>>>>>>>FMV_GAS_POS>>>>>>>>>>>>>>>"
        print p.myCoeffs.getFMV_GAS_POS()
        print ">>>>>>>>>>>>>>>>FMV_GAS_POS>>>>>>>>>>>>>>>"
        self.assertFalse(p.myCoeffs.hasGas("CO2"))
        self.assertFalse(p.myCoeffs.hasGas("O3"))
        self.assertFalse(p.myCoeffs.hasGas("N2O"))
        self.assertFalse(p.myCoeffs.hasGas("CO"))
        self.assertFalse(p.myCoeffs.hasGas("CH4"))

        # p.myProfile.display()
        p.ctrlCoherence()
        self.check_option(p)

        err = p.runDirect()

        self.assertEqual(err, 0)
        profileName = p.config.ENV[
            "RTTOV_GUI_PREFIX"] + "/PROFILES/standard54lev_allgas.H5"
        print ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>", profileName
        p.openProfile(profileName, 2)
        print ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>before save profile "
        p.saveProfile("/tmp/toto.h5")
        print ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>after save profile "
#        err=p.runDirect()
#
#        self.assertEqual(err,0)

if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.test_read_profile']
    unittest.main()
