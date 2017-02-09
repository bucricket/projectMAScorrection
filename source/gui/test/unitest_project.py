# -*- coding: utf-8 -*-
import unittest
import rttov
import h5py
import rmodel
import rttovgui_unittest_class


class Test(rttovgui_unittest_class.RttovGuiUnitTest):

    def setUp(self):
        file = "../rttov_tests/cldaer101lev_allgas.H5"
        self.p = rmodel.project.Project()
        profileName = self.p.config.ENV[
            "RTTOV_GUI_PROFILE_DIR"] + "/standard54lev_allgas.H5"
        self.p.openProfile(profileName)

    def tearDown(self):
        pass

    def test_read_profile(self):
        self.assertGreater(self.p.myProfile['BE'],  0.20)

    def test_run_pc_1(self):
        coefFile = self.p.config.ENV[
            'RTTOV_GUI_COEFF_DIR'] + "/rttov9pred101L/rtcoef_metop_2_iasi.H5"
        pcFile = self.p.config.ENV[
            'RTTOV_GUI_COEFF_DIR'] + "/pc/pccoef_metop_2_iasi.H5"
        self.p.myCoeffs.fileName["standard"] = coefFile
        self.p.myCoeffs.fileName["PC"] = pcFile
        err = self.p.loadCoefficients()
        self.assertEqual(err, 0)
        self.assertFalse(self.p.isPC())
        self.p.myOption["ADDPC"] = True
        self.assertTrue(self.p.isPC())

        self.p.myProfile["SKIN"]["SURFTYPE"] = 1
        err = self.p.runPC()
        self.assertEqual(err, 0)

    def test_run_pc_2(self):
        print "-------------------AIRS-------------------------"
        self.p.myProfile["SKIN"]["SURFTYPE"] = 1
        coefFile = self.p.config.ENV[
            'RTTOV_GUI_COEFF_DIR'] + "/rttov9pred101L/rtcoef_eos_2_airs.H5"
        pcFile = self.p.config.ENV[
            'RTTOV_GUI_COEFF_DIR'] + "/pc/pccoef_eos_2_airs.H5"
        self.p.myCoeffs.fileName["standard"] = coefFile
        self.p.myCoeffs.fileName["PC"] = pcFile
        print self.p.myCoeffs.fileName
        err = self.p.loadCoefficients()
        self.assertEqual(err, 0)
        err = self.p.runPC()
        self.assertEqual(err, 0)
        self.check_option(self.p)
        print "-------------------IASI NG-------------------------"
        coefFile = self.p.config.ENV[
         'RTTOV_GUI_COEFF_DIR'] + "/rttov9pred101L/rtcoef_metopsg_1_iasing.H5"
        pcFile = self.p.config.ENV[
            'RTTOV_GUI_COEFF_DIR'] + "/pc/pccoef_metopsg_1_iasing.H5"
        self.p.myCoeffs.fileName["standard"] = coefFile
        self.p.myCoeffs.fileName["PC"] = pcFile
        print self.p.myCoeffs.fileName
        err = self.p.loadCoefficients()
        self.assertEqual(err, 0)
        err = self.p.runPC()
        self.assertEqual(err, 0)
        self.check_option(self.p)

    def test_profile_without_option(self):
        print ">>>>>>>>>>>>>>>>>>>test_profile_without_option"
        profileName = self.p.config.ENV[
            "RTTOV_GUI_PREFIX"] + "/rttov_tests/profile_without_options.h5"
        err = self.p.openProfile(profileName)
        self.assertFalse(self.p.myOption["CO2_DATA"])

    def test_read_profile_ascii_1(self):
        print ">>>>>>>>test_read_profile_ascii_1"
        profileName = self.p.config.ENV[
            "RTTOV_GUI_PREFIX"] + "/rttov_tests/p_us76.py"
        err = self.p.openAsciiProfile(profileName)
        self.assertEqual(err, 0)

    def test_read_profile_ascii_2(self):
        print ">>>>>>>>test_read_profile_ascii_2"
        profileName = self.p.config.ENV[
            "RTTOV_GUI_PREFIX"] + "/rttov_tests/ECMWF83_prof83.py"
        err = self.p.openAsciiProfile(profileName)
        self.assertEqual(err, 0)

    def test_read_profile_ascii_3(self):
        print ">>>>>>>>test_read_profile_ascii_3"
        profileName = self.p.config.ENV[
            "RTTOV_GUI_PREFIX"] + "/rttov_tests/profile_without_options.h5"
        err = self.p.openAsciiProfile(profileName)
        self.assertEqual(err, 1)

    def test_pc(self):
        print ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>test pc"
        p = rmodel.project.Project()
        profileName = p.config.ENV[
            "RTTOV_GUI_PROFILE_DIR"] + "/standard54lev_allgas.H5"
        print profileName
        err = p.openProfile(profileName)
        self.assertEqual(err, 0)
        print "---------------------TEST PC ---------------------------------"
        coefFile = p.config.ENV['RTTOV_GUI_COEFF_DIR'] + \
            "/rttov9pred101L/rtcoef_metop_2_iasi.H5"
        pcFile = p.config.ENV['RTTOV_GUI_COEFF_DIR'] + \
            "/pc/pccoef_metop_2_iasi.H5"
        p.myCoeffs.fileName["standard"] = coefFile
        p.myCoeffs.fileName["PC"] = pcFile
        print p.myCoeffs.fileName
        err = p.loadCoefficients()
        self.assertEqual(err, 0)
        print "isPC ?", p.isPC()
        self.assertFalse(p.isPC())
        p.myOption["ADDPC"] = True
        p.ctrlCoherence()
        self.check_option(p)
        print "isPC ?", p.isPC()
        self.assertTrue(p.isPC())
        p.myProfile['SKIN']['SURFTYPE'] = 1
        print "-------------------runPC------------------------"
        p.myOption["ADDRADREC"] = True
        p.myOption["ADDSOLAR"] = True
        p.myOption.display()
        err = p.runPC()
        self.assertEqual(err, 0)
        self.assertFalse(p.myOption["ADDSOLAR"])
        self.check_option(p)
        print "-------------------AIRS-------------------------"
        p.dropCoefficients()
        coefFile = p.config.ENV['RTTOV_GUI_COEFF_DIR'] + \
            "/rttov9pred101L/rtcoef_eos_2_airs.H5"
        pcFile = p.config.ENV['RTTOV_GUI_COEFF_DIR'] + \
            "/pc/pccoef_eos_2_airs.H5"
        p.myCoeffs.fileName["standard"] = coefFile
        p.myCoeffs.fileName["PC"] = pcFile
        print p.myCoeffs.fileName
        err = p.loadCoefficients()
        self.assertEqual(err, 0)
        err = p.runPC()
        self.assertEqual(err, 0)
        err = p.runPCK()
        self.assertEqual(err, 0)
        self.check_option(p)
        print "-------------------IASI NG-------------------------"
        coefFile = p.config.ENV['RTTOV_GUI_COEFF_DIR'] + \
            "/rttov9pred101L/rtcoef_metopsg_1_iasing.H5"
        pcFile = p.config.ENV['RTTOV_GUI_COEFF_DIR'] + \
            "/pc/pccoef_metopsg_1_iasing.H5"
        p.myCoeffs.fileName["standard"] = coefFile
        p.myCoeffs.fileName["PC"] = pcFile
        print p.myCoeffs.fileName
        err = p.loadCoefficients()
        self.assertEqual(err, 0)
        err = p.runPC()
        self.assertEqual(err, 0)
        err = p.runPCK()
        self.assertEqual(err, 0)
        self.check_option(p)


if __name__ == "__main__":
    unittest.main()
