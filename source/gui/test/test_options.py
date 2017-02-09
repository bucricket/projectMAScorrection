'''
Created on May 7, 2014

@author: pascale
'''
import unittest
import rmodel
import logging
import glob
import rttovgui_unittest_class


class Test(rttovgui_unittest_class.RttovGuiUnitTest):

    def setUp(self):
        level_logging = logging.DEBUG
        self.p = rmodel.project.Project()

        logging.basicConfig(filename=(self.p.config.ENV['GUI_WRK_DIR'] +
                                      "/rttovgui_inittest.log"),
                            format=("[%(asctime)s] %(levelname)s "
                                    "[%(module)s:%(funcName)s:%"
                                    "(lineno)d] %(message)s"),
                            level=level_logging,
                            datefmt="%Y:%m:%d %H:%M:%S",
                            filemode="w")
        logging.info("start main controller")

        profileName = self.p.config.ENV[
            "RTTOV_GUI_PREFIX"] + "/rttov_tests/cldaer101lev_allgas.H5"
        self.p.openProfile(profileName)
        self.p.myProfile.display()


#    def tearDown(self):
#        pass

    def controleOptionsGas(self):

        for gas in ["CO2", "N2O", "CO", "CH4"]:
            print "controle gas", gas
            if self.p.myCoeffs.hasGas(gas) and self.p.myProfile.hasGas(gas):

                self.assertTrue(self.p.myOption.statusOption[gas + "_DATA"])
            else:
                self.assertFalse(self.p.myOption[gas + "_DATA"])
                self.assertFalse(self.p.myOption.statusOption[gas + "_DATA"])
        gas = "OZONE"
        if self.p.myCoeffs.hasGas("O3") and self.p.myProfile.hasGas("O3"):
            print "gas=OZONE"
            self.assertTrue(self.p.myOption.statusOption[gas + "_DATA"])
        else:
            print ("self.p.myProfile.hasGas('O3')",
                   self.p.myProfile.hasGas("O3"))
            print self.p.myProfile["O3"]
            print "self.p.myCoeff.hasGas('O3')", self.p.myCoeffs.hasGas("O3")
            self.assertFalse(self.p.myOption[gas + "_DATA"])
            self.assertFalse(self.p.myOption.statusOption[gas + "_DATA"])

    def controlePC(self):
        if (self.p.myOption["ADDPC"]):
            self.assertTrue(self.p.myOption.statusOption["ADDRADREC"])
            self.assertTrue(self.p.myOption.statusOption["IPCBND"])
            self.assertTrue(self.p.myOption.statusOption["IPCREG"])
        else:
            self.assertFalse(self.p.myOption.statusOption["ADDRADREC"])
            self.assertFalse(self.p.myOption.statusOption["IPCBND"])
            self.assertFalse(self.p.myOption.statusOption["IPCREG"])
            self.assertFalse(self.p.myOption["ADDRADREC"])
            self.assertFalse(self.p.myOption["ADDPC"])

    def test_solar1(self):
        self.p.myCoeffs.fileName["standard"] = self.p.config.ENV[
            "RTTOV_GUI_COEFF_DIR"] + "/rttov9pred54L/rtcoef_eos_2_modis.dat"
        err = self.p.loadCoefficients()

        self.p.ctrlCoherence()
        self.p.myOption.display()
        self.assertEqual(err, 0)
        self.assertTrue(self.p.myOption.statusOption["ADDSOLAR"])
        self.assertFalse(self.p.myOption["ADDSOLAR"])
        self.assertFalse(self.p.myOption["ADDCLOUDS"])
        self.assertFalse(self.p.myOption["ADDAEROSL"])
        print "OZONE?", self.p.myOption["OZONE_DATA"]
        self.controleOptionsGas()

        print "OZONE?", self.p.myOption["OZONE_DATA"]

        self.assertFalse(self.p.myOption["CO2_DATA"])
        self.assertFalse(self.p.myOption["N2O_DATA"])
        self.assertFalse(self.p.myOption["CO_DATA"])
        self.assertFalse(self.p.myOption.statusOption["ADDPC"])
        self.assertFalse(self.p.myOption["ADDPC"])
        print "OZONE?", self.p.myOption["OZONE_DATA"]
        self.p.myProfile.removeGas("O3")
        print "OZONE?", self.p.myOption["OZONE_DATA"]
        self.p.ctrlCoherence()
        print "OZONE?", self.p.myOption["OZONE_DATA"]
        self.assertFalse(self.p.myOption["OZONE_DATA"])
        self.assertFalse(self.p.myOption.statusOption["OZONE_DATA"])
        self.assertFalse(self.p.myProfile.hasGas("O3"))
        err = self.p.runDirect()
        self.assertEqual(err, 0)
        err = self.p.runK()
        self.assertEqual(err, 0)
        print "OZONE?", self.p.myOption["OZONE_DATA"]
        err = self.p.runDirect()
        self.assertEqual(err, 0)
        self.assertFalse(self.p.myOption["OZONE_DATA"])
        self.assertFalse(self.p.myOption.statusOption["OZONE_DATA"])
        err = self.p.dropCoefficients()
        self.assertEqual(err, 0)
        print "***************************************************test iasi *"
        self.p.myCoeffs.fileName["standard"] = self.p.config.ENV[
            "RTTOV_GUI_COEFF_DIR"] + "/rttov9pred101L/rtcoef_metop_2_iasi.H5"
        err = self.p.loadCoefficients()
        print "OZONE?", self.p.myOption["OZONE_DATA"]
        self.assertEqual(err, 0)
        self.p.ctrlCoherence()
        self.p.myOption.display()
        self.controleOptionsGas()
        self.assertFalse(self.p.myOption["OZONE_DATA"])
        self.assertFalse(self.p.myOption.statusOption["OZONE_DATA"])

        self.assertFalse(self.p.myOption["ADDSOLAR"])
        self.assertTrue(self.p.myOption.statusOption["ADDSOLAR"])

        self.assertFalse(self.p.myOption["ADDPC"])
        self.controlePC()

        print ("CH4?=", self.p.myOption["CH4_DATA"], "status?",
               self.p.myOption.statusOption["CH4_DATA"])
        self.p.myProfile.removeGas("CH4")

        self.p.ctrlCoherence()
        print ("CH4?=", self.p.myOption["CH4_DATA"], "status?",
               self.p.myOption.statusOption["CH4_DATA"])
        self.controleOptionsGas()

        self.assertFalse(self.p.myOption["CH4_DATA"])
        self.assertFalse(self.p.myOption.statusOption["CH4_DATA"])

        self.assertFalse(self.p.myOption["ADDPC"])
        self.assertTrue(self.p.myOption.statusOption["ADDSOLAR"])
        err = self.p.runDirect()
        self.assertEqual(err, 0)
        print ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>ok "
        self.p.myCoeffs.fileName["PC"] = self.p.config.ENV[
            "RTTOV_GUI_COEFF_DIR"] + "/pc/pccoef_metop_2_iasi.H5"
        err = self.p.loadCoefficients()
        self.assertEqual(err, 0)
        self.p.ctrlCoherence()

        self.assertTrue(self.p.myOption.statusOption["ADDPC"])
        self.controlePC()

        self.assertTrue(self.p.myOption.statusOption["ADDSOLAR"])

        self.assertFalse(self.p.myOption.statusOption["ADDCLOUDS"])
        self.assertFalse(self.p.myOption["ADDCLOUDS"])

        self.p.myOption["ADDPC"] = True
        self.p.ctrlCoherence()
        self.controlePC()

        self.assertFalse(self.p.myOption.statusOption["ADDSOLAR"])
        self.assertFalse(self.p.myOption["ADDSOLAR"])
        self.assertFalse(self.p.myOption.statusOption["ADDCLOUDS"])
        self.assertFalse(self.p.myOption["ADDCLOUDS"])
        err = self.p.runPC()
        self.assertEqual(err, 0)
        print">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>ok 1"
        self.p.myOption["ADDPC"] = False
        err = self.p.runDirect()
        print">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>ok 2"
        self.assertEqual(err, 0)
        self.p.ctrlCoherence()
        self.assertTrue(self.p.myOption.statusOption["ADDSOLAR"])
        self.controlePC()
        err = self.p.runDirect()
        self.assertEqual(err, 0)
        print">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>ok 3"
        self.p.myOption["ADDSOLAR"] = True
        self.p.ctrlCoherence()
        self.assertFalse(self.p.myOption.statusOption["ADDPC"])
        self.controlePC()

        err = self.p.runDirect()
        self.assertEqual(err, 0)
        print">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>ok 4"
        self.p.myOption["N2O_DATA"] = True
        self.p.ctrlCoherence()

        err = self.p.runDirect()
        self.assertEqual(err, 0)
        print">>>>>>>>>>>>>>>>>>>>>>>>ok 5"
        self.p.myOption["ADDPC"] = True

        self.p.ctrlCoherence()
        self.p.myOption["ADDSOLAR"] = True
        err = self.p.runPC()
        self.assertEqual(err, 0)

        print ">>>>>>>>>>>>>>>>>>>>>>>>>End test solar1"

    def test_removeO3(self):
        print ">>>>>>>>>>>>>>>>>>test_removeO3"

        self.p.myProfile.removeGas("O3")
        list_std = []
        for coefFile in glob.glob(self.p.config.ENV[
                            "RTTOV_GUI_COEFF_DIR"] + "/rttov7pred54L/*.h5"):
            list_std.append(coefFile)
        for coefFile in glob.glob(self.p.config.ENV[
                            "RTTOV_GUI_COEFF_DIR"] + "/rttov8pred54L/*.h5"):
            list_std.append(coefFile)
        for coefFile in glob.glob(self.p.config.ENV[
                            "RTTOV_GUI_COEFF_DIR"] + "/rttov9pred54L/*.h5"):
            list_std.append(coefFile)
        for coefFile in glob.glob(self.p.config.ENV[
                            "RTTOV_GUI_COEFF_DIR"] + "/rttov7pred101L/*.h5"):
            list_std.append(coefFile)
        for coefFile in glob.glob(self.p.config.ENV[
                            "RTTOV_GUI_COEFF_DIR"] + "/rttov8pred101L/*.h5"):
            list_std.append(coefFile)
        for coefFile in glob.glob(self.p.config.ENV[
                            "RTTOV_GUI_COEFF_DIR"] + "/rttov9pred101L/*.h5"):
            list_std.append(coefFile)
        for coefFile in glob.glob(self.p.config.ENV[
                            "RTTOV_GUI_COEFF_DIR"] + "/rttov7pred54L/*.dat"):
            list_std.append(coefFile)
        for coefFile in glob.glob(self.p.config.ENV[
                            "RTTOV_GUI_COEFF_DIR"] + "/rttov8pred54L/*.dat"):
            list_std.append(coefFile)
        for coefFile in glob.glob(self.p.config.ENV[
                            "RTTOV_GUI_COEFF_DIR"] + "/rttov9pred54L/*.dat"):
            list_std.append(coefFile)
        for coefFile in glob.glob(self.p.config.ENV[
                            "RTTOV_GUI_COEFF_DIR"] + "/rttov7pred101L/*.dat"):
            list_std.append(coefFile)
        for coefFile in glob.glob(self.p.config.ENV[
                            "RTTOV_GUI_COEFF_DIR"] + "/rttov8pred101L/*.dat"):
            list_std.append(coefFile)
        for coefFile in glob.glob(self.p.config.ENV[
                            "RTTOV_GUI_COEFF_DIR"] + "/rttov9pred101L/*.dat"):
            list_std.append(coefFile)

        for coefFile in list_std:
            print ">>>>>>>>>>>>>>>>>>>>>coefFile;", coefFile
            self.p.myCoeffs.fileName["standard"] = coefFile
            err = self.p.loadCoefficients()
            self.assertEqual(err, 0)
            print ">>>>>>>>>>>>>>>>>>>>>coefFile;", coefFile

            self.p.ctrlCoherence()
            self.check_option(self.p)
            err = self.p.runDirect()
            self.assertEqual(err, 0)
            err = self.p.runK()
            self.assertEqual(err, 0)
            self.check_option(self.p)
            print ">>>>>>>>>>>>>>>>>>>>>>OK", coefFile

        print ">>>>>>>>>>>>>>>>> end test remove O3"


if __name__ == "__main__":
    unittest.main()
