'''
Created on May 14, 2014

@author: pascale
'''
import unittest
import rmodel


class RttovGuiUnitTest(unittest.TestCase):

    def test_dummy(self):
        pass

    def check_option(self, p):

        if p is None:
            return
        if p.isPC():
            self.assertTrue(p.myOption["ADDPC"])
            self.assertTrue(p.myOption.statusOption["ADDPC"])
        else:
            self.assertFalse(p.myOption["ADDPC"])
            self.assertTrue(p.myOption.statusOption["DO_LAMBERTIAN"])

        if p.myCoeffs.hasPC():
            self.assertTrue(p.myOption.statusOption["ADDPC"])
        else:
            self.assertFalse(p.myOption.statusOption["ADDPC"])
            self.assertFalse(p.myOption["ADDPC"])

        if p.myCoeffs.isPCClouds() and p.myProfile.hasClouds():
            self.assertTrue(p.myOption.statusOption["ADDPC"])

        if p.myCoeffs.hasSolar():
            if not p.isPC():
                self.assertTrue(p.myOption.statusOption["ADDSOLAR"])
        else:
            self.assertFalse(p.myOption.statusOption["ADDSOLAR"])
            self.assertFalse(p.myOption["ADDSOLAR"])

        if not p.myCoeffs.isMW():
            self.assertFalse(p.myOption.statusOption["FASTEM_VERSION"])
            self.assertFalse(p.myOption.statusOption["CLW_DATA"])
            self.assertFalse(p.myOption["CLW_DATA"])

        if p.myCoeffs.isMW():
            self.assertTrue(p.myOption.statusOption["FASTEM_VERSION"])
            self.assertTrue(p.myOption.statusOption["SUPPLY_FOAM_FRACTION"])

            if p.myProfile["CLW"] is not None:
                self.assertTrue(p.myOption.statusOption["CLW_DATA"])
            self.assertTrue(p.myOption.statusOption["DO_LAMBERTIAN"])
            self.assertFalse(p.myOption.statusOption["ADDSOLAR"])
            self.assertFalse(p.myOption.statusOption["DO_NLTE_CORRECTION"])
            self.assertFalse(p.myOption.statusOption["ADDAEROSL"])
            self.assertFalse(p.myOption.statusOption["ADDCLOUDS"])
            self.assertFalse(p.myOption.statusOption["CLDSTR_THRESHOLD"])
            self.assertFalse(p.myOption.statusOption["OZONE_DATA"])
            self.assertFalse(p.myOption.statusOption["CO2_DATA"])
            self.assertFalse(p.myOption.statusOption["N2O_DATA"])
            self.assertFalse(p.myOption.statusOption["CO_DATA"])
            self.assertFalse(p.myOption.statusOption["CH4_DATA"])
            self.assertFalse(p.myOption.statusOption["ADDPC"])
            self.assertFalse(p.myOption.statusOption["IPCBND"])
            self.assertFalse(p.myOption.statusOption["IPCREG"])
            self.assertFalse(p.myOption.statusOption["ADDRADREC"])

            self.assertFalse(p.myOption["ADDSOLAR"])
            self.assertFalse(p.myOption["DO_NLTE_CORRECTION"])
            self.assertFalse(p.myOption["ADDAEROSL"])
            self.assertFalse(p.myOption["ADDCLOUDS"])

            self.assertFalse(p.myOption["OZONE_DATA"])
            self.assertFalse(p.myOption["CO2_DATA"])
            self.assertFalse(p.myOption["N2O_DATA"])
            self.assertFalse(p.myOption["CO_DATA"])
            self.assertFalse(p.myOption["CH4_DATA"])
            self.assertFalse(p.myOption["ADDPC"])
            self.assertFalse(p.myOption["ADDRADREC"])

        else:  # not MW
            self.assertFalse(p.myOption.statusOption["FASTEM_VERSION"])
            self.assertFalse(p.myOption.statusOption["CLW_DATA"])
            self.assertFalse(p.myOption["CLW_DATA"])

            if p.isPC():  # not MW not PC
                for gas in ("CO", "CO2", "N2O", "CH4"):
                    if p.myCoeffs.hasGas(gas) and p.myProfile.hasGas(gas):
                        self.assertFalse(p.myOption[gas + "_DATA"])
                        self.assertFalse(
                            p.myOption.statusOption[gas + "_DATA"])
                self.assertFalse(p.myOption["ADDAEROSL"])
                self.assertFalse(p.myOption.statusOption["ADDAEROSL"])
                self.assertFalse(p.myOption["DO_NLTE_CORRECTION"])
                self.assertFalse(p.myOption.statusOption["DO_NLTE_CORRECTION"])
            else:
                for gas in ("CO", "CO2", "N2O", "CH4"):
                    if p.myCoeffs.hasGas(gas) and p.myProfile.hasGas(gas):
                        self.assertTrue(p.myOption.statusOption[gas + "_DATA"])
                    else:
                        self.assertFalse(p.myOption[gas + "_DATA"])
                        self.assertFalse(
                            p.myOption.statusOption[gas + "_DATA"])
                if p.myCoeffs.hasAerosols() and p.myProfile.hasAerosols():
                    self.assertTrue(p.myOption.statusOption["ADDAEROSL"])
                else:
                    self.assertFalse(p.myOption["ADDAEROSL"])
                    self.assertFalse(p.myOption.statusOption["ADDAEROSL"])

            if p.myCoeffs.hasGas("O3") and p.myProfile.hasGas("O3"):
                self.assertTrue(p.myOption.statusOption["OZONE_DATA"])
            else:
                self.assertFalse(p.myOption["OZONE_DATA"])
                self.assertFalse(p.myOption.statusOption["OZONE_DATA"])
        if not p.myCoeffs.hasNLTE():
            self.assertFalse(p.myOption.statusOption["DO_NLTE_CORRECTION"])
            self.assertFalse(p.myOption["DO_NLTE_CORRECTION"])
        if p.myOption["ADDCLOUDS"]:
            self.assertTrue(p.myOption.statusOption["CLDSTR_THRESHOLD"])
            if p.myOption["ADDPC"]:
                self.assertFalse(p.myOption.statusOption["CLDSTR_SIMPLE"])
                self.assertFalse(p.myOption["CLDSTR_SIMPLE"])
            else:
                self.assertTrue(p.myOption.statusOption["CLDSTR_SIMPLE"])
        else:
            self.assertFalse(p.myOption.statusOption["CLDSTR_SIMPLE"])
            self.assertFalse(p.myOption["CLDSTR_SIMPLE"])
            self.assertFalse(p.myOption.statusOption["CLDSTR_THRESHOLD"])
        if not p.myOption.statusOption["ADDCLOUDS"]:
            self.assertFalse(p.myOption.statusOption["CLDSTR_SIMPLE"])
            self.assertFalse(p.myOption["CLDSTR_SIMPLE"])
            self.assertFalse(p.myOption.statusOption["CLDSTR_THRESHOLD"])


class test_test(RttovGuiUnitTest):

    def test_check_option(self):
        print "test rttovgui_unittest_class"
        profileName = "../rttov_tests/cldaer101lev_allgas.H5"
        p = rmodel.project.Project()
        profileName = p.config.ENV[
            "RTTOV_GUI_PROFILE_DIR"] + "/us76_43lev_allgas.H5"
        print "openProfile", profileName
        p.openProfile(profileName)

        coefFile = p.config.ENV['RTTOV_GUI_COEFF_DIR'] + \
            "/rttov9pred101L/rtcoef_metop_2_iasi.H5"
        pcFile = p.config.ENV['RTTOV_GUI_COEFF_DIR'] + \
            "/pc/pccoef_metop_2_iasi.H5"
        p.myCoeffs.fileName["standard"] = coefFile
        p.myCoeffs.fileName["PC"] = pcFile

        p.loadCoefficients()
        p.ctrlCoherence()
        self.check_option(p)

if __name__ == "__main__":
    unittest.main()
