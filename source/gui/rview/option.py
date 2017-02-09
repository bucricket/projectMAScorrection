# -*- coding: utf-8 -*-
import wx
import rmodel
import util
import wx.lib.agw.floatspin as FS
import os
import logging


class OptionView(util.GenericView):
    """ Option Dialog window of the application """
    helpPage = os.environ["RTTOV_GUI_PREFIX"] + "/doc/helpOptions.html"
    helpTitle = "Options help"

    def __init__(self, parent, project, controller=None):
        self.project = project
        self.myOptions = project.myOption
        self.myControler = controller
        self.optionsThemes = self.myOptions.optionsThemes
        self.optionsThemesListe = self.myOptions.optionsThemesListe
        self.comboChoice = self.myOptions.comboChoice
        util.GenericView.__init__(self, parent,  "Options")

        self.CreateMenuBar()
        self.SetSize((570, 660))
        self.SetMinSize((570, 660))
        self.SetTitle('OPTIONS')

        self.panel1 = wx.Panel(self, -1)
        self.masterSizer = wx.BoxSizer(wx.VERTICAL)
        self.panel1.SetSizerAndFit(self.masterSizer)

        self.cb = {}
        self.combo = {}
        self.boxes = {}
        self.boxsizers = {}
        self.gb = wx.GridBagSizer(hgap=5, vgap=5)

        mysizer = {}
        for i in range(6):
            self.boxes[self.optionsThemes[i]] = wx.StaticBox(
                self.panel1, label=self.optionsThemes[i])
            self.boxsizers[self.optionsThemes[i]] = wx.StaticBoxSizer(
                self.boxes[self.optionsThemes[i]], wx.VERTICAL)

            for parametre in self.optionsThemesListe[self.optionsThemes[i]]:
                if parametre in ("IPCREG", 'IPCBND', 'FASTEM_VERSION',
                                 'INTERP_MODE'):
                    mysizer[parametre] = wx.BoxSizer(wx.HORIZONTAL)
                    self.combo[parametre] = wx.ComboBox(
                                                self.panel1,
                                                choices=self.comboChoice[
                                                    parametre].values(),
                                                size=(190, 26))
                    mysizer[parametre].Add(self.combo[parametre], border=5)
                    mysizer[parametre].Add(wx.StaticText(
                        self.panel1, -1, parametre.swapcase()), border=5)
                    self.boxsizers[self.optionsThemes[i]].Add(
                        mysizer[parametre], flag=wx.LEFT | wx.TOP | wx.EXPAND,
                        border=5)
                else:
                    if (parametre == "CLDSTR_THRESHOLD"):
                        self.thresholdItem = FS.FloatSpin(self.panel1, -1,
                                                          min_val=-1,
                                                          max_val=1,
                                                          increment=0.001,
                                                          agwStyle=FS.FS_LEFT)
                        self.thresholdItem.SetFormat("%f")
                        self.thresholdItem.SetDigits(3)
                        sizer2 = wx.BoxSizer(wx.HORIZONTAL)
                        sizer2.Add(wx.StaticText(
                                        self.panel1, -1,
                                        "CLDSTR_THRESHOLD".swapcase()),
                                   border=5)
                        sizer2.Add(self.thresholdItem, border=5)
                        self.boxsizers[self.optionsThemes[i]].Add(
                            sizer2, flag=wx.LEFT | wx.TOP, border=5)
                    else:
                        self.cb[parametre] = wx.CheckBox(
                            self.panel1, label=parametre.swapcase())
                        self.boxsizers[self.optionsThemes[i]].Add(
                            self.cb[parametre], flag=wx.LEFT | wx.TOP,
                            border=5)

            if i < 2:
                if i == 1:
                    self.gb.Add(self.boxsizers[self.optionsThemes[i]],
                                pos=(i, 0), span=(3, 1), border=20,
                                flag=(wx.EXPAND | wx.ALIGN_CENTER_HORIZONTAL |
                                      wx.ALIGN_CENTER_VERTICAL))
                else:
                    self.gb.Add(self.boxsizers[self.optionsThemes[i]],
                                pos=(i, 0), border=20,
                                flag=(wx.EXPAND | wx.ALIGN_CENTER_HORIZONTAL |
                                      wx.ALIGN_CENTER_VERTICAL))
            else:

                self.gb.Add(self.boxsizers[self.optionsThemes[i]],
                            pos=(i - 2, 1), border=20, flag=wx.EXPAND)

        self.masterSizer.Add(self.gb, flag=wx.ALIGN_CENTER)
        self.masterSizer.Add((10, 10), flag=wx.EXPAND)
        self._CreateButtons()

        self.masterSizer.Add(self.btnSizer, flag=wx.ALIGN_CENTER)
        self.SetValueItems()
        self.sb = self.CreateStatusBar()
        self.Centre()
        self.Show(True)

        self.Bind(wx.EVT_CLOSE, self.OnClose)
        self.combo['IPCBND'].Bind(wx.EVT_COMBOBOX, self.UpdateIPCREG)
        self.cb['ADDPC'].Bind(wx.EVT_CHECKBOX, self.UpdatePC)
        self.cb['ADDRADREC'].Bind(wx.EVT_CHECKBOX, self.UpdatePCaddradrec)
        self.cb['ADDCLOUDS'].Bind(wx.EVT_CHECKBOX, self.UpdateClouds)
        self.UpdatePC(None)

        # hide none yet used option
        self.cb["CLDSTR_SIMPLE"].Hide()
        self.cb["USER_AER_OPT_PARAM"].Hide()
        self.cb["USER_CLD_OPT_PARAM"].Hide()

    def UpdateClouds(self, e):
        if self.cb["ADDCLOUDS"].GetValue():
            if not self.cb["ADDPC"].GetValue():
                self.cb["CLDSTR_SIMPLE"].Enable()
            else:
                self.cb["CLDSTR_SIMPLE"].Disable()
                self.cb["CLDSTR_SIMPLE"].SetValue(False)
            self.thresholdItem.Enable()
        else:
            self.cb["CLDSTR_SIMPLE"].Disable()
            self.thresholdItem.Disable()
            self.cb["CLDSTR_SIMPLE"].SetValue(False)

    def UpdatePCaddradrec(self, e):
        if self.cb["ADDRADREC"].GetValue():
            self.cb["SWITCHRAD"].Enable()
        else:
            self.cb["SWITCHRAD"].Disable()
            self.cb["SWITCHRAD"].SetValue(False)

    def UpdatePC(self, e):
        """ enable or disable pc options if addpc is checked or not """
        self.GetValuesItems()
        self.project.ctrlCoherence()
        self.SetValueItems()

    def UpdateIPCREG(self, e):
        bandValue = self.combo["IPCBND"].GetValue()
        band = [k for k, v in self.comboChoice["IPCBND"].iteritems()
                if v == bandValue][0]
        maxIPCREG = self.project.myCoeffs.getMaxIPCREG(band)
        self.combo["IPCREG"].SetItems(
            map(lambda x: str(x), range(1, maxIPCREG + 1)))
        myval = self.comboChoice['IPCREG'][self.myOptions['IPCREG']]
        self.combo["IPCREG"].SetValue(myval)

    def _MakeBinding(self):
        """ set the trivial Binding for the View """
        self.Bind(wx.EVT_BUTTON, self.OnCancel, self.cancelBtn)

    def OnItemFocus(self, e):  # TODO AREVOIR
        for (i, value) in self.cb.items():
            if (value == e.GetEventObject()):
                self.sb.PushStatusText(
                    self.myOptions[i + '_ATTRIBUTE'][u'COMMENT'], 1)

    def OnItemUnFocus(self, e):  # TODO AREVOIR
        for (i, value) in self.cb.items():
            if (value == e.GetEventObject()):
                self.sb.PushStatusText(
                    self.myOptions[i + '_ATTRIBUTE'][u'COMMENT'], 1)

    def _CreateButtons(self):

        self.btnSizer = wx.BoxSizer(wx.HORIZONTAL)
        self.cancelBtn = wx.Button(self.panel1, wx.ID_CANCEL, label="Revert")
        self.cancelBtn.SetHelpText("Revert to previous options")
        self.btnSizer.Add(
            self.cancelBtn, flag=wx.ALIGN_RIGHT | wx.RIGHT, border=10)

        self.applyBtn = wx.Button(self.panel1, wx.ID_OK, label="Apply")
        self.applyBtn.SetHelpText("Apply options")
        self.applyBtn.SetDefault()
        self.btnSizer.Add(self.applyBtn, flag=wx.ALIGN_RIGHT |
                          wx.RIGHT, border=10)

        # binding cancel button
        self.cancelBtn.Bind(wx.EVT_BUTTON, self.OnCancel)

    def SetValueItems(self):
        for i in self.myOptions.options_list_logical:
            logging.debug(str(
                i) + " " + str(self.myOptions.statusOption[i]) + " " +
                               str(self.cb[i].GetValue()))
            self.cb[i].Enable(self.myOptions.statusOption[i])
            self.cb[i].SetValue(self.myOptions[i])

        for param in ['IPCREG', 'FASTEM_VERSION', "INTERP_MODE", "IPCBND"]:
            logging.debug("SetValueItems " + param +
                          str(self.myOptions[param]) +
                          str(self.comboChoice[param]))

            myval = self.comboChoice[param][self.myOptions[param]]
            logging.debug("SetValueItems " + param +
                          str(self.myOptions[param]) + " myval " + str(myval))
            self.combo[param].Enable(self.myOptions.statusOption[param])
            self.combo[param].SetItems(
                self.myOptions.comboChoice[param].values())

            self.combo[param].SetValue(myval)

        self.thresholdItem.SetValue(self.myOptions['CLDSTR_THRESHOLD'])
        self.thresholdItem.Enable(
            self.myOptions.statusOption['CLDSTR_THRESHOLD'])

        if self.project.myCoeffs.loadCoeffs:
            if self.project.myCoeffs.hasPC():
                self.UpdateIPCREG(None)

    def SetOptions(self, options):
        self.myOptions = options
        self.SetValueItems()

    def GetOptions(self):
        return self.myOptions

    def GetValuesItems(self):
        for i in self.myOptions.options_list_logical:
            self.myOptions[i] = self.cb[i].GetValue()
        for i in ['IPCREG', 'FASTEM_VERSION', 'IPCBND', 'INTERP_MODE']:
            value = self.combo[i].GetValue()
            self.myOptions[i] = [k for k, v in self.comboChoice[i].iteritems()
                                 if v == value][0]

            logging.debug("value" + str(value) +
                          " result " + str(self.myOptions[i]))

        try:
            self.myOptions['CLDSTR_THRESHOLD'] = float(
                self.thresholdItem.GetValue())
        except ValueError:
            self.ShowErrorMessageDialogBox(
                var='CLDSTR_THRESHOLD', type='float')


#

    def ShowErrorMessageDialogBox(self, varName, varType):
        message = "variable " + varName + " must be of " + varType + " type"
        dlg = wx.MessageDialog(
            None, message, caption="Error", style=wx.ICON_ERROR)
        dlg.ShowModal()
        dlg.Destroy()

    def OnCancel(self, e):
        """ Cancel modifications made in the frame
            set all widget with initial option values"""
        self.SetValueItems()

    def OnApply(self, e):
        """ take value from the windows and save it in options"""
        self.GetValuesItems()

    def OnSave(self, e):
        """ take value from the windows and save it in options"""
        self.GetValuesItems()

    def OnClose(self, e):
        """ close the option windows"""
        if self.myControler is not None:
            self.myControler.optionView = None
        self.Destroy()

    def _initItem(self):
        self.applyItem = 0
        self.saveItem = 0

    def MenuData(self):
        """ define the data for the menu
        """
        return(("&File",  # File Menu
                ("Apply options", "Apply the options",
                 self.OnApply, "applyOptions", True),
                ("Save options",
                 "Apply and Save the options and the profile in a file",
                 self.OnSave, "saveOptions", True),
                ("", "", "", "", True),
                ('&Quit', 'Quit', self.OnQuit, "quit", True)),
               ("&Help",  # Help Menu
                ("About", "About screen", self.OnAbout, "about", True),
                ("&Help", "Help", self.OnHelpHTML, "help", True)))


if __name__ == "__main__":
    p = rmodel.project.Project()
    p.openProfile(p.config.ENV["RTTOV_GUI_PREFIX"] +
                  '/rttov_tests/cldaer101lev_allgas.H5')
    for option in p.myOption.statusOption.keys():
        p.myOption.statusOption[option] = True

    coefFile = p.config.ENV['RTTOV_GUI_COEFF_DIR'] + \
        "/rttov9pred101L/rtcoef_metop_2_iasi.H5"
    pcFile = p.config.ENV['RTTOV_GUI_COEFF_DIR'] + "/pc/pccoef_metop_2_iasi.H5"
    p.myCoeffs.fileName["standard"] = coefFile
    p.myCoeffs.fileName["PC"] = pcFile

    p.loadCoefficients()

    ex = wx.App()
    p.myOption.display()
    mv = OptionView(None, p)

    ex.MainLoop()
