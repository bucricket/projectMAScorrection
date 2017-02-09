# -*- coding: utf-8 -*-

from util import *
import wx
import rview


class OptionController (GenericController):
    """ controller for the Option windows """

    def __init__(self, parent):
        """ Constructor of the MainController
            define a Project object
            create the OptionView Object
            make the binding with the options view

        """
        self.theParentFrame = parent
        super(GenericController, self).__init__()
        self.optionView = None
        pub.subscribe(self.OptionsChangedListener, "Options CHANGED")
        self.controlerName = "OptionController"

    def OptionsChangedListener(self, msg):
        """ update the Options window if it exist """
        # print "OptionsChangedListener"
        if self.optionView is not None:
            self.optionView.SetOptions(self.project.myOption)

    def MakeBinding(self):
        """  Define the different binding of the application  """
        # binding for the MaineView (mv)

        self.mv.Bind(wx.EVT_MENU, self.OnOptions,
                     self.mv.items["optionsWindow"])

    def _MakeBindingOption(self):

        self.optionView.Bind(wx.EVT_MENU, self.applyOptions,
                             self.optionView.items['applyOptions'])
        self.optionView.Bind(
            wx.EVT_BUTTON, self.applyOptions, self.optionView.applyBtn)

    def OnOptions(self, e):

        if self.project.loadProfile:
            self.project.ctrlCoherence()
            self.ShowOptions(self.project)
        else:
            self.theParentFrame.WarmError(
                "You must open a Profile or create a new Profile ")

    def applyOptions(self, e):
        """ apply the options """
        # print self.optionView.myOptions["IPCREG"]
        self.optionView.OnApply(e)
        self.project.ctrlCoherence()
        self.project.updateOptions(self.optionView.myOptions)
        self.write("options applied")
        # print self.optionView.myOptions["IPCREG"]

    def ShowOptions(self, project):
        if not self.optionView:
            self.optionView = rview.option.OptionView(
                self.theParentFrame, project, controller=self)
            self._MakeBindingOption()
        else:
            if self.optionView.IsIconized():
                self.optionView.Restore()
