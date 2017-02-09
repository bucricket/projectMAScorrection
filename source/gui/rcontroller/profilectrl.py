# -*- coding: utf-8 -*-

from util import *

import rview
import wx
import logging


class ProfileController (GenericController):
    """ controller for the profile windows """

    def __init__(self, parent):
        """ Constructor of the ProfileController
            create the ProfileView Object
            make the binding with the profile view

        """
        self.theParentFrame = parent
        super(GenericController, self).__init__()
        # surface view
        self.profileView = None
        self.Redraw = True

        # subscribe to project changes
        pub.subscribe(self.ProfileChangedListener, "Profile CHANGED")
        self.controlerName = "ProfileController"

    def ProfileChangedListener(self, msg):
        if self.profileView:
            if self.Redraw:
                self.profileView.RePlotAll(self.project.myProfile)

    def OnProfile(self, e):
        """ create the profile view """
        if self.project.loadProfile:
            self.ShowProfile(self.project.myProfile)
        else:
            self.theParentFrame.WarmError("You must open a Profile ")

    def ShowProfile(self, profile):
        if not self.profileView:
            self.profileView = rview.profileframe.ProfileView(
                self.theParentFrame, profile, controler=self)
            self.MakeProfileBinding()
        else:
            if self.profileView.IsIconized():
                self.profileView.Restore()

    def applyProfileChange(self, e):
        """ apply the change to the profile """
        self.write("apply profile changes")
        aProfile = self.profileView.OnApplyChange(e)
        self.Redraw = False
        # we don't want to redraw everything in this case
        # (project.updateProfile send a message "Profile CHANGED")
        self.project.updateProfile(aProfile)
        self.Redraw = True

    def saveProfile(self, e):
        """ apply profile change and save profile+option """
        logging.debug("debug saveProfile surface controller")
        profile = self.profileView.OnApplyChange(e)
        self.Redraw = False
        self.project.updateProfile(profile)
        self.SaveProfile(e, self.profileView)
        self.Redraw = True

    def saveProfileAs(self, e):
        """ apply profile change and save profile+option with a new name """
        logging.debug("debug saveProfileAs surface controller")
        # re-initialize the name of the saved profile file to None
        self.project.savedProfileFileName = None
        self.saveProfile(e)

    def MakeProfileBinding(self):
        """ binding the binding for the profile view """
        self.profileView.Bind(
            wx.EVT_MENU, self.applyProfileChange,
            self.profileView.items['applyChange'])
        self.profileView.Bind(wx.EVT_MENU, self.saveProfileAs,
                              self.profileView.items['saveProfileAs'])
        self.profileView.Bind(wx.EVT_MENU, self.saveProfile,
                              self.profileView.items['saveProfile'])
