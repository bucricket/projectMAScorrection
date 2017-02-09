# -*- coding: utf-8 -*-

from util import *
import rview
import wx
import logging


class SurfaceController (GenericController):
    """ controller for the surface windows """

    def __init__(self, parent):
        """ Constructor of the SurfaceController
            create the SurfaceView Object
            make the binding with the surface view

        """
        self.theParentFrame = parent
        super(GenericController, self).__init__()
        # surface view
        self.surfaceView = None
        self.emisController = None
        self.reflController = None
        # subscribe to project changes
        pub.subscribe(self.ProfileLoadedListener, "Profile CHANGED")
        pub.subscribe(self.CoeffsLoadedListener, "Coefficients CHANGED")
        pub.subscribe(self.SurfaceChangedListener, "Surface CHANGED")
        self.controlerName = "SurfaceController"

    def SurfaceChangedListener(self, msg):
        logging.debug("SurfaceChangedListener")
        self.prepLog()
        if self.surfaceView:
            logging.debug("plotReflectanceEmissivity")
            self.surfaceView.SetProfileValues(self.project.myProfile)
            self.surfaceView.PlotReflectanceEmissivity(
                self.project.myEmissivity, self.project.myReflectance,
                self.project.myCoeffs.WAVENUMBERS)
            logging.debug(str(self.project.myReflectance))

    def ProfileLoadedListener(self, msg):
        """ update window of the surface window if it exists """
        self.prepLog()
        logging.debug("SurfaceController ProfileLoadedListener")
        if self.surfaceView:
            self.surfaceView.SetProfileValues(self.project.myProfile)

        logging.debug("useAtlas :" + str(self.project.useAtlas))
        needOpenAtlas = False
        if self.project.useAtlas:
            for param in self.project.ListeParamAtlas:
                if param == "DATE":
                    for i in range(3):
                        if self.project.paramAtlas[
                                param][i] != self.project.myProfile[param][i]:
                            needOpenAtlas = True
                elif self.project.paramAtlas[
                                param] != self.project.myProfile[param]:
                    needOpenAtlas = True
            self.write("Need Re-Open Atlas : " + str(needOpenAtlas))
            logging.debug("need to re-open Atlas")
            if needOpenAtlas:
                answer = self.theParentFrame.WarmQuestion(
                    "Profile has changed : do you want to re-open Atlas ?")
                if answer == wx.ID_YES:
                    self.project.openAtlas()

    def CoeffsLoadedListener(self, msg):
        """ update menus of the surface window if it exists """
        self.prepLog()
        if self.surfaceView:
            self.UpdateMenus()

    def UpdateMenus(self):
        """ allow all the menus items """
        self.surfaceView.EnableMenuItem('createEmissivityReflectance')
        self.surfaceView.EnableMenuItem('loadAtlas')
        self.surfaceView.EnableMenuItem('saveSurface')
        self.surfaceView.EnableMenuItem('saveSurfaceAs')
        self.surfaceView.EnableMenuItem('modifyEmissivity')
        self.surfaceView.EnableMenuItem('modifyReflectance')

    def UpdateMenusItemModifyEmissivityReflectance(self):
        pass

    def OnSurface(self, e):
        """ create the surface view """
        logging.debug("surfacectrl OnSurface")
        self.prepLog()
        if self.project.loadProfile:
            self.ShowSurface(self.project)

            if self.project.myCoeffs.loadCoeffs:
                self.UpdateMenus()
            if self.project.myEmissivity["EMIS_IN"] is not None:
                self.surfaceView.PlotReflectanceEmissivity(
                    self.project.myEmissivity, self.project.myReflectance,
                    self.project.myCoeffs.WAVENUMBERS)
        else:
            self.theParentFrame.WarmError("You must open a Profile ")

    def ShowSurface(self, project):
        logging.debug("surfacectrl showSurface")
        self.prepLog()
        if not self.surfaceView:
            self.surfaceView = rview.surface.SurfaceView(
                self.theParentFrame, project)
            if project.myEmissivity['EMIS_IN'] is not None:
                self.UpdateMenusItemModifyEmissivityReflectance()
            # create the 2 controllers for emissivity and reflectance edition
            self.emisController = ReflEmisController(self.surfaceView, "EMIS")
            self.reflController = ReflEmisController(self.surfaceView, "REFL")
            self.MakeSurfaceBinding()
        else:
            if self.surfaceView.IsIconized():
                self.surfaceView.Restore()

    def applyProfileChange(self, e):
        """ apply the change to the profile """
        self.prepLog()
        aProfile = self.surfaceView.OnApplyChange(e)
        self.project.updateProfile(aProfile)

    def openEmissivityReflectance(self, e):
        """ open an EmissityReflectance file """
        self.prepLog()
        surfaceFileName = self.surfaceView.OnOpenFile(
            e, self.project.config.ENV['RTTOV_GUI_PREFIX'])
        self.write("filename= " + surfaceFileName)
        if (surfaceFileName):
            try:
                self.project.openSurface(surfaceFileName)
            except IOError:
                self.write("Error while opening " + surfaceFileName)

    def createEmissivityReflectance(self, e):
        """ create EmissivityReflectance values with defaut values """
        self.prepLog()
        logging.debug("surface controller create EmissivityReflectance ")
        self.project.createSurface()

    def loadAtlasEmissivityReflectance(self, e):
        """ Load EmissivityReflectance values from Atlas"""
        self.prepLog()
        logging.debug("surface controller loadAtlasEmissivityReflectance ")
        if self.project.myProfile['SKIN']['SURFTYPE'] == 1:
            msg = "No atlas over sea !"
            self.surfaceView.WarmError(msg)
        else:
            if self.project.myCoeffs.hasSolar():
                msg = "No BRDF data for sea ice surface"
                if self.project.myProfile['SKIN']['SURFTYPE'] == 2:
                    self.surfaceView.WarmError(msg)
                else:
                    self.project.openAtlas()
            else:
                self.project.openAtlas()

    def modifyEmissivity(self, e):
        """ modify Emissivity values """
        self.emisController.OnGridView(e)

    def modifyReflectance(self, e):
        """ modify Emissivity values """
        self.reflController.OnGridView(e)

    def saveProfile(self, e):
        """ apply profile change and save profile+option """
        self.prepLog()
        logging.debug("debug saveProfile surface controller")
        # self.project.myProfile is maybe modified elsewhere (profile editor)
        # (but not necessary with the listener)
        profile = self.surfaceView.OnApplyChange(e, self.project.myProfile)
        self.project.updateProfile(profile)
        self.SaveProfile(e, self.surfaceView)

    def saveProfileAs(self, e):
        """ apply profile change and save profile+option with a new name """
        logging.debug("debug saveProfileAs surface controller")
        # re-initialize the name of the saved profile file to None
        self.project.savedProfileFileName = None
        self.saveProfile(e)

    def saveSurface(self, e):
        """ apply surface file default name """

        logging.debug("debug saveSurface surface controller")
        # self.project.myProfile is maybe modified elsewhere (profile editor)
        if self.project.savedSurfaceFileName is not None:
            self.project.saveSurface(self.project.savedSurfaceFileName)
        else:
            self.saveSurfaceAs(e)

    def saveSurfaceAs(self, e):
        """ save  surface file with a new name """
        self.prepLog()
        logging.debug("debug saveSurfaceAs surface controller")
        # re-initialize the name of the saved profile file to None
        fileName = self.surfaceView.OnSaveFileAs(e, title="Save Surface File")
        if fileName is not None:
            self.project.saveSurface(fileName)
            self.project.savedSurfaceFileName = fileName

    def MakeSurfaceBinding(self):
        # binding for the SurfaceView
        self.surfaceView.Bind(
            wx.EVT_MENU, self.applyProfileChange,
            self.surfaceView.items['applyChange'])
        self.surfaceView.Bind(
                wx.EVT_MENU, self.createEmissivityReflectance,
                self.surfaceView.items['createEmissivityReflectance'])
        self.surfaceView.Bind(
                wx.EVT_MENU, self.modifyEmissivity,
                self.surfaceView.items['modifyEmissivity'])
        self.surfaceView.Bind(
                wx.EVT_MENU, self.modifyReflectance,
                self.surfaceView.items['modifyReflectance'])
        self.surfaceView.Bind(wx.EVT_MENU, self.saveProfileAs,
                              self.surfaceView.items['saveProfileAs'])
        self.surfaceView.Bind(
            wx.EVT_MENU,
            self.loadAtlasEmissivityReflectance,
            self.surfaceView.items['loadAtlas'])
        self.surfaceView.Bind(wx.EVT_MENU, self.saveProfile,
                              self.surfaceView.items['saveProfile'])
        self.surfaceView.Bind(wx.EVT_MENU, self.saveSurface,
                              self.surfaceView.items['saveSurface'])
        self.surfaceView.Bind(wx.EVT_MENU, self.saveSurfaceAs,
                              self.surfaceView.items['saveSurfaceAs'])
        self.surfaceView.Bind(
            wx.EVT_BUTTON, self.applyProfileChange, self.surfaceView.applyBtn)
        self.surfaceView.Bind(wx.EVT_CLOSE, self.OnClose)

    def OnClose(self, Event):
        if self.surfaceView is not None:
            self.surfaceView.Destroy()
        self.surfaceView = None


class ReflEmisController(GenericController):
    """ controller for the Emissivity or Reflectance editor windows """

    def __init__(self, parent, dataName):
        """ Constructor of the Emissivity or Reflectance window
            create the GridView Object
            make the binding with the GridView

        """

        self.dataName = dataName
        self.theParentFrame = parent
        super(GenericController, self).__init__()
        self.controlerName = "ReflEmisController"
        self.gridView = None
        logging.debug("Create ReflEmisController")
        pub.subscribe(self.SurfaceChangedListener, "Surface CHANGED")

    def Close(self):
        if self.gridView is not None:
            self.gridView.Close()
            self.gridView = None

    def SurfaceChangedListener(self, msg):
        logging.debug("ReflEmisController" + self.dataName +
                      "SurfaceChangedListener")
        if self.gridView is not None:
            if self.dataName == 'REFL':
                logging.debug(str(self.project.myReflectance))
                logging.debug(
                    "call self.gridView.OnUpdateData"
                    "(self.project.myReflectance REFL")
                self.gridView.OnUpdateData(self.project.myReflectance, "REFL")
            else:
                logging.debug(
                    "call self.gridView.OnUpdateData("
                    "self.project.myReflectance EMIS")
                self.gridView.OnUpdateData(self.project.myEmissivity, "EMIS")
        else:
            logging.debug("no gridView")

    def OnGridView(self, e):
        if self.gridView is None:
            if self.dataName == "REFL":
                logging.debug("create GridView REFL")
                self.gridView = rview.surfedit.GridView(
                    self.theParentFrame, self.project.myReflectance,
                    "Surface BRDF", "REFL", controller=self)
                self.MakeBinding()
            else:
                logging.debug("create GridView EMIS")
                self.gridView = rview.surfedit.GridView(
                    self.theParentFrame, self.project.myEmissivity,
                    "Emissivity", "EMIS", controller=self)
                self.MakeBinding()
        else:
            pass

    def OnApplyChange(self, e):
        logging.debug("controller Emis/refl " +
                      self.dataName + " OnApplyChange")
        self.prepLog()
        data = self.gridView.OnApplyChange(e)
        if self.dataName == 'REFL':
            self.project.updateReflectance(data)
        else:
            self.project.updateEmissivity(data)

    def MakeBinding(self):
        # binding for the SurfaceView
        self.gridView.Bind(wx.EVT_MENU, self.OnApplyChange,
                           self.gridView.items['applyChange'])
