# -*- coding: utf-8 -*-

import wx
import rmodel
import h5py
import rttov
import util
from rview import layeritem
import locale
import numpy
import matplotlib
import sys
import copy
import logging
from profileframeutils import GenericPlotItemPanel, MyNotebook
from profileframeutils import kindOfItem, PlotItemPanelAll


locale.setlocale(locale.LC_ALL, '')


class PlotItemPanel(GenericPlotItemPanel):
    """ plot on a PlotPanel one curve """

    def __init__(self, parent, value, pression, theName, liste_item=None,
                 kind="GASES", xlegend="(mw/cm-1/ster/sq.m)/ppmv",
                 layerstyle=False, layer=None, yInPressions=True, tskin=None):
        edit = False
        GenericPlotItemPanel.__init__(self, parent, value, pression, theName,
                                      liste_item, kind, xlegend,
                                      edit, layerstyle, layer, yInPressions,
                                      tskin)
        self.SetTickSize(8)


class KProfileView (util.GenericViewRadio):
    """ Profile window of the application """
    helpTitle = "Help Profile"
    helpMessage = """
    Select and visualize a component profile on the right panel
    Click left button to modify the profile.
    Click left and drag a zone to zoom in.
    Click right button to zoom out.
    Apply your changes or save the profile
    for the next run of RTTOV.
        """

    def __init__(self, parent, profile, channel=1, baseProfile=None,
                 edit=False, yInPressions=True, runNumber=1):
        self.edit = edit
        logging.debug("profile" + str(profile['T'].shape))
        logging.debug("baseProfile" + str(baseProfile['T'].shape))
        self.myProfile = copy.deepcopy(profile)
        self.yInPressions = yInPressions
        if baseProfile is not None:
            self.myProfile['P'] = baseProfile['P']
            self._ComputeLayers(self.myProfile['P'])
        else:
            self.yInPressions = False

        self.my_list_cloud = []
        for layeritem in self.myProfile.cloud_list:
            self.my_list_cloud.append(layeritem)
        self.my_list_cloud.append("CFRAC")
        self.my_list_cloud.append("CLW")
        util.GenericView.__init__(self, parent,  "PROFILE")

        sizer = wx.BoxSizer(wx.HORIZONTAL)
        self.SetSizer(sizer)
        self.CreateMenuBar()
        if baseProfile is not None:
            self.items["ypressions"].Enable(True)
            self.items["ypressions"].Check(True)
            self.items["ylevels"].Check(False)
        else:

            self.items["ypressions"].Enable(False)
            self.items["ypressions"].Check(False)
            self.items["ylevels"].Check(True)
        self.SetSize((1450, 700))
        self.SetMinSize((1300, 700))
        self.SetTitle("run %d k profile channel %d" % (runNumber, channel))
        # panel 1 notebook with all curves (GASES, AEROSOLS, CLOUDS)
        self.panel1 = wx.Panel(self, -1, style=wx.BORDER_SIMPLE)
        self.panel1.SetSize((200, 500))
        sizer.Add(self.panel1, 1, wx.EXPAND)
        # creation of notebook for the panel 1
        self.nb_all = MyNotebook(self.panel1, isRightPage=False)

        sizer1 = wx.BoxSizer()
        sizer1.Add(self.nb_all, 1, wx.EXPAND)
        self.panel1.SetSizer(sizer1)

        # panel 2 notebook with one curve
        self.panel2 = wx.Panel(self, -1, style=wx.BORDER_SIMPLE)
        sizer.Add(self.panel2, 1, wx.EXPAND)
        # creation of the notebook for the panel 2
        self.nb = MyNotebook(self.panel2)
        self.axesDef = []

        # creation des graphiques
        self.Plot(self.myProfile)

        # create a second sizer for the notebook
        sizer2 = wx.BoxSizer()
        sizer2.Add(self.nb, 1, wx.EXPAND)
        self.panel2.SetSizer(sizer2)
        self.sb = self.CreateStatusBar()
        self.sb.SetBackgroundColour('WHITE')
        txt = ''
        self.sb.SetStatusText(txt)
        self.Centre()
        self.Show(True)

    def PlotLeft(self, profile=None):
        # plot panel 1 with all gas
        self.allGraphicsPages = {}
        self.allGraphicsPages['GASES'] = PlotItemPanelAll(
                                            self.nb_all,
                                            self.myProfile,
                                            kind='GASES',
                                            xlegendT=self.myProfile[
                                                      'T_ATTRIBUTE']['UNITS'],
                                            yInPressions=self.yInPressions,
                                            addTskin=True)
        self.allGraphicsPages['GASES'].SetTickSize(8)
        self.nb_all.AddPage(self.allGraphicsPages['GASES'], 'GASES')

        if self.myProfile.anyAerosol():
            self.allGraphicsPages['AEROSOLS'] = PlotItemPanelAll(
                                                    self.nb_all,
                                                    self.myProfile,
                                                    kind="AEROSOLS",
                                                    layer=self.layer,
                                                    xlegendT=self.myProfile[
                                                       'T_ATTRIBUTE']['UNITS'],
                                                    yInPressions=(
                                                        self.yInPressions))
            self.allGraphicsPages['AEROSOLS'].SetTickSize(8)
            self.nb_all.AddPage(self.allGraphicsPages['AEROSOLS'], 'AEROSOLS')

        if self.myProfile.anyCloud():
            self.allGraphicsPages['CLOUDS'] = PlotItemPanelAll(
                                                self.nb_all,
                                                self.myProfile,
                                                xlegendT=self.myProfile[
                                                       'T_ATTRIBUTE']['UNITS'],
                                                kind="CLOUDS",
                                                layer=self.layer,
                                                yInPressions=self.yInPressions)

            self.allGraphicsPages['CLOUDS'].SetTickSize(8)
            self.nb_all.AddPage(self.allGraphicsPages['CLOUDS'], 'CLOUDS')

    def Plot(self, profile=None):
        if profile is not None:
            self.myProfile = profile
            self._ComputeLayers(self.myProfile['P'])
        self.graphicPages = {}
        self.graphicPages['T'] = PlotItemPanel(self.nb,
                                               self.myProfile['T'],
                                               self.myProfile['P'],
                                               theName='T',
                                               xlegend=self.myProfile[
                                                   'T_ATTRIBUTE']['UNITS'],
                                               yInPressions=self.yInPressions,
                                               tskin=self.myProfile[
                                                                'SKIN']['T'])
        self.nb.AddPage(self.graphicPages['T'], 'T')

        for gas in self.myProfile.gas_list:
            if self.myProfile[gas] is not None:
                self.graphicPages[gas] = PlotItemPanel(self.nb,
                                                       self.myProfile[gas],
                                                       self.myProfile['P'],
                                                       theName=gas,
                                                       xlegend=self.myProfile[
                                                           gas + '_ATTRIBUTE'][
                                                                    'UNITS'],
                                                       yInPressions=(
                                                           self.yInPressions))
                self.nb.AddPage(self.graphicPages[gas], gas)

        if self.myProfile.anyAerosol():
            for aerosol in self.myProfile.aerosol_list:
                if self.myProfile[aerosol] is not None:
                    self.graphicPages[aerosol] = PlotItemPanel(
                                                    self.nb,
                                                    self.myProfile[aerosol],
                                                    self.myProfile['P'],
                                                    theName=aerosol,
                                                    kind="AEROSOLS",
                                                    xlegend=self.myProfile[
                                                       aerosol + '_ATTRIBUTE'][
                                                                    'UNITS'],
                                                    layer=self.layer,
                                                    yInPressions=(
                                                       self.yInPressions),
                                                    layerstyle=True)
                    self.nb.AddPage(self.graphicPages[aerosol], aerosol)

        if self.myProfile.anyCloud():
            for cloud in self.my_list_cloud:
                if self.myProfile[cloud] is not None:
                    if cloud == "CFRAC":
                        mylegend = "CFRAC"
                    else:
                        mylegend = "layer mean content (g/m3)"
                    self.graphicPages[cloud] = PlotItemPanel(
                                                    self.nb,
                                                    self.myProfile[cloud],
                                                    self.myProfile['P'],
                                                    theName=cloud,
                                                    kind="CLOUDS",
                                                    xlegend=mylegend,
                                                    layer=self.layer,
                                                    yInPressions=(
                                                            self.yInPressions),
                                                    layerstyle=True)
                    self.nb.AddPage(self.graphicPages[cloud], cloud)

        # delete empty graphicPages
        for key in self.graphicPages.keys():
            if self.myProfile[key] is None:
                del self.graphicPages[key]

        # plot panel 1 with all gas
        self.PlotLeft()

    def _ComputeLayers(self, pression):
        """ Compute the mean value of pression in a layer """
        foo = numpy.empty(pression.shape[0] - 1)
        for i in range(foo.shape[0]):
            foo[i] = (pression[i + 1] + pression[i]) / 2
        self.layer = foo

    def _MakeBinding(self):
        """ set the trivial Binding for the View """
        # binding cancel button

    def RePlotAll(self, profile=None):
        """ Plot the 2 panels with (new) profile
            (delete everything before redraw) """

        if profile is not None:
            self.myProfile = profile
            self._ComputeLayers(self.myProfile['P'])
        # remove all pages of the notebook
        self.nb.DeleteAllPages()
        self.nb_all.DeleteAllPages()
        self.Plot()

    def RePlotAllLeftPanel(self, profile=None):
        """ Plot the 2 panels with (new) profile
            (delete everything before redraw) """

        if profile is not None:
            self.myProfile = profile
            self._ComputeLayers(self.myProfile['P'])
        # remove all pages of the notebook
        self.nb_all.DeleteAllPages()
        self.PlotLeft()

    def addRightPage(self, item):
        """ add an new item page  """
        kind = kindOfItem[item]
        if kind == "GASES":
            myY = self.myProfile['P']
        else:
            myY = self.layer
        self.graphicPages[item] = PlotItemPanel(self.nb, self.myProfile[item],
                                                self.myProfile['P'],
                                                theName=item, kind=kind,
                                                xlegend=self.myProfile[
                                                        item + '_ATTRIBUTE'][
                                                                    'UNITS'],
                                                layer=myY,
                                                yInPressions=self.yInPression)
        self.nb.AddPage(self.graphicPages[item], item)
        if kind == "CLOUDS":
            if 'CFRAC' not in self.graphicPages:
                item = 'CFRAC'
                self.graphicPages[item] = PlotItemPanel(
                                            self.nb, self.myProfile[item],
                                            self.myProfile['P'], theName=item,
                                            kind=kind, xlegend="CFRAC",
                                            layer=myY,
                                            yInPressions=self.yInPression)
                self.nb.AddPage(self.graphicPages[item], item)
            if 'CLW' not in self.graphicPages and self.myProfile[
                                                            'CLW'] is not None:
                item = "CLW"
                self.graphicPages[item] = PlotItemPanel(
                                            self.nb, self.myProfile[item],
                                            self.myProfile['P'], theName=item,
                                            kind=kind, xlegend="CLW",
                                            layer=myY,
                                            yInPressions=self.yInPression)
                self.nb.AddPage(self.graphicPages[item], item)

    def OnClose(self, e):
        """ close the surface windows"""
        self.Close()

    def OnMouseMove(self, e):
        """ print x y value of the left plot in the status bar  """
        pass

    def OnYPressions(self, e):
        self.yInPressions = True
        self.RePlotAll()

    def OnYLevels(self, e):
        self.yInPressions = False
        self.RePlotAll()

    def MenuData(self):
        """ define the data for the menu
        """
        return(("&File",  # File Menu
                ('&Quit', 'Quit', self.OnQuit, "quit", True, False)),
               ("&Edit",  # Edit Menu
                ("Yaxis in pressure units", "put y in pressure units",
                 self.OnYPressions, "ypressions", False, True),
                ("Yaxis in level units", "put y in level unit",
                 self.OnYLevels, "ylevels", True, True)),

               ("&Help",  # Help Menu
                ("About", "About screen", self.OnAbout, "about", True, False)))


if __name__ == "__main__":

    print "version matplotlib :", matplotlib.__version__
    p = rmodel.project.Project()
    ex = wx.App()
    fh5 = sys.argv[1]
    print "f=", fh5
    kmat = rttov.kmatrix.Kmatrix()
    frad = h5py.File(fh5, 'r')
    h5 = frad['/']
    kmat.loadh5(h5)
    baseProfile = rttov.profile.Profile()

    frad.close()
    baseProfile = rmodel.project.OpenAProfile(fh5, 1)

    profile = kmat.getchanprof(1)
    print "tskin", profile['SKIN']['T']

    profile.display()
    print "T shape", profile['T'].shape
    print " NLEVELS : ", profile['NLEVELS']
    print profile['NLEVELS']
    print "P"
    print baseProfile["P"]
    print "T"
    print profile['T']
    print max(profile['T'])

    frame = KProfileView(None, profile, channel=1,
                         yInPressions=True, baseProfile=baseProfile)
    frame.Show()

    ex.MainLoop()

    print "loop"
    ex.MainLoop()
    # ex.MainLoop()
