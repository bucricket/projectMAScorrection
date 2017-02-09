# -*- coding: utf-8 -*-

import wx
import rmodel
import util
import locale
import wxmpl
import numpy
import matplotlib
import wx.lib.agw.floatspin as FS
from profileframeutils import Data, GenericPlotItemPanel, MyNotebook
from profileframeutils import PlotItemPanelAll
from util import kindOfItem, axesDef
import copy


locale.setlocale(locale.LC_ALL, '')


class ComputeDialog(wx.Dialog):

    def __init__(
            self, parent, ID, title, text="", pmin0=0, pmax0=1000, scale0=1.,
            offset0=0., increment=None, size=wx.DefaultSize,
            pos=wx.DefaultPosition,
            style=wx.DEFAULT_DIALOG_STYLE,
            useMetal=False, limits_scale=(0.01, 10), limits_offset=(-100, 100)
    ):

        # Instead of calling wx.Dialog.__init__ we precreate the dialog
        # so we can set an extra style that must be set before
        # creation, and then we create the GUI object using the Create
        # method.
        if increment is None:
            increment = 10.
        pre = wx.PreDialog()
        pre.SetExtraStyle(wx.DIALOG_EX_CONTEXTHELP)
        pre.Create(parent, ID, title, pos, size, style)

        # This next step is the most important, it turns this Python
        # object into the real wrapper of the dialog (instead of pre)
        # as far as the wxPython extension is concerned.
        self.PostCreate(pre)

        # This extra style can be set after the UI object has been created.
        if 'wxMac' in wx.PlatformInfo and useMetal:
            self.SetExtraStyle(wx.DIALOG_EX_METAL)

        # Now continue with the normal construction of the dialog
        # contents
        sizer = wx.BoxSizer(wx.VERTICAL)

        label = wx.StaticText(self, -1, text)
        label.SetHelpText("help text")
        sizer.Add(label, 0, wx.ALIGN_CENTRE | wx.ALL, 5)

        # box 1 pmin
        box = wx.BoxSizer(wx.HORIZONTAL)

        label = wx.StaticText(self, -1, "p min:")
        box.Add(label, 0, wx.ALIGN_CENTRE | wx.ALL, 5)
        self.pminFS = FS.FloatSpin(self, -1, min_val=0.00001, max_val=1150,
                                   digits=3, value=pmin0,
                                   increment=increment, agwStyle=FS.FS_LEFT)
        box.Add(self.pminFS, 1, wx.ALIGN_CENTRE | wx.ALL, 5)
        sizer.Add(box, 0, wx.GROW | wx.ALIGN_CENTER_VERTICAL | wx.ALL, 5)

        # box 2 pmax
        box = wx.BoxSizer(wx.HORIZONTAL)
        label = wx.StaticText(self, -1, "p max:")
        box.Add(label, 0, wx.ALIGN_CENTRE | wx.ALL, 5)
        self.pmaxFS = FS.FloatSpin(self, -1, min_val=0.00001, max_val=1150,
                                   digits=3, value=pmax0,
                                   increment=increment, agwStyle=FS.FS_LEFT)
        box.Add(self.pmaxFS, 1, wx.ALIGN_CENTRE | wx.ALL, 5)
        sizer.Add(box, 0, wx.GROW | wx.ALIGN_CENTER_VERTICAL | wx.ALL, 5)

        # box 3 scale
        box = wx.BoxSizer(wx.HORIZONTAL)
        label = wx.StaticText(self, -1, "scale:")
        box.Add(label, 0, wx.ALIGN_CENTRE | wx.ALL, 5)
        self.scale = FS.FloatSpin(self, -1, min_val=limits_scale[0],
                                  max_val=limits_scale[1], digits=2,
                                  value=scale0,
                                  increment=1, agwStyle=FS.FS_LEFT)
        box.Add(self.scale, 1, wx.ALIGN_CENTRE | wx.ALL, 5)

        # box 4 offset
        label = wx.StaticText(self, -1, "offset:")
        box.Add(label, 0, wx.ALIGN_CENTRE | wx.ALL, 5)
        self.offset = FS.FloatSpin(self, -1, min_val=limits_offset[0],
                                   max_val=limits_offset[1], digits=3,
                                   value=offset0,
                                   increment=1, agwStyle=FS.FS_LEFT)
        box.Add(self.offset, 1, wx.ALIGN_CENTRE | wx.ALL, 5)
        sizer.Add(box, 0, wx.GROW | wx.ALIGN_CENTER_VERTICAL | wx.ALL, 5)

        line = wx.StaticLine(self, -1, size=(20, -1), style=wx.LI_HORIZONTAL)
        sizer.Add(line, 0, wx.GROW | wx.ALIGN_CENTER_VERTICAL |
                  wx.RIGHT | wx.TOP, 5)

        # buttons
        btnsizer = wx.StdDialogButtonSizer()

        btn = wx.Button(self, wx.ID_OK)

        btn.SetDefault()
        btnsizer.AddButton(btn)
        btn = wx.Button(self, wx.ID_CANCEL)
        btnsizer.AddButton(btn)
        btnsizer.Realize()

        sizer.Add(btnsizer, 0, wx.ALIGN_CENTER_VERTICAL | wx.ALL, 5)

        self.SetSizer(sizer)
        sizer.Fit(self)
        self.Show()

    def GetSelections(self):

        return (self.pminFS.GetValue(), self.pmaxFS.GetValue(),
                self.scale.GetValue(), self.offset.GetValue())


class AxesDialog(wx.Dialog):

    def __init__(
            self, parent, ID, title, text="", xmin0=0, xmax0=0,
            scale="linear", increment=None, size=wx.DefaultSize,
            pos=wx.DefaultPosition,
            style=wx.DEFAULT_DIALOG_STYLE,
            useMetal=False,
    ):

        # Instead of calling wx.Dialog.__init__ we precreate the dialog
        # so we can set an extra style that must be set before
        # creation, and then we create the GUI object using the Create
        # method.
        if increment is None:
            increment = 10
        pre = wx.PreDialog()
        pre.SetExtraStyle(wx.DIALOG_EX_CONTEXTHELP)
        pre.Create(parent, ID, title, pos, size, style)

        # This next step is the most important, it turns this Python
        # object into the real wrapper of the dialog (instead of pre)
        # as far as the wxPython extension is concerned.
        self.PostCreate(pre)

        # This extra style can be set after the UI object has been created.
        if 'wxMac' in wx.PlatformInfo and useMetal:
            self.SetExtraStyle(wx.DIALOG_EX_METAL)

        # Now continue with the normal construction of the dialog
        # contents
        sizer = wx.BoxSizer(wx.VERTICAL)

        label = wx.StaticText(self, -1, text)
        label.SetHelpText("help text")
        sizer.Add(label, 0, wx.ALIGN_CENTRE | wx.ALL, 5)

        # box 1 xmin
        box = wx.BoxSizer(wx.HORIZONTAL)

        label = wx.StaticText(self, -1, "x min:")
        box.Add(label, 0, wx.ALIGN_CENTRE | wx.ALL, 5)
        self.xminFS = FS.FloatSpin(self, -1, min_val=0, max_val=1000000000,
                                   digits=2, value=xmin0,
                                   increment=increment, agwStyle=FS.FS_LEFT)
        box.Add(self.xminFS, 1, wx.ALIGN_CENTRE | wx.ALL, 5)
        sizer.Add(box, 0, wx.GROW | wx.ALIGN_CENTER_VERTICAL | wx.ALL, 5)

        # box 2 xmax
        box = wx.BoxSizer(wx.HORIZONTAL)
        label = wx.StaticText(self, -1, "x max:")
        box.Add(label, 0, wx.ALIGN_CENTRE | wx.ALL, 5)
        self.xmaxFS = FS.FloatSpin(self, -1, min_val=0, max_val=1000000000,
                                   digits=2, value=xmax0,
                                   increment=increment, agwStyle=FS.FS_LEFT)
        box.Add(self.xmaxFS, 1, wx.ALIGN_CENTRE | wx.ALL, 5)
        sizer.Add(box, 0, wx.GROW | wx.ALIGN_CENTER_VERTICAL | wx.ALL, 5)

        # buttons
        btnsizer = wx.StdDialogButtonSizer()

        btn = wx.Button(self, wx.ID_OK)

        btn.SetDefault()
        btnsizer.AddButton(btn)
        btn = wx.Button(self, wx.ID_CANCEL)
        btnsizer.AddButton(btn)
        btnsizer.Realize()

        sizer.Add(btnsizer, 0, wx.ALIGN_CENTER_VERTICAL | wx.ALL, 5)

        self.SetSizer(sizer)
        sizer.Fit(self)
        self.Show()

    def GetSelections(self):
        return (self.xminFS.GetValue(), self.xmaxFS.GetValue())


class PlotItemPanel(GenericPlotItemPanel):
    """ plot on a PlotPanel one curve """

    def __init__(self, parent, value, pression, theName, liste_item=None,
                 kind="GASES", xlegend="ppmv",
                 layerstyle=False, layer=None, yInPressions=True, tskin=None):
        edit = True
        GenericPlotItemPanel.__init__(self, parent, value, pression, theName,
                                      liste_item,
                                      kind=kind, xlegend=xlegend,
                                      edit=edit, layerstyle=layerstyle,
                                      layer=layer, yInPressions=yInPressions,
                                      tskin=tskin)

    def ComputeProfile(self, pmin, pmax, scale, offset):
        """ change the values of the profile between pmin and pmax :
                multiply by scale and add offset  """

        if not self.layerstyle:
            y = numpy.zeros(self.x.shape[0]) + self.x
            self.valueHistory.append(y)
            for i in range(self.x.shape[0]):
                p = self.pression[i]
                if p >= pmin and p <= pmax:
                    val = self.x[i]
                    self.x[i] = max(0.0001, scale * val + offset)
            self.data.setChanged(True)
            self.data.myUpdate(self.x, self.y)
            self.Update()

        else:
            y = numpy.zeros(self.xlayeritem.shape[0]) + self.xlayeritem
            self.valueHistory.append(y)
            for i in range(self.xlayeritem.shape[0]):
                p = self.ylayeritem[i]
                if p >= pmin and p <= pmax:
                    self.xlayeritem[i] = max(
                        0, scale * self.xlayeritem[i] + offset)
            self.data.setChanged(True)
            self.data.myUpdate(self.xlayeritem, self.ylayeritem)
            self.x = self.myLayeritem.getLayeritem(self.xlayeritem)
            self.Update()


class ProfileView(util.GenericViewRadio):
    """ Profile window of the application """
    helpTitle = "Help Profile"
    helpMessage = """
    Select and visualize a component profile on the right panel
    Click the middle button or the right button to modify the profile.
    Use the matplotlib toolbar to zoom.
    Apply your changes or save the profile
    for the next run of RTTOV.
        """

    def __init__(self, parent, profile, edit=True, controler=None):
        self.edit = edit
        self.myProfileRef = profile
        self.myProfile = copy.deepcopy(profile)
        self.my_list_cloud = []
        self.myControler = controler
        for cloud in self.myProfile.cloud_list:
            self.my_list_cloud.append(cloud)
        self.my_list_cloud.append("CFRAC")
        self.my_list_cloud.append("CLW")
        util.GenericView.__init__(self, parent,  "profile")

        sizer = wx.BoxSizer(wx.HORIZONTAL)
        self.SetSizer(sizer)
        self.CreateMenuBar()
        self.SetSize((1300, 700))
        self.SetMinSize((1300, 700))
        self.SetTitle('PROFILE EDITOR')
        self._ComputeLayers(self.myProfile['P'])
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

    def PlotLeft(self):
        # plot panel 1 with all gas
        self.allGraphicsPages = {}

        self.allGraphicsPages['GASES'] = PlotItemPanelAll(
                                    self.nb_all,
                                    self.myProfile,
                                    kind='GASES',
                                    xlegend=self.myProfile.gas_units_legend,
                                    xlegendT="Temperature (K)",
                                    tickSize=10)
        self.nb_all.AddPage(self.allGraphicsPages['GASES'], 'GASES')

        if self.myProfile.anyAerosol():
            self.allGraphicsPages['AEROSOLS'] = PlotItemPanelAll(
                                            self.nb_all,
                                            self.myProfile,
                                            kind="AEROSOLS",
                                            layer=self.layer,
                                            xlegend="number density (cm-3)",
                                            xlegendT="Temperature (K)",
                                            XinLog=True,
                                            tickSize=10)
            self.nb_all.AddPage(self.allGraphicsPages['AEROSOLS'], 'AEROSOLS')

        if self.myProfile.anyCloud():
            self.allGraphicsPages['CLOUDS'] = PlotItemPanelAll(
                                        self.nb_all,
                                        self.myProfile,
                                        kind="CLOUDS",
                                        layer=self.layer,
                                        xlegend="layer mean content (g/m3)",
                                        xlegendT="Temperature (K)",
                                        XinLog=True,
                                        tickSize=10)
            self.nb_all.AddPage(self.allGraphicsPages['CLOUDS'], 'CLOUDS')

    def Plot(self, profile=None):
        if profile is not None:
            self.myProfileRef = profile
            self.myProfile = copy.deepcopy(profile)
            self._ComputeLayers(self.myProfile['P'])
        self.graphicPages = {}
        self.graphicPages['T'] = PlotItemPanel(self.nb, self.myProfile['T'],
                                               self.myProfile['P'],
                                               theName='T',
                                               xlegend="temperature (K)")
        self.nb.AddPage(self.graphicPages['T'], 'T')
        self.graphicPages['T'].ConnectCanvasEVT_POINT(self.OnPoint)

        for gas in self.myProfile.gas_list:
            if self.myProfile[gas] is not None:
                self.graphicPages[gas] = PlotItemPanel(
                                    self.nb,
                                    self.myProfile[gas],
                                    self.myProfile['P'],
                                    xlegend=self.myProfile.gas_units_legend,
                                    theName=gas)
                self.nb.AddPage(self.graphicPages[gas], gas)
                self.graphicPages[gas].ConnectCanvasEVT_POINT(self.OnPoint)

        if self.myProfile.anyAerosol():
            for aerosol in self.myProfile.aerosol_list:
                if self.myProfile[aerosol] is not None:
                    self.graphicPages[aerosol] = PlotItemPanel(
                                            self.nb,
                                            self.myProfile[aerosol],
                                            self.myProfile['P'],
                                            theName=aerosol,
                                            kind="AEROSOLS",
                                            layerstyle=True,
                                            xlegend="number density (cm-3)",
                                            layer=self.layer)
                    self.nb.AddPage(self.graphicPages[aerosol], aerosol)
                    self.graphicPages[
                        aerosol].ConnectCanvasEVT_POINT(self.OnPoint)

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
                                                layerstyle=True,
                                                xlegend=mylegend,
                                                layer=self.layer)
                    self.nb.AddPage(self.graphicPages[cloud], cloud)
                    self.graphicPages[
                        cloud].ConnectCanvasEVT_POINT(self.OnPoint)

        # delete empty graphicPages
        for key in self.graphicPages.keys():
            if self.myProfile[key] is None:
                del self.graphicPages[key]

        # plot panel 1 with all gas
        self.PlotLeft()

    def OnPoint(self, e):
        """ update the graphic when receive event e """
        wx.BeginBusyCursor()
        wx.SafeYield(None, True)
        cp = self.nb.GetCurrentPage()
        hasDoSomething = cp.OnPoint(e)
        if hasDoSomething:
            self.UpdateAllGraphics()
        wx.EndBusyCursor()

    def UpdateAllGraphics(self):
        """ update the left panel after getting all information
            from the right panels"""
        # get the name of the RightPlotPanel
        cp = self.nb.GetCurrentPage()
        theName = cp.theName
        theKind = cp.kind
        self.GetProfile()
        # update the profile with curves from all rights panels

        self.allGraphicsPages[theKind].UpdateData(self.myProfile, theName)
        if theName == "T":
            for kind in ("AEROSOLS", "CLOUDS"):
                if kind in self.allGraphicsPages:
                    self.allGraphicsPages[kind].UpdateData(
                        self.myProfile, theName)

    def _ComputeLayers(self, pression):
        """ Compute the mean value of pression in a layer """
        foo = numpy.empty(pression.shape[0] - 1)
        for i in range(foo.shape[0]):
            foo[i] = (pression[i + 1] + pression[i]) / 2
        self.layer = foo

    def _MakeBinding(self):
        """ set the trivial Binding for the View """
        # binding cancel button

    def GetProfile(self):
        """ return the profile as modified by the edition frame """
        self.write("Get profile from the edition panels")
        self.myProfile["T"] = self.graphicPages["T"].GetItem()
        for gas in self.myProfile.gas_list:
            if self.myProfile[gas]is not None:
                self.myProfile[gas] = self.graphicPages[gas].GetItem()
        if self.myProfile.anyAerosol:
            for aerosol in self.myProfile.aerosol_list:
                if self.myProfile[aerosol] is not None:
                    self.myProfile[aerosol] = self.graphicPages[
                        aerosol].GetItem()

        if self.myProfile.anyCloud:
            for item in self.my_list_cloud:
                if self.myProfile[item] is not None:
                    self.myProfile[item] = self.graphicPages[item].GetItem()

        return self.myProfile

    def RePlotAll(self, profile=None):
        """ Plot the 2 panels with (new) profile (delete everything
            before redraw) """

        if profile is not None:
            self.myProfileRef = profile
            self.myProfile = copy.deepcopy(profile)
            self._ComputeLayers(self.myProfile['P'])
        # remove all pages of the notebook
        self.nb.DeleteAllPages()
        self.nb_all.DeleteAllPages()
        self.Plot()

    def RePlotAllLeftPanel(self, profile=None):
        """ Plot the left panel with (new) profile (delete everything
            before redraw) """

        if profile is not None:
            self.myProfileRef = profile
            self.myProfile = copy.deepcopy(profile)
            self._ComputeLayers(self.myProfile['P'])
        # remove all pages of the notebook
        self.nb_all.DeleteAllPages()
        self.PlotLeft()

    def deleteRightPage(self, item):
        """ delete the item page """
        # delete the notebook page
        # delete the graphic page
        for i in range(self.nb.GetPageCount()):
            if self.nb.GetPageText(i) == item:
                self.graphicPages[item].DisconnectCanvasEVT_POINT()
                self.nb.DeletePage(i)
                self.graphicPages.pop(item)
                break
        kind = kindOfItem[item]
        if kind == "CLOUDS" and not self.myProfile.anyCloud():
            if 'CFRAC' in self.graphicPages:
                self.deleteRightPage('CFRAC')
            if 'CLW' in self.graphicPages:
                self.deleteRightPage('CLW')

    def addRightPage(self, item):
        """ add an new item page  """
        kind = kindOfItem[item]
        if kind == "GASES":
            myY = self.myProfile['P']
            myLayerstyle = False
        else:
            myY = self.layer
            myLayerstyle = True
        self.graphicPages[item] = PlotItemPanel(
                        self.nb, self.myProfile[item], self.myProfile['P'],
                        theName=item, kind=kind,
                        layerstyle=myLayerstyle,
                        xlegend=self.myProfile.gas_units_legend, layer=myY)
        self.nb.AddPage(self.graphicPages[item], item)
        if kind == "CLOUDS":
            if 'CFRAC' not in self.graphicPages:
                item = 'CFRAC'
                self.graphicPages[item] = PlotItemPanel(
                        self.nb, self.myProfile[item], self.myProfile['P'],
                        theName=item, kind=kind,
                        layerstyle=True,
                        xlegend="CFRAC", layer=myY)
                self.nb.AddPage(self.graphicPages[item], item)
        self.graphicPages[item].ConnectCanvasEVT_POINT(self.OnPoint)

    def replotLeftPage(self):
        """ replot the left panel """
        self.RePlotAllLeftPanel()

    def OnUndo(self, e):
        pageName = self.nb.GetPage(self.nb.GetSelection()).theName
        self.graphicPages[pageName].OnUndo()
        self.UpdateAllGraphics()

    def OnRedo(self, e):
        pageName = self.nb.GetPage(self.nb.GetSelection()).theName
        self.graphicPages[pageName].OnRedo()
        self.UpdateAllGraphics()

    def OnClose(self, e):
        """ close the profile windows"""
        print ">>> Close profileView"
        if self.myControler is not None:
            self.myControler.profileView = None
        self.Close()

    def OnControlCFRAC(self):
        """ Control CFRAC versus Clouds and update the graphic if necessary """
        self.myProfile.ctrlCoherenceClouds()
        if 'CFRAC' in self.graphicPages:
            self.graphicPages['CFRAC'].UpdateData(self.myProfile['CFRAC'])

    def OnSaveProfile(self, e):
        """ return the name for a profile File """
        self.myProfile = self.GetProfile()
        self.OnControlCFRAC()
        fileSaved = self.OnSaveFile(e, "Save a profile")
        return fileSaved

    def OnApplyChange(self, e):
        """ get profile change from panels and return a profile """
        self.myProfile = self.GetProfile()
        self.OnControlCFRAC()
        return self.myProfile

    def OnSaveProfileAs(self, e):
        """ return the name for a profile File """
        self.myProfile = self.GetProfile()
        self.OnControlCFRAC()
        fileSaved = self.OnSaveFile(e, "Save a profile")
        return fileSaved

    def OnInsert(self, e):
        for item in self.graphicPages.keys():
            self.graphicPages[item].onInsert = True

    def OnRemove(self, e):
        for item in self.graphicPages.keys():
            self.graphicPages[item].onInsert = False

    def OnAddgas(self, e):
        self.changeItem(self.myProfile.gas_list, "no gas to add",
                        "Choose the gases to add to the profile",
                        self.myProfile.addGas, False)

    def OnAddAerosol(self, e):
        self.changeItem(self.myProfile.aerosol_list, "no aerosol to add",
                        "Choose the aerosols to add to the profile",
                        self.myProfile.addAerosol, False, kind="AEROSOLS")

    def OnAddCloud(self, e):
        self.changeItem(self.myProfile.cloud_list, "no layeritem to remove",
                        "Choose the clouds to add to the profile",
                        self.myProfile.addCloud, False, kind="CLOUDS")

    def OnRemovegas(self, e):
        self.changeItem(self.myProfile.gas_list, "no gas to remove",
                        "Choose the gases to remove from the profile",
                        self.myProfile.removeGas)

    def OnRemoveAerosol(self, e):
        self.changeItem(self.myProfile.aerosol_list, "no aerosol to remove",
                        "Choose the aerosols to remove from the profile",
                        self.myProfile.removeAerosol, kind="AEROSOLS")

    def OnRemoveCloud(self, e):
        self.changeItem(self.myProfile.cloud_list, "no layeritem to remove",
                        "Choose the clouds to remove from the profile",
                        self.myProfile.removeCloud, kind="CLOUDS")

    def changeItem(self, liste_item, message1, message2, action, remove=True,
                   kind="GASES"):
        """ perform the action on self.myProfile """
        wx.BeginBusyCursor()
        myList = []
        if remove:
            for item in liste_item:
                # cannot remove 'Q'
                if not (item == 'Q'):
                    if self.myProfile[item] is not None:
                        myList.append(item)
        else:
            for item in liste_item:
                if self.myProfile[item] is None:
                    myList.append(item)
        if myList == []:
            self.ShowInfoMessage(message1)
        else:
            list_to_change = self. ShowDialogList(message2, myList)
            if list_to_change is not None:
                for item in list_to_change:
                    action(item)
                    if remove:
                        self.deleteRightPage(item)
                    else:
                        self.addRightPage(item)
                # controle the coherence of clouds
                if kind == "CLOUDS":
                    print "control coherence layeritem recompute CFRAC"
                    self.OnControlCFRAC()
                self.replotLeftPage()
        wx.EndBusyCursor()

    def OnReplaceAerosolByClim(self, e):
        listItem = ["Continental clean", "Continental average",
                    "Continental polluted", "Urban",
                    "Desert", "Maritime clean", "Maritime polluted",
                    "Maritime tropical", "Arctic", "Antarctic"]
        selection = self.ShowDialogSingleChoice("Choose climatology", listItem)
        if selection is not None:
            self.write("replace aerosols by climatology " +
                       str(listItem[selection]))
            self.myProfile.replaceByAerosolClim(selection + 1)
        self.RePlotAll()

    def ShowDialogSingleChoice(self, label, listItem):
        selection = None
        dlg = wx.SingleChoiceDialog(self, label, "", listItem)

        if (dlg.ShowModal() == wx.ID_OK):
            selection = dlg.GetSelection()
            self.write("Selection:  " + str(selection))

        dlg.Destroy()
        return selection

    def ShowDialogList(self, label, listItem):
        strings = None
        dlg = wx.MultiChoiceDialog(self, label, "", listItem)

        if (dlg.ShowModal() == wx.ID_OK):
            selections = dlg.GetSelections()
            strings = [listItem[x] for x in selections]
            self.write("Selections:  " + str(strings))

        dlg.Destroy()
        return strings

    def OnMouseMove(self, e):
        """ print x y value of the left plot in the status bar  """
        pass

    def OnEditXAxe(self, e):
        self.GetProfile()
        mypage = self.nb.GetCurrentPage()
        item = mypage.theName
        (xmin, xmax) = mypage.axes.get_xlim()
        theScale = mypage.axes.get_xscale()
        if mypage.kind == "CLOUDS":
            theIncrement = 0.1
        else:
            theIncrement = 10
        dlg = AxesDialog(self, -1, "X axe edit",
                         text="enter X limits for " + item,
                         xmin0=xmin, xmax0=xmax, scale=theScale,
                         increment=theIncrement)

        if (dlg.ShowModal() == wx.ID_OK):
            (xmin, xmax) = dlg.GetSelections()
            self.graphicPages[item].SetXlimits(xmin=xmin, xmax=xmax)
            self.graphicPages[item].Update()
        dlg.Close()

    def OnChangeProfile(self, e):
        self.GetProfile()
        mypage = self.nb.GetCurrentPage()
        item = mypage.theName
        (pmax, pmin) = mypage.axes.get_ylim()
        min_offset = min(mypage.x)
        max_offset = max(mypage.x) * 2
        if mypage.kind == "CLOUDS":
            max_offset = 1
        dlg = ComputeDialog(self, -1, "change profile values",
                            text="enter value for " + item,
                            pmin0=pmin, pmax0=pmax,
                            limits_offset=(-min_offset, max_offset))
        if (dlg.ShowModal() == wx.ID_OK):
            (pmin, pmax, scale, offset) = dlg.GetSelections()
            self.write("update values for " + str(item))
            self.write("between pmin: %d  pmax: %d " % (pmin, pmax))
            self.write("between scale: %g  offset: %g" % (scale, offset))
            mypage.ComputeProfile(pmin, pmax, scale, offset)
            self.UpdateAllGraphics()
        dlg.Close()

    def MenuData(self):
        """ define the data for the menu
        """
        if self.edit:
            return(("&File",  # File Menu
                    ("Apply change", "Apply the profile ",
                     self.OnApplyChange, "applyChange", True, False),
                    ("Save profile", "Save the profile file",
                        self.OnSaveProfile, "saveProfile", True, False),
                    ("Save profile as", "Save the profile file",
                        self.OnSaveProfileAs, "saveProfileAs", True, False),
                    ('&Quit', 'Quit', self.OnClose, "quit", True, False)),
                   ("&Edit",  # Edit Menu
                    ("Undo", "Undo the last operation",
                     self.OnUndo, "undo", True, False),
                    ("Redo", "Redo the last operation",
                        self.OnRedo, "redo", True, False),
                    ("", "", "", "", True, False),
                    ("insert", "mode for add or modify points",
                        self.OnInsert, "insert", True, True),
                    ("remove", "mode for suppress points",
                        self.OnRemove, "remove", True, True),
                    ("", "", "", "", True, False),
                    ("change profile values", "chage the profile values",
                        self.OnChangeProfile, "changProfile", True, False),
                    ("edit x axe", "configure x axe",
                        self.OnEditXAxe, "edit x axe", True, False),

                    ("", "", "", "", True, False),
                    ("Add gas", "add a gas", self.OnAddgas,
                     "Add gas", True, False),
                    ("Remove gas", "remove a gas",
                        self.OnRemovegas, "Remove gas", True, False),
                    ("Add aerosol", "add an aerosol",
                        self.OnAddAerosol, "Add aerosol", True, False),
                    ("Remove aerosol", "remove an aerosol",
                        self.OnRemoveAerosol, "Remove aerosol", True, False),
                    ("Add cloud", "add a cloud",
                        self.OnAddCloud, "Add cloud", True, False),
                    ("Remove cloud", "remove a cloud",
                        self.OnRemoveCloud, "Remove cloud", True, False),
                    ("Replace Aerosol by clim", "Replace Aerosol by clim",
                     self.OnReplaceAerosolByClim, "Replace Aerosol by clim",
                     True, False)),
                   ("&Help",  # Help Menu
                    ("About", "About screen", self.OnAbout, "about",
                     True, False),
                    ("&Help", "Help", self.OnHelp, "help", True, False)))
        else:
            return(("&File",  # File Menu
                    ('&Quit', 'Quit', self.OnQuit, "quit", True, False)),
                   ("&Help",  # Help Menu
                    ("About", "About screen", self.OnAbout, "about",
                     True, False)))


if __name__ == "__main__":

    print "version matplotlib :", matplotlib.__version__
    p = rmodel.project.Project()

    p.openProfile(p.config.ENV["RTTOV_GUI_PREFIX"] +
                  '/rttov_tests/cldaer101lev_allgas.H5', 1)

    ex = wx.App()

    myProfileView = ProfileView(None, p.myProfile)
    myProfileView.Bind(wx.EVT_MENU, myProfileView.OnApplyChange,
                       myProfileView.items['applyChange'])

    profile = myProfileView.GetProfile()
    print "P"
    print profile['P']
    print "T"
    print profile["T"]
    print "loop"
    ex.MainLoop()
