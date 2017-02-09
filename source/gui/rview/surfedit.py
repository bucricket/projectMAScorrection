# -*- coding: utf-8 -*-

import util
try:
    import wx.grid
except ImportError:
    import sys
    sys.stderr.write('ERROR: wxPython is not installed\n')
    sys.exit(1)

import rmodel
import copy


class DataTable(wx.grid.PyGridTableBase):
    """ interface between the model and the view implements PyGridTableBase"""

    def __init__(self, data, dataName='EMIS'):
        """ init the data which could be either EMIS DATA ou REFL DATA """

        self.myData = copy.deepcopy(data)
        self.dataName = dataName
        wx.grid.PyGridTableBase.__init__(self)

        self.colsNames = {0: dataName + "_IN",
                          1: dataName + "_OUT",
                          2: "CALC" + dataName, }

    def GetAttr(self, row, col, kind):
        attr = wx.grid.GridCellAttr()
        if col == 2:
            attr.SetRenderer(wx.grid.GridCellBoolRenderer())
        else:
            if col == 1:
                attr.SetReadOnly()
                attr.SetBackgroundColour("light grey")
            else:
                attr.SetRenderer(wx.grid.GridCellFloatRenderer())
        return attr

    def GetNumberRows(self):
        name = self.dataName + '_IN'
        if (len(self.myData[name].shape) == 0):
            nbRows = 1
        else:
            nbRows = self.myData[name].shape[0]
        return nbRows

    def GetNumberCols(self):
        return 3

    def IsEmptyCell(self, row, col):
        return False

    def GetValue(self, row, col):
        value = self.myData[self.colsNames[col]][row]
        if col == 2:
            if value:
                return 1
            else:
                return 0
        else:
            return value

    def SetValue(self, row, col, value):
        if col == 2:
            if value == "0" or value == "" or value == "False":
                myvalue = False
            else:
                myvalue = True
        else:
            myvalue = value

        oldData = self.myData[self.colsNames[col]][row]
        try:
            self.myData[self.colsNames[col]][row] = myvalue
            if self.myData[self.colsNames[col]][row] > 1.:
                self.myData[self.colsNames[col]][row] = 1.
            if self.myData[self.colsNames[col]][row] < 0.:
                self.myData[self.colsNames[col]][row] = 0.
        except ValueError:
            self.myData[self.colsNames[col]][row] = oldData

    def GetColLabelValue(self, col):
        return self.colsNames[col]

    def UpdateData(self, data):
        self.myData = data


class GridView(util.GenericView):
    """ grid editor of the application for reflectance / emissivity """

    helpMessage = """

    This window allows you to modify reflectance or emissivity input values.
    Tip for check boxes : write 1 for check off , write 0 or nothing otherwise.

    """

    def __init__(self, parent, data, title, name, controller=None):
        self.myData = data
        util.GenericView.__init__(self, parent, title)
        self.myControler = controller
        sizer = wx.BoxSizer(wx.HORIZONTAL)
        self.SetSizer(sizer)
        self.CreateMenuBar()
        self.SetSize((400, 600))
        self.SetMinSize((400, 600))
        self.SetTitle(title)
        self.grid = wx.grid.Grid(self, -1)
        self.table = DataTable(data, name)
        self.grid.SetTable(self.table, True)
        self.grid.AutoSize()
        self.nbRow = self.table.GetNumberRows()
        sizer.Add(self.grid, 1, wx.EXPAND)
        self.Fit()
        self.sb = self.CreateStatusBar()
        self.Centre()
        self.Show(True)
        self.Bind(wx.EVT_CLOSE, self.OnQuit)

    def OnUpdateData(self, data, name):
        self.write("surfedit Update Data for " + name)

        if len(data[name + '_IN'].shape) == 0:
            newNbRow = 1
        else:
            newNbRow = data[name + '_IN'].shape[0]

        if newNbRow != self.nbRow:
            self.write(
                "surfedit nbRow has changed must make new DataTable for " +
                name)
            self.table = DataTable(data, name)
            self.grid.SetTable(self.table, True)
            self.nbRow = self.table.GetNumberRows()
            self.Refresh()
        else:
            self.table.UpdateData(data)
            self.Refresh()

    def OnDoubleClick(self, e):
        print "doubleClick"

    def OnApplyChange(self, e):
        self.myData = self.table.myData
        return self.myData

    def OnSaveSurface(self, e): pass

    def OnSaveSurfaceAs(self, e): pass

    def OnQuit(self, e):
        if self.myControler is not None:
            self.myControler.gridView = None
        self.Destroy()

    def MenuData(self):
        """ define the data for the menu
        """
        return(("&File",  # File Menu

                ("Apply change", "Apply change for the surface ",
                 self.OnApplyChange, "applyChange", True),
                ('&Quit', 'Quit', self.OnQuit, "quit", True)),
               ("&Help",  # Help Menu
                ("About", "About screen", self.OnAbout, "about", True),
                ("&Help", "Help", self.OnHelp, "help", True)))


if __name__ == "__main__":

    p = rmodel.project.Project()
    prefix = p.config.ENV['RTTOV_GUI_PREFIX']
    p.myCoeffs.nchannels = 1
    p.openSurface(prefix + "/rttov_tests/surface_atlas_1c.h5")
    emis = p.myEmissivity
    print "--emis---------"
    print emis
    print "-----------------"
    ex = wx.App()
    gv = GridView(None, emis, "Emissivity", 'EMIS')

    ex.MainLoop()
