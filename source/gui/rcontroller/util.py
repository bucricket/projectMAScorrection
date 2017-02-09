# -*- coding: utf-8 -*-

try:
    import wx

except ImportError:
    import sys
    sys.stderr.write('ERROR: check your python installation\n')
    sys.stderr.write('ERROR: wxPython is not installed\n')
    sys.exit(1)
try:
    import datetime
except ImportError:
    import sys
    sys.stderr.write('ERROR: datetime is not installed\n')
    sys.exit(1)


import logging
if wx.VERSION[1] == 8:
    from wx.lib.pubsub import Publisher as pub

else:

    # wxPython 2.9 :
    # wx.lib.pubsub: Pusub now defaults to the new "kwarg" version of the API.
    # In order to continue using the original "arg1" API you
    # will need to import wx.lib.pubsub.setuparg1 before importing any other
    # pubsub modules.

    import wx.lib.pubsub.setuparg1
    from wx.lib.pubsub import pub


class UndoRedo:
    """ object keep to list for options and profile
        an undo list every modification of an object profile or options
        is saved in the list
        when the list are too long they must be cleaned """

    forUndo = {"profiles": [], "options": []}
    forRedo = {"profiles": [], "options": []}

    def __init__(self):
        pass

    def clean(self, typeItem):
        pass  # TODO

    def saveforUndo(self, item, typeItem):
        self.forUndo[typeItem].append(item)
        self.clean()

    def undo(self, typeItem):
        item = self.forUndo[typeItem].pop()
        del self.forUndo[typeItem][-1]
        self.forRedo[typeItem].append(item)
        return item

    def redo(self, typeItem):
        item = self.forRedo[typeItem].pop()
        del self.forRedo[typeItem][-1]
        self.forUndo[typeItem].append(item)
        return item


class GenericController (object):
    """ define Generic actions for all controllers
        SaveProfile etc
        project is a static variable : known for all instances
        undoredo is also a static variable """
    project = None
    undoRedo = None

    def __init__(self):
        if (self.undoRedo is None):
            self.undoRedo = UndoRedo()
        pub.subscribe(self.OnChanged, "Options CHANGED")
        pub.subscribe(self.OnChanged, "Profile CHANGED")
        self.controlerName = ""

    def prepLog(self):
        """ before any call to project methods put the cursor
            of the TextEdit of the main window at the end
            in order to prevent messed output """

        mainW = wx.GetTopLevelWindows()[0]
        if mainW is not None:
            if mainW.log is not None:
                mainW.log.SetInsertionPointEnd()

    def OnChanged(self, msg):
        print "message", msg

    def write(self, msg):
        self.prepLog()
        logging.info(self.controlerName + " " + msg)
        print datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"), msg

    def SaveProfile(self, e, parentFrame):
        """ Save the Profile (Profile+Option)  """
        try:
            fileName = None
            if (self.project.savedProfileFileName) is not None:
                fileName = self.project.savedProfileFileName
            else:
                fileName = parentFrame.OnSaveFile(e)
                self.project.savedProfileFileName = fileName
            if(fileName):
                self.project.saveProfile(fileName)
                print fileName, "saved"
        except IOError:
            print "Error while saving the profile file ", fileName

    def SaveSurface(self, e, parentFrame):
        """ Save the Surface File (Emissivity+Reflectance)  """
        print ("SaveSurface GenericController")
        try:
            fileName = None
            if (self.project.savedSurfaceFileName):
                fileName = self.project.savedSurfaceFileName
            else:
                fileName = self.project.surfaceFileName
            if(fileName):
                self.project.saveSurface(fileName)
                print fileName, "saved"
        except IOError:
            print "Error while saving the surface file "
