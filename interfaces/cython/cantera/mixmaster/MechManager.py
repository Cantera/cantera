from cantera import *

import sys
if sys.version_info[0] == 3:
    from tkinter import *
else:
    from Tkinter import *

from .ControlPanel import ControlWindow
from .ControlPanel import make_menu, menuitem_state
#from Cantera.Examples.Tk import _mechdir
import os

# automatically-loaded mechanisms
_autoload = [
    (' GRI-Mech 3.0', 'gri30.cti'),
    (' Air', 'air.cti'),
    (' H/O/Ar', 'h2o2.cti')
    ]

def testit():
    pass

class MechManager(Frame):

    def __init__(self,master,app):
        Frame.__init__(self,master)
        #self.config(relief=GROOVE, bd=4)
        self.app = app
        self.master = master
        self.mechindx = IntVar()
        self.mechindx.set(1)

        #m = Label(self, text = 'Loaded Mechanisms')
        #m.grid(column=0,row=0)
#         m.bind('<Double-1>',self.show)
#         self.mechindx.set(0)
        self.mechanisms = []
        self.mlist = [ [] ]
        i = 1
        #for m in self.mechanisms:
        #    self.mlist.append((m[0], self.setMechanism, 'check', self.mechindx, i))
        #    i = i + 1
        #self.mlist.append([])

        self.mechmenu = make_menu('Mixtures', self, self.mlist)
        self.mechmenu.grid(row=0,column=0,sticky=W)

        self.mfr = None

    def addMechanism(self, name, mech):
        self.mechanisms.append((name, mech))
        il = len(self.mechanisms)
        self.mlist[-1] = (name, self.setMechanism, 'check', self.mechindx, il)
        self.mlist.append([])

        self.mechmenu = make_menu('Mixtures', self, self.mlist)
        self.mechindx.set(il)
        self.mechmenu.grid(row=0,column=0,sticky=W)


    def delMechanism(self, mech):
        self.mechanisms.remove(mech)
        self.show()

##     def show(self,event=None):
##         print 'show'
##         if self.mfr:
##             self.mfr.destroy()
##         self.mfr = Frame(self)
##         self.mfr.grid(row=1,column=0)
##         self.mfr.config(relief=GROOVE, bd=4)
##         Label(self.mfr,text='jkl').grid(row=0,column=0)
##         i = 0
##         for name, mech in self.mechanisms:
##             Radiobutton(self.mfr, text=name, variable=self.mechindx,
##                         value = i,
##                         command=self.setMechanism).grid(row=i,column=0)
##             i = i + 1
##         print 'end'


    def setMechanism(self, event=None):
        i = self.mechindx.get()
        self.app.mech = self.mechanisms[i-1][1]
        self.app.makeMix()
        self.app.makeWindows()
