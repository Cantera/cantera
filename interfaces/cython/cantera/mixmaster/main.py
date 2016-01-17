#############################################################################
#
#  MixMaster
#
#############################################################################

# options
_app_title = 'MixMaster'
_app_version = '1.0'


# functionality imports
import sys
if sys.version_info[0] == 3:
    from tkinter import *
    from tkinter import messagebox
    from tkinter.filedialog import askopenfilename
else:
    from Tkinter import *
    import tkMessageBox as messagebox
    from tkFileDialog import askopenfilename

import sys, os, string

# Cantera imports
from cantera import *
from numpy import zeros
from . import utilities

# local imports
from .TransportFrame import TransportFrame
from .CompositionFrame import MixtureFrame
from .ThermoFrame import ThermoFrame
from .ImportFrame import ImportFrame
from .DataFrame import DataFrame
from .KineticsFrame import SpeciesKineticsFrame, ReactionKineticsFrame, ReactionPathFrame
#from Edit import EditFrame
from .MechManager import MechManager, _autoload
from .UnitChooser import UnitVar
from .ControlPanel import ControlWindow
from .ControlPanel import make_menu, menuitem_state
from .Mix import Mix, Species

def testit():
    return


class MixMaster:

    def stop(self):
        sys.exit(0)

    def openmech(self):
        pathname = askopenfilename(filetypes=[("Cantera Input Files", "*.cti"),
                                              ("XML Files", "*.xml *.ctml"),
                                              ("All Files", "*.*")])
        if pathname:
            self.loadmech('', pathname)


    def loadmech(self, mechname, pathname, mw=1):

        p = os.path.normpath(os.path.dirname(pathname))
        self.fname = os.path.basename(pathname)
        ff = os.path.splitext(self.fname)

        try:
            self.mech = Solution(pathname)
            self.mechname = ff[0]

        except:
            utilities.handleError('could not create gas mixture object: '
                                  +ff[0]+'\n')
            self.mechname = 'Error'
            return

        self.makeMix()

        if not mechname:
            mechname = self.mechname

        self.mechframe.addMechanism(mechname, self.mech)
        if mw==1:
            self.makeWindows()


    def addWindow(self, name, w):
        """Add a new window, or replace an existing one."""
        wstate = ''
        if name in self._windows:
            try:
                wstate = self._windows[name].master.state()
                self._windows[name].master.destroy()
            except:
                pass
        else:
            wstate = 'withdrawn'
        self._windows[name] = w
        self._vis[name] = IntVar()
        if wstate == 'withdrawn':
            self._windows[name].master.withdraw()
        else:
            self._windows[name].show()



    def update(self):
        """Update all windows to reflect the current mixture state."""
        for w in self._windows.keys():
            try:
                m = self._windows[w].master
                if m.state() != 'withdrawn':
                    self._windows[w].show()
            except:
                pass
        self.thermo.showState()
        self.mixfr.show()



    def makeMix(self):
        self.mix = Mix(self.mech)
        nsp = self.mech.n_species
        self.species = []
        nm = self.mech.species_names

        for k in range(nsp):
            self.species.append(Species(self.mech, nm[k]))

        x = self.mech.X
        self.mix.setMoles(x)
        self.mix.set(temperature = self.mech.T,
                     pressure = self.mech.P)



    def __init__(self, master=None):
        """
        Create a new Application instance.
        Usually this is only called once.
        """
        if master:
            self.master = master
        else:
            t = Tk()
            self.master = t

        self._windows = {}
        self._vis = {}
        self.windows = []

        self.cwin = ControlWindow(_app_title,self.master)
        self.cwin.master.resizable(FALSE,FALSE)

        self.menubar = Frame(self.cwin, relief=GROOVE,bd=2)
        self.menubar.grid(row=0,column=0,sticky=N+W+E)

        self.mixfr = None
        self.thermo = None
        self.transport = None
        self.kinetics = None
        self.rxndata = None
        self.rxnpaths = None
        self.edit = None
        self.fname = None

        self.mechframe = MechManager(self.cwin, self)
        self.mechframe.grid(row=1,column=0,sticky=N+W)

        fileitems = [('Load Mixture...', self.openmech),
                     ('Import Mechanism File...',self.importfile),
                     'separator',
                     ('Load Data File...',self.showdata),
                     'separator',
                     ('Exit', self.stop),
                     []
                     ]
        self.filemenu = make_menu('File', self.menubar, fileitems)


        self.vtherm = IntVar()
        self.vcomp = IntVar()
        self.vtran = IntVar()
        self.vkin = IntVar()
        self.vrxn = IntVar()
        self.vrxn.set(0)
        self.vtherm.set(1)
        self.vedit = IntVar()


        dataitems = [(' Import Flame Data', testit),
                     (' Import CSV Data', testit),
                     []]
        #self.datamenu = make_menu('Data', self.menubar, dataitems)


        #toolitems = [(' Convert...', self.importfile),
        #             []]
        #self.toolmenu = make_menu('Tools', self.menubar, toolitems)


        w = [(' Thermodynamic State', self.showthermo, 'check', self.vtherm),
             (' Composition', self.showcomp, 'check', self.vcomp),
             'separator',
             (' Kinetics', self.showkinetics, 'check', self.vkin),
             (' Reactions...', self.showrxns),
             (' Reaction Paths...', self.showrpaths),
             []]

        self.viewmenu = make_menu('Windows', self.menubar, w)

        self.helpmenu = make_menu('Help', self.menubar,
                                  [('About '+_app_title+'...', self.aboutmix),
                                   ('About Cantera...', testit),
                                   []

                                   ])

        # load the preloaded mechanisms
        for m in _autoload:
            self.loadmech(m[0],m[1],0)

        self.makeWindows()
        self.addWindow('import',ImportFrame(self))

        self.vtherm.set(1)
        self.showthermo()
##         self.vcomp.set(1)
##         self.showcomp()

        self.master.iconify()
        self.master.update()
        self.master.deiconify()
        self.cwin.mainloop()


    def importfile(self):
        #self.vimport.set(1)
        w = self._windows['import']
        w.show()


    def makeWindows(self):
#        if self.mixfr:
        for w in self.windows:
            try:
                w.destroy()
            except:
                pass

        fr = [MixtureFrame, ThermoFrame, TransportFrame]

        self.mixfr = MixtureFrame(self.cwin, self)
        self.thermo = ThermoFrame(self.cwin, self)

#        self.transport = TransportFrame(self.cwin, self)
        self.kinetics = SpeciesKineticsFrame(self.cwin, self)

        self.addWindow('rxndata',ReactionKineticsFrame(self.vrxn, self))
        self.addWindow('rxnpaths',ReactionPathFrame(self))
        self.addWindow('dataset',DataFrame(None, self))


        #self.edit = EditFrame(t, self)

        self.windows = [self.mixfr,
                        self.thermo, self.transport,
                        self.kinetics]

        self.showthermo()
        self.showcomp()
        #self.showtransport()
        self.showkinetics()
        #self.showrxns()
        #self.showrpaths()
        #self.showdata()

        if self.mech:
            self.mechframe.grid(row=1,column=0)
        else:
            self.mechframe.grid_forget()
        #self.showedit()

    def show(self, frame, vis, row, col):
        if vis:
            frame.grid(row=row,column=col,sticky=N+E+S+W)
        else:
            frame.grid_forget()

    def showthermo(self):
        if self.thermo:
            self.show(self.thermo, self.vtherm.get(), 7, 0)

    def showcomp(self):
        if self.mixfr:
            self.show(self.mixfr, self.vcomp.get(), 8, 0)

    def showkinetics(self):
        if self.kinetics:
            self.show(self.kinetics, self.vkin.get(), 10, 0)

    def showrxns(self):
        self._windows['rxndata'].show()

    def showrpaths(self):
        self._windows['rxnpaths'].show()

    def showdata(self):
        self._windows['dataset'].browseForDatafile()

    def aboutmix(self):

        m = messagebox.showinfo(title = 'About MixMaster',
                                message = """
                     MixMaster

                    version """+_app_version+"""

written by:

Prof. David G. Goodwin
California Institute of Technology

copyright 2003
California Institute of Technology
     """)



if __name__ == "__main__":
    MixMaster()
