import sys

if sys.version_info[0] == 3:
    from tkinter import *
else:
    from Tkinter import *

from cantera import *
import numpy as np

from .SpeciesInfo import SpeciesInfo
#from KineticsFrame import KineticsFrame

_CUTOFF = 1.e-15
_ATOL = 1.e-15
_RTOL = 1.e-7

class CompFrame(Frame):
    def __init__(self,master):
        Frame.__init__(self,master)
        self.config(relief=FLAT, bd=4)
        self.top = self.master.top
        self.controls=Frame(self)
        self.hide = IntVar()
        self.hide.set(0)
        self.comp = IntVar()
        self.comp.set(0)
        self.controls.grid(column=1,row=0,sticky=W+E+N)
        self.makeControls()
        mf = self.master

    def makeControls(self):
        Radiobutton(self.controls,text='Moles',
                    variable=self.comp,value=0,
                    command=self.show).grid(column=0,row=0,sticky=W)
        Radiobutton(self.controls,text='Mass',
                    variable=self.comp,value=1,
                    command=self.show).grid(column=0,row=1,sticky=W)
        Radiobutton(self.controls,text='Concentration',
                    variable=self.comp,value=2,
                    command=self.show).grid(column=0,row=2,sticky=W)
        Button(self.controls,text='Clear',
               command=self.zero).grid(column=0,row=4,sticky=W+E)
        Button(self.controls,text='Normalize',
               command=self.norm).grid(column=0,row=5,sticky=W+E)
        Checkbutton(self.controls,text='Hide Missing\nSpecies',
                    variable=self.hide,onvalue=1,
                    offvalue=0,command=self.master.redo).grid(column=0,
                                                              row=3,
                                                              sticky=W)

    def norm(self):
        mf = self.master
        mf.update()

        data = mf.comp
        sum = 0.0
        for sp in data:
            sum += sp
        for i in range(len(mf.comp)):
            mf.comp[i] /= sum
        self.show()

    def set(self):
        c = self.comp.get()
        mix = self.top.mix
        mf = self.master
        g = mix.g
        if c == 0:
            mix.setMoles(mf.comp)

        elif c == 1:
            mix.setMass(mf.comp)

        elif c == 2:
            pass
        self.top.thermo.setState()
        self.top.kinetics.show()

    def show(self):
        mf = self.master
        mf.active = self
        c = self.comp.get()
        mix = self.top.mix
        g = mix.g
        if c == 0:
            mf.var.set("Moles")
            #mf.data = spdict(mix.g, mix.moles())
            mf.comp = mix.moles()

        elif c == 1:
            mf.var.set("Mass")
            #mf.data = spdict(mix.g,mix.mass())
            mf.comp = mix.mass()

        elif c == 2:
            mf.var.set("Concentration")
            mf.comp = g.concentrations
            #mf.data = spdict(mix,mix,mf.comp)

        for s in mf.variable.keys():
            try:
                k = g.species_index(s)
                if mf.comp[k] > _CUTOFF:
                    mf.variable[s].set(mf.comp[k])
                else:
                    mf.variable[s].set(0.0)
            except:
                pass



    def zero(self):
        mf = self.master
        mf.comp *= 0.0
        self.show()



class MixtureFrame(Frame):
    def __init__(self,master,top):
        Frame.__init__(self,master)
        self.config(relief=GROOVE, bd=4)
        self.top = top
        self.top.mixframe = self
        self.g = self.top.mix.g
        #self.scroll = Scrollbar(self)
        self.entries=Frame(self)
        #self.scroll.config(command=self.entries.xview)
        #self.scroll.grid(column=0,row=1)
        self.var = StringVar()
        self.var.set("Moles")
        self.comp = np.array(self.top.mix.moles())
        self.names = self.top.mix.speciesNames()
        self.nsp = len(self.names)
        #self.data = self.top.mix.moleDict()
        self.makeControls()
        self.makeEntries()
        self.entries.bind('<Double-l>',self.minimize)
        self.ctype = 0
        self.newcomp = 0

    def makeControls(self):
        self.c = CompFrame(self)
        #self.k = KineticsFrame(self)
        self.active = self.c
        self.c.grid(column=1,row=0,sticky=E+W+N+S)
        #self.k.grid(column=2,row=0,sticky=E+W+N+S)

    def update(self):
        self.newcomp = 0
        for s in self.variable.keys():
            k = self.g.species_index(s)
            current = self.comp[k]
            val = self.variable[s].get()
            dv = abs(val - current)
            if dv > _RTOL*abs(current) + _ATOL:
                self.comp[k] = val
                self.newcomp = 1

    def show(self):
        self.active.show()
##              for k in range(self.nsp):
##                      sp = self.names[k]
##                      if self.comp[k] > _CUTOFF:
##                              self.variable[sp].set(self.comp[k])
##                      else:
##                              self.variable[sp].set(0.0)

    def redo(self):
        self.update()
        self.entries.destroy()
        self.entries=Frame(self)
        self.makeEntries()

    def minimize(self,Event=None):
        self.c.hide.set(1)
        self.redo()
        self.c.grid_forget()
        self.entries.bind("<Double-1>",self.maximize)

    def maximize(self,Event=None):
        self.c.hide.set(0)
        self.redo()
        self.c.grid(column=1,row=0,sticky=E+W+N+S)
        self.entries.bind("<Double-1>",self.minimize)

    def up(self, x):
        self.update()
        if self.newcomp:
            self.c.set()
            self.c.show()
            self.top.update()
            #thermo.showState()
            #self.top.kinetics.show()

    def makeEntries(self):
        self.entries.grid(row=0,column=0,sticky=W+N+S+E)
        self.entries.config(relief=FLAT,bd=4)
        DATAKEYS = self.top.species
        self.variable = {}

        n=0
        ncol = 3
        col = 0
        row = 60

        equil = 0
        if self.top.thermo:
            equil = self.top.thermo.equil.get()

        for sp in DATAKEYS:
            s = sp # self.top.species[sp]
            k = s.index
            if row > 25:
                row = 0
                col = col + 2
                l = Label(self.entries,text='Species')
                l.grid(column=col,row=row,sticky=E+W)
                e1 = Entry(self.entries)
                e1.grid(column=col+1,row=row,sticky=E+W)
                e1['textvariable'] = self.var
                e1.config(state=DISABLED)
                e1.config(bg='lightyellow',relief=RIDGE)
                row = row + 1

            spname = s.name
            val = self.comp[k]
            if not self.c.hide.get() or val: showit = 1
            else:  showit = 0

            l=SpeciesInfo(self.entries,species=s,
                          text=spname,relief=FLAT,justify=RIGHT,
                          fg='darkblue')
            entry1 = Entry(self.entries)
            self.variable[spname] = DoubleVar()
            self.variable[spname].set(self.comp[k])
            entry1['textvariable']=self.variable[spname]
            entry1.bind('<Any-Leave>',self.up)
            if showit:
                l.grid(column= col ,row=row,sticky=E)
                entry1.grid(column=col+1,row=row)
                n=n+1
                row = row + 1
            if equil == 1:
                entry1.config(state=DISABLED,bg='lightgray')
##                 if self.c.hide.get():
##                  b=Button(self.entries,height=1,command=self.maximize)
##              else:
##                  b=Button(self.entries,command=self.minimize)
##                 b.grid(column=col,columnspan=2, row=row+1)
