

from cantera import *

import sys
if sys.version_info[0] == 3:
    from tkinter import *
else:
    from Tkinter import *

from .Units import temperature, pressure, density, specificEnergy, specificEntropy
from .UnitChooser import UnitVar
from .ThermoProp import ThermoProp
from .utilities import handleError

_PRESSURE = 1
_TEMPERATURE = 0
_DENSITY = 2
_INTENERGY = 3
_ENTHALPY = 4
_ENTROPY = 5

class ThermoFrame(Frame):
    def __init__(self,master,top):
        Frame.__init__(self,master)
        self.config(relief=GROOVE, bd=4)
        self.top = top
        self.mix = self.top.mix
        self.warn = 0
        self.internal = Frame(self)
        self.internal.pack(side=LEFT,anchor=N+W,padx=2,pady=2)
        self.controls=Frame(self.internal)
        self.controls.pack(side=LEFT,anchor=N+W,padx=4,pady=5)

        self.entries=Frame(self.internal)
        self.entries.pack(side=LEFT,anchor=N,padx=4,pady=2)
        self.makeEntries()
        self.makeControls()
        self.showState()

    def makeControls(self):
        Button(self.controls,text='Set State', width=15,
               command=self.setState).grid(column=0,row=0)
        self.equil = IntVar()
        self.equil.set(0)
        Button(self.controls,text='Equilibrate', width=15,
               command=self.eqset).grid(column=0,row=1)
##              Radiobutton(self.controls,text='Frozen',variable = self.equil,
##                             command=self.freeze,value=0).grid(column=0,row=2,sticky='W')
##              Radiobutton(self.controls,text='Equilibrium',
##                          variable=self.equil,
##                             command=self.eqset,value=1).grid(column=0,row=3,sticky='W')

    def eqset(self):
        self.equil.set(1)
        self.setState()
        self.equil.set(0)
        #if self.top.mixframe:
        #       self.top.mixframe.redo()

    def freeze(self):
        self.equil.set(0)
        if self.top.mixframe:
            self.top.mixframe.redo()

    def makeEntries(self):
        self.entries.pack()
        self.variable = {}
        self.prop = []
        props = ['Temperature', 'Pressure', 'Density',
                 'Internal Energy', 'Enthalpy', 'Entropy']
        units = [temperature, pressure, density, specificEnergy,
                 specificEnergy, specificEntropy]
        defaultunit = [0, 2, 0, 1, 1, 1]
        for i in range(len(props)):
            self.prop.append(ThermoProp(self.entries, self, i, props[i],
                                        0.0, units[i], defaultunit[i]))
            #self.prop[-1].entry.bind("<Any-Leave>",self.setState)
        self.last2 = self.prop[3]
        self.last1 = self.prop[2]
        self.prop[0].checked.set(1)
        self.prop[0].check()
        self.prop[1].checked.set(1)
        self.prop[1].check()
        self.showState()

    def checkTPBoxes(self):
        if not self.prop[0].isChecked():
            self.prop[0].checked.set(1)
            self.prop[0].check()
        if not self.prop[1].isChecked():
            self.prop[1].checked.set(1)
            self.prop[1].check()

    def showState(self):
        self.prop[_TEMPERATURE].set(self.mix.g.T)
        self.prop[_PRESSURE].set(self.mix.g.P)
        self.prop[_DENSITY].set(self.mix.g.density)
        self.prop[_INTENERGY].set(self.mix.g.int_energy_mass)
        self.prop[_ENTHALPY].set(self.mix.g.enthalpy_mass)
        self.prop[_ENTROPY].set(self.mix.g.entropy_mass)

    def setState(self,event=None):
        if event:
            self.warn = 0
        else:
            self.warn = 1
        self.top.mixfr.update()
        i = self.equil.get()
        optlist = ['frozen','equilibrium']
        opt = [optlist[i]]

        if self.prop[_PRESSURE].isChecked() \
           and self.prop[_TEMPERATURE].isChecked():
            self.mix.set(
                    temperature = self.prop[_TEMPERATURE].get(),
                    pressure = self.prop[_PRESSURE].get(),
                    equil=i)

        elif self.prop[_DENSITY].isChecked() \
             and self.prop[_TEMPERATURE].isChecked():
            self.mix.set(
                    temperature = self.prop[_TEMPERATURE].get(),
                    density = self.prop[_DENSITY].get(),
                    equil=i)

        elif self.prop[_ENTROPY].isChecked() \
             and self.prop[_PRESSURE].isChecked():
            self.mix.set(pressure = self.prop[_PRESSURE].get(),
                         entropy = self.prop[_ENTROPY].get(),
                         equil=i)

        elif self.prop[_ENTHALPY].isChecked() \
             and self.prop[_PRESSURE].isChecked():
            self.mix.set(pressure = self.prop[_PRESSURE].get(),
                         enthalpy = self.prop[_ENTHALPY].get(),
                         equil=i)

        elif self.prop[_INTENERGY].isChecked() \
             and self.prop[_DENSITY].isChecked():
            self.mix.set(density = self.prop[_DENSITY].get(),
                         intEnergy = self.prop[_INTENERGY].get(),
                         equil=i)
        else:
            if self.warn > 0:
                handleError("unsupported property pair")

        self.top.update()
