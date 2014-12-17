from __future__ import print_function
import sys

if sys.version_info[0] == 3:
    from tkinter import *
else:
    from Tkinter import *

from .ElementFrame import getElements
from .utilities import handleError
from cantera import *
from .config import *
from .SpeciesFrame import getSpecies

def testit():
    pass

class EditFrame(Frame):

    def redraw(self):
        try:
            self.eframe.destroy()
            self.sframe.destroy()
            self.rframe.destroy()
        except:
            pass
        self.addElementFrame()
        self.addSpeciesFrame()
        self.addReactionFrame()

    def __init__(self, master, app):
        Frame.__init__(self, master)
        self.mix = app.mix
        print(self.mix, dir(self.mix))
        self.app = app
        self.master = master
        self.master.title("Cantera Mechanism Editor")
        self.redraw()

    def addReactionFrame(self):
        self.rframe = Frame(self)
        self.rframe.config(relief=GROOVE,bd=4)
        self.rframe.grid(row=2,column=0,columnspan=10,sticky=E+W)
        b=Button(self.rframe,text='Reactions',command=testit)
        b.grid(column=5, row=0)

    def addElementFrame(self):
        self.eframe = Frame(self)
        self.eframe.config(relief=GROOVE,bd=4)
        self.eframe.grid(row=0,column=0,columnspan=10,sticky=E+W)
        self.element_labels = []
        n = 0
        for el in self.mix._mech.elementNames():
            x = Label(self.eframe,text=el,fg='darkblue')
            x.grid(column = n, row=0)
            self.element_labels.append(x)
            n = n + 1
        b=Button(self.eframe,text='Element',command=self.chooseElements, default=ACTIVE)
        b.grid(column=0, row=1, columnspan=10)


    def addSpeciesFrame(self):
        self.sframe = Frame(self)
        self.sframe.config(relief=GROOVE,bd=4)
        self.sframe.grid(row=1,column=0,columnspan=10,sticky=E+W)
        r = 0
        c = 0
        splist = self.app.species
        self.spcheck = []
        self.spec = []
        for i in range(self.app.mech.nSpecies()):
            self.spec.append(IntVar())
            self.spec[i].set(1)
            self.spcheck.append( Checkbutton(self.sframe,
                                             text=splist[i].name,
                                             variable=self.spec[i],
                                             onvalue = 1, offvalue = 0) )
            self.spcheck[i].grid(row = r, column = c, sticky = N+W)
            self.spcheck[i].bind("<Button-3>", self.editSpecies)
            c = c + 1
            if c > 4:
                c, r = 0, r + 1

    def getspecies(self):
        print(getSpecies(self.mix.speciesNames(),
                         self.mix.speciesNames()))

    def editSpecies(self, event=None):
        e = Toplevel(event.widget.master)
        w = event.widget
        txt = w.cget('text')
        sp = self.app.mix.species[txt]

        # name, etc.
        e1 = Frame(e, relief=FLAT)
        self.addEntry(e1,'Name',0,0,sp.name)
        self.addEntry(e1,'ID Tag',1,0,sp.id)
        self.addEntry(e1,'Phase',2,0,sp.phase)
        e1.grid(row=0,column=0)

        # elements
        elframe = Frame(e)
        elframe.grid(row=1,column=0)
        Label(elframe,text='Elemental Composition').grid(row=0,column=0,columnspan=2,sticky=E+W)

        i = 0
        for el in self.app.mech.elementNames():
            self.addEntry(elframe,el,i,0,self.mech.nAtoms(sp, el))
            i = i + 1

        # thermo
        thframe = Frame(e)
        thframe.grid(row=0,rowspan=2,column=1)
        thframe.config(relief=GROOVE,bd=4)
        i = 0
        Label(thframe,text='Thermodynamic Properties').grid(row=0,
              column=0, columnspan=4, sticky=E+W)
        if isinstance(sp.thermoParam(),NasaPolynomial):
            Label(thframe,text='Parametrization:').grid(row=1,column=1)
            self.addEntry(thframe,'',2,0,'NasaPolynomial')
            Label(thframe,text='Temperatures (min, mid, max):').grid(row=3,column=1)
            self.addEntry(thframe,'',4,0,str(sp.minTemp))
            self.addEntry(thframe,'',5,0,str(sp.midTemp))
            self.addEntry(thframe,'',6,0,str(sp.maxTemp))
        low = Frame(thframe)
        low.config(relief=GROOVE,bd=4)
        low.grid(row=1,rowspan=6,column=3,columnspan=2)
        Label(low,text='Coefficients for the Low\n Temperature Range').grid(row=0,column=0,columnspan=2,sticky=E+W)
        c = sp.thermoParam().coefficients(sp.minTemp)
        for j in range(7):
            self.addEntry(low,'a'+str(j),j+3,0,str(c[j]))
        high = Frame(thframe)
        high.config(relief=GROOVE,bd=4)
        high.grid(row=1,rowspan=6,column=5,columnspan=2)
        Label(high,text='Coefficients for the High\n Temperature Range').grid(row=0,column=0,columnspan=2,sticky=E+W)
        c = sp.thermoParam().coefficients(sp.maxTemp)
        for j in range(7):
            self.addEntry(high,'a'+str(j),j+3,0,str(c[j]))

        com = Frame(e)
        com.grid(row=10,column=0,columnspan=5)
        ok = Button(com,text='OK',default=ACTIVE)
        ok.grid(row=0,column=0)
        ok.bind('<1>',self.modifySpecies)
        Button(com,text='Cancel',command=e.destroy).grid(row=0,column=1)
        self.especies = e

    def modifySpecies(self,event=None):
        button = event.widget
        e = self.especies
        for fr in e.children.values():
            for item in fr.children.values():
                try:
                    print(item.cget('selection'))
                except:
                    pass
        e.destroy()

    def addEntry(self,master,name,row,column,text):
        if name:
            Label(master, text=name).grid(row=row, column=column)
        nm = Entry(master)
        nm.grid(row=row, column=column+1)
        nm.insert(END,text)

    def chooseElements(self):
        oldel = self.mix.g.elementNames()
        newel = getElements(self.mix.g.elementNames())
        removeList = []
        for el in oldel:
            if not el in newel:
                removeList.append(el)
        #self.app.mech.removeElements(removeList)
        addList = []
        for el in newel:
            if not el in oldel:
                addList.append(el)
        #self.app.mech.addElements(addList)
        try:
            self.redraw()
            self.app.makeWindows()
        except:
            handleError('Edit err')

        self.app.mix = IdealGasMixture(self.app.mech)
        self.mix = self.app.mix
        nn = self.mix.speciesList[0].name
        self.mix.set(temperature = 300.0, pressure = 101325.0, moles = {nn:1.0})
        for label in self.element_labels:
            label.destroy()
        self.element_labels = []
        n = 0
        for el in self.mix._mech.elementList():
            x = Label(self.eframe,text=el.symbol(),fg='darkblue')
            x.grid(column = n, row=0)
            self.element_labels.append(x)
            n = n + 1

        self.app.makeWindows()
