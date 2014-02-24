#
#  function getElements displays a periodic table, and returns a list of
#  the selected elements
#

from Tkinter import *
from types import *
import tkMessageBox

from Cantera import *

class SpeciesFrame(Frame):

    def __init__(self, master, speciesList = [], selected=[]):
        Frame.__init__(self,master)
        self.master = master
        self.control = Frame(self)
        self.species = {}
        for sp in speciesList:
            self.species[sp.name] = sp

        self.control.config(relief=GROOVE,bd=4)
        Button(self.control, text = 'Display',command=self.show).pack(fill=X,pady=3, padx=10)
        Button(self.control, text = 'Clear',command=self.clear).pack(fill=X,pady=3, padx=10)
        Button(self.control, text = '  OK  ',command=self.get).pack(side=BOTTOM,
                                                                    fill=X,pady=3, padx=10)
        Button(self.control, text = 'Cancel',command=self.master.quit).pack(side=BOTTOM,
                                                                            fill=X,pady=3, padx=10)
        self.entries = Frame(self)
        self.entries.pack(side=LEFT)
        self.control.pack(side=RIGHT,fill=Y)
        self.c = {}
        self.selected = selected
        n=0
        ncol = 8
        rw = 1
        col = 0
        list = self.species.values()
        list.sort()
        for sp in list:
            el = sp.name
            self.species[el] = Frame(self.entries)
            self.species[el].config(relief=GROOVE, bd=4, bg=self.color(el))
            self.c[el] = Button(self.species[el],text=el,bg=self.color(el),width=6,relief=FLAT)
            self.c[el].pack()
            self.c[el].bind("<Button-1>",self.setColors)
            self.species[el].grid(row= rw, column = col,sticky=W+N+E+S)
            col = col + 1
            if col > ncol:
                rw = rw + 1
                col = 0
        Label(self.entries,text='select the species to be included, and then press OK.\nTo view the properties of the selected species, press Display ').grid(row=0, column=2, columnspan=10, sticky=W)


    def select(self, el):
        self.c[el]['relief'] = RAISED
        self.c[el]['bg'] = self.color(el, sel=1)

    def deselect(self, el):
        self.c[el]['relief'] = FLAT
        self.c[el]['bg'] = self.color(el, sel=0)

    def selectSpecies(self,splist):
        for sp in splist:
            spname = sp.name
            self.select(spname)

    def setColors(self,event):
        el = event.widget['text']
        if event.widget['relief'] == RAISED:
            event.widget['relief'] = FLAT
            back = self.color(el, sel=0)
            fore = '#ffffff'
        elif event.widget['relief'] == FLAT:
            event.widget['relief'] = RAISED
            fore = '#000000'
            back = self.color(el, sel=1)
        event.widget['bg'] = back
        event.widget['fg'] = fore

    def color(self, el, sel=0):
        _normal = ['#88dddd','#005500','#dd8888']
        _selected = ['#aaffff','#88dd88','#ffaaaa']
        #row, column = _pos[el]
        if sel: list = _selected
        else: list = _normal
        return list[1]
        #if column < 3:
        #    return list[0]
        #elif column > 12:
        #    return list[1]
        #else:
        #    return list[2]

    def show(self):
        selected = []
        for sp in self.species.values():
            if self.c[sp.name]['relief'] == RAISED:
                selected.append(sp)
        #showElementProperties(selected)

    def get(self):
        self.selected = []
        for sp in self.species.values():
            if self.c[sp.name]['relief'] == RAISED:
                self.selected.append(sp)
        #self.master.quit()'
        self.master.destroy()

    def clear(self):
        for sp in self.species.values():
            self.c[sp]['bg'] = self.color(sp, sel=0)
            self.c[sp]['relief'] = FLAT

## class ElementPropertyFrame(Frame):
##     def __init__(self,master,ellist):
##         Frame.__init__(self,master)
##         n = 1
##         ellist.sort()
##         Label(self,text='Name').grid(column=0,row=0,sticky=W+S,padx=10,pady=10)
##         Label(self,text='Atomic \nNumber').grid(column=1,row=0,sticky=W+S,padx=10,pady=10)
##         Label(self,
##               text='Atomic \nWeight').grid(column=2,
##                                            row=0,
##                                            sticky=W+S,
##                                            padx=10,
##                                            pady=10)
##         for el in ellist:
##             Label(self,
##                   text=el.name).grid(column=0,
##                                      row=n,
##                                      sticky=W,
##                                      padx=10)
##             Label(self,
##                   text=`el.atomicNumber`).grid(column=1,
##                                                row=n,
##                                                sticky=W,
##                                                padx=10)
##             Label(self,
##                   text=`el.atomicWeight`).grid(column=2,
##                                                row=n,
##                                                sticky=W,
##                                                padx=10)
##             n = n + 1


# utility functions

def getSpecies(splist=[],selected=[]):
    master = Toplevel()
    master.title('Species')
    t = SpeciesFrame(master,splist,selected)
    if splist: t.selectSpecies(splist)
    t.pack()
    t.focus_set()
    t.grab_set()
    t.wait_window()
    try:
        master.destroy()
    except TclError:
        pass
    return t.selected


# display table of selected element properties in a window
def showElementProperties(ellist):
    m = Tk()
    m.title('Element Properties')
    elem = []
    ElementPropertyFrame(m, ellist).pack()


if __name__ == "__main__":
    print getSpecies()
