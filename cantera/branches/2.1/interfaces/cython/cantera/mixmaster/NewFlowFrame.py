from cantera import *

import sys
if sys.version_info[0] == 3:
    from tkinter import *
else:
    from Tkinter import *

from .SpeciesInfo import SpeciesInfo

_CUTOFF = 1.e-15
_ATOL = 1.e-15
_RTOL = 1.e-7

class NewFlowFrame(Frame):
    def __init__(self,master):
        Frame.__init__(self,master)
        self.config(relief=GROOVE, bd=4)
        self.app = self.master.app
        self.controls=Frame(self)
        self.hide = IntVar()
        self.hide.set(0)
        self.p = DoubleVar()
        #self.comp.set(1.0)
        self.controls.grid(column=1,row=0,sticky=W+E+N)
        #self.makeControls()
        mf = self.master

        e1 = Entry(self)
        e1.grid(column=0,row=0,sticky=E+W)
        e1['textvariable'] = self.p
        #e1.config(state=ENABLED)
        e1.config(relief=RIDGE)

##      def makeControls(self):
##              Radiobutton(self.controls,text='Moles',
##                          variable=self.comp,value=0,command=self.show).grid(column=0,row=0,sticky=W)
##              Radiobutton(self.controls,text='Mass',variable=self.comp,value=1,command=self.show).grid(column=0,row=1,sticky=W)
##              Radiobutton(self.controls,text='Concentration',variable=self.comp,value=2,command=self.show).grid(column=0,row=2,sticky=W)
##              Button(self.controls,text='Clear',command=self.zero).grid(column=0,row=4,sticky=W+E)
##              Button(self.controls,text='Normalize',command=self.norm).grid(column=0,row=5,sticky=W+E)
##              Checkbutton(self.controls,text='Hide Missing\nSpecies',
##                          variable=self.hide,onvalue=1,offvalue=0,command=self.master.redo).grid(column=0,row=3,sticky=W)


##      def makeControls(self):
##              self.c = CompFrame(self)
##              self.c.grid(column=1,row=0,sticky=E+W+N+S)

##      def redo(self):
##              self.update()
##              self.entries.destroy()
##              self.entries=Frame(self)
##              self.makeEntries()

##      def minimize(self,Event=None):
##          self.c.hide.set(1)
##          self.redo()
##          self.c.grid_forget()
##          self.entries.bind("<Double-1>",self.maximize)

##      def maximize(self,Event=None):
##          self.c.hide.set(0)
##          self.redo()
##          self.c.grid(column=1,row=0,sticky=E+W+N+S)
##          self.entries.bind("<Double-1>",self.minimize)


##      def makeEntries(self):
##              self.entries.grid(row=0,column=0,sticky=W+N+S+E)
##              self.entries.config(relief=GROOVE,bd=4)
##              DATAKEYS = self.top.species
##              self.variable = {}

##              n=0
##              ncol = 3
##              col = 0
##              row = 60

##              presbox =

##              for sp in DATAKEYS:
##                      s = sp # self.top.species[sp]
##                      k = s.index
##                      if row > 15:
##                              row = 0
##                              col = col + 2
##                              l = Label(self.entries,text='Species')
##                              l.grid(column=col,row=row,sticky=E+W)
##                              e1 = Entry(self.entries)
##                              e1.grid(column=col+1,row=row,sticky=E+W)
##                              e1['textvariable'] = self.var
##                              e1.config(state=DISABLED)
##                              e1.config(bg='lightyellow',relief=RIDGE)
##                              row = row + 1

##                      spname = s.name
##                      val = self.comp[k]
##                      if not self.c.hide.get() or val: showit = 1
##                      else:  showit = 0

##                      l=SpeciesInfo(self.entries,species=s,
##                                       text=spname,relief=FLAT,justify=RIGHT,
##                                    fg='darkblue')
##                      entry1 = Entry(self.entries)
##                      self.variable[spname] = DoubleVar()
##                      self.variable[spname].set(self.comp[k])
##                      entry1['textvariable']=self.variable[spname]
##                      entry1.bind('<Any-Leave>',self.up)
##                      if showit:
##                              l.grid(column= col ,row=row,sticky=E)
##                              entry1.grid(column=col+1,row=row)
##                              n=n+1
##                              row = row + 1
##                 if self.c.hide.get():
##                  b=Button(self.entries,height=1,command=self.maximize)
##              else:
##                  b=Button(self.entries,command=self.minimize)
##                 b.grid(column=col,columnspan=2, row=row+1)
