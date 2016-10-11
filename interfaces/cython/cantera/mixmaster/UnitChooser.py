# This file is part of Cantera. See License.txt in the top-level directory or
# at http://www.cantera.org/license.txt for license and copyright information.

import sys
if sys.version_info[0] == 3:
    from tkinter import *
else:
    from Tkinter import *

import re

class UnitVar(Frame):
    def __init__(self,master,unitmod,defaultunit=0):
        Frame.__init__(self,master)
        self.x = DoubleVar()
        self.xsi = 0.0
        self.x.set(0.0)
        self.unitmod = unitmod
        try:
            self.unitlist = self.unitmod.units
        except:
            self.unitlist = []
            unitlist=dir(self.unitmod)
            for it in unitlist:
                if it[0] != '_':
                    self.unitlist.append(it)
        self.v = Entry(self,textvariable=self.x)
        self.s = StringVar()
        tmp = re.sub('__',' / ',self.unitlist[defaultunit])
        self.s.set(tmp)
        self.conv = eval('self.unitmod.'+re.sub(' / ','__',self.s.get())).value
        self.u = Label(self)
        self.u.config(textvariable=self.s,fg='darkblue')
        self.u.bind('<Double-1>', self.select)
        self.u.bind('<Any-Enter>',self.highlight)
        self.u.bind('<Any-Leave>',self.nohighlight)
        self.v.grid(row=0,column=0)
        self.u.grid(row=0,column=1)

    def highlight(self, event=None):
        self.u.config(fg='yellow')

    def nohighlight(self, event=None):
        self.u.config(fg='darkblue')

    def select(self, event):
        self.new=Toplevel()
        self.new.title("Units")
        self.new.transient(self.master)
        self.new.bind("<Return>", self.finished,"+")

        r=0
        c=0
        for each in self.unitlist:
            if each[0] != '_' and each[:1] != '__' and each != 'SI':
                each = re.sub('__',' / ',each)
                Radiobutton(self.new,
                            text=each,
                            variable=self.u['textvariable'],
                            value=each,
                            command=self.update,
                            ).grid(column=c, row=r, sticky=W)
                r += 1
                if (r>10):
                    r=0
                    c += 1
                    r += 1

        b=Button(self.new,text='OK',command=self.finished, default=ACTIVE)
        b.grid(column=c, row=r)

        self.new.grab_set()
        self.new.focus_set()
        self.new.wait_window()

    def finished(self,event=None):
        self.new.destroy()

    def update(self):
        self.xsi = self.x.get() * self.conv
        self.conv = eval('self.unitmod.'+re.sub(' / ','__',self.s.get())).value
        self.x.set(self.xsi/self.conv)

    def get(self):
        self.xsi = self.x.get() * self.conv
        return self.xsi

    def set(self,value):
        self.xsi = value
        self.x.set(value/self.conv)
