import sys
if sys.version_info[0] == 3:
    from tkinter import *
else:
    from Tkinter import *

import re, math
from cantera import *
from .Units import temperature, specificEnergy, specificEntropy
from .UnitChooser import UnitVar
from .GraphFrame import Graph

def testit():
    pass

class SpeciesInfo(Label):
    def __init__(self,master,phase=None,species=None,**opt):
        Label.__init__(self,master,opt)
        self.sp = species
        self.phase = phase
        self.bind('<Double-1>', self.show)
        self.bind('<Button-3>', self.show)
        self.bind('<Any-Enter>', self.highlight)
        self.bind('<Any-Leave>', self.nohighlight)


    def highlight(self, event=None):
        self.config(fg='yellow')

    def nohighlight(self, event=None):
        self.config(fg='darkblue')

    def show(self, event):
        self.new=Toplevel()
        self.new.title(self.sp.symbol)
        #self.new.transient(self.master)
        self.new.bind("<Return>", self.update,"+")
        self.cpr = 0.0
        self.t = 0.0
        self.cpl = 0.0
        self.tl = 0.0
        self.cpp = [[(0.0, 0.0, 'red')]]

        # elemental composition
        self.eframe = Frame(self.new)
        self.eframe.config(relief=GROOVE,bd=4)
        self.eframe.grid(row=0,column=0,columnspan=10,sticky=E+W)
        r = 1
        Label(self.eframe,text='Atoms:')\
                               .grid(row=0,column=0,sticky=N+W)
        for el, c in self.sp.composition():
            Label(self.eframe,text=repr(int(c))+' '+el).grid(row=0,column=r)
            r = r + 1


        # thermodynamic properties
        self.thermo = Frame(self.new)
        self.thermo.config(relief=GROOVE,bd=4)
        self.thermo.grid(row=1,column=0,columnspan=10,sticky=N+E+W)
        Label(self.thermo,text = 'Standard Heat of Formation at 298 K: ').grid(row=0, column=0, sticky=W)
        Label(self.thermo,text = '%8.2f kJ/mol' % (self.sp.hf0*1.0e-6)).grid(row=0, column=1, sticky=W)
        Label(self.thermo,text = 'Molar Mass: ').grid(row=1, column=0, sticky=W)
        Label(self.thermo,text = self.sp.molecularWeight).grid(row=1, column=1, sticky=W)
        labels = ['Temperature', 'c_p', 'Enthalpy', 'Entropy']
        units = [temperature, specificEntropy, specificEnergy, specificEntropy]
        whichone = [0, 1, 1, 1]

        r = 2
        self.prop = []
        for prop in labels:
            Label(self.thermo,text=prop).grid(row=r,column=0,sticky=W)
            p = UnitVar(self.thermo,units[r-2],whichone[r-2])
            p.grid(row=r,column=1,sticky=W)
            p.v.config(state=DISABLED,bg='lightgray')
            self.prop.append(p)
            r = r + 1

        tmin = self.sp.minTemp
        tmax = self.sp.maxTemp
        cp = self.sp.cp_R(tmin)
        hh = self.sp.enthalpy_RT(tmin)
        ss = self.sp.entropy_R(tmin)

        self.prop[0].bind("<Any-Enter>", self.decouple)
        self.prop[0].bind("<Any-Leave>", self.update)
        self.prop[0].bind("<Key>", self.update)
        self.prop[0].v.config(state=NORMAL,bg='white')
        self.prop[0].set(300.0)

        self.graphs = Frame(self.new)
        self.graphs.config(relief=GROOVE,bd=4)
        self.graphs.grid(row=2,column=0,columnspan=10,sticky=E+W)

        self.cpdata = []
        self.hdata = []
        self.sdata = []
        t = tmin
        n = int((tmax - tmin)/100.0)
        while t <= tmax:
            self.cpdata.append((t,self.sp.cp_R(t)))
            self.hdata.append((t,self.sp.enthalpy_RT(t)))
            self.sdata.append((t,self.sp.entropy_R(t)))
            t = t + n

        # specific heat

        Label(self.graphs,text='c_p/R').grid(row=0,column=0,sticky=W+E)
        ymin, ymax, dtick = self.plotLimits(self.cpdata)
        self.cpg = Graph(self.graphs,'',tmin,tmax,ymin,ymax,
                         pixelX=150,pixelY=150)
        self.cpg.canvas.config(bg='white')
        self.cpg.grid(row=1,column=0,columnspan=2,sticky=W+E)
        self.ticks(ymin, ymax, dtick, tmin, tmax, self.cpg)

        # enthalpy
        Label(self.graphs,text='enthalpy/RT').grid(row=0,column=3,sticky=W+E)
        ymin, ymax, dtick = self.plotLimits(self.hdata)
        self.hg = Graph(self.graphs,'',tmin,tmax,ymin,ymax,
                        pixelX=150,pixelY=150)
        self.hg.canvas.config(bg='white')
        self.hg.grid(row=1,column=3,columnspan=2,sticky=W+E)
        self.ticks(ymin, ymax, dtick, tmin, tmax, self.hg)

        # entropy
        Label(self.graphs,text='entropy/R').grid(row=0,column=5,sticky=W+E)
        ymin, ymax, dtick = self.plotLimits(self.sdata)
        self.sg = Graph(self.graphs,'',tmin,tmax,ymin,ymax,
                        pixelX=150,pixelY=150)
        self.sg.canvas.config(bg='white')
        self.sg.grid(row=1,column=5,columnspan=2,sticky=W+E)
        self.ticks(ymin, ymax, dtick, tmin, tmax, self.sg)

        n = int((tmax - tmin)/100.0)
        t = tmin
        self.cpp = []

        for t, cp in self.cpdata:
            self.cpg.join([(t,cp,'red')])
        for t, h in self.hdata:
            self.hg.join([(t,h,'green')])
        for t, s in self.sdata:
            self.sg.join([(t,s,'blue')])

        self.cpdot = self.cpg.plot(tmin,cp,'red')
        self.hdot = self.hg.plot(tmin,hh,'green')
        self.sdot = self.sg.plot(tmin,ss,'blue')

        b=Button(self.new,text=' OK ',command=self.finished, default=ACTIVE)
        #ed=Button(self.new,text='Edit',command=testit)
        b.grid(column=0, row=4,sticky=W)
        #ed.grid(column=1,row=4,sticky=W)

        self.scfr = Frame(self.new)
        self.scfr.config(relief=GROOVE,bd=4)
        self.scfr.grid(row=3,column=0,columnspan=10,sticky=N+E+W)
        self.sc = Scale(self.scfr,command=self.update,variable = self.prop[0].x,
                   orient='horizontal',digits=7,length=400)
        self.sc.config(cnf={'from':tmin,'to':tmax})
        self.sc.bind('<Any-Enter>',self.couple)
        self.scfr.bind('<Any-Leave>',self.decouple)
        self.sc.grid(row=0,column=0,columnspan=10)

    def decouple(self,event=None):
        d = DoubleVar()
        xx = self.prop[0].get()
        d.set(xx)
        self.sc.config(variable = d)

    def couple(self,event=None):
        self.sc.config(variable = self.prop[0].x)
        #self.update()

    def update(self,event=None):
        try:
            tmp = self.prop[0].get()
            cnd = self.sp.cp_R(tmp)
            cc = cnd*GasConstant
            self.prop[1].set(cc)
            hnd = self.sp.enthalpy_RT(tmp)
            hh = hnd*tmp*GasConstant
            self.prop[2].set(hh)
            snd = self.sp.entropy_R(tmp)
            ss = snd*tmp*GasConstant
            self.prop[3].set(ss)


            self.cppoint = tmp, cnd
            self.hpoint = tmp, hnd
            self.spoint = tmp, snd
            if hasattr(self, 'cpdot'):
                self.cpg.delete(self.cpdot)
                self.cpdot = self.cpg.plot(self.cppoint[0], self.cppoint[1],'red')
                self.hg.delete(self.hdot)
                self.hdot = self.hg.plot(self.hpoint[0], self.hpoint[1],'green')
                self.sg.delete(self.sdot)
                self.sdot = self.sg.plot(self.spoint[0], self.spoint[1],'blue')
        except:
            pass

    def plotLimits(self, xy):
        ymax = -1.e10
        ymin = 1.e10
        for x, y in xy:
            if y > ymax: ymax = y
            if y < ymin: ymin = y

        dy = abs(ymax - ymin)
        if dy < 0.2*ymin:
            ymin = ymin*.9
            ymax = ymax*1.1
            dy = abs(ymax - ymin)
        else:
            ymin = ymin - 0.1*dy
            ymax = ymax + 0.1*dy
            dy = abs(ymax - ymin)

        p10 = math.floor(math.log10(0.1*dy))
        fctr = math.pow(10.0, p10)
        mm = [2.0, 2.5, 2.0]
        i = 0
        while dy/fctr > 5:
            fctr = mm[i % 3]*fctr
            i = i + 1
        ymin = fctr*math.floor(ymin/fctr)
        ymax = fctr*(math.floor(ymax/fctr + 1))
        return (ymin, ymax, fctr)

    def ticks(self, ymin, ymax, dtick, tmin, tmax, plot):
        ytick = ymin
        eps = 1.e-3
        while ytick <= ymax:
            if abs(ytick) < eps:
                plot.join([(tmin, ytick, 'gray')])
                plot.join([(tmax, ytick, 'gray')])
                plot.last_points = []
            else:
                plot.join([(tmin, ytick, 'gray')])
                plot.join([(tmin + 0.05*(tmax - tmin), ytick, 'gray')])
                plot.last_points = []
                plot.join([(2.0*tmax, ytick, 'gray')])
                plot.join([(tmax - 0.05*(tmax - tmin), ytick, 'gray')])
                plot.last_points = []

            ytick = ytick + dtick

    def finished(self,event=None):
        self.new.destroy()
