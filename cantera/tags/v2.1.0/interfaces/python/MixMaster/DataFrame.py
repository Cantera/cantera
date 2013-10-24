import os, math, string
from Tkinter import *
from Cantera import *
from Cantera.num import *
from Cantera import num
from tkFileDialog import askopenfilename
from GraphFrame import Graph
from DataGraph import DataGraph, plotLimits
from ControlPanel import make_menu


U_LOC = 1
V_LOC = 2
T_LOC = 3
P_LOC = 4
Y_LOC = 5

def testit(e = None):
    pass

class DataFrame(Frame):

    def __init__(self,master,top):
#           if master==None:
        self.master = Toplevel()
        self.master.protocol("WM_DELETE_WINDOW",self.hide)
        #else:
        #           self.master = master

        #self.vis = vis
        Frame.__init__(self,self.master)
        self.config(relief=GROOVE, bd=4)
        self.top = top
        self.mix = self.top.mix
        self.g = self.top.mix.g
        self.data = None
        self.zdata = None
        self.ydata = None
        self.plt = None
        #self.pltwhat = None
        self.datasets = []
        self.vars = []
        self.whichsoln = IntVar()
        self.loc = IntVar()
#        self.loc.set(1)
        self.lastloc = T_LOC # self.loc.get()
        self.datafile = StringVar()
        self.solnid = StringVar()
        self.gr = Frame(self)
        self.n = IntVar()

        self.scframe = Frame(self)
        self.sc = Scale(self.scframe, variable = self.n,
                        orient='horizontal',digits=0,
                        length=300,resolution=1,command=self.updateplot)
        self.sc.config(cnf={'from':0,'to':1})
        Label(self.scframe,text='Grid Point').grid(column=0,row=0)
        self.sc.grid(row=0,column=1)
        self.sc.bind('<ButtonRelease-1>',self.updateState)
        self.gr.grid(row=4,column=0,columnspan=10)

        self.grid(column=0,row=10)
        self.makeMenu()
        self.hide()


    def makeMenu(self):
        self.menubar = Frame(self, relief=GROOVE,bd=2)
        self.menubar.grid(row=0,column=0,sticky=N+W+E,columnspan=10)
        f = [('Open...',self.browseForDatafile)]
        #make_menu('File',self.menubar,items)
        make_menu('File',self.menubar,f)
        make_menu('Dataset',self.menubar,self.datasets)
        make_menu('Plot',self.menubar,self.vars)


    def browseForDatafile(self, e=None):
        pathname = askopenfilename(
                filetypes=[("Data Files", ("*.xml","*.csv","*.dat")),
                           ("All Files", "*.*")])
        if pathname:
            self.datafile.set(pathname)
        self.show()
        self.getSoln()

    def getSoln(self):
        fname = os.path.basename(self.datafile.get())
        ff = os.path.splitext(fname)
        self.datasets = []
        if len(ff) == 2 and (ff[1] == '.xml' or ff[1] == '.ctml'):

            x = XML.XML_Node('root',src=self.datafile.get())
            c = x.child('ctml')
            self.solns = c.children('simulation')
            if len(self.solns) > 1:
                i = 0
                for soln in self.solns:
                    self.datasets.append((soln['id'],self.pickSoln,
                                          'check',self.whichsoln,i))
                    i += 1
            self.solnid.set(self.solns[-1]['id'])
            self.soln = self.solns[-1]

            self.importData()

        elif len(ff) == 2 and (ff[1] == '.csv' or ff[1] == '.CSV'):
            self.importCSV()

        self.makeMenu()
        if self.loc.get() <= 0:
            self.loc.set(self.lastloc)


    def importCSV(self):
        self.lastloc = self.loc.get()
        if self.lastloc <= 0: self.lastloc = T_LOC
        self.vars = []
        self.zdata = None
        self.ydata = None
        if self.plt:
            self.plt.destroy()

        f = open(self.datafile.get(),'r')
        lines = f.readlines()
        vars = string.split(lines[0],',')
        nlines = len(lines)
        self.np = nlines - 1
        nv = len(vars)
        vv = []
        for n in range(nv):
            nm = vars[n].split()
            if n < nv - 1 or (len(nm) > 0 and nm[0].isalnum()):
                vv.append(nm[0])
            else:
                break
        nv = len(vv)
        vars = vv
        fdata = zeros((nv, self.np),'d')
        for n in range(self.np):
            v = string.split(lines[n+1],',')
            for j in range(nv):
                try:
                    fdata[j,n] = float(v[j])
                except:
                    fdata[j,n] = 0.0

        self.nsp = self.g.nSpecies()
        self.y = zeros(self.nsp,'d')
        self.data = zeros((self.nsp+6,self.np),'d')
        self.data[0,:] = fdata[0,:]
        self.label = ['-']*(self.nsp+6)
        self.label[0] = vars[0]
        w = []
        for n in range(1,nv-1):
            try:
                k = self.g.speciesIndex(vars[n])
            except:
                k = -1
            v2 = vars[n]
            if v2 == 'T':
                self.data[T_LOC,:] = fdata[n,:]
                self.label[T_LOC] = vars[n]
                w.append(('T', self.newplot, 'check', self.loc, T_LOC))
            elif v2 == 'P':
                self.data[P_LOC,:] = fdata[n,:]
                self.label[P_LOC] = vars[n]
                w.append((vars[n], self.newplot, 'check', self.loc, P_LOC))
            elif v2 == 'u':
                self.data[U_LOC,:] = fdata[n,:]
                self.label[U_LOC] = vars[n]
                w.append((vars[n], self.newplot, 'check', self.loc, U_LOC))
            elif v2 == 'V':
                self.data[V_LOC,:] = fdata[n,:]
                self.label[V_LOC] = vars[n]
                w.append((vars[n], self.newplot, 'check', self.loc, V_LOC))
            elif k >= 0:
                self.data[k+Y_LOC,:] = fdata[n,:]
                self.label[k+Y_LOC] = vars[n]
                w.append((vars[n], self.newplot, 'check', self.loc, k + Y_LOC))

        if self.data[P_LOC,0] == 0.0:
            self.data[P_LOC,:] = ones(self.np,'d')*OneAtm
            print 'Warning: no pressure data. P set to 1 atm.'

        self.sc.config(cnf={'from':0,'to':self.np-1})
        if self.loc.get() <= 0:
            self.loc.set(self.lastloc)
        self.updateplot()

        self.vars = w
        #self.makeMenu()
        self.scframe.grid(row=5,column=0,columnspan=10)


    def pickSoln(self):
        self.solnid.set(self.solns[self.whichsoln.get()]['id'])
        self.soln = self.solns[self.whichsoln.get()]
#        self.t.destroy()
        self.importData()


    def importData(self):

        self.lastloc = self.loc.get()
        if self.lastloc <= 0: self.lastloc = T_LOC
        self.vars = []
        self.zdata = None
        self.ydata = None
        if self.plt:
            self.plt.destroy()

        self.nsp = self.g.nSpecies()
        self.label = ['-']*(self.nsp + 6)

        self.y = zeros(self.nsp,'d')
        gdata = self.soln.child('flowfield/grid_data')
        xp = self.soln.child('flowfield').children('float')
        p = 0.0
        for x in xp:
            if x['title'] == 'pressure':
                p = float(x.value())
        fa = gdata.children('floatArray')
        self.np = int(fa[0]['size'])

        self.data = zeros((self.nsp+6,self.np),'d')
        w = []
        for f in fa:
            t = f['title']
            try:
                k = self.g.speciesIndex(t)
            except:
                k = -1
            v = XML.getFloatArray(f)
            if t == 'z' or t == 't':
                self.data[0,:] = v
                self.label[0] = t
            elif k >= 0:
                self.data[k + Y_LOC] = v
                self.label[k + Y_LOC] = t
                w.append((t, self.newplot, 'check', self.loc, k + Y_LOC))
            elif t == 'T':
                self.data[T_LOC,:] = v
                self.label[T_LOC] = t
                w.append((t, self.newplot, 'check', self.loc, T_LOC))
            elif t == 'u':
                self.data[U_LOC,:] = v
                self.label[U_LOC] = t
                w.append((t, self.newplot, 'check', self.loc, U_LOC))
            elif t == 'V':
                self.data[V_LOC,:] = v
                self.label[V_LOC] = t
                w.append((t, self.newplot, 'check', self.loc, V_LOC))

        self.data[P_LOC,:] = ones(self.np,'d')*p
        self.label[P_LOC] = 'P (Pa)'
        self.sc.config(cnf={'from':0,'to':self.np-1})
        if self.loc.get() <= 0:
            self.loc.set(self.lastloc)
        self.updateplot()

        self.vars = w
        self.scframe.grid(row=5,column=0,columnspan=10)


    def hide(self):
        #self.vis.set(0)
        self.master.withdraw()
        #if self.pltwhat: self.pltwhat.withdraw()

    def show(self, e=None):
        self.master.deiconify()

    def updateState(self, e=None):
        n = self.n.get()
        if self.plt: self.plt.update()

        for k in range(self.nsp):
            self.y[k] = self.data[k+Y_LOC,n]

        self.top.thermo.checkTPBoxes()
        self.mix.setMass(self.y)
        self.mix.set(temperature = self.data[T_LOC,n],
                     pressure = self.data[P_LOC,n])

        self.top.update()

    def newplot(self,e=0):
        loc = self.loc.get()
        self.zdata = self.data[0,:]
        self.ydata = self.data[loc,:]
        npts = len(self.zdata)

        ylog = 0
        if loc >= Y_LOC:
            for n in range(npts):
                if self.ydata[n] <= 0.0:
                    #print n, self.ydata[n]
                    self.ydata[n] = 1.0e-20
            self.ydata = num.log10(self.ydata)
            ylog = 1

        self.gdata = []
        zmin = self.zdata[0]
        zmax = self.zdata[-1]
        for n in range(npts):
            self.gdata.append((self.zdata[n],self.ydata[n]))

        ymin, ymax, dtick = plotLimits(self.ydata)
        if loc > 0:
            self.plt = DataGraph(self.gr,self.data, 0, loc,
                                 title='',
                                 label=(self.label[0],self.label[loc]),
                                 logscale=(0,ylog),
                                 pixelX=500,pixelY=400)
            self.plt.canvas.config(bg='white')
            self.plt.grid(row=1,column=0,columnspan=2,sticky=W+E)
            n = self.n.get()
            self.gdot = self.plt.plot(n,'red')

    def updateplot(self,event=None):
        if self.data == None: return

        if self.zdata == None:
            self.newplot()

        n = self.n.get()
        self.pnt = self.zdata[n], self.ydata[n]
        if hasattr(self, 'gdot'):
            self.plt.delete(self.gdot)
            self.gdot = self.plt.plot(n,'red')


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
