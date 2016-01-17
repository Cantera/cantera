from __future__ import print_function

import os, math, sys

if sys.version_info[0] == 3:
    from tkinter import *
else:
    from Tkinter import *

from cantera import *

from .SpeciesInfo import SpeciesInfo
import webbrowser

_CUTOFF = 1.e-15
_ATOL = 1.e-15
_RTOL = 1.e-7

def showsvg():
    f = open('_rp_svg.html','w')
    f.write('<embed src="rxnpath.svg" name="rxnpath" height=500\n')
    f.write('type="image/svg-xml" pluginspage="http://www.adobe.com/svg/viewer/install/">\n')
    f.close()
    webbrowser.open('file:///'+os.getcwd()+'/_rp_svg.html')

def showpng():
    f = open('_rp_png.html','w')
    f.write('<img src="rxnpath.png" height=500/>\n')
    f.close()
    webbrowser.open('file:///'+os.getcwd()+'/_rp_png.html')


class KineticsFrame(Frame):
    def __init__(self,master):
        Frame.__init__(self,master)
        self.config(relief=FLAT, bd=4)
        self.top = self.master.top
        self.controls=Frame(self)
        self.hide = IntVar()
        self.hide.set(0)
        self.comp = IntVar()
        self.comp.set(2)
        self.controls.grid(column=1,row=0,sticky=W+E+N)
        self.makeControls()
        mf = self.master

    def makeControls(self):
        Radiobutton(self.controls,text='Creation Rates',
                    variable=self.comp,value=0,
                    command=self.show).grid(column=0,row=0,sticky=W)
        Radiobutton(self.controls,text='Destruction Rates',
                    variable=self.comp,value=1,
                    command=self.show).grid(column=0,row=1,sticky=W)
        Radiobutton(self.controls,text='Net Production Rates',
                    variable=self.comp,value=2,
                    command=self.show).grid(column=0,row=2,sticky=W)

    def show(self):
        mf = self.master
        mf.active = self
        c = self.comp.get()
        mix = self.top.mix
        g = mix.g
        if c == 0:
            mf.var.set("Creation Rates")
            #mf.data = spdict(mix.g, mix.moles())
            mf.comp = g.creation_rates

        elif c == 1:
            mf.var.set("Destruction Rates")
            #mf.data = spdict(mix.g,mix.mass())
            mf.comp = g.destruction_rates

        elif c == 2:
            mf.var.set("Net Production Rates")
            mf.comp = g.net_production_rates
            #mf.data = spdict(mix,mix,mf.comp)

        for s in mf.variable.keys():
            try:
                k = g.species_index(s)
                if mf.comp[k] > _CUTOFF or -mf.comp[k] > _CUTOFF:
                    mf.variable[s].set(mf.comp[k])
                else:
                    mf.variable[s].set(0.0)
            except:
                pass

class SpeciesKineticsFrame(Frame):
    def __init__(self,master,top):
        Frame.__init__(self,master)
        self.config(relief=GROOVE, bd=4)
        self.top = top
        self.top.kinetics = self
        self.g = self.top.mix.g
        self.entries=Frame(self)
        self.var = StringVar()
        self.var.set("Net Production Rates")
        self.names = self.top.mix.speciesNames()
        self.nsp = len(self.names)
        self.comp = [0.0]*self.nsp
        self.makeControls()
        self.makeEntries()
        self.entries.bind('<Double-l>',self.minimize)
        self.ctype = 0

    def makeControls(self):
        self.c = KineticsFrame(self)
        #self.rr = ReactionKineticsFrame(self, self.top)
        self.c.grid(column=1,row=0,sticky=E+W+N+S)
        #self.rr.grid(column=0,row=1,sticky=E+W+N+S)

    def show(self):
        self.c.show()

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

    def makeEntries(self):
        self.entries.grid(row=0,column=0,sticky=W+N+S+E)
        self.entries.config(relief=FLAT,bd=4)
        DATAKEYS = self.top.species
        self.variable = {}

        n=0
        ncol = 3
        col = 0
        row = 60

        for sp in DATAKEYS:
            s = sp
            k = s.index
            if row > 15:
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
                entry1.config(state=DISABLED,bg='lightgray')


class ReactionKineticsFrame(Frame):
    def __init__(self,vis,top):
        self.master = Toplevel()
        self.master.protocol("WM_DELETE_WINDOW",self.hide)
        self.vis = vis
        Frame.__init__(self,self.master)
        self.config(relief=GROOVE, bd=4)
        self.top = top
        self.g = self.top.mix.g
        nr = self.g.n_reactions
        self.eqs=Text(self,width=40,height=30)
        self.data = []
        self.start = DoubleVar()
        if nr > 30:
            self.end = self.start.get()+30
        else:
            self.end = self.start.get()+nr

        for i in range(4):
            self.data.append(Text(self,width=15,height=30))

        for n in range(nr):
            s = self.g.reaction_equation(n)
            self.eqs.insert(END,s+'\n')
        self.eqs.grid(column=0,row=1,sticky=W+E+N)
        for i in range(4):
            self.data[i].grid(column=i+1,row=1,sticky=W+E+N)
        Label(self, text='Reaction').grid(column=0,row=0,sticky=W+E+N)
        Label(self, text='Fwd ROP').grid(column=1,row=0,sticky=W+E+N)
        Label(self, text='Rev ROP').grid(column=2,row=0,sticky=W+E+N)
        Label(self, text='Net ROP').grid(column=3,row=0,sticky=W+E+N)
        Label(self, text='Kp').grid(column=4,row=0,sticky=W+E+N)

        self.scfr = Frame(self)
        self.scfr.config(relief=GROOVE,bd=4)

##              self.sc = Scrollbar(self.scfr,command=self.show,
##                              variable = self.start,
##                              orient='horizontal',length=400)
        self.sc = Scale(self.scfr,command=self.show,
                        variable=self.start,
                        orient='vertical',length=400)
#               self.sc.config(cnf={'from':0,'to':nr},variable = self.start)
        #self.sc.bind('<Any-Enter>',self.couple)
        #self.scfr.bind('<Any-Leave>',self.decouple)
        self.sc.pack(side=RIGHT,fill=Y)
        self.scfr.grid(row=0,column=6,rowspan=10,sticky=N+E+W)
        self.grid(column=0,row=0)

        self.hide()

##      def decouple(self,event=None):
##              d = DoubleVar()
##              xx = self.start.get()
##              d.set(xx)
##              self.sc.config(variable = d)

##      def couple(self,event=None):
##              self.sc.config(variable = self.start)

    def hide(self):
#               self.vis.set(0)
        self.master.withdraw()

    def show(self,e=None,b=None,c=None):
        v = self.vis.get()
        print(e,b,c)
        #if v == 0:
        #       self.hide()
        #       return

        self.master.deiconify()
        nr = self.g.n_reactions
        frop = self.g.forward_rates_of_progress
        rrop = self.g.reverse_rates_of_progress
        kp = self.g.equilibrium_constants
        self.data[0].delete(1.0,END)
        self.data[1].delete(1.0,END)
        self.data[2].delete(1.0,END)
        self.data[3].delete(1.0,END)
        self.eqs.delete(1.0,END)

        n0 = int(self.start.get())
        nn = nr - n0
        if nn > 30: nn = 30
        for n in range(n0, nn+n0):
            s = '%12.5e \n' % (frop[n],)
            self.data[0].insert(END,s)
            s = '%12.5e \n' % (rrop[n],)
            self.data[1].insert(END,s)
            s = '%12.5e \n' % (frop[n] - rrop[n],)
            self.data[2].insert(END,s)
            s = '%12.5e \n' % (kp[n],)
            self.data[3].insert(END,s)
            self.eqs.insert(END, self.g.reaction_equation(n)+'\n')

class ReactionPathFrame(Frame):

    def __init__(self,top):
        self.master = Toplevel()
        self.master.protocol("WM_DELETE_WINDOW",self.hide)
        #self.vis = vis
        Frame.__init__(self,self.master)
        self.config(relief=GROOVE, bd=4)
        self.grid(column=0,row=0)
        self.top = top
        self.g = self.top.mix.g
        self.el = IntVar()
        self.el.set(0)
        self.thresh = DoubleVar()

        scframe = Frame(self)
        self.sc = Scale(scframe,variable = self.thresh,
                        orient='horizontal',digits=3,length=300,resolution=0.01)
        self.sc.config(cnf={'from':-6,'to':0})
        Label(scframe,text='log10 Threshold').grid(column=0,row=0)
        self.sc.grid(row=0,column=1,columnspan=10)
        self.sc.bind('<ButtonRelease-1>',self.show)
        scframe.grid(row=3,column=0,columnspan=10)

        enames = self.g.element_names
        self.nel = len(enames)

        i = 1
        eframe = Frame(self)
        Label(eframe,text='Element').grid(column=0,row=0,sticky=W)
        for e in enames:
            Radiobutton(eframe,text=e,
                        variable=self.el,value=i-1,
                        command=self.show).grid(column=i,row=0,sticky=W)
            i += 1
        eframe.grid(row=0,column=0)

        self.detailed = IntVar()
        Checkbutton(self, text = 'Show details', variable=self.detailed,
                 command=self.show).grid(column=1,row=0)
        self.net = IntVar()
        Checkbutton(self, text = 'Show net flux',
                    variable=self.net,
                    command=self.show).grid(column=2,row=0)
        self.local = StringVar()
        Label(self,text='Species').grid(column=1,row=1,sticky=E)
        sp = Entry(self, textvariable=self.local,
              width=15)
        sp.grid(column=2,row=1)
        sp.bind('<Any-Leave>',self.show)



        self.fmt = StringVar()
        self.fmt.set('svg')
        i = 1
        fmtframe = Frame(self)
        fmtframe.config(relief=GROOVE, bd=4)
        self.browser = IntVar()
        self.browser.set(0)
        Checkbutton(fmtframe, text = 'Display in Web Browser',
                    variable=self.browser,
                    command=self.show).grid(column=0,columnspan=6,row=0)
        Label(fmtframe,text='Format').grid(column=0,row=1,sticky=W)
        for e in ['svg', 'png', 'gif', 'jpg']:
            Radiobutton(fmtframe,text=e,
                        variable=self.fmt,value=e,
                        command=self.show).grid(column=i,row=1,sticky=W)
            i += 1
        fmtframe.grid(row=5,column=0,columnspan=10,sticky=E+W)

        self.cv = Canvas(self,relief=SUNKEN,bd=1)
        self.cv.grid(column=0,row=4,sticky=W+E+N,columnspan=10)

        pframe = Frame(self)
        pframe.config(relief=GROOVE, bd=4)
        self.dot = StringVar()
        self.dot.set('dot -Tgif rxnpath.dot > rxnpath.gif')
        Label(pframe,text='DOT command:').grid(column=0,row=0,sticky=W)
        Entry(pframe,width=60,textvariable=self.dot).grid(column=0,
                                                 row=1,sticky=W)
        pframe.grid(row=6,column=0,columnspan=10,sticky=E+W)

        self.thresh.set(-2.0)
        self.hide()

    def hide(self):
        #self.vis.set(0)
        self.master.withdraw()

    def show(self,e=None):

        self.master.deiconify()
        el = self.g.element_name(self.el.get())
        det = False
        if self.detailed.get() == 1: det = True
        flow = 'OneWayFlow'
        if self.net.get() == 1: flow = 'NetFlow'

        self.d = ReactionPathDiagram(self.g, el)
        self.d.arrow_width = -2
        self.d.flow_type = flow
        self.d.show_details = det
        self.d.threshold = math.pow(10.0, self.thresh.get())
        node = self.local.get()
        try:
            k = self.g.species_index(node)
            self.d.display_only(k)
        except:
            self.d.display_only(-1)

        self.d.write_dot('rxnpath.dot')

        if self.browser.get() == 1:
            fmt = self.fmt.get()
            os.system('dot -T'+fmt+' rxnpath.dot > rxnpath.'+fmt)
            if fmt == 'svg': showsvg()
            elif fmt == 'png': showpng()
            else:
                path = 'file:///'+os.getcwd()+'/rxnpath.'+fmt
                webbrowser.open(path)
            try:
                self.cv.delete(self.image)
            except:
                pass
            self.cv.configure(width=0, height=0)
        else:
            os.system(self.dot.get())
            self.rp = None
            try:
                self.cv.delete(self.image)
            except:
                pass
            try:
                self.rp = PhotoImage(file='rxnpath.gif')
                self.cv.configure(width=self.rp.width(),
                                  height=self.rp.height())

                self.image = self.cv.create_image(0,0,anchor=NW,
                                                  image=self.rp)
            except:
                pass
