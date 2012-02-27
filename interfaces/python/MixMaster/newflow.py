
from Tkinter import *
from tkFileDialog import askopenfilename
import tkMessageBox

from Cantera.gases import IdealGasMix
from Cantera import *

class NewFlowDialog:

    def __init__(self, parent):

        top = self.top = Toplevel(parent)

        pl = Label(top, text='Pressure')
        pl.grid(row = 0, column = 0)

        geom = Frame(top, bd=2, relief=GROOVE)
        geom.grid(row = 1, column = 0)
        lb = Listbox(geom)
        for item in ["One-Dimensional", "Stagnation"]:
            lb.insert(END, item)
        lb.grid(row = 0, column = 0)
        glb = Listbox(geom)
        for item in ["Axisymmetric","2D"]:
            glb.insert(END, item)
        glb.grid(row = 1, column = 0)

        # ------------- pressure input ----------------

        self.p = DoubleVar()
        self.pbox = Entry(top, textvariable = self.p)
        self.pbox.grid(row = 0, column = 1)



        # ------------- gas file name input -----------


        gasf  = Frame(top, bd=2, relief=GROOVE)
        gasf.grid(row = 4, column = 0, columnspan=2)
        gl = Label(gasf, text='Gas Mixture Specification')
        gl.grid(row = 0, column = 0)


        self.infile = StringVar()
        Label(gasf, text='Mixture Input File').grid(row = 1, column = 0)
        Entry(gasf, textvariable = self.infile).grid(row = 1, column = 1)
        Button(gasf, text='Browse..', command=self.getinfile).grid(row = 1,
                                                                  column = 2)

        self.spfile = StringVar()
        Label(gasf, text='Species Database').grid(row = 2, column = 0)
        Entry(gasf, textvariable = self.spfile).grid(row = 2, column = 1)
        Button(gasf, text='Browse..', command=self.getspfile).grid(row = 2,
                                                                  column = 2)


        self.trfile = StringVar()
        Label(gasf, text='Transport Database').grid(row = 3, column = 0)
        Entry(gasf, textvariable = self.trfile).grid(row = 3, column = 1)
        Button(gasf, text='Browse..', command=self.gettrfile).grid(row = 3,
                                                                  column = 2)


        # ------------- grid -------------------------

        gf = Frame(top, bd=2, relief=GROOVE)
        gf.grid(row = 5, column = 0, columnspan=2)

        gr = Label(gf, text='Initial Grid')
        gr.grid(row = 0, column = 0)

        self.zleft = DoubleVar()
        self.zright = DoubleVar()
        ll = Label(gf, text='Left boundary at ')
        rl = Label(gf, text='Right boundary at ')
        lbb = Entry(gf, textvariable = self.zleft)
        rbb = Entry(gf, textvariable = self.zright)
        ll.grid(row = 1, column = 0)
        rl.grid(row = 2, column = 0)
        lbb.grid(row = 1, column = 1)
        rbb.grid(row = 2, column = 1)

        ok = Button(top, text = 'OK', command=self.ok)
        ok.grid(row = 20, column = 20)



    def ok(self):
        p = self.p.get()
        try:
            infile = self.infile.get()
            spfile = self.spfile.get()
            trfile = self.trfile.get()
            if spfile and trfile:
                self.gas = IdealGasMix(import_file = infile,
                                       thermo_db = spfile,
                                       transport_db = trfile)
            elif spfile:
                self.gas = IdealGasMix(import_file = infile,
                                       thermo_db = spfile)
            else:
                self.gas = IdealGasMix(import_file = infile)

        except:
            tkMessageBox.showerror('Create Gas',
                                   'Error reading file %s. See log file for more information.' % infile)

        #self.flow = Flow1D(flow_type = ftype, flow_geom = fgeom,
        #                   pressure = p, grid = gr, gas = g)
        self.top.destroy()

    def getinfile(self):
        pathname = askopenfilename(filetypes=[
            ("Input Files", "*.xml *.inp"),
            ("All Files", "*.*")])
        self.infile.set(pathname)

    def getspfile(self):
        pathname = askopenfilename(filetypes=[
            ("Species Data Files", "*.xml *.dat"),
            ("All Files", "*.*")])
        self.spfile.set(pathname)

    def gettrfile(self):
        pathname = askopenfilename(filetypes=[
            ("Transport Data Files", "*.xml *.dat"),
            ("All Files", "*.*")])
        self.trfile.set(pathname)
