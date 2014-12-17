import os, math,sys

if sys.version_info[0] == 3:
    from tkinter import *
    from tkinter.filedialog import askopenfilename
else:
    from Tkinter import *
    from tkFileDialog import askopenfilename

from cantera import *
#from Cantera.ck2ctml import ck2ctml

class ImportFrame(Frame):
    def __init__(self,top):
        self.master = Toplevel()
        self.master.title('Convert and Import CK File')
        self.master.protocol("WM_DELETE_WINDOW",self.hide)

        Frame.__init__(self,self.master)
        self.config(relief=GROOVE, bd=4)
        self.top = top
        self.infile = StringVar()

        Label(self,text="Input File").grid(row=0,column=0)
        Entry(self, width=40,
              textvariable=self.infile).grid(column=1,row=0)
        Button(self, text='Browse',
               command=self.browseForInput).grid(row=0,column=2)

        self.thermo = StringVar()
        Label(self,text="Thermodynamic Database").grid(row=1,column=0)
        Entry(self, width=40,
              textvariable=self.thermo).grid(column=1,row=1)
        Button(self, text='Browse',
               command=self.browseForThermo).grid(row=1,column=2)


        self.transport = StringVar()
        Label(self,text="Transport Database").grid(row=2,column=0)
        Entry(self, width=40,
              textvariable=self.transport).grid(column=1,row=2)
        Button(self, text='Browse',
               command=self.browseForTransport).grid(row=2,column=2)

        bframe = Frame(self)
        bframe.config(relief=GROOVE, bd=1)
        bframe.grid(row=100,column=0)
        Button(bframe, text='OK', width=8, command=self.importfile).grid(row=0,column=0)
        self.grid(column=0,row=0)
        Button(bframe, text='Cancel', width=8, command=self.hide).grid(row=0,column=1)
        self.grid(column=0,row=0)
        self.hide()

    def browseForInput(self, e=None):
        pathname = askopenfilename(
                filetypes=[("Reaction Mechanism Files",
                            ("*.inp","*.mech","*.ck2")),
                            ("All Files", "*.*")])
        if pathname:
            self.infile.set(pathname)
        self.show()

    def browseForThermo(self, e=None):
        pathname = askopenfilename(
                filetypes=[("Thermodynamic Databases",
                            ("*.dat","*.inp","*.therm")),
                           ("All Files", "*.*")])
        if pathname:
            self.thermo.set(pathname)
        self.show()

    def browseForTransport(self, e=None):
        pathname = askopenfilename(
                filetypes=[("Transport Databases", "*.dat"),
                           ("All Files", "*.*")])
        if pathname:
            self.transport.set(pathname)
        self.show()


    def importfile(self):
        ckfile = self.infile.get()
        thermdb = self.thermo.get()
        trandb = self.transport.get()
        p = os.path.normpath(os.path.dirname(ckfile))
        fname = os.path.basename(ckfile)
        ff = os.path.splitext(fname)
        nm = ""
        if len(ff) > 1: nm = ff[0]
        else: nm = ff
        outfile = p+os.sep+nm+'.xml'
        try:
            print('not supported.')
            #ck2ctml(infile = ckfile, thermo = thermdb,
            #       transport = trandb, outfile = outfile,
            #       id = nm)
            self.hide()
            return

        except:
            print('Errors were encountered. See log file ck2ctml.log')
            self.hide()
            return

        self.top.loadmech(nm,outfile,1)
        self.hide()

##              cmd = 'ck2ctml -i '+ckfile+' -o '+outfile
##              if thermdb <> "":
##                      cmd += ' -t '+thermdb
##              if trandb <> "":
##                      cmd += ' -tr '+trandb
##              cmd += ' -id '+nm
##              ok = os.system(cmd)
##              if ok == 0:
##                      self.top.loadmech(nm,outfile,1)


    def hide(self):
        #self.vis.set(0)
        self.master.withdraw()

    def show(self):
        #v = self.vis.get()
        #if v == 0:
        #       self.hide()
        #       return

        self.master.deiconify()
