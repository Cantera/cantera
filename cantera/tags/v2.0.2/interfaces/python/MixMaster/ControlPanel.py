from types import *
from Tkinter import *
from ScrolledText import ScrolledText
#import datawindow
#import filewindow

def ff():
    print ' hi '

class ControlWindow(Frame):
    fncs = [ff]*10

    def __init__(self, title, master=None):
        self.app = master
        Frame.__init__(self,master)
        self.grid(row=0,column=0,sticky=E+W+N+S)
        self.master.title(title)


    def addButtons(self, label, funcs):
        self.buttonholder = Frame(self, relief=FLAT, bd=2)
        self.buttonholder.pack(side=TOP,anchor=W)
        b = Label(self.buttonholder,text=label)
        b.pack(side=LEFT,fill=X)
        for f in funcs:
            b=Button(self.buttonholder,
                     text=f[0],command=f[1], padx=1,pady=1)
            b.pack(side=LEFT,fill=X)

    def disableButtons(self, *buttons):
        for button in self.buttonholder.slaves():
            if (button.cget('text') in buttons):
                try:
                    button.config(state=DISABLED)
                except:
                    pass

    def enableButtons(self, *buttons):
        for button in self.buttonholder.slaves():
            if (button.cget('text') in buttons):
                try:
                    button.config(state=NORMAL)
                except:
                    pass



    def newFrame(self, label, var):
        fr = Frame(self, relief = RIDGE, bd = 2)
        fr.pack(side=TOP,fill=X)
        c = Checkbutton(fr, variable=var)
        c.pack(side = LEFT, fill = X)
        b = Label(fr,text=label,foreground="NavyBlue")
        b.pack(side=LEFT,fill=X)
        return fr

    ##creates a new Toplevel object
    ##options:  transient=<callback for window close>,
    ##                      placement=(<screen x-coord>, <screen y-coord>)
    def newWindow(self, master, title, **options):
        new = Toplevel(master)
        new.title(title)
        #new.config(takefocus=0)
        if 'transient' in options.keys():
            new.transient(master)
            if options['transient']:
                new.protocol('WM_DELETE_WINDOW', options['transient'])
        if 'placement' in options.keys():
            new.geometry("+%d+%d" % tuple(options['placement']))
        return new

    ##routes mouse and keyboard events to the window and
    ##waits for it to close before returning
    def makemodal(self, window):
        window.focus_set()
        window.grab_set()
        window.wait_window()
        return


    def PlotMenu(self, fr, label, funcs):
        filebutton = Menubutton(fr,text=label, padx=3,pady=1)
        filebutton.pack(side=LEFT)
        filemenu = Menu(filebutton,tearoff=TRUE)
        i = 0
        for f in funcs:
            filemenu.add_command(label=f[0], command=f[1])
            i = i + 1
        filebutton['menu']=filemenu
        return filemenu

def testevent(event):
    print 'event ',event.value

def make_menu(name, menubar, list):
    nc = len(name)
    button=Menubutton(menubar, text=name, width=nc+4, padx=3,pady=1)
    button.pack(side=LEFT)
    menu = Menu(button,tearoff=FALSE)
    m = menu
    i = 0
    for entry in list:
        i += 1
        if entry == 'separator':
            menu.add_separator({})
        elif type(entry)==ListType:
            for num in entry:
                menu.entryconfig(num,state=DISABLED)
        elif type(entry[1]) != ListType:
            if i == 20:
                i = 0
                submenu = Menu(button,tearoff=FALSE)
                m.add_cascade(label='More...',
                                 menu=submenu)
                m = submenu
            if len(entry) == 2 or entry[2] == 'command':
                m.add_command(label=entry[0],
                              command=entry[1])
            elif entry[2] == 'check':
                entry[3].set(0)
                if len(entry) >= 5: val = entry[4]
                else: val = 1
                m.add_checkbutton(label=entry[0],
                                     command=entry[1],
                                     variable = entry[3],
                                     onvalue=val)
        else:
            submenu=make_menu(entry[0], menu, entry[1])
            m.add_cascade(label=entry[0],
                             menu=submenu)
    button['menu']=menu
    return button

def menuitem_state(button, *statelist):
    for menu in button.children.keys():
        if isinstance(button.children[menu], Menu):
            for (commandnum, onoff) in statelist:
                if onoff==0:
                    button.children[menu].entryconfig(commandnum,state=DISABLED)
                if onoff==1:
                    button.children[menu].entryconfig(commandnum,state=NORMAL)
        else:
            pass

class ArgumentWindow(Toplevel):
    import tkMessageBox
    def __init__(self, sim, **options):
        Toplevel.__init__(self, sim.cwin)
        self.resizable(FALSE,FALSE)
        self.protocol("WM_DELETE_WINDOW", lambda:0)  #self.cancelled)
        self.transient(sim.cwin)
        if 'placement' in options.keys():
            self.geometry("+%d+%d" % tuple(options['placement']))
        self.title('Thermal Model Initialization')
        self.sim = sim

        self.make_options()

        buttonframe = Frame(self)
        buttonframe.pack(side=BOTTOM)
        b1=Button(buttonframe, text='OK', command=self.callback)
        b1.pack(side=LEFT)
        #b2=Button(buttonframe, text='Cancel', command=self.cancelled)
        #b2.pack(side=LEFT)
        self.bind("<Return>", self.callback)
        #self.bind("<Escape>", self.cancelled)

        self.initial_focus = self
        self.initial_focus.focus_set()
        #self.wait_window(self)

    def make_options(self):
        pass

        ###  must override this function    ###
        ###  with the entry forms           ###
        ###  be sure to use pack or a       ###
        ###      frame that is packed into self ###

    def getArguments(self):
        pass

        ###  must override this function  ###
        ###  with the validation checking ###
        ###  must return None if error,   ###
        ###  and non_null if ok                   ###


    def callback(self, event=None):
        g=self.getArguments()
        if not g:
            self.initial_focus.focus_set()
            return
        self.withdraw()
        self.update_idletasks()

        self.assign(g)
        self.cancelled()

    def assign(self, obj):
        pass

        ###  must override this function  ###
        ###  to do the assignment in sim  ###

    def cancelled(self,event=None):
        self.sim.cwin.focus_set()
        self.destroy()


if __name__=='__main__':
    t = Tk()
    ControlWindow(t).mainloop()
