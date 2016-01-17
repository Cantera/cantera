import sys
if sys.version_info[0] == 3:
    from tkinter import *
else:
    from Tkinter import *

class TransportFrame(Frame):
    def show(self, i, frame, row, col):
        if self.checked[i].get():
            frame.grid(row=row,column=col,sticky=N+E+S+W)
        else:
            frame.grid_forget()

    def showcomp(self):
        self.show(0, self.top.mixfr, 8, 0)

    def showthermo(self):
        self.show(1, self.top.thermo, 7, 0)

    def __init__(self,master,top):
        self.top = top
        self.c = []
        self.checked = []
        Frame.__init__(self,master)
        self.config(relief=GROOVE, bd=4)
        lbl = ['multicomponent', 'mixture-averaged']
        cmds = [self.showcomp, self.showthermo]
        for i in range(2):
            self.checked.append(IntVar())
            self.checked[i].set(0)
            self.c.append(Checkbutton(self,
                               text=lbl[i],
                               variable=self.checked[i],
                               onvalue=1,
                               offvalue=0,
                               command=cmds[i]
                               ))
            self.c[i].grid(column=i,row=0, sticky=W+N)
