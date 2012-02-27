

# functionality imports
from Tkinter import *

from Cantera.gui import menu, newflow

class App:
    def __init__(self, master):

        try:
            self.root = master.root
        except:
            self.root = master

        self.frame = Frame(master)
        self.frame.grid(row = 0, column = 0)

        self.makemenu(self.frame)

        self.quitbutton = Button(self.frame, text = "Quit",
                                 command = self.frame.quit)
        self.quitbutton.grid(row = 1, column = 0)

        self.newbutton = Button(self.frame, text = "New...",
                                command = self.notyet)
        self.newbutton.grid(row = 1, column = 1)

    def notyet(self):
        print 'not yet!'


    def newflow(self):
        n = newflow.NewFlowDialog(self.root)

    def makemenu(self,frame):
        self.menubar = Frame(frame, relief=FLAT, bd=0)
        self.menubar.grid(row = 0, column = 0)

        self.filemenu = menu.make_menu('File', self.menubar,
                                  [('New...', self.newflow),
                                   ('Open...', self.notyet),
                                   ('Save As...', self.notyet),
                                   'separator',
                                   ('Exit', frame.quit),
                                   []
                                   ])



root = Tk()

app = App(root)

root.mainloop()
