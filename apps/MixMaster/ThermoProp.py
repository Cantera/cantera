from Tkinter import *
from UnitChooser import UnitVar

class ThermoProp:
	def __init__(self, master, thermoframe, row, name, value, units, defaultunit=0):
		self.value = DoubleVar()
                self.thermoframe = thermoframe
		self.entry = UnitVar(master,units,defaultunit)
		self.entry.grid(column=1,row=row,sticky=W)
		self.entry.v.config(state=DISABLED,bg='lightgray')			
		self.checked=IntVar()
		self.checked.set(0)
		self.c=Checkbutton(master, 
			      text=name,
			      variable=self.checked, 
			      onvalue=1,
			      offvalue=0,
			      command=self.check
			      )
		self.c.grid(column=0,row=row, sticky=W+N)

	def check(self):
            self._check()
            self.thermoframe.last2.checked.set(0)
            self.thermoframe.last2._check()            
            self.thermoframe.last2 = self.thermoframe.last1
            self.thermoframe.last1 = self
                
	def _check(self):
		if self.isChecked():
			self.entry.v.config(state=NORMAL,bg='white')
		else:
			self.entry.v.config(state=DISABLED,bg='lightgray') 

	def isChecked(self):
		return self.checked.get()

	def set(self, value):
		self.entry.set(value)
                
	def get(self):
		return self.entry.get()
			

