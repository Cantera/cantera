# -*- coding: utf-8 -*-

#################################################################################################
## Copyright (C) 2017 The University Of Dayton. All Rights Reserved.
##
## No part of this program may be photocopied, transferred, or otherwise reproduced in machine
## or human readable form without the prior written consent of The University Of Dayton.
#################################################################################################
#################################################################################################
## THIS SOFTWARE IS PROVIDED "AS IS" AND ANY EXPRESSED OR IMPLIED WARRANTIES,
## INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND
## FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
## UNIVERSITY OF DAYTON OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
## INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
## LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA,
## OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
## LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
## NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE,
## EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.  THE UNIVERSITY OF DAYTON
## HAS NO OBLIGATION TO SUPPORT THE SOFTWARE.
#################################################################################################

"""
Please cite this work as follows:
Briones, A.M., Olding, R., Sykes, J.P., Rankin, B.A., McDevitt, K., Heyne, J.S., 
"Combustion Modeling Software Development, Verification and Validation," Proceedings of the 2018 ASME Power & Energy Conference,
PowerEnergy2018-7433, Lake Buena Vista, FL. 
"""
import tkinter
from tkinter import ttk
import numpy as np

# plot dialog class
class PlotDlg():

    # constructor - display plot dialog box
    def __init__(self, CVar, ZVar, cvarIndex, zvarIndex, keyList, name):
        self._dlg = tkinter.Tk()
        self._dlg.config(borderwidth=8)
        self._cvarNumIndex = tkinter.IntVar(self._dlg, value=len(CVar))
        self._zvarNumIndex = tkinter.IntVar(self._dlg, value=len(ZVar))
        self._cvarIndex = tkinter.IntVar(self._dlg, value=cvarIndex)
        self._zvarIndex = tkinter.IntVar(self._dlg, value=zvarIndex)
        self._variableName = tkinter.StringVar(self._dlg, value=name)
        self._plot = False
        self._CVar = CVar
        self._ZVar = ZVar
        self._dlg.title("Plot Selection")
        # column titles
        tkinter.Label(self._dlg, text="Min").grid( row=1, column=2)
        tkinter.Label(self._dlg, text="Max").grid( row=1, column=3)
        tkinter.Label(self._dlg, text="Indices").grid( row=1, column=4)
        tkinter.Label(self._dlg, text="Selection").grid( row=1, column=5)
        #ZVar row
        tkinter.Label(self._dlg, text="ZVar").grid( row=2, column=1, sticky="w")
        tkinter.Label(self._dlg, text=str(np.min(ZVar))).grid( row=2, column=2)
        tkinter.Label(self._dlg, text=str(np.max(ZVar))).grid( row=2, column=3)
        tkinter.Label(self._dlg, textvariable=self._zvarNumIndex).grid( row=2, column=4)
        tkinter.Entry(self._dlg, textvariable=self._zvarIndex, width=6).grid(row=2, column=5, pady=2, padx=4)
        # CVar row
        tkinter.Label(self._dlg, text="CVar").grid( row=3, column=1, sticky="w")
        tkinter.Label(self._dlg, text=str(np.min(CVar))).grid( row=3, column=2)
        tkinter.Label(self._dlg, text=str(np.max(CVar))).grid( row=3, column=3)
        tkinter.Label(self._dlg, textvariable=self._cvarNumIndex).grid( row=3, column=4)
        tkinter.Entry(self._dlg, textvariable=self._cvarIndex, width=6).grid(row=3, column=5, pady=2, padx=4)
        #variable selection
        tkinter.Label(self._dlg, text="Variable").grid( row=4, column=1, sticky="w")
        self._variable_combobox = ttk.Combobox(self._dlg, textvariable=self._variableName, values = keyList, width=20).grid(row=4,column=2,columnspan=3, pady=4) 
        # button row
        tkinter.Button(self._dlg, text = "Plot", command = self.RetryCallBack, width = 6).grid( row=5, column=3, pady=6)
        tkinter.Button(self._dlg, text = "Exit", command = self.ExitCallBack, width = 6).grid( row=5, column=5, pady=6)
        
        # center Dialog on screen
        self._dlg.withdraw()
        self._dlg.update_idletasks()  # Update "requested size" from geometry manager

        x = (self._dlg.winfo_screenwidth() - self._dlg.winfo_reqwidth()) / 2
        y = (self._dlg.winfo_screenheight() - self._dlg.winfo_reqheight()) / 2
        self._dlg.geometry("+%d+%d" % (x, y))

        # This seems to draw the window frame immediately, so only call deiconify()
        # after setting correct window position
        self._dlg.deiconify()

    @property
    def cvar_index(self):
        '''cvar index'''
        return self._cvarIndex.get()

    @property
    def zvar_index(self):
        '''zvar index'''
        return self._zvarIndex.get()

    @property
    def variable_name(self):
        '''name of variable selected'''
        return self._variableName.get()
    
    # wait for user to hit button
    def WaitResult(self):
        self._dlg.mainloop()
        return self._plot

    # retry button handler
    def RetryCallBack(self):
        self._plot = True
        if self._cvarIndex.get() < 0:
            self._cvarIndex.set(0)
        if self._cvarIndex.get() >= len(self._CVar):
            self._cvarIndex.set(len(self._CVar)-1)
        if self._zvarIndex.get() < 0:
            self._zvarIndex.set(0)
        if self._zvarIndex.get() >= len(self._ZVar):
            self._zvarIndex.set(len(self._ZVar)-1)
        self._dlg.quit()
        self._dlg.destroy()

    # exit buttonhandler
    def ExitCallBack(self):
        self._plot = False
        self._dlg.quit()
        self._dlg.destroy()
