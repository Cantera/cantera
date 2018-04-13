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

import os
import tkinter
import sys

###############################################################################
#
#  Defines class used to display the 'Hit Ok to continue' dialog
#
###############################################################################

class MessageDlg():

    # constructor - display retry dialog box
    def __init__(self, message):
        self._dlg = tkinter.Tk()
        self._dlg.config(borderwidth=8)
        self._dlg.title("Message Box")

        # property definitions
        self._message = tkinter.StringVar(self._dlg, value=message)

        # message text
        tkinter.Label(self._dlg, textvariable=self._message, padx=2).grid( row=1, column=1)
        tkinter.Label(self._dlg, text="Please hit OK to continue").grid(row=2, column=1)

        # retry and exit buttons
        tkinter.Button(self._dlg, text = "Ok", command = self.OkCallBack, width = 6).grid( row=3, column=1, pady=8)
    
    # wait for user to hit button
    def WaitResult(self):
        self._dlg.mainloop()

    # OK button handler
    def OkCallBack(self):
        self._dlg.quit()
        self._dlg.destroy()


# function that flushes stdout after each print
def PrintFlush(*arg):
    print(*arg)
    sys.stdout.flush()
