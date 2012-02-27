import string
import os, sys
import types, traceback

try:
    from Tkinter import Tk
    import tkMessageBox
    _hasTk = 1
except:
    _hasTk = 0


def write_CSV(f,x):
    """write list x to file f in comma-separated-value format."""
    for e in x:
        f.write( `e`+',')
    f.write('\n')


def _print_value(name, value, unitstr):
    print string.rjust(name, 15)+  \
          string.rjust('%10.5e' %value, 15) + ' ' + \
          string.ljust(unitstr,5),

def hasTk():
    try:
        import tkMessageBox
        return 1
    except:
        return 0

def handleError(message = '<error>', window = None,
                fatal = 0, warning = 0, options = None):
    if warning:
        tkMessageBox.showwarning(title = 'Warning', message = message,
                                 parent = window)
    else:
        m = tkMessageBox.showerror(title = 'Error', message = message,
                                   parent = window)
