import _cantera
from exceptions import *
import os

def ck2ctml(infile = '', thermo = '-', transport = '-', outfile = '', id = ''):
    if not infile:
        raise CanteraError('No input file specified')
    fname = os.path.basename(infile)
    ff = os.path.splitext(fname)
    if len(ff) == 2:
        mechname = ff[0]
    else:
        mechname = ff
    if not outfile:
        outfile = mechname + '.xml'
    if not id:
        id = mechname
    _cantera.ck2ctml(infile, thermo, transport, outfile, id)

    

