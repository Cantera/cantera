"""This module defines classes and functions used to model gas mixtures."""

from constants import *
from solution import Solution
from ck2ctml import ck2ctml
#import _cantera
import os

def IdealGasMix(src="", root=None, transport='None',
                thermo = "", trandb = ""):
    p = os.path.normpath(os.path.dirname(src))        
    fname = os.path.basename(src)
    ff = os.path.splitext(fname)
    nm = ""
    if len(ff) > 1:
        nm = ff[0]
        ext = ff[1]
    else:
        nm = ff
        ext = ''
    if ext <> '.xml' and ext <> '.XML' and ext <> '.ctml' and ext <> '.CTML':
        outfile = p+os.sep+nm+'.xml'
        ck2ctml(infile = src, outfile = outfile, thermo = thermo,
                transport = trandb, id = nm)
        return Solution(src=outfile, root=None, transport=transport)
    else:
        return Solution(src=src, root=root, transport=transport)

def GRI30(transport='None'):
    return Solution(src="gri30.xml#gri30_hw", transport=transport)

def Air():
    return Solution(src="air.xml#air")

def Argon():
    return Solution(src="argon.xml#argon")

## def H_O_AR(transport=None, chem = 1):
##     """
##     The hydrogen/oxygen/argon portion of GRI-Mech 3.0.

##     Parameters:
##     transport --- transport model (None, 'Mix,' or 'Multi'). Default: None
##     chem      --- chemistry disabled if chem = 0. Default: enabled.
##     """
##     if chem == 1:
##         return Solution(import_file='h2o2.inp', thermo_db="",
##                         eos=1, kmodel=1, trmodel=transport,
##                         transport_db='gri30_tran.dat',
##                         validate=0)
##     else:
##         return Solution(import_file='h2o2_noch.inp', thermo_db="", 
##                         eos=1, kmodel=1, trmodel=transport,
##                         transport_db='gri30_tran.dat',
##                         validate=0)        


## def GRI30(transport=None):
##     """
##     GRI-Mech 3.0.

##     Parameters:
##     transport --- transport model (None, 'Mix,' or 'Multi'). Default: None
##     """
##     return Solution(import_file='gri30.xml', thermo_db="",
##                     eos=1, kmodel=2, trmodel=transport, id="gri30",
##                     transport_db='gri30_tran.dat',
##                     validate=0)


## def Air(transport=None, T=300.0, P=OneAtm):
##     """
##     Edited version of GRI-Mech 3.0 containing only O/N/AR species. Initial
##     composition is set to 21% O2, 78% N2, 1% Ar.
    
##     Parameters:
##     transport --- transport model (None, 'Mix,' or 'Multi'). Default: None
##     T         --- temperature
##     P         --- pressure
##     """
##     gas = Solution(import_file='air.inp', thermo_db="nasathermo.dat",
##                    eos=1, kmodel=1, trmodel=transport, id="air",
##                    transport_db='gri30_tran.dat',
##                    validate=0)
##     gas.setState_TPX(T, P, 'O2:0.21, N2:0.78, Ar:0.01')
##     return gas

    
