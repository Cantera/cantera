"""Gas mixtures.

These functions all return instances of class Solution that represent
gas mixtures.

"""
# for pydoc
import solution, constants, ck2ctml

from constants import *
from Cantera.solution import Solution
from ck2ctml import ck2ctml
#import _cantera
import os

def IdealGasMix(src="", root=None, transport='None',
                thermo = "", trandb = ""):
    """Return a Solution object representing an ideal gas mixture.

    src       --- input file
    root      --- root of an XML tree containing the phase specification.
                  Specify src or root but not both.
    thermo    --- auxiliary thermo database
    transport --- transport model
    trandb    --- transport database
    """
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
##     if ext <> '.xml' and ext <> '.XML' and ext <> '.ctml' and ext <> '.CTML':
##         outfile = p+os.sep+nm+'.xml'
##         if ext == '.py':
##             from Cantera import pip
##             pip.process(fname)
##         else:
##             ck2ctml(infile = src, outfile = outfile, thermo = thermo,
##                     transport = trandb, id = nm)
##         return Solution(src=outfile, root=None, transport=transport)
##     else:
    return Solution(src=src, root=root, transport=transport)

def GRI30(transport='None'):
    """Return a Solution instance implementing reaction mechanism
    GRI-Mech 3.0."""
    return Solution(src="gri30.xml#gri30_hw", transport=transport)

def Air():
    """Return a Solution instance implementing the O/N/Ar portion of
    reaction mechanism GRI-Mech 3.0. The initial composition is set to
    that of air"""    
    return Solution(src="air.xml#air")

def Argon():
    """Return a Solution instance representing pure argon."""    
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

    
