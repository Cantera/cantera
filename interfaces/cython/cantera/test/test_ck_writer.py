from cantera import soln2ck
import os
import cantera as ct
from cantera import ck2cti
A=ct.Solution('gri30.cti')
soln2ck.write(A)
os.system('ck2cti --input=pym_gri30.inp')
new_solution = ct.Solution('pym_gri30.cti')
os.system('rm pym_gri30.inp')
os.system('rm ./pym_gri30.cti')
