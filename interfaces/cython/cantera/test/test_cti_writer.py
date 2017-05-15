from cantera import soln2cti
import os
import cantera as ct
A=ct.Solution('gri30.cti')
soln2cti.write(A)
solution=ct.Solution('pym_gri30.cti')
os.system('rm ./pym_gri30.cti')
