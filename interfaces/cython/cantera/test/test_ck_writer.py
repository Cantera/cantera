from cantera import soln2ck
import os
import test_mechanism_from_solution
import cantera as ct
from cantera import ck2cti

#use soln2ck to write a chemkin file
original_solution = ct.Solution('gri30.cti')
soln2ck.write(original_solution)

#use ck2cti to convert that chemkin file to .cti
os.system('ck2cti --input=pym_gri30.inp')

#read in now twice-converted file and compare it to the original
new_solution = ct.Solution('pym_gri30.cti')
test_mechanism_from_solution.test(original_solution, new_solution)

os.system('rm pym_gri30.inp')
os.system('rm ./pym_gri30.cti')
