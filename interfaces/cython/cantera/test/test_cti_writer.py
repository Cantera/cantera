from cantera import soln2cti
import test_mechanism_from_solution
import os
import cantera as ct

#use coln2cti to write to a cantera file format
original_solution=ct.Solution('gri30.cti')
soln2cti.write(original_solution)

#compare solution from new file to the original
new_solution=ct.Solution('pym_gri30.cti')
test_mechanism_from_solution.test(original_solution, new_solution)

os.system('rm ./pym_gri30.cti')
