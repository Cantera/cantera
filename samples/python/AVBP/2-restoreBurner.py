# Importing necessary libraries
import cantera as ct
import numpy as np
import matplotlib.pyplot as plt
import os
os.chdir(os.path.dirname(os.path.abspath(__file__)))

# Setting the mechanism name
mech = 'BISETTI'

# Creating gas object
#if cantera version is 2.5.0 or higher, use the following line instead
if ct.__version__ >= '2.5.0':
    gas = ct.Solution('./inputs/%s.yaml' % mech)
else:
    gas = ct.Solution('./inputs/%s.cti' % mech)
    
# Loading experimental data (temperature vs. height)
z, T = np.genfromtxt('./inputs/T_vs_x.dat').T

# Creating flame object of the same size as the experimental flame with 25 soot sections
f = ct.BurnerFlame(gas=gas, width=z[-1], sections = 25)

# Setting transport model
f.transport_model = 'Mix'

# Setting up soot calculation
f.soot_setup(precursors        =['A2'],
             retroaction       =True,
             condensation      =True,
             coagulation       =True,
             surface_growth    =True,
             oxidation         =True,
             #fractal_aggregates=True,
             radiation         =True,
             trash_section     =-1,
             #show_sections     = True
             )

# Restoring the non-sooting flame
if ct.__version__ >= '2.5.0':
    f.restore('./RESULTS/%s.yaml' % mech, 'non-sooting')
else:
    f.restore('./RESULTS/%s.xml' % mech, 'non-sooting')
f.radiation_enabled = True

# Solving the sooting flame
f.solve(1,refine_grid='disabled')

# # Saving the result
if ct.__version__ >= '2.5.0':
    f.save('./RESULTS/%s.yaml' % mech,'sooting',overwrite=True)
else:
    f.save('./RESULTS/%s.xml' % mech,'sooting')