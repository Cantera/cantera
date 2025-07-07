# Importing necessary libraries
import cantera as ct
import numpy as np
import os
os.chdir(os.path.dirname(os.path.abspath(__file__)))

# Setting the mechanism name
mech = 'BISETTI'

#Setting pressure
pressure = ct.one_atm

# Setting temperature
temperature = 298.00

# Setting inlet velocity
velocity = 6.73e-2

# Setting flame initial composition
composition = 'C2H4:14.08, O2:18.05, N2:67.87'

# Creating gas object
if ct.__version__ >= '2.5.0':
    gas = ct.Solution('./inputs/%s.yaml' % mech)
else:
    gas = ct.Solution('./inputs/%s.cti' % mech)

# Settng gas properties
gas.TPX = temperature, pressure, composition

# Loading experimental data (temperature vs. height)
z, T = np.genfromtxt('./inputs/T_vs_x.dat').T

# Creating flame object of the same size as the experimental flame
f = ct.BurnerFlame(gas=gas, width=z[-1])

# Setting the flame's mass flow rate
f.burner.mdot = velocity * gas.density

# Disabling energy equation
f.energy_enabled = False

#Imposing proposed temperature profile
f.flame.set_fixed_temp_profile(z/max(z), T)

# Setting the transport model to mixture-averaged
f.transport_model = 'mixture-averaged'

# Setting mesh refinement criteria
f.set_refine_criteria(ratio=2.0, slope=0.1, curve=0.1)

# Solving the problem
f.solve(1, 'refine')

# Showing the result
f.show()

# Saving the result
if ct.__version__ >= '2.5.0':
    f.save('./RESULTS/%s.yaml' % mech, 'non-sooting',overwrite=True)
else:
    f.save('./RESULTS/%s.xml' % mech, 'non-sooting')
