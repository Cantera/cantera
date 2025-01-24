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
f = ct.BurnerFlame(gas=gas, width=z[-1], sections=25)

# Setting up soot computation
f.soot_setup(precursors        =['A2'],
             #fractal_aggregates=True,
             trash_section     =-1
              )

# Restoring the sooting flame
if ct.__version__ >= '2.5.0':
    f.restore('./RESULTS/%s.yaml' % mech, 'sooting')
else:
    f.restore('./RESULTS/%s.xml' % mech, 'sooting')

# Plot volume fraction
# Particles under 2nm and in last section are excluded
plt.figure(0)
plt.plot(f.grid * 100, f.soot_fv(min=2e-7, last=-1))
plt.yscale('log')
plt.ylabel('Volume fraction')
plt.xlabel('HAB ($cm$)')

# Plot particles number density
# Particles under 2nm and in last section are excluded
plt.figure(1)
plt.plot(f.grid * 100, f.soot_Np(min=2e-7, last=-1))
plt.yscale('log')
plt.ylabel('Soot particles number density ($cm^{-3}$)')
plt.xlabel('HAB ($cm$)')

# Displaying graphs
plt.show()
