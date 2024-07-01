#%%
"""
A freely-propagating, premixed hydrogen flat flame with multicomponent
transport properties.

Requires: cantera >= 3.0
Keywords: combustion, 1D flow, premixed flame, multicomponent transport,
          saving output
"""
import sys, os
sys.path.append("C:/Users/pjsin/Documents/cantera/build/python")
from pathlib import Path
import cantera as ct
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
df = pd.read_csv('black.csv')
p_list = np.linspace(0,20,50)[1:]

# Simulation parameters
# p = ct.one_atm  # pressure [Pa]
Tin = 300.0  # unburned gas temperature [K]
reactants = 'H2:0.1071, O2:0.1785, He:0.7144'  # premixed gas composition
width = 0.03  # m
loglevel = 0  # amount of diagnostic output (0 to 8)

# Solution object used to compute mixture properties, set to the state of the
# upstream fuel-air mixture



plt.figure()
mbr = []
for p in p_list:

    gas = ct.Solution('test/data/alzuetamechanism.yaml')
    gas.TPX = Tin, p*ct.one_atm, reactants
    
    f = ct.FreeFlame(gas, width=width)
    f.set_refine_criteria(ratio=3, slope=0.06, curve=0.12)
    
    
    f.transport_model = 'mixture-averaged'
    f.solve(loglevel=loglevel, auto=True)
    
    mbr.append(f.velocity[0]*f.density[0] / 10) # g/cm2*s
    
plt.plot(p_list,mbr,label='simulation')
plt.scatter(df['Pressure'], df['Mass burning rate'], label='Burke')
plt.legend()
plt.xlabel('Pressure [atm]')
plt.ylabel('Mass Burning Rate [g/cm2*s]')
plt.show()    
    
    
    # print(f"mixture-averaged flamespeed = {f.velocity[0]:7f} m/s")
    
    
    
    # f.transport_model = 'multicomponent'
    # f.solve(loglevel)  # don't use 'auto' on subsequent solves
    # print(f"multicomponent flamespeed = {f.velocity[0]:7f} m/s")

# %%
