"""
This example generates a saturated steam table and plots the vapor dome. The
steam table corresponds to data typically found in thermodynamic text books
and uses the same customary units.

Requires: Cantera >= 2.5.0, matplotlib >= 2.0, pandas >= 1.1.0, numpy >= 1.12
"""

import cantera as ct
import pandas as pd
import numpy as np
from matplotlib import pyplot as plt

w = ct.Water()

# create colums
columns = ['T', 'P',
           'vf', 'vfg', 'vg',
           'uf', 'ufg', 'ug',
           'hf', 'hfg', 'hg',
           'sf', 'sfg', 'sg']

# temperatures correspond to Engineering Thermodynamics, Moran et al. (9th ed),
# Table A-2; additional data points are added close to the critical point;
# w.min_temp is equal to the triple point temperature
degc = np.hstack([np.array([w.min_temp - 273.15, 4, 5, 6, 8]),
                  np.arange(10, 37), np.array([38]),
                  np.arange(40, 100, 5), np.arange(100, 300, 10),
                  np.arange(300, 380, 20), np.arange(370, 374),
                  np.array([w.critical_temperature - 273.15])])

df = pd.DataFrame(0, index=np.arange(len(degc)), columns=columns)
df.T = degc

arr = ct.SolutionArray(w, len(degc))

# saturated vapor data
arr.TQ = degc + 273.15, 1
df.P = arr.P_sat / 1.e5
df.vg = arr.v
df.ug = arr.int_energy_mass / 1.e3
df.hg = arr.enthalpy_mass / 1.e3
df.sg = arr.entropy_mass / 1.e3

# saturated liquid data
arr.TQ = degc + 273.15, 0
df.vf = arr.v
df.uf = arr.int_energy_mass / 1.e3
df.hf = arr.enthalpy_mass / 1.e3
df.sf = arr.entropy_mass / 1.e3

# delta values
df.vfg = df.vg - df.vf
df.ufg = df.ug - df.uf
df.hfg = df.hg - df.hf
df.sfg = df.sg - df.sf

# reference state (triple point; liquid state)
w.TQ = w.min_temp, 0
uf0 = w.int_energy_mass / 1.e3
hf0 = w.enthalpy_mass / 1.e3
sf0 = w.entropy_mass / 1.e3
pv0 = w.P * w.v / 1.e3

# change reference state
df.ug -= uf0
df.uf -= uf0
df.hg -= hf0 - pv0
df.hf -= hf0 - pv0
df.sg -= sf0
df.sf -= sf0

# print and write saturated steam table to csv file
print(df)
df.to_csv('saturated_steam_T.csv', index=False)

# illustrate the vapor dome in a P-v diagram
plt.semilogx(df.vf.values, df.P.values, label='Saturated liquid')
plt.semilogx(df.vg.values, df.P.values, label='Saturated vapor')
plt.semilogx(df.vg.values[-1], df.P.values[-1], 'o', label='Critical point')
plt.xlabel(r'Specific volume - $v$ ($\mathrm{m^3/kg}$)')
plt.ylabel(r'Presssure - $P$ (bar)')
plt.legend()

# illustrate the vapor dome in a T-s diagram
plt.figure()
plt.plot(df.sf.values, df['T'].values, label='Saturated liquid')
plt.plot(df.sg.values, df['T'].values, label='Saturated vapor')
plt.plot(df.sg.values[-1], df['T'].values[-1], 'o', label='Critical point')
plt.xlabel(r'Specific entropy - $s$ ($\mathrm{kJ/kg-K}$)')
plt.ylabel(r'Temperature - $T$ (${}^\circ C$)')
plt.legend()

plt.show()
