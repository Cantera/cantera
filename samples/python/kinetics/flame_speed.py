"""
Simulate flame speeds across a range of equivalence ratios, and observe the impact of incorporating the
reduced-pressure linear mixture rule (LMR-R) in such calculations.

Here we will consider a mixture of NH3/air at 1 atm and 296 K, and compare results against the
experimental measurements of Ronney [1]

Two models are compared in this example:
(i)  A 2023 model of H2 and NH3 chemistry published by Alzueta et al. [2]
(ii) An adapted version of this model that has applied the reduced-pressure linear mixture
     rule (LMR-R) and ab initio third-body efficiencies [3]

References:
[1] P. D. Ronney, Effect of chemistry and transport properties on near-limit flames at microgravity,
    Combust. Sci. Tech. 59 (1988) 123 -- 141.
[2] M. U. Alzueta, I. Salas, H. Hashemi, P. Glarborg, CO-assisted NH3 oxidation, Combust. Flame 257
    (2023) 112438.
[3] P. J. Singal, J. Lee, L. Lei, R. L. Speth, M. P. Burke, Implementation of New Mixture Rules
    Has a Substantial Impact on Combustion Predictions for H2 and NH3, Proc. Combust. Inst. 40
    (2024).

Requires: cantera >= 3.1
Keywords: burning velocity, flame speed, equivalence ratio, mixture rule, LMR-R
"""

import sys, os
sys.path.append("C:/Users/pjsin/Documents/cantera/build/python")
import cantera as ct
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

plt.figure(figsize=(3.5,2.5))

models = {    
          'Original':'alzueta.yaml',
          'LMR-R':'alzueta_LMRR.yaml',
          }
colours = ["xkcd:grey",'xkcd:purple']

fuel_list = np.linspace(0.14,0.4,20) # mole fractions of fuel to simulate across
alpha = 1.0
a_st = 0.75
Tin = 296  # unburned gas temperature [K]
p=760

for k, m in enumerate(models):
  vel_list = []
  phi_list = []
  for j, fuel_frac in enumerate(fuel_list):
      gas = ct.Solution(list(models.values())[k])
      NH3 = alpha*fuel_frac
      H2 = (1-alpha)*fuel_frac
      ox_frac = 1 - fuel_frac # oxidizer fraction
      O2 = ox_frac*0.21
      N2 = ox_frac*0.79
      phi = np.divide(fuel_frac/O2,1/a_st)
      phi_list.append(phi)
      X = {'NH3':NH3,'H2':H2,'O2':O2,'N2':N2}
      gas.TPX = Tin, (p/760)*ct.one_atm, X
      f = ct.FreeFlame(gas, width=0.03)
      f.set_refine_criteria(ratio=3, slope=0.05, curve=0.05)
      f.transport_model = 'multicomponent'
      f.soret_enabled = True
      f.solve(loglevel=1, auto=True)
      vel_list.append(f.velocity[0] * 100) # cm/s
  plt.plot(phi_list, vel_list, color=colours[k],label=m)

expData = {
   'X_NH3': [16.3,16.4,17.0,18.0,19.0,20.0,21.9,24.0,26.0,28.5,29.0,30.0,31.0,31.5],
   'vel': [1.35,1.48,2.30,3.36,4.01,5.88,6.80,8.14,6.73,5.00,4.78,3.3,2.9,3.0]
}
X_NH3 = np.divide(expData['X_NH3'],100)
X_O2 = np.multiply(np.subtract(1,X_NH3), 0.21)
phi_data = np.divide(np.divide(X_NH3,X_O2),np.divide(4,3))
plt.plot(phi_data,expData['vel'],'o',fillstyle='none',color='k',label='Ronney')
plt.legend(fontsize=8,frameon=False, loc='lower right')
plt.ylabel(r'Burning velocity [cm $\rm s^{-1}$]')
plt.xlabel(r'Equivalence Ratio')
plt.xlim([0.6, 2.1])
plt.ylim([0, 12])
plt.show()