"""
Flame speed as a function of equivalence ratio
==============================================

Simulate flame speeds across a range of equivalence ratios, and observe the impact of
incorporating the reduced-pressure linear mixture rule (LMR-R) in such calculations.

Here we will consider a mixture of NH3/air (1 atm, 296 K) and compare results against
the experimental measurements of Ronney. [1] Two models are compared in this example:

1. A 2023 model of H2 and NH3 chemistry published by Alzueta et al. [2]
2. An adapted version of this model that has applied the reduced-pressure linear mixture
   rule (LMR-R) and ab initio third-body efficiencies. [3]

References:

    [1] P. D. Ronney, Effect of chemistry and transport properties on near-limit flames
    at microgravity, Combust. Sci. Tech. 59 (1988) 123 -- 141.

    [2] M. U. Alzueta, I. Salas, H. Hashemi, P. Glarborg, CO-assisted NH3 oxidation,
    Combust. Flame 257 (2023) 112438.

    [3] P. J. Singal, J. Lee, L. Lei, R. L. Speth, M. P. Burke, Implementation of New
    Mixture Rules Has a Substantial Impact on Combustion Predictions for H2 and NH3,
    Proc. Combust. Inst. 40 (2024).

Requires: cantera >= 3.1

.. tags:: flame speed, kinetics
"""

import cantera as ct
import matplotlib.pyplot as plt
import numpy as np

fig, ax = plt.subplots()
file = 'example_data/ammonia-CO-H2-Alzueta-2023.yaml'
models = {'Original': 'baseline', 'LMR-R': 'linear-Burke'}
colours = ["xkcd:grey",'xkcd:purple']
Tin = 296  # unburned gas temperature [K]
p=760  # pressure [torr]
n=16 # number of points to simulate
phi_list = np.linspace(0.6,2.0,n) # equivalence ratios to simulate across
for k, m in enumerate(models):
    vel_list = []
    gas = ct.Solution(file, name=models[m])
    for j, phi in enumerate(phi_list):
        gas.set_equivalence_ratio(phi, 'NH3', {'O2':1, 'N2': 3.76})
        gas.TP = Tin, (p/760)*ct.one_atm
        f = ct.FreeFlame(gas, width=0.03)
        f.set_refine_criteria(ratio=3, slope=0.06, curve=0.10)
        # f.transport_model = 'multicomponent' # optionally enable
        # f.soret_enabled = True  # optionally enable
        f.solve(loglevel=1, auto=True)
        vel_list.append(f.velocity[0] * 100) # cm/s
    ax.plot(phi_list, vel_list, color=colours[k], label=m)
expData = {
   'X_NH3': [16.3,16.4,17.0,18.0,19.0,20.0,21.9,24.0,26.0,28.5,29.0,30.0,31.0,31.5],
   'vel': [1.35,1.48,2.30,3.36,4.01,5.88,6.80,8.14,6.73,5.00,4.78,3.3,2.9,3.0]
}
X_NH3 = np.divide(expData['X_NH3'],100)
X_O2 = np.multiply(np.subtract(1,X_NH3), 0.21)
phi_data = np.divide(np.divide(X_NH3,X_O2),np.divide(4,3))
ax.plot(phi_data, expData['vel'], 'o', fillstyle='none', color='k', label='Ronney')
ax.legend(frameon=False, loc='upper right')
ax.set_ylabel(r'Burning velocity [cm $\rm s^{-1}$]')
ax.set_xlabel(r'Equivalence Ratio')
plt.show()