"""
Simulate ignition delay times for a shock tube across a range of initial temperatures,
and observe the impact of incorporating the reduced-pressure linear mixture rule (LMR-R)
in such calculations. 

Here we predict ignition delay times behind reflected shock waves for a mixture of
3% H2/1.5% O2/20% CO2/Ar at ~12 atm, and compare results against the experimental measurements
of Shao et al. [1]

Two models are compared in this example:
(i)  A 2023 model of H2 and NH3 chemistry published by Alzueta et al. [2]
(ii) An adapted version of this model that has applied the reduced-pressure linear mixture
     rule (LMR-R) and ab initio third-body efficiencies [3]

References:
[1] J. Shao, R. Choudhary, A. Susa, D. F. Davidson, R. K. Hanson, Shock tube study of the rate
    constants for H + O2 + M → HO2 + M (M = Ar, H2O, CO2, N2) at elevated pressures, Proc.
    Combust. Inst. 37 (2019) 145 -- 152.
[2] M. U. Alzueta, I. Salas, H. Hashemi, P. Glarborg, CO-assisted NH3 oxidation, Combust. Flame
    257 (2023) 112438.
[3] P. J. Singal, J. Lee, L. Lei, R. L. Speth, M. P. Burke, Implementation of New Mixture Rules
    Has a Substantial Impact on Combustion Predictions for H2 and NH3, Proc. Combust. Inst. 40
    (2024).

Requires: cantera >= 3.1
Keywords: ignition delay time, shock tube, mixture rule, LMR-R
"""

import sys, os
sys.path.append("C:/Users/pjsin/Documents/cantera/build/python")
import cantera as ct
import matplotlib.pyplot as plt
import time
import numpy as np

plt.figure()

models = {    
          'Original':'alzueta.yaml',
          'LMR-R':'alzueta_LMRR.yaml',
          }

colours = ["xkcd:grey",'xkcd:purple']

def ignitionDelay(states, species):
    i_ign = np.gradient(states(species).Y.T[0]).argmax()
    return states.t[i_ign]

expData = {
    'T': [1223, 1267.5, 1203, 1235, 1189, 1185.5, 1191, 1228, 1185, 1224, 1238],
    'IDT': [217,93,281,178,221,236,210,118,234,191,144]
}

T_list = np.linspace(1100,1300,9)
plt.semilogy(expData['T'],expData['IDT'],'o',fillstyle='none',linestyle='none',color='k',label='Shao et al.')
for k, m in enumerate(models):
    estimatedIgnitionDelayTimes = np.ones(len(T_list))
    estimatedIgnitionDelayTimes[:] = 0.05
    ignitionDelays_RG = np.zeros(len(T_list))
    for j, T in enumerate(T_list):
        gas = ct.Solution(list(models.values())[k])
        gas.TPX = T, 12*ct.one_atm, {'H2':0.03, 'O2':0.015, 'CO2':0.20, 'Ar':1-0.2-0.03-0.015}
        r = ct.Reactor(contents=gas)
        reactorNetwork = ct.ReactorNet([r])
        timeHistory = ct.SolutionArray(gas, extra=['t'])
        t0 = time.time()
        t = 0
        counter = 1
        while t < estimatedIgnitionDelayTimes[j]:
            t = reactorNetwork.step()
            if counter % 10 == 0:
                timeHistory.append(r.thermo.state, t=t)
            counter += 1
        tau = ignitionDelay(timeHistory, 'oh')
        t1 = time.time()
        ignitionDelays_RG[j] = tau
    plt.semilogy(T_list, 1e6*ignitionDelays_RG, color=colours[k], label=m)

plt.ylabel(r'Ignition delay [$\mathdefault{\mu s}$]')
plt.xlabel('Temperature [K]')
plt.title('Ignition delay times for 3% H$_2$/1.5% O$_2$/20% CO$_2$/Ar at ~12 atm',fontsize=10)

plt.legend(fontsize=10, frameon=False, loc='upper right')
plt.show()     