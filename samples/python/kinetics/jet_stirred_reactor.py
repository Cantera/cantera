"""
Simulate temperature profiles and species profiles in a jet-stirred reactor across a range of
initial temperatures, and observe the impact of incorporating the reduced-pressure linear mixture
rule (LMR-R) in such calculations. Code for transient simulations is also provided. 

Here we will consider a mixture of H2/O2/NH3/Ar (with 10% NH3) at 1.2 atm, and compare results
against the experimental measurements of Sabia et al. [1]

Two models are compared in this example:
(i)  A 2023 model of H2 and NH3 chemistry published by Alzueta et al. [2]
(ii) An adapted version of this model that has applied the reduced-pressure linear mixture
     rule (LMR-R) and ab initio third-body efficiencies [3]

References:
[1] P. Sabia, M. V. Manna, R. Ragucci, M. de Joannon, Mutual inhibition effect of hydrogen and ammonia
    in oxidation processes and the role of ammonia as “strong” collider in third-molecular reactions,
    Int. J. Hydrogen Energy 45 (2020) 32113 -- 32127.
[2] M. U. Alzueta, I. Salas, H. Hashemi, P. Glarborg, CO-assisted NH3 oxidation, Combust. Flame 257
    (2023) 112438.
[3] P. J. Singal, J. Lee, L. Lei, R. L. Speth, M. P. Burke, Implementation of New Mixture Rules
    Has a Substantial Impact on Combustion Predictions for H2 and NH3, Proc. Combust. Inst. 40
    (2024).

Requires: cantera >= 3.1
Keywords: jet-stirred reactor, species profile, temperature profile, mixture rule, LMR-R
"""

import sys, os
sys.path.append("C:/Users/pjsin/Documents/cantera/build/python")
import numpy as np
import pandas as pd
import time as time
import cantera as ct
import matplotlib.pyplot as plt

f, ax = plt.subplots(1, 3, figsize=(6.5, 2.5)) 
plt.subplots_adjust(wspace=0.6)

models = {    
          'Original':'alzueta.yaml',
          'LMR-R':'alzueta_LMRR.yaml',
          }

colours = ["xkcd:grey",'xkcd:purple']

T_list = np.linspace(800,1050,50) # [K]
reactorTemperature = 1000  # [K]
reactorPressure = 1.2  # [atm]
residenceTime = 0.5 # [s]
reactorVolume = 0.000113 # [m3]
reactorRadius = np.cbrt(reactorVolume*3/4/np.pi) # [m3]
reactorSurfaceArea = 4*np.pi*np.square(reactorRadius) # [m3]
pressureValveCoefficient = 2e-5
maxPressureRiseAllowed = 0.01
maxSimulationTime = 50  # [s]
heatTransferCoefficient = 79.5 # [W/m2/K]


def getTimeHistory(gas): 
    fuelAirMixtureTank = ct.Reservoir(gas)
    exhaust = ct.Reservoir(gas)
    env = ct.Reservoir(gas)
    stirredReactor = ct.IdealGasReactor(gas, energy='on', volume=reactorVolume)
    ct.MassFlowController(upstream=fuelAirMixtureTank,
                          downstream=stirredReactor,
                          mdot=stirredReactor.mass/residenceTime)
    ct.Valve(upstream=stirredReactor,
             downstream=exhaust,
             K=pressureValveCoefficient)
    ct.Wall(stirredReactor, env, A=reactorSurfaceArea, U=heatTransferCoefficient)
    reactorNetwork = ct.ReactorNet([stirredReactor])
    columnNames = ['pressure'] + [stirredReactor.component_name(item) for item in range(stirredReactor.n_vars)]
    timeHistory = pd.DataFrame(columns=columnNames)
    t = 0
    counter = 1
    while t < maxSimulationTime:
        t = reactorNetwork.step()
        if(counter%10 == 0):
            state = np.hstack([stirredReactor.thermo.P, stirredReactor.mass, 
                        stirredReactor.volume, stirredReactor.T, stirredReactor.thermo.X])
            timeHistory.loc[t] = state
        counter += 1
    return timeHistory

def getTemperatureDependence(gas):
    stirredReactor = ct.IdealGasReactor(gas, energy='on', volume=reactorVolume)
    columnNames = ['pressure'] + [stirredReactor.component_name(item) for item in range(stirredReactor.n_vars)]
    tempDependence = pd.DataFrame(columns=columnNames)
    for T in T_list:
        gas.TPX = T, reactorPressure*ct.one_atm, reactants
        fuelAirMixtureTank = ct.Reservoir(gas)
        exhaust = ct.Reservoir(gas)
        env = ct.Reservoir(gas)
        stirredReactor = ct.IdealGasReactor(gas, energy='on', volume=reactorVolume)
        ct.MassFlowController(upstream=fuelAirMixtureTank,
                          downstream=stirredReactor,
                          mdot=stirredReactor.mass/residenceTime)
        ct.Valve(upstream=stirredReactor,
                downstream=exhaust,
                K=pressureValveCoefficient)
        ct.Wall(stirredReactor, env, A=reactorSurfaceArea, U=heatTransferCoefficient)
        reactorNetwork = ct.ReactorNet([stirredReactor])
        t = 0
        while t < maxSimulationTime:
            t = reactorNetwork.step()
        state = np.hstack([stirredReactor.thermo.P, 
                        stirredReactor.mass, 
                        stirredReactor.volume, 
                        stirredReactor.T, 
                        stirredReactor.thermo.X])
        tempDependence.loc[T] = state
    return tempDependence

########################### TEMPERATURE DEPENDENCE ###################################
## Comment out this code if plotting time history instead of T-dependence
for k,m in enumerate(models):
    reactants = {'H2': 0.03, 'O2': 0.03, 'Ar': 0.846, 'NH3':0.094}
    gas = ct.Solution(list(models.values())[k])
    gas.TPX = reactorTemperature, reactorPressure*ct.one_atm, reactants
    tempDependence = getTemperatureDependence(gas)
    ax[0].plot(tempDependence.index, np.subtract(tempDependence['temperature'],tempDependence.index), color=colours[k],label=m) 
    ax[1].plot(tempDependence.index, tempDependence['O2']*100, color=colours[k])   
    ax[2].plot(tempDependence.index, tempDependence['H2']*100, color=colours[k])
    expData = {
        'T': [807,843,855,870,884,904,925,945,965,995,1018],
        'deltaT': [0.051,0.051,0.051,0.051,0.101,0.606,1.414,2.626,4.091,6.768,8.586],
        'X_O2': [3.076,3.053,3.050,3.037,3.024,3.015,2.966,2.924,2.794,2.597,2.261],
        'X_H2': [3.030,3.038,3.038,3.038,3.030,2.993,2.948,2.829,2.693,2.434,2.126]
    }
ax[0].plot(expData['T'],expData['deltaT'],'o',fillstyle='none',color='k',label="Sabia et al.")
ax[1].plot(expData['T'],expData['X_O2'],'o',fillstyle='none',color='k')
ax[2].plot(expData['T'],expData['X_H2'],'o',fillstyle='none',color='k')
ax[0].legend(fontsize=8,frameon=False,loc='upper left')
ax[0].set_ylabel('$\Delta$ T [K]')
ax[1].set_xlabel('Temperature [K]')
ax[1].set_ylabel('O$_2$ mole fraction [%]')
ax[2].set_ylabel('H$_2$ mole fraction [%]')
ax[0].set_xlim([780,1070])
ax[1].set_xlim([780,1070])
ax[2].set_xlim([780,1070])
plt.show()  

########################### TIME HISTORY ############################################
## Uncomment this code to plot the time history instead of T-dependence
# for k,m in enumerate(models):
#     reactants = {'H2': 0.03, 'O2': 0.03, 'Ar': 0.846, 'NH3':0.094}
#     gas = ct.Solution(list(models.values())[k])
#     gas.TPX = reactorTemperature, reactorPressure*ct.one_atm, reactants
#     timeHistory = getTimeHistory(gas)
#     ax[0].plot(timeHistory.index*1e3, np.subtract(timeHistory['temperature'],np.ones(len(timeHistory['temperature']))*reactorTemperature), color=colours[k])   
#     ax[1].plot(timeHistory.index*1e3, timeHistory['O2']*100, color=colours[k])   
#     ax[2].plot(timeHistory.index*1e3, timeHistory['H2']*100, color=colours[k])
# plt.show()  

   