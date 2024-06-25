#%%
from __future__ import division
from __future__ import print_function
import sys, os
# sys.path.append("C:/Users/pjsin/Documents/cantera/build/python")
import cantera as ct
import matplotlib.pyplot as plt
import pandas as pd
import time
import scipy
import scipy.optimize
from scipy.optimize import curve_fit
import numpy as np
from matplotlib import gridspec
import pandas as pd
import numpy as np
import time
# import cantera as ct
import os.path
from os import path
import matplotlib.pyplot as plt
plt.rcParams.update(plt.rcParamsDefault)
import matplotlib as mpl

mpl.rc('font',family='Times New Roman')
mpl.rcParams['mathtext.fontset'] = 'stix'

# mpl.rcParams['font.size'] = 10
# mpl.rcParams['xtick.labelsize'] = 6
# mpl.rcParams['ytick.labelsize'] = 6
from matplotlib.legend_handler import HandlerTuple
# plt.rcParams['axes.labelsize'] = 8
# mpl.rcParams['xtick.major.width'] = 0.5  # Width of major ticks on x-axis
# mpl.rcParams['ytick.major.width'] = 0.5  # Width of major ticks on y-axis
# mpl.rcParams['xtick.minor.width'] = 0.5  # Width of minor ticks on x-axis
# mpl.rcParams['ytick.minor.width'] = 0.5  # Width of minor ticks on y-axis
# mpl.rcParams['xtick.major.size'] = 2.5  # Length of major ticks on x-axis
# mpl.rcParams['ytick.major.size'] = 2.5  # Length of major ticks on y-axis
# mpl.rcParams['xtick.minor.size'] = 1.5  # Length of minor ticks on x-axis
# mpl.rcParams['ytick.minor.size'] = 1.5  # Length of minor ticks on y-axis

save_plots = True

# f, ax = plt.subplots(1, 2, figsize=(9, 5))

f = plt.figure(figsize=(9,7))


spec = gridspec.GridSpec(ncols=1, nrows=2, hspace=0, height_ratios=[3, 2])

# import matplotlib.ticker as ticker
# plt.subplots_adjust(hspace=0.3)
# plt.subplots_adjust()
# ax.yaxis.set_major_locator(ticker.MultipleLocator(5))
# ax[1].yaxis.set_major_locator(ticker.MultipleLocator(0.25))
# ax[2].yaxis.set_major_locator(ticker.MultipleLocator(0.5))
# ax.xaxis.set_major_locator(ticker.MultipleLocator(50))
# ax[1].xaxis.set_major_locator(ticker.MultipleLocator(50))
# ax[2].xaxis.set_major_locator(ticker.MultipleLocator(50))
# ax.xaxis.set_major_formatter(ticker.StrMethodFormatter("{x:.0f}"))
# ax.yaxis.set_major_formatter(ticker.StrMethodFormatter("{x:.0f}"))
# ax[1].xaxis.set_major_formatter(ticker.StrMethodFormatter("{x:.0f}"))
# ax[1].yaxis.set_major_formatter(ticker.StrMethodFormatter("{x:.1f}"))
# ax[2].xaxis.set_major_formatter(ticker.StrMethodFormatter("{x:.0f}"))
# ax[2].yaxis.set_major_formatter(ticker.StrMethodFormatter("{x:.1f}"))
# ax[0].annotate('(d)', xy=(0.95, 0.95), xycoords='axes fraction',ha='right', va='top')
# ax[1].annotate('(e)', xy=(0.95, 0.95), xycoords='axes fraction',ha='right', va='top')
# ax[2].annotate('(f)', xy=(0.95, 0.95), xycoords='axes fraction',ha='right', va='top')
lw=0.7
mw=0.5
msz=3.5
dpi=1000
lgdw=1
lgdfsz=6

name = 'JSR_net_rates_of_progress_2'
models = {
    
            'glarborg':'test/data/reduced_glarborg.cti'
        #   'LMR-R':"test/data/alzuetamechanism_LMRR.yaml", 
        #   'Alzueta':"test/data/alzuetamechanism.yaml",            
        #   'Ar':"test/data/alzuetamechanism_LMRR_allAR.yaml",
        #   r'H$_2$O':"test/data/alzuetamechanism_LMRR_allH2O.yaml",
          }

T_list = np.concatenate((np.linspace(800,829.9,10),np.linspace(830,844,100),np.linspace(848.8,870.9,100),np.linspace(871,1050,50)))
# T_list = np.linspace(800,1050,50)
P = 1.2
tau = 0.5

# NH3_con_list = [0.02, 0.05, 0.075, 0.10]
# H2O_con_list = [0, 0.05, 0.10, 0.15, 0.20]
diluent = 0.94
# H2Opercent_list = [0, 0.05, 0.10, 0.15, 0.20]
H2Opercent_list = [0.20]

colors = ['xkcd:purple',"xkcd:grey",'r','b']
lines =['-','--','-','-','-']


reactions = ['H + O2 <=> O + OH', 'H + O2 (+M) <=> HO2 (+M)', 
             'H + HO2 <=> H2 + O2', 'H + HO2 <=> 2 OH',
             'HO2 + OH <=> H2O + O2', 'HO2 + OH <=> H2O + O2']


##############################################################################################################################

reactorTemperature = 1000  # Kelvin
reactorPressure = P*ct.one_atm  # in atm. This equals 1.06 bars
residenceTime = tau  # s
reactorVolume = 0.000113 #30.5*(1e-2)**3  # m3
reactorRadius = np.cbrt(reactorVolume*3/4*np.pi) # m3
reactorSurfaceArea = 4*np.pi*np.square(reactorRadius) # m3
pressureValveCoefficient = 0.01
maxPressureRiseAllowed = 0.01
maxSimulationTime = 50  # seconds
heatTransferCoefficient = 7949.6
heatTransferCoefficient = 7.9496*2.2
tempDependence = []

##############################################################################################################################

for k,m in enumerate(models):
    for i, H2Opercent in enumerate(H2Opercent_list):
        H2O = diluent * H2Opercent
        Ar = diluent * (1-H2Opercent)
        # reactants = {'H2': 0.03, "H":1e-6, 'O2': 0.03, 'AR': Ar, 'H2O':H2O}    
        reactants = {'H2': 0.03, 'O2': 0.03, 'AR': Ar, 'H2O':H2O}    
        
        # for j, T in enumerate(T_list):
        
        concentrations = reactants
        gas = ct.Solution(list(models.values())[k])
        gas.TPX = reactorTemperature, reactorPressure, concentrations 

        fuelAirMixtureTank = ct.Reservoir(gas)
        exhaust = ct.Reservoir(gas)
        env = ct.Reservoir(gas)
        stirredReactor = ct.IdealGasReactor(gas, energy='on', volume=reactorVolume)
        massFlowController = ct.MassFlowController(upstream=fuelAirMixtureTank,
                                                    downstream=stirredReactor,
                                                    mdot=stirredReactor.mass/residenceTime)
        pressureRegulator = ct.Valve(upstream=stirredReactor,
                                    downstream=exhaust,
                                    K=pressureValveCoefficient)
        w2 = ct.Wall(stirredReactor, env, A=reactorSurfaceArea, U=heatTransferCoefficient)
        reactorNetwork = ct.ReactorNet([stirredReactor])
        columnNames = [stirredReactor.component_name(item) for item in range(stirredReactor.n_vars)]
        
        columnNames = ['pressure'] + columnNames
        timeHistory = pd.DataFrame(columns=columnNames)
        tic = time.time()
        # reactorNetwork.rtol = 1.0e-6
        # reactorNetwork.atol = 1.0e-15
        t = 0
        counter = 1
        while t < maxSimulationTime:
            t = reactorNetwork.step()
            if(counter%10 == 0):
                state = np.hstack([stirredReactor.thermo.P, stirredReactor.mass, 
                            stirredReactor.volume, stirredReactor.T, stirredReactor.thermo.X])
                timeHistory.loc[t] = state
            counter += 1
        toc = time.time()
        
        pressureDifferential = timeHistory['pressure'].max()-timeHistory['pressure'].min()
        if(abs(pressureDifferential/reactorPressure) > maxPressureRiseAllowed):
            print("WARNING: Non-trivial pressure rise in the reactor. Adjust K value in valve")
            
        tempDependence.append(pd.DataFrame(columns=timeHistory.columns))
        tempDependence[i].index.name = 'Temperature'
        
        net_rates=[]    
        inletConcentrations = concentrations        
        for j,T in enumerate(T_list): #temperature in T:
            reactorTemperature = T #temperature  # Kelvin
            gas.TPX = reactorTemperature, reactorPressure, inletConcentrations
            timeHistory = pd.DataFrame(columns=columnNames)
            fuelAirMixtureTank = ct.Reservoir(gas)
            exhaust = ct.Reservoir(gas)
            env = ct.Reservoir(gas)
            # gas.TPX = reactorTemperature, reactorPressure, concentrations
            stirredReactor = ct.IdealGasReactor(gas, energy='on', volume=reactorVolume)
            
            massFlowController = ct.MassFlowController(upstream=fuelAirMixtureTank,
                                                    downstream=stirredReactor,
                                                    mdot=stirredReactor.mass/residenceTime)
            pressureRegulator = ct.Valve(upstream=stirredReactor, 
                                        downstream=exhaust, 
                                        K=pressureValveCoefficient)
            w2 = ct.Wall(stirredReactor, env, A=reactorSurfaceArea, U=heatTransferCoefficient)
            reactorNetwork = ct.ReactorNet([stirredReactor])
            tic = time.time()
            t = 0
            while t < maxSimulationTime:
                t = reactorNetwork.step()
            state = np.hstack([stirredReactor.thermo.P, 
                            stirredReactor.mass, 
                            stirredReactor.volume, 
                            stirredReactor.T, 
                            stirredReactor.thermo.X])
            toc = time.time()
            concentrations = stirredReactor.thermo.X
            tempDependence[i].loc[T] = state
            net_rates.append(list(gas.net_rates_of_progress))
            
        net_rates = np.array(net_rates).T.tolist()
            
        ax1 = f.add_subplot(spec[0])
        ax1.plot(tempDependence[i].index, np.add(net_rates[gas.reaction_equations().index(reactions[0])],net_rates[gas.reaction_equations().index(reactions[1])]), color='r', linestyle='solid', linewidth=2, label='H + O2')
        ax1.plot(tempDependence[i].index, np.add(net_rates[gas.reaction_equations().index(reactions[2])],net_rates[gas.reaction_equations().index(reactions[3])]), color='b', linestyle='solid', linewidth=2, label='H + HO2')
        ax1.plot(tempDependence[i].index, np.add(net_rates[gas.reaction_equations().index(reactions[4])],net_rates[gas.reaction_equations().index(reactions[4])+1]), color='g', linestyle='solid', linewidth=2, label='OH + HO2')
        # ax1.plot(tempDependence[i].index, np.add(net_rates[gas.reaction_equations().index('H + O2 <=> O + OH')],net_rates[gas.reaction_equations().index('H + O2 (+M) <=> HO2 (+M)')]), color='r', linestyle='solid', linewidth=2, label='H + O2')
        # ax1.plot(tempDependence[i].index, np.add(net_rates[gas.reaction_equations().index('H + O2 <=> O + OH')],net_rates[gas.reaction_equations().index('H + O2 (+M) <=> HO2 (+M)')]), color='r', linestyle='solid', linewidth=2, label='H + O2')

        # ax1.set_xlabel('Temperature [K]')
        ax1.set_ylabel(r'Rate of Progress [kmol/m$^{3}$/s]')
        ax1.tick_params(axis='both',direction='in')
        ax1.legend(frameon='off')        
        
        ax2 = f.add_subplot(spec[1])
        # ax1.plot(tempDependence[i].index, net_rates[gas.reaction_equations().index(reactions[0])], color='r', linestyle='solid', linewidth=2, label=reactions[0])
        # ax1.plot(tempDependence[i].index, net_rates[gas.reaction_equations().index(reactions[1])], color='b', linestyle='solid', linewidth=2, label=reactions[1]) 
        # ax2.plot(tempDependence[i].index, np.multiply(np.divide(np.add(net_rates[gas.reaction_equations().index(reactions[2])],net_rates[gas.reaction_equations().index(reactions[3])]),np.multiply(2,net_rates[gas.reaction_equations().index(reactions[2])])),np.divide(np.multiply(2,net_rates[gas.reaction_equations().index(reactions[0])]),net_rates[gas.reaction_equations().index(reactions[1])])), color='k', linestyle='solid', linewidth=2)   
        # ax[1].plot(tempDependence[i].index, tempDependence[i]['O2']*100, color=colors[k], linestyle='solid', linewidth=lw, label=m)   
        # ax[2].plot(tempDependence[i].index, tempDependence[i]['H2']*100, color=colors[k], linestyle='solid', linewidth=lw, label=m) 
        ax2.set_xlabel('Temperature [K]')
        ax2.set_ylabel(r'Ratio')
        ax2.tick_params(axis='both',direction='in')
        # ax2.legend()         

# path="G:\\Mon disque\\Columbia\\Burke Lab\\01 Mixture Rules Project\\Graph Reading\\1 JSR H2O\\"
 
# T_20_data = pd.read_csv(path+'JSR_T_H2O_20_data.csv') 
# O2_20_data = pd.read_csv(path+'JSR_O2_H2O_20_data.csv') 
# H2_20_data = pd.read_csv(path+'JSR_H2_H2O_20_data.csv') 
# ax[0].plot(T_20_data.iloc[:, 0],T_20_data.iloc[:, 1],marker='o',fillstyle='none',linestyle='none',color='k',markersize=msz,markeredgewidth=mw, label="Sabia et al.")
# ax[1].plot(O2_20_data.iloc[:, 0],O2_20_data.iloc[:, 1],marker='o',fillstyle='none',linestyle='none',color='k',markersize=msz,markeredgewidth=mw,label="Sabia et al.")
# ax[2].plot(H2_20_data.iloc[:, 0],H2_20_data.iloc[:, 1],marker='o',fillstyle='none',linestyle='none',color='k',markersize=msz,markeredgewidth=mw, label="Sabia et al.")


# ax[0].legend(frameon=False)#,loc='lower right')
            
# ax[1].set_xlabel('Temperature [K]')
# ax[1].set_ylabel('O$_2$ mole fraction [%]')
# ax[1].tick_params(axis='both',direction='in')
# # ax[1].legend(frameon=False)#,loc='upper right')

# ax[2].set_xlabel('Temperature [K]')
# ax[2].set_ylabel('H$_2$ mole fraction [%]')
# ax[2].tick_params(axis='both',direction='in')



if save_plots == True:
    plt.savefig(name+'.pdf', dpi=1000, bbox_inches='tight')
    plt.savefig(name+'.png', dpi=dpi, bbox_inches='tight')
plt.show()     