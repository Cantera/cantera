#%%
import sys, os
sys.path.append("C:/Users/pjsin/Documents/cantera/build/python")
import cantera as ct
import matplotlib.pyplot as plt
import pandas as pd 
import time
import numpy as np
import matplotlib as mpl
mpl.rc('font',family='Times New Roman')
mpl.rcParams['mathtext.fontset'] = 'stix'

mpl.rcParams['font.size'] = 8
mpl.rcParams['xtick.labelsize'] = 7
mpl.rcParams['ytick.labelsize'] = 7
from matplotlib.legend_handler import HandlerTuple
plt.rcParams['axes.labelsize'] = 8
mpl.rcParams['xtick.major.width'] = 0.5  # Width of major ticks on x-axis
mpl.rcParams['ytick.major.width'] = 0.5  # Width of major ticks on y-axis
mpl.rcParams['xtick.minor.width'] = 0.5  # Width of minor ticks on x-axis
mpl.rcParams['ytick.minor.width'] = 0.5  # Width of minor ticks on y-axis
mpl.rcParams['xtick.major.size'] = 2.5  # Length of major ticks on x-axis
mpl.rcParams['ytick.major.size'] = 2.5  # Length of major ticks on y-axis
mpl.rcParams['xtick.minor.size'] = 1.5  # Length of minor ticks on x-axis
mpl.rcParams['ytick.minor.size'] = 1.5  # Length of minor ticks on y-axis

save_plots = True
figsize = (5.3333, 2.6) # (width, height)
fig, ax = plt.subplots(2, 2, figsize=figsize)

import matplotlib.ticker as ticker
plt.subplots_adjust(wspace=0.15)
ax[0,0].yaxis.set_major_locator(ticker.LogLocator(base=10.0, numticks=3))
ax[0,1].yaxis.set_major_locator(ticker.LogLocator(base=10.0, numticks=3))
ax[1,0].yaxis.set_major_locator(ticker.LogLocator(base=10.0, numticks=3))
ax[1,1].yaxis.set_major_locator(ticker.LogLocator(base=10.0, numticks=3))

ax[0,0].xaxis.set_major_locator(ticker.MultipleLocator(50))
ax[0,1].xaxis.set_major_locator(ticker.MultipleLocator(50))
ax[1,0].xaxis.set_major_locator(ticker.MultipleLocator(50))
ax[1,1].xaxis.set_major_locator(ticker.MultipleLocator(50))

ax[0,0].xaxis.set_major_formatter(ticker.StrMethodFormatter("{x:.0f}"))
ax[0,0].yaxis.set_major_formatter(ticker.StrMethodFormatter("{x:.0e}"))
ax[0,1].xaxis.set_major_formatter(ticker.StrMethodFormatter("{x:.0f}"))
ax[0,1].yaxis.set_major_formatter(ticker.StrMethodFormatter("{x:.1e}"))
ax[1,0].xaxis.set_major_formatter(ticker.StrMethodFormatter("{x:.0f}"))
ax[1,0].yaxis.set_major_formatter(ticker.StrMethodFormatter("{x:.1e}"))
ax[1,1].xaxis.set_major_formatter(ticker.StrMethodFormatter("{x:.0f}"))
ax[1,1].yaxis.set_major_formatter(ticker.StrMethodFormatter("{x:.1e}"))

lw=0.7
mw=0.5
msz=2.5
dpi=3000
lgdw=0.6
lgdfsz=5
gridsz=15


name = 'IDT_shao'
path=os.getcwd()

# models = {
#           'LMR-R':"test/data/alzuetamechanism_LMRR.yaml", 
#           'A priori':"test/data/alzuetamechanism.yaml",             
#           'Ar':"test/data/alzuetamechanism_LMRR_allAR.yaml",
#           r'H$_2$O':"test/data/alzuetamechanism_LMRR_allH2O.yaml",
#           }

# colors = ['xkcd:purple',"xkcd:grey",'r','b']

models = {
          'LMR-R':"test/data/alzuetamechanism_LMRR.yaml",          
          'Ar':"test/data/alzuetamechanism_LMRR_allAR.yaml",
          r'H$_2$O':"test/data/alzuetamechanism_LMRR_allH2O.yaml",
          }

colors = ['xkcd:purple','r','b']

def ignitionDelay(states, species):
    # i_ign = states(species).Y.argmax()
    i_ign = np.gradient(states(species).Y.T[0]).argmax()
    return states.t[i_ign]

################################################################################################

path='G:\\Mon disque\\Columbia\\Burke Lab\\01 Mixture Rules Project\\Graph Reading\\'
df = pd.read_csv(path+'Shao_IDT\\1.csv')
p_df = df['P']
T_df = df['T']
IDT_df = df['IDT']

T_list = np.linspace(1100,1300,gridsz)#[::-1]
for k, m in enumerate(models):
    estimatedIgnitionDelayTimes = np.ones(len(T_list))
    estimatedIgnitionDelayTimes[:] = 0.05
    ignitionDelays_RG = np.zeros(len(T_list))
    for j, T in enumerate(T_list):
        gas = ct.Solution(list(models.values())[k])
        gas.TPX = T, 17*ct.one_atm, {'H2':0.03, 'O2':0.015, 'Ar':1-0.03-0.015}
        r = ct.Reactor(contents=gas)
        reactorNetwork = ct.ReactorNet([r])
        timeHistory = ct.SolutionArray(gas, extra=['t'])
        t0 = time.time()
        t = 0
        counter = 1
        while t < estimatedIgnitionDelayTimes[j]:
            t = reactorNetwork.step()
            if counter % 1 == 0:
                timeHistory.append(r.thermo.state, t=t)
            counter += 1
        tau = ignitionDelay(timeHistory, 'oh')
        t1 = time.time()
        ignitionDelays_RG[j] = tau
    ax[0, 0].semilogy(T_list, 1e6*ignitionDelays_RG, '-', linestyle='solid', linewidth=lw, color=colors[k], label=m)
ax[0, 0].semilogy(T_df,IDT_df,'o',fillstyle='none',linestyle='none',color='k',markersize=msz,markeredgewidth=mw,label='Shao et al.')
ax[0,0].legend(fontsize=lgdfsz, frameon=False, loc='lower left',handlelength=lgdw)  
ax[0, 0].set_ylabel(r'Ignition delay [$\mathdefault{\mu s}$]')
# ax[0, 0].set_xlabel(r'Temperature [K]', fontsize=18)
ax[0, 0].tick_params(axis='both', direction="in")
ax[0, 0].tick_params(axis='both', which='minor', direction="in")#, bottom=False)
ax[0, 0].annotate('(a)', xy=(0.95, 0.9), xycoords='axes fraction', ha='right', va='top')

################################################################################################

df = pd.read_csv(path+'\\Shao_IDT\\2.csv')
p_df = df['P']
T_df = df['T']
IDT_df = df['IDT']
T_list = np.linspace(1100,1300,gridsz)#[::-1]
for k, m in enumerate(models):
    estimatedIgnitionDelayTimes = np.ones(len(T_list))
    estimatedIgnitionDelayTimes[:] = 0.05
    ignitionDelays_RG = np.zeros(len(T_list))
    for j, T in enumerate(T_list):
        gas = ct.Solution(list(models.values())[k])
        gas.TPX = T, 13*ct.one_atm, {'H2':0.03, 'O2':0.015, 'N2':1-0.03-0.015}
        r = ct.Reactor(contents=gas)
        reactorNetwork = ct.ReactorNet([r])
        timeHistory = ct.SolutionArray(gas, extra=['t'])
        t0 = time.time()
        t = 0
        counter = 1
        while t < estimatedIgnitionDelayTimes[j]:
            t = reactorNetwork.step()
            if counter % 1 == 0:
                timeHistory.append(r.thermo.state, t=t)
            counter += 1
        tau = ignitionDelay(timeHistory, 'oh')
        t1 = time.time()
        ignitionDelays_RG[j] = tau
    ax[0,1].semilogy(T_list, 1e6*ignitionDelays_RG, '-', linestyle='solid',linewidth=lw, color=colors[k], label=m)
ax[0,1].semilogy(T_df,IDT_df,'o',fillstyle='none',linestyle='none',color='k',markersize=msz,markeredgewidth=mw,label='Shao et al.')
# ax[1].legend(fontsize=15, frameon=False)#, loc='upper right')  
# ax[0,1].set_ylabel(r'Ignition delay [$\mathdefault{\mu s}$]', fontsize=18)
# ax[0,1].set_xlabel(r'Temperature [K]', fontsize=18)
# ax[0,1].legend(fontsize=10, frameon=False, loc='upper right')  
ax[0,1].tick_params(axis='both', direction="in")
ax[0,1].tick_params(axis='both', which='minor', direction="in")#, bottom=False)
ax[0,1].annotate('(b)', xy=(0.95, 0.9), xycoords='axes fraction', ha='right', va='top')

################################################################################################

df = pd.read_csv(path+'\\Shao_IDT\\3.csv')
p_df = df['P']
T_df = df['T']
IDT_df = df['IDT']
H2O_df = df['H2O']
T_list = np.linspace(1200,1400,gridsz)#[::-1]
for k, m in enumerate(models):
    estimatedIgnitionDelayTimes = np.ones(len(T_list))
    estimatedIgnitionDelayTimes[:] = 0.05
    ignitionDelays_RG = np.zeros(len(T_list))
    for j, T in enumerate(T_list):
        gas = ct.Solution(list(models.values())[k])
        gas.TPX = T, 15*ct.one_atm, {'H2':0.03, 'O2':0.015, 'H2O':0.09, 'Ar':1-0.03-0.015-0.09}
        r = ct.Reactor(contents=gas)
        reactorNetwork = ct.ReactorNet([r])
        timeHistory = ct.SolutionArray(gas, extra=['t'])
        t0 = time.time()
        t = 0
        counter = 1
        while t < estimatedIgnitionDelayTimes[j]:
            t = reactorNetwork.step()
            if counter % 1 == 0:
                timeHistory.append(r.thermo.state, t=t)
            counter += 1
        tau = ignitionDelay(timeHistory, 'oh')
        t1 = time.time()
        ignitionDelays_RG[j] = tau
    ax[1,0].semilogy(T_list, 1e6*ignitionDelays_RG, '-', linestyle='solid',linewidth=lw, color=colors[k], label=m)
ax[1,0].semilogy(T_df,IDT_df,'o',fillstyle='none',linestyle='none',color='k',markersize=msz,markeredgewidth=mw,label='Shao et al.')
# ax[2].legend(fontsize=15, frameon=False)#, loc='upper right')  
ax[1,0].set_ylabel(r'Ignition delay [$\mathdefault{\mu s}$]')
ax[1,0].set_xlabel(r'Temperature [K]')
ax[1,0].tick_params(axis='both', direction="in")
ax[1,0].tick_params(axis='both', which='minor', direction="in")#, bottom=False)
ax[1,0].annotate('(c)', xy=(0.95, 0.9), xycoords='axes fraction',ha='right', va='top')

################################################################################################

df = pd.read_csv(path+'\\Shao_IDT\\4.csv')
p_df = df['P']
T_df = df['T']
IDT_df = df['IDT']
T_list = np.linspace(1100,1300,gridsz)#[::-1]
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
            if counter % 1 == 0:
                timeHistory.append(r.thermo.state, t=t)
            counter += 1
        tau = ignitionDelay(timeHistory, 'oh')
        t1 = time.time()
        ignitionDelays_RG[j] = tau
    ax[1,1].semilogy(T_list, 1e6*ignitionDelays_RG, '-', linestyle='solid',linewidth=lw, color=colors[k], label=m)
ax[1,1].semilogy(T_df,IDT_df,'o',fillstyle='none',linestyle='none',color='k',markersize=msz,markeredgewidth=mw,label='Shao et al.')
# ax[3].legend(fontsize=15, frameon=False)#, loc='upper right')  
# ax[1,1].set_ylabel(r'Ignition delay [$\mathdefault{\mu s}$]', fontsize=18)
ax[1,1].set_xlabel(r'Temperature [K]')
ax[1,1].tick_params(axis='both', direction="in")
ax[1,1].tick_params(axis='both', which='minor', direction="in")#, bottom=False)
ax[1,1].annotate('(d)', xy=(0.95, 0.9), xycoords='axes fraction',ha='right', va='top')
# plt.subplots_adjust(wspace=0.4, hspace=0.4)

if save_plots == True:
    plt.savefig(name+'.pdf', dpi=1000, bbox_inches='tight')
    plt.savefig(name+'.png', dpi=dpi, bbox_inches='tight')
plt.show()     