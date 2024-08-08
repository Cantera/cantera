import sys, os
sys.path.append("C:/Users/pjsin/Documents/cantera/build/python")
import cantera as ct
import matplotlib.pyplot as plt
import pandas as pd 
import time
import numpy as np
import matplotlib as mpl
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('--figwidth', type=float, help="figwidth = ")
parser.add_argument('--figheight', type=float, help="figheight = ")
parser.add_argument('--fsz', type=float, help="mpl.rcParams['font.size'] = ", default=8)
parser.add_argument('--fszxtick', type=float, help="mpl.rcParams['xtick.labelsize'] = ", default=7)
parser.add_argument('--fszytick', type=float, help="mpl.rcParams['ytick.labelsize'] = ", default=7)
parser.add_argument('--fszaxlab', type=float, help="mpl.rcParams['axes.labelsize'] = ", default=8)
parser.add_argument('--lw', type=float, help="lw = ", default=0.7)
parser.add_argument('--mw', type=float, help="mw = ", default=0.5)
parser.add_argument('--msz', type=float, help="msz = ", default=2.5)
parser.add_argument('--lgdw', type=float, help="lgdw = ", default=0.6)
parser.add_argument('--lgdfsz', type=float, help="lgdw = ", default=5)
parser.add_argument('--gridsz', type=int, help="gridsz = ", default=10)
parser.add_argument('--dpi', type=int, help="dpi = ", default=1000)
args = parser.parse_args()
lw=args.lw
mw=args.mw
msz=args.msz
dpi=args.dpi
lgdw=args.lgdw
lgdfsz=args.lgdfsz
gridsz=args.gridsz

mpl.rc('font',family='Times New Roman')
mpl.rcParams['mathtext.fontset'] = 'stix'
mpl.rcParams['font.size'] = args.fsz
mpl.rcParams['xtick.labelsize'] = args.fszxtick
mpl.rcParams['ytick.labelsize'] = args.fszytick
from matplotlib.legend_handler import HandlerTuple
plt.rcParams['axes.labelsize'] = args.fszaxlab
mpl.rcParams['xtick.major.width'] = 0.5  # Width of major ticks on x-axis
mpl.rcParams['ytick.major.width'] = 0.5  # Width of major ticks on y-axis
mpl.rcParams['xtick.minor.width'] = 0.5  # Width of minor ticks on x-axis
mpl.rcParams['ytick.minor.width'] = 0.5  # Width of minor ticks on y-axis
mpl.rcParams['xtick.major.size'] = 2.5  # Length of major ticks on x-axis
mpl.rcParams['ytick.major.size'] = 2.5  # Length of major ticks on y-axis
mpl.rcParams['xtick.minor.size'] = 1.5  # Length of minor ticks on x-axis
mpl.rcParams['ytick.minor.size'] = 1.5  # Length of minor ticks on y-axis

f, ax = plt.subplots(1, 1, figsize=(args.figwidth, args.figheight)) 

# figsize = (6.5, 1.5) # (width, height)
# lw=0.7
# mw=0.5
# msz=2.5
# dpi=1000
# lgdw=0.6
# lgdfsz=5
# gridsz=10

save_plots = True

import matplotlib.ticker as ticker
plt.subplots_adjust(wspace=0.25)
ax.yaxis.set_major_locator(ticker.LogLocator(base=10.0, numticks=3))
ax.xaxis.set_major_locator(ticker.MultipleLocator(50))
ax.xaxis.set_major_formatter(ticker.StrMethodFormatter("{x:.0f}"))

name = 'IDT_keromnes'
path=os.getcwd()


def ignitionDelay(states, species):
    # i_ign = states(species).Y.argmax()
    i_ign = np.gradient(states(species).Y.T[0]).argmax()
    return states.t[i_ign]

################################################################################################
# path='G:\\Mon disque\\Columbia\\Burke Lab\\01 Mixture Rules Project\\Graph Reading\\'
# df = pd.read_csv(path+'10 IDT H2 CO (Keromnes)\\exp_50H2_50CO.csv')
# ax.semilogy(df['T'],df['IDT'],'o',fillstyle='none',linestyle='none',color='g',markersize=msz,markeredgewidth=mw,label='Exp 50H2/50CO', zorder=12)

# df = pd.read_csv(path+'10 IDT H2 CO (Keromnes)\\exp_85H2_85CO.csv')
# ax.semilogy(df['T'],df['IDT'],'o',fillstyle='none',linestyle='none',color='b',markersize=msz,markeredgewidth=mw,label='Exp 85H2/85CO', zorder=12)

# df = pd.read_csv(path+'10 IDT H2 CO (Keromnes)\\model_50H2_50CO.csv')
# ax.semilogy(df['T'],df['IDT'],'x',fillstyle='none',linestyle='none',color='g',markersize=msz,markeredgewidth=mw,label='Model 50H2/50CO', zorder=12)

# df = pd.read_csv(path+'10 IDT H2 CO (Keromnes)\\model_85H2_85CO.csv')
# ax.semilogy(df['T'],df['IDT'],'x',fillstyle='none',linestyle='none',color='b',markersize=msz,markeredgewidth=mw,label='Model 85H2/85CO', zorder=12)

#PLOT THE 50H2/50CO LINE
T_list = np.linspace(1000/8,1000/10,gridsz)#[::-1]
estimatedIgnitionDelayTimes = np.ones(len(T_list))
estimatedIgnitionDelayTimes[:] = 0.05
ignitionDelays_RG = np.zeros(len(T_list))
for j, T in enumerate(T_list):
    gas = ct.Solution("G:\\Mon disque\\Columbia\\Burke Lab\\07 Mechanisms\\01 Syngas (H2-CO)\\Keromnes-2013.yaml")
    # gas.TP = T, 16*ct.one_atm
    # fuel = {'H2':0.50, 'CO':0.50}
    # oxidizer = {'O2':1.0}
    # diluent = {'AR':0.50,'N2':0.50}
    # gas.set_equivalence_ratio(0.5, fuel, oxidizer, basis='mole', diluent=diluent, fraction='diluent:5/6')
    gas.TPX = T, 16*ct.one_atm, {'H2':0.0174, 'CO':0.0174, 'O2':0.0352, 'N2':0.4732, 'AR':0.4568}
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
ax.semilogy(np.divide(1000,T_list), 1e3*ignitionDelays_RG, '-', linestyle='solid', linewidth=lw, color="xkcd:grey", label="50H2/50CO")

#PLOT THE 85H2/15CO LINE
T_list = np.linspace(1000/8,1000/10,gridsz)#[::-1]
estimatedIgnitionDelayTimes = np.ones(len(T_list))
estimatedIgnitionDelayTimes[:] = 0.05
ignitionDelays_RG = np.zeros(len(T_list))
for j, T in enumerate(T_list):
    gas = ct.Solution("G:\\Mon disque\\Columbia\\Burke Lab\\07 Mechanisms\\01 Syngas (H2-CO)\\Keromnes-2013.yaml")
    # gas.TP = T, 16*ct.one_atm
    # fuel = {'H2':0.50, 'CO':0.50}
    # oxidizer = {'O2':1.0}
    # diluent = {'AR':0.50,'N2':0.50}
    # gas.set_equivalence_ratio(0.5, fuel, oxidizer, basis='mole', diluent=diluent, fraction='diluent:5/6')
    gas.TPX = T, 16*ct.one_atm, {'H2':0.0298, 'CO':0.0052, 'O2':0.0351, 'N2':0.4657, 'AR':0.4642}
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
ax.semilogy(np.divide(1000,T_list), 1e3*ignitionDelays_RG, '-', linestyle='solid', linewidth=lw, color='xkcd:purple', label="85H2/15CO")
    

ax.legend(fontsize=lgdfsz, frameon=False, loc='upper right',handlelength=lgdw)
ax.set_ylabel(r'Ignition delay [$\mathdefault{\mu s}$]')
ax.set_xlabel(r'Temperature [K]', fontsize=18)
ax.tick_params(axis='both', direction="in")
ax.tick_params(axis='both', which='minor', direction="in")#, bottom=False)

# plt.subplots_adjust(top=0.98)
if save_plots == True:
    plt.savefig('burkelab_SimScripts/figures/'+name+'.pdf', dpi=2000, bbox_inches='tight')
    plt.savefig('burkelab_SimScripts/figures/'+name+'.png', dpi=2000, bbox_inches='tight')
# plt.show()     