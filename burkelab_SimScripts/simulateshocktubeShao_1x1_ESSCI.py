from __future__ import division
from __future__ import print_function
import sys, os
sys.path.append("C:/Users/pjsin/Documents/cantera/build/python")
import cantera as ct
import matplotlib.pyplot as plt
import pandas as pd
import time
import scipy
import scipy.optimize
from scipy.optimize import curve_fit
import numpy as np
from matplotlib import gridspec
import argparse

import matplotlib as mpl
import sys, os
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

import matplotlib.ticker as ticker

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
save_plots = True
# figsize=(3.8,2)
f, ax = plt.subplots(1, 1, figsize=(args.figwidth, args.figheight)) 
name = 'ShockTubeSpeciesProfile_H2O' #os.path.splitext(os.path.basename(__file__))[0]

refSpecies='H2O'
X_H2O2 = 1163e-6
X_H2O = 1330e-6
X_O2 = 665e-6
X_CO2= 0.2*(1-X_H2O2-X_H2O-X_O2)
X_Ar = 1-X_CO2
def plotXvsTime(fname,pltlabel,pltcolour,lstyle='solid',zorder_value=10):
    # gas = ct.Solution('test/data/Burke_H2_ArBath.yaml')
    gas = ct.Solution(fname)
    gas.TPX = 1196, 2.127*101325, {'H2O2':X_H2O2, 'H2O':X_H2O, 'O2':X_O2, 'CO2':X_CO2, 'AR':X_Ar}
    r = ct.Reactor(contents=gas,energy="on")
    reactorNetwork = ct.ReactorNet([r]) # this will be the only reactor in the network
    timeHistory = ct.SolutionArray(gas, extra=['t'])
    estIgnitDelay = 1
    t = 0
    counter = 1
    while t < estIgnitDelay:
        t = reactorNetwork.step()
        if counter % 10 == 0:
            timeHistory.append(r.thermo.state, t=t)
        counter += 1
    tConv = 1e6 #time conversion factor (1e6 converts to microseconds)
    timeShift=0 # [seconds]
    shiftedTime = tConv*(timeHistory.t - timeShift)
    moleFrac = timeHistory(refSpecies).X 
    ax.plot(shiftedTime, moleFrac*100, color=pltcolour,label=pltlabel,linestyle=lstyle,linewidth=lw,zorder=zorder_value)

path="G:\\Mon disque\\Columbia\\Burke Lab\\01 Mixture Rules Project\\Graph Reading\\"
# path=os.getcwd()
def plotPoints(filename,mkr='none',mkrw='none',mkrsz='none',line='none',fill='none',colour='k',subplot='off',pltLabel="_hidden",zorder_value=10): 
    dataset = pd.read_csv(filename)
    ax.plot(dataset.iloc[:,0],dataset.iloc[:,1]*100,mkr,linewidth=0.7,fillstyle=fill,linestyle=line,color=colour,label=pltLabel,markersize=mkrsz,markeredgewidth=mkrw,zorder=zorder_value)

plotXvsTime("test/data/alzuetamechanism.yaml","Alzueta","xkcd:grey",zorder_value=90)
plotXvsTime("test\\data\\alzuetamechanism_epsNH3_T=300K.yaml",r"$\epsilon_{0,NH_3}(300K)$","orange",zorder_value=82)
plotXvsTime("test\\data\\alzuetamechanism_epsNH3_T=2000K.yaml",r"$\epsilon_{0,NH_3}(2000K)$","xkcd:teal",zorder_value=81)
plotXvsTime("test/data/alzuetamechanism_LMRR_allAR.yaml","Ar","r",zorder_value=70)
plotXvsTime("test/data/alzuetamechanism_LMRR_allH2O.yaml",r'$\rm H_2O$',"b",zorder_value=80)
plotXvsTime("test/data/alzuetamechanism_LMRR.yaml","LMR-R","xkcd:purple",zorder_value=100)
# plotPoints(path+'\\7 SP H2O X vs t (Shock Tube) (Shao)\\expData.csv',pltLabel='Shao et al.',line=':',colour='k')
plotPoints(path+'\\7 SP H2O X vs t (Shock Tube) (Shao)\\expData.csv',mkr='o',mkrsz=msz,pltLabel='Shao et al.',mkrw=mw,zorder_value=110)
# plotPoints(path+'\\7 SP H2O X vs t (Shock Tube) (Shao)\\troe_k0co2.csv',pltLabel='Troe et al.',line='solid',colour='g')
    
ax.legend(fontsize=lgdfsz,handlelength=lgdw, frameon=False, loc='lower right')  
ax.set_ylabel(r'$\rm H_2O$ mole fraction [%]')
ax.set_xlabel(r'Time [$\mathdefault{\mu s}$]')
ax.tick_params(axis='both', direction="in")#, labelsize=7)
ax.set_xlim([0.0001,299.999])
ax.set_ylim([0.120001,0.269999])

import matplotlib.ticker as ticker
ax.xaxis.set_major_locator(ticker.MultipleLocator(50))
ax.xaxis.set_major_formatter(ticker.StrMethodFormatter("{x:.0f}"))
ax.yaxis.set_major_locator(ticker.MultipleLocator(0.03))
ax.yaxis.set_major_formatter(ticker.StrMethodFormatter("{x:.2f}"))

if save_plots == True:
    plt.savefig('burkelab_SimScripts/figures/'+name+'_ESSCI.pdf', dpi=1000, bbox_inches='tight')
    plt.savefig('burkelab_SimScripts/figures/'+name+'_ESSCI.svg', dpi=1000, bbox_inches='tight')
# plt.show()     