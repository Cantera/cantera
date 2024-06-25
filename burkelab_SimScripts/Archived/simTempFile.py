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

hfont = {'fontname':'sans-serif','fontweight':550,'fontsize':10,'fontstretch':500}

# # REPRODUCING BURKE FIGURE 1

path="G:\\Mon disque\\Columbia\\Burke Lab\\01 Mixture Rules Project\\Graph Reading\\"

def plotMassBurnRate(fname,labell,colour):
    # p_list = np.linspace(0,20,50)[1:]
    p_list = np.linspace(0,20,50)[1:]
    # Simulation parameters
    # p = ct.one_atm  # pressure [Pa]
    Tin = 300.0  # unburned gas temperature [K]
    reactants = 'H2:0.1071, O2:0.1785, He:0.7144'  # premixed gas composition
    width = 0.03  # m
    loglevel = 2  # amount of diagnostic output (0 to 8)
    mbr = []
    for p in p_list:
        gas = ct.Solution(fname)
        gas.TPX = Tin, p*ct.one_atm, reactants
        f = ct.FreeFlame(gas, width=width)
        f.set_refine_criteria(ratio=3, slope=0.06, curve=0.12)
        # f.set_max_jac_age(50, 50)
        f.transport_model = 'mixture-averaged'
        # f.transport_model = 'multicomponent'
        # f.set_max_time_step(1.0e-5)
        # f.set_max_steps(2000)
        # f.set_sensitivity_threshold(1.0e-6)
        # f.rtol = 1.0e-10
        # f.atol = 1.0e-15
        # f.transport_model = 'multicomponent'
        f.solve(loglevel=loglevel, auto=True)
        mbr.append(f.velocity[0]*f.density[0] / 10) # g/cm2*s
    plt.plot(p_list,mbr,label=labell,color=colour)

def plotPoints(filename,mkr='none',line='none',fill='none',colour='k',subplot='off',pltLabel="_hidden"):
    dataset = pd.read_csv(filename)
    plt.plot(dataset.iloc[:,0],dataset.iloc[:,1],mkr,fillstyle=fill,linestyle=line,color=colour,label=pltLabel)

plt.figure()
plotPoints(path+'5 FS H2O (Burke)\\black.csv',pltLabel='Keromnes et al.',line='solid',colour='k')
plotPoints(path+'5 FS H2O (Burke)\\exp_pts.csv',pltLabel='Burke et al.',mkr='s',fill='full',colour='k')
plotMassBurnRate('test/data/alzuetamechanism_LMRR.yaml', "LMR-R","m")
# plotMassBurnRate('test/data/alzuetamechanism_LMRR_allAR.yaml', "AR","r")
# plotMassBurnRate('test/data/alzuetamechanism_LMRR_allH2O.yaml', "H2O","b")
plt.legend()
plt.xlabel('Pressure [atm]')
plt.ylabel('Mass Burning Rate [g/cm2*s]')
plt.xlim([0,20])
plt.ylim([0,0.05])
plt.show()    