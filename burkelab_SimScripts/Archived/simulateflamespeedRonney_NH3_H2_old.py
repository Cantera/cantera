
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

import sys, os
import matplotlib.pyplot as plt
import pandas as pd 
import numpy as np
import matplotlib as mpl
mpl.rc('font',family='Times New Roman')
mpl.rcParams['mathtext.fontset'] = 'stix'
from matplotlib.legend_handler import HandlerTuple

save_plots = True
fig, ax = plt.subplots(1, 1, figsize=(9, 5))
name = 'ronney_flamespeed'

models = {
          'LMR-R':"test/data/alzuetamechanism_LMRR.yaml", 
          'Ab initio':"test/data/alzuetamechanism.yaml",                 
          'Ar':"test/data/alzuetamechanism_LMRR_allAR.yaml",
          r'H$_2$O':"test/data/alzuetamechanism_LMRR_allH2O.yaml",
          }


import csv

# Function to save data to CSV
def save_to_csv(filename, data):
    with open(filename, 'w', newline='') as csvfile:
        writer = csv.writer(csvfile)
        writer.writerows(data)


# p_list = [50,100,250,760,1500]
p_list=[760]
NH3_list = np.linspace(0.175,0.32,5)

colors = ['xkcd:purple',"xkcd:grey",'r','b']
lines = ['-','-','--',':',(0, (5, 10)),(0, (3, 5, 1, 5, 1, 5))]

Tin = 300.0  # unburned gas temperature [K]
width = 0.03  # m
loglevel = 1  # amount of diagnostic output (0 to 8)

for k, m in enumerate(models):
    for i, p in enumerate(p_list):
        mbr = []
        phi_list = []
        for j, NH3 in enumerate(NH3_list):
            gas = ct.Solution(list(models.values())[k])
            ox_frac = 1 - NH3 # oxidizer fraction
            O2 = ox_frac*0.21
            N2 = ox_frac*0.79
            phi = np.divide(NH3/O2, 4/3)
            phi_list.append(phi)
            X = {'NH3':NH3, 'O2':O2, 'N2':N2}
            gas.TPX = Tin, (p/760)*ct.one_atm, X
            f = ct.FreeFlame(gas, width=width)
            f.set_refine_criteria(ratio=3, slope=0.06, curve=0.12)
            f.transport_model = 'mixture-averaged'
            # f.transport_model = 'multicomponent'
            f.solve(loglevel=loglevel, auto=True)
            mbr.append(f.velocity[0] * 100) # cm/s
        ax.plot(phi_list,mbr,label=m, color=colors[k],linestyle=lines[i])

        # Save phi_list and mbr to CSV
        csv_filename = f'{m}_{i}_data.csv'
        data = zip(phi_list, mbr)
        save_to_csv(csv_filename, data)
        
        
def plotPoints(fname, label, shape,color):
    dataset = pd.read_csv(fname)
    NH3_list = np.divide(dataset.iloc[:,0],100)
    ox_frac_list = np.subtract(1,NH3_list)
    O2_list = np.multiply(ox_frac_list, 0.21)
    phi_list = np.divide(np.divide(NH3_list,O2_list),np.divide(4,3))
    ax.plot(phi_list,dataset.iloc[:,1],shape,linestyle='none',color=color,label=label)
    
path="G:\\Mon disque\\Columbia\\Burke Lab\\01 Mixture Rules Project\\Graph Reading\\"
# # plotPoints(path+'\\6 FS NH3 (Stagni-Ronney)\\50torr.csv','50 torr','o','k')
# # plotPoints(path+'\\6 FS NH3 (Stagni-Ronney)\\100torr.csv','100 torr','^','k')
# # plotPoints(path+'\\6 FS NH3 (Stagni-Ronney)\\250torr.csv','250 torr','v','k')
# plotPoints(path+'\\6 FS NH3 (Stagni-Ronney)\\760torr.csv','Ronney','s','k')
# # plotPoints(path+'\\6 FS NH3 (Stagni-Ronney)\\1500torr.csv','1500 torr','D','k')

dataset=pd.read_csv(path+'\\6 FS NH3 (Stagni-Ronney)\\760torr.csv')
ax.plot(dataset.iloc[:,0],dataset.iloc[:,1]*100,marker='o',markersize=7,linewidth=3,fillstyle='none',linestyle='none',color='k',label='Ronney')

ax.legend(fontsize=15, frameon=False)#, loc='upper right')  
ax.set_ylabel(r'Burning velocity [cm $\rm s^{-1}$]', fontsize=18)
ax.set_xlabel(r'Equivalence Ratio', fontsize=18)
ax.tick_params(axis='both', direction="in", labelsize=15)
ax.tick_params(axis='both', which='minor', direction="in")

if save_plots == True:
    plt.savefig(name+'.pdf', dpi=1000, bbox_inches='tight')
    plt.savefig(name+'.png', dpi=1000, bbox_inches='tight')
plt.show()     