
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

mpl.rcParams['font.size'] = 10
mpl.rcParams['xtick.labelsize'] = 7
mpl.rcParams['ytick.labelsize'] = 7
from matplotlib.legend_handler import HandlerTuple
plt.rcParams['axes.labelsize'] = 10
mpl.rcParams['xtick.major.width'] = 0.5  # Width of major ticks on x-axis
mpl.rcParams['ytick.major.width'] = 0.5  # Width of major ticks on y-axis
mpl.rcParams['xtick.minor.width'] = 0.5  # Width of minor ticks on x-axis
mpl.rcParams['ytick.minor.width'] = 0.5  # Width of minor ticks on y-axis
mpl.rcParams['xtick.major.size'] = 2.5  # Length of major ticks on x-axis
mpl.rcParams['ytick.major.size'] = 2.5  # Length of major ticks on y-axis
mpl.rcParams['xtick.minor.size'] = 1.5  # Length of minor ticks on x-axis
mpl.rcParams['ytick.minor.size'] = 1.5  # Length of minor ticks on y-axis

save_plots = True
fig, ax = plt.subplots(1, 1, figsize=(3.8, 2))

import matplotlib.ticker as ticker


lw=0.7
mw=0.5
msz=3.5
dpi=1000
lgdw=0.6
lgdfsz=7
    
# ax[idxs[x]].yaxis.set_major_locator(ticker.MultipleLocator(5))
# ax[idxs[x]].xaxis.set_major_locator(ticker.MultipleLocator(50))
# ax[idxs[x]].xaxis.set_major_formatter(ticker.StrMethodFormatter("{x:.0f}"))
# ax[idxs[x]].yaxis.set_major_formatter(ticker.StrMethodFormatter("{x:.0f}"))

models = {
        'LMR-R':"test/data/alzuetamechanism_LMRR.yaml", 
        'Alzueta':"test/data/alzuetamechanism.yaml",                 
        'Ar':"test/data/alzuetamechanism_LMRR_allAR.yaml",
        r'H$_2$O':"test/data/alzuetamechanism_LMRR_allH2O.yaml",
        }
        
# def plotPoints(fname, label, shape,color,x):
#     dataset = pd.read_csv(fname)
#     NH3_list = np.divide(dataset.iloc[:,0],100)
#     ox_frac_list = np.subtract(1,NH3_list)
#     O2_list = np.multiply(ox_frac_list, 0.21)
#     phi_list = np.divide(np.divide(NH3_list,O2_list),np.divide(4,3))
#     ax[x].plot(phi_list,dataset.iloc[:,1],marker=shape,fillstyle='none',markersize=3.5,markeredgewidth=0.5,linestyle='none',color=color,label=label)
    
# path="G:\\Mon disque\\Columbia\\Burke Lab\\01 Mixture Rules Project\\Graph Reading\\"
# # plotPoints(path+'\\6 FS NH3 (Stagni-Ronney)\\50torr.csv','50 torr','o','k')
# # plotPoints(path+'\\6 FS NH3 (Stagni-Ronney)\\100torr.csv','100 torr','^','k')
# # plotPoints(path+'\\6 FS NH3 (Stagni-Ronney)\\250torr.csv','250 torr','v','k')
# plotPoints(path+'\\6 FS NH3 (Stagni-Ronney)\\760torr.csv','Ronney','o','k',x)
# # plotPoints(path+'\\6 FS NH3 (Stagni-Ronney)\\1500torr.csv','1500 torr','D','k')

# dataset=pd.read_csv(path+'\\6 FS NH3 (Stagni-Ronney)\\760torr.csv')
# ax.plot(dataset.iloc[:,0],dataset.iloc[:,1]*100,marker='o',markersize=7,linewidth=3,fillstyle='none',linestyle='none',color='k',label='Ronney')

path="C:\\Users\\pjsin\\Documents\\cantera\\burkelab_SimScripts\\RonneyResults_Mar04\\"

dataset=pd.read_csv(path+f'Ar_0_data_1alpha.csv')
ax.plot(dataset.iloc[2:,0],dataset.iloc[2:,1],linewidth=lw,color='r',label='Ar')

dataset=pd.read_csv(path+f'H2O_0_data_1alpha.csv')
ax.plot(dataset.iloc[2:,0],dataset.iloc[2:,1],linewidth=lw,color='b',label=r'H$_2$O')

dataset=pd.read_csv(path+f'Alzueta_0_data_1alpha.csv')
ax.plot(dataset.iloc[2:,0],dataset.iloc[2:,1],linewidth=lw,color="xkcd:grey",label='Alzueta')

dataset=pd.read_csv(path+f'LMR-R_0_data_1alpha.csv')
ax.plot(dataset.iloc[2:,0],dataset.iloc[2:,1],linewidth=lw,color='xkcd:purple',label='LMR-R')


path="G:\\Mon disque\\Columbia\\Burke Lab\\01 Mixture Rules Project\\Graph Reading\\"

dataset = pd.read_csv(path+'\\6 FS NH3 (Stagni-Ronney)\\760torr.csv')
NH3_list = np.divide(dataset.iloc[:,0],100)
ox_frac_list = np.subtract(1,NH3_list)
O2_list = np.multiply(ox_frac_list, 0.21)
phi_list = np.divide(np.divide(NH3_list,O2_list),np.divide(4,3))
ax.plot(phi_list,dataset.iloc[:,1],marker='o',fillstyle='none',markersize=msz,markeredgewidth=0.5,linestyle='none',color='k',label='Ronney')

# dataset = pd.read_csv(path+f'\\AlzuetaFig15\\1pt0_NH3.csv')
# ax.plot(dataset.iloc[:,0],dataset.iloc[:,1],marker='s',markersize=msz,markeredgewidth=0.5,linestyle='none',color='b',label='Alz Sim (Graph Read)')
ax.legend(fontsize=lgdfsz, frameon=False, loc='upper right') 

ax.set_ylabel(r'Burning velocity [cm $\rm s^{-1}$]')
ax.set_xlabel(r'Equivalence Ratio')
ax.tick_params(axis='both', direction="in")
ax.tick_params(axis='both', which='minor', direction="in")
ax.set_xlim([0.6, 1.8])

name = f'ronney_flamespeed_Mar04_3'
if save_plots == True:
    plt.savefig("C:\\Users\\pjsin\\Documents\\cantera\\burkelab_SimScripts\\figures\\Flame Speed Plots\\"+name+'.pdf', dpi=1000, bbox_inches='tight')
    plt.savefig("C:\\Users\\pjsin\\Documents\\cantera\\burkelab_SimScripts\\figures\\Flame Speed Plots\\"+name+'.png', dpi=dpi, bbox_inches='tight')

# plt.show()     