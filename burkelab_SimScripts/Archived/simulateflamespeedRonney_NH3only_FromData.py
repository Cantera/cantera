
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
parser.add_argument('--date', type=str, help="sim date = ",default='May28')
parser.add_argument('--slopeVal', type=float, help="slope value = ",default=-1)
parser.add_argument('--curveVal', type=float, help="curve value = ",default=-1)
parser.add_argument('--paper', type=str, help="paper = ",default='PCI')

args = parser.parse_args()
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
fig, ax = plt.subplots(1,1,figsize=(args.figwidth, args.figheight))
# plt.subplots_adjust(wspace=0.4, hspace=1)


lw=args.lw
mw=args.mw
msz=args.msz
dpi=args.dpi
lgdw=args.lgdw
lgdfsz=args.lgdfsz
date=args.date
fslope=args.slopeVal
fcurve=args.curveVal
import matplotlib.ticker as ticker
ax.xaxis.set_major_locator(ticker.MultipleLocator(0.2))
ax.xaxis.set_major_formatter(ticker.StrMethodFormatter("{x:.1f}"))
ax.yaxis.set_major_locator(ticker.MultipleLocator(2))
ax.yaxis.set_major_formatter(ticker.StrMethodFormatter("{x:.0f}"))



if fslope != -1:
    path="C:\\Users\\pjsin\\Documents\\cantera\\burkelab_SimScripts\\RonneyResults_"+date+f' (slope={fslope} curve={fcurve})\\'
else:
    path="C:\\Users\\pjsin\\Documents\\cantera\\burkelab_SimScripts\\RonneyResults_"+date+"\\"


dataset=pd.read_csv(path+f'Ar_0_data_1.0alpha.csv')
ax.plot(dataset.iloc[:,0],dataset.iloc[:,1],linewidth=lw,color='r',label='Ar')

dataset=pd.read_csv(path+f'H2O_0_data_1.0alpha.csv')
ax.plot(dataset.iloc[:,0],dataset.iloc[:,1],linewidth=lw,color='b',label=r'H$_2$O')

if args.paper == 'ESSCI':
    dataset=pd.read_csv(path+f'Alzueta_0_data_1.0alpha.csv')
    ax.plot(dataset.iloc[:,0],dataset.iloc[:,1],linewidth=lw,color="xkcd:grey",label='Alzueta')

dataset=pd.read_csv(path+f'LMR-R_0_data_1.0alpha.csv')
ax.plot(dataset.iloc[:,0],dataset.iloc[:,1],linewidth=lw,color='xkcd:purple',label='LMR-R')

# ax.set_title(f'{round(alpha*100)}% NH3/{round((1-alpha)*100)}% H2')

path="G:\\Mon disque\\Columbia\\Burke Lab\\01 Mixture Rules Project\\Graph Reading\\"

dataset = pd.read_csv(path+'\\6 FS NH3 (Stagni-Ronney)\\760torr.csv')
NH3_list = np.divide(dataset.iloc[:,0],100)
ox_frac_list = np.subtract(1,NH3_list)
O2_list = np.multiply(ox_frac_list, 0.21)
phi_list = np.divide(np.divide(NH3_list,O2_list),np.divide(4,3))
ax.plot(phi_list,dataset.iloc[:,1],marker='o',fillstyle='none',markersize=msz,markeredgewidth=0.5,linestyle='none',color='k',label='Ronney')
ax.legend(fontsize=lgdfsz, frameon=False, loc='upper right') 
ax.set_xlabel(r'Equivalence Ratio')
ax.set_ylabel(r'Burning velocity [cm $\rm s^{-1}$]')
ax.tick_params(axis='both', direction="in")
ax.tick_params(axis='both', which='minor', direction="in")
ax.set_xlim([0.6001, 1.7999])
ax.set_ylim([0.001, 11.9999])

if fslope != -1:
    name = f'ronney_flamespeed_1.0Alpha_'+date+f' (slope={fslope} curve={fcurve})'+f'_{args.paper}'
else:
    name = f'ronney_flamespeed_1.0Alpha_'+date+f'_{args.paper}'
    
if save_plots == True:
    plt.savefig("C:\\Users\\pjsin\\Documents\\cantera\\burkelab_SimScripts\\figures\\Flame Speed Plots\\"+name+'.pdf', dpi=1000, bbox_inches='tight')
    plt.savefig("C:\\Users\\pjsin\\Documents\\cantera\\burkelab_SimScripts\\figures\\Flame Speed Plots\\"+name+'.png', dpi=dpi, bbox_inches='tight')

# plt.show()     