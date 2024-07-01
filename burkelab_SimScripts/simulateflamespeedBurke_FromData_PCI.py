
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
import sys, os
import matplotlib.pyplot as plt
import pandas as pd 
import numpy as np
import matplotlib as mpl

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
parser.add_argument('--title', type=str, help="title = ",default='null')


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

import matplotlib.ticker as ticker
ax.xaxis.set_major_locator(ticker.MultipleLocator(4))
ax.xaxis.set_major_formatter(ticker.StrMethodFormatter("{x:.0f}"))
ax.yaxis.set_major_locator(ticker.MultipleLocator(0.03))
ax.yaxis.set_major_formatter(ticker.StrMethodFormatter("{x:.2f}"))

lw=args.lw
mw=args.mw
msz=args.msz
dpi=args.dpi
lgdw=args.lgdw
lgdfsz=args.lgdfsz
date=args.date
fslope=args.slopeVal
fcurve=args.curveVal

if fslope != -1:
    path="C:\\Users\\pjsin\\Documents\\cantera\\burkelab_SimScripts\\BurkeSongResults_"+date+f' (slope={fslope} curve={fcurve})\\'
else:
    path="C:\\Users\\pjsin\\Documents\\cantera\\burkelab_SimScripts\\BurkeSongResults_"+date+'\\'

# dataset=pd.read_csv(path+f'Alzueta_data.csv')
# ax.plot(dataset.iloc[:,0],dataset.iloc[:,1],linewidth=lw,color="xkcd:grey",label='Alzueta',zorder=80)

# dataset=pd.read_csv(path+f'Alzueta-300K_data.csv')
# ax.plot(dataset.iloc[:,0],dataset.iloc[:,1],linewidth=lw,color="cyan",label='$\epsilon_{NH_3}(300K)$',zorder=60)

# dataset=pd.read_csv(path+f'Alzueta-2000K_data.csv')
# ax.plot(dataset.iloc[:,0],dataset.iloc[:,1],linewidth=lw,color="xkcd:teal",label='$\epsilon_{NH_3}(2000K)$',zorder=70)

dataset=pd.read_csv(path+f'Ar_data.csv')
ax.plot(dataset.iloc[:,0],dataset.iloc[:,1],linewidth=lw,color='r',label='Ar',zorder=20)

dataset=pd.read_csv(path+f'H2O_data.csv')
ax.plot(dataset.iloc[:,0],dataset.iloc[:,1],linewidth=lw,color='b',label=r'H$_2$O',zorder=300)

dataset=pd.read_csv(path+f'LMR-R_data.csv')
ax.plot(dataset.iloc[:,0],dataset.iloc[:,1],linewidth=lw,color='xkcd:purple',label='LMR-R',zorder=100)

# if args.title != 'null':
#     ax.set_title(args.title)

path="G:\\Mon disque\\Columbia\\Burke Lab\\01 Mixture Rules Project\\Graph Reading\\"

# dataset = pd.read_csv(path+'\\5 FS H2O (Burke)\\black.csv')
# pressures = dataset.iloc[:,0]
# mbr_list = dataset.iloc[:,1]
# ax.plot(pressures,mbr_list,marker='^',fillstyle='none',markersize=1,markeredgewidth=mw,linestyle='none',color='b',label='Burke Model',zorder=1)
dataset = pd.read_csv(path+'\\5 FS H2O (Burke)\\exp_pts.csv',header=None)
pressures = dataset.iloc[:,0]
mbr_list = dataset.iloc[:,1]
ax.plot(pressures,mbr_list,marker='o',fillstyle='none',markersize=msz,markeredgewidth=mw,linestyle='none',color='k',label='Burke/Song',zorder=100)
ax.legend(fontsize=lgdfsz, frameon=False, loc='right', handlelength=lgdw) 

ax.set_ylabel(r'Mass burning rate [g $\rm cm^{-2}$ $\rm s^{-1}$]')
ax.set_xlabel(r'Pressure [atm]')
ax.tick_params(axis='both', direction="in")
ax.tick_params(axis='both', which='minor', direction="in")
ax.set_xlim([0.001, 19.999])
ax.set_ylim([-0.005, 0.1299])


if fslope != -1:
    name = f"burkesong_flamespeed_"+date+f' (slope={fslope} curve={fcurve})_PCI'
else:
    name = f"burkesong_flamespeed_"+date+"_PCI"

path="C:\\Users\\pjsin\\Documents\\cantera\\burkelab_SimScripts\\figures\\Flame Speed Plots\\"
if save_plots == True:
    plt.savefig(path+name+'.pdf', dpi=500, bbox_inches='tight')
    plt.savefig(path+name+'.svg', dpi=500, bbox_inches='tight')

# plt.show()     