
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
fig, ax = plt.subplots(2,3,figsize=(args.figwidth, args.figheight))
# plt.subplots_adjust(wspace=0.4, hspace=1)
import matplotlib.ticker as ticker


if args.title != 'null':
    fig.text(0.5, 1, args.title, ha='center', va='center')

lw=args.lw
mw=args.mw
msz=args.msz
dpi=args.dpi
lgdw=args.lgdw
lgdfsz=args.lgdfsz
date=args.date
fslope=args.slopeVal
fcurve=args.curveVal

alpha_list=[1.0,0.8,0.6,0.4,0.2,0.0]
idxs=[(0,0),(0,1),(0,2),(1,0),(1,1),(1,2)]
# alpha_list=[1]
# idxs=[(0,0)]
# pct=["1pt0","0pt8","0pt6","0pt4","0pt2"]
plt.subplots_adjust(wspace=0.3,hspace=0.3)
for x, alpha in enumerate(alpha_list):
    ax[idxs[x]].xaxis.set_major_locator(ticker.MultipleLocator(0.2))
    ax[idxs[x]].xaxis.set_major_formatter(ticker.StrMethodFormatter("{x:.1f}"))
    ax[idxs[x]].yaxis.set_major_formatter(ticker.StrMethodFormatter("{x:.0f}"))

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

    if fslope != -1:
        path="C:\\Users\\pjsin\\Documents\\cantera\\burkelab_SimScripts\\RonneyResults_"+date+f' (slope={fslope} curve={fcurve})\\'
    else:
        path="C:\\Users\\pjsin\\Documents\\cantera\\burkelab_SimScripts\\RonneyResults_"+date+"\\"
    
    dataset=pd.read_csv(path+f'Alzueta_0_data_{alpha}alpha.csv')
    ax[idxs[x]].plot(dataset.iloc[:,0],dataset.iloc[:,1],linewidth=lw,color="xkcd:grey",label='Alzueta')

    dataset=pd.read_csv(path+f'Ar_0_data_{alpha}alpha.csv')
    ax[idxs[x]].plot(dataset.iloc[:,0],dataset.iloc[:,1],linewidth=lw,color='r',label='Ar')

    dataset=pd.read_csv(path+f'H2O_0_data_{alpha}alpha.csv')
    ax[idxs[x]].plot(dataset.iloc[:,0],dataset.iloc[:,1],linewidth=lw,color='b',label=r'H$_2$O')

    dataset=pd.read_csv(path+f'LMR-R_0_data_{alpha}alpha.csv')
    ax[idxs[x]].plot(dataset.iloc[:,0],dataset.iloc[:,1],linewidth=lw,color='xkcd:purple',label='LMR-R')

    ax[idxs[x]].set_title(f'{round(alpha*100)}% NH3/{round((1-alpha)*100)}% H2')

    path="G:\\Mon disque\\Columbia\\Burke Lab\\01 Mixture Rules Project\\Graph Reading\\"

    if x==0:
        dataset = pd.read_csv(path+'\\6 FS NH3 (Stagni-Ronney)\\760torr.csv')
        NH3_list = np.divide(dataset.iloc[:,0],100)
        ox_frac_list = np.subtract(1,NH3_list)
        O2_list = np.multiply(ox_frac_list, 0.21)
        phi_list = np.divide(np.divide(NH3_list,O2_list),np.divide(4,3))
        ax[idxs[x]].plot(phi_list,dataset.iloc[:,1],marker='o',fillstyle='none',markersize=msz,markeredgewidth=0.5,linestyle='none',color='k',label='Ronney')

        dataset = pd.read_csv(path+f'\\AlzuetaFig15\\1pt0_NH3.csv')
        ax[idxs[x]].plot(dataset.iloc[:,0],dataset.iloc[:,1],marker='s',markersize=msz,markeredgewidth=0.5,linestyle='none',color='b',label='Alz Sim (Graph Read)')
        ax[idxs[x]].legend(fontsize=lgdfsz, frameon=False, loc='lower right') 
    
    if x==2:
        dataset = pd.read_csv(path+f'\\Han\\han_0pt6_NH3.csv')
        ax[idxs[x]].plot(dataset.iloc[:,0],dataset.iloc[:,1],marker='s',fillstyle='none',markersize=msz,markeredgewidth=0.5,linestyle='none',color='k',label='Han')
        dataset = pd.read_csv(path+f'\\Wang\\wang_0pt6_NH3.csv')
        ax[idxs[x]].plot(dataset.iloc[:,0],dataset.iloc[:,1],marker='x',fillstyle='none',markersize=msz,markeredgewidth=0.5,linestyle='none',color='k',label='Wang')
        ax[idxs[x]].legend(fontsize=lgdfsz, frameon=False, loc='lower right') 

    if x==3:
        dataset = pd.read_csv(path+f'\\Wang\\wang_0pt4_NH3.csv')
        ax[idxs[x]].plot(dataset.iloc[:,0],dataset.iloc[:,1],marker='x',fillstyle='none',markersize=msz,markeredgewidth=0.5,linestyle='none',color='k',label='Wang')
        ax[idxs[x]].legend(fontsize=lgdfsz, frameon=False, loc='lower right') 
    if x>2:
        ax[idxs[x]].set_xlabel(r'Equivalence Ratio')
    if x==0 or x==3:
        ax[idxs[x]].set_ylabel(r'Burning velocity [cm $\rm s^{-1}$]')
    ax[idxs[x]].tick_params(axis='both', direction="in")
    ax[idxs[x]].tick_params(axis='both', which='minor', direction="in")
    ax[idxs[x]].set_xlim([0.6001,1.7999])

if fslope != -1:
    name = f'ronney_flamespeed_allAlpha_'+date+f' (slope={fslope} curve={fcurve})'
else:
    name = f'ronney_flamespeed_allAlpha_'+date
    
if save_plots == True:
    plt.savefig("C:\\Users\\pjsin\\Documents\\cantera\\burkelab_SimScripts\\figures\\Flame Speed Plots\\"+name+'.pdf', dpi=1000, bbox_inches='tight')
    plt.savefig("C:\\Users\\pjsin\\Documents\\cantera\\burkelab_SimScripts\\figures\\Flame Speed Plots\\"+name+'.png', dpi=dpi, bbox_inches='tight')

# plt.show()     