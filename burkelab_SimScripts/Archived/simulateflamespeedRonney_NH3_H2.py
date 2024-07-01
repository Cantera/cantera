# To run file: 
# python "C:\\Users\\pjsin\\Documents\\cantera\\burkelab_SimScripts\\simulateflamespeedRonney_NH3_H2.py"

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
import csv
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('--gridsz', type=int, help="gridsz = ", default=10)
parser.add_argument('--date', type=str, help="sim date = ",default='May28')
parser.add_argument('--slopeVal', type=float, help="slope value = ")
parser.add_argument('--curveVal', type=float, help="curve value = ")
parser.add_argument('--transport', type=str, help="flame transport = ")
args = parser.parse_args()

############# CHANGE THESE ####################################################################################
gridsz = args.gridsz
date=args.date

# fuel_list = np.linspace(0.14,0.5,gridsz) #fuel mole fractions
fuel_list = np.linspace(0.14,0.4,gridsz)
# alpha_list = [1.0,0.8,0.6,0.4,0.2,0.0]
# a_st = [0.75,0.7,0.65,0.6,0.55,0.5]
alpha_list = [1.0,0.6]
a_st = [0.75,0.65]
# alpha_list = [1]
# a_st = [0.75]
fratio=3
fslope=args.slopeVal #should be low enough that the results don't depend on the value. 0.02 for both is a good place to start. Try 0.01 and 0.05 and see if there are any differences
fcurve=args.curveVal
ftransport=args.transport # 'multicomponent' or 'mixture-averaged'
models = {    
          'Alzueta':"C:\\Users\\pjsin\\Documents\\cantera\\test\\data\\alzuetamechanism.yaml",  
          'Alzueta-300K':"C:\\Users\\pjsin\\Documents\\cantera\\test\\data\\alzuetamechanism_epsNH3_T=300K.yaml",  
          'Alzueta-2000K':"C:\\Users\\pjsin\\Documents\\cantera\\test\\data\\alzuetamechanism_epsNH3_T=2000K.yaml",            
          'Ar':"C:\\Users\\pjsin\\Documents\\cantera\\test\\data\\alzuetamechanism_LMRR_allAR.yaml",
          r'H$_2$O':"C:\\Users\\pjsin\\Documents\\cantera\\test\\data\\alzuetamechanism_LMRR_allH2O.yaml",
          'LMR-R':"C:\\Users\\pjsin\\Documents\\cantera\\test\\data\\alzuetamechanism_LMRR.yaml", 
          }
###############################################################################################################

# p_list = [50,100,250,760,1500]
p_list=[760]
Tin = 300.0  # unburned gas temperature [K]
width = 0.03  # m
loglevel = 1  # amount of diagnostic output (0 to 8)

def save_to_csv(filename, data):
    with open(filename, 'w', newline='') as csvfile:
        writer = csv.writer(csvfile)
        writer.writerows(data)

for x, alpha in enumerate(alpha_list):
    for k, m in enumerate(models):
        for i, p in enumerate(p_list):
            mbr = []
            phi_list = []
            for j, fuel_frac in enumerate(fuel_list):
                gas = ct.Solution(list(models.values())[k])
                NH3 = alpha*fuel_frac
                H2 = (1-alpha)*fuel_frac
                
                ox_frac = 1 - fuel_frac # oxidizer fraction
                O2 = ox_frac*0.21
                N2 = ox_frac*0.79
                phi = np.divide(fuel_frac/O2,1/a_st[x]) # THIS STEP MUST BE REVIEWED
                phi_list.append(phi)
                X = {'NH3':NH3,'H2':H2,'O2':O2,'N2':N2}
                gas.TPX = Tin, (p/760)*ct.one_atm, X
                f = ct.FreeFlame(gas, width=width)
                f.set_refine_criteria(ratio=fratio, slope=fslope, curve=fcurve)
                f.transport_model = ftransport
                f.soret_enabled = True
                f.solve(loglevel=loglevel, auto=True)
                mbr.append(f.velocity[0] * 100) # cm/s

            # Save phi_list and mbr to CSV
            path=f'C:\\Users\\pjsin\\Documents\\cantera\\burkelab_SimScripts\\RonneyResults_'+date+f' (slope={fslope} curve={fcurve})'
            os.makedirs(path,exist_ok=True)
            csv_filename =path+f'\\{m}_{i}_data_{alpha}alpha.csv'
            data = zip(phi_list, mbr)
            save_to_csv(csv_filename, data)