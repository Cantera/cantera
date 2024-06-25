#REACTION 16: H+O2(+M)=HO2(=M)

#%%
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


# %%
#STEP 1: GET RATE CONSTANTS FROM TROE PARAMETERS
def getKfromTroe(P_list,T_list,rxn,ref_collider):
    file = 'C:\\Users\\pjsin\\Documents\\cantera\\test\\data\\sandbox.yaml'
    reaction = rxn
    gas = ct.Solution(file)
    def getK(Temp,Pres,X) :
        gas.TPX = Temp,Pres,X
        k = gas.forward_rate_constants[gas.reaction_equations().index(reaction)]
        return k
    plt.figure()
    k_P=[]
    for P in P_list:
        k_T = []
        for T in T_list:
            k_T.append(getK(T,P,{ref_collider:1}))
        k_P.append(k_T)
        plt.plot(T_list,np.log(k_T),label=str(P/101325)+" atm")
    plt.legend(bbox_to_anchor=(1.05, 1.0), loc='upper left')
    plt.title("Troe Rates for "+rxn)
    plt.xlabel("Temperature")
    plt.ylabel("logK [units still undetermined]")
    plt.show()
    return k_P

P_list=[0.0001*101325,0.001*101325,0.01*101325,0.1*101325,1*101325,10*101325,100*101325,1000*101325,10000*101325]
T_list=np.linspace(200,2000,50)
k_P = getKfromTroe(P_list,T_list,'H2O2 (+M) <=> 2 OH (+M)','AR')

# %% STEP 2: PERFORM A PLOG FIT FOR AR COLLIDER
def plogFit(P_list,T_list,k_P,rxn):
    def arrhenius(T, A, n, Ea):
        return np.log(A) + n*np.log(T)+ (-Ea/(1.987*T))
    plt.figure(figsize=(8,6))
    # dataset = pd.read_csv(fname)
    for i,p in enumerate(P_list):
        k_data = k_P[i] #rate constants across a temperature range, for a given P
        popt, pcov = curve_fit(arrhenius, T_list, np.log(k_data),maxfev = 2000)
        print(("- {P: %.3e atm, A: %.5e, b: %.5e, Ea: %.5e}")%(p/101325, popt[0],popt[1],popt[2]))
        lnk_fit = arrhenius(T_list,popt[0],popt[1],popt[2])
        plt.plot(T_list,lnk_fit,label=str(p/101325) + ' fit',linestyle='solid')
        plt.scatter(T_list,lnk_fit,marker='x',label=str(p/101325) + ' data')
    plt.title('PLOG Fit for '+rxn)
    plt.xlabel('Temperature [K]')
    plt.ylabel('ln ( k/[M] [cm^6/molecule^2*s] )')
    plt.legend(bbox_to_anchor=(1.05, 1.0), loc='upper left')
    plt.show()
plogFit(P_list,T_list,k_P,'H2O2 (+M) <=> 2 OH (+M)')

# %% STEP 3: COMPUTE ARRHENIUS PARAMETERS FOR COLLIDERS USING COLLISION EFFICIENCIES FROM JASPER
from scipy.optimize import least_squares
def plot_ratefit(temperatures, rate_constants, pltcolour, labell): #3-parameter fit (A, b, Ea)
    def arrhenius_rate(T, A, beta, Ea):
        # R = 8.314  # Gas constant in J/(mol K)
        R = 1.987 # cal/molK
        return A * T**beta * np.exp(-Ea / (R * T))
    def fit_function(params, T, ln_rate_constants):
        A, beta, Ea = params
        return np.log(arrhenius_rate(T, A, beta, Ea)) - ln_rate_constants
    initial_guess = [3, 0.5, 50.0]  
    result = least_squares(fit_function, initial_guess, args=(temperatures, np.log(rate_constants)))
    A_fit, beta_fit, Ea_fit = result.x
    T_range = np.linspace(200, 2000, 100)
    fit_curve = np.exp(np.log(arrhenius_rate(T_range, A_fit, beta_fit, Ea_fit)))
    print(('- collider: "%s"\n  low-P-rate-constant: {A: %.5e, b: %.5e, Ea: %.5e}')%(labell.upper(),A_fit, beta_fit, Ea_fit))
    plt.plot(T_range, np.log(fit_curve), label=labell, color=pltcolour)
    plt.scatter(temperatures, np.log(rate_constants), label=None, color=pltcolour)

# OBTAIN ARRHENIUS PARAMETERS THAT DESCRIBE THEIR EIG0 FUNCTION (RELATIVE TO EIG0 FOR AR)
pltcolours = ['k', 'g', 'b', 'r', 'olive', 'orange', 'brown', 'indigo']
plt.figure()
plot_ratefit(np.array([300,1000,2000]), np.array([1.50,1.58,1.63]), pltcolours[1], "N2")
plot_ratefit(np.array([300,1000,2000]), np.array([5.22,4.14,3.27]), pltcolours[5], "CO2")
plot_ratefit(np.array([300,1000,2000]), np.array([6.83,12.0,16.3]), pltcolours[0], "H2O2") 
plot_ratefit(np.array([300,1000,2000]), np.array([5.51,10.2,13.3]), pltcolours[2], "H2O") 
plt.xlabel('Temperature (K)')
plt.ylabel('Rate Constant')
plt.legend()
plt.title("Arrhenius fits for R25: H2O2 (+M) <=> OH + OH (+M)")
plt.show()

# %%
# SPECIAL CASE 2: Assigning H2O parameters to ALL species
speciesList=['NH3', 'NO2', 'NO', 'N2O', 'H2', 'O2', 'O3', 'H', 'O', 'OH', 'HO2', 'H2O', 'H2O2', 'CO',
    'CO2', 'HOCO', 'CH2O', 'HCO', 'NH2', 'NH', 'N', 'NNH', 'N2H4', 'N2H3', 'tHNNH', 'cHNNH', 'H2NN',
    'NH2OH', 'H2NO', 'HNOH', 'HNO', 'HON', 'HONO', 'HNO2', 'NO3', 'HONO2', 'H2NCO', 'HNCO', 'NCO',
    'AR', 'HE', 'N2']
#H2O2 (+M) <=> OH + OH (+M)
for species in speciesList:
    print(("- collider: '%s'\n  low-P-rate-constant: {A: 1.36377e+00, b: 3.06592e-01, Ea: 2.10079e+02}")%(species))
# %%
