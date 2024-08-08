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

# PLOG table already given w.r.t. Ar, so no need to calculate one. Skip to Step 3.

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
plt.figure() #using Jasper efficiencies
plot_ratefit(np.array([300, 1000, 2000]), np.array([3.15,2.47,2.25]), pltcolours[1], "N2")
plot_ratefit(np.array([300, 1000, 2000]), np.array([1.55,1.48,1.70]), pltcolours[2], "O2")
plot_ratefit(np.array([300, 1000, 2000]), np.array([11.2,13.4,14.3]), pltcolours[3], "CO2")
plot_ratefit(np.array([300, 1000, 2000]), np.array([13.9,20.0,22.2]), pltcolours[4], "NH3")
plot_ratefit(np.array([300, 1000, 2000]), np.array([9.94,13.3,14.3]), pltcolours[0], "CH4")
plot_ratefit(np.array([300, 1000, 2000]), np.array([14.0,23.6,27.9]), pltcolours[5], "H2O")
plt.xlabel('Temperature (K)')
plt.ylabel('Rate Constant')
plt.legend()
plt.title("Arrhenius fits for R72: NH3 (+M) <=> H + NH2 (+M)")
plt.show()

# %%
# SPECIAL CASE 2: Assigning H2O parameters to ALL species
speciesList=['NH3', 'NO2', 'NO', 'N2O', 'H2', 'O2', 'O3', 'H', 'O', 'OH', 'HO2', 'H2O', 'H2O2', 'CO',
    'CO2', 'HOCO', 'CH2O', 'HCO', 'NH2', 'NH', 'N', 'NNH', 'N2H4', 'N2H3', 'tHNNH', 'cHNNH', 'H2NN',
    'NH2OH', 'H2NO', 'HNOH', 'HNO', 'HON', 'HONO', 'HNO2', 'NO3', 'HONO2', 'H2NCO', 'HNCO', 'NCO',
    'AR', 'HE', 'N2']
#NH3 (+M) <=> H + NH2 (+M)
for species in speciesList:
    print(("- collider: '%s'\n  low-P-rate-constant: {A: 1.14560e+01, b: 1.27501e-01, Ea: 3.13959e+02}")%(species))

# %%
