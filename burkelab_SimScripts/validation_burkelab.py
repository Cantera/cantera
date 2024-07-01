
#Validation script
import sys, os
sys.path.append("C:/Users/pjsin/Documents/cantera/build/python")
import cantera as bklabct
import numpy as np
import matplotlib.pyplot as plt
# bklabct.print_stack_trace_on_segfault()

file = 'C:\\Users\\pjsin\\Documents\\cantera\\test\\data\\kineticsfromscratch_LMRtest.yaml'
reactions = ['H + O2 (+M) <=> HO2 (+M)']
gas = bklabct.Solution(file)
#Temp = np.linspace(750,2500,50)
Temp=[1000]
# Pres = np.logspace(-2,2,5)
Pres = [101325] # units: Pa
# Pres=[1e-6,1e-5,1e-4,0.001,0.01,0.1, 1, 10, 100, 1000, 10000, 100000, 1e6, 1e7, 1e8]
for i, R in enumerate(reactions):
    k_list=[]
    for j, P in enumerate(Pres):
        temp_list = []
        for k,T in enumerate(Temp):
            gas.TPX = T,P,{'Ar':0,'H2O':1}
            # gas.TPX = T,P,{'H2O':1.0}
            # gas.TPX = T,P,{'Ar':1.0}
            # print(gas.TPX)
            rc = gas.forward_rate_constants[gas.reaction_equations().index(R)]
            temp_list.append(rc)
        k_list.append(temp_list)  
        print("%.5e" % k_list[j][0])

    # plt.figure()
    # plt.title(R)
    # for j,P in enumerate(Pres):
    #     plt.plot(Temp,k_list[j],label=str(P)+' atm') 
    #     #plt.semilogy(Temp,k_list[j],label=str(P)+' atm')    
    # plt.legend()
    # plt.xlabel('Temperature [K]')
    # plt.ylabel('k')
    
    # plt.savefig("reaction1",bbox_inches="tight")
    # plt.show()
    
# 1.0810363605840141e-13, 2.1620727211680282e-13, 3.243109081752042e-13, 4.324145442336056e-13,

# print()

# # %%
# T=1000 #K
# P=101325 #Pa
# def k (dict, T):
#     if 'eig0' in dict:
#         A = dict['eig0']['A']
#         b = dict['eig0']['b']
#         Ea = dict['eig0']['Ea']
#     if 'plog' in dict:
#         A = dict['plog']['A']
#         b = dict['plog']['b']
#         Ea = dict['plog']['Ea']
#     R = 1.987 # cal/molK
#     return A*T**(b)*np.exp(-Ea/R/T)
# AR = {'X':0.5, 'eig0':{'A': 2.20621e-02, 'b': 4.74036e-01, 'Ea': -1.13148e+02}}
# H2O = {'X':0.5, 'eig0':{'A': 1.04529e-01, 'b': 5.50787e-01, 'Ea': -2.32675e+02}}

# eig0_mix = AR['X'] * k(AR,T) + H2O['X'] * k(H2O,T)

# def Peff(eig0,eig0_mix,P,species):
#     Peff = np.exp(np.log(P)+np.log(eig0_mix)-np.log(eig0))/101325
#     label = "Peff_"+species+" [atm] = "
#     print(label, Peff)
#     return Peff
# Peff_AR = Peff(k(AR,T),eig0_mix,P,"AR")
# Peff_H2O = Peff(k(H2O,T),eig0_mix,P,"H2O")
# #%%

# def interpolate(p, p1, p2, y1, y2):
#     """Interpolates to find y at a specific p using points (p1, y1) and (p2, y2)."""
#     return y1 + (y2 - y1) * ((p - p1) / (p2 - p1))

# # Given data points
# data = {
#     1.000e-04: {"A": 5.30514e+12, "b": -2.80725, "Ea": 499.267},
#     1.000e-03: {"A": 5.25581e+13, "b": -2.80630, "Ea": 499.946}
# }

# # Target pressure
# p_target = 10**(-3.5)

# # Interpolate A, b, Ea
# A_interpolated = interpolate(p_target, 1.000e-04, 1.000e-03, data[1.000e-04]["A"], data[1.000e-03]["A"])
# b_interpolated = interpolate(p_target, 1.000e-04, 1.000e-03, data[1.000e-04]["b"], data[1.000e-03]["b"])
# Ea_interpolated = interpolate(p_target, 1.000e-04, 1.000e-03, data[1.000e-04]["Ea"], data[1.000e-03]["Ea"])


# AR = {'X':0.5, 'plog':{'A': 2.20621e-02, 'b': 4.74036e-01, 'Ea': -1.13148e+02}}
# H2O = {'X':0.5, 'plog':{'A': 1.04529e-01, 'b': 5.50787e-01, 'Ea': -2.32675e+02}}
