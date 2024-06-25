#Validation script
import sys, os
sys.path.append("C:/Users/pjsin/Documents/cantera/build/python")
import cantera as bklabct
import numpy as np
import matplotlib.pyplot as plt
bklabct.print_stack_trace_on_segfault()

file = 'C:\\Users\\pjsin\\Documents\\cantera\\test\\data\\kineticsfromscratch_LMRtest.yaml'
reactions = ['H + O2 (+M) <=> HO2 (+M)']
# reactions = ['HO2 (+M) <=> H + O2 (+M)']
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
            gas.TPX = T,P,{'Ar':0.5,'H2O':0.5}
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

