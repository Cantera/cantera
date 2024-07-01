#%%
#Validation script
import sys, os
sys.path.append("C:/Users/pjsin/Documents/cantera/build/python")
import matplotlib.pyplot as plt
import cantera as bklabct
import numpy as np

# bklabct.print_stack_trace_on_segfault()
#%%
# def printRateConstant(Temp,Pres,X) :
#     # file = 'test/data/kineticsfromscratch_LMRtest.yaml'
#     # file = 'test/data/alzuetamechanism_LMRR.yaml'
#     file = 'test/data/sandbox.yaml'
#     # reactions = ['H + O2 (+M) <=> HO2 (+M)']
#     reactions = ['NH2 + NH2 (+M) <=> N2H4 (+M)']
#     gas = bklabct.Solution(file)
#     for i, R in enumerate(reactions):
#         k_list=[]
#         for j, P in enumerate(Pres):
#             temp_list = []
#             for k,T in enumerate(Temp):
#                 gas.TPX = T,P,X
#                 rc = gas.forward_rate_constants[gas.reaction_equations().index(R)]
#                 temp_list.append(rc)
#             k_list.append(temp_list)
#             print(("%.5e      %s")%(k_list[j][0],str(X)))


# X_list = [{'AR':1.0,'H2O':0.0},{'AR':0.7,'H2O':0.3},{'AR':0.0,'H2O':1.0},
#           {'AR':1.0,'N2':0.0},{'AR':0.7,'N2':0.3},{'AR':0.0,'N2':1.0},
#           {'AR':1.0,'H2':0.0},{'AR':0.7,'H2':0.3},{'AR':0.0,'H2':1.0},
#           {'AR':1.0,'CO2':0.0},{'AR':0.7,'CO2':0.3},{'AR':0.0,'CO2':1.0},
#           {'AR':1.0,'NH3':0.0},{'AR':0.7,'NH3':0.3},{'AR':0.0,'NH3':1.0},
#           {'AR':1.0,'H2O2':0.0},{'AR':0.7,'H2O2':0.3},{'AR':0.0,'H2O2':1.0},
#           {'HE':1.0,'AR':0.0}]
# for X in X_list:
#     printRateConstant([2000],[2*101325],X)


# def plotK(X,rxn):
#     P=10/760
#     Temp=np.linspace(400,4000)
#     file = 'test/data/alzuetamechanism_LMRR.yaml'
#     gas = bklabct.Solution(file)
#     k_list=[]
#     for k,T in enumerate(Temp):
#         gas.TPX = T,P,X
#         rc = gas.forward_rate_constants[gas.reaction_equations().index(rxn)]
#         k_list.append(rc)
#     plt.plot(Temp,np.log(k_list),label=str(X))

# # reactions=['H2O (+M) <=> H + OH (+M)','H + O2 (+M) <=> HO2 (+M)',
# #            'H2O2 (+M) <=> OH + OH (+M)','NH3 (+M) <=> NH2 + H (+M)',
# #            '2 NH2 (+M) <=> N2H4 (+M)','HNO (+M) <=> H + NO (+M)']
  
# plt.figure()
# rxn='HNO (+M) <=> H + NO (+M)'
# plotK({'AR':1},rxn)
# plotK({'N2':1},rxn)
# plotK({'H2O':1},rxn)
# plt.legend(bbox_to_anchor=(1.05, 1.0), loc='upper left')
# plt.ylabel("logK")
# plt.xlabel("T [K]")
# plt.title("Collider eig0 vals for "+rxn)
# plt.show()

# plt.figure()
# rxn='2 NH2 (+M) <=> N2H4 (+M)'
# plotK({'AR':1},rxn)
# plotK({'N2':1},rxn)
# plotK({'O2':1},rxn)
# plotK({'NH3':1},rxn)
# plotK({'H2O':1},rxn)
# plt.legend(bbox_to_anchor=(1.05, 1.0), loc='upper left')
# plt.ylabel("logK")
# plt.xlabel("T [K]")
# plt.title("Collider eig0 vals for "+rxn)
# plt.show()

# plt.figure()
# rxn='NH3 (+M) <=> H + NH2 (+M)'
# plotK({'AR':1},rxn)
# plotK({'N2':1},rxn)
# plotK({'O2':1},rxn)
# plotK({'CO2':1},rxn)
# plotK({'NH3':1},rxn)
# # plotK({'CH4':1},rxn) #CH4 commmented-out in yaml bc not in species list
# plotK({'H2O':1},rxn)
# plt.legend(bbox_to_anchor=(1.05, 1.0), loc='upper left')
# plt.ylabel("logK")
# plt.xlabel("T [K]")
# plt.title("Collider eig0 vals for "+rxn)
# plt.show()

# plt.figure()
# rxn='H2O2 (+M) <=> 2 OH (+M)'
# plotK({'AR':1},rxn)
# plotK({'N2':1},rxn)
# plotK({'CO2':1},rxn)
# plotK({'H2O2':1},rxn)
# plotK({'H2O':1},rxn)
# plt.legend(bbox_to_anchor=(1.05, 1.0), loc='upper left')
# plt.ylabel("logK")
# plt.xlabel("T [K]")
# plt.title("Collider eig0 vals for "+rxn)
# plt.show()

# plt.figure()
# rxn='H + O2 (+M) <=> HO2 (+M)'
# plotK({'AR':1},rxn)
# plotK({'HE':1},rxn)
# plotK({'N2':1},rxn)
# plotK({'H2':1},rxn)
# plotK({'CO2':1},rxn)
# plotK({'NH3':1},rxn)
# plotK({'H2O':1},rxn)
# plt.legend(bbox_to_anchor=(1.05, 1.0), loc='upper left')
# plt.ylabel("logK")
# plt.xlabel("T [K]")
# plt.title("Collider eig0 vals for "+rxn)
# plt.show()

# plt.figure()
# rxn='H2O (+M) <=> H + OH (+M)'
# plotK({'N2':1},rxn)
# plotK({'AR':1},rxn)
# plotK({'H2O':1},rxn)
# plt.legend(bbox_to_anchor=(1.05, 1.0), loc='upper left')
# plt.ylabel("logK")
# plt.xlabel("T [K]")
# plt.title("Collider eig0 vals for "+rxn)
# plt.show()

# # 'H2O (+M) <=> H + OH (+M)' #reaction 11
# # 'H + O2 (+M) <=> HO2 (+M)' #reaction 13
# # 'H2O2 (+M) <=> OH + OH (+M)' #reaction 22
# # 'NH3 (+M) <=> NH2 + H (+M)' #reaction 69
# # '2 NH2 (+M) <=> N2H4 (+M)' #reaction 87
# # 'HNO (+M) <=> H + NO (+M)' #reaction 164
# # %%


def getK(T,P,X,rxn):
    file = 'test/data/alzuetamechanism_LMRR.yaml'
    gas = bklabct.Solution(file)
    gas.TPX = T,P,X
    rc = gas.forward_rate_constants[gas.reaction_equations().index(rxn)]
    print(rc)
    # print(gas.species_names)
# rxn1='H2O2 (+M) <=> 2 OH (+M)'
# rxn2='H2O2 (+M) <=> H + HO2 (+M)'
# print("Test 1")
# getK(1000,1e-3*101325,{"H2O":1}, rxn1)
# getK(1000,1e-3*101325,{"H2O":1},rxn2)
# print("Test 2")
# getK(1000,101325,{"H2O":1}, rxn1)
# getK(1000,101325,{"H2O":1},rxn2)
# print("Test 3")
# getK(1000,1e-1*101325,{"H2O":1}, rxn1)
# getK(1000,1e-1*101325,{"H2O":1},rxn2)
# print("Test 4")
# getK(1000,1e-3*101325,{"H2O":0.5,'AR':0.5}, rxn1)
# getK(1000,1e-3*101325,{"H2O":0.5, 'NH3':0.5},rxn2)
# print("Test 5")
# getK(1000,101325,{"NH3":1}, rxn1)
# getK(1000,101325,{"HE":1},rxn2)
# print("Test 6")
# getK(1000,1e-1*101325,{"AR":1}, rxn1)
# getK(1000,1e-1*101325,{"HE":5,"NH3":0.5},rxn2)


def getK(T,P,X,rxn):
    # file = 'test/data/alzuetamechanism_LMRR.yaml'
    file = 'test/data/kineticsfromscratch_LMRtest_V2.yaml'
    gas = bklabct.Solution(file)
    gas.TPX = T,P,X
    rc = gas.forward_rate_constants[gas.reaction_equations().index(rxn)]
    print(rc)
    # print(gas.species_names)
rxn1='2 OH (+ M) <=> H2O2 (+ M)'
rxn2='H + O2 <=> HO2'
rxn3='HO2 <=> OH + O'
print("Test 1")
getK(1000,1e-3*101325,{"H2O":1}, rxn1)
getK(1000,1e-3*101325,{"H2O":1},rxn2)
getK(1000,1e-3*101325,{"H2O":1},rxn3)
# print("Test 2")
# getK(1000,101325,{"H2O":1}, rxn1)
# getK(1000,101325,{"H2O":1},rxn2)
# print("Test 3")
# getK(1000,1e-1*101325,{"H2O":1}, rxn1)
# getK(1000,1e-1*101325,{"H2O":1},rxn2)
# print("Test 4")
# getK(1000,1e-3*101325,{"H2O":0.5,'AR':0.5}, rxn1)
# getK(1000,1e-3*101325,{"H2O":0.5, 'NH3':0.5},rxn2)
# print("Test 5")
# getK(1000,101325,{"NH3":1}, rxn1)
# getK(1000,101325,{"HE":1},rxn2)
# print("Test 6")
# getK(1000,1e-1*101325,{"AR":1}, rxn1)
# getK(1000,1e-1*101325,{"HE":5,"NH3":0.5},rxn2)



# rxn='HNO (+M) <=> H + NO (+M)'
# print(rxn)
# getK(1000,101325,{'AR':1},rxn) #P=1atm
# getK(1000,101325*1.71,{'N2':1},rxn)
# getK(1000,101325*8.15,{'H2O':1},rxn)

# rxn='2 NH2 (+M) <=> N2H4 (+M)'
# print(rxn)
# getK(1000,101325,{'AR':1},rxn)
# getK(1000,101325*1.69,{'N2':1},rxn)
# getK(1000,101325*1.17,{'O2':1},rxn)
# getK(1000,101325*8.25,{'NH3':1},rxn)
# getK(1000,101325*7.98,{'H2O':1},rxn)

# rxn='NH3 (+M) <=> H + NH2 (+M)'
# print(rxn)
# getK(1000,101325,{'AR':1},rxn)
# getK(1000,101325*2.47,{'N2':1},rxn)
# getK(1000,101325*1.48,{'O2':1},rxn)
# getK(1000,101325*13.4,{'CO2':1},rxn)
# getK(1000,101325*20.0,{'NH3':1},rxn)
# getK(1000,101325*23.6,{'H2O':1},rxn)

# rxn='H2O2 (+M) <=> 2 OH (+M)'
# print(rxn)
# getK(1000,101325,{'AR':1},rxn)
# getK(1000,101325/1.58,{'N2':1},rxn)
# getK(1000,101325/4.14,{'CO2':1},rxn)
# getK(1000,101325/12.0,{'H2O2':1},rxn)
# getK(1000,101325/10.2,{'H2O':1},rxn)

# rxn='H + O2 (+M) <=> HO2 (+M)'
# print(rxn)
# getK(1000,101325,{'AR':1},rxn)
# getK(1000,101325,{'HE':1},rxn)
# getK(1000,101325/1.58,{'N2':1},rxn)
# getK(1000,101325/3.07,{'H2':1},rxn)
# getK(1000,101325/8.94,{'CO2':1},rxn)
# getK(1000,101325/17.9,{'NH3':1},rxn)
# getK(1000,101325/22.2,{'H2O':1},rxn)

# rxn='H2O (+M) <=> H + OH (+M)'
# print(rxn)
# getK(1000,101325,{'N2':1},rxn)
# getK(1000,101325/(1/1.62),{'AR':1},rxn)
# getK(1000,101325/(8.55/1.62),{'H2O':1},rxn)

# %%

gas=bklabct.Solution('test/data/alzuetamechanism_LMRR.yaml')
for i, reaction in enumerate(gas.reactions()):
    print(f'Reaction {i + 1}: {reaction.equation}')


# %%

#H + OH (+M) is a possible candidate. We're having trouble with H, so its prob one of the non ammonia reactions
#Jasper's 3rd body efficiencies work in either direction. But the rate constant data is direction-dependent

# def plotK(fname,Tlist,P,X,rxn,colour,style,label,direction):
#     kvals=[]
#     for T in Tlist: 
#         gas = bklabct.Solution(fname)
#         gas.TPX = T,P*101325,X
#         if direction=='f':
#             rc = gas.forward_rate_constants[gas.reaction_equations().index(rxn)]
#         if direction =='b':
#             rc= gas.reverse_rate_constants[gas.reaction_equations().index(rxn)]
#         kvals.append(rc)
#     plt.semilogy(Tlist,kvals,color=colour,linestyle=style,label=label)

import numpy as np

def plotK(fname,Tlist,P,X,rxn,colour,style,label,direction):
    kvals=[]
    for T in Tlist: 
        gas = bklabct.Solution(fname)
        gas.TPX = T,P*101325,X
        if direction=='f':
            rc = gas.forward_rate_constants[gas.reaction_equations().index(rxn)]
            rc = rc*1000 # cm3/mol*s
        if direction =='b':
            rc= gas.reverse_rate_constants[gas.reaction_equations().index(rxn)] # m6/kmol2*s
            rc = rc*1e6 # cm6/mol2*s
            R = 82.5 # cm3*atm/mol*K
            C = np.divide(P,np.multiply(T,R)) # mol/cm3
            rc = np.multiply(rc,C) # cm3/mol*s


        kvals.append(rc)
    plt.semilogy(Tlist,kvals,color=colour,linestyle=style,label=label)

Tlist=np.linspace(200,2000)
rxns1 = ['H + OH (+M) <=> H2O (+M)', 'H + O2 (+M) <=> HO2 (+M)', 'H2O2 (+M) <=> 2 OH (+M)','NH3 (+M) <=> H + NH2 (+M)', '2 NH2 (+M) <=> N2H4 (+M)', 'HNO (+M) <=> H + NO (+M)']
rxns2 = ['H + OH (+M) <=> H2O (+M)', 'H + O2 (+M) <=> HO2 (+M)', 'H2O2 (+M) <=> 2 OH (+M)','H + NH2 (+M) <=> NH3 (+M)', '2 NH2 (+M) <=> N2H4 (+M)', 'H + NO (+M) <=> HNO (+M)']
# rxns1 = ['H + OH <=> H2O', 'H + O2 <=> HO2']
# rxns2 = ['H + OH (+M) <=> H2O (+M)','H + O2 (+M) <=> HO2 (+M)']
rxns1 = ['H + OH <=> H2O']
rxns2 = ['H2O + M <=> H + OH + M']
dirs1=['f']
dirs2=['b']


# dirs1=['f','f','f','f','f','f']
# dirs2=['f','f','f','b','f','b']

# 'NH3 (+M) <=> H + NH2 (+M)', 
# kmols, m3, kg, s
for i in range(len(rxns1)):
    plt.figure()
    plotK('test/data/sandbox2.yaml',Tlist,1,{'H':1},rxns1[i],'b','solid','NEW: '+rxns1[i],dirs1[i])
    plotK('test/data/alzuetamechanism.yaml',Tlist,1,{'H':1},rxns2[i],'k','solid','OLD: '+rxns2[i],dirs2[i])
    plt.legend()
    plt.show()

# # %%

# def plotK(fname,Tlist,P,X,rxn,colour,style,label,direction,div,divval):
#     kvals=[]
#     for T in Tlist: 
#         gas = bklabct.Solution(fname)
#         gas.TPX = T,P*101325,X
#         if div=='y':
#             if direction=='f':
#                 rc = gas.forward_rate_constants[gas.reaction_equations().index(rxn)]*divval
#             if direction =='b':
#                 rc= gas.reverse_rate_constants[gas.reaction_equations().index(rxn)]*divval
#             kvals.append(rc)
#         else:
#             if direction=='f':
#                 rc = gas.forward_rate_constants[gas.reaction_equations().index(rxn)]
#             if direction =='b':
#                 rc= gas.reverse_rate_constants[gas.reaction_equations().index(rxn)]
#             kvals.append(rc)            
#     plt.semilogy(Tlist,kvals,color=colour,linestyle=style,label=label)

# Tlist=np.linspace(200,2000)
# rxns1 = ['H + OH (+M) <=> H2O (+M)', 'H + O2 (+M) <=> HO2 (+M)']
# rxns2 = ['H2O + M <=> H + OH + M', 'H + O2 (+M) <=> HO2 (+M)']

# dirs1=['f','f']
# dirs2=['b','f']

# divvals=[1000,1/6.02e23*1000]

# # 'NH3 (+M) <=> H + NH2 (+M)', 
# # kmols, m3, kg, s
# for i in range(len(rxns1)):
#     plt.figure()
#     plotK('test/data/alzuetamechanism_LMRR.yaml',Tlist,1e-2,{'H':1},rxns1[i],'b','solid','NEW: '+rxns1[i],dirs1[i],'y',divvals[i])
#     plotK('test/data/alzuetamechanism.yaml',Tlist,1e-2,{'H':1},rxns2[i],'k','solid','OLD: '+rxns2[i],dirs2[i],'n',divvals[i])
#     plt.legend()
#     plt.show()
# # %%

# %%
