"""Conversion factors to SI (m, kg, kmol, s)"""

from constants import Avogadro, GasConstant

kmol = 1.0
mol  = 1.e-3
molecule = kmol/Avogadro

m    = 1.0
cm   = 0.01
mm   = 0.001

m2   = 1.0
cm2  = 1.e-4
mm2  = 1.e-6
A2   = 1.e-20

m3   = 1.0
cm3  = 1.e-6
mm3  = 1.e-9

J    = 1.0
kJ   = 1000.0
cal  = 4.184
kcal = 4184.0

K    = 1.0

kJ_per_mol  = kJ/mol
cal_per_mol = cal/mol
kcal_per_mol = kcal/mol

mol_per_cm2 = mol/cm/cm
molecule_per_cm2 = molecule/cm/cm

# pressure
Pa   = 1.0
kPa  = 1000.0
atm  = 1.01325e5
bar  = 1.0e5
torr = atm/760.0

# mass
kg   = 1.0
gm   = 1000.0

# mass flux
kg_per_m2_per_s  = 1.0
g_per_cm2_per_s  = gm/(cm*cm)

_lengthdict = {'m':m, 'cm':cm, 'mm':mm}
def length(u):
    return _lengthdict[u]

_moldict = {'kmol':kmol, 'mol':mol, 'molecule':molecule}
def mole(u):
    return _moldict[u]

_eadict = {'kJ_per_mol':kJ_per_mol,
           'kcal_per_mol':kcal_per_mol,
           'cal_per_mol':cal_per_mol,
           'K':GasConstant}

def actEnergy(u):
    return _eadict[u]/GasConstant
