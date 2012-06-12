"""
Physical Constants

These values are the same as those in the C++ header file ct_defs.h in
the Cantera kernel.
"""

import math

#: One atmosphere in Pascals
OneAtm = 101325.0

#: The ideal gas constant in J/kmo-K
GasConstant = 8314.4621

#: Avogadro's Number, /kmol
Avogadro = 6.02214129e26

#: The ideal gas constant in cal/mol-K
GasConst_cal_mol_K = GasConstant / 4184.0

#: Boltzmann-s constant
Boltzmann = GasConstant / Avogadro

#: The Stefan-Boltzmann constant, W/m^2K^4
StefanBoltz = 5.670373e-8

#: The charge on an electron (C)
ElectronCharge = 1.602176565e-19

#: The mass of an electron (kg)
ElectronMass = 9.10938291e-31

Pi = math.pi

#: Faraday's constant, C/kmol
Faraday = ElectronCharge * Avogadro

#: Planck's constant (J/s)
Planck = 6.62607009e-34

#: Speed of Light (m/s).
lightSpeed = 299792458.0

#: Permeability of free space :math:`\mu_0` in N/A^2.
permeability_0 = 4.0e-7*Pi ## N/A^2

#: Permittivity of free space
epsilon_0 = 1.0 / (lightSpeed*lightSpeed*permeability_0) ## Farads/m = C^2/N/m^2
