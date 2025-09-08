# This file is part of Cantera. See License.txt in the top-level directory or
# at https://cantera.org/license.txt for license and copyright information.

#cython: language_level=3
#distutils: language=c++

from .constants cimport *

#: Avogadro's constant, 1/kmol
avogadro = CxxAvogadro

#: The ideal gas constant in J/kmol/K
gas_constant = CxxGasConstant

#: One atmosphere in Pascals
one_atm = CxxOneAtm

#: Boltzmann constant
boltzmann = CxxBoltzmann

#: Planck constant (J/s)
planck = CxxPlanck

#: The Stefan-Boltzmann constant, W/m²/K⁴
stefan_boltzmann = CxxStefanBoltz

#: The charge on an electron (C)
electron_charge = CxxElectronCharge

#: The mass of an electron (kg)
electron_mass = CxxElectronMass

#: Faraday constant, C/kmol
faraday = CxxFaraday

#: Speed of Light (m/s).
light_speed = CxxLightSpeed

#: Permeability of free space :math:`\mu_0` in N/A².
permeability_0 = CxxPermeability_0

#: Permittivity of free space (Farads/m = C²/N/m²)
epsilon_0 = CxxEpsilon_0
