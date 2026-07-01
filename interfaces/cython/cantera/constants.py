# This file is part of Cantera. See License.txt in the top-level directory or
# at https://cantera.org/license.txt for license and copyright information.

# distutils: language = c++
# cython: language_level=3
# pyright: reportUndefinedVariable=false

#: Avogadro's constant, 1/kmol
avogadro: float = CxxAvogadro

#: The ideal gas constant in J/kmol/K
gas_constant: float = CxxGasConstant

#: One atmosphere in Pascals
one_atm: float = CxxOneAtm

#: Boltzmann constant
boltzmann: float = CxxBoltzmann

#: Planck constant (J/s)
planck: float = CxxPlanck

#: The Stefan-Boltzmann constant, W/m²/K⁴
stefan_boltzmann: float = CxxStefanBoltz

#: The charge on an electron (C)
electron_charge: float = CxxElectronCharge

#: The mass of an electron (kg)
electron_mass: float = CxxElectronMass

#: Faraday constant, C/kmol
faraday: float = CxxFaraday

#: Speed of Light (m/s).
light_speed: float = CxxLightSpeed

#: Permeability of free space :math:`\mu_0` in N/A².
permeability_0: float = CxxPermeability_0

#: Permittivity of free space (Farads/m = C²/N/m²)
epsilon_0: float = CxxEpsilon_0
