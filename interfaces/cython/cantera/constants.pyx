# This file is part of Cantera. See License.txt in the top-level directory or
# at https://cantera.org/license.txt for license and copyright information.

#cython: language_level=3
#distutils: language=c++

from .constants cimport *

#: Avogadro's Number, /kmol
avogadro: float = CxxAvogadro

#: The ideal gas constant in J/kmol-K
gas_constant: float = CxxGasConstant

#: One atmosphere in Pascals
one_atm: float = CxxOneAtm

#: Boltzmann constant
boltzmann: float = CxxBoltzmann

#: Planck constant (J/s)
planck: float = CxxPlanck

#: The Stefan-Boltzmann constant, W/m^2K^4
stefan_boltzmann: float = CxxStefanBoltz

#: The charge on an electron (C)
electron_charge: float = CxxElectronCharge

#: The mass of an electron (kg)
electron_mass: float = CxxElectronMass

#: Faraday constant, C/kmol
faraday: float = CxxFaraday

#: Speed of Light (m/s).
light_speed: float = CxxLightSpeed

#: Permeability of free space :math:`\mu_0` in N/A^2.
permeability_0: float = CxxPermeability_0

#: Permittivity of free space (Farads/m = C^2/N/m^2)
epsilon_0: float = CxxEpsilon_0
