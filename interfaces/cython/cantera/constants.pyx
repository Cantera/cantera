# This file is part of Cantera. See License.txt in the top-level directory or
# at https://cantera.org/license.txt for license and copyright information.

#cython: language_level=3
#distutils: language=c++

cdef extern from "cantera/base/ct_defs.h" namespace "Cantera":
    cdef double CxxAvogadro "Cantera::Avogadro"
    cdef double CxxGasConstant "Cantera::GasConstant"
    cdef double CxxOneAtm "Cantera::OneAtm"
    cdef double CxxBoltzmann "Cantera::Boltzmann"
    cdef double CxxPlanck "Cantera::Planck"
    cdef double CxxStefanBoltz "Cantera::StefanBoltz"
    cdef double CxxFaraday "Cantera::Faraday"
    cdef double CxxElectronCharge "Cantera::ElectronCharge"
    cdef double CxxElectronMass "Cantera::ElectronMass"
    cdef double CxxLightSpeed "Cantera::lightSpeed"
    cdef double CxxPermeability_0 "Cantera::permeability_0"
    cdef double CxxEpsilon_0 "Cantera::epsilon_0"

#: Avogadro's Number, /kmol
avogadro = CxxAvogadro

#: The ideal gas constant in J/kmol-K
gas_constant = CxxGasConstant

#: One atmosphere in Pascals
one_atm = CxxOneAtm

#: Boltzmann constant
boltzmann = CxxBoltzmann

#: Planck constant (J/s)
planck = CxxPlanck

#: The Stefan-Boltzmann constant, W/m^2K^4
stefan_boltzmann = CxxStefanBoltz

#: The charge on an electron (C)
electron_charge = CxxElectronCharge

#: The mass of an electron (kg)
electron_mass = CxxElectronMass

#: Faraday constant, C/kmol
faraday = CxxFaraday

#: Speed of Light (m/s).
light_speed = CxxLightSpeed

#: Permeability of free space :math:`\mu_0` in N/A^2.
permeability_0 = CxxPermeability_0

#: Permittivity of free space (Farads/m = C^2/N/m^2)
epsilon_0 = CxxEpsilon_0
