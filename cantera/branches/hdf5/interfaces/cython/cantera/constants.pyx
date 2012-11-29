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
Avogadro = CxxAvogadro

#: The ideal gas constant in J/kmo-K
GasConstant = CxxGasConstant

#: One atmosphere in Pascals
OneAtm = CxxOneAtm

#: Boltzmann constant
Boltzmann = CxxBoltzmann

#: Planck constant (J/s)
Planck = CxxPlanck

#: The Stefan-Boltzmann constant, W/m^2K^4
StefanBoltz = CxxStefanBoltz

#: The charge on an electron (C)
ElectronCharge = CxxElectronCharge

#: The mass of an electron (kg)
ElectronMass = CxxElectronMass

#: Faraday constant, C/kmol
Faraday = CxxFaraday

#: Speed of Light (m/s).
lightSpeed = CxxLightSpeed

#: Permeability of free space :math:`\mu_0` in N/A^2.
permeability_0 = CxxPermeability_0

#: Permittivity of free space (Farads/m = C^2/N/m^2)
epsilon_0 = CxxEpsilon_0
