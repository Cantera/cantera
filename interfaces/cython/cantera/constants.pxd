# This file is part of Cantera. See License.txt in the top-level directory or
# at https://cantera.org/license.txt for license and copyright information.

#cython: language_level=3
#distutils: language = c++

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
