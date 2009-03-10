"""
Physical Constants

These values are the same as those in the C++ header file ct_defs.h in
the Cantera kernel.
"""

import math

## One atmosphere in Pascals
OneAtm = 101325.0

## The ideal gas constant in J/kmo-K
GasConstant = 8314.47215

## Avogadro's Number, /kmol
Avogadro = 6.02214179e26

## The ideal gas constant in cal/mol-K
GasConst_cal_mol_K = 1.987

## Boltzmann-s constant
Boltzmann = GasConstant / Avogadro

## The Stefan-Boltzmann constant, W/m^2K^4
StefanBoltz = 5.6704004e-8

## The charge on an electron (C)
ElectronCharge = 1.60217648740e-19

## The mass of an electron (kg)
ElectronMass = 9.1093821545e-31

Pi = 3.1415926

## Faraday's constant, C/kmol
Faraday = ElectronCharge * Avogadro

## Planck's constant (J/s)
Planck = 6.6262e-34

## Permittivity of free space
epsilon_0 = 8.85417817e-12  ## Farads/m = C^2/N/m^2

## Permeability of free space \f$ \mu_0 \f$ in N/A^2.
permeability_0 = 4.0e-7*Pi; ## N/A^2

## Speed of Light (m/s).
lightSpeed = 1.0/math.sqrt(epsilon_0 * permeability_0);
