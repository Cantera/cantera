! This module is the only 'public' one - that is, the only one visible in
! an application program. It's primary purpose is to provide generic
! procedure names that map to specific procedures depending on the
! argument types.
!
! @warning  The Fortran API is an experimental part of %Cantera and
!   may be changed or removed without notice.

! This file is part of Cantera. See License.txt in the top-level directory or
! at https://cantera.org/license.txt for license and copyright information.

MODULE CANTERA

  USE cantera_thermo
  USE cantera_kinetics
  USE cantera_transport
  USE cantera_funcs

  INTERFACE addCanteraDirectory
     MODULE PROCEDURE ctfunc_addCanteraDirectory
  END INTERFACE addCanteraDirectory

  INTERFACE advanceCoverages
     MODULE PROCEDURE ctkin_advanceCoverages
  END INTERFACE advanceCoverages

  INTERFACE chemPotentials
     MODULE PROCEDURE ctthermo_chemPotentials
  END INTERFACE chemPotentials

  INTERFACE cp_mass
     MODULE PROCEDURE ctthermo_cp_mass
  END INTERFACE cp_mass

  INTERFACE cp_mole
     MODULE PROCEDURE ctthermo_cp_mole
  END INTERFACE cp_mole

  INTERFACE cv_mass
     MODULE PROCEDURE ctthermo_cv_mass
  END INTERFACE cv_mass

  INTERFACE cv_mole
     MODULE PROCEDURE ctthermo_cv_mole
  END INTERFACE cv_mole

  INTERFACE density
     MODULE PROCEDURE ctthermo_density
  END INTERFACE density

  INTERFACE elementIndex
     MODULE PROCEDURE ctthermo_elementIndex
  END INTERFACE elementIndex

  INTERFACE enthalpy_mass
     MODULE PROCEDURE ctthermo_enthalpy_mass
  END INTERFACE enthalpy_mass

  INTERFACE enthalpy_mole
     MODULE PROCEDURE ctthermo_enthalpy_mole
  END INTERFACE enthalpy_mole

  INTERFACE entropy_mass
     MODULE PROCEDURE ctthermo_entropy_mass
  END INTERFACE entropy_mass

  INTERFACE entropy_mole
     MODULE PROCEDURE ctthermo_entropy_mole
  END INTERFACE entropy_mole

  INTERFACE getEosType
     MODULE PROCEDURE ctthermo_getEosType
  END INTERFACE getEosType

  INTERFACE equilibrate
     MODULE PROCEDURE ctthermo_equilibrate
  END INTERFACE equilibrate

  INTERFACE getAtomicWeights
     MODULE PROCEDURE ctthermo_getAtomicWeights
  END INTERFACE getAtomicWeights

  INTERFACE getBinDiffCoeffs
     MODULE PROCEDURE ctrans_getBinDiffCoeffs
  END INTERFACE getBinDiffCoeffs

  INTERFACE getCanteraError
     MODULE PROCEDURE ctfunc_getCanteraError
  END INTERFACE getCanteraError

  INTERFACE turnOffExitOnError
     MODULE PROCEDURE ctfunc_TurnOffExitOnError
  END INTERFACE turnOffExitOnError

  INTERFACE getCp_R
     MODULE PROCEDURE ctthermo_getCp_R
  END INTERFACE getCp_R

  INTERFACE getCreationRates
     MODULE PROCEDURE ctkin_getCreationRates
  END INTERFACE getCreationRates

  INTERFACE getDestructionRates
     MODULE PROCEDURE ctkin_getDestructionRates
  END INTERFACE getDestructionRates

  INTERFACE getElementName
     MODULE PROCEDURE ctthermo_getElementName
  END INTERFACE getElementName

  INTERFACE getEnthalpies_RT
     MODULE PROCEDURE ctthermo_getEnthalpies_RT
  END INTERFACE getEnthalpies_RT

  INTERFACE getGibbs_RT
     MODULE PROCEDURE ctthermo_getGibbs_RT
  END INTERFACE getGibbs_RT

  INTERFACE getEntropies_R
     MODULE PROCEDURE ctthermo_getEntropies_R
  END INTERFACE getEntropies_R

  INTERFACE getEquilibriumConstants
     MODULE PROCEDURE ctkin_getEquilibriumConstants
  END INTERFACE getEquilibriumConstants

  INTERFACE getFwdRatesOfProgress
     MODULE PROCEDURE ctkin_getFwdRatesOfProgress
  END INTERFACE getFwdRatesOfProgress

  INTERFACE getMassFractions
     MODULE PROCEDURE ctthermo_getMassFractions
  END INTERFACE getMassFractions

  INTERFACE getMixDiffCoeffs
     MODULE PROCEDURE ctrans_getMixDiffCoeffs
  END INTERFACE getMixDiffCoeffs

  INTERFACE getMixDiffCoeffsMole
     MODULE PROCEDURE ctrans_getMixDiffCoeffsMole
  END INTERFACE getMixDiffCoeffsMole

  INTERFACE getMixDiffCoeffsMass
     MODULE PROCEDURE ctrans_getMixDiffCoeffsMass
  END INTERFACE getMixDiffCoeffsMass

  INTERFACE getMoleFractions
     MODULE PROCEDURE ctthermo_getMoleFractions
  END INTERFACE getMoleFractions

  INTERFACE getMolecularWeights
     MODULE PROCEDURE ctthermo_getMolecularWeights
  END INTERFACE getMolecularWeights

  INTERFACE getMultiDiffCoeffs
     MODULE PROCEDURE ctrans_getMultiDiffCoeffs
  END INTERFACE getMultiDiffCoeffs

  INTERFACE getNetProductionRates
     MODULE PROCEDURE ctkin_getNetProductionRates
  END INTERFACE getNetProductionRates

  INTERFACE getNetRatesOfProgress
     MODULE PROCEDURE ctkin_getNetRatesOfProgress
  END INTERFACE getNetRatesOfProgress

  INTERFACE getPartialMolarIntEnergies
     MODULE PROCEDURE ctthermo_getPartialMolarIntEnerg_R
  END INTERFACE getPartialMolarIntEnergies

  INTERFACE getPartialMolarCp
     MODULE PROCEDURE ctthermo_getPartialMolarCp
  END INTERFACE getPartialMolarCp

  INTERFACE getPartialMolarEnthalpies
     MODULE PROCEDURE ctthermo_getPartialMolarEnthalpies
  END INTERFACE getPartialMolarEnthalpies

  INTERFACE getReactionString
     MODULE PROCEDURE ctkin_getReactionString
  END INTERFACE getReactionString

  INTERFACE getRevRatesOfProgress
     MODULE PROCEDURE ctkin_getRevRatesOfProgress
  END INTERFACE getRevRatesOfProgress

  INTERFACE getSpeciesName
     MODULE PROCEDURE ctthermo_getSpeciesName
  END INTERFACE getSpeciesName

  INTERFACE getThermalDiffCoeffs
     MODULE PROCEDURE ctrans_getThermalDiffCoeffs
  END INTERFACE getThermalDiffCoeffs

  INTERFACE gibbs_mass
     MODULE PROCEDURE ctthermo_gibbs_mass
  END INTERFACE gibbs_mass

  INTERFACE gibbs_mole
     MODULE PROCEDURE ctthermo_gibbs_mole
  END INTERFACE gibbs_mole

  INTERFACE importPhase
     MODULE PROCEDURE ctfunc_importPhase
  END INTERFACE importPhase

  INTERFACE importInterface
     MODULE PROCEDURE ctfunc_importInterface
  END INTERFACE importInterface

  INTERFACE intEnergy_mass
     MODULE PROCEDURE ctthermo_intEnergy_mass
  END INTERFACE intEnergy_mass

  INTERFACE intEnergy_mole
     MODULE PROCEDURE ctthermo_intEnergy_mole
  END INTERFACE intEnergy_mole

  INTERFACE isReversible
     MODULE PROCEDURE ctkin_isReversible
  END INTERFACE isReversible

  INTERFACE kineticsSpeciesIndex
     MODULE PROCEDURE ctkin_kineticsSpeciesIndex
  END INTERFACE kineticsSpeciesIndex

  INTERFACE kineticsStart
     MODULE PROCEDURE ctkin_kineticsStart
  END INTERFACE kineticsStart

  INTERFACE getKineticsType
     MODULE PROCEDURE ctkin_getKineticsType
  END INTERFACE getKineticsType

  INTERFACE massFraction
     MODULE PROCEDURE ctthermo_massFraction
  END INTERFACE massFraction

  INTERFACE maxTemp
     MODULE PROCEDURE ctthermo_maxTemp
  END INTERFACE maxTemp

  INTERFACE meanMolecularWeight
     MODULE PROCEDURE ctthermo_meanMolecularWeight
  END INTERFACE meanMolecularWeight

  INTERFACE minTemp
     MODULE PROCEDURE ctthermo_minTemp
  END INTERFACE minTemp

  INTERFACE molarDensity
     MODULE PROCEDURE ctthermo_molarDensity
  END INTERFACE molarDensity

  INTERFACE moleFraction
     MODULE PROCEDURE ctthermo_moleFraction
  END INTERFACE moleFraction

  INTERFACE multiplier
     MODULE PROCEDURE ctkin_multiplier
  END INTERFACE multiplier

  INTERFACE nAtoms
     MODULE PROCEDURE ctthermo_nAtoms
  END INTERFACE nAtoms

  INTERFACE nElements
     MODULE PROCEDURE ctthermo_nElements
  END INTERFACE nElements

  INTERFACE nReactions
     MODULE PROCEDURE ctkin_nReactions
  END INTERFACE nReactions

  INTERFACE nSpecies
     MODULE PROCEDURE ctthermo_nSpecies
  END INTERFACE nSpecies

  INTERFACE nTotalSpecies
     MODULE PROCEDURE ctkin_nTotalSpecies
  END INTERFACE nTotalSpecies

  INTERFACE nPhases
     MODULE PROCEDURE ctkin_nPhases
  END INTERFACE nPhases

  INTERFACE phaseIndex
     MODULE PROCEDURE ctkin_phaseIndex
  END INTERFACE phaseIndex

  INTERFACE phase_report
     MODULE PROCEDURE ctfunc_phase_report
  END INTERFACE phase_report

  INTERFACE pressure
     MODULE PROCEDURE ctthermo_pressure
  END INTERFACE pressure

  INTERFACE productStoichCoeff
     MODULE PROCEDURE ctkin_productStoichCoeff
  END INTERFACE productStoichCoeff

  INTERFACE reactantStoichCoeff
     MODULE PROCEDURE ctkin_reactantStoichCoeff
  END INTERFACE reactantStoichCoeff

  INTERFACE getReactionType
     MODULE PROCEDURE ctkin_getReactionType
  END INTERFACE getReactionType

  INTERFACE refPressure
     MODULE PROCEDURE ctthermo_refPressure
  END INTERFACE refPressure

  INTERFACE setDensity
     MODULE PROCEDURE ctthermo_setDensity
  END INTERFACE setDensity

  INTERFACE setMassFractions
     MODULE PROCEDURE ctthermo_setMassFractions
  END INTERFACE setMassFractions

  INTERFACE setMassFractionsByName
     MODULE PROCEDURE ctthermo_setMassFractionsByName
  END INTERFACE setMassFractionsByName

  INTERFACE setMoleFractions
     MODULE PROCEDURE ctthermo_setMoleFractions
  END INTERFACE setMoleFractions

  INTERFACE setMoleFractionsByName
     MODULE PROCEDURE ctthermo_setMoleFractionsByName
  END INTERFACE setMoleFractionsByName

  INTERFACE setMultiplier
     MODULE PROCEDURE ctkin_setMultiplier
  END INTERFACE setMultiplier

  INTERFACE setPressure
     MODULE PROCEDURE ctthermo_setPressure
  END INTERFACE setPressure

  INTERFACE setState_HP
     MODULE PROCEDURE ctthermo_setState_HP
  END INTERFACE setState_HP

  INTERFACE setState_SP
     MODULE PROCEDURE ctthermo_setState_SP
  END INTERFACE setState_SP

  INTERFACE setState_SV
     MODULE PROCEDURE ctthermo_setState_SV
  END INTERFACE setState_SV

  INTERFACE setState_TPX
     MODULE PROCEDURE ctthermo_setState_TPX
     MODULE PROCEDURE ctstring_setState_TPX
  END INTERFACE setState_TPX

  INTERFACE setState_TRX
     MODULE PROCEDURE ctthermo_setState_TRX
     MODULE PROCEDURE ctstring_setState_TRX
  END INTERFACE setState_TRX

  INTERFACE setState_TRY
     MODULE PROCEDURE ctthermo_setState_TRY
     MODULE PROCEDURE ctstring_setState_TRY
  END INTERFACE setState_TRY

  INTERFACE setState_TPY
     MODULE PROCEDURE ctthermo_setState_TPY
     MODULE PROCEDURE ctstring_setState_TPY
  END INTERFACE setState_TPY

  INTERFACE setState_UV
     MODULE PROCEDURE ctthermo_setState_UV
  END INTERFACE setState_UV

  INTERFACE setTemperature
     MODULE PROCEDURE ctthermo_setTemperature
  END INTERFACE setTemperature

  INTERFACE speciesIndex
     MODULE PROCEDURE ctthermo_speciesIndex
  END INTERFACE speciesIndex

  INTERFACE temperature
     MODULE PROCEDURE ctthermo_temperature
  END INTERFACE temperature

  INTERFACE electricalConductivity
     MODULE PROCEDURE ctrans_electricalConductivity
  END INTERFACE electricalConductivity

  INTERFACE thermalConductivity
     MODULE PROCEDURE ctrans_thermalConductivity
  END INTERFACE thermalConductivity

  INTERFACE viscosity
     MODULE PROCEDURE ctrans_viscosity
  END INTERFACE viscosity

END MODULE CANTERA
