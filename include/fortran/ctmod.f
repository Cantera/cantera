! ------------------------------------------------------------------------
! 
! Copyright 2002 California Institute of Technology
! 
! Version 1.2.0.1
! Mon Jan  7 07:13:18 2002
! This file generated automatically.
! Cantera Fortran 90 Interface module
! ------------------------------------------------------------------------

      module Cantera

      use cttypes
      use ctmixture
      use cttransport
      use ctthermo
      use ctreactor
      use ctfdev
      use ctutils

      double precision, parameter :: OneAtm = 1.01325D5
      double precision, parameter :: Avogadro = 6.022136736D23
      double precision, parameter :: GasConstant = 8314.0D0
      double precision, parameter :: StefanBoltz = 5.67D-8

      interface assignment (=)
         module procedure mix_copy
         module procedure reac_copy
      end interface

      interface addDirectory
      module procedure util_addDirectory
      end interface

      interface addElement
      module procedure mix_addElement
      end interface

      interface advance
      module procedure reac_advance
      end interface

      interface atomicWeight
      module procedure mix_atwt
      end interface

      interface charge
      module procedure mix_charge
      end interface

      interface contents
      module procedure reac_contents
      end interface

      interface copy
      module procedure mix_copy
      module procedure reac_copy
      module procedure fdev_copy
      end interface

      interface cp_mass
      module procedure mix_cp_mass
      end interface

      interface cp_mole
      module procedure mix_cp_mole
      end interface

      interface critPressure
      module procedure mix_critpres
      end interface

      interface critTemperature
      module procedure mix_crittemp
      end interface

      interface cv_mass
      module procedure mix_cv_mass
      end interface

      interface cv_mole
      module procedure mix_cv_mole
      end interface

      interface delete
      module procedure mix_delete
      module procedure trans_delete
      end interface

      interface density
      module procedure mix_density
      module procedure reac_density
      end interface

      interface disableChemistry
      module procedure reac_chemoff
      end interface

      interface downstream
      module procedure fdev_downstream
      end interface

      interface elementIndex
      module procedure mix_eindex
      end interface

      interface enableChemistry
      module procedure reac_chemon
      end interface

      interface enthalpy_mass
      module procedure mix_enthalpy_mass
      module procedure reac_enthalpy_mass
      end interface

      interface enthalpy_mole
      module procedure mix_enthalpy_mole
      end interface

      interface entropy_mass
      module procedure mix_entropy_mass
      end interface

      interface entropy_mole
      module procedure mix_entropy_mole
      end interface

      interface equationOfState
      module procedure mix_get_eos
      end interface

      interface equilibrate
      module procedure mix_equilib
      end interface

      interface getChemPotentials_RT
      module procedure mix_gchempot
      end interface

      interface getConcentrations
      module procedure mix_getconc
      end interface

      interface getCp_R
      module procedure mix_gcp_r
      end interface

      interface getCreationRates
      module procedure mix_get_cdot
      end interface

      interface getDestructionRates
      module procedure mix_get_ddot
      end interface

      interface getElementNames
      module procedure mix_getElementNames
      end interface

      interface getEnthalpy_RT
      module procedure mix_gh_rt
      end interface

      interface getEntropy_R
      module procedure mix_gs_r
      end interface

      interface getEquilibriumConstants
      module procedure mix_get_kc
      end interface

      interface getFwdRatesOfProgress
      module procedure mix_get_fwdrop
      end interface

      interface getGains
      module procedure fdev_getGains
      end interface

      interface getGibbs_RT
      module procedure mix_ggibbs_rt
      end interface

      interface getMassFractions
      module procedure mix_gety
      end interface

      interface getMoleFractions
      module procedure mix_getx
      end interface

      interface getMolecularWeights
      module procedure mix_gmolwts
      end interface

      interface getMultiDiffCoeffs
      module procedure mix_gmultidiff
      end interface

      interface getNetProductionRates
      module procedure mix_get_wdot
      end interface

      interface getNetRatesOfProgress
      module procedure mix_get_netrop
      end interface

      interface getReactionString
      module procedure mix_getrxnstr
      end interface

      interface getRevRatesOfProgress
      module procedure mix_get_revrop
      end interface

      interface getSpeciesFluxes
      module procedure mix_gflux
      end interface

      interface getSpeciesNames
      module procedure mix_getspnm
      end interface

      interface getSpeciesViscosities
      module procedure mix_gspvisc
      end interface

      interface getThermalDiffCoeffs
      module procedure mix_gtdiff
      end interface

      interface gibbs_mass
      module procedure mix_gibbs_mass
      end interface

      interface gibbs_mole
      module procedure mix_gibbs_mole
      end interface

      interface install
      module procedure fdev_install
      end interface

      interface intEnergy_mass
      module procedure mix_umass
      module procedure reac_umass
      end interface

      interface intEnergy_mole
      module procedure mix_umole
      end interface

      interface mass
      module procedure reac_mass
      end interface

      interface massFlowRate
      module procedure fdev_mfrate
      end interface

      interface massFraction
      module procedure mix_ybyname
      end interface

      interface maxError
      module procedure fdev_maxError
      end interface

      interface maxTemp
      module procedure mix_maxTemp
      end interface

      interface meanMolecularWeight
      module procedure mix_meanmw
      end interface

      interface mean_X
      module procedure mix_mean_X
      end interface

      interface mean_Y
      module procedure mix_mean_Y
      end interface

      interface minTemp
      module procedure mix_minTemp
      end interface

      interface molarDensity
      module procedure mix_moldens
      end interface

      interface moleFraction
      module procedure mix_xbyname
      end interface

      interface molecularWeight
      module procedure mix_molwt
      end interface

      interface nAtoms
      module procedure mix_nAtoms
      end interface

      interface nElements
      module procedure mix_nElements
      end interface

      interface nReactions
      module procedure mix_nReactions
      end interface

      interface nSpecies
      module procedure mix_nSpecies
      end interface

      interface netStoichCoeff
      module procedure mix_nstoich
      end interface

      interface potentialEnergy
      module procedure mix_pe
      end interface

      interface pressure
      module procedure mix_pressure
      module procedure reac_pressure
      end interface

      interface printSummary
      module procedure util_printSummary
      end interface

      interface productStoichCoeff
      module procedure mix_pstoich
      end interface

      interface reactantStoichCoeff
      module procedure mix_rstoich
      end interface

      interface ready
      module procedure mix_ready
      module procedure fdev_ready
      end interface

      interface refPressure
      module procedure mix_refp
      end interface

      interface reset
      module procedure fdev_reset
      end interface

      interface residenceTime
      module procedure reac_residenceTime
      end interface

      interface restoreState
      module procedure mix_restoreState
      end interface

      interface satPressure
      module procedure mix_satpres
      end interface

      interface satTemperature
      module procedure mix_sattemp
      end interface

      interface saveState
      module procedure mix_saveState
      end interface

      interface setArea
      module procedure reac_setArea
      end interface

      interface setConcentrations
      module procedure mix_sconc
      end interface

      interface setDensity
      module procedure mix_setDensity
      end interface

      interface setEmissivity
      module procedure reac_setEmissivity
      end interface

      interface setEquationOfState
      module procedure mix_set_eos
      end interface

      interface setExtPressure
      module procedure reac_setepr
      end interface

      interface setExtRadTemp
      module procedure reac_setExtRadTemp
      end interface

      interface setExtTemp
      module procedure reac_setExtTemp
      end interface

      interface setGains
      module procedure fdev_setGains
      end interface

      interface setHeatTransferCoeff
      module procedure reac_seth
      end interface

      interface setInitialTime
      module procedure reac_setitm
      end interface

      interface setInitialVolume
      module procedure reac_setivol
      end interface

      interface setMassFractions
      module procedure mix_sety
      end interface

      interface setMassFractions_NoNorm
      module procedure mix_synonorm
      end interface

      interface setMaxStep
      module procedure reac_setMaxStep
      end interface

      interface setMixture
      module procedure reac_insmix
      end interface

      interface setMoleFractions
      module procedure mix_setx
      end interface

      interface setMoleFractions_NoNorm
      module procedure mix_sxnonorm
      end interface

      interface setOptions
      module procedure trans_setopt
      end interface

      interface setPotentialEnergy
      module procedure mix_setpe_k
      end interface

      interface setPressure
      module procedure mix_setPressure
      end interface

      interface setSetpoint
      module procedure fdev_setspnt
      end interface

      interface setState_HP
      module procedure mix_setState_HP
      end interface

      interface setState_PX
      module procedure mix_setState_PX
      end interface

      interface setState_PY
      module procedure mix_setState_PY
      end interface

      interface setState_RX
      module procedure mix_setState_RX
      end interface

      interface setState_RY
      module procedure mix_setState_RY
      end interface

      interface setState_SP
      module procedure mix_setState_SP
      end interface

      interface setState_SV
      module procedure mix_setState_SV
      end interface

      interface setState_TP
      module procedure mix_setState_TP
      end interface

      interface setState_TPX
      module procedure mix_setTPXstr
      module procedure mix_setTPXarray
      end interface

      interface setState_TPY
      module procedure mix_setTPYstr
      module procedure mix_setTPYarray
      end interface

      interface setState_TR
      module procedure mix_setState_TR
      end interface

      interface setState_TRX
      module procedure mix_setState_TRX
      end interface

      interface setState_TRY
      module procedure mix_setState_TRY
      end interface

      interface setState_TX
      module procedure mix_setState_TX
      end interface

      interface setState_TY
      module procedure mix_setState_TY
      end interface

      interface setState_UV
      module procedure mix_setState_UV
      end interface

      interface setTemperature
      module procedure mix_settemp
      end interface

      interface setTransport
      module procedure mix_set_trans
      end interface

      interface setVDotCoeff
      module procedure reac_setVDotCoeff
      end interface

      interface speciesIndex
      module procedure mix_speciesIndex
      end interface

      interface sum_xlogQ
      module procedure mix_sum_xlogQ
      end interface

      interface temperature
      module procedure mix_temperature
      module procedure reac_temperature
      end interface

      interface thermalConductivity
      module procedure mix_tcon
      end interface

      interface time
      module procedure reac_time
      end interface

      interface transportMgr
      module procedure mix_transportMgr
      end interface

      interface update
      module procedure fdev_update
      end interface

      interface upstream
      module procedure fdev_upstream
      end interface

      interface viscosity
      module procedure mix_visc
      end interface

      interface volume
      module procedure reac_volume
      end interface

      end module
