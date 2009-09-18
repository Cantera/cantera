echo on
cd ..\..\..\

if not exist build\include\cantera\kernel mkdir build\include\cantera\kernel

cd Cantera\src\thermo


copy AdsorbateThermo.h          ..\..\..\build\include\cantera\kernel
copy ConstCpPoly.h          ..\..\..\build\include\cantera\kernel
copy Constituents.h          ..\..\..\build\include\cantera\kernel
copy Crystal.h          ..\..\..\build\include\cantera\kernel
copy DebyeHuckel.h          ..\..\..\build\include\cantera\kernel
copy EdgePhase.h          ..\..\..\build\include\cantera\kernel
copy electrolytes.h          ..\..\..\build\include\cantera\kernel
copy Elements.h          ..\..\..\build\include\cantera\kernel
copy GeneralSpeciesThermo.h          ..\..\..\build\include\cantera\kernel
copy GibbsExcessVPSSTP.h          ..\..\..\build\include\cantera\kernel
copy HMWSoln.h          ..\..\..\build\include\cantera\kernel
copy IdealGasPhase.h          ..\..\..\build\include\cantera\kernel
copy IdealMolalSoln.h          ..\..\..\build\include\cantera\kernel
copy IdealSolidSolnPhase.h          ..\..\..\build\include\cantera\kernel
copy IdealSolnGasVPSS.h          ..\..\..\build\include\cantera\kernel
copy LatticePhase.h          ..\..\..\build\include\cantera\kernel
copy LatticeSolidPhase.h          ..\..\..\build\include\cantera\kernel
copy MargulesVPSSTP.h          ..\..\..\build\include\cantera\kernel
copy MetalPhase.h          ..\..\..\build\include\cantera\kernel
copy MineralEQ3.h          ..\..\..\build\include\cantera\kernel
copy mix_defs.h          ..\..\..\build\include\cantera\kernel
copy MolalityVPSSTP.h          ..\..\..\build\include\cantera\kernel
copy Mu0Poly.h          ..\..\..\build\include\cantera\kernel
copy Nasa9Poly1.h          ..\..\..\build\include\cantera\kernel
copy Nasa9PolyMultiTempRegion.h          ..\..\..\build\include\cantera\kernel
copy NasaPoly1.h          ..\..\..\build\include\cantera\kernel
copy NasaPoly2.h          ..\..\..\build\include\cantera\kernel
copy NasaThermo.h          ..\..\..\build\include\cantera\kernel
copy PDSS.h          ..\..\..\build\include\cantera\kernel
copy PDSS_ConstVol.h          ..\..\..\build\include\cantera\kernel
copy PDSS_HKFT.h          ..\..\..\build\include\cantera\kernel
copy PDSS_IdealGas.h          ..\..\..\build\include\cantera\kernel
copy PDSS_Water.h          ..\..\..\build\include\cantera\kernel
copy Phase.h          ..\..\..\build\include\cantera\kernel
copy PseudoBinaryVPSSTP.h          ..\..\..\build\include\cantera\kernel
copy PureFluidPhase.h          ..\..\..\build\include\cantera\kernel
copy SemiconductorPhase.h          ..\..\..\build\include\cantera\kernel
copy ShomateThermo.h          ..\..\..\build\include\cantera\kernel
copy SimpleThermo.h          ..\..\..\build\include\cantera\kernel
copy SingleSPeciesTP.h          ..\..\..\build\include\cantera\kernel
copy SpeciesThermo.h          ..\..\..\build\include\cantera\kernel
copy SpeciesThermoFactory.h          ..\..\..\build\include\cantera\kernel
copy SpeciesThermoInterpType.h          ..\..\..\build\include\cantera\kernel
copy SpeciesThermoMgr.h          ..\..\..\build\include\cantera\kernel
copy SpeciesThermoTypes.h          ..\..\..\build\include\cantera\kernel
copy State.h          ..\..\..\build\include\cantera\kernel
copy StoichSubstance.h          ..\..\..\build\include\cantera\kernel
copy StoichSubstanceSSTP.h          ..\..\..\build\include\cantera\kernel
copy SurfPhase.h          ..\..\..\build\include\cantera\kernel
copy ThermoFactory.h          ..\..\..\build\include\cantera\kernel
copy ThermoPhase.h          ..\..\..\build\include\cantera\kernel
copy VPSSMgr.h          ..\..\..\build\include\cantera\kernel
copy VPSSMGr_ConstVol.h          ..\..\..\build\include\cantera\kernel
copy VPSSMgr_General.h          ..\..\..\build\include\cantera\kernel
copy VPSSMgr_IdealGas.h          ..\..\..\build\include\cantera\kernel
copy VPSSMgr_types.h          ..\..\..\build\include\cantera\kernel
copy VPSSMgr_Water_ConstVol.h          ..\..\..\build\include\cantera\kernel
copy VPSSMgr_Water_HKFT.h          ..\..\..\build\include\cantera\kernel
copy VPSSMgrFactory.h          ..\..\..\build\include\cantera\kernel
copy VPStandardStateTP.h          ..\..\..\build\include\cantera\kernel
copy WaterProps.h          ..\..\..\build\include\cantera\kernel
copy WaterPropsIAPWS.h          ..\..\..\build\include\cantera\kernel
copy WaterPropsIAPWSphi.h          ..\..\..\build\include\cantera\kernel
copy WaterSSTP.h          ..\..\..\build\include\cantera\kernel

cd ..\..\..\win32\vc9\thermo
echo off
echo 'ok'
