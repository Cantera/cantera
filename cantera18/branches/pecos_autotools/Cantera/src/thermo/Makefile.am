h_sources  = State.h Elements.h Constituents.h Phase.h mix_defs.h  \
             ThermoPhase.h IdealGasPhase.h ConstDensityThermo.h    \
             SpeciesThermoFactory.h ThermoFactory.h                \
             NasaPoly1.h NasaPoly2.h NasaThermo.h Nasa9Poly1.h     \
	     Nasa9PolyMultiTempRegion.h WaterPropsIAPWS.h          \
             ShomateThermo.h ShomatePoly.h ConstCpPoly.h           \
             SimpleThermo.h SpeciesThermoMgr.h                     \
             SpeciesThermoInterpType.h VPSSMgr_Water_ConstVol.h    \
             GeneralSpeciesThermo.h Mu0Poly.h                      \
             speciesThermoTypes.h SpeciesThermo.h SurfPhase.h      \
             EdgePhase.h IdealSolnGasVPSS.h MolalityVPSSTP.h       \
             VPSSMgr.h  VPSSMgrFactory.h  VPSSMgr_General.h        \
             VPSSMgr_IdealGas.h VPSSMgr_ConstVol.h PDSS_ConstVol.h \
	     PDSS_IdealGas.h PDSS_SSVol.h DebyeHuckel.h            \
	     VPStandardStateTP.h PDSS_Water.h WaterProps.h         \
	     WaterPropsIAPWSphi.h VPSSMgr_Water_HKFT.h PDSS_HKFT.h \
             PDSS_IonsFromNeutral.h IonsFromNeutralVPSSTP.h        \
             GibbsExcessVPSSTP.h LatticePhase.h IdealMolalSoln.h   \
	     PDSS.h HMWSoln.h WaterSSTP.h MetalSHEelectrons.h      \
	     IdealSolidSolnPhase.h StoichSubstanceSSTP.h           \
             SingleSpeciesTP.h MineralEQ3.h PseudoBinaryVPSSTP.h   \
             MargulesVPSSTP.h StoichSubstance.h electrolytes.h     \
             PureFluidPhase.h LatticeSolidPhase.h MetalPhase.h     \
	     VPSSMgr_types.h PecosGasPhase.h StatMech.h

cc_sources = State.cpp Elements.cpp Constituents.cpp Phase.cpp             \
             ThermoPhase.cpp IdealGasPhase.cpp ConstDensityThermo.cpp      \
             SpeciesThermoFactory.cpp ConstCpPoly.cpp Nasa9Poly1.cpp       \
	     Nasa9PolyMultiTempRegion.cpp PDSS_Water.cpp PDSS_HKFT.cpp     \
             Mu0Poly.cpp GeneralSpeciesThermo.cpp SurfPhase.cpp            \
             ThermoFactory.cpp phasereport.cpp SpeciesThermoInterpType.cpp \
             VPSSMgr.cpp  VPSSMgrFactory.cpp  VPSSMgr_General.cpp          \
	     IdealSolnGasVPSS.cpp MolalityVPSSTP.cpp VPStandardStateTP.cpp \
             VPSSMgr_IdealGas.cpp VPSSMgr_ConstVol.cpp PDSS_ConstVol.cpp   \
	     PDSS_IdealGas.cpp PDSS_SSVol.cpp DebyeHuckel.cpp PDSS.cpp     \
	     WaterProps.cpp WaterPropsIAPWS.cpp	WaterPropsIAPWSphi.cpp     \
             VPSSMgr_Water_HKFT.cpp VPSSMgr_Water_ConstVol.cpp             \
             PDSS_IonsFromNeutral.cpp IonsFromNeutralVPSSTP.cpp            \
             GibbsExcessVPSSTP.cpp LatticePhase.cpp IdealMolalSoln.cpp     \
	     HMWSoln.cpp HMWSoln_input.cpp WaterSSTP.cpp                   \
             MetalSHEelectrons.cpp                                         \
             IdealSolidSolnPhase.cpp LatticeSolidPhase.cpp                 \
             SingleSpeciesTP.cpp MineralEQ3.cpp                            \
             PseudoBinaryVPSSTP.cpp MargulesVPSSTP.cpp                     \
	     StoichSubstanceSSTP.cpp PureFluidPhase.cpp StatMech.cpp       \
	     StoichSubstance.cpp PecosGasPhase.cpp

AM_CPPFLAGS             = -I../base
AM_CXXFLAGS             = $(AM_CPPFLAGS)

lib_LTLIBRARIES         = $(top_builddir)/build/lib/libthermo.la
library_includedir      = $(includedir)

#-----------------------
# Cantera thermo C/C++ library
#-----------------------

__top_builddir__build_lib_libthermo_la_LDFLAGS      = $(all_libraries) -release $(GENERIC_RELEASE)
__top_builddir__build_lib_libthermo_la_SOURCES      = $(cc_sources) $(h_sources)

__top_builddir__build_lib_libthermo_include_HEADERS = $(h_sources)
__top_builddir__build_lib_libthermo_includedir = $(prefix)/include

# header file accumulation
all:
	@echo copying headers $<
	cp -p *.h $(top_builddir)/build/include/
