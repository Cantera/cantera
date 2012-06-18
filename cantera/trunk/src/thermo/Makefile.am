h_sources  = 
cc_sources = Elements.cpp Phase.cpp RedlichKisterVPSSTP.cpp                \
             ThermoPhase.cpp IdealGasPhase.cpp ConstDensityThermo.cpp      \
             SpeciesThermoFactory.cpp ConstCpPoly.cpp Nasa9Poly1.cpp       \
	     Nasa9PolyMultiTempRegion.cpp PDSS_Water.cpp PDSS_HKFT.cpp     \
             Mu0Poly.cpp GeneralSpeciesThermo.cpp SurfPhase.cpp            \
             ThermoFactory.cpp SpeciesThermoInterpType.cpp                 \
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
	     StoichSubstanceSSTP.cpp PureFluidPhase.cpp                    \
	     StoichSubstance.cpp 
#PecosGasPhase.cpp

AM_CPPFLAGS             = -I$(top_builddir)/include
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
