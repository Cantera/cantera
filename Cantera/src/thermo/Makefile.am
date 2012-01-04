h_sources  = State.h Elements.h Constituents.h Phase.h mix_defs.h  \
             ThermoPhase.h IdealGasPhase.h ConstDensityThermo.h    \
             SpeciesThermoFactory.h ThermoFactory.h                \
             NasaPoly1.h NasaPoly2.h NasaThermo.h Nasa9Poly1.h     \
	     Nasa9PolyMultiTempRegion.h                            \
             ShomateThermo.h ShomatePoly.h ConstCpPoly.h           \
             SimpleThermo.h SpeciesThermoMgr.h                     \
             SpeciesThermoInterpType.h                             \
             GeneralSpeciesThermo.h Mu0Poly.h                      \
             speciesThermoTypes.h SpeciesThermo.h SurfPhase.h      \
             EdgePhase.h IdealSolnGasVPSS.h                        \
             VPSSMgr.h  VPSSMgrFactory.h  VPSSMgr_General.h        \
             VPSSMgr_IdealGas.h VPSSMgr_ConstVol.h PDSS_ConstVol.h \
	     PDSS_IdealGas.h PDSS_SSVol.h

cc_sources = State.cpp Elements.cpp Constituents.cpp Phase.cpp             \
             ThermoPhase.cpp IdealGasPhase.cpp ConstDensityThermo.cpp      \
             SpeciesThermoFactory.cpp ConstCpPoly.cpp Nasa9Poly1.cpp       \
	     Nasa9PolyMultiTempRegion.cpp                                  \
             Mu0Poly.cpp GeneralSpeciesThermo.cpp SurfPhase.cpp            \
             ThermoFactory.cpp phasereport.cpp SpeciesThermoInterpType.cpp \
             VPSSMgr.cpp  VPSSMgrFactory.cpp  VPSSMgr_General.cpp          \
	     IdealSolnGasVPSS.cpp                                          \
             VPSSMgr_IdealGas.cpp VPSSMgr_ConstVol.cpp PDSS_ConstVol.cpp   \
	     PDSS_IdealGas.cpp PDSS_SSVol.cpp 

INCLUDES                = -I../base -I../../../
AM_CXXFLAGS             = $(INCLUDES)

lib_LTLIBRARIES         = $(top_builddir)/build/lib/libthermo.la
library_includedir      = $(includedir)

#-----------------------
# Cantera thermo C/C++ library
#-----------------------

__top_builddir__build_lib_libthermo_la_LDFLAGS      = $(all_libraries) -release $(GENERIC_RELEASE)
__top_builddir__build_lib_libthermo_la_SOURCES      = $(cc_sources) $(h_sources)
