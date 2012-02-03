h_sources  = importKinetics.h GRI_30_Kinetics.h KineticsFactory.h      \
             Kinetics.h GasKinetics.h                                  \
             FalloffFactory.h  ReactionStoichMgr.h reaction_defs.h     \
             FalloffMgr.h ThirdBodyMgr.h RateCoeffMgr.h ReactionData.h \
             RxnRates.h Enhanced3BConc.h StoichManager.h solveSP.h     \
             AqueousKinetics.h InterfaceKinetics.h ImplicitSurfChem.h  \
             Group.h

cc_sources = importKinetics.cpp GRI_30_Kinetics.cpp KineticsFactory.cpp  \
             GasKinetics.cpp AqueousKinetics.cpp FalloffFactory.cpp      \
	     ReactionStoichMgr.cpp Kinetics.cpp solveSP.cpp              \
	     ReactionPath.cpp InterfaceKinetics.cpp ImplicitSurfChem.cpp \
             Group.cpp

AM_CPPFLAGS             = -I$(top_builddir)/build/include/ -I$(top_builddir)/build/include/cantera -I$(top_builddir)/build/include/cantera/kernel
AM_CXXFLAGS             = $(AM_CPPFLAGS)

lib_LTLIBRARIES         = $(top_builddir)/build/lib/libkinetics.la
library_includedir      = $(includedir)

#-----------------------
# Cantera Kinetics C/C++ library
#-----------------------

__top_builddir__build_lib_libkinetics_la_LDFLAGS      = $(all_libraries) -release $(GENERIC_RELEASE)
__top_builddir__build_lib_libkinetics_la_SOURCES      = $(cc_sources) $(h_sources)