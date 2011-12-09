h_sources  = importKinetics.h GRI_30_Kinetics.h KineticsFactory.h \
             Kinetics.h GasKinetics.h \
             FalloffFactory.h  ReactionStoichMgr.h reaction_defs.h \
             FalloffMgr.h ThirdBodyMgr.h RateCoeffMgr.h ReactionData.h \
             RxnRates.h Enhanced3BConc.h StoichManager.h solveSP.h \
             AqueousKinetics.h

cc_sources = importKinetics.cpp GRI_30_Kinetics.cpp KineticsFactory.cpp \
             GasKinetics.cpp AqueousKinetics.cpp FalloffFactory.cpp 
	     ReactionStoichMgr.cpp Kinetics.cpp solveSP.cpp

INCLUDES                = -I../base -I../thermo -I../../../
AM_CXXFLAGS             = $(INCLUDES)

lib_LTLIBRARIES         = libkinetics.la
library_includedir      = $(includedir)

#-----------------------
# Cantera Kinetics C/C++ library
#-----------------------

libkinetics_la_LDFLAGS      = $(all_libraries) -release $(GENERIC_RELEASE)
libkinetics_la_SOURCES      = $(cc_sources) $(h_sources)