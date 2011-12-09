h_sources  = ChemEquil.h MultiPhase.h MultiPhaseEquil.h equil.h \
	     PropertyCalculator.h

cc_sources = BasisOptimize.cpp ChemEquil.cpp MultiPhase.cpp \
             MultiPhaseEquil.cpp equilibrate.cpp

INCLUDES    = -I../base -I../thermo -I../numerics -I../../../ 
INCLUDES   += -I../../../ext/f2c_libs/
AM_CXXFLAGS = $(INCLUDES)

lib_LTLIBRARIES         = libequil.la
library_includedir      = $(includedir)

#-----------------------
# Cantera Converters C/C++ library
#-----------------------

libequil_la_LDFLAGS      = $(all_libraries) -release $(GENERIC_RELEASE)
libequil_la_SOURCES      = $(cc_sources) $(h_sources)