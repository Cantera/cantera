h_sources  = ChemEquil.h MultiPhase.h MultiPhaseEquil.h equil.h \
	     PropertyCalculator.h

cc_sources = BasisOptimize.cpp ChemEquil.cpp MultiPhase.cpp \
             MultiPhaseEquil.cpp equilibrate.cpp

AM_CPPFLAGS    = -I../base -I../thermo -I../numerics -I../../../ 
AM_CPPFLAGS   += -I../../../ext/f2c_libs/
AM_CXXFLAGS = $(AM_CPPFLAGS)

lib_LTLIBRARIES         = $(top_builddir)/build/lib/libequil.la
library_includedir      = $(includedir)

#-----------------------
# Cantera Converters C/C++ library
#-----------------------

__top_builddir__build_lib_libequil_la_LDFLAGS      = $(all_libraries) -release $(GENERIC_RELEASE)
__top_builddir__build_lib_libequil_la_SOURCES      = $(cc_sources) $(h_sources)