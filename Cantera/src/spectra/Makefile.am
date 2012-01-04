h_sources  = DiatomicMolecule.h LineBroadener.h Nuclei.h \
             rotor.h spectralUtilities.h

cc_sources = LineBroadener.cpp rotor.cpp spectralUtilities.cpp

INCLUDES                = -I../base -I../../../
AM_CXXFLAGS             = $(INCLUDES)

lib_LTLIBRARIES         = $(top_builddir)/build/lib/libctspectra.la
library_includedir      = $(includedir)

#-----------------------
# Cantera Spectra C/C++ library
#-----------------------

__top_builddir__build_lib_libctspectra_la_LDFLAGS      = $(all_libraries) -release $(GENERIC_RELEASE)
__top_builddir__build_lib_libctspectra_la_SOURCES      = $(cc_sources) $(h_sources)