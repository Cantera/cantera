h_sources  = DiatomicMolecule.h LineBroadener.h Nuclei.h \
             rotor.h spectralUtilities.h

cc_sources = LineBroadener.cpp rotor.cpp spectralUtilities.cpp

INCLUDES                = -I../base -I../../../
AM_CXXFLAGS             = $(INCLUDES)

lib_LTLIBRARIES         = libctspectra.la
library_includedir      = $(includedir)

#-----------------------
# Cantera Spectra C/C++ library
#-----------------------

libctspectra_la_LDFLAGS      = $(all_libraries) -release $(GENERIC_RELEASE)
libctspectra_la_SOURCES      = $(cc_sources) $(h_sources)