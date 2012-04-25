h_sources  = 
cc_sources = LineBroadener.cpp rotor.cpp spectralUtilities.cpp

AM_CPPFLAGS             = -I$(top_builddir)/include
AM_CXXFLAGS             = $(AM_CPPFLAGS)

lib_LTLIBRARIES         = $(top_builddir)/build/lib/libctspectra.la
library_includedir      = $(includedir)

#-----------------------
# Cantera Spectra C/C++ library
#-----------------------

__top_builddir__build_lib_libctspectra_la_LDFLAGS      = $(all_libraries) -release $(GENERIC_RELEASE)
__top_builddir__build_lib_libctspectra_la_SOURCES      = $(cc_sources) $(h_sources)
