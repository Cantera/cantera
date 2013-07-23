
h_sources  = $(top_builddir)/include/cantera/spectra/LineBroadener.h \
             $(top_builddir)/include/cantera/spectra/rotor.h
#             $(top_builddir)/include/cantera/spectra/spectralUtilities.h

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

__top_builddir__build_lib_libctspectra_include_HEADERS = $(h_sources)
__top_builddir__build_lib_libctspectra_includedir = $(prefix)/include/cantera/spectra
