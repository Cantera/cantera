cc_sources              = cxxutils.cpp  importPhase.cpp
AM_CPPFLAGS             = -I../../src/thermo -I../../src/base
AM_CXXFLAGS             = $(AM_CPPFLAGS)

lib_LTLIBRARIES         = $(top_builddir)/build/lib/libctcxx.la
library_includedir      = $(top_builddir)/build/include
library_include_HEADERS = $(h_sources)

#-----------------------
# Cantera Converters C/C++ library
#-----------------------

__top_builddir__build_lib_libctcxx_la_LDFLAGS      = $(all_libraries) -release $(GENERIC_RELEASE)
__top_builddir__build_lib_libctcxx_la_SOURCES      = $(cc_sources) $(h_sources)

CLEANFILES = *.o
