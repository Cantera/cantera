cc_sources              = user.cpp

AM_CPPFLAGS             = -I$(top_builddir)/build/include/ -I$(top_builddir)/build/include/cantera
AM_CXXFLAGS             = $(AM_CPPFLAGS)

lib_LTLIBRARIES         = $(top_builddir)/build/lib/libuser.la
library_includedir      = $(top_builddir)/build/include
library_include_HEADERS = $(h_sources)

#-----------------------
# Cantera Converters C/C++ library
#-----------------------

__top_builddir__build_lib_libuser_la_LDFLAGS      = $(all_libraries) -release $(GENERIC_RELEASE)
__top_builddir__build_lib_libuser_la_SOURCES      = $(cc_sources) $(h_sources)

CLEANFILES = *.o
