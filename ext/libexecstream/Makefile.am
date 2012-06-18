cc_sources              = exec-stream.cpp
h_sources               = exec-stream.h

AM_CPPFLAGS             = -I.
AM_CXXFLAGS             = $(AM_CPPFLAGS)
AM_FCFLAGS              = $(AM_CPPFLAGS)

lib_LTLIBRARIES         = $(top_builddir)/build/lib/libexecstream.la
library_includedir      = $(top_builddir)/include
library_include_HEADERS = $(h_sources)

#-----------------------
# Execstream library
#-----------------------

__top_builddir__build_lib_libexecstream_la_LDFLAGS      = $(all_libraries) -release $(GENERIC_RELEASE)
__top_builddir__build_lib_libexecstream_la_SOURCES      = $(cc_sources) $(h_sources)

CLEANFILES = *.o
