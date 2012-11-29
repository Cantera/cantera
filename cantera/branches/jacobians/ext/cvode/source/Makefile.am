cc_sources              = band.c cvband.c cvbandpre.c cvdense.c cvdiag.c   \
                          cvode.c cvspgmr.c dense.c iterativ.c llnlmath.c  \
	                  nvector.c spgmr.c

h_sources               = 
AM_CPPFLAGS             = -I../include
AM_CXXFLAGS             = $(AM_CPPFLAGS)
AM_FCFLAGS              = $(AM_CPPFLAGS)

lib_LTLIBRARIES         = $(top_builddir)/build/lib/libcvode.la
library_includedir      = $(top_builddir)/include
library_include_HEADERS = $(h_sources)

#-----------------------
# Cantera Converters C/C++ library
#-----------------------

__top_builddir__build_lib_libcvode_la_LDFLAGS      = $(all_libraries) -release $(GENERIC_RELEASE)
__top_builddir__build_lib_libcvode_la_SOURCES      = $(fc_sources) $(cc_sources)

__top_builddir__build_lib_libcvode_include_HEADERS = $(h_sources)
__top_builddir__build_lib_libcvode_includedir = $(prefix)/include

CLEANFILES = *.o
