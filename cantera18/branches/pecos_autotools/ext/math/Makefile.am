cc_sources              = mach.cpp printstring.c

fc_sources              = ddaspk.f dgbfa.f dgbsl.f            \
		          dgefa.f dgesl.f dp1vlu.f dpcoef.f   \
			  dpolft.f fdump.f j4save.f pcoef.f   \
			  polfit.f pvalue.f xercnt.f xerhlt.f \
			  xermsg.f xerprn.f xersve.f xgetua.f 

AM_CPPFLAGS             = 
AM_CXXFLAGS             = $(AM_CPPFLAGS)
AM_FCFLAGS              = $(AM_CPPFLAGS)

lib_LTLIBRARIES         = $(top_builddir)/build/lib/libctmath.la
library_includedir      = $(top_builddir)/build/include
library_include_HEADERS = $(h_sources)

#-----------------------
# Cantera Converters C/C++ library
#-----------------------

__top_builddir__build_lib_libctmath_la_LDFLAGS      = $(all_libraries) -release $(GENERIC_RELEASE)
__top_builddir__build_lib_libctmath_la_SOURCES      = $(fc_sources) $(cc_sources)

CLEANFILES = *.o
