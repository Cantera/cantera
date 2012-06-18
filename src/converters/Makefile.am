h_sources  = 
cc_sources = atomicWeightDB.cpp CKReader.cpp  \
	     Reaction.cpp thermoFunctions.cpp ck2ct.cpp CKParser.cpp \
	     ckr_utils.cpp NASA9Parser.cpp Species.cpp writelog.cpp

AM_CPPFLAGS             = -I$(top_builddir)/include
AM_CXXFLAGS             = $(AM_CPPFLAGS)

lib_LTLIBRARIES         = $(top_builddir)/build/lib/libconverters.la
library_includedir      = $(includedir)

#-----------------------
# Cantera Converters C/C++ library
#-----------------------

__top_builddir__build_lib_libconverters_la_LDFLAGS      = $(all_libraries) -release $(GENERIC_RELEASE)
__top_builddir__build_lib_libconverters_la_SOURCES      = $(cc_sources) $(h_sources)

__top_builddir__build_lib_libconverters_include_HEADERS = $(h_sources)
__top_builddir__build_lib_libconverters_includedir = $(prefix)/include
