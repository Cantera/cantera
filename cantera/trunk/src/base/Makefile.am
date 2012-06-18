h_sources               = LogPrintCtrl.h application.h units.h

cc_sources              = ct2ctml.cpp ctml.cpp plots.cpp                   \
			  stringUtils.cpp xml.cpp clockWC.cpp              \
	                  PrintCtrl.cpp LogPrintCtrl.cpp mdp_allo.cpp      \
			  checkFinite.cpp application.cpp ctexceptions.cpp \
	                  global.cpp

AM_LDFLAGS              = -L$(top_builddir)/build/lib/ 
LIBS                    = -lexecstream

AM_CPPFLAGS             = -I$(top_builddir) -I$(top_builddir)/include -I$(top_builddir)/ext/libexecstream
AM_CXXFLAGS             = $(AM_CPPFLAGS)

lib_LTLIBRARIES         = $(top_builddir)/build/lib/libctbase.la
library_includedir      = $(top_builddir)/build/include

otherincludedir         = $(top_builddir)/build/include/cantera/kernel
otherinclude_HEADERS    = $(h_sources)

#-----------------------
# Cantera Converters C/C++ library
#-----------------------

__top_builddir__build_lib_libctbase_la_LDFLAGS      = $(all_libraries) -release $(GENERIC_RELEASE)
__top_builddir__build_lib_libctbase_la_SOURCES      = $(cc_sources) $(h_sources)

__top_builddir__build_lib_libctbase_include_HEADERS = $(h_sources)
__top_builddir__build_lib_libctbase_includedir = $(prefix)/include

CLEANFILES = *.o
