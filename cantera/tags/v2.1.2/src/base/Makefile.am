h_sources               = $(top_builddir)/include/cantera/base/utilities.h \
                          $(top_builddir)/include/cantera/base/ct_defs.h \
                          $(top_builddir)/include/cantera/base/ctexceptions.h \
                          $(top_builddir)/include/cantera/base/logger.h  \
                          $(top_builddir)/include/cantera/base/ctml.h \
                          $(top_builddir)/include/cantera/base/plots.h \
                          $(top_builddir)/include/cantera/base/stringUtils.h \
                          $(top_builddir)/include/cantera/base/xml.h  \
			  $(top_builddir)/include/cantera/base/Array.h \
                          $(top_builddir)/include/cantera/base/vec_functions.h \
                          $(top_builddir)/include/cantera/base/global.h \
                          $(top_builddir)/include/cantera/base/XML_Writer.h  \
                          $(top_builddir)/include/cantera/base/FactoryBase.h \
                          $(top_builddir)/include/cantera/base/clockWC.h \
	                  $(top_builddir)/include/cantera/base/PrintCtrl.h \
			  LogPrintCtrl.h units.h application.h

cc_sources              = ct2ctml.cpp ctml.cpp plots.cpp                   \
			  stringUtils.cpp xml.cpp clockWC.cpp              \
	                  PrintCtrl.cpp LogPrintCtrl.cpp       \
			  checkFinite.cpp application.cpp ctexceptions.cpp \
	                  global.cpp

AM_LDFLAGS              = -L$(top_builddir)/build/lib/ 
LIBS                    = -L$(top_builddir)/build/lib/ -lexecstream

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
__top_builddir__build_lib_libctbase_includedir = $(prefix)/include/cantera/base

CLEANFILES = *.o
