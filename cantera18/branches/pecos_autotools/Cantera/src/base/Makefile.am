h_sources               = utilities.h ct_defs.h ctexceptions.h logger.h  \
                          ctml.h plots.h stringUtils.h xml.h             \
			  Array.h vec_functions.h global.h XML_Writer.h  \
                          FactoryBase.h clockWC.h PrintCtrl.h            \
	 		  LogPrintCtrl.h mdp_allo.h units.h

cc_sources              = ct2ctml.cpp ctml.cpp misc.cpp plots.cpp     \
			  stringUtils.cpp xml.cpp clockWC.cpp         \
	                  PrintCtrl.cpp LogPrintCtrl.cpp mdp_allo.cpp \
			  checkFinite.cpp	

AM_CPPFLAGS             = -I$(top_builddir)
AM_CXXFLAGS             = $(AM_CPPFLAGS)

lib_LTLIBRARIES         = $(top_builddir)/build/lib/libctbase.la
library_includedir      = $(top_builddir)/build/include
library_include_HEADERS = utilities.h

otherincludedir         = $(top_builddir)/build/include/cantera/kernel
otherinclude_HEADERS    = $(h_sources)

#-----------------------
# Cantera Converters C/C++ library
#-----------------------

__top_builddir__build_lib_libctbase_la_LDFLAGS      = $(all_libraries) -release $(GENERIC_RELEASE)
__top_builddir__build_lib_libctbase_la_SOURCES      = $(cc_sources) $(h_sources)

__top_builddir__build_lib_libctbase_include_HEADERS = $(h_sources)
__top_builddir__build_lib_libctbase_includedir = $(prefix)/include

nodist___top_builddir__build_lib_libctbase_include_HEADERS = $(top_srcdir)/config.h

# header file accumulation -- only copy if file updated!
all:
	@echo copying headers $<
	cp -p *.h $(top_builddir)/build/include/

CLEANFILES = *.o
