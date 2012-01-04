h_sources               = utilities.h ct_defs.h ctexceptions.h logger.h  \
                          ctml.h plots.h stringUtils.h xml.h             \
			  Array.h vec_functions.h global.h XML_Writer.h  \
                          FactoryBase.h clockWC.h PrintCtrl.h            \
	 		  LogPrintCtrl.h  mdp_allo.h config.h   

cc_sources              = ct2ctml.cpp ctml.cpp misc.cpp plots.cpp     \
			  stringUtils.cpp xml.cpp clockWC.cpp         \
	                  PrintCtrl.cpp LogPrintCtrl.cpp mdp_allo.cpp \
			  checkFinite.cpp	

INCLUDES                = 
AM_CXXFLAGS             = $(INCLUDES)

lib_LTLIBRARIES         = $(top_builddir)/build/lib/libctbase.la
library_includedir      = $(includedir)

#-----------------------
# Cantera Converters C/C++ library
#-----------------------

__top_builddir__build_lib_libctbase_la_LDFLAGS      = $(all_libraries) -release $(GENERIC_RELEASE)
__top_builddir__build_lib_libctbase_la_SOURCES      = $(cc_sources) $(h_sources)
