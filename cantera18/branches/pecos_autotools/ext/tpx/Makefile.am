h_sources               = Methane.h Nitrogen.h Oxygen.h Water.h Hydrogen.h RedlichKwong.h \
		          CarbonDioxide.h Heptane.h lk.h Sub.h utils.h HFC134a.h

cc_sources              = Methane.cpp Nitrogen.cpp Oxygen.cpp Water.cpp Hydrogen.cpp RedlichKwong.cpp \
		          CarbonDioxide.cpp Heptane.cpp lk.cpp Sub.cpp utils.cpp HFC134a.cpp

INCLUDES                = 
AM_CXXFLAGS             = $(INCLUDES)

lib_LTLIBRARIES         = $(top_builddir)/build/lib/libtpx.la
library_includedir      = $(top_builddir)/build/include
library_include_HEADERS = $(h_sources)

#-----------------------
# Cantera Converters C/C++ library
#-----------------------

__top_builddir__build_lib_libtpx_la_LDFLAGS      = $(all_libraries) -release $(GENERIC_RELEASE)
__top_builddir__build_lib_libtpx_la_SOURCES      = $(cc_sources) $(h_sources)

CLEANFILES = *.o
