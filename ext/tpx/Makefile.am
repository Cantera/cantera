h_sources               = Methane.h Nitrogen.h Oxygen.h Water.h Hydrogen.h RedlichKwong.h \
		          CarbonDioxide.h Heptane.h lk.h Sub.h utils.h HFC134a.h

cc_sources              = Methane.cpp Nitrogen.cpp Oxygen.cpp Water.cpp Hydrogen.cpp RedlichKwong.cpp \
		          CarbonDioxide.cpp Heptane.cpp lk.cpp Sub.cpp utils.cpp HFC134a.cpp

AM_CPPFLAGS             = 
AM_CXXFLAGS             = $(AM_CPPFLAGS)

lib_LTLIBRARIES         = $(top_builddir)/build/lib/libtpx.la
library_includedir      = $(top_builddir)/build/include
library_include_HEADERS = $(h_sources)

#-----------------------
# Cantera Converters C/C++ library
#-----------------------

__top_builddir__build_lib_libtpx_la_LDFLAGS      = $(all_libraries) -release $(GENERIC_RELEASE)
__top_builddir__build_lib_libtpx_la_SOURCES      = $(cc_sources) $(h_sources)

__top_builddir__build_lib_libtpx_include_HEADERS = $(h_sources)
__top_builddir__build_lib_libtpx_includedir = $(prefix)/include


CLEANFILES = *.o

# header file accumulation
all:
	@echo copying headers $<
	cp -p *.h $(top_builddir)/build/include/
