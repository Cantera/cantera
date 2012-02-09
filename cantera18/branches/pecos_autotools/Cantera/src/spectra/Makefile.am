h_sources  = DiatomicMolecule.h LineBroadener.h Nuclei.h \
             rotor.h spectralUtilities.h

cc_sources = LineBroadener.cpp rotor.cpp spectralUtilities.cpp

AM_CPPFLAGS             = -I../base
AM_CXXFLAGS             = $(AM_CPPFLAGS)

lib_LTLIBRARIES         = $(top_builddir)/build/lib/libctspectra.la
library_includedir      = $(includedir)

#-----------------------
# Cantera Spectra C/C++ library
#-----------------------

__top_builddir__build_lib_libctspectra_la_LDFLAGS      = $(all_libraries) -release $(GENERIC_RELEASE)
__top_builddir__build_lib_libctspectra_la_SOURCES      = $(cc_sources) $(h_sources)


# header file accumulation
all:
	@echo copying headers $<
	cp -p $(h_sources) $(top_builddir)/build/include/
