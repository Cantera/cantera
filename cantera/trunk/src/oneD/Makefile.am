cc_sources = MultiJac.cpp MultiNewton.cpp newton_utils.cpp OneDim.cpp \
             StFlow.cpp boundaries1D.cpp refine.cpp Sim1D.cpp Domain1D.cpp

h_sources  = 
AM_CPPFLAGS             = -I$(top_builddir)/include
AM_CXXFLAGS             = $(AM_CPPFLAGS) -fPIC

lib_LTLIBRARIES         = $(top_builddir)/build/lib/liboneD.la
library_includedir      = $(includedir)

#-----------------------
# Cantera OneD C/C++ library
#-----------------------

__top_builddir__build_lib_liboneD_la_LDFLAGS      = $(all_libraries) -release $(GENERIC_RELEASE)
__top_builddir__build_lib_liboneD_la_SOURCES      = $(cc_sources) $(h_sources)

__top_builddir__build_lib_liboneD_include_HEADERS = $(h_sources)
__top_builddir__build_lib_liboneD_includedir = $(prefix)/include
