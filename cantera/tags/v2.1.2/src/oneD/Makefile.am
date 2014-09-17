h_sources  = $(top_builddir)/include/cantera/oneD/Inlet1D.h \
             $(top_builddir)/include/cantera/oneD/MultiJac.h \
             $(top_builddir)/include/cantera/oneD/Sim1D.h \
             $(top_builddir)/include/cantera/oneD/StFlow.h \
             $(top_builddir)/include/cantera/oneD/Domain1D.h \
             $(top_builddir)/include/cantera/oneD/MultiNewton.h \
             $(top_builddir)/include/cantera/oneD/OneDim.h \
             $(top_builddir)/include/cantera/oneD/refine.h

#             $(top_builddir)/include/cantera/oneD/Surf1D.h  
#             $(top_builddir)/include/cantera/oneD/Resid1D.h
#             $(top_builddir)/include/cantera/oneD/Solid1D.h

cc_sources = MultiJac.cpp MultiNewton.cpp newton_utils.cpp OneDim.cpp \
             StFlow.cpp boundaries1D.cpp refine.cpp Sim1D.cpp Domain1D.cpp

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
__top_builddir__build_lib_liboneD_includedir = $(prefix)/include/cantera/oneD
