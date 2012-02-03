cc_sources = MultiJac.cpp MultiNewton.cpp newton_utils.cpp OneDim.cpp \
             StFlow.cpp boundaries1D.cpp refine.cpp Sim1D.cpp Domain1D.cpp

h_sources  = Inlet1D.h MultiJac.h Sim1D.h StFlow.h \
	     Surf1D.h  Domain1D.h MultiNewton.h OneDim.h \
	     Resid1D.h Solid1D.h refine.h

AM_CPPFLAGS             = -I../base -I../thermo -I../kinetics -I../transport -I../numerics -I../../../
AM_CXXFLAGS             = $(AM_CPPFLAGS) -fPIC

lib_LTLIBRARIES         = $(top_builddir)/build/lib/liboneD.la
library_includedir      = $(includedir)

#-----------------------
# Cantera OneD C/C++ library
#-----------------------

__top_builddir__build_lib_liboneD_la_LDFLAGS      = $(all_libraries) -release $(GENERIC_RELEASE)
__top_builddir__build_lib_liboneD_la_SOURCES      = $(cc_sources) $(h_sources)


# header file accumulation
all:
	@echo copying headers $<
	cp $(h_sources) $(top_builddir)/build/include/
