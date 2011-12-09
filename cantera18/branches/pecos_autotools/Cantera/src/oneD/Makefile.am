h_sources  = Domain1D.h  Inlet1D.h  MultiJac.h  MultiNewton.h  \
	     OneDim.h  refine.h  Resid1D.h  Sim1D.h  Solid1D.h \
	     StFlow.h  Surf1D.h

cc_sources = boundaries1D.cpp  Domain1D.cpp  MultiJac.cpp       \
	     MultiNewton.cpp  newton_utils.cpp  oneD_files.cpp  \
	     OneDim.cpp  refine.cpp  Sim1D.cpp  Solid1D.cpp     \
	     StFlow.cpp

INCLUDES                = -I../base -I../thermo -I../kinetics -I../transport -I../numerics -I../../../
AM_CXXFLAGS             = $(INCLUDES) -fPIC

lib_LTLIBRARIES         = liboneD.la
library_includedir      = $(includedir)

#-----------------------
# Cantera OneD C/C++ library
#-----------------------

liboneD_la_LDFLAGS      = $(all_libraries) -release $(GENERIC_RELEASE)
liboneD_la_SOURCES      = $(cc_sources) $(h_sources)