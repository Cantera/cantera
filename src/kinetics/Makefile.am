h_sources  = 
cc_sources = AqueousKinetics.cpp FalloffFactory.cpp GRI_30_Kinetics.cpp \
             GasKinetics.cpp Group.cpp ImplicitSurfChem.cpp InterfaceKinetics.cpp \
             Kinetics.cpp KineticsFactory.cpp ReactionPath.cpp ReactionStoichMgr.cpp \
             importKinetics.cpp solveSP.cpp


AM_CPPFLAGS             = -I$(top_builddir)/include/
AM_CXXFLAGS             = $(AM_CPPFLAGS)

lib_LTLIBRARIES         = $(top_builddir)/build/lib/libkinetics.la
library_includedir      = $(includedir)

#-----------------------
# Cantera Kinetics C/C++ library
#-----------------------

__top_builddir__build_lib_libkinetics_la_LDFLAGS      = $(all_libraries) -release $(GENERIC_RELEASE)
__top_builddir__build_lib_libkinetics_la_SOURCES      = $(cc_sources) $(h_sources)

__top_builddir__build_lib_libkinetics_include_HEADERS = $(h_sources)
__top_builddir__build_lib_libkinetics_includedir = $(prefix)/include