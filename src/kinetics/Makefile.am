
h_sources  = $(top_builddir)/include/cantera/kinetics/importKinetics.h \
             $(top_builddir)/include/cantera/kinetics/GRI_30_Kinetics.h \
             $(top_builddir)/include/cantera/kinetics/KineticsFactory.h \
             $(top_builddir)/include/cantera/kinetics/Kinetics.h \
             $(top_builddir)/include/cantera/kinetics/GasKinetics.h  \
             $(top_builddir)/include/cantera/kinetics/FalloffFactory.h \
             $(top_builddir)/include/cantera/kinetics/ReactionStoichMgr.h \
             $(top_builddir)/include/cantera/kinetics/reaction_defs.h \
             $(top_builddir)/include/cantera/kinetics/FalloffMgr.h  \
             $(top_builddir)/include/cantera/kinetics/ThirdBodyMgr.h  \
             $(top_builddir)/include/cantera/kinetics/RateCoeffMgr.h \
             $(top_builddir)/include/cantera/kinetics/ReactionData.h \
             $(top_builddir)/include/cantera/kinetics/RxnRates.h \
             $(top_builddir)/include/cantera/kinetics/Enhanced3BConc.h \
             $(top_builddir)/include/cantera/kinetics/StoichManager.h \
             $(top_builddir)/include/cantera/kinetics/AqueousKinetics.h \
             $(top_builddir)/include/cantera/kinetics/InterfaceKinetics.h  \
             $(top_builddir)/include/cantera/kinetics/Group.h \
             $(top_builddir)/include/cantera/kinetics/ReactionPath.h \
             $(top_builddir)/include/cantera/kinetics/EdgeKinetics.h

#             $(top_builddir)/include/cantera/kinetics/solveSP.h 
#             $(top_builddir)/include/cantera/kinetics/ImplicitSurfChem.h 

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
__top_builddir__build_lib_libkinetics_includedir = $(prefix)/include/cantera/base