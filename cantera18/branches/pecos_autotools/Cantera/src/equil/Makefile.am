h_sources = vcs_internal.h vcs_VolPhase.h vcs_solve.h vcs_prob.h \
            vcs_IntStarStar.h vcs_DoubleStarStar.h vcs_defs.h    \
            vcs_MultiPhaseEquil.h vcs_Exception.h                \
            vcs_SpeciesProperties.h vcs_species_thermo.h

cc_sources = BasisOptimize.cpp vcs_elem.cpp vcs_IntStarStar.cpp             \
             vcs_prob.cpp vcs_solve.cpp vcs_VolPhase.cpp ChemEquil.cpp      \
             vcs_elem_rearrange.cpp vcs_linmaxc.cpp vcs_rearrange.cpp       \
	     vcs_solve_TP.cpp equilibrate.cpp vcs_equilibrate.cpp           \
	     vcs_MultiPhaseEquil.cpp vcs_report.cpp                         \
	     vcs_SpeciesProperties.cpp MultiPhase.cpp vcs_Exception.cpp     \
	     vcs_nondim.cpp vcs_root1d.cpp vcs_species_thermo.cpp           \
	     MultiPhaseEquil.cpp vcs_Gibbs.cpp vcs_phaseStability.cpp       \
             vcs_rxnadj.cpp vcs_TP.cpp vcs_DoubleStarStar.cpp vcs_inest.cpp \
	     vcs_prep.cpp vcs_setMolesLinProg.cpp vcs_util.cpp

cc_sources += vcs_dbocls.c vcs_dbols.c vcs_dbolsm.c vcs_dmout.c            \
              vcs_dvout.c vcs_ivout.c vcs_xerror.c                          

INCLUDES    = -I../base -I../thermo -I../numerics -I../../../ 
INCLUDES   += -I../../../ext/f2c_libs/
AM_CXXFLAGS = $(INCLUDES)

lib_LTLIBRARIES         = libequil.la
library_includedir      = $(includedir)

#-----------------------
# Cantera Converters C/C++ library
#-----------------------

libequil_la_LDFLAGS      = $(all_libraries) -release $(GENERIC_RELEASE)
libequil_la_SOURCES      = $(cc_sources) $(h_sources)