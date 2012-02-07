h_sources  = ChemEquil.h       PropertyCalculator.h  vcs_internal.h  \
             equil.h            vcs_defs.h             vcs_IntStarStar.h      vcs_SpeciesProperties.h \
             MultiPhaseEquil.h  vcs_DoubleStarStar.h   vcs_MultiPhaseEquil.h  vcs_species_thermo.h \
             MultiPhase.h       vcs_Exception.h        vcs_prob.h             vcs_VolPhase.h \
	     vcs_solve.h

#ChemEquil.h MultiPhase.h MultiPhaseEquil.h equil.h        \
	     PropertyCalculator.h vcs_equilibrate.h                    \
             vcs_MultiPhaseEquil.h vcs_VolPhase.h vcs_DoubleStarStar.h \
             vcs_prob.h vcs_SpeciesProperties.h vcs_internal.h         \
             vcs_species_thermo.h vcs_TP.h vcs_util.h vcs_VolPhase.h   \
             vcs_MultiPhaseEquil.h

cc_sources = BasisOptimize.cpp       vcs_Gibbs.cpp            vcs_root1d.cpp \
	     ChemEquil.cpp           vcs_rxnadj.cpp           vcs_TP.cpp \
             equilibrate.cpp         vcs_IntStarStar.cpp      vcs_setMolesLinProg.cpp \
             MultiPhase.cpp          vcs_MultiPhaseEquil.cpp  vcs_solve.cpp \
             MultiPhaseEquil.cpp     vcs_nondim.cpp           vcs_solve_TP.cpp \
             vcs_DoubleStarStar.cpp  vcs_phaseStability.cpp   vcs_SpeciesProperties.cpp \
             vcs_elem.cpp            vcs_prep.cpp             vcs_species_thermo.cpp \
             vcs_elem_rearrange.cpp  vcs_prob.cpp             vcs_inest.cpp \
             vcs_equilibrate.cpp     vcs_rearrange.cpp        vcs_util.cpp \
             vcs_Exception.cpp       vcs_report.cpp           vcs_VolPhase.cpp \
	     vcs_linmaxc.cpp vcs_dbocls.c vcs_dbolsm.c vcs_dbols.c vcs_dvout.c \
	     vcs_ivout.c vcs_dmout.c vcs_xerror.c

#BasisOptimize.cpp ChemEquil.cpp MultiPhase.cpp            \
             MultiPhaseEquil.cpp equilibrate.cpp vcs_equilibrate.cpp   \
	     vcs_MultiPhaseEquil.cpp vcs_VolPhase.cpp                  \
             vcs_DoubleStarStar.cpp vcs_prob.cpp                       \
	     vcs_SpeciesProperties.cpp vcs_MultiPhaseEquil.cpp         \
             vcs_species_thermo.cpp vcs_TP.cpp vcs_util.cpp            \
             vcs_VolPhase.cpp

AM_CPPFLAGS    = -I$(top_builddir)/build/include
AM_CXXFLAGS = $(AM_CPPFLAGS)

lib_LTLIBRARIES         = $(top_builddir)/build/lib/libequil.la
library_includedir      = $(includedir)

#-----------------------
# Cantera Converters C/C++ library
#-----------------------

__top_builddir__build_lib_libequil_la_LDFLAGS      = $(all_libraries) -release $(GENERIC_RELEASE)
__top_builddir__build_lib_libequil_la_SOURCES      = $(cc_sources) $(h_sources)

# header file accumulation
all:
	@echo copying headers $<
	cp -u *.h $(top_builddir)/build/include/
