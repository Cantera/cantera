h_sources = vcs_internal.h vcs_VolPhase.h vcs_solve.h vcs_prob.h \
            vcs_IntStarStar.h    vcs_DoubleStarStar.h vcs_defs.h \
            vcs_MultiPhaseEquil.h  vcs_Exception.h  \
            vcs_SpeciesProperties.h vcs_species_thermo.h

cc_sources = vcs_solve_TP.o vcs_VolPhase.o vcs_solve.o vcs_prob.o \
             vcs_TP.o vcs_report.o vcs_util.o \
             vcs_IntStarStar.o  vcs_DoubleStarStar.o  vcs_elem.o \
             vcs_elem_rearrange.o   vcs_MultiPhaseEquil.o \
             vcs_nondim.o vcs_Exception.o \
             vcs_inest.o vcs_rearrange.o \
             vcs_root1d.o vcs_rxnadj.o \
             vcs_SpeciesProperties.o vcs_equilibrate.o \
             vcs_prep.o vcs_species_thermo.o vcs_Gibbs.o vcs_phaseStabil

INCLUDES                = = -I../base -I../thermo -I../numerics
AM_CXXFLAGS             = $(INCLUDES)

lib_LTLIBRARIES         = libequil.la
library_includedir      = $(includedir)

#-----------------------
# Cantera Converters C/C++ library
#-----------------------

libequil_la_LDFLAGS      = $(all_libraries) -release $(GENERIC_RELEASE)
libequil_la_SOURCES      = $(cc_sources) $(h_sources)

# Basic Cantera Thermodynamics Object Files
EQUIL_OBJ = BasisOptimize.o ChemEquil.o MultiPhase.o MultiPhaseEquil.o equilibrate.o

EQUIL_H = ChemEquil.h MultiPhase.h MultiPhaseEquil.h equil.h PropertyCalculator.h
