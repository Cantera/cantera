h_sources  = 
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
	     vcs_dbocls.c vcs_dbolsm.c vcs_dbols.c vcs_dvout.c \
	     vcs_ivout.c vcs_dmout.c vcs_xerror.c

AM_CPPFLAGS    = -I$(top_builddir)/include -I$(top_builddir)/ext/f2c_libs/
AM_CXXFLAGS = $(AM_CPPFLAGS)

lib_LTLIBRARIES         = $(top_builddir)/build/lib/libequil.la
library_includedir      = $(includedir)

#-----------------------
# Cantera Converters C/C++ library
#-----------------------

__top_builddir__build_lib_libequil_la_LDFLAGS      = $(all_libraries) -release $(GENERIC_RELEASE)
__top_builddir__build_lib_libequil_la_SOURCES      = $(cc_sources) $(h_sources)

__top_builddir__build_lib_libequil_include_HEADERS = $(h_sources)
__top_builddir__build_lib_libequil_includedir = $(prefix)/include
