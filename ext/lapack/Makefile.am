fc_sources = dbdsqr.f dgbtrf.f dgbtf2.f dgbtrs.f dgbsv.f  \
             dgebd2.f dgebrd.f dgelq2.f dgelqf.f dgelss.f \
             dgeqr2.f dgeqrf.f dgetf2.f dgetrf.f dgetri.f \
             dgetrs.f dlabad.f dlabrd.f dlacpy.f dlamch.f \
             dlange.f dlapy2.f dlarf.f dlarfb.f dlarfg.f  \
             dlarft.f dlartg.f dlas2.f dlascl.f dlaset.f  \
             dlasq1.f dlasq2.f dlasq3.f dlasq4.f dlasr.f  \
             dlasrt.f dlassq.f dlasv2.f dlaswp.f dorg2r.f \
             dorgbr.f dorgl2.f dorglq.f dorgqr.f dorm2r.f \
             dormbr.f dorml2.f dormlq.f dormqr.f drscl.f  \
             ilaenv.f 

INCLUDES                = 
AM_CXXFLAGS             = $(INCLUDES)
AM_FCFLAGS              = $(INCLUDES)

lib_LTLIBRARIES         = $(top_builddir)/build/lib/libctlapack.la
library_includedir      = $(top_builddir)/build/include
library_include_HEADERS = $(h_sources)

#-----------------------
# Cantera Converters C/C++ library
#-----------------------

__top_builddir__build_lib_libctlapack_la_LDFLAGS      = $(all_libraries) -release $(GENERIC_RELEASE)
__top_builddir__build_lib_libctlapack_la_SOURCES      = $(fc_sources) $(cc_sources)

CLEANFILES = *.o
