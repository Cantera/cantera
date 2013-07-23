cc_sources = dbdsqr.c dgbtrf.c dgbtf2.c dgbtrs.c dgbsv.c dgebd2.c dgebrd.c  \
             dgelq2.c dgelqf.c dgelss.c dgeqr2.c dgeqrf.c dgetf2.c dgetrf.c \
             dgetri.c dgetrs.c dlabad.c dlabrd.c dlacpy.c dlamch.c dlange.c \
             dlapy2.c dlarf.c dlarfb.c dlarfg.c dlarft.c dlartg.c dlas2.c   \
             dlascl.c dlaset.c dlasq1.c dlasq2.c dlasq3.c dlasq4.c dlasq5.c \
             dlasq6.c dlasr.c dlasrt.c dlassq.c dlasv2.c dlaswp.c dorg2r.c  \
             dorgbr.c dorgl2.c dorglq.c dorgqr.c dorm2r.c dormbr.c dorml2.c \
             dormlq.c dormqr.c drscl.c dtrtri.c dtrti2.c ieeeck.c ilaenv.c  \
	     dpotrs.c dtrtrs.c dgbcon.c dgecon.c dlatrs.c dlatbs.c dlacon.c \
             dtrcon.c dlantr.c dpotrf.c dpotf2.c

h_sources               = blaswrap.h

AM_CPPFLAGS             = -I../f2c_libs
AM_CXXFLAGS             = $(AM_CPPFLAGS)
AM_FCFLAGS              = $(AM_CPPFLAGS)

lib_LTLIBRARIES         = $(top_builddir)/build/lib/libctlapack.la
library_includedir      = $(top_builddir)/include
library_include_HEADERS = $(h_sources)

#-----------------------
# Cantera Converters C/C++ library
#-----------------------

__top_builddir__build_lib_libctlapack_la_LDFLAGS      = $(all_libraries) -release $(GENERIC_RELEASE)
__top_builddir__build_lib_libctlapack_la_SOURCES      = $(cc_sources) $(h_sources)

__top_builddir__build_lib_libctlapack_include_HEADERS = $(h_sources)
__top_builddir__build_lib_libctlapack_includedir = $(prefix)/include


CLEANFILES = *.o

# header file accumulation
#all:
#	@echo copying headers $<
#	cp -p *.h $(top_builddir)/build/include/
