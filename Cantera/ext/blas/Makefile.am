fc_sources              = dasum.f dcabs1.f ddot.f dgemm.f dger.f drot.f drotm.f \
                          dsbmv.f dsdot.f dspr2.f dswap.f dsymv.f dsyr2k.f      \
                          dsyrk.f dtbsv.f dtpsv.f dtrmv.f dtrsv.f dznrm2.f      \
                          idamax.f izamax.f xerbla.f daxpy.f dcopy.f dgbmv.f    \
                          dgemv.f dnrm2.f drotg.f drotmg.f dscal.f dspmv.f      \
                          dspr.f dsymm.f dsyr2.f dsyr.f dtbmv.f dtpmv.f dtrmm.f \
                          dtrsm.f dzasum.f icamax.f isamax.f lsame.f

AM_CPPFLAGS             = 
AM_CXXFLAGS             = $(AM_CPPFLAGS)
AM_FCFLAGS              = $(AM_CPPFLAGS)

lib_LTLIBRARIES         = $(top_builddir)/build/lib/libctblas.la
library_includedir      = $(top_builddir)/build/include
library_include_HEADERS = $(h_sources)

#-----------------------
# Cantera Converters C/C++ library
#-----------------------

__top_builddir__build_lib_libctblas_la_LDFLAGS      = $(all_libraries) -release $(GENERIC_RELEASE)
__top_builddir__build_lib_libctblas_la_SOURCES      = $(fc_sources) $(cc_sources)

CLEANFILES = *.o
