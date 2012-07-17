cc_sources = \
dasum.c \
daxpy.c \
dcabs1.c \
dcopy.c \
ddot.c \
dgbmv.c \
dgemm.c \
dgemv.c \
dger.c \
dnrm2.c \
drot.c \
drotg.c \
drotm.c \
drotmg.c \
dsbmv.c \
dscal.c \
dsdot.c \
dspmv.c \
dspr.c \
dspr2.c \
dswap.c \
dsymm.c \
dsymv.c \
dsyr.c \
dsyr2.c \
dsyr2k.c \
dsyrk.c \
dtbmv.c \
dtbsv.c \
dtpmv.c \
dtpsv.c \
dtrmm.c \
dtrmv.c \
dtrsm.c \
dtrsv.c \
dzasum.c \
dznrm2.c \
idamax.c \
lsame.c \
xerbla.c 

h_sources               = blaswrap.h

AM_CPPFLAGS             = -I../f2c_libs
AM_CXXFLAGS             = $(AM_CPPFLAGS)
AM_FCFLAGS              = $(AM_CPPFLAGS)

lib_LTLIBRARIES         = $(top_builddir)/build/lib/libctblas.la
library_includedir      = $(top_builddir)/include
library_include_HEADERS = $(h_sources)

#-----------------------
# Cantera Converters C/C++ library
#-----------------------

__top_builddir__build_lib_libctblas_la_LDFLAGS      = $(all_libraries) -release $(GENERIC_RELEASE)
__top_builddir__build_lib_libctblas_la_SOURCES      = $(cc_sources) $(h_sources)

__top_builddir__build_lib_libctblas_include_HEADERS = $(h_sources)
__top_builddir__build_lib_libctblas_includedir = $(prefix)/include

CLEANFILES = *.o
