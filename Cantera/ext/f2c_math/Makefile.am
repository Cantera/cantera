cc_sources = \
mach.cpp \
ddaspk.c \
dgbefa.c \
dgbsl.c \
dgefa.c \
dgesl.c \
dp1vlu.c \
dpcoef.c \
dpolft.c \
fdump.c \
j4save.c \
pcoef.c \
polfit.c \
pvalue.c \
xercnt.c \
xerhlt.c \
xermsg.c \
xerprn.c \
xersve.c \
xgetua.c \
printstring.c

h_sources               = cblas.h gmres.h mkl_cblas.h

AM_CPPFLAGS             = -I../f2c_libs
AM_CXXFLAGS             = $(AM_CPPFLAGS)
AM_FCFLAGS              = $(AM_CPPFLAGS)

lib_LTLIBRARIES         = $(top_builddir)/build/lib/libctmath.la
library_includedir      = $(top_builddir)/build/include
library_include_HEADERS = $(h_sources)

#-----------------------
# Cantera Converters C/C++ library
#-----------------------

__top_builddir__build_lib_libctmath_la_LDFLAGS      = $(all_libraries) -release $(GENERIC_RELEASE)
__top_builddir__build_lib_libctmath_la_SOURCES      = $(cc_sources) $(h_sources)

CLEANFILES = *.o

# header file accumulation
all:
	@echo copying headers $<
	cp *.h $(top_builddir)/build/include/
