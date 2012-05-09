h_sources         = Cabinet.h clib_defs.h ctbdry.h ctfunc.h ct.h ctmultiphase.h \
                    ctnum.h ctonedim.h ctreactor.h ctrpath.h ctstagn.h ctsurf.h \
                    ctxml.h Storage.h

cc_sources        = ctbdry.cpp ct.cpp ctfunc.cpp ctmultiphase.cpp     \
                    ctonedim.cpp ctreactor.cpp ctrpath.cpp            \
                    ctsurf.cpp ctxml.cpp Storage.cpp

INC                     = -I./ -I$(top_builddir)/build/include/
AM_CPPFLAGS             = $(INC)
AM_CXXFLAGS             = $(AM_CPPFLAGS)
AM_FCFLAGS              = $(AM_CPPFLAGS)

LINK                    = -luser -loneD -lzeroD -lequil -lkinetics -ltransport -lthermo
LINK                   += -lctnumerics -lctmath -ltpx -lctspectra -lconverters -lctbase -lcvode
LINK                   += -lctlapack -lctblas -lctf2c -lctcxx -lm  -lctf2c -lstdc++

AM_LDFLAGS              = -L$(top_builddir)/build/lib/

lib_LTLIBRARIES         = $(top_builddir)/build/lib/libclib.la
library_includedir      = $(top_builddir)/build/include
library_include_HEADERS = $(h_sources)

#-----------------------
# Cantera Converters C/C++ library
#-----------------------

__top_builddir__build_lib_libclib_la_LDFLAGS      = $(all_libraries) -release $(GENERIC_RELEASE)
__top_builddir__build_lib_libclib_la_SOURCES      = $(cc_sources) $(h_sources)

__top_builddir__build_lib_libclib_include_HEADERS = $(h_sources)
__top_builddir__build_lib_libclib_includedir      = $(prefix)/include

# header file accumulation
all-local:
	@echo copying headers $<
	cp -p *.h $(top_builddir)/build/include/
