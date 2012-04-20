h_sources  = AqueousTransport.h LiquidTransportData.h L_matrix.h  \
	     MultiTransport.h TransportBase.h WaterTransport.h    \
             DustyGasTransport.h LiquidTransport.h MixTransport.h \
	     SimpleTransport.h  TransportFactory.h FtnTransport.h \
	     LiquidTransportParams.h  MMCollisionInt.h  	  \
             SolidTransport.h TransportParams.h PecosTransport.h

cc_sources = AqueousTransport.cpp LiquidTransport.cpp MMCollisionInt.cpp \
             SimpleTransport.cpp  TransportBase.cpp WaterTransport.cpp   \
             DustyGasTransport.cpp  MixTransport.cpp MultiTransport.cpp  \
	     SolidTransport.cpp TransportFactory.cpp PecosTransport.cpp

AM_CPPFLAGS             = -I../base -I../thermo -I../numerics
AM_CXXFLAGS             = $(AM_CPPFLAGS)

lib_LTLIBRARIES         = $(top_builddir)/build/lib/libtransport.la
library_includedir      = $(includedir)

#-----------------------
# Cantera Transport C/C++ library
#-----------------------

__top_builddir__build_lib_libtransport_la_LDFLAGS      = $(all_libraries) -release $(GENERIC_RELEASE)
__top_builddir__build_lib_libtransport_la_SOURCES      = $(cc_sources) $(h_sources)

__top_builddir__build_lib_libtransport_include_HEADERS = $(h_sources)
__top_builddir__build_lib_libtransport_includedir = $(prefix)/include


# header file accumulation
all:
	@echo copying headers $<
	cp -p *.h $(top_builddir)/build/include/
