h_sources  = 
cc_sources = AqueousTransport.cpp LiquidTransport.cpp MMCollisionInt.cpp \
             SimpleTransport.cpp  TransportBase.cpp WaterTransport.cpp   \
             DustyGasTransport.cpp  MixTransport.cpp MultiTransport.cpp  \
	     SolidTransport.cpp TransportFactory.cpp

AM_CPPFLAGS             = -I$(top_builddir)/include
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
