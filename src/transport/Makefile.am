h_sources  = $(top_builddir)/include/cantera/transport/AqueousTransport.h \
             $(top_builddir)/include/cantera/transport/LiquidTransportData.h \
             $(top_builddir)/include/cantera/transport/MultiTransport.h \
             $(top_builddir)/include/cantera/transport/TransportBase.h \
             $(top_builddir)/include/cantera/transport/WaterTransport.h    \
             $(top_builddir)/include/cantera/transport/DustyGasTransport.h \
             $(top_builddir)/include/cantera/transport/LiquidTransport.h \
             $(top_builddir)/include/cantera/transport/MixTransport.h \
             $(top_builddir)/include/cantera/transport/SimpleTransport.h \
             $(top_builddir)/include/cantera/transport/TransportFactory.h \
             $(top_builddir)/include/cantera/transport/LiquidTransportParams.h \
             $(top_builddir)/include/cantera/transport/SolidTransport.h \
             $(top_builddir)/include/cantera/transport/TransportParams.h \
             $(top_builddir)/include/cantera/transport/PecosTransport.h

#             $(top_builddir)/include/cantera/transport/L_matrix.h  
#             $(top_builddir)/include/cantera/transport/FtnTransport.h
#             $(top_builddir)/include/cantera/transport/MMCollisionInt.h 
cc_sources = AqueousTransport.cpp LiquidTransport.cpp MMCollisionInt.cpp \
             SimpleTransport.cpp  TransportBase.cpp WaterTransport.cpp   \
             DustyGasTransport.cpp  MixTransport.cpp MultiTransport.cpp  \
	     SolidTransport.cpp TransportFactory.cpp LTPspecies.cpp      \
	     LiquidTranInteraction.cpp LiquidTransportData.cpp           \
             LiquidTransportParams.cpp TortuosityBase.cpp                \
             TortuosityBruggeman.cpp TortuosityMaxwell.cpp               \
             TortuosityPercolation.cpp TransportParams.cpp               \
	     GasTransport.cpp PecosTransport.cpp

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
__top_builddir__build_lib_libtransport_includedir = $(prefix)/include/cantera/transport
