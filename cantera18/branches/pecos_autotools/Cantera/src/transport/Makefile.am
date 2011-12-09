h_sources  = AqueousTransport.h LiquidTransportData.h L_matrix.h  \
	     MultiTransport.h TransportBase.h WaterTransport.h    \
             DustyGasTransport.h LiquidTransport.h MixTransport.h \
	     SimpleTransport.h  TransportFactory.h FtnTransport.h \
	     LiquidTransportParams.h  MMCollisionInt.h  	  \
             SolidTransport.h   TransportParams.h

cc_sources = AqueousTransport.cpp LiquidTransport.cpp MMCollisionInt.cpp \
             SimpleTransport.cpp  TransportBase.cpp WaterTransport.cpp   \
             DustyGasTransport.cpp  MixTransport.cpp MultiTransport.cpp  \
	     SolidTransport.cpp   TransportFactory.cpp 

INCLUDES                = -I../base -I../thermo -I../numerics -I../../../
AM_CXXFLAGS             = $(INCLUDES)

lib_LTLIBRARIES         = libtransport.la
library_includedir      = $(includedir)

#-----------------------
# Cantera Transport C/C++ library
#-----------------------

libtransport_la_LDFLAGS      = $(all_libraries) -release $(GENERIC_RELEASE)
libtransport_la_SOURCES      = $(cc_sources) $(h_sources)