SET (TRANSPORT_SRCS  TransportFactory.cpp MultiTransport.cpp MixTransport.cpp
                                      MMCollisionInt.cpp SolidTransport.cpp DustyGasTransport.cpp)

INCLUDE_DIRECTORIES (${PROJECT_SOURCE_DIR})
INCLUDE_DIRECTORIES (${PROJECT_SOURCE_DIR}/Cantera/src/base)
INCLUDE_DIRECTORIES (${PROJECT_SOURCE_DIR}/Cantera/src/thermo)
INCLUDE_DIRECTORIES (${PROJECT_SOURCE_DIR}/Cantera/src/numerics)

ADD_LIBRARY(transport ${TRANSPORT_SRCS})

TARGET_LINK_LIBRARIES (transport numerics thermo ctbase)

SET (TRANSPORT_H   TransportFactory.h MultiTransport.h MixTransport.h 
       		       MMCollisionInt.h  SolidTransport.h DustyGasTransport.h 
           	       TransportBase.h L_matrix.h TransportParams.h )
INSTALL_FILES(/include/cantera/kernel FILES ${TRANSPORT_H})


