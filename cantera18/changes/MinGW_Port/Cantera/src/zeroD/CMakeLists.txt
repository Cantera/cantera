SET (ZEROD_SRCS Reactor.cpp ReactorBase.cpp FlowDevice.cpp Wall.cpp ReactorNet.cpp
      FlowReactor.cpp ConstPressureReactor.cpp ReactorFactory.cpp )

INCLUDE_DIRECTORIES (${PROJECT_SOURCE_DIR})
INCLUDE_DIRECTORIES (${PROJECT_SOURCE_DIR}/Cantera/src/base)
INCLUDE_DIRECTORIES (${PROJECT_SOURCE_DIR}/Cantera/src/thermo)
INCLUDE_DIRECTORIES (${PROJECT_SOURCE_DIR}/Cantera/src/numerics)
INCLUDE_DIRECTORIES (${PROJECT_SOURCE_DIR}/Cantera/src/kinetics)

ADD_LIBRARY(zeroD  ${ZEROD_SRCS} )

SET (ZEROD_H Reactor.h ReactorBase.h FlowDevice.h Wall.h ReactorNet.h
          flowControllers.h Reservoir.h FlowReactor.h
          ConstPressureReactor.h ReactorFactory.h )
INSTALL_FILES(/include/cantera/kernel FILES ${ZEROD_H})

