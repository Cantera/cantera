SET (ONED_SRCS  MultiJac.cpp MultiNewton.cpp newton_utils.cpp OneDim.cpp
          StFlow.cpp boundaries1D.cpp refine.cpp Sim1D.cpp Domain1D.cpp )

INCLUDE_DIRECTORIES (${PROJECT_SOURCE_DIR})
INCLUDE_DIRECTORIES (${PROJECT_SOURCE_DIR}/Cantera/src/base)
INCLUDE_DIRECTORIES (${PROJECT_SOURCE_DIR}/Cantera/src/thermo)
INCLUDE_DIRECTORIES (${PROJECT_SOURCE_DIR}/Cantera/src/numerics)
INCLUDE_DIRECTORIES (${PROJECT_SOURCE_DIR}/Cantera/src/kinetics)
INCLUDE_DIRECTORIES (${PROJECT_SOURCE_DIR}/Cantera/src/transport)

ADD_LIBRARY(oneD  ${ONED_SRCS} )

SET (ONED_H  Inlet1D.h   MultiJac.h  Sim1D.h     StFlow.h
          Surf1D.h    Domain1D.h  MultiNewton.h   OneDim.h
          Resid1D.h   Solid1D.h   refine.h)

INSTALL_FILES(/include/cantera/kernel FILES ${ONED_H})

