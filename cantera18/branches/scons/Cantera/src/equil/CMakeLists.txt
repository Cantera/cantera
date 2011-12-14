SET (EQUIL_SRCS BasisOptimize.cpp ChemEquil.cpp MultiPhase.cpp MultiPhaseEquil.cpp
       equilibrate.cpp )

INCLUDE_DIRECTORIES (${PROJECT_SOURCE_DIR})
INCLUDE_DIRECTORIES (${PROJECT_SOURCE_DIR}/Cantera/src/base)
INCLUDE_DIRECTORIES (${PROJECT_SOURCE_DIR}/Cantera/src/thermo)
INCLUDE_DIRECTORIES (${PROJECT_SOURCE_DIR}/Cantera/src/numerics)

ADD_LIBRARY(equil  ${EQUIL_SRCS} )
TARGET_LINK_LIBRARIES( equil numerics thermo ctbase)

SET(EQUIL_H  ChemEquil.h MultiPhase.h MultiPhaseEquil.h equil.h PropertyCalculator.h)
INSTALL_FILES(/include/cantera/kernel FILES ${EQUIL_H})

