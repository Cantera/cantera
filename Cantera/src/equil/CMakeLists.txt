SET (EQUIL_SRCS BasisOptimize.cpp ChemEquil.cpp MultiPhase.cpp MultiPhaseEquil.cpp
       equilibrate.cpp )

INCLUDE_DIRECTORIES (${PROJECT_SOURCE_DIR}/Cantera/src/base)
INCLUDE_DIRECTORIES (${PROJECT_SOURCE_DIR}/Cantera/src/thermo)
INCLUDE_DIRECTORIES (${PROJECT_SOURCE_DIR}/Cantera/src/numerics)

SET (LIBRARY_OUTPUT_PATH 
  ${PROJECT_BINARY_DIR}/build/lib/${CMAKE_SYSTEM_PROCESSOR}-${CMAKE_SYSTEM_VERSION} CACHE PATH "Single directory for all libraries"
  )

ADD_LIBRARY(equil  ${EQUIL_SRCS} )
TARGET_LINK_LIBRARIES( equil thermo ctbase)

SET(EQUIL_H  ChemEquil.h MultiPhase.h MultiPhaseEquil.h equil.h PropertyCalculator.h)


