SET (SPECTRA_SRCS  rotor.cpp LineBroadener.cpp spectralUtilities.cpp)

INCLUDE_DIRECTORIES (${PROJECT_SOURCE_DIR})
INCLUDE_DIRECTORIES (${PROJECT_SOURCE_DIR}/Cantera/src/base)
INCLUDE_DIRECTORIES (${PROJECT_SOURCE_DIR}/Cantera/src/thermo)
INCLUDE_DIRECTORIES (${PROJECT_SOURCE_DIR}/Cantera/src/numerics)
INCLUDE_DIRECTORIES (${PROJECT_SOURCE_DIR}/Cantera/src/transport)

ADD_LIBRARY(spectra  ${SPECTRA_SRCS} )

SET (SPECTRA_H  rotor.h LineBroadener.h Nuclei.h spectralUtilities.h)
INSTALL_FILES(/include/cantera/kernel FILES ${SPECTRA_H})

