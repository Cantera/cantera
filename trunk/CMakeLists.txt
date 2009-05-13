cmake_minimum_required (VERSION 2.6)

PROJECT (Cantera)

SET (PROJECT_BINARY_DIR  ${PROJECT_SOURCE_DIR}/build )

CONFIGURE_FILE  (${PROJECT_SOURCE_DIR}/config/config.h.in
   ${PROJECT_BINARY_DIR}/Headers/config.h )

ADD_SUBDIRECTORY (BackEnd)
ADD_SUBDIRECTORY (FrontEnds)


SET (CMAKE_INSTALL_PREFIX /usr/local/)

INSTALL_FILES (/include  ${PROJECT_BINARY_DIR}/Headers/Cantera.h)
