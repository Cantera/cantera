cmake_minimum_required (VERSION 2.6)

CONFIGURE_FILE (${PROJECT_SOURCE_DIR}/config/config.h.in  ${PROJECT_BINARY_DIR}/Headers/config.h)

ADD_SUBDIRECTORY (BackEnd)
ADD_SUBDIRECTORY (FrontEnds)

				
