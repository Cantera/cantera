cmake_minimum_required (VERSION 2.6)

PROJECT (Cantera)

<<<<<<< .mine
ADD_EXECUTABLE (starting start.cpp)
=======
SET (PROJECT_BINARY_DIR  ${PROJECT_SOURCE_DIR}/build )

CONFIGURE_FILE (${PROJECT_SOURCE_DIR}/config/config.h.in
   ${PROJECT_BINARY_DIR}/Headers/config.h )

>>>>>>> .r78
ADD_SUBDIRECTORY (BackEnd)
ADD_SUBDIRECTORY (FrontEnds)
ADD_EXECUTABLE (testit testit.cpp )

<<<<<<< .mine
=======
SET (CMAKE_INSTALL_PREFIX /usr/local)

>>>>>>> .r78
