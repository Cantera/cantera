PROJECT (Cantera)

#----------------------------
#     user-configurable settings
#----------------------------
INCLUDE (config.cmake)


#operating system
SET(os ${CMAKE_SYSTEM_NAME})
if (os STREQUAL "Darwin")
    SET(DARWIN 1)
endif (os STREQUAL "Darwin")

if (PYTHON_CMD STREQUAL "default")
   INCLUDE( FindPythonInterp )
   SET( PYTHON_EXE ${PYTHON_EXECUTABLE} )
else (PYTHON_CMD STREQUAL "default")
  SET( PYTHON_EXE ${PYTHON_CMD} )
endif (PYTHON_CMD STREQUAL "default")

INCLUDE(cmake/thermo.cmake)

INCLUDE(cmake/fortran.cmake)

# configuration
CONFIGURE_FILE (
    ${PROJECT_SOURCE_DIR}/config.h_cmake.in
    ${PROJECT_BINARY_DIR}/config.h )

include_directories(${PROJECT_BINARY_DIR})

set(LIBRARY_OUTPUT_PATH ${PROJECT_BINARY_DIR}/build/lib/${CMAKE_SYSTEM_PROCESSOR}-${CMAKE_SYSTEM_VERSION})

add_subdirectory(Cantera)


