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

# configuration
CONFIGURE_FILE (
    ${PROJECT_SOURCE_DIR}/config.h_cmake.in
    ${PROJECT_BINARY_DIR}/config.h )

INCLUDE_DIRECTORIES (${PROJECT_BINARY_DIR})
INCLUDE_DIRECTORIES (${PROJECT_BINARY_DIR}/build/include/cantera)

SET (LIBRARY_OUTPUT_PATH 
  ${PROJECT_BINARY_DIR}/build/lib/${CMAKE_SYSTEM_PROCESSOR}-${CMAKE_SYSTEM_VERSION} CACHE PATH "Single directory for all libraries"
  )

SET (EXECUTABLE_OUTPUT_PATH 
  ${PROJECT_BINARY_DIR}/build/bin/${CMAKE_SYSTEM_PROCESSOR}-${CMAKE_SYSTEM_VERSION} CACHE PATH "Single directory for all executables"
  )

MARK_AS_ADVANCED (
    LIBRARY_OUTPUT_PATH
    EXECUTABLE_OUTPUT_PATH
   )

ADD_SUBDIRECTORY (ext)
ADD_SUBDIRECTORY (Cantera)

SET ( CANTERA_INCLUDE_DIR ${PROJECT_SOURCE_DIR}/Cantera/cxx/include )

SET(CMAKE_INSTALL_PREFIX /Applications/Cantera/include )
INSTALL_FILES (  /cantera  .h  ${CANTERA_CXX_HEADERS} )
