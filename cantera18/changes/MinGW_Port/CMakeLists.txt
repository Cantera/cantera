PROJECT (Cantera)

#----------------------------
#     user-configurable settings
#----------------------------
INCLUDE (config_cantera.cmake)


#----------------------------
#     Python
#----------------------------     
#if (PYTHON_CMD STREQUAL "default")
INCLUDE( FindPythonInterp )
INCLUDE( FindPythonLibs)
SET( PYTHON_EXE ${PYTHON_EXECUTABLE} )
#else (PYTHON_CMD STREQUAL "default")
#  SET( PYTHON_EXE ${PYTHON_CMD} )
#endif (PYTHON_CMD STREQUAL "default")


#---------------------------------------
#      configuration
#---------------------------------------
CONFIGURE_FILE (
    ${PROJECT_SOURCE_DIR}/config.h_cmake.in
    ${PROJECT_BINARY_DIR}/config.h )


#----------------------------------------
#      output paths
#----------------------------------------
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

INSTALL_FILES(/include/cantera/kernel FILES config.h)
INSTALL_FILES(/include/cantera FILES config.h)
ADD_SUBDIRECTORY (ext)
ADD_SUBDIRECTORY (Cantera)

SET(CMAKE_INSTALL_PREFIX /Applications/Cantera )
INSTALL_FILES (  /include/cantera  .h  ${CANTERA_CXX_HEADERS} )
