SET (CLIB_SRCS ct.cpp Storage.cpp ctsurf.cpp ctrpath.cpp ctreactor.cpp
                        ctfunc.cpp ctxml.cpp ctonedim.cpp ctmultiphase.cpp )

INCLUDE_DIRECTORIES (${PROJECT_SOURCE_DIR}/Cantera/src/base)
INCLUDE_DIRECTORIES (${PROJECT_SOURCE_DIR}/Cantera/src/thermo)
INCLUDE_DIRECTORIES (${PROJECT_SOURCE_DIR}/Cantera/src/equil)
INCLUDE_DIRECTORIES (${PROJECT_SOURCE_DIR}/Cantera/src/numerics)
INCLUDE_DIRECTORIES (${PROJECT_SOURCE_DIR}/Cantera/src/kinetics)
INCLUDE_DIRECTORIES (${PROJECT_SOURCE_DIR}/Cantera/src/transport)
INCLUDE_DIRECTORIES (${PROJECT_SOURCE_DIR}/Cantera/src/converters)
INCLUDE_DIRECTORIES (${PROJECT_SOURCE_DIR}/Cantera/src/zeroD)
INCLUDE_DIRECTORIES (${PROJECT_SOURCE_DIR}/Cantera/src/oneD)

SET (LIBRARY_OUTPUT_PATH 
  ${PROJECT_BINARY_DIR}/build/lib/${CMAKE_SYSTEM_PROCESSOR}-${CMAKE_SYSTEM_VERSION} CACHE PATH "Single directory for all libraries"
  )

ADD_LIBRARY(clib ${CLIB_SRCS})

SET (CLIB_H  ct.h Storage.h ctsurf.h ctrpath.h ctreactor.h
                        ctfunc.h ctxml.h ctonedim.h ctmultiphase.h )

TARGET_LINK_LIBRARIES (clib oneD zeroD equil kinetics thermo ctbase)

#             @phase_header_files@)
