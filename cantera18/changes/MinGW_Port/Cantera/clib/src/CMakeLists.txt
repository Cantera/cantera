SET (CLIB_SRCS ct.cpp Storage.cpp ctsurf.cpp ctrpath.cpp ctreactor.cpp
                        ctfunc.cpp ctxml.cpp ctonedim.cpp ctmultiphase.cpp )

INCLUDE_DIRECTORIES (${PROJECT_SOURCE_DIR})
INCLUDE_DIRECTORIES (${PROJECT_SOURCE_DIR}/Cantera/src/base)
INCLUDE_DIRECTORIES (${PROJECT_SOURCE_DIR}/Cantera/src/thermo)
INCLUDE_DIRECTORIES (${PROJECT_SOURCE_DIR}/Cantera/src/equil)
INCLUDE_DIRECTORIES (${PROJECT_SOURCE_DIR}/Cantera/src/numerics)
INCLUDE_DIRECTORIES (${PROJECT_SOURCE_DIR}/Cantera/src/kinetics)
INCLUDE_DIRECTORIES (${PROJECT_SOURCE_DIR}/Cantera/src/transport)
INCLUDE_DIRECTORIES (${PROJECT_SOURCE_DIR}/Cantera/src/converters)
INCLUDE_DIRECTORIES (${PROJECT_SOURCE_DIR}/Cantera/src/zeroD)
INCLUDE_DIRECTORIES (${PROJECT_SOURCE_DIR}/Cantera/src/oneD)

# INCLUDE_DIRECTORIES (${PROJECT_BINARY_DIR}/build/include/cantera)

ADD_LIBRARY ( clib ${CLIB_SRCS} )

TARGET_LINK_LIBRARIES (clib oneD zeroD equil kinetics thermo ctbase)

SET (CLIB_H  ct.h Storage.h ctsurf.h ctrpath.h ctreactor.h
                        ctfunc.h ctxml.h ctonedim.h ctmultiphase.h )



