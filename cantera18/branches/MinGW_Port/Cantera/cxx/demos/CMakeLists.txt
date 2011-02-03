SET (COMBUST_SRCS combustor.cpp)

INCLUDE_DIRECTORIES (${PROJECT_BINARY_DIR}/build/include)
INCLUDE_DIRECTORIES (${PROJECT_BINARY_DIR})

ADD_EXECUTABLE(combustor.x ${COMBUST_SRCS})
TARGET_LINK_LIBRARIES ( combustor.x zeroD kinetics thermo ctbase tpx cvode )

ADD_EXECUTABLE(flamespeed.x flamespeed.cpp )
TARGET_LINK_LIBRARIES ( flamespeed.x oneD equil transport kinetics numerics 
                                         thermo ctbase tpx ctmath ctlapack ctblas ctf2c )


