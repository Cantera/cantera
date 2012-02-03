SET ( TPX_SRCS  HFC134a.cpp Heptane.cpp Hydrogen.cpp Methane.cpp 
        Nitrogen.cpp  Oxygen.cpp  RedlichKwong.cpp   CarbonDioxide.cpp 
        Sub.cpp Water.cpp utils.cpp lk.cpp )

INCLUDE_DIRECTORIES ( ${PROJECT_BINARY_DIR} ) 
#IINCLUDE_DIRECTORIES ( ${PROJECT_BINARY_DIR}/build/include/cantera ) 

ADD_LIBRARY( tpx ${TPX_SRCS}  )

