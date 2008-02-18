
INCLUDE_DIRECTORIES (${PROJECT_SOURCE_DIR})

SET (CTBASE_SRCS misc.cpp ct2ctml.cpp ctml.cpp 
                             plots.cpp stringUtils.cpp xml.cpp clockWC.cpp)
ADD_LIBRARY(ctbase ${CTBASE_SRCS})

SET (CTBASE_H global.h ctml.h 
                       ct_defs.h ctexceptions.h logger.h XML_Writer.h 
                       ctml.h plots.h stringUtils.h xml.h utilities.h 
                       Array.h vec_functions.h global.h FactoryBase.h clockWC.h )
INSTALL_FILES(/include/cantera/kernel FILES ${CTBASE_H})
