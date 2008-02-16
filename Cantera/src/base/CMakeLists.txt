SET (CTBASE_SRCS misc.cpp ct2ctml.cpp ctml.cpp 
                             plots.cpp stringUtils.cpp xml.cpp clockWC.cpp)
ADD_LIBRARY(ctbase ${CTBASE_SRCS})

SET (CTBASE_H global.h ct2ctml.h ctml.h 
                        ct_defs.h ctexceptions.h logger.h XML_Writer.h 
                        ctml.h plots.h stringUtils.h xml.h config.h utilities.h 
                       Array.h vec_functions.h global.h FactoryBase.h clockWC.h )
INSTALL_FILES(/include/cantera/kernel FILES ${CTBASE_H})
