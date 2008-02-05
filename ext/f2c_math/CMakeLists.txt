SET ( CTMATH_SRCS  mach.cpp  ddaspk.c  dgbefa.c  dgbsl.c  dgefa.c  dgesl.c 
         dp1vlu.c dpcoef.c  dpolft.c fdump.c j4save.c pcoef.c polfit.c pvalue.c 
         xercnt.c xerhlt.c xermsg.c xerprn.c xersve.c  xgetua.c printstring.c)

INCLUDE_DIRECTORIES ( ${PROJECT_BINARY_DIR} ) 
INCLUDE_DIRECTORIES ( ${PROJECT_BINARY_DIR}/ext/f2c_libs ) 
INCLUDE_DIRECTORIES ( ${PROJECT_BINARY_DIR}/build/include/cantera ) 

ADD_LIBRARY( ctmath ${CTMATH_SRCS} )

