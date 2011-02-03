SET (F2C_BLAS_SRCS  dasum.c  daxpy.c  dcabs1.c  dcopy.c  ddot.c  
dgbmv.c  dgemm.c  dgemv.c  dger.c  dnrm2.c  drot.c  drotg.c  
drotm.c  drotmg.c  dsbmv.c  dscal.c  dsdot.c  dspmv.c  dspr.c  dspr2.c  
dswap.c  dsymm.c  dsymv.c  dsyr.c  dsyr2.c  dsyr2k.c  dsyrk.c  dtbmv.c  
dtbsv.c  dtpmv.c  dtpsv.c  dtrmm.c  dtrmv.c  dtrsm.c  dtrsv.c  
dzasum.c  dznrm2.c  idamax.c  lsame.c  xerbla.c )

INCLUDE_DIRECTORIES ( ${PROJECT_SOURCE_DIR}/ext/f2c_libs )
ADD_LIBRARY(ctblas ${F2C_BLAS_SRCS})


