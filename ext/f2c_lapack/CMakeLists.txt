
SET (F2C_LAPACK_SRCS  dgbtrf.c  dgbtrs.c  dgbtf2.c dgbsv.c  
   dgebd2.c   dgebrd.c   dgelq2.c   dgelqf.c   dgelss.c   dgeqr2.c   dgeqrf.c   dgetf2.c  
   dgetrf.c   dgetri.c   dgetrs.c   dlabad.c   dlabrd.c   dlacpy.c   dlamch.c  dlange.c  
   dlapy2.c  dlarf.c  dlarfb.c  dlarfg.c  dlarft.c  dlartg.c  dlas2.c  dlascl.c  dlaset.c  
   dlasq1.c  dlasq2.c  dlasq3.c  dlasq4.c  dlasq5.c  dlasq6.c  dlasr.c  dlasrt.c  
   dlassq.c   dlasv2.c  dlaswp.c  dorg2r.c  dorgbr.c  dorgl2.c  dorglq.c  dorgqr.c  
   dorm2r.c  dormbr.c  dorml2.c  dormlq.c  dormqr.c  drscl.c  dtrtri.c  
   dtrti2.c  ieeeck.c  ilaenv.c  )

ADD_LIBRARY(ctlapack ${F2C_LAPACK_SRCS})



