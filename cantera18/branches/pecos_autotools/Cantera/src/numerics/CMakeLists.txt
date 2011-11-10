SET (NUMERICS_SRCS  DenseMatrix.cpp funcs.cpp Func1.cpp 
                 ODE_integrators.cpp  BandMatrix.cpp DAE_solvers.cpp 
                 sort.cpp  )

INCLUDE_DIRECTORIES (${PROJECT_SOURCE_DIR})
INCLUDE_DIRECTORIES (${PROJECT_SOURCE_DIR}/Cantera/src/base)

ADD_LIBRARY(numerics ${NUMERICS_SRCS})

TARGET_LINK_LIBRARIES (numerics ctbase cvode)

SET(NUMERICS_H ArrayViewer.h CVodeInt.h CVodesIntegrator.h DenseMatrix.h 
                 funcs.h ctlapack.h Func1.h FuncEval.h 
                 polyfit.h
                 BandMatrix.h Integrator.h DAE_Solver.h ResidEval.h sort.h)
INSTALL_FILES( /include/cantera/kernel FILES ${NUMERICS_H})
