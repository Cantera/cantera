h_sources  = ArrayViewer.h DenseMatrix.h \
             funcs.h ctlapack.h Func1.h FuncEval.h \
             polyfit.h\
             BandMatrix.h Integrator.h DAE_Solver.h ResidEval.h sort.h \
             SquareMatrix.h ResidJacEval.h NonlinearSolver.h

cc_sources = DenseMatrix.cpp funcs.cpp Func1.cpp ODE_integrators.cpp \
	     BandMatrix.cpp DAE_solvers.cpp funcs.cpp sort.cpp       \
	     SquareMatrix.cpp ResidJacEval.cpp NonlinearSolver.cpp

INCLUDES                = -I../base -I../../../
AM_CXXFLAGS             = $(INCLUDES)

lib_LTLIBRARIES         = libctnumerics.la
library_includedir      = $(includedir)

#-----------------------
# Cantera numerics C/C++ library
#-----------------------

libctnumerics_la_LDFLAGS      = $(all_libraries) -release $(GENERIC_RELEASE)
libctnumerics_la_SOURCES      = $(cc_sources) $(h_sources)