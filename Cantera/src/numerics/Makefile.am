h_sources  = ArrayViewer.h DenseMatrix.h funcs.h ctlapack.h Func1.h \
             FuncEval.h polyfit.h \
             BandMatrix.h Integrator.h DAE_Solver.h ResidEval.h sort.h \
             SquareMatrix.h ResidJacEval.h NonlinearSolver.h CVodeInt.h

cc_sources = DenseMatrix.cpp funcs.cpp Func1.cpp ODE_integrators.cpp \
	     BandMatrix.cpp DAE_solvers.cpp sort.cpp                 \
	     SquareMatrix.cpp ResidJacEval.cpp NonlinearSolver.cpp   \
	     CVodeInt.cpp

AM_CPPFLAGS             = -I$(top_builddir)/build/include/ -I$(top_builddir)/Cantera/ext/cvode
AM_CXXFLAGS             = $(AM_CPPFLAGS)

lib_LTLIBRARIES         = $(top_builddir)/build/lib/libctnumerics.la
library_includedir      = $(includedir)

#-----------------------
# Cantera numerics C/C++ library
#-----------------------

__top_builddir__build_lib_libctnumerics_la_LDFLAGS      = $(all_libraries) -release $(GENERIC_RELEASE)
__top_builddir__build_lib_libctnumerics_la_SOURCES      = $(cc_sources) $(h_sources)


# header file accumulation
all:
	@echo copying headers $<
	cp *.h $(top_builddir)/build/include/
