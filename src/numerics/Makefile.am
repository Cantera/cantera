h_sources  = $(top_builddir)/include/cantera/numerics/DenseMatrix.h \
             $(top_builddir)/include/cantera/numerics/funcs.h \
             $(top_builddir)/include/cantera/numerics/ctlapack.h \
             $(top_builddir)/include/cantera/numerics/Func1.h \
             $(top_builddir)/include/cantera/numerics/FuncEval.h \
             $(top_builddir)/include/cantera/numerics/polyfit.h \
             $(top_builddir)/include/cantera/numerics/BandMatrix.h \
             $(top_builddir)/include/cantera/numerics/Integrator.h \
             $(top_builddir)/include/cantera/numerics/DAE_Solver.h \
             $(top_builddir)/include/cantera/numerics/ResidEval.h \
             $(top_builddir)/include/cantera/numerics/SquareMatrix.h \
             $(top_builddir)/include/cantera/numerics/ResidJacEval.h \
             $(top_builddir)/include/cantera/numerics/NonlinearSolver.h 

# $(top_builddir)/include/cantera/numerics/ArrayViewer.h
#             $(top_builddir)/include/cantera/numerics/sort.h

#            $(top_builddir)/include/cantera/numerics/CVodeInt.h

cc_sources = BEulerInt.cpp DenseMatrix.cpp funcs.cpp Func1.cpp ODE_integrators.cpp \
	     BandMatrix.cpp DAE_solvers.cpp CVodeInt.cpp GeneralMatrix.cpp   \
	     SquareMatrix.cpp ResidJacEval.cpp NonlinearSolver.cpp IDA_Solver.cpp \
	     RootFind.cpp solveProb.cpp 

AM_CPPFLAGS             = -I$(top_builddir)/include/ -I$(top_builddir)
AM_CXXFLAGS             = $(AM_CPPFLAGS)

lib_LTLIBRARIES         = $(top_builddir)/build/lib/libctnumerics.la
library_includedir      = $(includedir)

#-----------------------
# Cantera Numerics C/C++ library
#-----------------------

__top_builddir__build_lib_libctnumerics_la_LDFLAGS      = $(all_libraries) -release $(GENERIC_RELEASE)
__top_builddir__build_lib_libctnumerics_la_SOURCES      = $(cc_sources) $(h_sources)

__top_builddir__build_lib_libctnumerics_include_HEADERS = $(h_sources)
__top_builddir__build_lib_libctnumerics_includedir = $(prefix)/include/cantera/numerics
