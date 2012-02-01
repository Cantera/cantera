h_sources  = example_utils.h

cc_sources = examples.cpp kinetics_example1.cpp kinetics_example2.cpp \
	     kinetics_example3.cpp  equil_example1.cpp                \
	     transport_example1.cpp transport_example2.cpp            \
	     rxnpath_example1.cpp 

LINK                 = -lkinetics -ltransport -lthermo -lctnumerics
LINK                += -lctbase -ltpx -lctmath -lconverters 
LINK                += -luser -lcvode -lstdc++ -lctlapack 
LINK                += -lctblas -lctf2c -lctcxx

#LIBS                 = -L$(top_builddir)/build/lib/ $(LINK)
LIBS                 = -L/workspace/src/cantera/pecos_autotools/build/lib/ $(LINK)
AM_CPPFLAGS          = -I. -I/workspace/src/cantera/pecos_autotools/build/include/ -I/workspace/src/cantera/pecos_autotools/build/include/cantera -I/workspace/src/cantera/pecos_autotools/build/include/cantera/kernel
AM_CXXFLAGS          = $(AM_CPPFLAGS)

bin_PROGRAMS         = cxx_examples
library_includedir   = $(includedir) -I/workspace/src/cantera/pecos_autotools/build/include/

#-----------------------
# Cantera cxx examples
#-----------------------

cxx_examples_LDFLAGS = $(all_libraries) -release $(GENERIC_RELEASE)
cxx_examples_SOURCES = $(cc_sources) $(h_sources)
