h_sources  = example_utils.h

cc_sources = examples.cpp kinetics_example1.cpp kinetics_example2.cpp \
	     kinetics_example3.cpp  equil_example1.cpp                \
	     transport_example1.cpp transport_example2.cpp            \
	     rxnpath_example1.cpp 

LINK                 = -lkinetics -ltransport -lthermo -lctnumerics
LINK                += -lctbase -ltpx -lctmath -lconverters 
LINK                += -luser -lcvode -lstdc++ -lctlapack 
LINK                += -lctblas
# -lctf2c -lctcxx -luser -lstdc++

LIBS                 = -L$(top_builddir)/build/lib/ $(LINK)

AM_CPPFLAGS          = -I../base -I$(top_builddir)/build/include/ -I$(top_builddir)/build/include/cantera
AM_CXXFLAGS          = $(AM_CPPFLAGS)

bin_PROGRAMS         = cxx_examples
library_includedir   = $(includedir)

#-----------------------
# Cantera cxx examples
#-----------------------

cxx_examples_LDFLAGS = $(all_libraries) -release $(GENERIC_RELEASE)
cxx_examples_SOURCES = $(cc_sources) $(h_sources)
