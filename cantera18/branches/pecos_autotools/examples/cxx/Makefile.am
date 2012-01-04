h_sources  = example_utils.h

cc_sources = equil_example1.cpp examples.cpp flame1.cpp kinetics_example1.cpp \
	     kinetics_example2.cpp kinetics_example3.cpp rxnpath_example1.cpp \
	     transport_example1.cpp transport_example2.cpp

INCLUDES                = -I../base -I../../build/include/
AM_CXXFLAGS             = $(INCLUDES)

bin_PROGRAMS            = cxx_examples
library_includedir      = $(includedir)

#-----------------------
# Cantera cxx examples
#-----------------------

cxx_examples_LDFLAGS      = $(all_libraries) -release $(GENERIC_RELEASE)
cxx_examples_SOURCES      = $(cc_sources) $(h_sources)
