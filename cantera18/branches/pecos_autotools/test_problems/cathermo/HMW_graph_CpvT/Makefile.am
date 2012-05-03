cc_sources = HMW_graph_CpvT.cpp sortAlgorithms.cpp
h_sources  = sortAlgorithms.h  TemperatureTable.h

INC                  = -I$(top_builddir)/build/include/
AM_CPPFLAGS          = $(INC)
AM_CXXFLAGS          = $(AM_CPPFLAGS)

LINK                 = -lctcxx -luser -loneD -lzeroD -lequil -lkinetics -ltransport -lthermo 
LINK                += -lctnumerics -lctmath -ltpx -lctspectra -lconverters -lctbase -lcvode 
LINK                += -lctlapack -lctblas -lctf2c -lm -lstdc++
AM_LDFLAGS           = -L$(top_builddir)/build/lib/ 
LIBS                 =  $(LINK)

bin_PROGRAMS         = HMW_graph_CpvT
library_includedir   = $(INC)
#library_include_HEADERS = $(h_sources)

#-----------------------
# Cantera DH graph test
#-----------------------

HMW_graph_CpvT_SOURCES = $(cc_sources)

TESTS_ENVIRONMENT = 
TESTS             = runtest

EXTRA_DIST = runtest