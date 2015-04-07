cc_sources = HMW_graph_VvT.cpp

INC                  = -I$(top_builddir)/include/ -I$(top_builddir)/test_problems/shared
AM_CPPFLAGS          = $(INC)
AM_CXXFLAGS          = $(AM_CPPFLAGS)

LINK                 = -loneD -lzeroD -lequil -lkinetics -ltransport -lthermo 
LINK                += -lctnumerics -lctmath -ltpx -lctspectra -lctbase -lcvode 
LINK                += -lctlapack -lctblas -lctf2c -lm -lstdc++ -lexecstream -pthread
LIBS                 =  $(LINK)
AM_LDFLAGS           = -L$(top_builddir)/build/lib/ 

bin_PROGRAMS         = HMW_graph_VvT
library_includedir   = $(INC)

#-----------------------
# Cantera DH graph test
#-----------------------

HMW_graph_VvT_SOURCES = $(cc_sources)

TESTS_ENVIRONMENT = 
TESTS             = runtest
