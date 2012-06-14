INC                  = -I. -I$(top_builddir)/include/cantera -I$(top_builddir)/test_problems/shared
AM_CPPFLAGS          = $(INC)
AM_CXXFLAGS          = $(AM_CPPFLAGS)

#LINK                 = -luser -loneD -lzeroD -lequil -lkinetics -ltransport -lthermo 
LINK                 = -loneD -lzeroD -lequil -lkinetics -ltransport -lthermo 
LINK                += -lctnumerics -lctmath -ltpx -lctspectra -lconverters -lctbase -lcvode 
LINK                += -lctlapack -lctblas -lctf2c -lm -lstdc++
AM_LDFLAGS           = -L$(top_builddir)/build/lib/ 
LIBS                 =  $(LINK)

bin_PROGRAMS         = DH_graph_1
library_includedir   = $(INC)

#-----------------------
# Cantera DH graph test
#-----------------------

DH_graph_1_SOURCES = DH_graph_1.cpp

TESTS_ENVIRONMENT = 
TESTS             = runtest
