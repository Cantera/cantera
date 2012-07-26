cc_sources = runDiamond.cpp

INC                  = -I. -I$(top_builddir)/include/
AM_CPPFLAGS          = $(INC)
AM_CXXFLAGS          = $(AM_CPPFLAGS)

LINK                 = -loneD -lzeroD -lequil -lkinetics -ltransport -lthermo 
LINK                += -lctnumerics -lctmath -ltpx -lctspectra -lctbase -lcvode 
LINK                += -lctlapack -lctblas -lctf2c -lm  -lctf2c -lstdc++ -lexecstream -pthread
AM_LDFLAGS           = -L$(top_builddir)/build/lib/ 
LIBS                 =  $(LINK)

bin_PROGRAMS         = runDiamond
library_includedir   = $(INC)

#-----------------------
# Cantera DH graph test
#-----------------------

runDiamond_SOURCES = $(cc_sources)

#TESTS_ENVIRONMENT = 
#TESTS             = runtest
