cc_sources = PecosTransport.cpp

INC                  = -I. -I$(top_builddir)/include/ -g
AM_CPPFLAGS          = $(INC)
AM_CXXFLAGS          = $(AM_CPPFLAGS)

LINK                 = -loneD -lzeroD -lequil -lkinetics -ltransport -lthermo 
LINK                += -lctnumerics -lctmath -ltpx -lctspectra -lctbase -lcvode 
LINK                += -lctlapack -lctblas -lctf2c -lm -lstdc++ -lexecstream -pthread
AM_LDFLAGS           = -L$(top_builddir)/build/lib/ 
LIBS                 =  $(LINK)

bin_PROGRAMS         = PecosTransport
library_includedir   = $(INC)

#-----------------------
# Cantera Pecos Transport test
#-----------------------

PecosTransport_SOURCES = $(cc_sources)


TESTS_ENVIRONMENT = 
TESTS             = runtest

EXTRA_DIST = runtest