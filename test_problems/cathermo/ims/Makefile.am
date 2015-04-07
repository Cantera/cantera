cc_sources = IMSTester.cpp

INC                  = -I$(top_builddir)/include/
AM_CPPFLAGS          = $(INC)
AM_CXXFLAGS          = $(AM_CPPFLAGS)

LINK                 = -loneD -lzeroD -lequil -lkinetics -ltransport -lthermo 
LINK                += -lctnumerics -lctmath -ltpx -lctspectra -lctbase -lcvode 
LINK                += -lctlapack -lctblas -lctf2c -lm -lstdc++ -lexecstream -pthread
AM_LDFLAGS           = -L$(top_builddir)/build/lib/
LIBS                 =  $(LINK)

bin_PROGRAMS         = IMSTester
library_includedir   = $(INC)

#-----------------------
# Cantera DH graph test
#-----------------------

IMSTester_SOURCES = $(cc_sources)

TESTS_ENVIRONMENT = 
TESTS             = runtest
