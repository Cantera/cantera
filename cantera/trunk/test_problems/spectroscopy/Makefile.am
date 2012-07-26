cc_sources = spectratest.cpp

INC                  = -I. -I$(top_builddir)/include/ -I$(top_builddir)/src/
AM_CPPFLAGS          = $(INC)
AM_CXXFLAGS          = $(AM_CPPFLAGS)

LINK                 = -loneD -lzeroD -lequil -lkinetics -ltransport -lthermo 
LINK                += -lctnumerics -lctmath -ltpx -lctspectra -lctbase -lcvode 
LINK                += -lctlapack -lctblas -lctf2c -lm  -lctf2c -lstdc++ -lexecstream -pthread
AM_LDFLAGS           = -L$(top_builddir)/build/lib/ 
LIBS                 =  $(LINK)

bin_PROGRAMS         = spectroscopy
library_includedir   = $(INC)

#-----------------------
# Cantera DH graph test
#-----------------------

spectroscopy_SOURCES = $(cc_sources)

TESTS_ENVIRONMENT = 
TESTS             = spectroscopy
