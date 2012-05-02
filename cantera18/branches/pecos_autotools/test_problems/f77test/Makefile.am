fc_sources = isentropic.f 
#ctlib.f
cc_sources = demo_ftnlib.cpp

INC                  = -I. -I$(top_builddir)/build/include/
AM_FCFLAGS           = -fno-second-underscore
AM_CPPFLAGS          = $(INC)
AM_CXXFLAGS          = $(AM_CPPFLAGS)

LINK                 = -luser -loneD -lzeroD -lequil -lkinetics -ltransport -lthermo 
LINK                += -lctnumerics -lctmath -ltpx -lctspectra -lconverters -lctbase -lcvode 
LINK                += -lctlapack -lctblas -lctf2c -lctcxx -lm  -lctf2c -lstdc++ -lgfortran
AM_LDFLAGS           = -L$(top_builddir)/build/lib/ 
LIBS                 =  $(LINK)

bin_PROGRAMS         = isentropic ctlib
library_includedir   = $(INC)

#-----------------------
# Cantera f77 test
#-----------------------

isentropic_SOURCES = $(cc_sources) $(fc_sources)
ctlib_SOURCES      = $(cc_sources) ctlib.f



TESTS_ENVIRONMENT = 
TESTS             = runtest

EXTRA_DIST = runtest