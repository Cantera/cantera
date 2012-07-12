cc_sources = statmech_test.cpp

INC                  = -I. -I$(top_builddir)/build/include/ -g
AM_CPPFLAGS          = $(INC)
AM_CXXFLAGS          = $(AM_CPPFLAGS)

LINK                 = -lctcxx -luser -loneD -lzeroD -lequil -lkinetics -ltransport -lthermo 
LINK                += -lctnumerics -lctmath -ltpx -lctspectra -lconverters -lctbase -lcvode 
LINK                += -lctlapack -lctblas -lctf2c -lm -lstdc++
AM_LDFLAGS           = -L$(top_builddir)/build/lib/ 
LIBS                 =  $(LINK)

bin_PROGRAMS         = statmech_test
library_includedir   = $(INC)

#-----------------------
# Cantera DH graph test
#-----------------------

statmech_test_SOURCES = $(cc_sources)


TESTS_ENVIRONMENT = 
TESTS             = runtest_stat

EXTRA_DIST = runtest_stat  