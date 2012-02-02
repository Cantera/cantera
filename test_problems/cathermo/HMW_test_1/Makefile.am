cc_sources = HMW_test_1.cpp

INC                  = -I. -I$(top_builddir)/build/include/ -I$(top_builddir)/build/include/cantera -I$(top_builddir)/build/include/cantera/kernel
AM_CPPFLAGS          = $(INC)
AM_CXXFLAGS          = $(AM_CPPFLAGS)

LINK                 = -lctcxx -luser -loneD -lzeroD -lequil -lkinetics -ltransport -lthermo 
LINK                += -lctnumerics -lctmath -ltpx -lctspectra -lconverters -lctbase -lcvode 
LINK                += -lctlapack -lctblas -lctf2c -lm -lstdc++
AM_LDFLAGS           = -L/workspace/src/cantera/pecos_autotools/build/lib/ 
LIBS                 =  $(LINK)

bin_PROGRAMS         = HMW_test_1
library_includedir   = $(INC)

#-----------------------
# Cantera DH graph test
#-----------------------

HMW_test_1_SOURCES = $(cc_sources)

TESTS_ENVIRONMENT = 
TESTS             = runtest
