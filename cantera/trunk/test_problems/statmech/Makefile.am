INC                  = -I. -I$(top_builddir)/include/
AM_CPPFLAGS          = $(INC)
AM_CXXFLAGS          = $(AM_CPPFLAGS)

LINK                 = -loneD -lzeroD -lequil -lkinetics -ltransport -lthermo 
LINK                += -lctnumerics -lctmath -ltpx -lctspectra -lctbase -lcvode 
LINK                += -lctlapack -lctblas -lctf2c -lm  -lctf2c -lstdc++ -lexecstream -pthread
AM_LDFLAGS           = -L$(top_builddir)/build/lib/ 
LIBS                 =  $(LINK)

check_PROGRAMS         = statmech_test \
                         statmech_test_poly \
                         statmech_test_Fe   \
	                 statmech_transport

library_includedir   = $(INC)

#-----------------------
# Cantera DH graph test
#-----------------------

statmech_test_SOURCES      = statmech_test.cpp
statmech_test_poly_SOURCES = statmech_test_poly.cpp
statmech_test_Fe_SOURCES   = statmech_test_Fe.cpp
statmech_transport_SOURCES = statmech_transport.cpp

TESTS_ENVIRONMENT = 	
TESTS             = statmech_test \
                    statmech_test_Fe \
	            statmech_transport
#                    statmech_test_poly

EXTRA_DIST = runtest_stat 