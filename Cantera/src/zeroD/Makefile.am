h_sources  = Reactor.h ReactorBase.h FlowDevice.h Wall.h ReactorNet.h \
             flowControllers.h PID_Controller.h Reservoir.h FlowReactor.h \
             ConstPressureReactor.h ReactorFactory.h

cc_sources = Reactor.cpp ReactorBase.cpp FlowDevice.cpp Wall.cpp \
             ReactorNet.cpp FlowReactor.cpp ConstPressureReactor.cpp \
             ReactorFactory.cpp

INCLUDES                = -I../base -I../thermo -I../kinetics -I../numerics -I../../../
AM_CXXFLAGS             = $(INCLUDES)

lib_LTLIBRARIES         = libzeroD.la
library_includedir      = $(includedir)

#-----------------------
# Cantera OneD C/C++ library
#-----------------------

libzeroD_la_LDFLAGS      = $(all_libraries) -release $(GENERIC_RELEASE)
libzeroD_la_SOURCES      = $(cc_sources) $(h_sources)