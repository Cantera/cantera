h_sources  = Reactor.h ReactorBase.h FlowDevice.h Wall.h ReactorNet.h \
             flowControllers.h PID_Controller.h Reservoir.h FlowReactor.h \
             ConstPressureReactor.h ReactorFactory.h

cc_sources = Reactor.cpp ReactorBase.cpp FlowDevice.cpp Wall.cpp \
             ReactorNet.cpp FlowReactor.cpp ConstPressureReactor.cpp \
             ReactorFactory.cpp

AM_CPPFLAGS             = -I../base -I../thermo -I../kinetics -I../numerics
AM_CXXFLAGS             = $(AM_CPPFLAGS)

lib_LTLIBRARIES         = $(top_builddir)/build/lib/libzeroD.la
library_includedir      = $(includedir)

#-----------------------
# Cantera OneD C/C++ library
#-----------------------

__top_builddir__build_lib_libzeroD_la_LDFLAGS      = $(all_libraries) -release $(GENERIC_RELEASE)
__top_builddir__build_lib_libzeroD_la_SOURCES      = $(cc_sources) $(h_sources)

__top_builddir__build_lib_libzeroD_include_HEADERS = $(h_sources)
__top_builddir__build_lib_libzeroD_includedir = $(prefix)/include

# header file accumulation
all:
	@echo copying headers $<
	cp -p $(h_sources) $(top_builddir)/build/include/
