h_sources  = $(top_builddir)/include/cantera/zeroD/Reactor.h \
             $(top_builddir)/include/cantera/zeroD/ReactorBase.h \
             $(top_builddir)/include/cantera/zeroD/FlowDevice.h \
             $(top_builddir)/include/cantera/zeroD/Wall.h \
             $(top_builddir)/include/cantera/zeroD/ReactorNet.h \
             $(top_builddir)/include/cantera/zeroD/flowControllers.h \
             $(top_builddir)/include/cantera/zeroD/Reservoir.h \
             $(top_builddir)/include/cantera/zeroD/FlowReactor.h \
             $(top_builddir)/include/cantera/zeroD/ConstPressureReactor.h \
             $(top_builddir)/include/cantera/zeroD/ReactorFactory.h

cc_sources = Reactor.cpp ReactorBase.cpp FlowDevice.cpp Wall.cpp \
             ReactorNet.cpp FlowReactor.cpp ConstPressureReactor.cpp \
             ReactorFactory.cpp

AM_CPPFLAGS             = -I$(top_builddir)/include
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
