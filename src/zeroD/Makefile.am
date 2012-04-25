h_sources  = 
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
