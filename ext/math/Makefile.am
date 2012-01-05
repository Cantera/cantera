h_sources               = 
cc_sources              = 

INCLUDES                = 
AM_CXXFLAGS             = $(INCLUDES)

lib_LTLIBRARIES         = $(top_builddir)/build/lib/libctbase.la
library_includedir      = $(top_builddir)/build/include
library_include_HEADERS = $(h_sources)

otherincludedir         = $(top_builddir)/build/include/cantera/kernel
otherinclude_HEADERS    = $(h_sources)

#-----------------------
# Cantera Converters C/C++ library
#-----------------------

__top_builddir__build_lib_libctbase_la_LDFLAGS      = $(all_libraries) -release $(GENERIC_RELEASE)
__top_builddir__build_lib_libctbase_la_SOURCES      = $(cc_sources) $(h_sources)

CLEANFILES = *.o
