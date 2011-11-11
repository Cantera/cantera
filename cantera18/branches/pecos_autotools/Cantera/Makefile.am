# $Id: $

# will need to add python and matlab 
SUBDIRS                 = src

h_sources               = cantera.h
cc_sources              = cantera.cpp

INCLUDES                = 
AM_FCFLAGS              = $(INCLUDES)
AM_CFLAGS               = $(INCLUDES)
AM_CPPFLAGS             = $(INCLUDES) 

lib_LTLIBRARIES         = libcantera.la
library_includedir      = $(includedir)
library_include_HEADERS = cantera.h
libcantera_la_LDFLAGS   = $(all_libraries) -release $(GENERIC_RELEASE)
libcantera_la_SOURCES   = $(cc_sources) $(h_sources)


