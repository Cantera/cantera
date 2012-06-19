AM_CPPFLAGS             = -I. -I$(top_builddir)/ext/f2c_libs/
AM_CXXFLAGS             = $(AM_CPPFLAGS)
AM_FCFLAGS              = $(AM_CPPFLAGS)

lib_LTLIBRARIES         = 
library_includedir      = 
library_include_HEADERS = $(h_sources)

bin_PROGRAMS         = arithchk

#-----------------------
# Cantera Arithchk 
#-----------------------

arithchk_SOURCES = arithchk.c

# header file creation
all:
	@echo creating arithchk.h
	./arithchk > ../arith.h

CLEANFILES = *.o
