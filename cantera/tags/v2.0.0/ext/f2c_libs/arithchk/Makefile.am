AM_CPPFLAGS             = -I. -I../
AM_CXXFLAGS             = $(AM_CPPFLAGS) -Wno-deprecated-declarations -pthread
AM_CFLAGS               = $(AM_CPPFLAGS) -Wno-deprecated-declarations -pthread -finline-functions -Wno-inline -g -Wno-format-security -DNO_FPINIT -DNDEBUG
AM_FCFLAGS              = $(AM_CPPFLAGS)
AM_LDFLAGS              = -L$(top_builddir)/build/lib/ -lm

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
