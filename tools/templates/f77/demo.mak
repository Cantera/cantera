#!/bin/sh

# This Makefile builds a Fortran 77 application that uses Cantera.  By
# default, the main program file is 'demo.f,' which prints out some
# properties of a reacting gas mixture. It uses the library
# 'demo_ftnlib.cpp,' which contains Fortran-callable functions that
# are implemented with C++ calls to Cantera.  

# To build program 'demo', simply type 'make', or 'make -f <this
# file>' if this file is named something other than 'Makefile.'  

# Once you have verified that the demo runs, edit this file to replace
# object file 'demo.o' with your own object file or files. You may
# continue to use 'demo_ftnlib' if it serves your needs, or else
# replace it with a different interface library.


#------------------------  edit this block ---------------------------------

# the name of the executable program to be created
PROG_NAME = demo

# the object files to be linked together. 
OBJS = demo.o demo_ftnlib.o

# additional flags to be passed to the linker. If your program
# requires other external libraries, put them here
LINK_OPTIONS = 

#---------------------------------------------------------------------------
# You probably don't need to edit anything below.


# the Fortran compiler
FORT = g77

# Fortran compile flags  
FORT_FLAGS = -O2 -fno-second-underscore 

# Fortran libraries
FORT_LIBS = -lg2c -lgcc

# the C++ compiler
CXX = g++

# C++ compile flags
CXX_FLAGS = -O0 -Wall

# external libraries
EXT_LIBS =  -loneD -lzeroD -ltransport -lconverters -lcantera -lrecipes -lcvode -lctlapack -lctmath -lctblas 

# the directory where the Cantera libraries are located
CANTERA_LIBDIR=/usr/local/lib/cantera

# the directory where Cantera include files may be found.
CANTERA_INCDIR=/usr/local/include/cantera

# flags passed to the C++ compiler/linker for the linking step
LCXX_FLAGS = -L$(CANTERA_LIBDIR) -O0 -Wall

# how to compile C++ source files to object files
.cpp.o:
	$(CXX) -c $< -I$(CANTERA_INCDIR) $(CXX_FLAGS)

# how to compile Fortran source files to object files
.f.o: 
	$(FORT) -c $< $(FORT_FLAGS)

PROGRAM = $(PROG_NAME)$(EXE_EXT)

DEPENDS = $(OBJS:.o=.d)

all: $(PROGRAM)

$(PROGRAM): $(OBJS)
	$(CXX) -o $(PROGRAM) $(OBJS) $(LCXX_FLAGS) $(CANTERA_LIBS) $(LINK_OPTIONS) $(EXT_LIBS)  $(FORT_LIBS)

%.d:
	g++ -MM $*.cpp > $*.d

clean:
	$(RM) $(OBJS) $(PROGRAM)

depends: $(DEPENDS)
	cat *.d > .depends
	$(RM) $(DEPENDS) 

TAGS: 
	etags *.h *.cpp

ifeq ($(wildcard .depends), .depends)
include .depends
endif






