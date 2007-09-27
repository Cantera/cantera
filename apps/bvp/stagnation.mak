#!/bin/sh

# This Makefile builds a C++ application that uses Cantera.  By
# default, the main program file is 'demo.cpp,' which prints out some
# properties of a reacting gas mixture. 

# To build program 'demo', simply type 'make', or 'make -f <this
# file>' if this file is named something other than 'Makefile.'  

# Once you have verified that the demo runs, edit this file to replace
# object file 'demo.o' with your own object file or files. 


#------------------------  edit this block ---------------------------------

# the name of the executable program to be created
PROG_NAME = stagnation.x

# the object files to be linked together. 
OBJS = AxiStagnBVP.o

# additional flags to be passed to the linker. If your program
# requires other external libraries, put them here
LINK_OPTIONS = -L/usr/local/lib -framework Accelerate 

#---------------------------------------------------------------------------
# You probably don't need to edit anything below.

# the C++ compiler
CXX = g++

# C++ compile flags
CXX_FLAGS = -O3 -Wall -fPIC 

# external libraries
EXT_LIBS =  -luser -loneD -lzeroD -lequil -lkinetics -ltransport -lthermo -lctnumerics -lcvode -lctbase  -lctmath -ltpx -lctf2c -lconverters -lctcxx

# Ending C++ linking libraries
LCXX_END_LIBS = -lm

# the directory where the Cantera libraries are located
CANTERA_LIBDIR=/Applications/Cantera/lib

# the directory where Cantera include files may be found.
CANTERA_INCDIR=/Applications/Cantera/include

# flags passed to the C++ compiler/linker for the linking step
LCXXFLAGS = -L$(CANTERA_LIBDIR)  -O3 -Wall -fPIC

# how to compile C++ source files to object files
.cpp.o:
	$(CXX) -c $< -I$(CANTERA_INCDIR) $(CXX_FLAGS)

PROGRAM = $(PROG_NAME)$(EXE_EXT)

DEPENDS = $(OBJS:.o=.d)

all: $(PROGRAM)

$(PROGRAM): $(OBJS)
	$(CXX) -o $(PROGRAM) $(OBJS) $(LCXXFLAGS)\
        $(CANTERA_LIBS) $(LINK_OPTIONS) $(EXT_LIBS) \
          $(LCXX_END_LIBS)

%.d:
	g++ -MM -I$(CANTERA_INCDIR) $*.cpp > $*.d

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






