#!/bin/sh

############################################################################
#
#  Makefile to compile and link a C++ application to 
#  Cantera.
#
#############################################################################

# addition to suffixes
.SUFFIXES : .d

# the name of the executable program to be created
PROG_NAME = air11
#PROG_NAME = air11

# the object files to be linked together. List those generated from Fortran
# and from C/C++ separately
OBJS = air11.o
#OBJS = air11.o

# additional flags to be passed to the linker. If your program
# requires other external libraries, put them here
LINK_OPTIONS = 

#############################################################################
# These links are to Cantera's install space
CANTERA_INCROOT = /home/andrea/WORKSPACE/cantera_lib/include
#
# Bring in the Cantera includes through the .mak file
#
include $(CANTERA_INCROOT)/cantera/Cantera.mak

# Check to see whether we are in the msvc++ environment
os_is_win = 0


# Foqrtran libraries
FORT_LIBS = -lgfortran -lm -lgcc_s -lquadmath
#FORT_LIBS =  -L/home/mpanesi/APPLICATIONS/gcc/gcc-4.6.3-install/lib/gcc/x86_64-unknown-linux-gnu/4.6.3 -L/home/mpanesi/APPLICATIONS/gcc/gcc-4.6.3-install/lib/gcc/x86_64-unknown-linux-gnu/4.6.3/../../../../lib64 -L/lib/../lib64 -L/usr/lib/../lib64 -L/home/mpanesi/APPLICATIONS/gcc/gcc-4.6.3-install/lib/gcc/x86_64-unknown-linux-gnu/4.6.3/../../.. -lgfortran -lm -lgcc_s -lquadmath


# the C++ compiler
CXX = g++
#
# The directory where Cantera include files may be found.
#
INCLUDE_DIRS = -I../../cantera/src -I$(CANTERA_INCROOT)
#
# C++ compile flags
CXX_FLAGS = -g -O0 -Wall   $(INCLUDE_DIRS) -DUSE_VCSNONIDEAL
#CXX_FLAGS = -O3 -Wall   $(INCLUDE_DIRS) -DUSE_VCSNONIDEAL
#
# Ending C++ linking libraries
LCXX_END_LIBS = -lm -lm -lstdc++ 
#LCXX_END_LIBS = -lm -lm -lstdc++ 

# flags passed to the C++ compiler/linker for the linking step
LCXX_FLAGS = -O3 -Wall  

# How to compile C++ source files to object files
.cpp.o:
	$(CXX) $(CXX_FLAGS) -c $< 

# How to compile the dependency file
.cpp.d:
	g++ -MM  $(CXX_FLAGS) $*.cpp > $*.d

# List of dependency files to be created
DEPENDS=$(OBJS:.o=.d)

# Program Name
PROGRAM = $(PROG_NAME)$(EXE_EXT)

all: $(PROGRAM) .depends

$(PROGRAM): $(OBJS) $(CANTERA_CORE_LIBS_DEP)
	$(CXX) -o $(PROGRAM) $(OBJS) $(LCXX_FLAGS) $(LINK_OPTIONS) \
                  $(CANTERA_TOTAL_LIBS) $(LCXX_END_LIBS)

# depends target -> forces recalculation of dependencies
depends:
	$(RM) *.d .depends
	make .depends

.depends: $(DEPENDS)
	cat *.d > .depends

# Do the test -> For the windows vc++ environment, we have to skip checking on
#                whether the program is uptodate, because we don't utilize make
#                in that environment to build programs.
test:
ifeq ($(os_is_win), 1)
else
	make $(PROGRAM)
endif
	./runtest

clean:
	$(RM) $(OBJS) *.o $(PROGRAM) $(DEPENDS) .depends
	(if test -d SunWS_cache ; then \
           $(RM) -rf SunWS_cache ; \
         fi )

