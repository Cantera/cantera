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
PROG_NAME = wtWater

# the object files to be linked together. List those generated from Fortran
# and from C/C++ separately
OBJS = wtWater.o

# additional flags to be passed to the linker. If your program
# requires other external libraries, put them here
LINK_OPTIONS =  

#############################################################################

# Check to see whether we are in the msvc++ environment
os_is_win = 0

# Fortran libraries
ifeq (0, 0)
FORT_LIBS = 
else
FORT_LIBS =
endif

# the C++ compiler
CXX = /home/sntools/extras/compilers/gcc-4.2.4/bin/gcc

# C++ compile flags
CXX_FLAGS = -g  -Wall -DDEBUG_HKM -DDEBUG_HKM_EPEQUIL  

# Ending C++ linking libraries
LCXX_END_LIBS = -lctf2c -lgfortran  -lm -lstdc++

# the directory where the Cantera libraries are located
CANTERA_LIBDIR=/ascldap/users/hkmoffa/Cantera/gc/canteraDevelop/build/lib/x86_64-unknown-linux-gnu

# required Cantera libraries
CANTERA_LIBS =  -luser -loneD -lzeroD -lequil -lVCSnonideal -lkinetics -ltransport -lthermo -lctnumerics -lctmath -ltpx -lctspectra -lconverters -lctbase -lsundials_cvodes -lsundials_nvecserial -lctlapack -lctblas -lctf2c -lctcxx

# Dependencies for CANTERA_LIBS
CANTERA_LIBS_DEP =  /ascldap/users/hkmoffa/Cantera/gc/canteraDevelop/build/lib/x86_64-unknown-linux-gnu/libuser.a /ascldap/users/hkmoffa/Cantera/gc/canteraDevelop/build/lib/x86_64-unknown-linux-gnu/liboneD.a /ascldap/users/hkmoffa/Cantera/gc/canteraDevelop/build/lib/x86_64-unknown-linux-gnu/libzeroD.a /ascldap/users/hkmoffa/Cantera/gc/canteraDevelop/build/lib/x86_64-unknown-linux-gnu/libequil.a /ascldap/users/hkmoffa/Cantera/gc/canteraDevelop/build/lib/x86_64-unknown-linux-gnu/libVCSnonideal.a /ascldap/users/hkmoffa/Cantera/gc/canteraDevelop/build/lib/x86_64-unknown-linux-gnu/libkinetics.a /ascldap/users/hkmoffa/Cantera/gc/canteraDevelop/build/lib/x86_64-unknown-linux-gnu/libtransport.a /ascldap/users/hkmoffa/Cantera/gc/canteraDevelop/build/lib/x86_64-unknown-linux-gnu/libthermo.a /ascldap/users/hkmoffa/Cantera/gc/canteraDevelop/build/lib/x86_64-unknown-linux-gnu/libctnumerics.a /ascldap/users/hkmoffa/Cantera/gc/canteraDevelop/build/lib/x86_64-unknown-linux-gnu/libctmath.a /ascldap/users/hkmoffa/Cantera/gc/canteraDevelop/build/lib/x86_64-unknown-linux-gnu/libtpx.a /ascldap/users/hkmoffa/Cantera/gc/canteraDevelop/build/lib/x86_64-unknown-linux-gnu/libctspectra.a /ascldap/users/hkmoffa/Cantera/gc/canteraDevelop/build/lib/x86_64-unknown-linux-gnu/libconverters.a /ascldap/users/hkmoffa/Cantera/gc/canteraDevelop/build/lib/x86_64-unknown-linux-gnu/libctbase.a /ascldap/users/hkmoffa/Cantera/gc/canteraDevelop/build/lib/x86_64-unknown-linux-gnu/libctf2c.a $(CANTERA_LIBDIR)/libctcxx.a

# the directory where Cantera include files may be found.
CANTERA_INCDIR=/ascldap/users/hkmoffa/Cantera/gc/canteraDevelop/build/include/cantera
INCLUDES=-I$(CANTERA_INCDIR) -I$(CANTERA_INCDIR)/kernel

# flags passed to the C++ compiler/linker for the linking step
LCXX_FLAGS = -L$(CANTERA_LIBDIR)  -L/ascldap/users/hkmoffa/Cantera/gc/canteraDevelop/build/lib/x86_64-unknown-linux-gnu -L/home/hkmoffa/arch/linux64/sundials/lib -g  -Wall -DDEBUG_HKM -DDEBUG_HKM_EPEQUIL  

# How to compile C++ source files to object files
.cpp.o:
	$(CXX) -c $< $(INCLUDES) $(CXX_FLAGS)

# How to compile the dependency file
.cpp.d:
	/home/sntools/extras/compilers/gcc-4.2.4/bin/g++ -MM $(INCLUDES) $(CXX_FLAGS) $*.cpp > $*.d

# List of dependency files to be created
DEPENDS=$(OBJS:.o=.d)

# Program Name
PROGRAM = $(PROG_NAME)$(EXE_EXT)

all: $(PROGRAM) .depends

$(PROGRAM): $(OBJS) $(CANTERA_LIBS_DEP)
	$(CXX) -o $(PROGRAM) $(OBJS) $(LCXX_FLAGS) $(LINK_OPTIONS) \
                  $(CANTERA_LIBS)  $(FORT_LIBS) \
                  $(LCXX_END_LIBS)

# Add an additional target for stability:
$(OBJS):  $(CANTERA_LIBDIR)/libctbase.a $(CANTERA_LIBDIR)/libthermo.a

# depends target -> forces recalculation of dependencies
depends:
	gmake .depends

.depends: $(DEPENDS)
	cat $(DEPENDS) > .depends

# Do the test -> For the windows vc++ environment, we have to skip checking on
#                whether the program is uptodate, because we don't utilize make
#                in that environment to build programs.
test:
ifeq ($(os_is_win), 1)
else
	@ gmake -s $(PROGRAM)
endif
	@ ./runtest

clean:
	$(RM) $(OBJS) $(PROGRAM) $(DEPENDS) .depends
	../../../bin/rm_cvsignore
	(if test -d SunWS_cache ; then \
           $(RM) -rf SunWS_cache ; \
         fi )

ifeq ($(wildcard .depends), .depends)
include .depends
endif

