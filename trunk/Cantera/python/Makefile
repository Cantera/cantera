#/bin/sh
###############################################################
#  $Author: hkmoffa $
#  $Date: 2008/01/01 19:30:21 $
#  $Revision: 1.26 $
#
#  Copyright 2001 California Institute of Technology
#  See file License.txt for licensing information
#
###############################################################

have_python_site_package_topdir=0
python_site_package_topdir=/Applications/Cantera

CANTERA_LIBDIR= /Users/dgg/dv/sf/cantera/build/lib/i386-apple-darwin9.6.0
LIB_DEPS =  /Users/dgg/dv/sf/cantera/build/lib/i386-apple-darwin9.6.0/libuser.a /Users/dgg/dv/sf/cantera/build/lib/i386-apple-darwin9.6.0/liboneD.a /Users/dgg/dv/sf/cantera/build/lib/i386-apple-darwin9.6.0/libzeroD.a /Users/dgg/dv/sf/cantera/build/lib/i386-apple-darwin9.6.0/libequil.a /Users/dgg/dv/sf/cantera/build/lib/i386-apple-darwin9.6.0/libkinetics.a /Users/dgg/dv/sf/cantera/build/lib/i386-apple-darwin9.6.0/libtransport.a /Users/dgg/dv/sf/cantera/build/lib/i386-apple-darwin9.6.0/libthermo.a /Users/dgg/dv/sf/cantera/build/lib/i386-apple-darwin9.6.0/libctnumerics.a /Users/dgg/dv/sf/cantera/build/lib/i386-apple-darwin9.6.0/libctmath.a /Users/dgg/dv/sf/cantera/build/lib/i386-apple-darwin9.6.0/libtpx.a /Users/dgg/dv/sf/cantera/build/lib/i386-apple-darwin9.6.0/libctspectra.a /Users/dgg/dv/sf/cantera/build/lib/i386-apple-darwin9.6.0/libconverters.a /Users/dgg/dv/sf/cantera/build/lib/i386-apple-darwin9.6.0/libctbase.a /Users/dgg/dv/sf/cantera/build/lib/i386-apple-darwin9.6.0/libctf2c.a

WIN_LIB_DEPS = $(CANTERA_LIBDIR)/cantera.lib $(CANTERA_LIBDIR)/zeroD.lib \
           $(CANTERA_LIBDIR)/oneD.lib \
           $(CANTERA_LIBDIR)/transport.lib
SRCS = src/ctphase_methods.cpp \
      src/ctthermo_methods.cpp \
      src/ctkinetics_methods.cpp \
      src/cttransport_methods.cpp \
      src/ctxml_methods.cpp \
      src/ctfuncs.cpp \
      src/ctsurf_methods.cpp \
      src/ctbndry_methods.cpp \
      src/ctrpath_methods.cpp \
      src/ctreactor_methods.cpp \
      src/ctfunc_methods.cpp \
      src/ctonedim_methods.cpp \
      src/methods.h

os_is_win = 0
use_clib_dll = 0

all: _build 

win: _winbuild

#
# HKM -> Python seems to want to compile C++ code with the CC compiler.
#        The fix is to assign the CC compiler to the CXX compiler. 
#        the double quotes are needed if the CXX variable is multi-token.
#
_build: $(SRCS) $(LIB_DEPS) Makefile setup.py
	touch src/pycantera.cpp 
	/bin/rm -f _build
	(CXX="g++"; export CXX; CC="g++"; export CC; /usr/bin/python setup.py build)
	echo 'ok' > _build

#
# HKM -> If clib is built as a dll, it needs to be copied by hand to the
#        site package directory.
#
ifeq ($(use_clib_dll), 1)
CLIB_DLL=../../build/lib/i686-pc-win32/clib.dll
CLIB_EXP=../../build/lib/i686-pc-win32/clib.exp
else
CLIB_DLL=
CLIB_EXP=
endif

_winbuild: $(SRCS) $(WIN_LIB_DEPS) $(CLIB_DLL)
	touch src/pycantera.cpp 
	(/usr/bin/python setup.py build)
ifeq ($(use_clib_dll), 1)
ifeq ($(have_python_site_package_topdir), 1)
	(/usr/bin/install -c -m 755 $(CLIB_DLL) $(CLIB_EXP)  \
            $(python_site_package_topdir)/lib/site-packages/Cantera)
endif
endif
	echo 'ok' > _winbuild

minbuild:
	(/usr/bin/python setup.py build)

install:
ifeq (0,1)
	(/usr/bin/python setup.py install --prefix="/Applications/Cantera")
else
	(/usr/bin/python setup.py install)
endif

clean:
	/usr/bin/python setup.py clean
	rm -f _build; rm -f _winbuild
	cd src; rm -f *.o

depends:
	echo '-'

test:
	cd examples; /usr/bin/python reactor1.py &> reactor1.out
	cd examples; /usr/bin/python flame1.py  &> flame1.out
	cd examples; /usr/bin/python diamond.py &> diamond.out
	cd examples; /usr/bin/python critProperties.py  &> critProperties.out

# end of file


