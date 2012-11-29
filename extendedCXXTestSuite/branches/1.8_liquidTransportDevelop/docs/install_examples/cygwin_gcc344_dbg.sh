#!/bin/sh
#
#
CANTERA_ROOT=/cygdrive/c/cygwin_env/cantera-1.8_liquidTransportDevelop
export CANTERA_ROOT
#
CANTERA_INSTALL_DIR=/cygdrive/c/cygwin_env/arch/cygwin/cantera-1.8_develop
export CANTERA_INSTALL_DIR
#
#CADS_ROOT=/home/hkmoffa/Cantera/cads_develop
#export CADS_ROOT
#
PYTHON_DIR=/usr/bin/python
export PYTHON_DIR
#
#  BUILD_WITH_F2C is optional on linux. 
#
BUILD_WITH_F2C='n'
#
#BUILD_WITH_MPI='y'
#
#
#
SET_PYTHON_SITE_PACKAGE_TOPDIR=y
export SET_PYTHON_SITE_PACKAGE_TOPDIR
#
PYTHON_SITE_PACKAGE_TOPDIR=/usr/local/python-2.3.5
export PYTHON_SITE_PACKAGE_TOPDIR
#
AFLAGS='-m32 -Wa,--32'
export AFLAGS
#
OPT='-g'
export OPT
#
do_purify='0'
if test "$do_purify" = 1
then
  PURIFY='purify'
  CXX=g++
  CC=cc
  F77=g77
else
  CXX=g++
  CC=cc
  F77=g77
  CXX_DEPENDS='CC -MM'
fi
export CXX
export F77
export CC
export PURIFY
export CXX_DEPENDS
#
CXXFLAGS="$AFLAGS $OPT -Wall -DUSE_VCSNONIDEAL"
export CXXFLAGS
#
FFLAGS="$AFLAGS $OPT"
export FFLAGS
#
LDFLAGS="$AFLAGS $OPT"
export LDFLAGS
#
LCXX_END_LIBS="-lm -lstdc++"
export LCXX_END_LIBS
#
#  there is a make on cygwin within the MATLAB commands. DO NOT USE THIS VERSION
#
MAKE=make
export MAKE
#
./preconfig

