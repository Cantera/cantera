#!/bin/sh
# ------------------------------------------------------------------------------
#
#     Example Cantera pre-preconfig configuration script
#
#  Platform:            cygwin
#  Compiler:            gcc
#  python:              v. 2.5.2 (full installation)
#  matlab:              no
#  f2c:                 yes
#
#
# Specific places where the User must custimize the installation
# is indiciated with the USER_INPUT_NEEDED lines.
#
#
#  Specify the installation directory here:
#  USER_INPUT_NEEDED
#
# -----------------------------------------------------------------------------
#
CANTERA_CONFIG_PREFIX=/cygdrive/c/cygwin_env/arch/cygwin/cantera-1.8_develop
export CANTERA_CONFIG_PREFIX

SET_PYTHON_SITE_PACKAGE_TOPDIR=y
export SET_PYTHON_SITE_PACKAGE_TOPDIR

PYTHON_SITE_PACKAGE_TOPDIR=$CANTERA_CONFIG_PREFIX
export PYTHON_SITE_PACKAGE_TOPDIR

PYTHON_CMD=/usr/bin/python
export PYTHON_CMD

#
# The full python package works within cygwin. Minimal is needed at the least
#
PYTHON_PACKAGE='full'
#PYTHON_PACKAGE='minimal'
export PYTHON_PACKAGE

WITH_IDEAL_SOLUTIONS="y"
export WITH_IDEAL_SOLUTIONS

WITH_ELECTROLYTES="y"
export WITH_ELECTROLYTES

WITH_VCSNONIDEAL="y"
export WITH_VCSNONIDEAL

#
#  Matlab on the pc is supported through the windows vc++ 8.0 compilation
#  environment
#
BUILD_MATLAB_TOOLBOX="n"
export BUILD_MATLAB_TOOLBOX

NUMARRAY_HOME=''
export NUMARRAY_HOME
#
# Turns on extra debugging
#
DEBUG_MODE='y'
export DEBUG_MODE
#
#  Numeric must have already been installed in the
#              /usr/lib/python2.5/site-packages
#  directory before the Cantera installation. If not, you can always
#  do a minimal python installation and forego forming Cantera's python interface.
#
USE_NUMERIC="y"
export USE_NUMERIC
#
#  We use F2C with the cygwin installation. However, the g77 interface should 
#  work too.
#
BUILD_WITH_F2C="y"
export BUILD_WITH_F2C

BITCOMPILE="32"
export BITCOMPILE
#
#  We use the default gcc compilers that come with cygwin. Currently, this
#  is gcc v. 3.4.4.
#
CXX='g++'
export CXX
#
CC='gcc'
export CC
#
F77='g77'
export F77
#
# Debug options
#
# -DDEBUG_HKM
# -DDEBUG_BASISOPTIMIZE
# -DDEBUG_CHEMEQUIL
#
CXXFLAGS="-g -Wall "
export CXXFLAGS

CFLAGS=$CXXFLAGS
export CFLAGS

FFLAGS="-g "
export FFLAGS

LDFLAGS='-g'
export LDFLAGS

LCXX_END_LIBS="-lm -lstdc++"
export LCXX_END_LIBS

# The gcc compiler in cygwin will put a .exe extension on anyway. So,
# we might as well define this field as .exe. cygwin will automagically
# equate executable files with no extension as having an .exe extension.
# This is very mysterious. However, the mystery goes away if you define
# these files to have a .exe suffix from the getgo.
EXE_EXT=".exe"
export EXE_EXT

EXTRA_LINK=""
export EXTRA_LINK

MAKE=make
export MAKE
#
# Specify the SUNDIALS option
#
USE_SUNDIALS='n'
export USE_SUNDIALS
#
# Specify where to find the sundials installation directories
#
SUNDIALS_HOME='/cygdrive/c/cygwin_env/libraries/sundials'
export SUNDIALS_HOME
#
#  Ok, fire off the main preconfig script
#
./preconfig

