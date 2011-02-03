#!/bin/sh
# ------------------------------------------------------------------------------
#
#     Example Cantera pre-preconfig configuration script
#       10/20/09
#
#  Platform:            CYGWIN_NT-5.1 1.7.1(0.218/5/3)
#  Compiler:            gcc v. 4.3.4
#  Optimization:        Debug version
#  python:              v. 2.5.2 (full installation)
#  matlab:              no
#  f2c:                 yes
#  num python:          numpy, which is distributed with cygwin
#
#
# ---------------------------------------------------------------------------
# The only place where the User must custimize this installation
# is located at the top here:
#
#   Specify the installation directory here:
#
Cantera_Install_Dir='/cygdrive/c/cygwin_env/arch/cygwin/cantera-1.8.0'
#
#   Specify the Sundials installation directory here (leave it as a null script
#   if you don't want to link sundials in
#
# Sundials_Home=''
Sundials_Home='/cygdrive/c/cygwin_env/arch/cygwin/sundials-2.3.0_dbg'
#
#   Specify the Sundials version number (either 2.2, 2.3, or 2.4)
Sundials_Version='2.3'
# Sundials_Version=2.4'
# -----------------------------------------------------------------------------
#
CANTERA_CONFIG_PREFIX=$Cantera_Install_Dir
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
#  environment. Therefore, we turn it off here.
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
#  Numpy must have already been installed in the
#              /usr/lib/python2.5/site-packages
#  directory before the Cantera installation. If not, you can always
#  do a minimal python installation and forego forming Cantera's python interface.
#
USE_NUMPY="y"
export USE_NUMPY
#
#      Specify where to find the include files within the numpy distribution
#
NUMPY_INC_DIR="/usr/lib/python2.5/site-packages/numpy/core/include"
export NUMPY_INC_DIR
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
if test -n "$Sundials_Home"
then
  USE_SUNDIALS='y'
  SUNDIALS_VERSION=$Sundials_Version
else
  USE_SUNDIALS='n'
fi
export USE_SUNDIALS
export SUNDIALS_VERSION
#
# Specify where to find the sundials installation directories
#
SUNDIALS_HOME=$Sundials_Home
export SUNDIALS_HOME
#
#  Ok, fire off the main preconfig script
#
./preconfig

