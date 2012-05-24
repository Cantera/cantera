# SYNOPSIS
#
#   Summarizes configuration settings.
#
#   AX_SUMMARIZE_CONFIG([, ACTION-IF-FOUND [, ACTION-IF-NOT-FOUND]]])
#
# DESCRIPTION
#
#   Outputs a summary of relevant configuration settings.
#

AC_DEFUN([AX_SUMMARIZE_CONFIG],
[

echo
echo '----------------------------------- SUMMARY -----------------------------------'
echo
echo Package version.................   : $PACKAGE-$VERSION
echo
echo C++ compiler....................   : $CXX
echo C++ compiler flags..............   : $CXXFLAGS
echo C compiler......................   : $CC
echo C compiler flags................   : $CFLAGS
echo Fortran compiler................   : $FC
echo Fortran compiler flags..........   : $FCFLAGS
echo Install dir.....................   : $prefix
echo Build user......................   : $USER
echo Build host......................   : $BUILD_HOST
echo Configure date..................   : $BUILD_DATE
echo Build architecture..............   : $BUILD_ARCH
echo SVN revision number.............   : $BUILD_VERSION
echo
echo
echo '-------------------------------------------------------------------------------'

echo
echo Configure complete, now type \'make\' and then \'make install\'.
echo
echo To verify your verification library, type \'make check\'
echo to run a suite of regression tests.

])
