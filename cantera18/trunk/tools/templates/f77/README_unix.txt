
Unix / linux / cygwin / Mac OS X Build Procedure
================================================

This demo program can be built using any compatible C++ and Fortran 77
compilers. A Makefile 'demo.mak' compatible with the GNU make utility
is provided. This file is produced when Cantera is built, and
specifies the same compilers and options used to build Cantera
itself. The include and library directories are also set to those
locations where these files were installed during the Cantera
installation procedure. If any of these parameters should be changed,
simply edit demo.mak before running it.

To build the demo program, type:

make -f demo.mak

To build your own application, simply edit demo.mak to specify your
source files.
