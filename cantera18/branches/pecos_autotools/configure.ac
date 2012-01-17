# -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ(2.61)
AC_INIT([cantera], [1.8.0], [nick@ices.utexas.edu])
AC_CONFIG_HEADER(config.h)
AC_CONFIG_AUX_DIR([build-aux])
AM_INIT_AUTOMAKE

# snarf and provide versioning numbers
AX_SPLIT_VERSION
GENERIC_MAJOR_VERSION=$AX_MAJOR_VERSION
GENERIC_MINOR_VERSION=$AX_MINOR_VERSION
GENERIC_MICRO_VERSION=$AX_POINT_VERSION
GENERIC_VERSION=$GENERIC_MAJOR_VERSION.$GENERIC_MINOR_VERSION.$GENERIC_MICRO_VERSION
GENERIC_RELEASE=$GENERIC_MAJOR_VERSION.$GENERIC_MINOR_VERSION

AC_SUBST(GENERIC_VERSION)
AC_SUBST(GENERIC_MAJOR_VERSION)
AC_SUBST(GENERIC_MINOR_VERSION)
AC_SUBST(GENERIC_MICRO_VERSION)
AC_SUBST(GENERIC_RELEASE)

PACKAGE_DESCRIPTION="An object-oriented software toolkit for chemical kinetics, thermodynamics, and transport processes."
AC_SUBST([PACKAGE_DESCRIPTION])
PACKAGE_URL="http://code.google.com/p/cantera/"
AC_SUBST([PACKAGE_URL])

# ------------------------------
# Checks for required programs
# ------------------------------

AC_PROG_CC
AC_PROG_CXX
AC_PROG_FC
AC_PROG_F77

AC_PROG_LIBTOOL
AM_SANITY_CHECK

# ----------------
# Doxygen support 
# ---------------- 

DX_HTML_FEATURE(ON)
DX_CHM_FEATURE(OFF)
DX_CHI_FEATURE(OFF)

DX_MAN_FEATURE(OFF)
DX_RTF_FEATURE(OFF)
DX_XML_FEATURE(OFF)
DX_PDF_FEATURE(ON)
DX_PS_FEATURE(OFF)
DX_INIT_DOXYGEN(cantera, doxygen/cantera.dox, docs) 

# Generate Output Files
AC_OUTPUT(Makefile examples/Makefile                                       \
		   examples/cxx/Makefile                                   \
                   test_problems/Makefile                                  \
		   doxygen/Makefile                                        \
		   ext/Makefile                                            \
		   ext/tpx/Makefile                                        \
		   ext/math/Makefile                                       \
		   ext/blas/Makefile                                       \
		   ext/cvode/Makefile                                      \
		   ext/lapack/Makefile                                     \
		   Cantera/Makefile                                        \
                   Cantera/user/Makefile                                   \
                   Cantera/src/Makefile                                    \
                   Cantera/src/base/Makefile                               \
                   Cantera/src/converters/Makefile                         \
		   Cantera/src/kinetics/Makefile                           \
                   Cantera/src/numerics/Makefile Cantera/src/oneD/Makefile \
                   Cantera/src/thermo/Makefile                             \
		   Cantera/src/transport/Makefile                          \
                   Cantera/src/spectra/Makefile Cantera/src/zeroD/Makefile \
                   Cantera/src/equil/Makefile)

# FINAL SUMMARY
AX_SUMMARIZE_CONFIG 
