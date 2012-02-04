# -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ(2.61)
AC_INIT([cantera], [1.8.0], [nick@ices.utexas.edu])
AC_CONFIG_HEADER(config.h)
AC_CONFIG_AUX_DIR([build-aux])
#AX_ENABLE_BUILDDIR
AM_INIT_AUTOMAKE(-Wall -Werror)

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

# add trailing underscore
AC_DEFINE([LAPACK_FTN_TRAILING_UNDERSCORE],1,[Define to add underscore after fortran functions])
LAPACK_FTN_TRAILING_UNDERSCORE=1

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

#                   examples/Makefile                                       \
#		   examples/cxx/Makefile                                   \

AC_OUTPUT(Makefile                                                         \
                   test_problems/Makefile                                  \
                   test_problems/ChemEquil_gri_matrix/Makefile             \
                   test_problems/cathermo/Makefile                         \
                   test_problems/cathermo/DH_graph_1/Makefile              \
                   test_problems/cathermo/DH_graph_acommon/Makefile        \
                   test_problems/cathermo/DH_graph_bdotak/Makefile         \
                   test_problems/cathermo/DH_graph_NM/Makefile             \
                   test_problems/cathermo/DH_graph_Pitzer/Makefile         \
                   test_problems/cathermo/HMW_dupl_test/Makefile           \
                   test_problems/cathermo/HMW_graph_CpvT/Makefile          \
                   test_problems/cathermo/HMW_graph_HvT/Makefile           \
                   test_problems/cathermo/HMW_graph_GvI/Makefile           \
                   test_problems/cathermo/HMW_graph_GvT/Makefile           \
                   test_problems/cathermo/HMW_graph_VvT/Makefile           \
                   test_problems/cathermo/HMW_test_1/Makefile              \
                   test_problems/cathermo/HMW_test_3/Makefile              \
                   test_problems/cathermo/ims/Makefile                     \
                   test_problems/cathermo/issp/Makefile                    \
                   test_problems/cathermo/stoichSubSSTP/Makefile           \
                   test_problems/cathermo/testIAPWS/Makefile               \
                   test_problems/cathermo/testIAPWSPres/Makefile           \
                   test_problems/cathermo/testIAPWSTripP/Makefile          \
                   test_problems/cathermo/testWaterPDSS/Makefile           \
                   test_problems/cathermo/testWaterTP/Makefile             \
                   test_problems/cathermo/VPissp/Makefile                  \
                   test_problems/cathermo/wtWater/Makefile                 \
		   doxygen/Makefile                                        \
		   Cantera/ext/Makefile                                    \
		   Cantera/ext/tpx/Makefile                                \
		   Cantera/ext/cvode/Makefile                              \
		   Cantera/ext/f2c_blas/Makefile                           \
   		   Cantera/ext/f2c_lapack/Makefile                         \
                   Cantera/ext/f2c_libs/Makefile                           \
		   Cantera/ext/f2c_math/Makefile                           \
		   Cantera/Makefile                                        \
                   Cantera/user/Makefile                                   \
                   Cantera/src/Makefile                                    \
                   Cantera/cxx/Makefile                                    \
                   Cantera/cxx/src/Makefile                                \
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

#  ext/lapack/Makefile                                     \


