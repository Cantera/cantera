

#  WORK IN PROGRESS
#  This configuration file is not yet functional; use the "preconfig" 
#  script instead.


#######################################################################
#
#           Cantera Configuration File
#
#  Edit this file to control how Cantera is built. Parameters can be set
#  here, or alternatively environment variables may be set before calling
#  this script.
#
#  The default configuration uses GNU compilers (gcc/g++/g77) and
#  builds as much of Cantera and its language interfaces as it can
#  (e.g. if MATLAB is installed on your system, the MATLAB toolbox
#  will be built automatically, otherwise it will be skipped. On linux
#  or Mac OS X, this default configuration should work, and most
#  likely you don't need to edit this file at all - just run it.
#  
#  NOTE: if you DO make changes to this file, save it with another name
#  so that it will not be overwritten if you update the source
#  distribution. 

#######################################################################

#----------------------------------------------------------------------
#         Language Interfaces
#----------------------------------------------------------------------
#
# Cantera has several programming language interfaces. Select the ones
# you want to build. The default is to try to build all language
# interfaces.
#
#
#----------------- Python --------------------------------------------
# 
# In addition to being one of the supported language interfaces,
# Python is used internally by Cantera, both in the build process and
# at run time (to process .cti input files). Therefore, you generally
# need to have Python on your system; if you don't, first install it
# from http://www.python.org before proceeding with the installation
# of Cantera.
#
# If you plan to work in Python, or you want to use the graphical
# MixMaster application, then you need the full Cantera Python
# Package. If, on the other hand, you will only use Cantera from some
# other language (e.g. MATLAB or Fortran 90/95) and only need Python
# to process .cti files, then you only need a minimal subset of the
# package (actually, only one file).

# Set PYTHON_PACKAGE to one of these four strings:
#    full      install everything needed to use Cantera from Python
#    minimal   install only enough to process .cti files
#    none      Don't install  or run any Python scripts during the 
#              build process
#    default   try to do a full installation, but fall back to a minimal
#              one in case of errors
OPTION(CANTERA_BUILD_PYTHON_PACKAGE "Build the Python Package?" ON)

SET(CANTERA_PYTHON_PACKAGE_TYPE "default" CACHE LIST "full or minimal")

# Cantera needs to know where to find the Python interpreter.  If
# PYTHON_CMD is set to "default", then cmake will look for the Python
# interpreter
SET( PYTHON_CMD "default")

# The Cantera Python interface can be built with either the numarray
# or Numeric packages. Set this to "y" to use Numeric, or anything
# else to use numarray. Using numarray is preferred.
SET( USE_NUMERIC "default")

# If numarray was installed using the --home option, set this to the
# home directory for numarray.
SET( NUMARRAY_HOME "$HOME/python_packages")

SET( CANTERA_VERSION "1.7.1")

SET( WITH_PURE_FLUIDS  1 )

SET( WITH_LATTICE_SOLID 0 )

SET( WITH_METAL 1 )

SET( WITH_STOICH_SUBSTANCE 1 )

SET ( WITH_IDEAL_SOLUTIONS 1 )

SET (WITH_ELECTROLYTES 1 )

