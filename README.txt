        
                             C A N T E R A

                              release 1.5

                               12/14/2004

      Copyright (c) 2001-2004 California Institute of Technology



License information
===================

See the file "License.txt" for information on the terms & conditions
for usage, and a DISCLAIMER OF ALL WARRANTIES.

All trademarks referenced herein are property of their respective
holders.


Web sites
==========

The main Cantera web site is http://www.cantera.org. This primarily serves
as a gateway to the other two web sites:

1. The Cantera User's Group. http://groups.yahoo.com/groups/cantera.   
   This site has a message board, and some miscellaneous files and utilities. 

2. The Cantera Sourceforge site. Distribution of the Cantera source code is
   done using Sourceforge. The site is http://sourceforge.net/projects/cantera.


Installing a Binary Version of Cantera
======================================

Binary installers are available for the Windows and Mac platforms. If
you want to install from one of these, download the appropriate
installer from the Cantera Sourceforge site and run it. This is the
simplest option if you want a standard installation, and plan to
primarily use Cantera from Python or MATLAB. 



Building Cantera from the source code
=====================================


1) Unix/linux/cygwin/Mac OSX build procedure
--------------------------------------------

Run the 'configure' script to build the Makefiles. 

By default, 'make install' will install under '/usr/local' on
linux/unix, and '/Applications/Cantera' on OS X. If you want to
install Cantera somewhere else, run 'configure' with the 'prefix'
option. For example, to install in directory 'cantera' within your
home directory:

configure --prefix=$HOME/cantera

If necessary, edit 'configure' to set options appropriate for your
system before running it.

After running 'configure', type:

make
make install

The last one may need to be run as super-user.

To test the installation, type

make test

After running 'make install', a script 'setup_cantera' will be written
to your home directory. Run it using 'source'
to configure the environment before using Cantera. You may want to add the line

source ~/setup_cantera

to your login script.

The build process requires a 'make' utility compatible with GNU
'make'.  If this has a different name on your system, define
environment variable MAKE to the name (e.g. 'gmake') before running
'configure'.

This procedure also builds the Python and MATLAB interfaces if 
your system is configured to use them. The requirements are:
    -- Python 2.x + NumPy for the Python interface
    -- MATLAB 13 or 14 for the MATLAB toolbox.
If either is missing or an error occurs, the interface is not installed.

Note that your C++ compiler must be compatible with the compiler used
to compile Python or MATLAB. For MATLAB 14, this means you will need
to build Cantera with g++ 3.x, while for MATLAB 13, you will need to
use g++ 2.95.


2) Windows Build Procedure
--------------------------

Cantera can be built under Windows using Visual C++ .NET. See the
document "cantera-vc7" at the Sourceforge site under "Cantera
Documentation / Building and Installing" for more details.


Configuring Matlab
--------------------

The Matlab toolbox uses one compiled MEX program written in C++.
Before you can build it, Matlab needs to be configured
for your compiler. In Matlab type:

mex -setup

and enter the number for the compiler you wish to use.

The Matlab toolbox is built automatically by the Cantera build
process, but can also be built manually.  To build the MEX file needed
for the Matlab toolbox, within Matlab go to to the 'cantera' directory
containing the toolbox and type 'buildux' on unix/linux/Mac OS X, or
'buildwin' on Windows.


Configuring Python
---------------------

Before you can build the Python interface, you need to have Python 2.0
or greater, and the 'numarray' or 'Numeric' package must be
installed. Python is available at www.python.org, and numarray and
Numeric are available through SourceForge (project 'NumPy')


Customizing
-----------

Before running configure, the following environment variables may be set:

MAKE                     set to 'make' utility compatible with GNU make

CXX                      C++ compiler

F77                      Fortran 77 compiler

PYTHON_CMD               Python interpreter to use with Cantera 
                         (default: 'python')

MATLAB_CMD               Matlab command (default: 'matlab')

Additional customization can be done by editing the configure script.

