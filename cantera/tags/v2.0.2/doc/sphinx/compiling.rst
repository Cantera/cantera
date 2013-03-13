*************************
Cantera Compilation Guide
*************************

.. toctree::
   :hidden:

   SCons Configuration Options <configuring>

This guide contains instructions for compiling Cantera on the following
operating systems:

* Linux

  * Ubuntu 10.04 LTS (Lucid Lynx) or newer
  * Debian 5.0 (Lenny) or newer

* Windows Vista or Windows 7 (32-bit or 64-bit versions)
* OS X 10.6 (Snow Leopard) or newer

In addition to the above operating systems, Cantera should work on any
Unix-like system where the necessary prerequisites are available, but some
additional configuration may be required.

Installation Prerequisites
==========================

Linux
-----

* For Ubuntu or Debian users, the following packages should be installed using
  your choice of package manager::

      g++ python scons libboost-all-dev libsundials-serial-dev subversion

* Building the python module also requires::

      python-dev python-numpy python-numpy-dev

* Building the Fortran interface also requires gfortran or another supported
  Fortran compiler.
* Users of other distributions should install the equivalent packages, which
  may have slightly different names.

Windows
-------

There are a number of requirements for the versions of software to install
depending on which interfaces (Python, Matlab) you want to build and what
architecture (32-bit or 64-bit) you want to use. See :ref:`sec-dependencies` for
the full list of dependencies.

* If you want to build the Python module, you must use the same version of the
  Microsoft compiler as was used to compile Python. For current versions of
  Python (2.6 and 2.7) this means that you must use Visual Studio 2008 or the
  equivalent version of the Windows SDK (see link below).
* The build process will produce a Python module compatible with the version of
  Python used for the compilation. To generate different modules for other
  versions of Python, you will need to install those versions of Python and
  recompile.
* If you want to build the Matlab toolbox and you have a 64-bit copy of
  Windows, by default you will be using a 64-bit copy of Matlab, and therefore
  you need to compile Cantera in 64-bit mode. For simplicity, it is highly
  recommended that you use a 64-bit version of Python to handle this
  automatically.
* There is no 64-bit installer for SCons under Windows, so you will need to
  download the ZIP version. After extracting it, start a command prompt in the
  unzipped folder and run::

     python setup.py install

* It is generally helpful to have SCons and Python in your PATH. This can
  usually be accomplished by adding the top-level Python directory
  (e.g. ``C:\Python27``) to your PATH. This is accessible from::

      Control Panel > System and Security > System > Advanced System Settings > Environment Variables

OS X
----

* Download and install Xcode from the
  `Apple Developer site <https://developer.apple.com/xcode/index.php>`_
* Cantera can be compiled with the command line tools that ship with either
  Xcode 3.x or Xcode 4.x.
* If you don't have numpy version >= 1.3, you can install a recent version with::

    sudo easy_install -U numpy

* If you want to build Cantera with Fortran 90 support, download gfortran from::

    http://gcc.gnu.org/wiki/GFortranBinaries#MacOS

* Download scons-2.1.0.tar.gz from scons.org and extract the contents. Install with either::

      sudo python setup.py install

  to install for all users, or::

     python setup.py install --user

  to install to a location in your home directory.

Downloading the Cantera source code
===================================

* Option 1: Download the most recent source tarball from the `Cantera Google
  Code Site <http://code.google.com/p/cantera/downloads/list>`_ and extract the
  contents.

* Option 2: Check out the code using Subversion::

    svn checkout http://cantera.googlecode.com/svn/cantera/trunk/ cantera

* Option 3: Check out the code using Git::

    git svn clone --std-layout http://cantera.googlecode.com/svn/cantera cantera

Determine configuration options
===============================

General
-------

* run ``scons help`` to see a list all configuration options for Cantera, or
  see :ref:`scons-config`.
* If the prerequisites are installed in standard locations, the default values
  should work.
* If you installed Sundials to a non-standard location (e.g. the libraries
  aren't in /usr/lib), you will need to specify the options::

    sundials_include=/path/to/sundials/include
    sundials_libdir=/path/to/sundials/lib

* If you want to build the Matlab toolbox, you will need to specify the path
  to the Matlab installation, e.g.::

    matlab_path=/opt/MATLAB/R2011a
    matlab_path="C:\Program Files\MATLAB\R2011a"
    matlab_path=/Applications/MATLAB_R2011a.app

  The above paths are typical defaults on Linux, Windows, and OS X,
  respectively.
* SCons saves configuration options specified on the command line in the file
  \b cantera.conf in the root directory of the source tree, so generally it is
  not necessary to respecify configuration options when rebuilding Cantera. To
  unset a previously set configuration option, either remove the corresponding
  line from cantera.conf or use the syntax::

    option_name=


Windows (MSVC)
--------------

* In Windows there aren't any proper default locations for many of the packages
  that Cantera depends on, so you will need to specify these paths explicitly.
* Remember to put double quotes around any paths with spaces in them, e.g.
  "C:\Program Files".
* By default, SCons attempts to use the same architecture and version of the
  Microsoft compiler as was used to compile Python, typically Visual Studio
  2008 or the equivalent version of the Windows SDK. If you aren't building the
  Python module, you can override this with the configuration options
  ``target_arch`` and ``msvc_version``.

.. note::

    The ``cantera.conf`` file uses the backslash character ``\`` as an escape
    character. When modifying this file, backslashes in paths need to be escaped
    like this: ``boost_inc_dir = 'C:\\Program Files (x86)\\boost\\include'``
    This does not apply to paths specified on the command line. Alternatively,
    you can use forward slashes in paths.


Windows (MinGW)
---------------

* To compile with MinGW, use the SCons command line option::

    toolchain=mingw

* The version of MinGW from http://www.mingw.org is 32-bit only, and therefore
  cannot be used to build a 64-bit Python module.

OS X
----

* The available compilers to compile Cantera will depend on the version of
  Xcode that is installed.

  * If Xcode 3 is installed, you can use either GCC by leaving the ``CC`` and
    ``CXX`` options unspecified, or setting them to::

      CC=gcc CXX=g++

    You can also use LLVM with the GCC frontend by specifying::

      CC=llvm-gcc CXX=llvm-g++

  * If Xcode 4 is installed, then you can either use LLVM-GCC as above or
    Clang by specifying::

      CC=clang CXX=clang++

* The Accelerate framework provides optimized versions of BLAS and LAPACK, so
  the ``blas_lapack_libs`` option should generally be left unspecified.

Intel Compilers
---------------
* Before compiling Cantera, you may need to set up the appropriate environment
  variables for the Intel compiler suite, e.g.::

    source /opt/intel/bin/compilervars.sh intel64

* For the Intel compiler to work with SCons, these environment variables need
  to be passed through SCons by using the command line option::

    env_vars=all

* If you want to use the Intel MKL versions of BLAS and LAPACK, you will need
  to provide additional options. The following are typically correct on
  64-bit Linux systems::

    blas_lapack_libs=mkl_rt blas_lapack_dir=$(MKLROOT)/lib/intel64

  Your final SCons call might then look something like::

    scons build env_vars=all CC=icc CXX=icpc F90=ifort F77=ifort blas_lapack_libs=mkl_rt blas_lapack_dir=$(MKLROOT)/lib/intel64

  When installing Cantera after building with the Intel compiler, the normal
  method of using ``sudo`` to install Cantera will not work because ``sudo``
  does not pass the environment variables needed by the Intel compiler.
  Instead, you will need to do something like::

    scons build ...
    sudo -s
    source /path/to/compilervars.sh intel64
    scons install
    exit

Compile Cantera & Test
======================

* Run scons with the list of desired configuration options, e.g.::

    scons build optimize=n blas_lapack_libs=blas,lapack prefix=/opt/cantera

* If Cantera compiles successfully, you should see a message that looks like::

    *******************************************************
    Compilation completed successfully.

    - To run the test suite, type 'scons test'.
    - To install, type '[sudo] scons install'.
    *******************************************************

* If you do not see this message, check the output for errors to see what went
  wrong.
* Cantera has a series of tests that can be run with the command::

    scons test

* When the tests finish, you should see a summary indicating the number of
  tests that passed and failed.

* If you have tests that fail, try looking at the following to determine the
  source of the error:

    * Messages printed to the console while running scons test
    * Output files generated by the tests

MinGW Compilation problems
--------------------------

* If you get a compiler error while compiling some of the "f2c" code, then your
  version of MinGW has a problem with the order of its internal include paths,
  such that it sees the GCC float.h before its own special version. To fix this
  problem edit the GCC float.h located at (roughly)::

    c:\MinGW\lib\gcc\mingw32\4.6.1\include\float.h

  and add the following just before the end (before the final #endif)

  .. code-block:: c++

      #ifndef _MINGW_FLOAT_H_
      #include_next <float.h>
      #endif

.. _sec-dependencies:

Software used by Cantera
========================

This section lists the versions of third-party software that are required to
build and use Cantera.

Compilers
---------

You must have one of the following C++ compilers installed on your system. A
Fortran compiler is required only if you plan to use Cantera from a Fortran
program.

* GNU compilers (C/C++/Fortran)

  * Known to work with version 4.6; Expected to work with version >= 4.3

* Clang/LLVM (C/C++)

  * Known to work with version 2.9
  * This is the version included with Apple Xcode 4.x

* Intel compilers (C/C++/Fortran)

  * Known to work with version 11.0 and 12.1; Expected to work with
    versions >= 11.0

* Microsoft compilers (C/C++)

  * Windows SDK: http://www.microsoft.com/download/en/details.aspx?id=3138
    This is equivalent to Visual Studio 2008. It is a free download.
  * Known to work with version 9.0 (Visual Studio 2008) and version 10.0
    (Visual Studio 2010).

* MinGW (C/C++/Fortran)

  * http://www.mingw.org
  * Known to work with version 4.6.
  * Supported versions of MinGW should be the same as the supported versions of
    GCC.

Other Required Software
-----------------------

* Subversion

  * For Windows: http://tortoisesvn.net/downloads.html
  * Known to work with versions >= 1.6

* SCons:

  * http://www.scons.org/download.php
  * Known to work with SCons 2.1.0; Expected to work with versions >= 1.0.0

* Python:

  * http://python.org/download/
  * Known to work with 2.6 and 2.7; Expected to work with versions >= 2.5
  * Does not yet work with 3.x

* Boost

  * http://www.boost.org/users/download/
  * Known to work with version 1.46; Expected to work with versions >= 1.40
  * Only the "header-only" portions of Boost are required. Cantera does not
    currently depend on any of the compiled Boost libraries.

Optional Programs
-----------------

* Numpy

  * Required to build the Cantera Python module.
  * http://sourceforge.net/projects/numpy/
  * Known to work with versions 1.3 and 1.5; Expected to work with version >= 1.1
  * Test suite requires version >= 1.3

* Matlab

  * Required to build the Cantera Matlab toolbox.
  * Known to work with 2009a, 2010a, and 2011b. Expected to work with
    versions >= 2009a.

* Sundials

  * Required to enable some features such as sensitivity analysis.
  * Strongly recommended if using reactor network or 1D simulation capabilities.
  * https://computation.llnl.gov/casc/sundials/download/download.html
  * Known to work with versions 2.4 and 2.5; Support for versions 2.3
    and 2.2 is deprecated.
  * To use Sundials with Cantera, you may need to compile it with the
    ``-fPIC`` flag. You can specify this flag when configuring Sundials::

          configure --with-cflags=-fPIC

* `Windows Installer XML (WiX) toolset <http://wix.sourceforge.net/>`_

  * Required to build MSI installers on Windows.
  * Known to work with version 3.5.

* `Distribute <http://pypi.python.org/pypi/distribute>`_ (Python)

  * Provides the ``easy_install`` command which can be used to install most of
    the other Python modules.

* Packages required for building Sphinx documentation

  * `Sphinx <http://sphinx.pocoo.org/>`_ (install with ``easy_install -U Sphinx``)
  * `Pygments <http://pygments.org/>`_ (install with ``easy_install -U pygments``)
  * `pyparsing <http://sourceforge.net/projects/pyparsing/>`_ (install with ``easy_install -U pyparsing``)

* `Doxygen <http://www.stack.nl/~dimitri/doxygen/>`_

  * Required for building the C++ API Documentation
