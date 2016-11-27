.. _sec-compiling:

*************************
Cantera Compilation Guide
*************************

.. toctree::
   :hidden:

   SCons Configuration Options <configuring>

This guide contains instructions for compiling Cantera on the following
operating systems:

* Linux

  * Ubuntu 12.04 LTS (Precise Pangolin) or newer; 16.04 LTS (Xenial Xerus) or
    newer is recommended
  * Debian 7.0 (Wheezy) or newer; 8.0 (Jessie) or newer is recommended
  * Fedora 20 or newer; 22 or newer is recommended

* Windows 7 or newer (32-bit or 64-bit versions)
* OS X 10.9 (Mavericks) or newer; 10.10 (Yosemite) or newer is recommended

In addition to the above operating systems, Cantera should work on any
Unix-like system where the necessary prerequisites are available, but some
additional configuration may be required.

Installation Prerequisites
==========================

Linux
-----

* For Ubuntu or Debian users, the following packages should be installed using
  your choice of package manager::

      g++ python scons libsundials-serial-dev libboost-dev

* Building the python module also requires::

      cython python-dev python-numpy python-numpy-dev python-setuptools

* Building the python3 module requires the following libraries::

      python3 python3-setuptools cython3 python3-numpy

* For Fedora (version 22 or higher) users, the following packages should
  be installed via the package manager::

      gcc-c++ python scons sundials-devel blas-devel lapack-devel boost-devel

  If your Fedora version is lower than 22, the `sundials-devel` package is not
  available, and you should build Sundials from source.

  Python module required packages are::

      python-setuptools python-devel Cython numpy

* Checking out the source code from version control requires Git (install
  ``git``).

* The minimum compatible Cython version is 0.23. If your distribution does not
  contain a suitable version, you may be able to install a more recent version
  using `pip`.

* Building the Fortran interface also requires `gfortran` or another supported
  Fortran compiler.
* Users of other distributions should install the equivalent packages, which
  may have slightly different names.

Windows
-------

There are a number of requirements for the versions of software to install
depending on which interfaces (Python, Matlab) you want to build and what
architecture (32-bit or 64-bit) you want to use. See :ref:`sec-dependencies` for
the full list of dependencies.

* The build process will produce a Python module compatible with the version of
  Python used for the compilation. To generate different modules for other
  versions of Python, you will need to install those versions of Python and
  recompile.
* If you want to build the Matlab toolbox and you have a 64-bit copy of
  Windows, by default you will be using a 64-bit copy of Matlab, and therefore
  you need to compile Cantera in 64-bit mode. For simplicity, it is highly
  recommended that you use a 64-bit version of Python to handle this
  automatically.

* It is generally helpful to have SCons and Python in your PATH. This can
  usually be accomplished by adding the top-level Python directory
  (e.g. ``C:\Python27``) to your PATH. This is accessible from::

      Control Panel > System and Security > System > Advanced System Settings > Environment Variables
* In order to use SCons to install Cantera to a system folder (e.g. ``C:\Program
  Files\Cantera``) you must run the ``scons install`` command in a command
  prompt that has been launched by selecting the *run as administrator* option.

OS X
----

* Download and install Xcode from the App Store

* From a Terminal, run::

      sudo xcode-select --install

  and agree to the Xcode license agreement

* OS X frequently includes out-of-date versions of Numpy (the 1.8.0rc1 version
  included with Sierra has known problems with Cantera). It is recommended to
  install a more recent version (e.g. using ``pip``) before compiling Cantera.

* If you want to build Cantera with Fortran 90 support, download gfortran from::

    http://gcc.gnu.org/wiki/GFortranBinaries#MacOS

* Download scons-2.x.y.tar.gz from scons.org and extract the contents. Install with either::

      sudo python setup.py install

  to install for all users, or::

     python setup.py install --user

  to install to a location in your home directory.

Downloading the Cantera source code
===================================

Stable Release
--------------

* **Option 1**: Check out the code using Git::

    git clone --recursive https://github.com/Cantera/cantera.git
    cd cantera

  Then, check out the tag of the most recent stable version::

    git checkout tags/v2.3.0

  A list of all the tags can be shown by::

    git tag --list

* **Option 2**: Download the most recent source tarball from `Github
  <https://github.com/Cantera/cantera/releases>`_ and extract the
  contents.

Beta Release
------------

* Check out the code using Git::

    git clone --recursive https://github.com/Cantera/cantera.git
    cd cantera

  Then pick either **Option 1** or **Option 2** below.

* **Option 1**: Check out the tag with the most recent beta release::

    git checkout tags/v2.3.0b1

  Note that the most recent beta version might be older than the most recent
  stable release. A list of all the tags, including stable and beta versions can
  be shown by::

    git tag --list

* **Option 2**: Check out the branch with all the bug fixes leading to the
  next minor release of the stable version::

    git checkout 2.3

  This branch has all the work on the 2.3.x version of the software.

Development Version
-------------------

* Check out the code using Git::

    git clone --recursive https://github.com/Cantera/cantera.git
    cd cantera

  Note that by default, the ``master`` branch is checked out, containing all of
  the feature updates and bug fixes to the code since the previous stable
  release. The master branch is usually an "alpha" release, corresponding to the
  ``a`` in the version number, and does not usually get a tag.

* Update an existing clone of the Git repo::

    cd /path/to/cantera
    git fetch
    git rebase origin/master
    git submodule update --init --recursive

Determine configuration options
===============================

General
-------

* run ``scons help`` to see a list all configuration options for Cantera, or
  see :ref:`scons-config`.
* Configuration options are specified as additional arguments to the ``scons``
  command, e.g.::

    scons build -j4 blas_lapack_libs=lapack,blas

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
  **cantera.conf** in the root directory of the source tree, so generally it is
  not necessary to respecify configuration options when rebuilding Cantera. To
  unset a previously set configuration option, either remove the corresponding
  line from cantera.conf or use the syntax::

    option_name=

* Sometimes, changes in your environment can cause SCons's configuration tests
  (e.g. checking for libraries or compiler capabilities) to unexpectedly fail.
  To force SCons to re-run these tests rather than trusting the cached results,
  run scons with the option ``--config=force``.


Python Module
-------------

The Cantera Python module is implemented using Cython, and as such building the
Cantera Python module requires the Cython package for Python.

The Python module is compatible with the following Python versions: 2.7
and 3.2 - 3.5.

Building for Python 2
.....................

By default, SCons will attempt to build the Cython-based Python module for
Python 2, if both Numpy and Cython are installed.

Building for Python 3
.....................

If SCons detects a Python 3 interpreter installed in a default location
(i.e. ``python3`` is on the path), it will try to build the Python module
for Python 3. The following SCons options control how the Python 3 module is
built::

    python3_package=[y|n]
    python3_cmd=/path/to/python3/interpreter
    python3_array_home=/path/to/numpy
    python3_prefix=/path/to/cantera/module

Note that even when building the Python 3 Cantera module, you should still use
Python 2 with SCons, as SCons does not currently support Python 3.

Windows (MSVC)
--------------

* In Windows there aren't any proper default locations for many of the packages
  that Cantera depends on, so you will need to specify these paths explicitly.
* Remember to put double quotes around any paths with spaces in them, e.g.
  "C:\Program Files".
* By default, SCons attempts to use the same architecture as the copy of Python
  that is running SCons, and the most recent installed version of the Visual
  Studio compiler. If you aren't building the Python module, you can override
  this with the configuration options ``target_arch`` and ``msvc_version``.

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
  cannot be used to build a 64-bit Python module. Versions of MinGW that provide
  a 64-bit compiler are available from http://mingw-w64.sourceforge.net/ .

OS X
----

* The Accelerate framework is automatically used to provide optimized versions
  of BLAS and LAPACK, so the ``blas_lapack_libs`` option should generally be
  left unspecified.

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

    scons build env_vars=all CC=icc CXX=icpc FORTRAN=ifort blas_lapack_libs=mkl_rt blas_lapack_dir=$(MKLROOT)/lib/intel64

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

Building Documentation
----------------------

* To build the Cantera HTML documentation, run the commands::

    scons doxygen
    scons sphinx

  or append the options `sphinx_docs=y` and `doxygen_docs=y` to the build
  command, e.g.::

    scons build doxygen_docs=y sphinx_docs=y

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

  * Known to work with version 4.8; Expected to work with version >= 4.6

* Clang/LLVM (C/C++)

  * Known to work with versions 3.5 and 3.8. Expected to work with version
    >= 3.1.
  * Works with the versions included with Xcode 5.1 and Xcode 6.1.

* Intel compilers (C/C++/Fortran)

  * Known to work with version 14.0.

* Microsoft compilers (C/C++)

  * Known to work with versions 12.0 (Visual Studio 2013) and 14.0 (Visual
    Studio 2015).

* MinGW (C/C++/Fortran)

  * http://mingw-w64.sourceforge.net/ (64-bit and 32-bit)
  * Known to work with Mingw-w64 3.0, which provides GCC 4.8. Expected to work
    with any version that provides a supported version of GCC and includes C++11
    thread support.

Other Required Software
-----------------------

* SCons:

  * http://scons.org/tag/releases.html
  * Linux & OS X: Known to work with SCons 2.4.1; Expected to work with versions >= 1.0.0
  * Version 2.3.6 or newer is required to use Visual Studio 2015.

* Python:

  * http://python.org/download/
  * Known to work with 2.7 and 3.5. Expected to work with versions >= 3.3.
  * The Cython module supports Python 2.7 and 3.x. However, SCons requires
    Python 2, so compilation of the Python 3 module requires two Python
    installations.

* Boost

  * http://www.boost.org/users/download/
  * Known to work with version 1.54; Expected to work with versions >= 1.48
  * Only the "header-only" portions of Boost are required. Cantera does not
    currently depend on any of the compiled Boost libraries.

* Sundials

  * If Sundials is not installed, it will be automatically downloaded and the
    necessary portions will be compiled and installed with Cantera.
  * https://computation.llnl.gov/casc/sundials/download/download.html
  * Known to work with versions 2.4, 2.5, 2.6, and 2.7.
  * To use Sundials with Cantera on a Linux/Unix system, it must be compiled
    with the ``-fPIC`` flag. You can specify this flag when configuring
    Sundials (2.4 or 2.5)::

          configure --with-cflags=-fPIC

    or Sundials 2.6 or 2.7::

          cmake -DCMAKE_C_FLAGS=-fPIC <other command-line options>

  .. note:: If you are compiling Sundials 2.5.0 on Windows using CMake, you need
            to edit the ``CMakeLists.txt`` file first and change the lines::

              SET(PACKAGE_STRING "SUNDIALS 2.4.0")
              SET(PACKAGE_VERSION "2.4.0")

            to read::

              SET(PACKAGE_STRING "SUNDIALS 2.5.0")
              SET(PACKAGE_VERSION "2.5.0")

            instead, so that Cantera can correctly identify the version of
            Sundials.

* Eigen

  * If Eigen is not installed, it will be automatically downloaded and installed
    with Cantera.
  * http://eigen.tuxfamily.org/
  * Known to work with version 3.2.8.

* fmt

  * If fmt (previously known as cppformat) is not installed, it will be
    automatically downloaded and the necessary portions will be compiled and
    installed with Cantera.
  * http://fmtlib.net/latest/index.html
  * Version 3.0.1 or newer is required.

* Google Test

  * If Google Test is not installed, it will be automatically downloaded and the
    necessary portions will be compiled as part of the Cantera build process.
  * https://github.com/google/googletest
  * Known to work with version 1.7.0.

Optional Programs
-----------------

* Numpy

  * Required to build the Cantera Python module, and to run significant portions
    of the test suite.
  * http://sourceforge.net/projects/numpy/
  * Known to work with versions 1.7-1.11; Expected to work with version >= 1.4

* `Cython <http://cython.org/>`_

  * Required version >=0.23 installed for Python 2.7 to build the Python module
    for both Python 2.7 and Python 3.x.

* `3to2 <http://pypi.python.org/pypi/3to2>`_

  * Used to convert Cython examples to Python 2 syntax.
  * Known to work with version 1.0

* Matlab

  * Required to build the Cantera Matlab toolbox.
  * Known to work with 2009a and 2014b. Expected to work with versions >= 2009a.

* `Windows Installer XML (WiX) toolset <http://wixtoolset.org/>`_

  * Required to build MSI installers on Windows.
  * Known to work with versions 3.5 and 3.8.

* `Pip <https://pip.pypa.io/en/stable/installing>`_ (Python)

  * Provides the ``pip`` command which can be used to install most of
    the other Python modules.

* Packages required for building Sphinx documentation

  * `Sphinx <http://sphinx.pocoo.org/>`_ (install with ``pip install --upgrade sphinx``)
  * `Pygments <http://pygments.org/>`_ (install with ``pip install --upgrade pygments``)
  * `pyparsing <http://sourceforge.net/projects/pyparsing/>`_ (install with ``pip install --upgrade pyparsing``)
  * `doxylink <http://pypi.python.org/pypi/sphinxcontrib-doxylink/>`_ (install with ``pip install --upgrade sphinxcontrib-doxylink``)
  * `matlabdomain <https://pypi.python.org/pypi/sphinxcontrib-matlabdomain>`_ (install with ``pip install sphinxcontrib-matlabdomain``)

* `Doxygen <http://www.stack.nl/~dimitri/doxygen/>`_

  * Required for building the C++ API Documentation
  * Version 1.8 or newer is recommended.
