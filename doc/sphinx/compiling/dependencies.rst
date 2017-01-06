
.. _sec-dependencies:

Software used by Cantera
========================

This section lists the versions of third-party software that are required to
build and use Cantera.

Compilers
---------

You must have one of the following C++ compilers installed on your system. A
Fortran compiler is required only if you plan to build the Fortran module.

* GNU compilers (C/C++/Fortran)

  * Known to work with version 4.8; Expected to work with version >= 4.6

* Clang/LLVM (C/C++)

  * Known to work with versions 3.5 and 3.8. Expected to work with version
    >= 3.1.
  * Works with the version included with Xcode 8.2.1.

* Intel compilers (C/C++/Fortran)

  * Known to work with version 14.0.

* Microsoft compilers (C/C++)

  * Known to work with versions 12.0 (Visual Studio 2013) and 14.0 (Visual
    Studio 2015).

* MinGW (C/C++/Fortran)

  * http://mingw-w64.sourceforge.net/ (64-bit and 32-bit)
  * http://tdm-gcc.tdragon.net/ (64-bit and 32-bit)
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

* SUNDIALS

  * If SUNDIALS is not installed, it will be automatically downloaded and the
    necessary portions will be compiled and installed with Cantera.
  * https://computation.llnl.gov/casc/sundials/download/download.html
  * Known to work with versions 2.4, 2.5, 2.6, and 2.7.
  * To use SUNDIALS with Cantera on a Linux/Unix system, it must be compiled
    with the ``-fPIC`` flag. You can specify this flag when configuring
    SUNDIALS (2.4 or 2.5)::

          configure --with-cflags=-fPIC

    or SUNDIALS 2.6 or 2.7::

          cmake -DCMAKE_C_FLAGS=-fPIC <other command-line options>

  .. note:: If you are compiling SUNDIALS 2.5.0 on Windows using CMake, you need
            to edit the ``CMakeLists.txt`` file first and change the lines::

              SET(PACKAGE_STRING "SUNDIALS 2.4.0")
              SET(PACKAGE_VERSION "2.4.0")

            to read::

              SET(PACKAGE_STRING "SUNDIALS 2.5.0")
              SET(PACKAGE_VERSION "2.5.0")

            instead, so that Cantera can correctly identify the version of
            SUNDIALS.

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

* `Numpy <http://www.numpy.org/>`_

  * Required to build the Cantera Python module, and to run significant portions
    of the test suite.
  * Known to work with versions 1.7-1.11; Expected to work with version >= 1.4

* `Cython <http://cython.org/>`_

  * Required version >=0.23 installed for Python 2.7 to build the Python module
    for both Python 2.7 and Python 3.x.

* `3to2 <http://pypi.python.org/pypi/3to2>`_

  * Used to convert Python examples to Python 2 syntax.
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
