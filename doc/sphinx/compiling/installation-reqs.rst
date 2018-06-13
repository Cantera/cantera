
.. contents::
   :local:

.. _sec-installation-reqs:

Installation Prerequisites
==========================

.. _sec-linux:

Linux
-----

General Notes
^^^^^^^^^^^^^

* To download the source code, installing ``git`` is highly recommended.

* SCons is only available for Python 2, so building the Python 3 module requires
  two installations of Python (one of Python 2 and one of Python 3), even if you
  do not intend to build the Python 2 module.

* The following instructions use the system-installed versions of Python, but
  alternate installations such as the Anaconda distribution of Python can be
  used as well.

* Cython is only required to be installed for the version of Python that also
  has SCons installed; following the instructions below will install Cython for
  the version of Python 2 installed in the system directories. The minimum
  compatible Cython version is 0.23. If your distribution does not contain a
  suitable version, you may be able to install a more recent version using
  ``pip``.

* Users of other distributions should install the equivalent packages, which
  may have slightly different names.

* In addition to the below operating systems, Cantera should work on any
  Unix-like system where the necessary prerequisites are available, but some
  additional configuration may be required.

.. _sec-ubuntu-debian-reqs:

Ubuntu & Debian
^^^^^^^^^^^^^^^

* Ubuntu 12.04 LTS (Precise Pangolin) or newer is required; 16.04 LTS (Xenial Xerus)
  or newer is recommended

* Debian 7.0 (Wheezy) or newer; 8.0 (Jessie) or newer is recommended

* The following packages must be installed to build any of the Cantera modules using
  your choice of package manager::

      g++ python scons libboost-dev

* In addition to the general packages, building the Python 2 module also requires::

      cython python-dev python-numpy python-numpy-dev python-setuptools

* In addition to the general packages, building the Python 3 module also requires::

      cython python3 python3-dev python3-setuptools python3-numpy

* In addition to the general packages, building the Fortran module also requires::

      gfortran

* In addition to the general packages, building the MATLAB toolbox also requires:

  * MATLAB version later than 2009a

    * Typically installed to::

        /opt/MATLAB/R20YYn

      where ``YY`` is a two digit year and ``n`` is either ``a`` or ``b``

.. _sec-fedora-reqs:

Fedora & RHEL
^^^^^^^^^^^^^

* The following packages must be installed to build any of the Cantera modules using
  your choice of package manager::

      gcc-c++ python scons boost-devel

* In addition to the general packages, building the Python 2 module also requires::

      python-setuptools python-devel Cython numpy

* In addition to the general packages, building the Python 3 module also requires::

      python3 python3-setuptools python3-devel Cython python3-numpy

* In addition to the general packages, building the Fortran module also requires::

      gcc-gfortran

* In addition to the general packages, building the MATLAB toolbox also requires:

  * MATLAB version later than 2009a

    * Typically installed to::

        /opt/MATLAB/R20YYn

      where ``YY`` is a two digit year and ``n`` is either ``a`` or ``b``

.. _sec-opensuse-reqs:

OpenSUSE & SUSE Linux Enterprise
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

* OpenSUSE 13.2 or newer; Leap 42.2 or newer recommended

* The following packages must be installed to build any of the Cantera modules using
  your choice of package manager::

      gcc-c++ python scons boost-devel

* In addition to the general packages, building the Python 2 module also requires::

      python-Cython python-devel python-numpy python-numpy-devel python-setuptools

* In addition to the general packages, building the Python 3 module also requires::

      python-Cython python3 python3-devel python3-setuptools python3-numpy python3-numpy-devel

* In addition to the general packages, building the Fortran module also requires::

      gcc-fortran

* In addition to the general packages, building the MATLAB toolbox also requires:

  * MATLAB version later than 2009a

    * Typically installed to::

        /opt/MATLAB/R20YYn

      where ``YY`` is a two digit year and ``n`` is either ``a`` or ``b``

.. _sec-windows:

Windows
-------

General Notes
^^^^^^^^^^^^^

* SCons is only available for Python 2, so building the Python 3 module requires
  two installations of Python (one of Python 2 and one of Python 3), even if you
  do not intend to build the Python 2 module.

* The build process will produce a Python module compatible with the version of
  Python used for the compilation. To generate different modules for other
  versions of Python, you will need to install those versions of Python and
  recompile.

* The following instructions use the versions of Python downloaded from
  https://www.python.org/downloads, but alternate installations such as the
  Anaconda distribution of Python can be used as well.

* If you want to build the Matlab toolbox and you have a 64-bit copy of Windows,
  by default you will be using a 64-bit copy of Matlab, and therefore you need
  to compile Cantera in 64-bit mode. For simplicity, it is highly recommended
  that you use a 64-bit version of Python to handle this automatically. Note
  that the default download from the Python website
  (https://www.python.org) is for a 32-bit installer, and you will
  need to select the 64-bit installer specifically.

* It is generally helpful to have SCons and Python in your ``PATH`` environment
  variable. This can be done by checking the appropriate box during the
  installation of Python or can be accomplished by adding the top-level Python
  directory and the ``Scripts`` subdirectory (e.g.,
  ``C:\Python27;C:\Python27\Scripts``) to your ``PATH``. The dialog to change
  the ``PATH`` is accessible from::

      Control Panel > System and Security > System > Advanced System Settings > Environment Variables

  Make sure that the installation of Python that has SCons comes first on your
  ``PATH``.

* In order to use SCons to install Cantera to a system folder (e.g. ``C:\Program
  Files\Cantera``) you must run the ``scons install`` command in a command
  prompt that has been launched by selecting the *Run as Administrator* option.

.. _sec-windows-reqs:

Windows Requirements
^^^^^^^^^^^^^^^^^^^^^^^

* Windows 7 or later; either 32-bit or 64-bit

* To build any of the Cantera modules, you will need to install

  * Python 2.7

    * https://www.python.org/downloads/

    * Be sure to choose the appropriate architecture for your system - either
      32-bit or 64-bit

    * When installing, make sure to choose the option to add to your ``PATH``

  * SCons

    * https://pypi.python.org/pypi/SCons

    * Be sure to choose the appropriate architecture for your system - either
      32-bit or 64-bit

  * One of the following supported compilers

    * Microsoft compilers

      * https://www.visualstudio.com/downloads/

      * Known to work with Visual Studio 2013 (MSVC 12.0) and Visual Studio 2015
        (MSVC 14.0)

    * MinGW compilers

      * http://mingw-w64.org/

      * http://tdm-gcc.tdragon.net/

      * Known to work with Mingw-w64 3.0, which provides GCC 4.8. Expected to
        work with any version that provides a supported version of GCC and
        includes C++11 thread support.

      * The version of MinGW from http://www.mingw.org/ cannot be used to build
        Cantera. Users must use MinGW-w64 or TDM-GCC.

  * The Boost headers

    * http://www.boost.org/doc/libs/1_63_0/more/getting_started/windows.html#get-boost

    * It is not necessary to compile the Boost libraries since Cantera only uses
      the headers from Boost

* In addition to the general software, building the Python 2 module also requires

  * Pip

    * Pip should be distributed with Python version 2.7.9 and higher.
      If you are using an older version of Python, see
      `these instructions to install pip <http://stackoverflow.com/a/12476379>`_

    * Most packages will be downloaded as Wheel (``*.whl``) files. To install
      these files, type::

          pip install C:\Path\to\downloaded\file\package-file-name.whl

  * Cython

    * http://www.lfd.uci.edu/~gohlke/pythonlibs/#cython

    * Download the ``*.whl`` file for your Python architecture (32-bit or 64-bit)
      and Python 2.7 (indicated by ``cp27`` in the file name).

    * Cython must be installed in the version of Python that has SCons installed

  * NumPy

    * http://www.lfd.uci.edu/~gohlke/pythonlibs/#numpy

    * Download the ``*.whl`` file for your Python architecture (32-bit or 64-bit)
      and Python 2.7 (indicated by ``cp27`` in the file name).

* In addition to the general software, building the Python 3 module also requires

  * Python 3

    * https://www.python.org/downloads/

    * Cantera supports Python 3.3 and higher

    * Be sure to choose the appropriate architecture for your system - either
      32-bit or 64-bit

    * Be careful that the installation of Python 3 does not come before Python 2
      on your ``PATH`` environment variable

  * Pip

    * Pip should be distributed with Python version 3.4 and higher.
      If you are using an older version of Python, see
      `these instructions to install pip <http://stackoverflow.com/a/12476379>`_

    * Most packages will be downloaded as Wheel (``*.whl``) files. To install
      these files, type::

          pip3 install C:\Path\to\downloaded\file\package-file-name.whl

  * Cython

    * http://www.lfd.uci.edu/~gohlke/pythonlibs/#cython

    * Download the ``*.whl`` file for your Python architecture (32-bit or 64-bit)
      and Python 2.7 (indicated by ``cp27`` in the file name).

    * Cython must be installed in the version of Python that has SCons installed

  * NumPy

    * http://www.lfd.uci.edu/~gohlke/pythonlibs/#numpy

    * Download the ``*.whl`` file for your Python architecture (32-bit or 64-bit)
      and Python 3.x (indicated by ``cp3x`` in the file name, where x matches
      your version of Python).

* In addition to the general software, building the MATLAB toolbox also requires:

  * MATLAB version later than 2009a

    * Typically installed to::

        C:\Program Files\MATLAB\R20YYn

      where ``YY`` is a two digit year and ``n`` is either ``a`` or ``b``/Applications/MATLAB_R2011a.app

.. _sec-macos:

OS X & macOS
------------

General Notes
^^^^^^^^^^^^^

* It is not recommended to use the system-installed version of Python to build
  Cantera. Instead, the following instructions use Homebrew to install a
  separate copy of Python, independent from the system Python.

* To download the source code, installing ``git`` via HomeBrew is highly recommended.

* SCons is only available for Python 2, so building the Python 3 module requires
  two installations of Python (one of Python 2 and one of Python 3), even if you
  do not intend to build the Python 2 module.

* Cython is only required to be installed for the version of Python that also
  has SCons installed; following the instructions below will install Cython for
  the version of Python 2 installed in the system directories. The minimum
  compatible Cython version is 0.23.

.. _sec-mac-os-reqs:

OS X & macOS Requirements
^^^^^^^^^^^^^^^^^^^^^^^^^

* OS X 10.9 (Mavericks) or newer required; 10.10 (Yosemite) or newer is recommended

* To build any of the Cantera modules, you will need to install

  * Xcode

    * Download and install from the App Store

    * From a Terminal, run::

        sudo xcode-select --install

      and agree to the Xcode license agreement

  * Homebrew

    * http://brew.sh

    * From a Terminal, run::

        /usr/bin/ruby -e "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/master/install)"

  * Once Homebrew is installed, the rest of the dependencies can be installed with::

      brew install python scons boost

* In addition to the general software, building the Python 2 module also requires::

    pip install cython numpy

* In addition to the general software, building the Python 3 module also requires::

    brew install python3
    pip install cython
    pip3 install numpy

  Note that Cython should be installed into the version of Python that has SCons
  installed.

* In addition to the general software, building the Fortran module also requires::

    brew install gcc

* In addition to the general software, building the MATLAB toolbox also requires:

  * MATLAB version later than 2009a

    * Typically installed to::

        /Applications/MATLAB_R20YYn.app

      where ``YY`` is a two digit year and ``n`` is either ``a`` or ``b``
