.. _sec-install:

******************
Installing Cantera
******************

.. contents::
   :local:
   :depth: 2

.. _sec-install-conda:

Conda
=====

`Anaconda <https://www.continuum.io/downloads>`_ and `Miniconda
<http://conda.pydata.org/miniconda.html>`_ are Python distributions for which
Cantera is available through the `conda` package manager. Both distributions are
available for Linux, OS X, and Windows. The base Anaconda distribution includes
a large number of Python packages that are widely used in scientific
applications. Miniconda is a minimal distribution, where all of the packages
available in Anaconda can be installed using the package manager. Note that
installing Cantera using conda will only provide the Cantera Python module. If
you want to use the other Cantera interfaces, see the OS-specific installation
options below.

For more details on how to use conda, see the `conda documentation
<http://conda.pydata.org/docs/intro.html>`_.

**Option 1: Create a new environment for Cantera**

If you have just installed Anaconda or Miniconda, the following instructions
will create a conda environment where you can use Cantera. For this example, the
environment is named ``spam``. From the command line, run::

    conda create -n spam -c cantera cantera ipython matplotlib

This will create an environment with Cantera, IPython, Matplotlib, and all their
dependencies installed. Although conda can install a large set of packages by
default, it is also possible to install packages such as Cantera that are
maintained independently. These additional channels from which packages may be
obtained are specified by adding the ``-c`` option in the ``install`` or
``create`` commands. In this case, we want to install Cantera from the
``cantera`` channel, so we add ``-c cantera`` and to tell conda to look at the
``cantera`` channel in addition to the default channels.

If you are running Linux or OS X, you can then activate this environment by
running::

    source activate spam

If you are running Windows, the equivalent command is::

    activate spam

**Option 2: Install Cantera in an existing environment**

First, activate your environment (assumed to be named ``baked_beans``; if you've
forgotten the name of the conda environment you wanted to use, the command
``conda env list`` can help). For Linux and OS X, this is done by running::

    source activate baked_beans

For Windows users, the command is::

    activate baked_beans

Then, install Cantera by running::

    conda install -c cantera cantera

**Option 3: Install the development version of Cantera**

To install a recent development snapshot (i.e. an alpha or beta version) of
Cantera in an existing environment, run::

    conda install -c cantera/label/dev cantera

If you later want to revert back to the stable version, first remove and then
reinstall Cantera::

    conda remove cantera
    conda install -c cantera cantera

.. _sec-install-win:

Windows
=======

Windows installers are provided for stable versions of Cantera. These
installation instructions are for Cantera 2.3.0. Use these installers if you
want to work with a copy of Python downloaded from `Python.org
<https://www.python.org/>`_. If you are using Anaconda / Miniconda, see the
directions :ref:`above <sec-install-conda>`.

1. **Choose your Python version and architecture**

   - On Windows, Installers are provided for Python 2.7, Python 3.4, and Python
     3.5. Python 3.5 is recommended unless you need to use legacy code that does
     not work with Python 3. You can install multiple Cantera Python modules
     simultaneously.

   - Cantera supports both 32- and 64- bit Python installations.

   - You need choose the matching Cantera installer for your Python version and
     machine architecture.

   - The rest of these instructions will refer to your chosen version of Python
     as *X.Y*.

   - If you are using Matlab, you must use the same architecture for Cantera and
     Matlab. Matlab defaults to 64-bit if you are running a 64-bit operating
     system.

2. **Install Python**

   - Go to `python.org <https://www.python.org/>`_.

     - *64-bit*: Download the most recent "Windows X86-64 MSI Installer" for
       Python *X.Y*.
     - *32-bit*: Download the most recent "Windows x86 MSI Installer" for
       Python *X.Y*.

   - Run the installer. The default installation options should be fine.

   - Python is required in order to work with `.cti` input files even if you are
     not using the Python interface to Cantera.

   - Cantera can also be used with alternative Python distributions such as the
     Enthought `Canopy <https://www.enthought.com/products/canopy/>`_
     distribution. These distributions will generally be based on the 64-bit
     version of Python 2.7, and will include Numpy as well as many other
     packages useful for scientific users.

3. **Install the Visual C++ Redistributable for Visual Studio 2015**

   - If you are using Python 3.5, you can skip this step as this will have
     already been installed when you installed Python.

   - Go to the `Microsoft Visual C++ Redistributable Download Page
     <https://www.microsoft.com/en-us/download/details.aspx?id=48145>`_.

     - *64-bit*: Download ``vc_redist.x64.exe``
     - *32-bit*: Download ``vc_redist.x86.exe``

   - Run the installer.

   - If this package is not installed, you will encounter the following error
     when importing the `cantera` module::

         ImportError: DLL load failed: The specified module could not be found.

4. **Install Numpy and optional Python packages**

   - Go to the `Unofficial Windows Binaries for Python Extension Packages page
     <http://www.lfd.uci.edu/~gohlke/pythonlibs/#numpy>`_.

   - Download the most recent release (distributed as a "wheel" archive) of the
     1.x series for Python *X.Y* that matches your Python architecture. In the
     filename, the digits after "cp" indicate the Python version, e.g.
     ``numpy‑1.11.2+mkl‑cp35‑none‑win_amd64.whl`` is the installer for 64-bit
     Python 3.5. The Windows installers for Cantera 2.3.0 require Numpy 1.10 or
     newer.

   - From an administrative command prompt, install the downloaded wheel using
     pip, e.g.::

         c:\python35\scripts\pip.exe install "%USERPROFILE%\Downloads\numpy‑1.11.2+mkl‑cp35‑none‑win_amd64.whl"

   - If you plan on using Cantera from Python, you may also want to install
     IPython (an advanced interactive Python interpreter) and Matplotlib (a
     plotting library), which are also available from the above link (note that
     you may also need to download additional dependencies for each of these
     packages). Matplotlib is required to run some of the Python examples.

5. **Remove old versions of Cantera**

   - Use The Windows "Add/Remove Programs" interface

   - Remove both the main Cantera package and the Python module.

   - The Python module will be listed as "Python *X.Y* Cantera ..."

6. **Install Cantera**

   - Go to the `Cantera Releases <https://github.com/Cantera/cantera/releases>`_
     page.

     - *64-bit*: Download **Cantera-2.3.0-x64.msi** and
       **Cantera-Python-2.3.0-x64-pyX.Y.msi**.
     - *32-bit*: Download **Cantera-2.3.0-x86.msi** and
       **Cantera-Python-2.3.0-x86-pyX.Y.msi**.

   - If you are only using the Python module, you do not need to download and
     install the base package.

   - Run the installer(s).

7. **Configure Matlab** (optional)

   - Set the environment variable ``PYTHON_CMD``

     - From the *Start* menu (Windows 7) or the *Start* screen (Windows 8) type
       "edit environment" and select "Edit environment variables for your
       account".
     - Add a *New* variable with ``PYTHON_CMD`` as the *name* and the full path
       to the Python executable (e.g. ``C:\python35\python.exe``) as the
       *value*.
     - Setting ``PYTHON_CMD`` is not necessary if the path to ``python.exe`` is
       in your ``PATH`` (which can be set from the same configuration dialog).

   - Launch Matlab

   - Go to *File->Set Path...*

   - Select *Add with Subfolders*

   - Browse to the folder ``C:\Program Files\Cantera\matlab\toolbox``

   - Select *Save*, then *Close*.

8. **Test the installation**

   - Python::

         import cantera
         gas = cantera.Solution('gri30.cti')
         h2o = cantera.PureFluid('liquidvapor.cti', 'water')

   - Matlab::

         gas = IdealGasMix('gri30.cti')
         h2o = Solution('liquidvapor.cti','water')

.. _sec-install-osx:

Mac OS X
========

Cantera can be installed on OS X using either Homebrew, MacPorts, or Anaconda /
Miniconda. If you are using Anaconda / Miniconda, see the directions
:ref:`above <sec-install-conda>`. With Homebrew, the current stable, or
development version of Cantera can be installed, and both the Python 2.7 and
Python 3.x modules are available, as well as the Matlab toolbox. The MacPorts
portfile supports the current stable version of Cantera and builds the Python
2.7 module.

Homebrew
---------
These instructions have been tested on Mac OS X 10.9 (Mavericks) with Xcode 5.1
and Mac OS X 10.10 (Yosemite) with Xcode 6.1. If you've used Homebrew before,
you can skip any steps which have already been completed.

1. **Install Xcode and Homebrew**

   - Install Xcode from the App Store

   - From a Terminal, run::

         sudo xcode-select --install
         sudo xcodebuild -license

     and agree to the Xcode license agreement.

   - Install `Homebrew <http://brew.sh/>`_ by running the following command in a
     Terminal::

         ruby -e "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/master/install)"

2. **Set up the compilation environment**

   - Run the following commands::

         brew tap homebrew/science
         brew update
         brew install python scons

   - Verify that your path is set up to use Homebrew's version of Python by
     running::

         which python

     If this command does not print ``/usr/local/bin/python``, add the following
     to ``~/.bash_profile`` (creating this file if it doesn't already exist; you
     can use the command line editor ``nano`` to edit this file)::

         export PATH=/usr/local/bin:$PATH

     and then run::

         source ~/.bash_profile

   - Install Python packages required to compile Cantera by running::

         pip install cython numpy

     Note that these packages are required even if you do not plan on using the
     Cantera Python 2 module.

   - If you plan on using Cantera from Python, you may also want to install
     IPython (an advanced interactive Python interpreter) and Matplotlib (a
     plotting library). Matplotlib is required to run some of the Python
     examples::

         pip install ipython matplotlib

   - If you want to build the Cantera Python 3 module, run::

         brew install python3
         pip3 install numpy cython

     and, optionally::

         pip3 install ipython matplotlib

3. **Compile and install Cantera**

   * To compile and install Cantera using the default configuration, run::

         brew install cantera

   * The following options are supported:

     ``--HEAD``
         Installs the current development version of Cantera.

     ``--with-python3``
         Install the Python 3 module.

     ``--with-matlab=/Applications/MATLAB_R2014a.app/``
         Installs the Matlab toolbox (with the path modified to match your
         installed Matlab version)

     ``--without-sundials``
         Do not use an external SUNDIALS version to build Cantera. Users
         choosing this option will not be able to run sensitivity analysis
         of Reactor Networks, but it may prevent errors when installing
         the Matlab toolbox.

     ``--without-check``
         NOT RECOMMENDED! Disable automatic testing of Cantera during the
         installation process.

   * These options are specified as additional arguments to the ``brew install``
     command, e.g.::

         brew install cantera --HEAD --with-python3

   * If you are installing the Matlab toolbox, the recommended command is::

         brew install cantera --with-matlab=/Applications/MATLAB_R2014a.app/ --without-sundials

   * If something goes wrong with the Homebrew install, re-run the command with
     the ``-v`` flag to get more verbose output that may help identify the
     source of the problem::

         brew install -v cantera

   * If Homebrew claims that it can't find a formula named ``cantera``, you may
     be able to fix it by running the commands::

         brew doctor
         brew tap --repair

4. **Test Cantera Installation (Python)**

   * The Python examples will be installed in::

         /usr/local/lib/pythonX.Y/site-packages/cantera/examples/

     where ``X.Y`` is your Python version, e.g. ``2.7``.

   * You may find it convenient to copy the examples to your Desktop::

         cp -r /usr/local/lib/python2.7/site-packages/cantera/examples ~/Desktop/cantera_examples

   * To run an example::

         cd cantera_examples/reactors
         python reactor1.py

5. **Test Cantera Installation (Matlab)**

   * The Matlab toolbox, if enabled, will be installed in::

         /usr/local/lib/cantera/matlab

   * To use the Cantera Matlab toolbox, run the following commands in Matlab
     (each time you start Matlab), or add them to a ``startup.m`` file located
     in ``/Users/$USER/Documents/MATLAB``, where ``$USER`` is your username::

         addpath(genpath('/usr/local/lib/cantera/matlab'))
         setenv('PYTHON_CMD', '/usr/local/bin/python')

   * The Matlab examples will be installed in::

         /usr/local/share/cantera/samples/matlab

   * You may find it convenient to copy the examples to your user directory::

         cp -r /usr/local/share/cantera/samples/matlab ~/Documents/MATLAB/cantera_examples

MacPorts
--------

If you have MacPorts installed (see https://www.macports.org/install.php), you
can install Cantera by executing::

    sudo port install cantera

from the command line. All dependencies will be installed automatically.

MacPorts installs its own Python interpreter. Be sure to be actually using it by
checking::

    sudo port select python python27

.. _sec-install-ubuntu:

Ubuntu
======

Ubuntu packages are provided for recent versions of Ubuntu using a Personal
Package Archive (PPA). As of Cantera 2.3.0, packages are available for Ubuntu
Ubuntu 14.04 LTS (Trusty Tahr), Ubuntu 16.04 (Xenial Xerus), and Ubuntu 16.10
(Yakkety Yak). To see which Ubuntu releases and Cantera versions are currently
available, visit https://launchpad.net/~speth/+archive/ubuntu/cantera

The available packages are:

- ``cantera-python`` - The Cantera Python module for Python 2.

- ``cantera-python3`` - The Cantera Python module for Python 3.

- ``cantera-dev`` - Libraries and header files for compiling your own C++ and
  Fortran 90 programs that use Cantera.

To add the Cantera PPA::

    sudo aptitude install python-software-properties
    sudo apt-add-repository ppa:speth/cantera
    sudo aptitude update

To install all of the Cantera packages::

    sudo aptitude install cantera-python cantera-python3 cantera-dev

or install whichever subset you need by adjusting the above command.

If you plan on using Cantera from Python, you may also want to install IPython
(an advanced interactive Python interpreter) and Matplotlib (a plotting
library), which are also available from the above link. Matplotlib is required
to run some of the Python examples. For Python 2, these packages can be
installed with::

    pip2 install ipython matplotlib

And for Python 3, these packages can be installed with::

    pip3 install ipython matplotlib

You may need to install ``pip`` first; instructions can be found on the
`pip installation instructions.
<https://pip.pypa.io/en/latest/installing.html#install-pip>`_
You may need to have superuser access to install packages into the system
directories. Alternatively, you can add ``--user`` after ``pip install`` but
before the package names to install into your local user directory. An
alternative method is to use the Ubuntu repositories, but these tend to
be very out of date. For Python 2, the command is::

    sudo aptitude install ipython python-matplotlib

And for Python 3, these packages can be installed with::

    sudo aptitude install ipython3 python3-matplotlib
