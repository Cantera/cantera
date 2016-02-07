.. _sec-install:

******************
Installing Cantera
******************

.. contents::
   :local:
   :depth: 2

.. _sec-install-win:

Windows
=======

Windows installers are provided for stable versions of Cantera. These
installation instructions are for Cantera 2.2.0.

1. **Choose your Python version and architecture**

   - On Windows, Cantera supports Python 2.7, Python 3.3, Python 3.4, and Python
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

   - Cantera can also be used with alternative Python distributions such as
     `Anaconda <https://store.continuum.io/cshop/anaconda/>`_ or the Enthought
     `Canopy <https://www.enthought.com/products/canopy/>`_ distribution. These
     distributions will generally be based on the 64-bit version of Python 2.7,
     and will include Numpy as well as many other packages useful for scientific
     users.

3. **Install Numpy and optional Python packages**

   - Go to the `Unofficial Windows Binaries for Python Extension Packages page
     <http://www.lfd.uci.edu/~gohlke/pythonlibs/#numpy>`_.

   - Download the most recent release (distributed as a "wheel" archive) of the
     1.x series for Python *X.Y* that matches your Python architecture. The
     binaries for Cantera 2.2.0 require Numpy 1.8.0 or newer, e.g. In the
     filename, the digits after "cp" indicate the Python version, e.g.
     ``numpy‑1.8.2+mkl‑cp35‑none‑win_amd64.whl`` is the installer for 64-bit
     Python 3.5.

   - From an administrative command prompt, install the downloaded wheel using
     pip, e.g.::

         c:\python35\scripts\pip.exe install "%USERPROFILE%\Downloads\numpy‑1.8.2+mkl‑cp35‑none‑win_amd64.whl"

   - If you plan on using Cantera from Python, you may also want to install
     IPython (an advanced interactive Python interpreter) and Matplotlib (a
     plotting library), which are also available from the above link (note that
     you may also need to download additional dependencies for each of these
     packages). Matplotlib is required to run some of the Python examples.

4. **Remove old versions of Cantera**

   - Use The Windows "Add/Remove Programs" interface

   - Remove both the main Cantera package and the Python module.

   - The Python module will be listed as "Python *X.Y* Cantera ..."

5. **Install Cantera**

   - Go to the `Cantera Releases <https://github.com/Cantera/cantera/releases>`_
     page.

     - *64-bit*: Download **Cantera-2.2.0-x64.msi** and
       **Cantera-Python-2.2.0-x64-pyX.Y.msi**.
     - *32-bit*: Download **Cantera-2.2.0-x86.msi** and
       **Cantera-Python-2.2.0-x86-pyX.Y.msi**.

   - If you are only using the Python module, you do not need to download and
     install the base package.

   - Run the installer(s).

6. **Configure Matlab** (optional)

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

7. **Test the installation**

   - Python::

         import cantera
         gas = cantera.Solution('gri30.cti')
         h2o = cantera.PureFluid('liquidvapor.cti', 'water')

   - Matlab::

         gas = IdealGasMix('gri30.cti')
         h2o = importPhase('liquidvapor.cti','water')

.. _sec-install-osx:

Mac OS X
========

Cantera can be installed on OS X using either Homebrew or MacPorts. With
Homebrew, the current stable, maintenance, or development versions of Cantera
can be installed, and both the Python 2.7 and Python 3.x modules are available,
as well as the Matlab toolbox. The MacPorts portfile supports the current stable
version of Cantera and builds the Python 2.7 module.

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
Package Archive (PPA). As of Cantera 2.2.1, packages are available for Ubuntu
Ubuntu 14.04 LTS (Trusty Tahr), Ubuntu 15.04 (Vivid Vervet), and Ubuntu 15.10
(Wily Werewolf). To see which Ubuntu releases and Cantera versions are currently
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
