.. _sec-install:

******************
Installing Cantera
******************

.. contents::
   :local:
   :depth: 1

.. _sec-install-win:

Windows
=======

Windows installers are provided for stable versions of Cantera. These
installation instructions are for Cantera 2.1.1.

1. **Choose your Python version and architecture**

   - On Windows, Cantera supports Python 2.7 and Python 3.3. Python 3.3 is
     recommended unless you need to use legacy code that does not work with
     Python 3. You can install both Cantera Python modules simultaneously.

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
       Python *X.Y* (i.e. prefer 3.3.5 to 3.3.4, but not 3.4.1).
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

3. **Install Numpy**

   - Go to the `Unofficial Windows Binaries for Python Extension Packages page
     <http://www.lfd.uci.edu/~gohlke/pythonlibs/#numpy>`_.

   - Download the most recent release of the 1.x series for Python *X.Y* that
     matches your Python architecture. The binaries for Cantera 2.1.1 require
     Numpy 1.8.0 or newer.

   - Run the installer.

4. **Remove old versions of Cantera**

   - Use The Windows "Add/Remove Programs" interface

   - Remove both the main Cantera package and the Python module.

   - The Python module will be listed as "Python *X.Y* Cantera ..."

5. **Install Cantera**

   - Go to the `Cantera Downloads
     <https://sourceforge.net/projects/cantera/files/cantera/2.1.1/>`_ page.

     - *64-bit*: Download **Cantera-2.1.1-x64.msi** and
       **Cantera-Python-2.1.1-x64-pyX.Y.msi**.
     - *32-bit*: Download **Cantera-2.1.1-x86.msi** and
       **Cantera-Python-2.1.1-x86-pyX.Y.msi**.

   - If you are only using the Python module, you do not need to download and
     install the base package.

   - Run the installer(s).

6. **Configure Matlab** (optional)

   - Set the environment variable ``PYTHON_CMD``

     - From the *Start* menu (Windows 7) or the *Start* screen (Windows 8) type
       "edit environment" and select "Edit environment variables for your
       account".
     - Add a *New* variable with ``PYTHON_CMD`` as the *name* and the full path
       to the Python executable (e.g. ``C:\python27\python.exe``) as the
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

The easiest way to install Cantera on OS X is by using Homebrew. These
instructions have been tested on Mac OS X 10.9 (Mavericks) with Xcode 5.1.

Prerequisites
-------------

If you've used Homebrew before, you may have already completed some of these
steps and can skip them.

- Install Xcode from the App store

- From a Terminal, run::

      sudo xcodebuild

  and agree to the Xcode license agreement

- Install `Homebrew <http://brew.sh/>`_

- Run the following commands:

      brew tap homebrew/science
      brew update
      brew install python scons sundials

- Put ``/usr/local/bin`` at the front of your path, e.g. add the following to
  ``~/.bash_profile`` (creating this file if it doesn't already exist)::

      export PATH=/usr/local/bin:$PATH

- Run::

      source ~/.bash_profile

- If you want to build the Cantera Python 2 module, run::

      pip install cython numpy

- If you want to build the Cantera Python 3 module, run::

      brew install python3
      pip3 install numpy cython

Installing Cantera
------------------

The installation command for Cantera supports several options:

- To install Cantera with additional patches that will be included in the next
  maintenance release, use the flag: ``--devel``

- To Install the current development version of Cantera, use the flag:
  ``--HEAD``

- To install the Matlab toolbox, use the flag
  ``--with-matlab=/Applications/MATLAB_R2014a.app/`` (with the version modified
  to match your installed Matlab version)

Install Cantera by adding the desired options to the ``brew install`` command,
e.g.::

    brew install cantera --devel --with-matlab=/Applications/MATLAB_R2014a.app/
