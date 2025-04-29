.. Cantera

|cantera|

|doi| |codecov| |ci| |release|


What is Cantera?
================

Cantera is an open-source collection of object-oriented software tools for
problems involving chemical kinetics, thermodynamics, and transport processes.
Among other things, it can be used to:

* Evaluate thermodynamic and transport properties of mixtures
* Compute chemical equilibrium
* Evaluate species chemical production rates
* Conduct kinetics simulations with large reaction mechanisms
* Simulate one-dimensional flames
* Conduct reaction path analysis
* Create process simulations using networks of stirred reactors
* Model non-ideal fluids

Cantera can be used from Python and Matlab, or in applications written in C++
and Fortran 90. A number of `examples of Cantera's capabilities
<https://github.com/Cantera/cantera-jupyter>`_ are available in the form of
Jupyter notebooks. These examples can be tried interactively, in the cloud by
using the following MyBinder link:

.. image:: https://mybinder.org/badge.svg
    :target: https://mybinder.org/repo/cantera/cantera-jupyter

Installation
============

To ensure a smooth installation, we highly recommend setting up Cantera within a virtual environment.
Follow the steps below for a **local** installation:

0. Clone the Cantera repository::

    cd /Users/$USER/Documents
    git clone git@github.com:cerfacs/cantera.git

1. Create a virtual environment in your home directory::

    cd $HOME
    python3 -m venv env_cantera
    source ~/env_cantera/bin/activate

2. Install specific python packages::

    pip install numpy==1.26.4
    pip install scipy==1.13

3. (optional) Create an alias to load ``env_cantera``::

    echo "alias load_cantera_env='source ~/env_cantera/bin/activate'" >> ~/.zshrc
    source ~/.zshrc

4. Navigate to the Cantera directory::

    cd /Users/$USER/Documents/cantera-avbp
    
5. Switch to the right branch and run the installation script::

    git checkout dev_cantera-avbp
    python3 install_cantera.py
    

6. When prompted about NFS, answer ``no``.
7. Wait for the compilation and installation process to complete successfully.

If everything goes well, you should see the following message::

  ********************************************************************************
  To use this newly installed Cantera version, update your environment variables 
  by adding the following lines to your.bashrc (or equivalent):

  #cantera-avbp-3.0
  export PYTHONPATH=/Users/$USER/cantera-avbp/INSTALL_DIR/lib/python3.9/site-packages:$PYTHONPATH 
  export PKG_CONFIG_PATH=/Users/$USER/cantera-avbp/INSTALL_DIR/lib/pkgconfig:$PKG_CONFIG_PATH
  export LD_LIBRARY_PATH=/Users/$USER/cantera-avbp/INSTALL_DIR/lib:$LD_LIBRARY_PATH
  export PATH=/Users/$USER/cantera-avbp/INSTALL_DIR/bin:$PATH
  export PYTHON_CMD=/Users/$USER/canavbp3/bin/python3

  #Only if you dont already have a custom lib folder:
  export CUSTOM_LIB=/Users/$USER/cantera-avbp/INSTALL_DIR/mech_lib
  export LD_LIBRARY_PATH=$CUSTOM_LIB:$LD_LIBRARY_PATH

  #Required for MacOS:
  export DYLD_LIBRARY_PATH=$CUSTOM_LIB

  DONT FORGET TO SOURCE YOUR .bashrc!

  ********************************************************************************


Copy and paste the provided lines into your ``.zshrc`` or ``.bashrc`` file and **source it**.

For **KRAKEN** users, use the following commands to install Cantera:

1. Load the necessary modules::

    module purge && module load python/3.9.5
    
2. Install ``ruamel.yaml`` package (not installed by default on KRAKEN)::
  
    pip3 install --user ruamel.yaml

3. Create a virtual environment in your home folder and activate it::

    cd $HOME
    python3 -m venv --system-site-packages env_cantera
    source ~/env_cantera/bin/activate

4. (optional) Create an alias to load ``env_cantera``::

    echo "alias load_cantera_env='source ~/env_cantera/bin/activate'" >> ~/.bashrc
    source ~/.bashrc

5. Clone cantera-avbp repository (e.g. in your ``/scratch``), navigate to it and run the installation script::

    git clone https://nitrox.cerfacs.fr/cantera/cantera-avbp.git /scratch/cfd/$USER
    cd /scratch/cfd/$USER/cantera-avbp
    git checkout 3.0
    python3 install_cantera.py

When prompted about NFS, answer 'yes' and wait for the compilation and installation to finish.

If everything goes well, you will see a similar message as above, but with an additional question::


  We suggest adding the alias 'load_cantera' to purge and load necessary modules in ~/.bashrc
  Do you want to add the alias to ~/.bashrc? (yes/no)


Answer ``yes`` to this question, then copy and paste the environment variables you need into your ``.bashrc`` file and **source it**.

Finally, load Cantera with the following command::

  load_cantera


Test if everything works using a sample script, for example::

  python3 samples/python/AVBP/ARC.py


Documentation
=============

The `documentation <https://cantera.org/documentation>`_
offers a number of starting points:

- `Python tutorial
  <https://cantera.org/tutorials/python-tutorial.html>`_
- `Application Examples in Python
  <https://cantera.org/examples/jupyter/index.html>`_
- `A guide to Cantera's input file format
  <https://cantera.org/tutorials/input-files.html>`_
- `Information about the Cantera community
  <https://cantera.org/community.html>`_
- `Affiliated packages
  <https://cantera.org/affiliated-packages.html>`_

`Documentation for the development version of Cantera
<https://cantera.org/documentation/dev-docs.html>`_ is also available.

Code of Conduct
===============

.. image:: https://img.shields.io/badge/Contributor%20Covenant-2.0-4baaaa.svg
    :alt: conduct
    :target: https://www.contributor-covenant.org/version/2/0/code_of_conduct/

In order to have a more open and welcoming community, Cantera adheres to a
`code of conduct <CODE_OF_CONDUCT.md>`_ adapted from the `Contributor Covenant
code of conduct <https://contributor-covenant.org/>`_.

Please adhere to this code of conduct in any interactions you have in the
Cantera community. It is strictly enforced on all official Cantera
repositories, websites, users' group, and other resources. If you encounter
someone violating these terms, please `contact the code of conduct team
<mailto:conduct@cantera.org>`_ (`@speth <https://github.com/speth>`_,
`@bryanwweber <https://github.com/bryanwweber>`_, and `@kyleniemeyer
<https://github.com/kyleniemeyer>`_) and we will address it as soon as
possible.

Development Site
================

The current development version is 3.0.0. The current stable version is
3.0.0. The `latest Cantera source code <https://github.com/Cantera/cantera>`_,
the `issue tracker <https://github.com/Cantera/cantera/issues>`_ for bugs and
enhancement requests, `downloads of Cantera releases and binary installers
<https://github.com/Cantera/cantera/releases>`_ , and the `Cantera wiki
<https://github.com/Cantera/cantera/wiki>`_ can all be found on Github.

Users' Group
============

The `Cantera Users' Group <https://groups.google.com/group/cantera-users>`_ is a
message board/mailing list for discussions amongst Cantera users.

Continuous Integration Status
=============================

|ci|

NumFOCUS
========

Cantera is a fiscally-sponsored project of `NumFOCUS <https://numfocus.org>`__,
a non-profit dedicated to supporting the open source scientific computing
community. Please consider `making a donation
<https://numfocus.salsalabs.org/donate-to-cantera/index.html>`__ to support the
development of Cantera through NumFOCUS.

.. image:: https://img.shields.io/badge/powered%20by-NumFOCUS-orange.svg?style=flat&colorA=E1523D&colorB=007D8A
    :target: https://numfocus.salsalabs.org/donate-to-cantera/index.html
    :alt: Powered by NumFOCUS

.. |cantera| image:: https://cantera.org/assets/img/cantera-logo.png
    :target: https://cantera.org
    :alt: cantera logo
    :width: 675px
    :align: middle

.. |ci| image:: https://github.com/Cantera/cantera/workflows/CI/badge.svg
    :target: https://github.com/Cantera/cantera/actions?query=workflow%3ACI+event%3Apush

.. |doi| image:: https://zenodo.org/badge/DOI/10.5281/zenodo.8137090.svg
   :target: https://doi.org/10.5281/zenodo.8137090

.. |codecov| image:: https://img.shields.io/codecov/c/github/Cantera/cantera/main.svg
   :target: https://codecov.io/gh/Cantera/cantera?branch=main

.. |release| image:: https://img.shields.io/github/release/cantera/cantera.svg
   :target: https://github.com/Cantera/cantera/releases
   :alt: GitHub release

.. |pip| image:: https://img.shields.io/pypi/v/cantera
   :target: https://pypi.org/project/Cantera/

.. |anaconda| image:: https://img.shields.io/conda/v/cantera/cantera
   :target: https://anaconda.org/Cantera/cantera

.. |conda-forge| image:: https://img.shields.io/conda/v/conda-forge/cantera
   :target: https://anaconda.org/conda-forge/cantera
