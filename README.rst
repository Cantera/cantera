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

|pip| |anaconda| |conda-forge|

`Installation instructions for the current release of Cantera
<https://cantera.org/install/index.html>`_ are available from the main `Cantera
documentation site <https://cantera.org>`_.

- The Python module can also be installed using pip on Windows, macOS, and Linux.

- Conda packages containing the Cantera Python and Matlab modules are also available
  for Windows, macOS, and Linux.

- Additional installation packages are provided for Windows, macOS, Ubuntu, Fedora,
  RHEL, Gentoo, and FreeBSD.

- For other platforms, or for users wishing to install a development version of
  Cantera, `compilation instructions <https://cantera.org/install/compiling-install.html#sec-compiling>`_
  are also available.

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
    :target: https://www.contributor-covenant.org/version/2/0/code-of-conduct.html

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

The current development version is 2.6.0. The current stable version is
2.6.0. The `latest Cantera source code <https://github.com/Cantera/cantera>`_,
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

.. |doi| image:: https://zenodo.org/badge/DOI/10.5281/zenodo.6387882.svg
   :target: https://doi.org/10.5281/zenodo.6387882

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