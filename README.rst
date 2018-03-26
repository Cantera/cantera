
.. Cantera

|cantera|

|doi| |codecov| |travisci| |appveyor| |release|


What is Cantera?
================

Cantera is an open-source collection of object-oriented software tools for
problems involving chemical kinetics, thermodynamics, and transport
processes. Among other things, it can be used to:

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
using the following Binder link:

.. image:: http://mybinder.org/badge.svg
    :target: http://mybinder.org:/repo/cantera/cantera-jupyter

Installation
============

`Installation instructions for the current release of Cantera
<http://cantera.github.io/docs/sphinx/html/install.html>`_ are available from
the main `Cantera documentation site
<http://cantera.github.io/docs/sphinx/html/index.html>`_. Installers are
provided for Windows (MSI packages), Mac OS X (through Homebrew), and
Ubuntu. Anaconda packages containing the Cantera Python module are also
available for Windows, OS X, and Linux.

.. image:: https://anaconda.org/cantera/cantera/badges/installer/conda.svg
    :target: https://anaconda.org/Cantera/cantera

For other platforms, or for users wishing to install a development version of
Cantera, `compilation instructions
<http://cantera.github.io/docs/sphinx/html/compiling/index.html>`_ are also available.

Documentation
=============

The `documentation <http://cantera.github.io/docs/sphinx/html/index.html>`_
offers a number of starting points:

- `Python tutorial
  <http://cantera.github.io/docs/sphinx/html/cython/tutorial.html>`_
- `Application Examples in Python
  <https://github.com/Cantera/cantera-jupyter#cantera-jupyter>`_
- `A guide to Cantera's input file format
  <http://cantera.github.io/docs/sphinx/html/cti/index.html>`_
- `A list of frequently asked questions
  <http://cantera.github.io/docs/sphinx/html/faq.html>`_

`Documentation for the development version of Cantera
<http://cantera.github.com/dev-docs/sphinx/html/index.html>`_ is also available.

Code of Conduct
===============

.. image:: https://img.shields.io/badge/code%20of%20conduct-contributor%20covenant-green.svg?style=flat-square
    :alt: conduct
    :target: http://contributor-covenant.org/version/1/4/

In order to have a more open and welcoming community, Cantera adheres to a
`code of conduct <CODE_OF_CONDUCT.md>`_ adapted from the `Contributor Covenent
code of conduct <http://contributor-covenant.org/>`_.

Please adhere to this code of conduct in any interactions you have in the
Cantera community. It is strictly enforced on all official Cantera
repositories, websites, users' group, and other resources.
If you encounter someone violating these terms, please
`contact the code of conduct team <mailto:conduct@cantera.org>`_
(`@speth <https://github.com/speth>`_,
`@bryanwweber <https://github.com/bryanwweber>`_, and
`@kyleniemeyer <https://github.com/kyleniemeyer>`_)
and we will address it as soon as possible.

Development Site
================

The current development version is 2.4.0b1. The `latest Cantera source code <https://github.com/Cantera/cantera>`_, the
`issue tracker <https://github.com/Cantera/cantera/issues>`_ for bugs and
enhancement requests, `downloads of Cantera releases and binary installers
<https://github.com/Cantera/cantera/releases>`_ , and the `Cantera wiki
<https://github.com/Cantera/cantera/wiki>`_ can all be found on Github.

Users' Group
============

The `Cantera Users' Group <http://groups.google.com/group/cantera-users>`_ is a
message board / mailing list for discussions amongst Cantera users.

Cantera Gitter Chat
============

.. image:: https://badges.gitter.im/org.png
   :target: https://gitter.im/Cantera/Lobby


The `Cantera Gitter Chat <https://gitter.im/Cantera/Lobby>`_ is a public chat client that is linked to users' Github account. The developers do not closely monitor the discussion, so *any* discussion at all of Cantera functionality such as how to use certain function calls, syntax problems, input files, etc. should be directed the User's Group. All conversations in the Gitter room will be covered under the Cantera Code of Conduct, so please be nice.

The chat room is a place to strengthen and develop the Cantera community, discuss tangentially-related topics such as how to model the underlying physics of a problem , share cool applications you’ve developed, etc. 

Summary:

“How do I perform this Cantera function call?” --> User's Group

"What do I do with the variables that a Cantera function call returns?” --> Chat


Continuous Integration Status
=============================

==============  ============  ===================
Platform        Site          Status
==============  ============  ===================
Linux & OS X    Travis CI     |travisci|
Windows x64     Appveyor      |appveyor|
==============  ============  ===================


NumFOCUS
========

Cantera is a fiscally-sponsored project of NumFOCUS, a non-profit dedicated to
supporting the open source scientific computing community.

.. image:: https://img.shields.io/badge/powered%20by-NumFOCUS-orange.svg?style=flat&colorA=E1523D&colorB=007D8A
    :target: https://numfocus.org
    :alt: Powered by NumFOCUS







.. |cantera| image:: http://cantera.github.io/docs/sphinx/html/_static/cantera-logo.png
    :target: http://cantera.github.io/docs/sphinx/html/index.html#
    :alt: cantera logo
    :width: 675px
    :align: middle
   
.. |travisci| image:: https://travis-ci.org/Cantera/cantera.svg?branch=master
    :target: https://travis-ci.org/Cantera/cantera

.. |appveyor| image:: https://ci.appveyor.com/api/projects/status/auhd35qn9cdmkpoj?svg=true
    :target: https://ci.appveyor.com/project/Cantera/cantera

.. |doi| image:: https://zenodo.org/badge/DOI/10.5281/zenodo.170284.svg
   :target: https://doi.org/10.5281/zenodo.1174508
    
.. |codecov| image:: https://img.shields.io/codecov/c/github/Cantera/cantera/master.svg
   :target: https://codecov.io/gh/Cantera/cantera?branch=master

.. |release| image:: https://img.shields.io/github/release/cantera/cantera.svg
   :target: https://github.com/Cantera/cantera/releases
   :alt: GitHub release
