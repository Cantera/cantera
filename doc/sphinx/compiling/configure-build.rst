
.. _sec-determine-config:

Determine configuration options
===============================

* Run ``scons help`` to see a list all configuration options for Cantera, or
  see :ref:`scons-config`.

* Configuration options are specified as additional arguments to the ``scons``
  command, e.g.::

    scons command option=value

  where ``scons`` is the program that manages the build steps, and ``command``
  is most commonly one of

    * ``build``
    * ``test``
    * ``clean``

  Other commands are possible, and are explained in :ref:`sec-build-commands`.

* SCons saves configuration options specified on the command line in the file
  ``cantera.conf`` in the root directory of the source tree, so generally it is
  not necessary to respecify configuration options when rebuilding Cantera. To
  unset a previously set configuration option, either remove the corresponding
  line from ``cantera.conf`` or use the syntax::

    scons command option_name=

* Sometimes, changes in your environment can cause SCons's configuration tests
  (e.g., checking for libraries or compiler capabilities) to unexpectedly fail.
  To force SCons to re-run these tests rather than trusting the cached results,
  run scons with the option ``--config=force``.

* The following lists of options are not complete, they show only some commonly
  used options. The entire list of options can be found in :ref:`scons-config`.

Common Options
^^^^^^^^^^^^^^^

* :ref:`blas_lapack_libs <blas-lapack-libs>`

  * On OS X, the Accelerate framework is automatically used to provide
    optimized versions of BLAS and LAPACK, so the ``blas_lapack_libs``
    option should generally be left unspecified.

* :ref:`blas_lapack_dir <blas-lapack-dir>`
* :ref:`boost_inc_dir <boost-inc-dir>`
* :ref:`debug <debug>`
* :ref:`optimize <optimize>`
* :ref:`prefix <prefix>`
* :ref:`sundials_include <sundials-include>`
* :ref:`sundials_libdir <sundials-libdir>`

Python 2 Module Options
^^^^^^^^^^^^^^^^^^^^^^^

By default, SCons will attempt to build the Cython-based Python module for
Python 2, if both Numpy and Cython are installed. The following options control
how the Python 2 module is built:

* :ref:`python_cmd <python-cmd>`
* :ref:`python_package <python-package>`
* :ref:`python_prefix <python-prefix>`

Python 3 Module Options
^^^^^^^^^^^^^^^^^^^^^^^

If SCons detects a Python 3 interpreter installed in a default location
(i.e., ``python3`` is on the ``PATH`` environment variable) or
``python3_package`` is ``y``, SCons will try to build the Python module
for Python 3. The following SCons options control how the Python 3 module is
built:

* :ref:`python3_cmd <python3-cmd>`
* :ref:`python3_package <python3-package>`
* :ref:`python3_prefix <python3-prefix>`

Note that even when building the Python 3 Cantera module, you should still use
Python 2 with SCons, as SCons does not currently support Python 3.

Windows Only Options
^^^^^^^^^^^^^^^^^^^^

.. note::

    The ``cantera.conf`` file uses the backslash character ``\`` as an escape
    character. When modifying this file, backslashes in paths need to be escaped
    like this: ``boost_inc_dir = 'C:\\Program Files (x86)\\boost\\include'``
    This does not apply to paths specified on the command line. Alternatively,
    you can use forward slashes (``/``) in paths.

* In Windows there aren't any proper default locations for many of the packages
  that Cantera depends on, so you will need to specify these paths explicitly.

* Remember to put double quotes around any paths with spaces in them, e.g.
  ``"C:\Program Files"``.

* By default, SCons attempts to use the same architecture as the copy of Python
  that is running SCons, and the most recent installed version of the Visual
  Studio compiler. If you aren't building the Python module, you can override
  this with the configuration options ``target_arch`` and ``msvc_version``.

* To compile with MinGW, specify the :ref:`toolchain <toolchain>` option::

    toolchain=mingw

* :ref:`msvc_version <msvc-version>`
* :ref:`target_arch <target-arch>`
* :ref:`toolchain <toolchain>`

MATLAB Toolbox Options
^^^^^^^^^^^^^^^^^^^^^^

Building the MATLAB toolbox requires an installed copy of MATLAB, and the path
to the directory where MATLAB is installed must be specified using the following
option:

* :ref:`matlab_path <matlab-path>`

Fortran Module Options
^^^^^^^^^^^^^^^^^^^^^^

Building the Fortran module requires a compatible Fortran comiler. SCons will
attempt to find a compatible compiler by default in the ``PATH`` environment
variable. The following options control how the Fortran module is built:

* :ref:`f90_interface <f90-interface>`
* :ref:`FORTRAN <FORTRAN>`

Documentation Options
^^^^^^^^^^^^^^^^^^^^^

The following options control if the documentation is built:

* :ref:`doxygen_docs <doxygen-docs>`
* :ref:`sphinx_docs <sphinx-docs>`

Less Common Options
^^^^^^^^^^^^^^^^^^^

* :ref:`CC <CC>`
* :ref:`CXX <CXX>`
* :ref:`env_vars <env-vars>`
* :ref:`layout <layout>`
* :ref:`VERBOSE <VERBOSE>`

.. _sec-build-commands:

Build Commands
==============

The following options are possible as commands to SCons, i.e., the first
argument after ``scons``::

    scons command

* ``scons help``
    Print a description of user-specifiable options.

* ``scons build``
    Compile Cantera and the language interfaces using
    default options.

* ``scons clean``
    Delete files created while building Cantera.

* ``[sudo] scons install``
    Install Cantera.

* ``[sudo] scons uninstall``
    Uninstall Cantera.

* ``scons test``
    Run all tests which did not previously pass or for which the
    results may have changed.

* ``scons test-reset``
    Reset the passing status of all tests.

* ``scons test-clean``
    Delete files created while running the tests.

* ``scons test-help``
    List available tests.

* ``scons test-NAME``
    Run the test named "NAME".

* ``scons <command> dump``
    Dump the state of the SCons environment to the
    screen instead of doing ``<command>``, e.g.
    ``scons build dump``. For debugging purposes.

* ``scons samples``
    Compile the C++ and Fortran samples.

* ``scons msi``
    Build a Windows installer (.msi) for Cantera.

* ``scons sphinx``
    Build the Sphinx documentation

* ``scons doxygen``
    Build the Doxygen documentation

Compile Cantera & Test
======================

* Run SCons with the list of desired configuration options::

    scons build ...

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

    * Messages printed to the console while running ``scons test``
    * Output files generated by the tests

Building Documentation
^^^^^^^^^^^^^^^^^^^^^^

* To build the Cantera HTML documentation, run the commands::

    scons doxygen
    scons sphinx

  or append the options ``sphinx_docs=y`` and ``doxygen_docs=y`` to the build
  command, e.g.::

    scons build doxygen_docs=y sphinx_docs=y
