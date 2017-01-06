
.. _sec-special-compiling-cases:

***********************
Special Compiling Cases
***********************

This guide explains some of the less common ways to build Cantera

.. contents::
   :local:

.. _sec-intel-compilers:

Intel Compilers
===============

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

* When installing Cantera after building with the Intel compiler, the normal
  method of using ``sudo`` to install Cantera to the system default directories
  will not work because ``sudo`` does not pass the environment variables needed
  by the Intel compiler. Instead, you will need to do something like::

    scons build ...
    sudo -s
    source /path/to/compilervars.sh intel64
    scons install
    exit

  Another option is to set the :ref:`prefix <prefix>` option to a directory
  for which you have write permissions, and specify the ``USER`` value to the
  :ref:`python_prefix <python-prefix>` or :ref:`python3_prefix <python3-prefix>`
  option.
