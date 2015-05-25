
******************************
Compiling Cantera C++ Programs
******************************

In general, it should be possible to use Cantera with any build system by
specifying the appropriate header and library paths, and specifying the required
libraries when linking. It is also necessary to specify the paths for libraries
used by Cantera, e.g. Sundials, BLAS, and LAPACK.

pkg-config
==========

On systems where the ``pkg-config`` program is installed, it can be used to
determine the correct compiler and linker flags for use with Cantera. For
example:

.. code-block:: bash

   g++ myProgram.cpp -o myProgram $(pkg-config --cflags --libs cantera)

It can also be used to populate variables in a Makefile:

.. code-block:: make

    CFLAGS += $(shell pkg-config --cflags cantera)
    LIBS += $(shell pkg-config --libs cantera)

Or in an SConstruct file::

    env.ParseConfig("pkg-config --cflags --libs cantera")

Note that ``pkg-config`` will work only if it can find the ``cantera.pc``
file. If Cantera's libraries are not installed in a standard location such as
``/usr/lib`` or ``/usr/local/lib``, you may need to set the ``PKG_CONFIG_PATH``
environment variable appropriately before using ``pkg-config``.

SCons
=====

SCons is a multi-platform, Python-based build system. It is the build system
used to compile Cantera. The description of how to build a project is contained
in a file named ``SConstruct``. The ``SConstruct`` file is actually a Python
script, which makes it very straightforward to add functionality to a
SCons-based build system.

A typical ``SConstruct`` file for compiling a program that uses Cantera might
look like this::

    env = Environment()

    env.Append(CCFLAGS='-g',
               CPPPATH=['/usr/local/cantera/include', 
                        '/usr/local/sundials/include'],
	       LIBS=['cantera', 'sundials_cvodes', 'sundials_ida', 
                     'sundials_nvecserial', 'lapack', 'blas'],
               LIBPATH=['/usr/local/cantera/lib',
                        '/usr/local/sundials/lib'],
               LINKFLAGS=['-g', '-pthread'])

    sample = env.Program('sample', 'sample.cpp')
    Default(sample)

This script establishes what SCons refers to as a "construction environment"
named ``env``, and sets the header (``CPPPATH``) and library (``LIBPATH``) paths
to include the directories containing the Cantera headers and libraries, as well
as libraries that Cantera depends on, such as Sundials, BLAS, and LAPACK. Then,
a program named ``sample`` is compiled using the single source file
``sample.cpp``.

Several other example ``SConstruct`` files are included with the C++ examples
contained in the ``samples`` subdirectory of the Cantera installation directory.

For more information on SCons, see the `SCons Wiki <http://scons.org/wiki/>`_
and the `SCons homepage <http://www.scons.org>`_.

Make
====

Cantera is distributed with an "include Makefile" that can be used with
Make-based build systems. This file ``Cantera.mak`` is located in the
``samples`` subdirectory of the Cantera installation directory. To use it, add a
line referencing this file to the top of your Makefile::

    include path/to/Cantera.mak

The path specified should be the relative path from the ``Makefile`` to
``Cantera.mak``. This file defines several variables which can be used in your
Makefile. The following is an example ``Makefile`` that uses the definitions
contained in ``Cantera.mak``:

.. code-block:: makefile

    include ../../Cantera.mak

    CC=gcc
    CXX=g++
    RM=rm -f
    CCFLAGS=-g
    CPPFLAGS=$(CANTERA_INCLUDES)
    LDFLAGS=
    LDLIBS=$(CANTERA_LIBS)

    SRCS=sample.cpp
    OBJS=$(subst .cpp,.o,$(SRCS))

    all: sample

    kinetics1: $(OBJS)
	    $(CXX) $(LDFLAGS) -o sample $(OBJS) $(LDLIBS)

    clean:
	    $(RM) $(OBJS)

    dist-clean: clean
	    $(RM) *~

