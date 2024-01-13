(sec-compiling-cplusplus)=
# Compiling Cantera C++ Applications

This guide shows you how to build C++ programs that use Cantera's features.

:::{seealso}
This section is about compiling applications that use Cantera; instructions to compile
the Cantera library itself are [over here](sec-compiling).
:::

## Build Systems

In general, it should be possible to use Cantera with any build system by specifying the
necessary header and library paths, as well as the appropriate compiler and linker
options. The correct options and paths depend on your system configuration and the
options that were used to compile Cantera. Cantera is installed with examples and
configuration data that will help with determining the correct options for some common
build systems, namely CMake and SCons.

### CMake

[CMake](https://cmake.org/) is a multi-platform build system that uses a high-level
project description to generate platform-specific build scripts (for example, on Linux,
CMake will generate Makefiles, and on Windows, it can generate Visual Studio `.sln`
files). The configuration file for a CMake project is called `CMakeLists.txt`. A typical
`CMakeLists.txt` file for compiling a program that uses Cantera might look like this:

```cmake
cmake_minimum_required(VERSION 3.1)
project (sample)

set(CMAKE_VERBOSE_MAKEFILE ON)
set(CMAKE_CXX_STANDARD 17)

find_package(Threads REQUIRED)

include_directories("/opt/cantera/include")
link_directories("/opt/cantera/lib")

add_executable(sample sample.cpp)
target_link_libraries(sample cantera_shared Threads::Threads)
```

Several example `CMakeLists.txt` files are included with the C++ examples contained in
the `samples/cxx` subdirectory of the Cantera installation directory, which have the
paths and lists of libraries correctly configured for the system on which they are
installed.

To build a program using CMake on Linux or macOS, run the following commands from the
directory containing the `CMakeLists.txt` file:

```bash
mkdir build
cd build
cmake ..
cmake --build .
```

This will create an executable named `sample` in the `build` directory.

To build a program using CMake on Windows, run the following commands from the directory
containing the `CMakeLists.txt` file:

```bash
mkdir build
cd build
cmake ..
cmake --build . --config Release
```

This will create an executable named `sample.exe` in the `build\Release` directory.

### SCons

[SCons](https://scons.org/) is a multi-platform, Python-based build system. It is the
build system used to compile Cantera. The description of how to build a project is
contained in a file named `SConstruct`. The `SConstruct` file is actually a Python
script, which makes it very straightforward to add functionality to a SCons-based build
system.

A typical `SConstruct` file for compiling a program that uses Cantera might look like
this:

```python
env = Environment()

env.Append(CCFLAGS='-g -std=c++17',
           CPPPATH=['/usr/local/cantera/include'],
           LIBS=['cantera_shared'],
           LIBPATH=['/usr/local/cantera/lib'],
           RPATH=['/usr/local/cantera/lib']
           LINKFLAGS=['-g', '-pthread'])

sample = env.Program('sample', 'sample.cpp')
Default(sample)
```

This script establishes what SCons refers to as a "construction environment" named
`env`, and sets the header (`CPPPATH`) and library (`LIBPATH`) paths to include the
directories containing the Cantera headers and libraries. Then, a program named `sample`
is compiled using the single source file `sample.cpp`.

To determine the appropriate settings for your system, take a look at one of the
pre-configured `SConstruct` files that are provided with the C++ examples contained in
the `samples/cxx` subdirectory of the Cantera installation directory.

To build a program using SCons, simply run the following command from a shell in the
directory containing the `SConstruct` file:

```bash
scons
```

:::{tip}
If you installed SCons using Conda, you may need to activate the appropriate Conda
environment so that the `scons` command will be on your path. On Windows, you may need
to run this command from a shell with the appropriate Visual Studio environment
variables set. This can be done either by starting the shell using the *Developer
Command Prompt for VS 20xx* shortcut in the Start menu, or by running the batch file:

```bat
C:\Program Files\Visual Studio 2022\VC\Auxiliary\Build\vcvars64.bat
```

in an existing shell, where the path specified will depend on the version and
installation path of Visual Studio.
:::

:::{seealso}
For more information on SCons, see the [SCons Wiki](https://github.com/SCons/scons/wiki/)
and the [SCons homepage](https://www.scons.org).
:::

### pkg-config

On systems where the `pkg-config` program is installed, it can be used to determine the
correct compiler and linker flags for use with Cantera. For example:

```bash
g++ myProgram.cpp -o myProgram $(pkg-config --cflags --libs cantera)
```

It can also be used to populate variables in a Makefile:

```make
CFLAGS += $(shell pkg-config --cflags cantera)
LIBS += $(shell pkg-config --libs cantera)
```

Or in an SConstruct file:

```python
env.ParseConfig("pkg-config --cflags --libs cantera")
```

````{tip}
`pkg-config` will work only if it can find the `cantera.pc` file. If Cantera's libraries
are not installed in a standard location such as `/usr/lib` or `/usr/local/lib`, you may
need to set the `PKG_CONFIG_PATH` environment variable appropriately before using
`pkg-config`, for example by running:

```bash
export PKG_CONFIG_PATH=$PKG_CONFIG_PATH:/path/to/lib/pkgconfig
```

where `/path/to/lib` should be replaced by Cantera's library installation path.
````

`pkg-config` is available through system package managers for most Linux distributions,
and can be installed using Conda or Homebrew on macOS.

## Libraries & Library Paths

### Choosing Libraries During Compilation

Applications can be linked to either the Cantera static library or dynamically linked to
the Cantera shared library. Dynamic linking is recommended generally, and required to
enable features such as the use of `ExtensibleRate` objects. The pre-configured
`CMakelists.txt` and `SConstruct` files included with the Cantera examples are set up to
use dynamic linking.

#### The Cantera Library

If Cantera was compiled with the `renamed_shared_libraries=y` option, then you can link
to the Cantera shared library by specifying the library name `cantera_shared` or to the
static library by specifying the library name `cantera`. If Cantera was compiled with
the `renamed_shared_libraries=n` option, then you can link to the shared library by
specifying the library named `cantera`.

The `renamed_shared_libraries=y` option is the default if you compiled Cantera yourself,
or if you installed packages for Windows. Cantera packages for Conda and Ubuntu use the
setting `renamed_shared_libraries=n`.

#### Additional Dependencies

If you link to the Cantera shared library, you only need to link to that and any of your
program's direct dependencies. You do not need to link to any of Cantera's dependencies
unless your program also uses them directly. One unexpected direct dependency your
program may have is on the `fmt` library, due to its use in C++ templates in Cantera.

If you link to the Cantera static library, you will also need to specify all of
Cantera's library dependencies when linking your program, as well as the directories
containing these libraries (if they are not in standard search directories).

### Runtime Library Paths

Your operating system needs to be able to find the shared library dependencies of your
program when it is run. This process is dependent primarily on your operating system.

#### Linux & macOS

If you linked to the Cantera shared library, you will need to provide the information
needed to find the Cantera library; the Cantera library then contains the information
needed to find its own dependencies such as SUNDIALS, LAPACK, and yaml-cpp. If you
linked to the Cantera static library, your program depends directly on Cantera's
dependencies instead, and you need to provide the information on where to find these
dependencies when you run your program.

There are several options for specifying library search paths:

1. Specify the "rpath" when compiling and linking your program. This is done with the
   compiler option for GCC/Clang `-Wl,-rpath,/path/to/libdir`, where `/path/to/libdir`
   is the directory containing the Cantera shared library. The build scripts provided
   with Cantera's examples are configured to use this option.
2. If the libraries are installed into a standard system location, such as `/usr/lib`
   or `/usr/local/lib` on Linux, they should be found automatically.
3. Set the `LD_LIBRARY_PATH` (Linux) or `DYLD_LIBRARY_PATH` (macOS) environment variable
   before running your program. For example, on Linux, use the command:

```bash
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/path/to/libdir
```

#### Windows

On Windows, all shared library (DLL) dependencies need to be on the `PATH`. You can add
the Cantera library directory to the `PATH` temporarily, for a single command prompt
session, by running a command like:

```bat
set PATH=%PATH%;C:\Program Files\Cantera\bin
```

where the path added depends on where you installed Cantera.
