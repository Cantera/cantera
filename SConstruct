"""
SCons build script for Cantera

Basic usage:

    'scons help' - Show this help message.

    'scons build' - Compile Cantera and the language interfaces using
                    default options.

    'scons clean' - Delete files created while building Cantera.

    'scons install' - Install Cantera.

    'scons uninstall' - Uninstall Cantera.

    'scons test' - Run all tests which did not previously pass or for which the
                   results may have changed.

    'scons test-reset' - Reset the passing status of all tests.

    'scons test-clean' - Delete files created while running the tests.

    'scons test-help' - List available tests.

    'scons test-NAME' - Run the test named "NAME".

    'scons samples' - Compile the C++ and Fortran samples.

    'scons msi' - Build a Windows installer (.msi) for Cantera.

    'scons sphinx' - Build the Sphinx documentation

    'scons doxygen' - Build the Doxygen documentation

Additional command options:

    'scons help --options' - Print a description of user-specifiable options.

    'scons help --list-options' - Print formatted list of available options.

    'scons help --option=<opt>' - Print the description of a specific option
                                  with name <opt>, for example
                                  'scons help --option=prefix'

    'scons <command> -j<X>' - Short form of '--jobs=<X>'. Run the <command> using '<X>'
                              parallel jobs, for example 'scons build -j4'

    'scons <command> -s' - Short form of '--silent' (or '--quiet'). Suppresses output.
"""
# Note that 'scons help' supports additional command options that are intended for
# internal use (debugging or reST parsing of options) and thus are not listed above:
#  --defaults ... list default values for all supported platforms
#  --restructured-text ... format configuration as reST
#  --dev ... add '-dev' to reST output
#  --output=<fname> ... send output to file (reST only)
#
# Other features not listed above:
# 'scons sdist' - Build PyPI packages.
# 'scons <command> dump' - Dump the state of the SCons environment to the
#                          screen instead of doing <command>, for example
#                          'scons build dump'. For debugging purposes.

# This f-string is deliberately here to trigger a SyntaxError when
# SConstruct is parsed by Python 2. This seems to be the most robust
# and simplest option that will reliably trigger an error in Python 2
# and provide actionable feedback for users.
f"""
Cantera must be built using Python 3.8 or higher. You can invoke SCons by executing
    python3 `which scons`
followed by any desired options.
"""

from pathlib import Path
import sys
import os
import platform
import subprocess
import re
import json
import textwrap
try:
    from packaging.version import parse as parse_version
except ImportError:
    from pkg_resources import parse_version
import SCons

# ensure that Python and SCons versions are sufficient for the build process
EnsurePythonVersion(3, 7)
EnsureSConsVersion(3, 0, 0)

from buildutils import *

if not COMMAND_LINE_TARGETS:
    # Print usage help
    logger.error("Missing command argument: type 'scons help' for information.")
    sys.exit(1)

if os.name not in ["nt", "posix"]:
    logger.error(f"Unrecognized operating system {os.name!r}")
    sys.exit(1)

valid_commands = ("build", "clean", "install", "uninstall",
                  "help", "msi", "samples", "sphinx", "doxygen", "dump",
                  "sdist")

# set default logging level
if GetOption("silent"):
    logger.logger.setLevel("ERROR")
else:
    logger.logger.setLevel("WARNING")

for command in COMMAND_LINE_TARGETS:
    if command not in valid_commands and not command.startswith('test'):
        logger.error(f"Unrecognized command line target: {command!r}")
        sys.exit(1)

    # update default logging level
    if command in ["build", "dump"] and not GetOption("silent"):
        logger.logger.setLevel("INFO")

if "clean" in COMMAND_LINE_TARGETS:
    remove_directory("build")
    remove_directory("stage")
    remove_directory(".sconf_temp")
    remove_directory("test/work")
    remove_file(".sconsign.dblite")
    remove_file("include/cantera/base/config.h")
    remove_file("src/extensions/PythonExtensionManager.os")
    remove_file("src/extensions/delegator.h")
    remove_file("src/pch/system.h.gch")
    remove_directory("include/cantera/ext")
    remove_file("config.log")
    remove_directory("doc/sphinx/matlab/examples")
    remove_file("doc/sphinx/matlab/examples.rst")
    for name in Path("doc/sphinx/matlab/").glob("**/*.rst"):
        if name.name != "index.rst":
            remove_file(name)
    remove_directory("doc/sphinx/cython/examples")
    remove_file("doc/sphinx/cython/examples.rst")
    for name in Path(".").glob("*.msi"):
        remove_file(name)
    for name in Path("site_scons").glob("**/*.pyc"):
        remove_file(name)
    remove_file("interfaces/matlab/toolbox/cantera_shared.dll")
    remove_file("interfaces/matlab/Contents.m")
    remove_file("interfaces/matlab/ctpath.m")
    for name in Path("interfaces/matlab/toolbox").glob("ctmethods.*"):
        remove_file(name)

    logger.status("Done removing output files.", print_level=False)

    if COMMAND_LINE_TARGETS == ["clean"]:
        # Just exit if there's nothing else to do
        sys.exit(0)
    else:
        Alias("clean", [])

if "test-clean" in COMMAND_LINE_TARGETS:
    remove_directory("build/test")
    remove_directory("test/work")
    remove_directory("build/python_local")

python_version = f"{sys.version_info.major}.{sys.version_info.minor}"
logger.info(
    f"SCons {SCons.__version__} is using the following Python interpreter:\n"
    f"    {sys.executable} (Python {python_version})", print_level=False)

# ******************************************
# *** Specify defaults for SCons options ***
# ******************************************

windows_options = [
    Option(
        "msvc_version",
        """Version of Visual Studio to use. The default is the newest installed version.
        Note that since multiple MSVC toolsets can be installed for a single version of
        Visual Studio, you probably want to use ``msvc_toolset_version`` unless you
        specifically installed multiple versions of Visual Studio. Windows MSVC only.
        """,
        ""),
    Option(
        "msvc_toolset_version",
        """Version of the MSVC toolset to use. The default is the default version for
        the given ``msvc_version``. Note that the toolset selected here must be
        installed in the MSVC version selected by ``msvc_version``. The default
        toolsets associated with various Visual Studio versions are:

        * '14.1' ('14.1x'): Visual Studio 2017
        * '14.2' ('14.2x'): Visual Studio 2019
        * '14.3' ('14.3x'): Visual Studio 2022.

        For version numbers in parentheses, 'x' is a placeholder for a minor version
        number. Windows MSVC only.""",
        ""),
    EnumOption(
        "target_arch",
        """Target architecture. The default is the same architecture as the
           installed version of Python. Windows only.""",
        {"Windows": "amd64"},
        ("amd64", "x86")),
    EnumOption(
        "toolchain",
        """The preferred compiler toolchain. If MSVC is not on the path but
           'g++' is on the path, 'mingw' is used as a backup. Windows only.""",
        {"Windows": "msvc"},
        ("msvc", "mingw", "intel")),
]

config_options = [
    Option(
        "AR",
        "The archiver to use.",
        "${AR}"),
    Option(
        "CXX",
        "The C++ compiler to use.",
        "${CXX}"),
    Option(
        "cxx_flags",
        """Compiler flags passed to the C++ compiler only. Separate multiple
           options with spaces, for example, "cxx_flags='-g -Wextra -O3 --std=c++14'"
           """,
        {
            "cl": "/EHsc /std:c++17",
            "default": "-std=c++17"
        }),
    Option(
        "CC",
        "The C compiler to use. This is only used to compile CVODE.",
        "${CC}"),
    Option(
        "cc_flags",
        """Compiler flags passed to both the C and C++ compilers, regardless of
           optimization level.""",
        {
            "cl": "/MD /nologo /D_SCL_SECURE_NO_WARNINGS /D_CRT_SECURE_NO_WARNINGS",
            "clang": "-fcolor-diagnostics",
            "default": "",
        }),
    PathOption(
        "prefix",
        r"""Set this to the directory where Cantera should be installed. If the Python
           executable found during compilation is managed by 'conda', the installation
           'prefix' defaults to the corresponding environment and the 'conda' layout
           will be used for installation (specifying any of the options 'prefix',
           'python_prefix', 'python_cmd', or 'layout' will override this default). On
           Windows systems, '$ProgramFiles' typically refers to "C:\Program Files".""",
        {"Windows": r"$ProgramFiles\Cantera", "default": "/usr/local"},
        PathVariable.PathAccept),
    PathOption(
        "libdirname",
        """Set this to the directory where Cantera libraries should be installed.
           Some distributions (for example, Fedora/RHEL) use 'lib64' instead of 'lib'
           on 64-bit systems or could use some other library directory name instead of
           'lib', depending on architecture and profile (for example, Gentoo 'libx32'
           on x32 profile). If the user didn't set the 'libdirname' configuration
           variable, set it to the default value 'lib'""",
        "lib", PathVariable.PathAccept),
    EnumOption(
        "python_package",
        """If you plan to work in Python, then you need the 'full' Cantera Python
           package. If, on the other hand, you will only use Cantera from some
           other language (for example, MATLAB or Fortran 90/95) and only need Python
           to process YAML files, then you only need a 'minimal' subset of the
           package and Cython and NumPy are not necessary. The 'none' option
           doesn't install any components of the Python interface. The default
           behavior is to build the full Python module for whichever version of
           Python is running SCons if the required prerequisites (NumPy and
           Cython) are installed. Note: 'y' is a synonym for 'full' and 'n'
           is a synonym for 'none'.""",
        "default", ("full", "minimal", "none", "n", "y", "default")),
    BoolOption(
        "python_sdist",
        """Setting this option to True builds the Python sdist.""",
        False),
    PathOption(
        "python_cmd",
        """Cantera needs to know where to find the Python interpreter. If
           'PYTHON_CMD' is not set, then the configuration process will use the
           same Python interpreter being used by SCons.""",
        "${PYTHON_CMD}",
        PathVariable.PathAccept),
    PathOption(
        "python_prefix",
        """Use this option if you want to install the Cantera Python package to
           an alternate location. On Unix-like systems, the default is the same
           as the 'prefix' option. If the 'python_prefix' option is set to
           the empty string or the 'prefix' option is not set, then the package
           will be installed to the system default 'site-packages' directory.
           To install to the current user's 'site-packages' directory, use
           'python_prefix=USER'.""",
        {"default": ""},
        PathVariable.PathAccept),
    EnumOption(
        "matlab_toolbox",
        """This variable controls whether the MATLAB toolbox will be built. If
           set to 'y', you will also need to set the value of the 'matlab_path'
           variable. If set to 'default', the MATLAB toolbox will be built if
           'matlab_path' is set.""",
        "default", ("y", "n", "default")),
    PathOption(
        "matlab_path",
        """Path to the MATLAB install directory. This should be the directory
           containing the 'extern', 'bin', etc. subdirectories. Typical values
           are: "C:\\Program Files\\MATLAB\\R2021a" on Windows,
           "/Applications/MATLAB_R2021a.app" on macOS, or
           "/opt/MATLAB/R2021a" on Linux.""",
        "", PathVariable.PathAccept),
    EnumOption(
        "f90_interface",
        """This variable controls whether the Fortran 90/95 interface will be
           built. If set to 'default', the builder will look for a compatible
           Fortran compiler in the 'PATH' environment variable, and compile
           the Fortran 90 interface if one is found.""",
        "default", ("y", "n", "default")),
    PathOption(
        "FORTRAN",
        """The Fortran (90) compiler. If unspecified, the builder will look for a
           compatible compiler (pgfortran, gfortran, ifort, ifx, g95) in the 'PATH'
           environment variable. Used only for compiling the Fortran 90 interface.""",
        "", PathVariable.PathAccept),
    Option(
        "FORTRANFLAGS",
        "Compilation options for the Fortran (90) compiler.",
        "-O3"),
    BoolOption(
        "coverage",
        """Enable collection of code coverage information with gcov.
           Available only when compiling with gcc.""",
        False),
    BoolOption(
        "doxygen_docs",
        "Build HTML documentation for the C++ interface using Doxygen.",
        False),
    BoolOption(
        "sphinx_docs",
        "Build HTML documentation for Cantera using Sphinx.",
        False),
    PathOption(
        "sphinx_cmd",
        "Command to use for building the Sphinx documentation.",
        "sphinx-build", PathVariable.PathAccept),
    Option(
        "sphinx_options",
        """Options passed to the 'sphinx_cmd' command line. Separate multiple
           options with spaces, for example, "-W --keep-going".""",
        "-W --keep-going"),
    EnumOption(
        "system_eigen",
        """Select whether to use Eigen from a system installation ('y'), from a
           Git submodule ('n'), or to decide automatically ('default'). If Eigen
           is not installed directly into a system include directory, for example, it
           is installed in '/opt/include/eigen3/Eigen', then you will need to add
           '/opt/include/eigen3' to 'extra_inc_dirs'.
           """,
        "default", ("default", "y", "n")),
    EnumOption(
        "system_fmt",
        """Select whether to use the fmt library from a system installation
           ('y'), from a Git submodule ('n'), or to decide automatically
           ('default'). If you do not want to use the Git submodule and fmt
           is not installed directly into system include and library
           directories, then you will need to add those directories to
           'extra_inc_dirs' and 'extra_lib_dirs'. This installation of fmt
           must include the shared version of the library, for example,
           'libfmt.so'.""",
        "default", ("default", "y", "n")),
    EnumOption(
        "hdf_support",
        """Select whether to support HDF5 container files natively ('y'), disable HDF5
           support ('n'), or to decide automatically based on the system configuration
           ('default'). Native HDF5 support uses the HDF5 library as well as the
           header-only HighFive C++ wrapper (see option 'system_highfive'). Specifying
           'hdf_include' or 'hdf_libdir' changes the default to 'y'.""",
        "default", ("default", "y", "n")),
    PathOption(
        "hdf_include",
        """The directory where the HDF5 header files are installed. This should be the
           directory that contains files 'H5Version.h' and 'H5Public.h', amongst others.
           Not needed if the headers are installed in a standard location, for example,
           '/usr/include'.""",
        "", PathVariable.PathAccept),
    PathOption(
        "hdf_libdir",
        """The directory where the HDF5 libraries are installed. Not needed if the
           libraries are installed in a standard location, for example, '/usr/lib'.""",
        "", PathVariable.PathAccept),
    EnumOption(
        "system_highfive",
        """Select whether to use HighFive from a system installation ('y'), from a
           Git submodule ('n'), or to decide automatically ('default'). If HighFive
           is not installed directly into a system include directory, for example, it
           is installed in '/opt/include/HighFive', then you will need to add
           '/opt/include/HighFive' to 'extra_inc_dirs'.""",
        "default", ("default", "y", "n")),
    EnumOption(
        "system_yamlcpp",
        """Select whether to use the yaml-cpp library from a system installation
           ('y'), from a Git submodule ('n'), or to decide automatically
           ('default'). If yaml-cpp is not installed directly into system
           include and library directories, then you will need to add those
           directories to 'extra_inc_dirs' and 'extra_lib_dirs'.""",
        "default", ("default", "y", "n")),
    EnumOption(
        "system_sundials",
        """Select whether to use SUNDIALS from a system installation ('y'), from
           a Git submodule ('n'), or to decide automatically ('default').
           Specifying 'sundials_include' or 'sundials_libdir' changes the
           default to 'y'.""",
        "default", ("default", "y", "n")),
    PathOption(
        "sundials_include",
        """The directory where the SUNDIALS header files are installed. This
           should be the directory that contains the "cvodes", "nvector", etc.
           subdirectories. Not needed if the headers are installed in a
           standard location, for example, '/usr/include'.""",
        "", PathVariable.PathAccept),
    PathOption(
        "sundials_libdir",
        """The directory where the SUNDIALS static libraries are installed.
           Not needed if the libraries are installed in a standard location,
           for example, '/usr/lib'.""",
        "", PathVariable.PathAccept),
    EnumOption(
        "system_blas_lapack",
        """Select whether to use BLAS/LAPACK from a system installation ('y'), use
           Eigen linear algebra support ('n'), or to decide automatically based on
           libraries detected on the system ('default'). Specifying 'blas_lapack_libs'
           or 'blas_lapack_dir' changes the default to 'y', whereas installing the
           Matlab toolbox changes the default to 'n'. On macOS, the 'default' option
           uses the Accelerate framework, whereas on other operating systems the
           preferred option depends on the CPU manufacturer. In general, OpenBLAS
           ('openblas') is prioritized over standard libraries ('lapack,blas'), with
           Eigen being used if no suitable BLAS/LAPACK libraries are detected. On Intel
           CPU's, MKL (Windows: 'mkl_rt' / Linux: 'mkl_rt,dl') has the highest priority,
           followed by the other options. Note that Eigen is required whether or not
           BLAS/LAPACK libraries are used.""",
        "default", ("default", "y", "n")),
    Option(
        "blas_lapack_libs",
        """Cantera can use BLAS and LAPACK libraries installed on your system if you
           have optimized versions available (see option 'system_blas_lapack'). To use
           specific versions of BLAS and LAPACK, set 'blas_lapack_libs' to the the list
           of libraries that should be passed to the linker, separated by commas, for
           example, "lapack,blas" or "lapack,f77blas,cblas,atlas".""",
        ""),
    PathOption(
        "blas_lapack_dir",
        """Directory containing the libraries specified by 'blas_lapack_libs'. Not
           needed if the libraries are installed in a standard location, for example,
           '/usr/lib'.""",
        "", PathVariable.PathAccept),
    BoolOption(
        "lapack_ftn_trailing_underscore",
        """Controls whether the LAPACK functions have a trailing underscore
           in the Fortran libraries.""",
        True),
    BoolOption(
        "lapack_ftn_string_len_at_end",
        """Controls whether the LAPACK functions have the string length
           argument at the end of the argument list ('yes') or after
           each argument ('no') in the Fortran libraries.""",
        True),
    EnumOption(
        "googletest",
        """Select whether to use gtest/gmock from system
           installation ('system'), from a Git submodule ('submodule'), to decide
           automatically ('default') or don't look for gtest/gmock ('none')
           and don't run tests that depend on gtest/gmock.""",
        "default", ("default", "system", "submodule", "none")),
    Option(
        "env_vars",
        """Environment variables to propagate through to SCons. Either the
           string 'all' or a comma separated list of variable names, for example,
           'LD_LIBRARY_PATH,HOME'.""",
        "PATH,LD_LIBRARY_PATH,DYLD_LIBRARY_PATH,PYTHONPATH,USERPROFILE"),
    BoolOption(
        "use_pch",
        "Use a precompiled-header to speed up compilation",
        {"icc": False, "default": True}),
    Option(
        "pch_flags",
        "Compiler flags when using precompiled-header.",
        {
            "cl": "/FIpch/system.h",
            "gcc": "-include src/pch/system.h",
            "icx": "-include-pch src/pch/system.h.gch",
            "clang": "-include-pch src/pch/system.h.gch",
            "default": "",
        }),
    Option(
        "thread_flags",
        "Compiler and linker flags for POSIX multithreading support.",
        {"Windows": "", "macOS": "", "default": "-pthread"}),
    BoolOption(
        "optimize",
        """Enable extra compiler optimizations specified by the
           'optimize_flags' variable, instead of the flags specified by the
           'no_optimize_flags' variable.""",
        True),
    Option(
        "optimize_flags",
        "Additional compiler flags passed to the C/C++ compiler when 'optimize=yes'.",
        {
            "cl": "/O2",
            "icc": "-O3 -fp-model precise",
            "icx": "-O3 -fp-model precise", # cannot assume finite math
            "gcc": "-O3 -Wno-inline",
            "default": "-O3",
        }),
    Option(
        "no_optimize_flags",
        "Additional compiler flags passed to the C/C++ compiler when 'optimize=no'.",
        {"cl": "/Od /Ob0", "default": "-O0"}),
    BoolOption(
        "debug",
        "Enable compiler debugging symbols.",
        True),
    Option(
        "debug_flags",
        "Additional compiler flags passed to the C/C++ compiler when 'debug=yes'.",
        {"cl": "/Zi /Fd${TARGET}.pdb", "default": "-g"}),
    Option(
        "no_debug_flags",
        "Additional compiler flags passed to the C/C++ compiler when 'debug=no'.",
        ""),
    Option(
        "debug_linker_flags",
        "Additional options passed to the linker when 'debug=yes'.",
        {"cl": "/DEBUG", "default": ""}),
    Option(
        "no_debug_linker_flags",
        "Additional options passed to the linker when 'debug=no'.",
        ""),
    Option(
        "warning_flags",
        """Additional compiler flags passed to the C/C++ compiler to enable
           extra warnings. Used only when compiling source code that is part
           of Cantera (for example, excluding code in the 'ext' directory).""",
        {
            "cl": "/W3",
            "default": "-Wall",
        }),
    Option(
        "extra_inc_dirs",
        """Additional directories to search for header files, with multiple
           directories separated by colons (*nix, macOS) or semicolons (Windows).
           If an active 'conda' environment is detected, the corresponding include
           path is automatically added.""",
        ""),
    Option(
        "extra_lib_dirs",
        """Additional directories to search for libraries, with multiple
           directories separated by colons (*nix, macOS) or semicolons (Windows).
           If an active 'conda' environment is detected, the corresponding library
           path is automatically added.""",
        ""),
    PathOption(
        "boost_inc_dir",
        """Location of the Boost header files. Not needed if the headers are
           installed in a standard location, for example, '/usr/include'.""",
        "",
        PathVariable.PathAccept),
    PathOption(
        "stage_dir",
        """Directory relative to the Cantera source directory to be
           used as a staging area for building for example, a Debian
           package. If specified, 'scons install' will install files
           to 'stage_dir/prefix/...'.""",
        "",
        PathVariable.PathAccept),
    BoolOption(
        "VERBOSE",
        """Create verbose output about what SCons is doing. Deprecated in Cantera 3.0
           and to be removed thereafter; replaceable by 'logging=debug'.
           """,
        False),
    EnumOption(
        "logging",
        """Select logging level for SCons output. By default, logging messages use
           the 'info' level for 'scons build' and the 'warning' level for all other
           commands. In case the SCons option '--silent' is passed, all messages below
           the 'error' level are suppressed.""",
        "default", ("debug", "info", "warning", "error", "default")),
    Option(
        "gtest_flags",
        """Additional options passed to each GTest test suite, for example,
           '--gtest_filter=*pattern*'. Separate multiple options with spaces.""",
        ""),
    BoolOption(
        "renamed_shared_libraries",
        """If this option is turned on, the shared libraries that are created
           will be renamed to have a '_shared' extension added to their base name.
           If not, the base names will be the same as the static libraries.
           In some cases this simplifies subsequent linking environments with
           static libraries and avoids a bug with using valgrind with
           the '-static' linking flag.""",
        True),
    BoolOption(
        "versioned_shared_library",
        """If enabled, create a versioned shared library, with symlinks to the
           more generic library name, for example, 'libcantera_shared.so.2.5.0' as the
           actual library and 'libcantera_shared.so' and 'libcantera_shared.so.2'
           as symlinks.""",
        {"mingw": False, "cl": False, "default": True}),
    BoolOption(
        "use_rpath_linkage",
        """If enabled, link to all shared libraries using 'rpath', that is, a fixed
           run-time search path for dynamic library loading.""",
        True),
    Option(
        "openmp_flag",
        """Compiler flags used for multiprocessing (only used to generate sample build
           scripts).""",
        {
            "cl": "/openmp",
            "icc": "-qopenmp",
            "icx": "-qopenmp",
            "apple-clang": "-Xpreprocessor -fopenmp",
            "default": "-fopenmp",
        }),
    EnumOption(
        "layout",
        """The layout of the directory structure. 'standard' installs files to
           several subdirectories under 'prefix', for example, 'prefix/bin',
           'prefix/include/cantera', 'prefix/lib' etc. This layout is best used in
           conjunction with "prefix='/usr/local'". 'compact' puts all installed files
           in the subdirectory defined by 'prefix'. This layout is best with a prefix
           like '/opt/cantera'. If the Python executable found during compilation is
           managed by 'conda', the layout will default to 'conda' irrespective of
           operating system. For the 'conda' layout, the Python package as well as all
           libraries and header files are installed into the active 'conda' environment.
           Input data, samples, and other files are installed in the 'shared/cantera'
           subdirectory of the active 'conda' environment.""",
        {"Windows": "compact", "default": "standard"},
        ("standard", "compact", "conda")),
    BoolOption(
        "package_build",
        """Used in combination with packaging tools (example: 'conda-build'). If
           enabled, the installed package will be independent from host and build
           environments, with all external library and include paths removed. Packaged
           C++ and Fortran samples assume that users will compile with local SDKs, which
           should be backwards compatible with the tools used for the build process.
        """,
        False),
    BoolOption(
        "fast_fail_tests",
        "If enabled, tests will exit at the first failure.",
        False),
    BoolOption(
        "skip_slow_tests",
        """If enabled, skip a subset of tests that are known to have long runtimes.
           Skipping these may be desirable when running with options that cause tests
           to run slowly, like disabling optimization or activating code profiling.""",
        False),
    BoolOption(
        "show_long_tests",
        "If enabled, duration of slowest tests will be shown.",
        False),
    BoolOption(
        "verbose_tests",
        "If enabled, verbose test output will be shown.",
        False),
]

config = Configuration()

if "help" in COMMAND_LINE_TARGETS:
    AddOption(
        "--options", dest="options",
        action="store_true", help="Print description of available options")
    AddOption(
        "--list-options", dest="list",
        action="store_true", help="List available options")
    AddOption(
        "--restructured-text", dest="rest",
        action="store_true", help="Format defaults as reST")
    AddOption(
        "--option", dest="option", nargs=1, type="string",
        action="store", help="Output help for specific option")
    AddOption(
        "--defaults", dest="defaults",
        action="store_true", help="All defaults (CLI only)")
    AddOption(
        "--dev", dest="dev",
        action="store_true", help="Append -dev (reST only)")
    AddOption(
        "--output", dest="output", nargs=1, type="string",
        action="store", help="Output file (reST only)")

    list = GetOption("list")
    rest = GetOption("rest")
    defaults = GetOption("defaults") is not None
    options = GetOption("options")
    option = GetOption("option")

    if not (list or rest or defaults or options or option):
        # show basic help information
        logger.status(__doc__, print_level=False)
        sys.exit(0)

    if defaults or rest or list:

        config.add(windows_options)
        config.add(config_options)

        if list:
            # show formatted list of options
            logger.status("\nConfiguration options for building Cantera:\n", print_level=False)
            logger.status(config.list_options() + "\n", print_level=False)
            sys.exit(0)

        if defaults:
            try:
                # print default values: if option is None, show description for all
                # available options, otherwise show description for specified option
                logger.status(config.help(option), print_level=False)
                sys.exit(0)
            except KeyError as err:
                message = "Run 'scons help --list-options' to see available options."
                logger.error(f"{err}.\n{message}")
                sys.exit(1)

        dev = GetOption("dev") is not None
        try:
            # format default values as reST: if option is None, all descriptions are
            # rendered, otherwise only the description of specified option is shown
            message = config.to_rest(option, dev=dev)
        except KeyError as err:
            message = "Run 'scons help --list-options' to see available options."
            logger.error(f"{err}.\n{message}")
            sys.exit(1)

        output = GetOption("output")
        if output:
            # write output to file
            output_file = Path(output).with_suffix(".rst")
            with open(output_file, "w+") as fid:
                fid.write(message)

            logger.status(f"Done writing output options to {output_file!r}.",
                          print_level=False)

        else:
            logger.status(message, print_level=False)

        sys.exit(0)

# **************************************
# *** Read user-configurable options ***
# **************************************

opts = Variables("cantera.conf")

extraEnvArgs = {}

if os.name == "nt":
    config.add(windows_options)
    config.add(config_options)

    config["prefix"].default = (Path(os.environ["ProgramFiles"]) / "Cantera").as_posix()
    config.select("Windows")

    # On Windows, target the same architecture as the current copy of Python,
    # unless the user specified another option.
    if "64 bit" not in sys.version:
        config["target_arch"].default = "x86"

    opts.AddVariables(*config.to_scons(("msvc_version", "msvc_toolset_version", "target_arch")))

    windows_compiler_env = Environment()
    opts.Update(windows_compiler_env)

    # Make an educated guess about the right default compiler
    if which("g++") and not which("cl.exe"):
        config["toolchain"].default = "mingw"

    if windows_compiler_env["msvc_version"] or windows_compiler_env["msvc_toolset_version"]:
        config["toolchain"].default = "msvc"

    opts.AddVariables(*config.to_scons("toolchain"))
    opts.Update(windows_compiler_env)

    if windows_compiler_env["toolchain"] == "msvc":
        toolchain = ["default"]
        if windows_compiler_env["msvc_version"]:
            extraEnvArgs["MSVC_VERSION"] = windows_compiler_env["msvc_version"]
        if windows_compiler_env["msvc_toolset_version"]:
            extraEnvArgs["MSVC_TOOLSET_VERSION"] = windows_compiler_env["msvc_toolset_version"]
        msvc_version = (windows_compiler_env["msvc_version"] or
                        windows_compiler_env.get("MSVC_VERSION"))
        logger.info(f"Compiling with MSVC version {msvc_version}", print_level=False)
        msvc_toolset = (windows_compiler_env["msvc_toolset_version"] or
                        windows_compiler_env.get("MSVC_TOOLSET_VERSION") or
                        f"{msvc_version} (default)")
        logger.info(f"Compiling with MSVC toolset {msvc_toolset}", print_level=False)

    elif windows_compiler_env["toolchain"] == "mingw":
        toolchain = ["mingw", "f90"]
        extraEnvArgs["F77"] = None
        # Next line fixes https://github.com/SCons/scons/issues/2683
        extraEnvArgs["WINDOWS_INSERT_DEF"] = 1

    elif windows_compiler_env["toolchain"] == "intel":
        toolchain = ["intelc"] # note: untested

    extraEnvArgs["TARGET_ARCH"] = windows_compiler_env["target_arch"]
    logger.info(f"Compiling for architecture: {windows_compiler_env['target_arch']}",
        print_level=False)
    logger.info(f"Compiling using the following toolchain(s): {repr(toolchain)}",
        print_level=False)
else:
    config.add(config_options)
    toolchain = ["default"]

env = Environment(tools=toolchain+["textfile", "subst", "recursiveInstall", "UnitsInterfaceBuilder", "wix", "gch"],
                  ENV={"PATH": os.environ["PATH"]},
                  toolchain=toolchain,
                  **extraEnvArgs)

env["OS"] = platform.system()
env["OS_BITS"] = int(platform.architecture()[0][:2])
if "cygwin" in env["OS"].lower():
    logger.error(f"Error: Operating system {os.name!r} is no longer supported.")
    sys.exit(1)

if "FRAMEWORKS" not in env:
    env["FRAMEWORKS"] = []

# Needed for Matlab to source ~/.matlab7rc.sh
if "HOME" in os.environ:
    env["ENV"]["HOME"] = os.environ["HOME"]

if os.name == "nt":
    env["INSTALL_MANPAGES"] = False

    # Fixes a linker error in Windows
    if "TMP" in os.environ:
        env["ENV"]["TMP"] = os.environ["TMP"]

    # Fixes issues with Python subprocesses. See http://bugs.python.org/issue13524
    env["ENV"]["SystemRoot"] = os.environ["SystemRoot"]

    # Fix an issue with Unicode sneaking into the environment on Windows
    for key,val in env["ENV"].items():
        env["ENV"][key] = str(val)

else:
    env["INSTALL_MANPAGES"] = True

add_RegressionTest(env)

opts.AddVariables(*config.to_scons(["AR", "CC", "CXX"], env=env))
opts.Update(env)

# Check if this is actually Apple's clang on macOS
env["using_apple_clang"] = False
if env["OS"] == "Darwin":
    result = subprocess.check_output([env.subst("$CC"), "--version"]).decode("utf-8")
    if "clang" in result.lower() and ("Xcode" in result or "Apple" in result):
        env["using_apple_clang"] = True
        config.select("apple-clang")

if "gcc" in env.subst("$CC") or "gnu-cc" in env.subst("$CC"):
    config.select("gcc")

elif env["CC"] == "cl": # Visual Studio
    config.select("cl")

elif "icc" in env.subst("$CC"):
    logger.warning("Support for the deprecated Intel compiler suite (icc/icpc) "
                   "will be removed after Cantera 3.0.\nConsider using the new "
                   "LLVM-based Intel oneAPI compilers (icx/icpx) instead.")
    config.select("icc")

elif "icx" in env.subst("$CC"):
    config.select("icx")

elif "clang" in env.subst("$CC"):
    config.select("clang")

else:
    logger.error(f"Unrecognized C compiler {env['CC']!r}")
    sys.exit(1)

if env["OS"] == "Windows":
    config.select("Windows")
elif env["OS"] == "Darwin":
    config.select("macOS")

# SHLIBVERSION fails with MinGW: http://scons.tigris.org/issues/show_bug.cgi?id=3035
if (env["toolchain"] == "mingw"):
    config.select("mingw")

config.select("default")
config["python_cmd"].default = sys.executable

opts.AddVariables(*config.to_scons())
opts.Update(env)
opts.Save('cantera.conf', env)

# Expand ~/ and environment variables used in cantera.conf (variables used on
# the command line will be expanded by the shell)
for option in opts.keys():
    original = env[option]
    if isinstance(original, str):
        modified = os.path.expandvars(os.path.expanduser(env[option]))
        if original != modified:
            logger.info(f"Expanding {original!r} to {modified!r}")
            env[option] = modified

if "help" in COMMAND_LINE_TARGETS:
    option = GetOption("option")
    try:
        # print configuration: if option is None, description is shown for all
        # options; otherwise description is shown for specified option
        logger.status(config.help(option, env=env), print_level=False)
        sys.exit(0)
    except KeyError as err:
        message = "Run 'scons help --list-options' to see available options."
        logger.error(f"{err}.\n{message}")
        sys.exit(1)

if 'doxygen' in COMMAND_LINE_TARGETS:
    env['doxygen_docs'] = True
if 'sphinx' in COMMAND_LINE_TARGETS:
    env['sphinx_docs'] = True
if "sdist" in COMMAND_LINE_TARGETS:
    env["python_sdist"] = True

    if env["python_package"] == "default":
        logger.info("'sdist' target was specified. Setting 'python_package' to none.")
        env["python_package"] = "none"
    elif env["python_package"] in ("full", "y"):
        logger.error("'sdist' target was specified. Cannot also build Python package.")
        sys.exit(1)
    for ext_dependency in ("sundials", "fmt", "yamlcpp", "eigen", "highfive"):
        if env[f"system_{ext_dependency}"] == "y":
            logger.error(f"'sdist' target was specified. Cannot use 'system_{ext_dependency}'.")
            sys.exit(1)
        else:
            env[f"system_{ext_dependency}"] = "n"

    logger.info("'sdist' target was specified. Setting 'use_pch' to False.")
    env["use_pch"] = False

for arg in ARGUMENTS:
    if arg not in config:
        logger.error(f"Encountered unexpected command line option: {arg!r}")
        sys.exit(1)

env["cantera_version"] = "3.0.0"
# For use where pre-release tags are not permitted (MSI, sonames)
env['cantera_pure_version'] = re.match(r'(\d+\.\d+\.\d+)', env['cantera_version']).group(0)
env['cantera_short_version'] = re.match(r'(\d+\.\d+)', env['cantera_version']).group(0)

def get_processor_name():
    """Check processor name"""
    # adapted from:
    # https://stackoverflow.com/questions/4842448/getting-processor-information-in-python
    if platform.system() == "Windows":
        return platform.processor().strip()
    elif platform.system() == "Darwin":
        os.environ['PATH'] = os.environ['PATH'] + os.pathsep + '/usr/sbin'
        command ="sysctl -n machdep.cpu.brand_string"
        return subprocess.check_output(command, shell=True).decode().strip()
    elif platform.system() == "Linux":
        command = "cat /proc/cpuinfo"
        all_info = subprocess.check_output(command, shell=True).decode().strip()
        for line in all_info.split("\n"):
            if "model name" in line:
                return re.sub(".*model name.*:", "", line, 1).strip()
    return ""

env["CPU"] = get_processor_name()
logger.info(f"Compiling on {env['CPU']!r}")

try:
    env["git_commit"] = get_command_output("git", "rev-parse", "--short", "HEAD")
    logger.info(f"Building Cantera from git commit {env['git_commit']!r}")
except (subprocess.CalledProcessError, FileNotFoundError):
    env["git_commit"] = "unknown"

# Print values of all build options:
# the (updated) "cantera.conf" combines all options that were specified by the user
cantera_conf = Path("cantera.conf").read_text()
logger.info("Configuration variables read from 'cantera.conf' and command line:")
logger.info(textwrap.indent(cantera_conf, "    "), print_level=False)

# ********************************************
# *** Configure system-specific properties ***
# ********************************************

loglevel = env["logging"]
if loglevel != "default":
    logger.logger.setLevel(loglevel.upper())

if env["VERBOSE"]:
    # @todo: Remove after Cantera 3.0
    logger.warning("Option 'VERBOSE' is deprecated: replaceable by 'logging=debug'")
    logger.logger.setLevel("DEBUG")

# Copy in external environment variables
if env['env_vars'] == 'all':
    env['ENV'].update(os.environ)
    if 'PYTHONHOME' in env['ENV']:
        del env['ENV']['PYTHONHOME']
elif env['env_vars']:
    for name in env['env_vars'].split(','):
        if name in os.environ:
            if name == 'PATH':
                env.AppendENVPath('PATH', os.environ['PATH'])
            else:
                env['ENV'][name] = os.environ[name]
            logger.debug(f"Propagating environment variable {name}={env['ENV'][name]}")
        elif name not in config["env_vars"].default.split(','):
            logger.warning(f"Failed to propagate environment variable {name!r}\n"
                           "Edit cantera.conf or the build command line to fix this.")

env["extra_inc_dirs"] = [Path(d).as_posix() for d in env["extra_inc_dirs"].split(os.pathsep) if d]
env["extra_lib_dirs"] = [Path(d).as_posix() for d in env["extra_lib_dirs"].split(os.pathsep) if d]

# Add conda library/include paths (if applicable) to extra
conda_prefix = os.environ.get("CONDA_PREFIX")
if conda_prefix is not None:
    conda_prefix = Path(conda_prefix)
    if os.name == "nt":
        conda_inc_dir = (conda_prefix / "Library" / "include").as_posix()
        conda_lib_dir = (conda_prefix / "Library" / env["libdirname"]).as_posix()
    else:
        conda_inc_dir = (conda_prefix / "include").as_posix()
        conda_lib_dir = (conda_prefix / env["libdirname"]).as_posix()
    env["extra_inc_dirs"].append(conda_inc_dir)
    env["extra_lib_dirs"].append(conda_lib_dir)
    logger.info(f"Adding conda include and library paths: {conda_prefix}")

add_system_include(env, env['extra_inc_dirs'])
env.Append(LIBPATH=env['extra_lib_dirs'])

if env['use_rpath_linkage']:
    env.Append(RPATH=env['extra_lib_dirs'])

if env['CC'] == 'cl':
    # embed manifest file
    env['LINKCOM'] = [env['LINKCOM'],
                      'if exist ${TARGET}.manifest mt.exe -nologo -manifest ${TARGET}.manifest -outputresource:$TARGET;1']
    env['SHLINKCOM'] = [env['SHLINKCOM'],
                        'if exist ${TARGET}.manifest mt.exe -nologo -manifest ${TARGET}.manifest -outputresource:$TARGET;2']
    env['FORTRAN_LINK'] = 'link'
else:
    env['FORTRAN_LINK'] = '$FORTRAN'

if env['boost_inc_dir']:
    env["boost_inc_dir"] = Path(env["boost_inc_dir"]).as_posix()
    add_system_include(env, env['boost_inc_dir'])

if env['blas_lapack_dir']:
    env["blas_lapack_dir"] = Path(env["blas_lapack_dir"]).as_posix()
    env.Append(LIBPATH=[env['blas_lapack_dir']])
    if env['use_rpath_linkage']:
        env.Append(RPATH=env['blas_lapack_dir'])

if env['system_sundials'] in ('y','default'):
    if env['sundials_include']:
        env["sundials_include"] = Path(env["sundials_include"]).as_posix()
        add_system_include(env, env['sundials_include'])
        env['system_sundials'] = 'y'
    if env['sundials_libdir']:
        env["sundials_libdir"] = Path(env["sundials_libdir"]).as_posix()
        env.Append(LIBPATH=[env['sundials_libdir']])
        env['system_sundials'] = 'y'
        if env['use_rpath_linkage']:
            env.Append(RPATH=env['sundials_libdir'])

# BLAS / LAPACK configuration
if env['blas_lapack_libs'] != '':
    env['blas_lapack_libs'] = env['blas_lapack_libs'].split(',')
    env['use_lapack'] = True
else:
    env['blas_lapack_libs'] = []
    env['use_lapack'] = False

# ************************************
# *** Compiler Configuration Tests ***
# ************************************

def CheckStatement(context, function, includes=""):
    context.Message('Checking for %s... ' % function)
    src = """
%(include)s
int main(int argc, char** argv) {
    %(func)s;
    return 0;
}
""" % {'func':function, 'include':includes}
    result = context.TryCompile(src, '.cpp')
    context.Result(result)
    return result

# Set up compiler options before running configuration tests
env.Append(CXXFLAGS=listify(env['cxx_flags']))
env['CCFLAGS'] = listify(env['cc_flags']) + listify(env['thread_flags'])
env['LINKFLAGS'] += listify(env['thread_flags'])
env['CPPDEFINES'] = {}

env['warning_flags'] = listify(env['warning_flags'])
env["pch_flags"] = listify(env["pch_flags"])
env["openmp_flag"] = listify(env["openmp_flag"])

if env['optimize']:
    env['CCFLAGS'] += listify(env['optimize_flags'])
    env.Append(CPPDEFINES=['NDEBUG'])
else:
    env['CCFLAGS'] += listify(env['no_optimize_flags'])

if env['debug']:
    env['CCFLAGS'] += listify(env['debug_flags'])
    env['LINKFLAGS'] += listify(env['debug_linker_flags'])
else:
    env['CCFLAGS'] += listify(env['no_debug_flags'])
    env['LINKFLAGS'] += listify(env['no_debug_linker_flags'])

if env['coverage']:
    if  'gcc' in env.subst('$CC') or 'clang' in env.subst('$CC'):
        env.Append(CCFLAGS=['-fprofile-arcs', '-ftest-coverage'])
        env.Append(LINKFLAGS=['-fprofile-arcs', '-ftest-coverage'])

    else:
        logger.error("Coverage testing is only available with GCC.")
        exit(0)

def config_error(message):
    if env["logging"].lower() == "debug":
        logger.error(message)
        debug_message = [
            f"\n{' Contents of config.log: ':*^80}\n",
            open("config.log").read().strip(),
            f"\n{' End of config.log ':*^80}",
        ]
        logger.debug("\n".join(debug_message), print_level=False)
    else:
        error_message = [message]
        error_message.append("\nSee 'config.log' for details.")
        logger.error("\n".join(error_message))
    sys.exit(1)

conf = Configure(env, custom_tests={'CheckStatement': CheckStatement})
env = conf.env  # Retain updates to `env` after the end of the `Configure` context

# First, a sanity check:
if not conf.CheckCXXHeader("cmath", "<>"):
    config_error(
        "The C++ compiler is not correctly configured (incomplete include paths)."
    )

def get_expression_value(includes, expression, defines=()):
    s = ['#define ' + d for d in defines]
    s.extend('#include ' + i for i in includes)
    s.extend(('#define Q(x) #x',
              '#define QUOTE(x) Q(x)',
              '#include <iostream>',
              'int main(int argc, char** argv) {',
              '    std::cout << %s << std::endl;' % expression,
              '    return 0;',
              '}\n'))
    return '\n'.join(s)

# Check that libraries link correctly
cmath_check_source = get_expression_value(["<cmath>"], "cos(0. * argc)")
retcode, cmath_works = conf.TryRun(cmath_check_source, ".cpp")
if cmath_works.strip() != "1":
    config_error(
        "The C++ compiler is not correctly configured (failed at linking stage)."
    )

# Check that NaN is treated correctly
nan_check_source = get_expression_value(["<cmath>"], "std::isnan(NAN + argc)")
retcode, nan_works = conf.TryRun(nan_check_source, ".cpp")
if nan_works.strip() != "1":
    config_error(
        "Cantera requires a working implementation of 'std::isnan'.\n"
        "If you have specified '-ffast-math' or equivalent as an optimization option,\n"
        "either remove this option or add the '-fno-finite-math-only option'."
    )

def split_version(version):
    """Split integer version into version string."""
    version = divmod(float(version.strip()), 10000)
    (fmt_maj, (fmt_min, fmt_pat)) = version[0], divmod(version[1], 100)
    return f"{fmt_maj:.0f}.{fmt_min:.0f}.{fmt_pat:.0f}"

# Check for fmt library and checkout submodule if needed
if env['system_fmt'] in ('y', 'default'):
    fmt_version_source = get_expression_value(
        ['<fmt/format.h>'], 'FMT_VERSION', ['FMT_HEADER_ONLY'])
    retcode, fmt_lib_version = conf.TryRun(fmt_version_source, '.cpp')
    if retcode:
        fmt_lib_version = split_version(fmt_lib_version)
        fmt_min_version = "6.1.2"
        if parse_version(fmt_lib_version) < parse_version(fmt_min_version):
            if env['system_fmt'] == 'y':
                config_error(
                    f"System fmt version {fmt_lib_version} is not supported;"
                    f"version {fmt_min_version} or higher is required.")
                logger.info(
                    f"System fmt version {fmt_lib_version} is not supported; "
                    "using private installation instead.")
        else:
            env['system_fmt'] = True
            logger.info(f"Using system installation of fmt library.")
    elif env['system_fmt'] == 'y':
        config_error('Expected system installation of fmt library, but it '
            'could not be found.')

if env['system_fmt'] in ('n', 'default'):
    if not os.path.exists('ext/fmt/include/fmt/ostream.h'):
        if not os.path.exists('.git'):
            config_error('fmt is missing. Install source in ext/fmt.')

        try:
            code = subprocess.call(['git','submodule','update','--init',
                                    '--recursive','ext/fmt'])
        except Exception:
            code = -1
        if code:
            config_error('fmt submodule checkout failed.\n'
                         'Try manually checking out the submodule with:\n\n'
                         '    git submodule update --init --recursive ext/fmt\n')

    fmt_version_source = get_expression_value(
        ['"../ext/fmt/include/fmt/format.h"'], 'FMT_VERSION', ['FMT_HEADER_ONLY'])
    retcode, fmt_lib_version = conf.TryRun(fmt_version_source, '.cpp')
    if not retcode:
        config_error('Expected private installation of fmt library, but it is '
            'not configured correctly.')

    fmt_lib_version = split_version(fmt_lib_version)
    env['system_fmt'] = False
    logger.info(f"Using private installation of fmt library.")

if env["OS"] == "Windows" and parse_version(fmt_lib_version) < parse_version("8.0.0"):
    # Workaround for symbols not exported on Windows in older fmt versions
    env.Append(CPPDEFINES={"FMT_HEADER_ONLY": 1})

logger.info(f"Using fmt version {fmt_lib_version}")

# Check for yaml-cpp library and checkout submodule if needed
if env['system_yamlcpp'] in ('y', 'default'):
    # We need the Mark() function, which was added in version 0.5.3
    if conf.CheckStatement('YAML::Node().Mark()', '#include "yaml-cpp/yaml.h"'):
        env['system_yamlcpp'] = True
        logger.info("Using system installation of yaml-cpp library.")

    elif env['system_yamlcpp'] == 'y':
        config_error("Expected system installation of yaml-cpp library, but it "
            "could not be found or it is too old (0.6 or newer is required).")

if env['system_yamlcpp'] in ('n', 'default'):
    env['system_yamlcpp'] = False
    logger.info("Using private installation of yaml-cpp library.")
    if not os.path.exists('ext/yaml-cpp/include/yaml-cpp/yaml.h'):
        if not os.path.exists('.git'):
            config_error('yaml-cpp is missing. Install source in ext/yaml-cpp.')

        try:
            code = subprocess.call(['git', 'submodule', 'update', '--init',
                                    '--recursive', 'ext/yaml-cpp'])
        except Exception:
            code = -1
        if code:
            config_error('yaml-cpp submodule checkout failed.\n'
                         'Try manually checking out the submodule with:\n\n'
                         '    git submodule update --init --recursive ext/yaml-cpp\n')

# Check for googletest and checkout submodule if needed
if env['googletest'] in ('system', 'default'):
    has_gtest = conf.CheckCXXHeader('gtest/gtest.h', '""')
    has_gmock = conf.CheckCXXHeader('gmock/gmock.h', '""')
    if has_gtest and has_gmock:
        env['googletest'] = 'system'
        logger.info("Using system installation of Googletest")
    elif env['googletest'] == 'system':
        config_error('Expected system installation of Googletest-1.8.0, but it '
                     'could not be found.')

if env['googletest'] in ('submodule', 'default'):
    env['googletest'] = 'submodule'
    has_gtest = os.path.exists('ext/googletest/googletest/include/gtest/gtest.h')
    has_gmock = os.path.exists('ext/googletest/googlemock/include/gmock/gmock.h')
    if not (has_gtest and has_gmock):
        if not os.path.exists('.git'):
            config_error('Googletest is missing. Install source in ext/googletest.')

        try:
            code = subprocess.call(['git','submodule','update','--init',
                                    '--recursive','ext/googletest'])
        except Exception:
            code = -1
        if code:
            config_error('Googletest not found and submodule checkout failed.\n'
                         'Try manually checking out the submodule with:\n\n'
                         '    git submodule update --init --recursive ext/googletest\n')
    logger.info("Using Googletest from Git submodule")

if env['googletest'] == 'none':
    logger.info("Not using Googletest -- unable to run complete test suite")

# Check for Eigen and checkout submodule if needed
if env["system_eigen"] in ("y", "default"):
    if conf.CheckCXXHeader("eigen3/Eigen/Dense", "<>"):
        env["system_eigen"] = True
        env["system_eigen_prefixed"] = True
        logger.info("Using system installation of Eigen.")
        eigen_include = "<eigen3/Eigen/Core>"
    elif conf.CheckCXXHeader("Eigen/Dense", "<>"):
        env["system_eigen"] = True
        env["system_eigen_prefixed"] = False
        logger.info("Using system installation of Eigen.")
        eigen_include = "<Eigen/Core>"
    elif env["system_eigen"] == "y":
        config_error("Expected system installation of Eigen, but it "
                     "could not be found.")

if env["system_eigen"] in ("n", "default"):
    env["system_eigen"] = False
    logger.info("Using private installation of Eigen.")
    if not os.path.exists("ext/eigen/Eigen/Dense"):
        if not os.path.exists(".git"):
            config_error("Eigen is missing. Install Eigen in ext/eigen.")

        try:
            code = subprocess.call(["git","submodule","update","--init",
                                    "--recursive","ext/eigen"])
        except Exception:
            code = -1
        if code:
            config_error("Eigen not found and submodule checkout failed.\n"
                         "Try manually checking out the submodule with:\n\n"
                         "    git submodule update --init --recursive ext/eigen\n")
    eigen_include = '"../ext/eigen/Eigen/Core"'

eigen_versions = 'QUOTE(EIGEN_WORLD_VERSION) "." QUOTE(EIGEN_MAJOR_VERSION) "." QUOTE(EIGEN_MINOR_VERSION)'
eigen_version_source = get_expression_value([eigen_include], eigen_versions)
retcode, eigen_lib_version = conf.TryRun(eigen_version_source, ".cpp")
env["EIGEN_LIB_VERSION"] = eigen_lib_version.strip()
logger.info(f"Found Eigen version {env['EIGEN_LIB_VERSION']}")

# Determine which standard library to link to when using Fortran to
# compile code that links to Cantera
if conf.CheckDeclaration('__GLIBCXX__', '#include <iostream>', 'C++'):
    env['cxx_stdlib'] = ['stdc++']
elif conf.CheckDeclaration('_LIBCPP_VERSION', '#include <iostream>', 'C++'):
    env['cxx_stdlib'] = ['c++']
else:
    env['cxx_stdlib'] = []

env['HAS_CLANG'] = conf.CheckDeclaration('__clang__', '', 'C++')
if not env["using_apple_clang"]:
    # This checks for these three libraries in order and stops when it finds the
    # first success. Intel = iomp5, LLVM/clang = omp, GCC = gomp. Since gomp is
    # likely to be installed on the system even if other compilers are installed
    # or in use, it needs to go last in the check.
    env['HAS_OPENMP'] = conf.CheckLibWithHeader(
        ["iomp5", "omp", "gomp"], "omp.h", language="C++"
    )
else:
    env["HAS_OPENMP"] = False
    logger.info("Not checking for OpenMP support due to using XCode compiler.")

boost_version_source = get_expression_value(['<boost/version.hpp>'], 'BOOST_LIB_VERSION')
retcode, boost_lib_version = conf.TryRun(boost_version_source, '.cpp')
env['BOOST_LIB_VERSION'] = '.'.join(boost_lib_version.strip().split('_'))
if not env['BOOST_LIB_VERSION']:
    config_error("Boost could not be found. Install Boost headers or set"
                 " 'boost_inc_dir' to point to the boost headers.")
else:
    logger.info(f"Found Boost version {env['BOOST_LIB_VERSION']}")
if parse_version(env['BOOST_LIB_VERSION']) < parse_version("1.70"):
    # Boost.DLL with std::filesystem (making it header-only) is available in Boost 1.70
    # or newer
    config_error("Cantera requires Boost version 1.70 or newer.")

# check BLAS/LAPACK installations
if env["system_blas_lapack"] == "n":
    env["blas_lapack_libs"] = []

elif env["blas_lapack_libs"] or env["blas_lapack_dir"]:
    if env["system_blas_lapack"] == "default":
        env["system_blas_lapack"] = "y"
    for lib in env["blas_lapack_libs"]:
        if not conf.CheckLib(lib, autoadd=False):
            config_error(f"Library {lib!r} could not be found.")

elif env["matlab_path"] != "" and env["matlab_toolbox"] in {"default", "y"}:
    # MATLAB provides the mwlapack and mwblas libraries in matlabroot/extern/lib.
    if env["system_blas_lapack"] == "default":
        env["system_blas_lapack"] = "n"

elif env["system_blas_lapack"] == "default":
    # auto-detect versions
    if env["OS"] == "Darwin":
        # Use macOS Accelerate framework by default
        blas_lapack_order = []
        env["use_lapack"] = True
        env.Append(FRAMEWORKS=["Accelerate"])
    elif "intel" in env["CPU"].lower():
        if env["OS"] == "Windows":
            blas_lapack_order = [["mkl_rt"], ["openblas"], ["lapack", "blas"]]
        else:
            blas_lapack_order = [["mkl_rt", "dl"], ["openblas"], ["lapack", "blas"]]
    else:
        # MKL is known to have deliberately sub-optimal performance on non-Intel
        # (i.e. AMD) processors
        blas_lapack_order = [["openblas"], ["lapack", "blas"]]
    for lib in blas_lapack_order:
        if all(conf.CheckLib(l, autoadd=False) for l in lib):
            env["blas_lapack_libs"] = lib
            env["use_lapack"] = True
            break
    if not env["blas_lapack_libs"] and env["OS"] != "Darwin":
        logger.info("No system BLAS/LAPACK libraries detected.")

if "mkl_rt" in env["blas_lapack_libs"]:
    mkl_version = textwrap.dedent("""\
        #include <iostream>
        #include "mkl.h"
        int main(int argc, char** argv) {
            MKLVersion Version;
            mkl_get_version(&Version);
            std::cout << Version.MajorVersion << "." << Version.MinorVersion << "."
                << Version.UpdateVersion << " for " << Version.Platform;
            return 0;
        }
    """)
    retcode, mkl_version = conf.TryRun(mkl_version, ".cpp")
    if retcode:
        logger.info(f"Using MKL {mkl_version}")
    else:
        logger.warning("Failed to determine MKL version.")

elif "openblas" in env["blas_lapack_libs"]:
    openblas_version_source = get_expression_value(["<openblas_config.h>"], "OPENBLAS_VERSION")
    retcode, openblas_version = conf.TryRun(openblas_version_source, ".cpp")
    if retcode:
        logger.info(f"Using {openblas_version.strip()}")
    else:
        logger.warning("Failed to determine OpenBLAS version.")

import SCons.Conftest, SCons.SConf
context = SCons.SConf.CheckContext(conf)

cvode_checks = [
    "CVodeCreate(CV_BDF, CV_NEWTON);",  # Sundials <= 3.2
    "CVodeCreate(CV_BDF);",  # Sundials>=4.0,<6.0
    "SUNContext ctx; SUNContext_Create(0, &ctx);"  # Sundials>=6.0
]
for cvode_call in cvode_checks:
    ret = SCons.Conftest.CheckLib(context,
                                  ['sundials_cvodes'],
                                  header='#include "cvodes/cvodes.h"',
                                  language='C++',
                                  call=cvode_call,
                                  autoadd=False,
                                  extra_libs=env['blas_lapack_libs'])
    # CheckLib returns False to indicate success
    if not ret:
        if env['system_sundials'] == 'default':
            env['system_sundials'] = 'y'
        break

# Execute if the cycle ends without 'break'
else:
    if env['system_sundials'] == 'default':
        env['system_sundials'] = 'n'
    elif env['system_sundials'] == 'y':
        config_error('Expected system installation of Sundials, but it could '
                     'not be found.')

# Checkout Sundials submodule if needed
if (env['system_sundials'] == 'n' and
    not os.path.exists('ext/sundials/include/cvodes/cvodes.h')):
    if not os.path.exists('.git'):
        config_error('Sundials is missing. Install source in ext/sundials.')

    try:
        code = subprocess.call(['git','submodule','update','--init',
                                '--recursive','ext/sundials'])
    except Exception:
        code = -1
    if code:
        config_error('Sundials not found and submodule checkout failed.\n'
                     'Try manually checking out the submodule with:\n\n'
                     '    git submodule update --init --recursive ext/sundials\n')

env['NEED_LIBM'] = not conf.CheckLibWithHeader(None, 'math.h', 'C',
                                               'double x; log(x);', False)
env['LIBM'] = ['m'] if env['NEED_LIBM'] else []

if env['system_sundials'] == 'y':
    for subdir in ('sundials', 'nvector', 'cvodes', 'idas', 'sunlinsol', 'sunmatrix',
                   'sunnonlinsol'):
        remove_directory('include/cantera/ext/' + subdir)

    # Determine Sundials version
    sundials_version_source = get_expression_value(['"sundials/sundials_config.h"'],
                                                   'QUOTE(SUNDIALS_VERSION)')
    retcode, sundials_version = conf.TryRun(sundials_version_source, '.cpp')
    if retcode == 0:
        config_error("Failed to determine Sundials version.")
    sundials_version = sundials_version.strip(' "\n')

    # Ignore the minor version, for example 2.4.x -> 2.4
    env['sundials_version'] = '.'.join(sundials_version.split('.')[:2])
    sundials_ver = parse_version(env['sundials_version'])
    if sundials_ver < parse_version("3.0") or sundials_ver >= parse_version("7.0"):
        logger.error(f"Sundials version {env['sundials_version']!r} is not supported.")
        sys.exit(1)
    elif sundials_ver > parse_version("6.6.0"):
        logger.warning(f"Sundials version {env['sundials_version']!r} has not been tested.")

    logger.info(f"Using system installation of Sundials version {sundials_version!r}.")

    # Determine whether or not Sundials was built with BLAS/LAPACK
    if sundials_ver < parse_version('5.5'):
        # In Sundials 2.6-5.5, SUNDIALS_BLAS_LAPACK is either defined or undefined
        env['has_sundials_lapack'] = conf.CheckDeclaration('SUNDIALS_BLAS_LAPACK',
                '#include "sundials/sundials_config.h"', 'C++')
    else:
        # In Sundials 5.5 and higher, two defines are included specific to the
        # SUNLINSOL packages indicating whether SUNDIALS has been built with LAPACK
        lapackband = conf.CheckDeclaration(
            "SUNDIALS_SUNLINSOL_LAPACKBAND",
            '#include "sundials/sundials_config.h"',
            "C++",
        )
        lapackdense = conf.CheckDeclaration(
            "SUNDIALS_SUNLINSOL_LAPACKDENSE",
            '#include "sundials/sundials_config.h"',
            "C++",
        )
        env["has_sundials_lapack"] = lapackband and lapackdense

    # In the case where a user is trying to link Cantera to an external BLAS/LAPACK
    # library, but Sundials was configured without this support, print a Warning.
    if not env['has_sundials_lapack'] and env['use_lapack']:
        logger.warning("External BLAS/LAPACK has been specified for Cantera "
                       "but Sundials was built without this support.")
else: # env['system_sundials'] == 'n'
    logger.info("Using private installation of Sundials version 5.3.")
    env['sundials_version'] = '5.3'
    env['has_sundials_lapack'] = int(env['use_lapack'])

if env["hdf_include"]:
    env["hdf_include"] = Path(env["hdf_include"]).as_posix()
    add_system_include(env, env["hdf_include"])
    env["hdf_support"] = "y"
if env["hdf_libdir"]:
    env["hdf_libdir"] = Path(env["hdf_libdir"]).as_posix()
    env.Append(LIBPATH=[env["hdf_libdir"]])
    env["hdf_support"] = "y"
    if env["use_rpath_linkage"]:
        env.Append(RPATH=env["hdf_libdir"])

if env["hdf_support"] == "n":
    env["use_hdf5"] = False
else:
    env["use_hdf5"] = conf.CheckLib("hdf5", autoadd=False)
    if not env["use_hdf5"] and env["hdf_support"] == "y":
        config_error("HDF5 support has been specified but libraries were not found.")

if env["use_hdf5"] and env["system_highfive"] in ("y", "default"):

    if conf.CheckLibWithHeader(
            "hdf5", "highfive/H5File.hpp", language="C++", autoadd=False):

        highfive_include = "<highfive/H5Version.hpp>"
        h5_version_source = get_expression_value(
            [highfive_include], "QUOTE(HIGHFIVE_VERSION)")
        retcode, h5_lib_version = conf.TryRun(h5_version_source, ".cpp")
        if retcode:
            if parse_version(h5_lib_version) < parse_version("2.5"):
                if env["system_highfive"] == "y":
                    config_error(
                        f"System HighFive version {h5_lib_version} is not "
                        "supported; version 2.5 or higher is required.")
                logger.info(
                    f"System HighFive version {h5_lib_version} is not supported.")
            else:
                env["system_highfive"] = True
            env["HIGHFIVE_VERSION"] = h5_lib_version.strip()
        else:
            config_error("Detected invalid HighFive configuration.")

        highfive_include = "<highfive/H5DataType.hpp>"
        logger.info("Using system installation of HighFive library.")

    elif env["system_highfive"] == "y":
        config_error(
            "Expected system installation of HighFive library, but it is either "
            "corrupted or could not be found.")

if env["use_hdf5"] and env["system_highfive"] in ("n", "default"):
    env["system_highfive"] = False
    if not Path("ext/HighFive/include").is_dir():
        if not os.path.exists(".git"):
            config_error("HighFive is missing. Install HighFive in ext/HighFive.")

        try:
            code = subprocess.call(["git", "submodule", "update", "--init",
                                    "--recursive", "ext/HighFive"])
        except Exception:
            code = -1
        if code:
            config_error("HighFive not found and submodule checkout failed.\n"
                        "Try manually checking out the submodule with:\n\n"
                        "    git submodule update --init --recursive ext/HighFive\n")

    def highfive_version(cmake_lists):
        """Read highfive version from CMakeLists.txt"""
        h5_version = Path(cmake_lists).read_text()
        h5_version = [line for line in h5_version.split("\n")
                      if line.startswith("project(HighFive")]
        return re.search('[0-9]+\.[0-9]+\.[0-9]+', h5_version[0]).group(0)

    env["HIGHFIVE_VERSION"] = highfive_version("ext/HighFive/CMakeLists.txt")
    highfive_include = '"../ext/HighFive/include/highfive/H5DataType.hpp"'
    logger.info(f"Using private installation of HighFive library.")

if env["use_hdf5"]:
    if conf.CheckStatement("HighFive::details::Boolean::HighFiveTrue",
                           f"#include {highfive_include}"):
        env["highfive_boolean"] = True

    hdf_version = textwrap.dedent("""\
        #include <iostream>
        #include "H5public.h"
        int main(int argc, char** argv) {
            std::cout << H5_VERS_MAJOR << "." << H5_VERS_MINOR << "." << H5_VERS_RELEASE;
            return 0;
        }
    """)
    retcode, hdf_version = conf.TryRun(hdf_version, ".cpp")
    if retcode:
        env["HDF_VERSION"] = hdf_version
        logger.info(
            f"Using HighFive version {env['HIGHFIVE_VERSION']} "
            f"for HDF5 {env['HDF_VERSION']}")
    else:
        logger.warning("Failed to determine HDF5 version.")

def set_fortran(pattern, value):
    # Set compiler / flags for all Fortran versions to be the same
    for version in ("FORTRAN", "F77", "F90", "F95", "F03", "F08"):
        env[pattern.format(version)] = value

# Try to find a working Fortran compiler:
def check_fortran(compiler, expected=False):
    hello_world = '''
program main
   write(*,'(a)') 'Hello, world!'
end program main
    '''
    if which(compiler):
        set_fortran("{}", compiler)
        success, output = conf.TryRun(hello_world, '.f90')
        if success and 'Hello, world!' in output:
            return True
        else:
            logger.warning(f"Unable to use {compiler!r} to compile the Fortran "
                           "interface. See config.log for details.")
            return False
    elif expected:
        logger.error(f"Could not find specified Fortran compiler: {compiler!r}")
        sys.exit(1)

    return False

set_fortran("{}FLAGS", env["FORTRANFLAGS"])

if env["using_apple_clang"] and env["f90_interface"] == "default":
    env["f90_interface"] = "n"

if env['f90_interface'] in ('y','default'):
    foundF90 = False
    if env['FORTRAN']:
        foundF90 = check_fortran(env['FORTRAN'], True)

    for compiler in ("pgfortran", "gfortran", "ifort", "ifx", "g95"):
        if foundF90:
            break
        foundF90 = check_fortran(compiler)

    if foundF90:
        logger.info(f"Using {env['FORTRAN']!r} to build the Fortran 90 interface")
        env['f90_interface'] = 'y'
    else:
        if env['f90_interface'] == 'y':
            logger.error("Could not find a suitable Fortran compiler to build the Fortran 90 interface.")
            sys.exit(1)
        else:
            env['f90_interface'] = 'n'
            env['FORTRAN'] = ''
            logger.info("Skipping compilation of the Fortran 90 interface.")

if 'pgfortran' in env['FORTRAN']:
    env['FORTRANMODDIRPREFIX'] = '-module '
elif 'gfortran' in env['FORTRAN']:
    env['FORTRANMODDIRPREFIX'] = '-J'
elif 'g95' in env['FORTRAN']:
    env['FORTRANMODDIRPREFIX'] = '-fmod='
elif 'ifort' in env['FORTRAN']:
    env['FORTRANMODDIRPREFIX'] = '-module '
elif "ifx" in env["FORTRAN"]:
    env["FORTRANMODDIRPREFIX"] = "-module "

set_fortran("{}", env["FORTRAN"])
set_fortran("SH{}", env["FORTRAN"])
env['FORTRANMODDIR'] = '${TARGET.dir}'

env = conf.Finish()

debug_message = [
    f"\n{' begin config.log ':-^80}\n",
    open("config.log").read().strip(),
    f"\n{' end config.log ':-^80}\n",
]
logger.debug("\n".join(debug_message), print_level=False)

env['python_cmd_esc'] = quoted(env['python_cmd'])

# Python Package Settings
python_min_version = parse_version("3.8")
# Newest Python version not supported/tested by Cantera
python_max_p1_version = parse_version("3.12")
# The string is used to set python_requires in setup.cfg.in
env["py_requires_ver_str"] = f">={python_min_version}"
if env["python_sdist"] or env["package_build"]:
    env["py_requires_ver_str"] += f",<{python_max_p1_version}"
cython_min_version = parse_version("0.29.31")
numpy_min_version = parse_version('1.12.0')

# We choose ruamel.yaml 0.15.34 as the minimum version
# since it is the highest version available in the Ubuntu
# 18.04 repositories and seems to work. Older versions such as
# 0.13.14 on CentOS7 and 0.10.23 on Ubuntu 16.04 raise an exception
# that they are missing the RoundTripRepresenter
ruamel_min_version = parse_version('0.15.34')

# Minimum pytest version assumed based on Ubuntu 20.04
pytest_min_version = parse_version("4.6.9")

# Pytest is required only to test the Python module
check_for_pytest = False
if env['python_package'] != 'none':
    check_for_pytest = "test" in COMMAND_LINE_TARGETS or any(
        target.startswith("test-python") for target in COMMAND_LINE_TARGETS
    )

# Check for the minimum ruamel.yaml version, 0.15.34, at install and test
# time. The check happens at install and test time because ruamel.yaml is
# only required to run the Python interface, not to build it.
check_for_ruamel_yaml = check_for_pytest or any(
    target in COMMAND_LINE_TARGETS
    for target in ["install", "test"]
)

if env['python_package'] == 'y':
    env['python_package'] = 'full'  # Allow 'y' as a synonym for 'full'
elif env['python_package'] == 'n':
  env['python_package'] = 'none'  # Allow 'n' as a synonym for 'none'

env['install_python_action'] = ''
env['python_module_loc'] = ''
env["python_module"] = None
env["ct_pyscriptdir"] = "<not installed>"

def check_module(name):
    return textwrap.dedent(f"""\
        try:
            import {name}
            versions["{name}"] = {name}.__version__
        except ImportError as {name}_err:
            err += str({name}_err) + "\\n"
    """)

if env['python_package'] != 'none':
    # Test to see if we can import numpy and Cython
    script = textwrap.dedent("""\
        import sys
        import json
        versions = {}
        versions["python"] = "{v.major}.{v.minor}".format(v=sys.version_info)
        err = ""
    """)
    script += check_module("numpy")
    script += check_module("Cython")

    if check_for_ruamel_yaml:
        script += textwrap.dedent("""\
            try:
                from ruamel import yaml
                versions["ruamel.yaml"] = yaml.__version__
            except ImportError as ru_err:
                try:
                    import ruamel_yaml as yaml
                    versions["ruamel.yaml"] = yaml.__version__
                except ImportError as ru_err_2:
                    err += str(ru_err) + "\\n"
                    err += str(ru_err_2) + "\\n"
        """)
    if check_for_pytest:
        script += check_module("pytest")

    script += textwrap.dedent("""\
        print("versions:", json.dumps(versions))
        if err:
            print(err)
    """)

    try:
        info = get_command_output(env["python_cmd"], "-c", script).splitlines()
    except OSError as err:
        logger.debug(f"Error checking for Python:\n{err}")
        warn_no_python = True
    except subprocess.CalledProcessError as err:
        logger.debug(f"Error checking for Python:\n{err} {err.output}")
        warn_no_python = True
    else:
        warn_no_python = False
        for line in info:
            if line.startswith("versions:"):
                versions = {
                    k: parse_version(v)
                    for k, v in json.loads(line.split(maxsplit=1)[1]).items()
                }
                break

        if check_for_ruamel_yaml:
            ruamel_yaml_version = versions.get("ruamel.yaml")
            if not ruamel_yaml_version:
                logger.error(
                    f"ruamel.yaml was not found. {ruamel_min_version} or newer "
                    "is required.")
                sys.exit(1)
            elif ruamel_yaml_version < ruamel_min_version:
                logger.error(
                    "ruamel.yaml is an incompatible version: Found "
                    f"{ruamel_yaml_version}, but {ruamel_min_version} or newer "
                    "is required.")
                sys.exit(1)

    python_version = versions["python"]
    if warn_no_python:
        if env["python_package"] == "default":
            logger.warning(
                "Not building the Python package because the Python interpreter "
                f"{env['python_cmd']!r} could not be found.")
            env["python_package"] = "none"
        else:
            logger.error(
                f"Could not execute the Python interpreter {env['python_cmd']!r}")
            sys.exit(1)
    elif python_version < python_min_version:
        if env["python_package"] in ("minimal", "full", "default"):
            logger.error(
                f"Python version is incompatible. Found {python_version} but "
                f"{python_min_version} or newer is required. In order to install "
                "Cantera without Python support, specify 'python_package=none'.")
            sys.exit(1)
    elif python_version >= python_max_p1_version:
        if env["python_package"] in ("minimal", "full", "default"):
            if env["python_sdist"] or env["package_build"]:
                # An error is triggered as the pip wheel cannot be built for
                # 'scons build' due to safeguards in setup.cfg(.in) files.
                logger.error(
                    f"Python version is incompatible. Found {python_version}, which "
                    "is not yet supported for package builds.")
                sys.exit(1)
            logger.warning(
                f"Python {python_version} is not supported for Cantera "
                f"{env['cantera_version']}. Python versions {python_max_p1_version} and "
                "newer are untested and may result in unexpected behavior. Proceed "
                "with caution.")
    elif env["python_package"] == "minimal":
        # If the minimal package was specified, no further checking needs to be done
        logger.info(f"Building the minimal Python package for Python {python_version}")
    else:

        if len(info) > 1:
            msg = ["Unexpected output while checking Python dependency versions:"]
            msg.extend(line for line in info if not line.startswith("versions:"))
            logger.warning("\n| ".join(msg))

        warn_no_full_package = False
        numpy_version = versions.get("numpy")
        if not numpy_version:
            logger.info("NumPy not found.")
            warn_no_full_package = True
        elif numpy_version < numpy_min_version:
            logger.warning(
                f"NumPy is an incompatible version: Found {numpy_version} but "
                f"{numpy_min_version} or newer is required.")
            warn_no_full_package = True
        else:
            logger.info(f"Using NumPy version {numpy_version}")

        cython_version = versions.get("Cython")
        if not cython_version:
            logger.info("Cython not found.")
            warn_no_full_package = True
        elif cython_version < cython_min_version:
            logger.warning(
                f"Cython is an incompatible version: Found {cython_version} but "
                f"{cython_min_version} or newer is required.")
            warn_no_full_package = True
        elif cython_version < parse_version("3.0.0"):
            logger.info(
                f"Using Cython version {cython_version} (uses legacy NumPy API)")
            env["numpy_1_7_API"] = True
        else:
            logger.info(f"Using Cython version {cython_version}")

        pytest_version = versions.get("pytest")
        if not check_for_pytest:
            pass
        elif not pytest_version:
            logger.error(
                f"pytest was not found. {pytest_min_version} or newer "
                "is required.")
            sys.exit(1)
        elif pytest_version < pytest_min_version:
            logger.error(
                "pytest is an incompatible version: Found "
                f"{pytest_version}, but {pytest_min_version} or newer "
                "is required.")
            sys.exit(1)

        if warn_no_full_package:
            msg = ("Unable to build the full Python package because compatible "
                   "versions of Numpy and/or Cython could not be found.")
            if env["python_package"] == "default":
                logger.warning(msg)
                logger.info(
                    f"Building the minimal Python package for Python {python_version}")
                env["python_package"] = "minimal"
            else:
                logger.error(msg)
                sys.exit(1)
        else:
            logger.info(
                f"Building the full Python package for Python {python_version}")
            env["python_package"] = "full"

if env["python_package"] == "full" and env["OS"] == "Darwin":
    # We need to know the macOS deployment target in advance to be able to determine
    # the name of the wheel file for the Python module. If this is not specified by the
    # MACOSX_DEPLOYMENT_TARGET environment variable, get the value from the Python
    # installation and use that.
    mac_target = env["ENV"].get("MACOSX_DEPLOYMENT_TARGET", None)
    if not mac_target:
        info = get_command_output(
            env["python_cmd"],
            "-c",
            "import sysconfig; print(sysconfig.get_platform())"
        )
        mac_target = info.split("-")[1]
        if parse_version(mac_target) < parse_version('10.15'):
            # macOS 10.15 is the minimum version with C++17 support
            mac_target = '10.15'

    env["ENV"]["MACOSX_DEPLOYMENT_TARGET"] = mac_target
    logger.info(f"MACOSX_DEPLOYMENT_TARGET = {mac_target}")

# Matlab Toolbox settings
if env["matlab_path"] != "" and env["matlab_toolbox"] == "default":
    env["matlab_toolbox"] = "y"

if env["matlab_toolbox"] == "y":
    matlab_path = env["matlab_path"]
    if matlab_path == "":
        logger.error(
            "Unable to build the Matlab toolbox because 'matlab_path' "
            "has not been set.")
        sys.exit(1)

    if env["blas_lapack_libs"] or env["system_blas_lapack"] == "y":
        logger.error(
            "The Matlab toolbox is incompatible with external BLAS and LAPACK "
            "libraries. Unset 'blas_lapack_libs' (for example, 'scons build "
            "blas_lapack_libs=') and/or 'system_blas_lapack' in order to build the "
            "Matlab toolbox, or set 'matlab_toolbox=n' to use the specified BLAS/"
            "LAPACK libraries and skip building the Matlab toolbox.")
        sys.exit(1)

    if env["system_sundials"] == "y":
        logger.error(
            "The Matlab toolbox is incompatible with external SUNDIALS "
            "libraries. Set system_sundials to no (for example, 'scons build "
            "system_sundials=n') in order to build the Matlab "
            "toolbox, or set 'matlab_toolbox=n' to use the specified "
            "SUNDIALS libraries and skip building the Matlab toolbox.")
        sys.exit(1)

    matlab_path = Path(matlab_path)
    env["matlab_path"] = matlab_path.as_posix()
    if not matlab_path.is_dir() and (matlab_path / "extern").is_dir():
        logger.error(
            f"Path set for 'matlab_path' is not correct. Path was {matlab_path!r}")
        sys.exit(1)


# **********************************************
# *** Set additional configuration variables ***
# **********************************************

# Identify options selected either on command line or in cantera.conf
selected_options = set(line.split("=")[0].strip()
    for line in cantera_conf.splitlines())

# Always set the stage directory before building an MSI installer
if "msi" in COMMAND_LINE_TARGETS:
    COMMAND_LINE_TARGETS.append("install")
    env["stage_dir"] = "stage"
    env["prefix"] = "."
    selected_options.add("prefix")
    selected_options.add("stage_dir")
    env["python_package"] = "none"

env["default_prefix"] = True
if "prefix" in selected_options:
    env["default_prefix"] = False

# Use posix-style paths on all platforms
env["prefix"] = Path(env["prefix"]).as_posix()

# Check whether Cantera should be installed into a conda environment
if conda_prefix is not None and sys.executable.startswith(str(conda_prefix)):
    # use conda layout unless any 'blocking' options were specified
    blocking_options = {"layout", "prefix", "python_prefix", "python_cmd"}
    if not selected_options & blocking_options:
        env["layout"] = "conda"
        # Directories where things will be after actually being installed. These
        # variables are the ones that are used to populate header files, scripts,
        # etc.
        conda_prefix = Path(conda_prefix)
        if "stage_dir" in selected_options:
            env["prefix"] = (conda_prefix.relative_to(conda_prefix.parents[2])).as_posix()
        else:
            env["prefix"] = conda_prefix.resolve().as_posix()
        logger.info(
            f"Using conda environment as default 'prefix': {env['prefix']}")
elif env["layout"] == "conda" and not env["package_build"]:
    logger.error("Layout option 'conda' requires a conda environment.")
    sys.exit(1)

prefix = Path(env["prefix"])

if env["layout"] == "conda" and os.name == "nt":
    env["ct_libdir"] = (prefix / "Library" / "lib").as_posix()
    env["ct_shlibdir"] = (prefix / "Library" / "bin").as_posix()
    env["ct_bindir"] = (prefix / "Scripts").as_posix()
    env["ct_python_bindir"] = (prefix / "Scripts").as_posix()
    env["ct_incdir"] = (prefix / "Library" / "include" / "cantera").as_posix()
    env["ct_incroot"] = (prefix / "Library" / "include").as_posix()
else:
    if "stage_dir" not in selected_options:
        prefix = prefix.resolve()
        env["prefix"] = prefix.as_posix()
    env["ct_libdir"] = (prefix / env["libdirname"]).as_posix()
    env["ct_bindir"] = (prefix / "bin").as_posix()

    # On Windows, the search path for DLLs is the "bin" dir. "lib" dirs are used only for
    # static and "import" libraries
    if env["OS"] == "Windows":
        env["ct_shlibdir"] = env["ct_bindir"]
    else:
        env["ct_shlibdir"] = env["ct_libdir"]

    env["ct_python_bindir"] = (prefix / "bin").as_posix()
    env["ct_incdir"] = (prefix / "include" / "cantera").as_posix()
    env["ct_incroot"] = (prefix / "include").as_posix()

env["ct_installroot"] = env["prefix"]

if env["layout"] == "compact":
    env["ct_datadir"] = (prefix / "data").as_posix()
    env["ct_sampledir"] = (prefix / "samples").as_posix()
    env["ct_docdir"] = (prefix / "doc").as_posix()
    env["ct_mandir"] = (prefix / "man1").as_posix()
    env["ct_matlab_dir"] = (prefix / "matlab" / "toolbox").as_posix()
else:
    env["ct_datadir"] = (prefix / "share" / "cantera" / "data").as_posix()
    env["ct_sampledir"] = (prefix  / "share" / "cantera" / "samples").as_posix()
    env["ct_docdir"] = (prefix / "share" / "cantera" / "doc").as_posix()
    env["ct_mandir"] = (prefix / "share" / "man" / "man1").as_posix()
    if env["layout"] == "conda":
        env["ct_matlab_dir"] = (
            prefix / "share" / "cantera" / "matlab" / "toolbox").as_posix()
    else:
        env["ct_matlab_dir"] = (
            prefix / env["libdirname"] / "cantera" / "matlab" / "toolbox").as_posix()


addInstallActions = ('install' in COMMAND_LINE_TARGETS or
                     'uninstall' in COMMAND_LINE_TARGETS)

# Directories where things will be staged for package creation. These
# variables should always be used by the Install(...) targets
if env["stage_dir"]:
    stage_prefix = Path(env["prefix"])
    # Strip the root off the prefix if it's absolute
    if stage_prefix.is_absolute():
        stage_prefix = Path(*stage_prefix.parts[1:])

    instRoot = (Path.cwd().joinpath(env["stage_dir"], stage_prefix)).as_posix()
else:
    instRoot = env["prefix"]

# Prevent setting Cantera installation path to source directory
if os.path.abspath(instRoot) == Dir('.').abspath:
    logger.error("cannot install Cantera into source directory.")
    sys.exit(1)

env["inst_root"] = instRoot
locations = ["libdir", "shlibdir", "bindir", "python_bindir", "incdir", "incroot",
    "matlab_dir", "datadir", "sampledir", "docdir", "mandir"]
for loc in locations:
    if env["prefix"] == ".":
        env[f"inst_{loc}"] = (Path(instRoot) / env[f"ct_{loc}"]).as_posix()
    else:
        env[f"inst_{loc}"] = env[f"ct_{loc}"].replace(env["ct_installroot"], instRoot)

if env['use_rpath_linkage']:
    env.Append(RPATH=env['ct_libdir'])

# **************************************
# *** Set options needed in config.h ***
# **************************************

configh = {}

configh['CANTERA_VERSION'] = quoted(env['cantera_version'])
configh['CANTERA_SHORT_VERSION'] = quoted(env['cantera_short_version'])

# Conditional defines
def cdefine(definevar, configvar, comp=True, value=1):
    if env.get(configvar) == comp:
        configh[definevar] = value
    else:
        configh[definevar] = None

if env['OS'] == 'Solaris' or env['HAS_CLANG']:
    env["RPATHPREFIX"] = "-Wl,-rpath,"

if env["OS"] == "Darwin" and env["use_rpath_linkage"] and not env.subst("$__RPATH"):
    # SCons doesn't want to specify RPATH on macOS, so circumvent that behavior by
    # specifying this directly as part of LINKFLAGS
    env.Append(LINKFLAGS=[env.subst(f'$RPATHPREFIX{x}$RPATHSUFFIX')
                          for x in env['RPATH']])

configh['CT_SUNDIALS_VERSION'] = env['sundials_version'].replace('.','')

if env.get('has_sundials_lapack') and env['use_lapack']:
    configh['CT_SUNDIALS_USE_LAPACK'] = 1
else:
    configh['CT_SUNDIALS_USE_LAPACK'] = 0

cdefine('LAPACK_FTN_STRING_LEN_AT_END', 'lapack_ftn_string_len_at_end')
cdefine('LAPACK_FTN_TRAILING_UNDERSCORE', 'lapack_ftn_trailing_underscore')
cdefine('FTN_TRAILING_UNDERSCORE', 'lapack_ftn_trailing_underscore')
cdefine('CT_USE_LAPACK', 'use_lapack')
cdefine("CT_USE_HDF5", "use_hdf5")
cdefine("CT_USE_SYSTEM_HIGHFIVE", "system_highfive")
cdefine("CT_USE_HIGHFIVE_BOOLEAN", "highfive_boolean")
cdefine("CT_USE_SYSTEM_EIGEN", "system_eigen")
cdefine("CT_USE_SYSTEM_EIGEN_PREFIXED", "system_eigen_prefixed")
cdefine('CT_USE_SYSTEM_FMT', 'system_fmt')
cdefine('CT_USE_SYSTEM_YAMLCPP', 'system_yamlcpp')

config_h_build = env.Command('build/src/config.h.build',
                             'include/cantera/base/config.h.in',
                       ConfigBuilder(configh))
# This separate copy operation, which SCons will skip if config.h.build is
# unmodified, prevents unnecessary rebuilds of the precompiled header
config_h = env.Command('include/cantera/base/config.h',
                       'build/src/config.h.build',
                       Copy('$TARGET', '$SOURCE'))
env.AlwaysBuild(config_h_build)
env['config_h_target'] = config_h

# *********************
# *** Build Cantera ***
# *********************

# Some options to speed up SCons
env.SetOption('max_drift', 2)
env.SetOption('implicit_cache', True)

buildTargets = []
env['build_targets'] = buildTargets
libraryTargets = [] # objects that go in the Cantera library
installTargets = []
sampleTargets = []

def build(targets):
    """ Wrapper to add target to list of build targets """
    buildTargets.extend(targets)
    return targets

def buildSample(*args, **kwargs):
    """ Wrapper to add target to list of samples """
    targets = args[0](*args[1:], **kwargs)
    sampleTargets.extend(targets)
    return targets

def install(*args, **kwargs):
    """ Wrapper to add target to list of install targets """
    if not addInstallActions:
        return
    if len(args) == 2:
        inst = env.Install(*args, **kwargs)
    else:
        inst = args[0](*args[1:], **kwargs)

    installTargets.extend(inst)
    return inst


env.SConsignFile()

env.Prepend(CPPPATH=[],
            LIBPATH=[Dir('build/lib')])

# Add build/lib to library search path in order to find Cantera shared library
if env['OS'] == 'Windows':
    env.PrependENVPath('PATH', Dir('#build/lib').abspath)
elif env['OS'] == 'Darwin':
    env.PrependENVPath('DYLD_LIBRARY_PATH', Dir('#build/lib').abspath)
else:
    env.PrependENVPath('LD_LIBRARY_PATH', Dir('#build/lib').abspath)

for yaml in multi_glob(env, "data", "yaml"):
    dest = Path() / "build" / "data" / yaml.name
    build(env.Command(str(dest), yaml.path, Copy("$TARGET", "$SOURCE")))

if addInstallActions:
    # Put headers in place
    install(env.RecursiveInstall, '$inst_incdir', 'include/cantera')

    # Data files
    install(env.RecursiveInstall, '$inst_datadir', 'build/data')


if env['system_sundials'] == 'y':
    env['sundials_libs'] = ['sundials_cvodes', 'sundials_idas', 'sundials_nvecserial']
    if env['use_lapack']:
        if env.get('has_sundials_lapack'):
            env['sundials_libs'].extend(('sundials_sunlinsollapackdense',
                                         'sundials_sunlinsollapackband'))
        else:
            env['sundials_libs'].extend(('sundials_sunlinsoldense',
                                         'sundials_sunlinsolband'))
else:
    env['sundials_libs'] = []

# List of shared libraries needed to link to Cantera
if env["renamed_shared_libraries"]:
    env["cantera_shared_libs"] = ["cantera_shared"]
else:
    env["cantera_shared_libs"] = ["cantera"]

# External libraries to link to
env["external_libs"] = []
env["external_libs"].extend(env["sundials_libs"])

if env["OS"] == "Linux":
    env["external_libs"].append("dl")

if env["use_hdf5"]:
    env["external_libs"].append("hdf5")

if env["system_fmt"]:
    env["external_libs"].append("fmt")
    # Usually need to link to fmt directly because of templated/inlined code that calls
    # fmt
    env["cantera_shared_libs"].append("fmt")

if env["system_yamlcpp"]:
    env["external_libs"].append("yaml-cpp")

if env["blas_lapack_libs"]:
    env["external_libs"].extend(env["blas_lapack_libs"])

# List of static libraries needed to link to Cantera
env["cantera_libs"] = ["cantera"] + env["external_libs"]

# Add targets from the SConscript files in the various subdirectories
Export('env', 'build', 'libraryTargets', 'install', 'buildSample', "configh")

# ext needs to come before src so that libraryTargets is fully populated
VariantDir('build/ext', 'ext', duplicate=0)
SConscript('build/ext/SConscript')

# Fortran needs to come before src so that libraryTargets is fully populated
if env['f90_interface'] == 'y':
    VariantDir('build/src/fortran/', 'src/fortran', duplicate=1)
    SConscript('build/src/fortran/SConscript')

VariantDir('build/src', 'src', duplicate=0)
SConscript('build/src/SConscript')

if env["python_package"] == "full":
    VariantDir("build/python", "interfaces/cython", duplicate=True)
    SConscript("build/python/SConscript")
elif env["python_package"] == "minimal":
    VariantDir("build/python_minimal", "interfaces/python_minimal", duplicate=True)
    SConscript("build/python_minimal/SConscript")

if env['CC'] != 'cl':
    VariantDir('build/platform', 'platform/posix', duplicate=0)
    SConscript('build/platform/SConscript')

if env['matlab_toolbox'] == 'y':
    SConscript('build/src/matlab/SConscript')

if env['doxygen_docs'] or env['sphinx_docs']:
    SConscript('doc/SConscript')

if env["python_sdist"]:
    VariantDir("build/python_sdist", "interfaces/python_sdist", duplicate=1)
    SConscript("interfaces/python_sdist/SConscript", variant_dir="build/python_sdist")

# Sample programs (also used from test_problems/SConscript)
VariantDir('build/samples', 'samples', duplicate=0)
sampledir_excludes = ['\\.o$', '^~$', '\\.in', 'SConscript']
SConscript('build/samples/cxx/SConscript')
SConscript('build/samples/clib/SConscript')

# Install C++ / C samples
install(env.RecursiveInstall, '$inst_sampledir/cxx',
        'samples/cxx', exclude=sampledir_excludes)
install(env.RecursiveInstall, '$inst_sampledir/clib',
        'samples/clib', exclude=sampledir_excludes)

if env['f90_interface'] == 'y':
    SConscript('build/samples/f77/SConscript')
    SConscript('build/samples/f90/SConscript')

    # install F90 / F77 samples
    install(env.RecursiveInstall, '$inst_sampledir/f77',
            'samples/f77', sampledir_excludes)
    install(env.RecursiveInstall, '$inst_sampledir/f90',
            'samples/f90', sampledir_excludes)

### Meta-targets ###
build_samples = Alias('samples', sampleTargets)

def postBuildMessage(target, source, env):
    build_message = [
        f"\n{' Compilation completed successfully ':*^80}\n",
        "- To run the test suite, type 'scons test'.",
        "- To list available tests, type 'scons test-help'.",
    ]
    if env["googletest"] == "none":
        build_message.append("  WARNING: You set the 'googletest' to 'none' "
            "and all its tests will be skipped.")
    build_message.append("- To install, type 'scons install'.")
    if os.name == 'nt':
        build_message.append("- To create a Windows MSI installer, type 'scons msi'.")
    build_message.append(f"\n{'*' * 80}\n")

    logger.status("\n".join(build_message), print_level=False)

finish_build = env.Command('finish_build', [], postBuildMessage)
env.Depends(finish_build, buildTargets)
build_cantera = Alias('build', finish_build)

Default('build')

def postInstallMessage(target, source, env):
    env_dict = env.Dictionary()
    locations = {
        "library files": "ct_libdir",
        "C++ headers": "ct_incroot",
        "samples": "ct_sampledir",
        "data files": "ct_datadir",
        "input file converters": "ct_pyscriptdir",
    }
    install_message = ["File locations:\n"]
    locations_message = "  {name:<28}{location}"
    for name, location in locations.items():
        install_message.append(locations_message.format(
            name=name, location=env_dict[location]
        ))

    if env["sphinx_docs"] or env["doxygen_docs"]:
        name = "HTML documentation"
        install_message.append(locations_message.format(
            name="HTML documentation", location=env_dict["inst_docdir"]
        ))

    if env["python_package"] == "full":
        env["python_example_loc"] = (Path(env["ct_sampledir"]) / "python").as_posix()
        install_message.append(locations_message.format(
            name="Python package", location=env_dict["python_module_loc"]
        ))
        install_message.append(locations_message.format(
            name="Python examples", location=env_dict["python_example_loc"]
        ))

    if env["matlab_toolbox"] == "y":
        env["matlab_sample_loc"] = (Path(env["ct_sampledir"]) / "matlab").as_posix()
        env["matlab_ctpath_loc"] = (Path(env["ct_matlab_dir"]) / "ctpath.m").as_posix()
        install_message.append(locations_message.format(
            name="Matlab toolbox", location=env_dict["ct_matlab_dir"]
        ))
        install_message.append(locations_message.format(
            name="Matlab samples", location=env_dict["matlab_sample_loc"]
        ))
        install_message.append(textwrap.dedent("""
            An m-file to set the correct matlab path for Cantera is at:

              {matlab_ctpath_loc!s}
        """.format(**env_dict)))

    status = f" Cantera {env['cantera_version']} has been successfully installed "
    install_message = [
        f"\n{status:*^80}\n",
        "\n".join(install_message),
        f"\n{'*' * 80}\n",
    ]
    logger.status("\n".join(install_message), print_level=False)

finish_install = env.Command("finish_install", [], postInstallMessage)
env.Depends(finish_install, installTargets)
install_cantera = Alias('install', finish_install)

### Uninstallation
def getParentDirs(path, top=True):
    head,tail = os.path.split(path)
    if head == os.path.abspath(env['prefix']):
        return [path]
    elif not tail:
        if head.endswith(os.sep):
            return []
        else:
            return [head]
    elif top:
        return getParentDirs(head, False)
    else:
        return getParentDirs(head, False) + [path]

# Files installed by SCons
allfiles = FindInstalledFiles()

# After removing files (which SCons keeps track of),
# remove any empty directories (which SCons doesn't track)
def removeDirectories(target, source, env):
    # Get all directories where files are installed
    alldirs = set()
    for f in allfiles:
        alldirs.update(getParentDirs(f.path))
    if env['layout'] == 'compact':
        alldirs.add(os.path.abspath(env['prefix']))
    # Sort in order of decreasing directory length so that empty subdirectories
    # will be removed before their parents are checked.
    alldirs = sorted(alldirs, key=lambda x: -len(x))

    # Don't remove directories that probably existed before installation,
    # even if they are empty
    keepDirs = ['local/share', 'local/lib', 'local/include', 'local/bin',
                'man/man1', 'dist-packages', 'site-packages']
    for d in alldirs:
        if any(d.endswith(k) for k in keepDirs):
            continue
        if os.path.isdir(d) and not os.listdir(d):
            os.rmdir(d)

uninstall = env.Command("uninstall", None, Delete(allfiles))
env.AddPostAction(uninstall, Action(removeDirectories))
if env["python_package"] == "full":
    env.AddPostAction(uninstall, Action("$python_cmd_esc -m pip uninstall -y Cantera"))
elif env["python_package"] == "minimal":
    env.AddPostAction(uninstall, Action("$python_cmd_esc -m pip uninstall -y Cantera_minimal"))

### Windows MSI Installer ###
if 'msi' in COMMAND_LINE_TARGETS:
    def build_wxs(target, source, env):
        import wxsgen
        wxs = wxsgen.WxsGenerator(env['stage_dir'],
                                  short_version=env['cantera_short_version'],
                                  full_version=env['cantera_pure_version'],
                                  x64=env['TARGET_ARCH']=='amd64',
                                  includeMatlab=env['matlab_toolbox']=='y')
        wxs.make_wxs(str(target[0]))

    wxs_target = env.Command('build/wix/cantera.wxs', [], build_wxs)
    env.AlwaysBuild(wxs_target)

    env.Append(WIXLIGHTFLAGS=['-ext', 'WixUIExtension'])
    msi_target = env.WiX('cantera.msi', ['build/wix/cantera.wxs'])
    env.Depends(wxs_target, installTargets)
    env.Depends(msi_target, wxs_target)
    build_msi = Alias('msi', msi_target)

### Tests ###
if any(target.startswith('test') for target in COMMAND_LINE_TARGETS):
    env['testNames'] = []
    env['test_results'] = env.Command('test_results', [], test_results.print_report)

    if env['python_package'] == 'none':
        # copy scripts from the full Cython module
        # (skipping 'yaml2ck', which depends on the full Python module)
        test_py_int = env.Command('#build/python_local/cantera/__init__.py',
                                  '#interfaces/python_minimal/cantera/__init__.py',
                                  Copy('$TARGET', '$SOURCE'))
        for script in ["ck2yaml", "ctml2yaml", "cti2yaml"]:
            s = env.Command('#build/python_local/cantera/{}.py'.format(script),
                            '#interfaces/cython/cantera/{}.py'.format(script),
                            Copy('$TARGET', '$SOURCE'))
            env.Depends(test_py_int, s)

        env.Depends(env['test_results'], test_py_int)

        env['python_cmd'] = sys.executable
        env.PrependENVPath('PYTHONPATH', Dir('build/python_local').abspath)
    else:
        env.PrependENVPath('PYTHONPATH', Dir('build/python').abspath)

    env['ENV']['PYTHON_CMD'] = env.subst('$python_cmd')

    # Tests written using the gtest framework, the Python unittest module,
    # or the Matlab xunit package.
    VariantDir('build/test', 'test', duplicate=0)
    SConscript('build/test/SConscript')

    # Regression tests
    SConscript('test_problems/SConscript')

    if 'test-help' in COMMAND_LINE_TARGETS:
        print('\n*** Available tests ***\n')
        for name in env['testNames']:
            print('test-%s' % name)
        sys.exit(0)

    Alias('test', env['test_results'])

### Dump (debugging SCons)
if 'dump' in COMMAND_LINE_TARGETS:
    import pprint
    # Typical usage: 'scons build dump'
    print('os.environ:\n', pprint.pprint(dict(os.environ)))
    print('env.Dump():\n', env.Dump())
    sys.exit(0)
