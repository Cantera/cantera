"""
SCons build script for Cantera

Basic usage:
    'scons help' - print a description of user-specifiable options.

    'scons build' - Compile Cantera and the language interfaces using
                    default options.

    'scons clean' - Delete files created while building Cantera.

    '[sudo] scons install' - Install Cantera.

    '[sudo] scons uninstall' - Uninstall Cantera.

    'scons test' - Run all tests which did not previously pass or for which the
                   results may have changed.

    'scons test-reset' - Reset the passing status of all tests.

    'scons test-clean' - Delete files created while running the tests.

    'scons test-help' - List available tests.

    'scons test-NAME' - Run the test named "NAME".

    'scons <command> dump' - Dump the state of the SCons environment to the
                             screen instead of doing <command>, e.g.
                             'scons build dump'. For debugging purposes.

    'scons samples' - Compile the C++ and Fortran samples.

    'scons msi' - Build a Windows installer (.msi) for Cantera.

    'scons sphinx' - Build the Sphinx documentation

    'scons doxygen' - Build the Doxygen documentation
"""

from __future__ import print_function
from buildutils import *

if not COMMAND_LINE_TARGETS:
    # Print usage help
    print(__doc__)
    sys.exit(0)

valid_commands = ('build','clean','install','uninstall',
                  'help','msi','samples','sphinx','doxygen','dump')

for command in COMMAND_LINE_TARGETS:
    if command not in valid_commands and not command.startswith('test'):
        print('ERROR: unrecognized command line target: %r' % command)
        sys.exit(0)

extraEnvArgs = {}

if 'clean' in COMMAND_LINE_TARGETS:
    removeDirectory('build')
    removeDirectory('stage')
    removeDirectory('.sconf_temp')
    removeDirectory('test/work')
    removeFile('.sconsign.dblite')
    removeFile('include/cantera/base/config.h')
    removeFile('src/pch/system.h.gch')
    removeDirectory('include/cantera/ext')
    removeFile('interfaces/cython/cantera/_cantera.cpp')
    removeFile('interfaces/cython/cantera/_cantera.h')
    removeFile('interfaces/cython/setup.py')
    removeFile('interfaces/python_minimal/setup.py')
    removeFile('config.log')
    removeDirectory('doc/sphinx/matlab/examples')
    removeFile('doc/sphinx/matlab/examples.rst')
    for name in os.listdir('doc/sphinx/matlab/'):
        if name.endswith('.rst') and name != 'index.rst':
            removeFile('doc/sphinx/matlab/' + name)
    removeDirectory('doc/sphinx/cython/examples')
    removeFile('doc/sphinx/cython/examples.rst')
    removeDirectory('interfaces/cython/Cantera.egg-info')
    removeDirectory('interfaces/python_minimal/Cantera_minimal.egg-info')
    for name in os.listdir('interfaces/cython/cantera/data/'):
        if name != '__init__.py':
            removeFile('interfaces/cython/cantera/data/' + name)
    removeDirectory('interfaces/cython/cantera/test/data/test_subdir')
    for name in os.listdir('interfaces/cython/cantera/test/data/'):
        if name != '__init__.py':
            removeFile('interfaces/cython/cantera/test/data/' + name)
    for name in os.listdir('.'):
        if name.endswith('.msi'):
            removeFile(name)
    for name in os.listdir('site_scons/'):
        if name.endswith('.pyc'):
            removeFile('site_scons/' + name)
    for name in os.listdir('site_scons/site_tools/'):
        if name.endswith('.pyc'):
            removeFile('site_scons/site_tools/' + name)
    for name in os.listdir('interfaces/python_minimal/cantera'):
        if name != '__init__.py':
            removeFile('interfaces/python_minimal/cantera/' + name)
    removeFile('interfaces/matlab/toolbox/cantera_shared.dll')
    removeFile('interfaces/matlab/Contents.m')
    removeFile('interfaces/matlab/ctpath.m')
    for name in os.listdir('interfaces/matlab/toolbox'):
        if name.startswith('ctmethods.'):
            removeFile('interfaces/matlab/toolbox/' + name)

    print('Done removing output files.')

    if COMMAND_LINE_TARGETS == ['clean']:
        # Just exit if there's nothing else to do
        sys.exit(0)
    else:
        Alias('clean', [])

if 'test-clean' in COMMAND_LINE_TARGETS:
    removeDirectory('build/test')
    removeDirectory('test/work')
    removeDirectory('build/python_local')

# ******************************************************
# *** Set system-dependent defaults for some options ***
# ******************************************************

print('INFO: SCons is using the following Python interpreter:', sys.executable)

opts = Variables('cantera.conf')

windows_compiler_options = []
if os.name == 'nt':
    # On Windows, target the same architecture as the current copy of Python,
    # unless the user specified another option.
    if '64 bit' in sys.version:
        target_arch = 'amd64'
    else:
        target_arch = 'x86'

    # Make an educated guess about the right default compiler
    if which('g++') and not which('cl.exe'):
        defaultToolchain = 'mingw'
    else:
        defaultToolchain = 'msvc'

    windows_compiler_options.extend([
        ('msvc_version',
         """Version of Visual Studio to use. The default is the newest
            installed version. Specify '12.0' for Visual Studio 2013 or '14.0'
            for Visual Studio 2015.""",
         ''),
        EnumVariable(
            'target_arch',
            """Target architecture. The default is the same architecture as the
               installed version of Python.""",
            target_arch, ('amd64', 'x86'))
    ])
    opts.AddVariables(*windows_compiler_options)

    pickCompilerEnv = Environment()
    opts.Update(pickCompilerEnv)

    if pickCompilerEnv['msvc_version']:
        defaultToolchain = 'msvc'

    windows_compiler_options.append(EnumVariable(
        'toolchain',
        """The preferred compiler toolchain.""",
        defaultToolchain, ('msvc', 'mingw', 'intel')))
    opts.AddVariables(windows_compiler_options[-1])
    opts.Update(pickCompilerEnv)

    if pickCompilerEnv['toolchain'] == 'msvc':
        toolchain = ['default']
        if pickCompilerEnv['msvc_version']:
            extraEnvArgs['MSVC_VERSION'] = pickCompilerEnv['msvc_version']
        print('INFO: Compiling with MSVC', (pickCompilerEnv['msvc_version'] or
                                            pickCompilerEnv['MSVC_VERSION']))

    elif pickCompilerEnv['toolchain'] == 'mingw':
        toolchain = ['mingw', 'f90']
        extraEnvArgs['F77'] = None
        # Next line fixes http://scons.tigris.org/issues/show_bug.cgi?id=2683
        extraEnvArgs['WINDOWS_INSERT_DEF'] = 1

    elif pickCompilerEnv['toolchain'] == 'intel':
        toolchain = ['intelc'] # note: untested

    extraEnvArgs['TARGET_ARCH'] = pickCompilerEnv['target_arch']
    print('INFO: Compiling for architecture:', pickCompilerEnv['target_arch'])
    print('INFO: Compiling using the following toolchain(s):', repr(toolchain))

else:
    toolchain = ['default']

env = Environment(tools=toolchain+['textfile', 'subst', 'recursiveInstall', 'wix', 'gch'],
                  ENV={'PATH': os.environ['PATH']},
                  toolchain=toolchain,
                  **extraEnvArgs)

env['OS'] = platform.system()
env['OS_BITS'] = int(platform.architecture()[0][:2])
if 'cygwin' in env['OS'].lower():
    env['OS'] = 'Cygwin' # remove Windows version suffix

# Fixes a linker error in Windows
if os.name == 'nt' and 'TMP' in os.environ:
    env['ENV']['TMP'] = os.environ['TMP']

# Fixes issues with Python subprocesses. See http://bugs.python.org/issue13524
if os.name == 'nt':
    env['ENV']['SystemRoot'] = os.environ['SystemRoot']

# Needed for Matlab to source ~/.matlab7rc.sh
if 'HOME' in os.environ:
    env['ENV']['HOME'] = os.environ['HOME']

# Fix an issue with Unicode sneaking into the environment on Windows
if os.name == 'nt':
    for key,val in env['ENV'].items():
        env['ENV'][key] = str(val)

if 'FRAMEWORKS' not in env:
    env['FRAMEWORKS'] = []

add_RegressionTest(env)

class defaults: pass

if os.name == 'posix':
    defaults.prefix = '/usr/local'
    defaults.boostIncDir = ''
    env['INSTALL_MANPAGES'] = True
elif os.name == 'nt':
    defaults.prefix = pjoin(os.environ['ProgramFiles'], 'Cantera')
    defaults.boostIncDir = ''
    env['INSTALL_MANPAGES'] = False
else:
    print("Error: Unrecognized operating system '%s'" % os.name)
    sys.exit(1)

compiler_options = [
    ('CXX',
     """The C++ compiler to use.""",
     env['CXX']),
    ('CC',
     """The C compiler to use. This is only used to compile CVODE.""",
     env['CC'])]
opts.AddVariables(*compiler_options)
opts.Update(env)

defaults.cxxFlags = ''
defaults.ccFlags = ''
defaults.noOptimizeCcFlags = '-O0'
defaults.optimizeCcFlags = '-O3'
defaults.debugCcFlags = '-g'
defaults.noDebugCcFlags = ''
defaults.debugLinkFlags = ''
defaults.noDebugLinkFlags = ''
defaults.warningFlags = '-Wall'
defaults.buildPch = False
env['pch_flags'] = []
env['openmp_flag'] = ['-fopenmp'] # used to generate sample build scripts

env['using_apple_clang'] = False
# Check if this is actually Apple's clang on macOS
if env['OS'] == 'Darwin':
    result = subprocess.check_output([env.subst('$CC'), '--version']).decode('utf-8')
    if 'clang' in result.lower() and 'Xcode' in result:
        env['using_apple_clang'] = True
        env['openmp_flag'].insert(0, '-Xpreprocessor')

if 'gcc' in env.subst('$CC') or 'gnu-cc' in env.subst('$CC'):
    defaults.optimizeCcFlags += ' -Wno-inline'
    if env['OS'] == 'Cygwin':
        # See http://stackoverflow.com/questions/18784112
        defaults.cxxFlags = '-std=gnu++0x'
    else:
        defaults.cxxFlags = '-std=c++0x'
    defaults.buildPch = True
    env['pch_flags'] = ['-include', 'src/pch/system.h']

elif env['CC'] == 'cl': # Visual Studio
    defaults.cxxFlags = ['/EHsc']
    defaults.ccFlags = ['/MD', '/nologo',
                        '/D_SCL_SECURE_NO_WARNINGS', '/D_CRT_SECURE_NO_WARNINGS']
    defaults.debugCcFlags = '/Zi /Fd${TARGET}.pdb'
    defaults.noOptimizeCcFlags = '/Od /Ob0'
    defaults.optimizeCcFlags = '/O2'
    defaults.debugLinkFlags = '/DEBUG'
    defaults.warningFlags = '/W3'
    defaults.buildPch = True
    env['pch_flags'] = ['/FIpch/system.h']
    env['openmp_flag'] = ['/openmp']

elif 'icc' in env.subst('$CC'):
    defaults.cxxFlags = '-std=c++0x'
    defaults.ccFlags = '-vec-report0 -diag-disable 1478'
    defaults.warningFlags = '-Wcheck'
    env['openmp_flag'] = ['-openmp']

elif 'clang' in env.subst('$CC'):
    defaults.ccFlags = '-fcolor-diagnostics'
    defaults.cxxFlags = '-std=c++11'
    defaults.buildPch = True
    env['pch_flags'] = ['-include-pch', 'src/pch/system.h.gch']

else:
    print("WARNING: Unrecognized C compiler '%s'" % env['CC'])

if env['OS'] in ('Windows', 'Darwin'):
    defaults.threadFlags = ''
else:
    defaults.threadFlags = '-pthread'

# InstallVersionedLib only fully functional in SCons >= 2.4.0
# SHLIBVERSION fails with MinGW: http://scons.tigris.org/issues/show_bug.cgi?id=3035
if (env['toolchain'] == 'mingw'
    or StrictVersion(SCons.__version__) < StrictVersion('2.4.0')):
    defaults.versionedSharedLibrary = False
else:
    defaults.versionedSharedLibrary = True

defaults.fsLayout = 'compact' if env['OS'] == 'Windows' else 'standard'
defaults.env_vars = 'PATH,LD_LIBRARY_PATH,PYTHONPATH'

defaults.python_prefix = '$prefix' if env['OS'] != 'Windows' else ''

# Transform lists into strings to keep cantera.conf clean
for key,value in defaults.__dict__.items():
    if isinstance(value, (list, tuple)):
        setattr(defaults, key, ' '.join(value))

# **************************************
# *** Read user-configurable options ***
# **************************************

config_options = [
    PathVariable(
        'prefix',
        'Set this to the directory where Cantera should be installed.',
        defaults.prefix, PathVariable.PathAccept),
    PathVariable(
        'libdirname',
        """Set this to the directory where Cantera libraries should be installed.
           Some distributions (e.g. Fedora/RHEL) use 'lib64' instead of 'lib' on 64-bit systems
           or could use some other library directory name instead of 'lib' depends
           on architecture and profile (e.g. Gentoo 'libx32' on x32 profile).
           If user didn't set 'libdirname' configuration variable set it to default value 'lib'""",
        'lib', PathVariable.PathAccept),
    EnumVariable(
        'python_package',
        """If you plan to work in Python, then you need the ``full`` Cantera Python
           package. If, on the other hand, you will only use Cantera from some
           other language (e.g. MATLAB or Fortran 90/95) and only need Python
           to process CTI files, then you only need a ``minimal`` subset of the
           package and Cython and NumPy are not necessary. The ``none`` option
           doesn't install any components of the Python interface. The default
           behavior is to build the full Python module for whichever version of
           Python is running SCons if the required prerequisites (NumPy and
           Cython) are installed. Note: ``y`` is a synonym for ``full`` and ``n``
           is a synonym for ``none``.""",
        'default', ('new', 'full', 'minimal', 'none', 'n', 'y', 'default')),
    PathVariable(
        'python_cmd',
        """Cantera needs to know where to find the Python interpreter. If
           PYTHON_CMD is not set, then the configuration process will use the
           same Python interpreter being used by SCons.""",
        sys.executable, PathVariable.PathAccept),
    PathVariable(
        'python_prefix',
        """Use this option if you want to install the Cantera Python package to
           an alternate location. On Unix-like systems, the default is the same
           as the 'prefix' option. If the 'python_prefix' option is set to
           the empty string or the 'prefix' option is not set, then the package
           will be installed to the system default 'site-packages' directory.
           To install to the current user's 'site-packages' directory, use
           'python_prefix=USER'.""",
        defaults.python_prefix, PathVariable.PathAccept),
    EnumVariable(
        'python3_package',
        """Deprecated synonym for the 'python_package' option. Will be overridden
        if 'python_package' is set.""",
        'default', ('full', 'minimal', 'y', 'n', 'none', 'default')),
    PathVariable(
        'python3_cmd',
        """Deprecated synonym for the 'python_cmd' option. Will be overridden
           if 'python_cmd' is set.""",
        sys.executable, PathVariable.PathAccept),
    PathVariable(
        'python3_prefix',
        """Deprecated synonym for the 'python_prefix' option. Will be overridden
           if 'python_prefix' is set.""",
        defaults.python_prefix, PathVariable.PathAccept),
     EnumVariable(
        'matlab_toolbox',
        """This variable controls whether the MATLAB toolbox will be built. If
           set to 'y', you will also need to set the value of the 'matlab_path'
           variable. If set to 'default', the MATLAB toolbox will be built if
           'matlab_path' is set.""",
        'default', ('y', 'n', 'default')),
    PathVariable(
        'matlab_path',
        """Path to the MATLAB install directory. This should be the directory
           containing the 'extern', 'bin', etc. subdirectories. Typical values
           are: "C:/Program Files/MATLAB/R2011a" on Windows,
           "/Applications/MATLAB_R2011a.app" on OS X, or
           "/opt/MATLAB/R2011a" on Linux.""",
        '', PathVariable.PathAccept),
    EnumVariable(
        'f90_interface',
        """This variable controls whether the Fortran 90/95 interface will be
           built. If set to 'default', the builder will look for a compatible
           Fortran compiler in the 'PATH' environment variable, and compile
           the Fortran 90 interface if one is found.""",
        'default', ('y', 'n', 'default')),
    PathVariable(
        'FORTRAN',
        """The Fortran (90) compiler. If unspecified, the builder will look for
           a compatible compiler (pgfortran, gfortran, ifort, g95) in the 'PATH' environment
           variable. Used only for compiling the Fortran 90 interface.""",
        '', PathVariable.PathAccept),
    ('FORTRANFLAGS',
     'Compilation options for the Fortran (90) compiler.',
     '-O3'),
    BoolVariable(
        'coverage',
        """Enable collection of code coverage information with gcov.
           Available only when compiling with gcc.""",
        False),
    BoolVariable(
        'doxygen_docs',
        """Build HTML documentation for the C++ interface using Doxygen.""",
        False),
    BoolVariable(
        'sphinx_docs',
        """Build HTML documentation for Cantera using Sphinx.""",
        False),
    PathVariable(
        'sphinx_cmd',
        """Command to use for building the Sphinx documentation.""",
        'sphinx-build', PathVariable.PathAccept),
    EnumVariable(
        'system_eigen',
        """Select whether to use Eigen from a system installation ('y'), from a
           Git submodule ('n'), or to decide automatically ('default'). If Eigen
           is not installed directly into a system include directory, e.g. it is
           installed in '/opt/include/eigen3/Eigen', then you will need to add
           '/opt/include/eigen3' to 'extra_inc_dirs'.
           """,
        'default', ('default', 'y', 'n')),
    EnumVariable(
        'system_fmt',
        """Select whether to use the fmt library from a system installation
           ('y'), from a Git submodule ('n'), or to decide automatically
           ('default').""",
        'default', ('default', 'y', 'n')),
    EnumVariable(
        'system_yamlcpp',
        """Select whether to use the yaml-cpp library from a system installation
           ('y'), from a Git submodule ('n'), or to decide automatically
           ('default'). If yaml-cpp is not installed directly into system
           include and library directories, then you will need to add those
           directories to 'extra_inc_dirs' and 'extra_lib_dirs'.""",
        'default', ('default', 'y', 'n')),
    EnumVariable(
        'system_sundials',
        """Select whether to use SUNDIALS from a system installation ('y'), from
           a Git submodule ('n'), or to decide automatically ('default').
           Specifying 'sundials_include' or 'sundials_libdir' changes the
           default to 'y'.""",
        'default', ('default', 'y', 'n')),
    PathVariable(
        'sundials_include',
        """The directory where the SUNDIALS header files are installed. This
           should be the directory that contains the "cvodes", "nvector", etc.
           subdirectories. Not needed if the headers are installed in a
           standard location, e.g., '/usr/include'.""",
        '', PathVariable.PathAccept),
    PathVariable(
        'sundials_libdir',
        """The directory where the SUNDIALS static libraries are installed.
           Not needed if the libraries are installed in a standard location,
           e.g., '/usr/lib'.""",
        '', PathVariable.PathAccept),
    (
        'blas_lapack_libs',
        """Cantera can use BLAS and LAPACK libraries available on your system if
           you have optimized versions available (e.g., Intel MKL). Otherwise,
           Cantera will use Eigen for linear algebra support. To use BLAS
           and LAPACK, set 'blas_lapack_libs' to the the list of libraries
           that should be passed to the linker, separated by commas, e.g.,
           "lapack,blas" or "lapack,f77blas,cblas,atlas". Eigen is required
           whether or not BLAS/LAPACK are used.""",
        ''),
    PathVariable(
        'blas_lapack_dir',
        """Directory containing the libraries specified by 'blas_lapack_libs'. Not
           needed if the libraries are installed in a standard location, e.g.
           ``/usr/lib``.""",
        '', PathVariable.PathAccept),
    EnumVariable(
        'lapack_names',
        """Set depending on whether the procedure names in the specified
           libraries are lowercase or uppercase. If you don't know, run 'nm' on
           the library file (e.g., 'nm libblas.a').""",
        'lower', ('lower','upper')),
    BoolVariable(
        'lapack_ftn_trailing_underscore',
        """Controls whether the LAPACK functions have a trailing underscore
           in the Fortran libraries.""",
        True),
    BoolVariable(
        'lapack_ftn_string_len_at_end',
        """Controls whether the LAPACK functions have the string length
           argument at the end of the argument list ('yes') or after
           each argument ('no') in the Fortran libraries.""",
        True),
    EnumVariable(
        'system_googletest',
        """Select whether to use gtest/gmock from system
           installation ('y'), from a Git submodule ('n'), or to decide
           automatically ('default'). Deprecated option, please use 'googletest' instead.
           This option is suppressed by 'googletest' option.""",
        'unspecified', ('unspecified','default', 'y', 'n')),
    EnumVariable(
        'googletest',
        """Select whether to use gtest/gmock from system
           installation ('system'), from a Git submodule ('submodule'), to decide
           automatically ('default') or don't look for gtest/gmock ('none')
           and don't run tests that depend on gtest/gmock. If this option is
           set then it suppresses the deprecated 'system_googletest' option.""",
        'default', ('default', 'system', 'submodule', 'none')),
    (
        'env_vars',
        """Environment variables to propagate through to SCons. Either the
           string "all" or a comma separated list of variable names, e.g.
           'LD_LIBRARY_PATH,HOME'.""",
        defaults.env_vars),
    BoolVariable(
        'use_pch',
        """Use a precompiled-header to speed up compilation""",
        defaults.buildPch),
    (
        'cxx_flags',
        """Compiler flags passed to the C++ compiler only. Separate multiple
           options with spaces, e.g., "cxx_flags='-g -Wextra -O3 --std=c++11'"
           """,
        defaults.cxxFlags),
    (
        'cc_flags',
        """Compiler flags passed to both the C and C++ compilers, regardless of optimization level.""",
        defaults.ccFlags),
    (
        'thread_flags',
        """Compiler and linker flags for POSIX multithreading support.""",
        defaults.threadFlags),
    BoolVariable(
        'optimize',
        """Enable extra compiler optimizations specified by the
           'optimize_flags' variable, instead of the flags specified by the
           'no_optimize_flags' variable.""",
        True),
    (
        'optimize_flags',
        """Additional compiler flags passed to the C/C++ compiler when 'optimize=yes'.""",
        defaults.optimizeCcFlags),
    (
        'no_optimize_flags',
        """Additional compiler flags passed to the C/C++ compiler when 'optimize=no'.""",
        defaults.noOptimizeCcFlags),
    BoolVariable(
        'debug',
        """Enable compiler debugging symbols.""",
        True),
    (
        'debug_flags',
        """Additional compiler flags passed to the C/C++ compiler when 'debug=yes'.""",
        defaults.debugCcFlags),
    (
        'no_debug_flags',
        """Additional compiler flags passed to the C/C++ compiler when 'debug=no'.""",
        defaults.noDebugCcFlags),
    (
        'debug_linker_flags',
        """Additional options passed to the linker when 'debug=yes'.""",
        defaults.debugLinkFlags),
    (
        'no_debug_linker_flags',
        """Additional options passed to the linker when 'debug=no'.""",
        defaults.noDebugLinkFlags),
    (
        'warning_flags',
        """Additional compiler flags passed to the C/C++ compiler to enable
           extra warnings. Used only when compiling source code that is part
           of Cantera (e.g. excluding code in the 'ext' directory).""",
        defaults.warningFlags),
    (
        'extra_inc_dirs',
        """Additional directories to search for header files, with multiple
        directories separated by colons (*nix, macOS) or semicolons (Windows)""",
        ''),
    (
        'extra_lib_dirs',
        """Additional directories to search for libraries, with multiple
        directories separated by colons (*nix, macOS) or semicolons (Windows)""",
        ''),
    PathVariable(
        'boost_inc_dir',
        """Location of the Boost header files. Not needed if the headers are
           installed in a standard location, e.g. '/usr/include'.""",
        defaults.boostIncDir, PathVariable.PathAccept),
    PathVariable(
        'stage_dir',
        """Directory relative to the Cantera source directory to be
           used as a staging area for building e.g. a Debian
           package. If specified, 'scons install' will install files
           to 'stage_dir/prefix/...'.""",
        '',
        PathVariable.PathAccept),
    BoolVariable(
        'VERBOSE',
        """Create verbose output about what SCons is doing.""",
        False),
    (
        'gtest_flags',
        """Additional options passed to each GTest test suite, e.g.
           '--gtest_filter=*pattern*'. Separate multiple options with spaces.""",
        ''),
    BoolVariable(
        'renamed_shared_libraries',
        """If this option is turned on, the shared libraries that are created
           will be renamed to have a '_shared' extension added to their base name.
           If not, the base names will be the same as the static libraries.
           In some cases this simplifies subsequent linking environments with
           static libraries and avoids a bug with using valgrind with
           the '-static' linking flag.""",
        True),
    BoolVariable(
        'versioned_shared_library',
        """If enabled, create a versioned shared library, with symlinks to the
           more generic library name, e.g. 'libcantera_shared.so.2.5.0' as the
           actual library and 'libcantera_shared.so' and 'libcantera_shared.so.2'
           as symlinks.
           """,
        defaults.versionedSharedLibrary),
    BoolVariable(
        'use_rpath_linkage',
        """If enabled, link to all shared libraries using 'rpath', i.e., a fixed
           run-time search path for dynamic library loading.""",
        True),
    EnumVariable(
        'layout',
        """The layout of the directory structure. 'standard' installs files to
           several subdirectories under 'prefix', e.g. 'prefix/bin',
           'prefix/include/cantera', 'prefix/lib' etc. This layout is best used in
           conjunction with "prefix'='/usr/local'". 'compact' puts all installed
           files in the subdirectory defined by 'prefix'. This layout is best
           with a prefix like '/opt/cantera'. 'debian' installs to the stage
           directory in a layout used for generating Debian packages.""",
        defaults.fsLayout, ('standard','compact','debian')),
    BoolVariable(
        "fast_fail_tests",
        """If enabled, tests will exit at the first failure.""",
        False,
    )
]

opts.AddVariables(*config_options)
opts.Update(env)
opts.Save('cantera.conf', env)

# Expand ~/ and environment variables used in cantera.conf (variables used on
# the command line will be expanded by the shell)
for option in opts.keys():
    original = env[option]
    if isinstance(original, str):
        modified = os.path.expandvars(os.path.expanduser(env[option]))
        if original != modified:
            print('INFO: Expanding {!r} to {!r}'.format(original, modified))
            env[option] = modified

if 'help' in COMMAND_LINE_TARGETS:
    ### Print help about configuration options and exit.
    print("""
        **************************************************
        *   Configuration options for building Cantera   *
        **************************************************

The following options can be passed to SCons to customize the Cantera
build process. They should be given in the form:

    scons build option1=value1 option2=value2

Variables set in this way will be stored in the 'cantera.conf' file and reused
automatically on subsequent invocations of scons. Alternatively, the
configuration options can be entered directly into 'cantera.conf' before
running 'scons build'. The format of this file is:

    option1 = 'value1'
    option2 = 'value2'

        **************************************************
""")

    for opt in opts.options:
        print('\n'.join(formatOption(env, opt)))
    sys.exit(0)

if 'doxygen' in COMMAND_LINE_TARGETS:
    env['doxygen_docs'] = True
if 'sphinx' in COMMAND_LINE_TARGETS:
    env['sphinx_docs'] = True

valid_arguments = (set(opt[0] for opt in windows_compiler_options) |
                   set(opt[0] for opt in compiler_options) |
                   set(opt[0] for opt in config_options))
for arg in ARGUMENTS:
    if arg not in valid_arguments:
        print('Encountered unexpected command line argument: %r' % arg)
        sys.exit(1)

# Require a StrictVersion-compatible version
env['cantera_version'] = "2.5.0a4"
ctversion = StrictVersion(env['cantera_version'])
# For use where pre-release tags are not permitted (MSI, sonames)
env['cantera_pure_version'] = '.'.join(str(x) for x in ctversion.version)
env['cantera_short_version'] = '.'.join(str(x) for x in ctversion.version[:2])

try:
    env['git_commit'] = getCommandOutput('git', 'rev-parse', '--short', 'HEAD')
except Exception:
    env['git_commit'] = 'unknown'

# Print values of all build options:
print("Configuration variables read from 'cantera.conf' and command line:")
for line in open('cantera.conf'):
    print('   ', line.strip())
print()

# ********************************************
# *** Configure system-specific properties ***
# ********************************************

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
            if env['VERBOSE']:
                print('Propagating environment variable {0}={1}'.format(name, env['ENV'][name]))
        elif name not in defaults.env_vars.split(','):
            print('WARNING: failed to propagate environment variable', repr(name))
            print('         Edit cantera.conf or the build command line to fix this.')

# @todo: Remove this Warning after Cantera 2.5
if os.pathsep == ';':
    for dirs in (env['extra_inc_dirs'], env['extra_lib_dirs']):
        if re.search(r':\w:', dirs):
            print('ERROR: Multiple entries in "extra_inc_dirs" and "extra_lib_dirs" '
                  'should be separated by semicolons (;) on Windows. Use of OS-specific '
                  'path separator introduced in Cantera 2.5.')
            sys.exit(1)

env['extra_inc_dirs'] = [d for d in env['extra_inc_dirs'].split(os.pathsep) if d]
env['extra_lib_dirs'] = [d for d in env['extra_lib_dirs'].split(os.pathsep) if d]

env.Append(CPPPATH=env['extra_inc_dirs'],
           LIBPATH=env['extra_lib_dirs'])

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
    env.Append(CPPPATH=env['boost_inc_dir'])

if env['blas_lapack_dir']:
    env.Append(LIBPATH=[env['blas_lapack_dir']])
    if env['use_rpath_linkage']:
        env.Append(RPATH=env['blas_lapack_dir'])

if env['system_sundials'] in ('y','default'):
    if env['sundials_include']:
        env.Append(CPPPATH=[env['sundials_include']])
        env['system_sundials'] = 'y'
    if env['sundials_libdir']:
        env.Append(LIBPATH=[env['sundials_libdir']])
        env['system_sundials'] = 'y'
        if env['use_rpath_linkage']:
            env.Append(RPATH=env['sundials_libdir'])

# BLAS / LAPACK configuration
if env['blas_lapack_libs'] != '':
    env['blas_lapack_libs'] = env['blas_lapack_libs'].split(',')
    env['use_lapack'] = True
elif env['OS'] == 'Darwin':
    env['blas_lapack_libs'] = []
    env['use_lapack'] = True
    env.Append(FRAMEWORKS=['Accelerate'])
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

conf = Configure(env, custom_tests={'CheckStatement': CheckStatement})

# Set up compiler options before running configuration tests
env['CXXFLAGS'] = listify(env['cxx_flags'])
env['CCFLAGS'] = listify(env['cc_flags']) + listify(env['thread_flags'])
env['LINKFLAGS'] += listify(env['thread_flags'])
env['CPPDEFINES'] = {}

env['warning_flags'] = listify(env['warning_flags'])

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
        print('Error: coverage testing is only available with GCC.')
        exit(0)

if env['toolchain'] == 'mingw':
    env.Append(LINKFLAGS=['-static-libgcc', '-static-libstdc++'])

def config_error(message):
    print('ERROR:', message)
    if env['VERBOSE']:
        print('*' * 25, 'Contents of config.log:', '*' * 25)
        print(open('config.log').read())
        print('*' * 28, 'End of config.log', '*' * 28)
    else:
        print("See 'config.log' for details.")
    sys.exit(1)

# First, a sanity check:
if not conf.CheckCXXHeader('cmath', '<>'):
    config_error('The C++ compiler is not correctly configured.')


def get_expression_value(includes, expression):
    s = ['#include ' + i for i in includes]
    s.extend(('#define Q(x) #x',
              '#define QUOTE(x) Q(x)',
              '#include <iostream>',
              '#ifndef SUNDIALS_PACKAGE_VERSION', # name change in Sundials >= 3.0
              '#define SUNDIALS_PACKAGE_VERSION SUNDIALS_VERSION',
              '#endif',
              'int main(int argc, char** argv) {',
              '    std::cout << %s << std::endl;' % expression,
              '    return 0;',
              '}\n'))
    return '\n'.join(s)

# Check for fmt library and checkout submodule if needed
# Test for 'ostream.h' to ensure that version >= 3.0.0 is available
if env['system_fmt'] in ('y', 'default'):
    if conf.CheckCXXHeader('fmt/ostream.h', '""'):
        env['system_fmt'] = True
        print("""INFO: Using system installation of fmt library.""")

    elif env['system_fmt'] == 'y':
        config_error('Expected system installation of fmt library, but it '
            'could not be found.')

if env['system_fmt'] in ('n', 'default'):
    env['system_fmt'] = False
    print("""INFO: Using private installation of fmt library.""")
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

fmt_include = '<fmt/format.h>' if env['system_fmt'] else '"../ext/fmt/include/fmt/format.h"'
fmt_version_source = get_expression_value([fmt_include], 'FMT_VERSION')
retcode, fmt_lib_version = conf.TryRun(fmt_version_source, '.cpp')
try:
    fmt_lib_version = divmod(float(fmt_lib_version.strip()), 10000)
    (fmt_maj, (fmt_min, fmt_pat)) = fmt_lib_version[0], divmod(fmt_lib_version[1], 100)
    env['FMT_VERSION'] = '{major:.0f}.{minor:.0f}.{patch:.0f}'.format(major=fmt_maj, minor=fmt_min, patch=fmt_pat)
    print('INFO: Found fmt version {}'.format(env['FMT_VERSION']))
except ValueError:
    env['FMT_VERSION'] = '0.0.0'
    print('INFO: Could not find version of fmt')

# Convert setting of the deprecated system_googletest option
if env['system_googletest'] != 'unspecified':
    print("WARNING: The 'system_googletest' option is deprecated. "
        "Use the 'googletest' option instead.")
if env['system_googletest'] == 'default':
    env['googletest'] = 'default'
elif env['system_googletest'] == 'y':
    env['googletest'] = 'system'
elif env['system_googletest'] == 'n':
    env['googletest'] = 'submodule'

# Check for yaml-cpp library and checkout submodule if needed
if env['system_yamlcpp'] in ('y', 'default'):
    if conf.CheckCXXHeader('yaml-cpp/yaml.h', '""'):
        env['system_yamlcpp'] = True
        print("""INFO: Using system installation of yaml-cpp library.""")

    elif env['system_yamlcpp'] == 'y':
        config_error('Expected system installation of yaml-cpp library, but it '
            'could not be found.')

if env['system_yamlcpp'] in ('n', 'default'):
    env['system_yamlcpp'] = False
    print("""INFO: Using private installation of yaml-cpp library.""")
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
        print("""INFO: Using system installation of Googletest""")
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
    print("""INFO: Using Googletest from Git submodule""")

if env['googletest'] == 'none':
    print("""INFO: Not using Googletest -- unable to run complete test suite""")

# Check for Eigen and checkout submodule if needed
if env['system_eigen'] in ('y', 'default'):
    if conf.CheckCXXHeader('Eigen/Dense', '<>'):
        env['system_eigen'] = True
        print("""INFO: Using system installation of Eigen.""")
    elif env['system_eigen'] == 'y':
        config_error('Expected system installation of Eigen, but it '
                     'could not be found.')

if env['system_eigen'] in ('n', 'default'):
    env['system_eigen'] = False
    print("""INFO: Using private installation of Eigen.""")
    if not os.path.exists('ext/eigen/Eigen/Dense'):
        if not os.path.exists('.git'):
            config_error('Eigen is missing. Install Eigen in ext/eigen.')

        try:
            code = subprocess.call(['git','submodule','update','--init',
                                    '--recursive','ext/eigen'])
        except Exception:
            code = -1
        if code:
            config_error('Eigen not found and submodule checkout failed.\n'
                         'Try manually checking out the submodule with:\n\n'
                         '    git submodule update --init --recursive ext/eigen\n')

eigen_include = '<Eigen/Core>' if env['system_eigen'] else '"../ext/eigen/Eigen/Core"'
eigen_versions = 'QUOTE(EIGEN_WORLD_VERSION) "." QUOTE(EIGEN_MAJOR_VERSION) "." QUOTE(EIGEN_MINOR_VERSION)'
eigen_version_source = get_expression_value([eigen_include], eigen_versions)
retcode, eigen_lib_version = conf.TryRun(eigen_version_source, '.cpp')
env['EIGEN_LIB_VERSION'] = eigen_lib_version.strip()
print('INFO: Found Eigen version {}'.format(env['EIGEN_LIB_VERSION']))

# Determine which standard library to link to when using Fortran to
# compile code that links to Cantera
if conf.CheckDeclaration('__GLIBCXX__', '#include <iostream>', 'C++'):
    env['cxx_stdlib'] = ['stdc++']
elif conf.CheckDeclaration('_LIBCPP_VERSION', '#include <iostream>', 'C++'):
    env['cxx_stdlib'] = ['c++']
else:
    env['cxx_stdlib'] = []

env['HAS_CLANG'] = conf.CheckDeclaration('__clang__', '', 'C++')

boost_version_source = get_expression_value(['<boost/version.hpp>'], 'BOOST_LIB_VERSION')
retcode, boost_lib_version = conf.TryRun(boost_version_source, '.cpp')
env['BOOST_LIB_VERSION'] = '.'.join(boost_lib_version.strip().split('_'))
print('INFO: Found Boost version {0}'.format(env['BOOST_LIB_VERSION']))
if not env['BOOST_LIB_VERSION']:
    config_error("Boost could not be found. Install Boost headers or set"
                 " 'boost_inc_dir' to point to the boost headers.")
# demangle is availble in Boost 1.55 or newer
env['has_demangle'] = conf.CheckDeclaration("boost::core::demangle",
                                '#include <boost/core/demangle.hpp>', 'C++')

import SCons.Conftest, SCons.SConf
context = SCons.SConf.CheckContext(conf)

# Check initially for Sundials<=3.2 and then for Sundials>=4.0
for cvode_call in ['CVodeCreate(CV_BDF, CV_NEWTON);','CVodeCreate(CV_BDF);']:
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
    for subdir in ('sundials', 'nvector', 'cvodes', 'ida', 'sunlinsol', 'sunmatrix'):
        removeDirectory('include/cantera/ext/'+subdir)

    # Determine Sundials version
    sundials_version_source = get_expression_value(['"sundials/sundials_config.h"'],
                                                   'QUOTE(SUNDIALS_PACKAGE_VERSION)')
    retcode, sundials_version = conf.TryRun(sundials_version_source, '.cpp')
    if retcode == 0:
        config_error("Failed to determine Sundials version.")
    sundials_version = sundials_version.strip(' "\n')

    # Ignore the minor version, e.g. 2.4.x -> 2.4
    env['sundials_version'] = '.'.join(sundials_version.split('.')[:2])
    sundials_ver = LooseVersion(env['sundials_version'])
    if sundials_ver < LooseVersion('2.4') or sundials_ver >= LooseVersion('6.0'):
        print("""ERROR: Sundials version %r is not supported.""" % env['sundials_version'])
        sys.exit(1)
    elif sundials_ver > LooseVersion('5.3'):
        print("WARNING: Sundials version %r has not been tested." % env['sundials_version'])

    print("""INFO: Using system installation of Sundials version %s.""" % sundials_version)

    #Determine whether or not Sundials was built with BLAS/LAPACK
    if sundials_ver < LooseVersion('2.6'):
        # In Sundials 2.4 / 2.5, SUNDIALS_BLAS_LAPACK is either 0 or 1
        sundials_blas_lapack = get_expression_value(['"sundials/sundials_config.h"'],
                                                       'SUNDIALS_BLAS_LAPACK')
        retcode, has_sundials_lapack = conf.TryRun(sundials_blas_lapack, '.cpp')
        if retcode == 0:
            config_error("Failed to determine Sundials BLAS/LAPACK.")
        env['has_sundials_lapack'] = int(has_sundials_lapack.strip())
    else:
        # In Sundials 2.6, SUNDIALS_BLAS_LAPACK is either defined or undefined
        env['has_sundials_lapack'] = conf.CheckDeclaration('SUNDIALS_BLAS_LAPACK',
                '#include "sundials/sundials_config.h"', 'C++')

    # In the case where a user is trying to link Cantera to an external BLAS/LAPACK
    # library, but Sundials was configured without this support, print a Warning.
    if not env['has_sundials_lapack'] and env['use_lapack']:
        print('WARNING: External BLAS/LAPACK has been specified for Cantera '
              'but Sundials was built without this support.')
else: # env['system_sundials'] == 'n'
    print("""INFO: Using private installation of Sundials version 5.3.""")
    env['sundials_version'] = '5.3'
    env['has_sundials_lapack'] = int(env['use_lapack'])


# Try to find a working Fortran compiler:
def check_fortran(compiler, expected=False):
    hello_world = '''
program main
   write(*,'(a)') 'Hello, world!'
end program main
    '''
    if which(compiler) is not None:
        env['F77'] = env['F90'] = env['F95'] = env['F03'] = env['FORTRAN'] = compiler
        success, output = conf.TryRun(hello_world, '.f90')
        if success and 'Hello, world!' in output:
            return True
        else:
            print("WARNING: Unable to use '%s' to compile the Fortran "
                  "interface. See config.log for details." % compiler)
            return False
    elif expected:
        print("ERROR: Couldn't find specified Fortran compiler: '%s'" % compiler)
        sys.exit(1)

    return False

env['F77FLAGS'] = env['F90FLAGS'] = env['F95FLAGS'] = env['F03FLAGS'] = env['FORTRANFLAGS']

if env['f90_interface'] in ('y','default'):
    foundF90 = False
    if env['FORTRAN']:
        foundF90 = check_fortran(env['FORTRAN'], True)

    for compiler in ('pgfortran', 'gfortran', 'ifort', 'g95'):
        if foundF90:
            break
        foundF90 = check_fortran(compiler)

    if foundF90:
        print("INFO: Using '%s' to build the Fortran 90 interface" % env['FORTRAN'])
        env['f90_interface'] = 'y'
    else:
        if env['f90_interface'] == 'y':
            print("ERROR: Couldn't find a suitable Fortran compiler to build the Fortran 90 interface.")
            sys.exit(1)
        else:
            env['f90_interface'] = 'n'
            env['FORTRAN'] = ''
            print("INFO: Skipping compilation of the Fortran 90 interface.")

if 'pgfortran' in env['FORTRAN']:
    env['FORTRANMODDIRPREFIX'] = '-module '
elif 'gfortran' in env['FORTRAN']:
    env['FORTRANMODDIRPREFIX'] = '-J'
elif 'g95' in env['FORTRAN']:
    env['FORTRANMODDIRPREFIX'] = '-fmod='
elif 'ifort' in env['FORTRAN']:
    env['FORTRANMODDIRPREFIX'] = '-module '

env['F77'] = env['F90'] = env['F95'] = env['F03'] = env['FORTRAN']

env['FORTRANMODDIR'] = '${TARGET.dir}'

env = conf.Finish()

if env['VERBOSE']:
    print('-------------------- begin config.log --------------------')
    print(open('config.log').read())
    print('--------------------- end config.log ---------------------')

env['python_cmd_esc'] = quoted(env['python_cmd'])

# Python Package Settings
python_min_version = LooseVersion('3.5')
# Note: cython_min_version is redefined below if the Python version is 3.8 or higher
cython_min_version = LooseVersion('0.23')
numpy_min_test_version = LooseVersion('1.8.1')

# Handle the version-specific Python package options
python_options = (('package', 'default'),
                  ('cmd', sys.executable),
                  ('prefix', defaults.python_prefix))
for option, default in python_options:
    if env['python3_' + option] != default:
        print("WARNING: Option 'python3_{0}' is deprecated."
            " Use 'python_{0}' instead.".format(option))
        if env['python_' + option] == default:
            env['python_' + option] = env['python3_' + option]
        else:
            print("WARNING: Ignoring option 'python3_{0}'"
                " because 'python_{0}' was also set.".format(option))

if env['python_package'] == 'new':
    print("WARNING: The 'new' option for the Python package is "
          "deprecated and will be removed in the future. Use "
          "'full' instead.")
    env['python_package'] = 'full'  # Allow 'new' as a synonym for 'full'
elif env['python_package'] == 'y':
    env['python_package'] = 'full'  # Allow 'y' as a synonym for 'full'
elif env['python_package'] == 'n':
  env['python_package'] = 'none'  # Allow 'n' as a synonym for 'none'

env['install_python_action'] = ''
env['python_module_loc'] = ''

if env['python_package'] != 'none':
    # Test to see if we can import numpy and Cython
    script = textwrap.dedent("""\
        import sys
        print('{v.major}.{v.minor}'.format(v=sys.version_info))
        err = ''
        try:
            import numpy
            print(numpy.__version__)
        except ImportError as np_err:
            print('0.0.0')
            err += str(np_err) + '\\n'
        try:
            import Cython
            print(Cython.__version__)
        except ImportError as cython_err:
            print('0.0.0')
            err += str(cython_err) + '\\n'
        if err:
            print(err)
    """)

    try:
        info = getCommandOutput(env['python_cmd'], '-c', script).splitlines()
    except OSError as err:
        if env['VERBOSE']:
            print('Error checking for Python:')
            print(err)
        warn_no_python = True
    except subprocess.CalledProcessError as err:
        if env['VERBOSE']:
            print('Error checking for Python:')
            print(err, err.output)
        warn_no_python = True
    else:
        warn_no_python = False
        python_version = LooseVersion(info[0])
        numpy_version = LooseVersion(info[1])
        cython_version = LooseVersion(info[2])

    if warn_no_python:
        if env['python_package'] == 'default':
            print('WARNING: Not building the Python package because the Python '
                  'interpreter {!r} could not be found'.format(env['python_cmd']))
            env['python_package'] = 'none'
        else:
            print('ERROR: Could not execute the Python interpreter {!r}'.format(
                env['python_cmd']))
            sys.exit(1)
    elif env['python_package'] == 'minimal':
        print('INFO: Building the minimal Python package for Python {}'.format(python_version))
    else:
        warn_no_full_package = False
        if len(info) > 3:
            print("WARNING: Unexpected output while checking Python / Numpy / Cython versions:")
            print('| ' + '\n| '.join(info[3:]))

        if python_version < python_min_version:
            print("WARNING: Python version is incompatible with the full Python module: "
                "Found {0} but {1} or newer is required".format(
                python_version, python_min_version))
            warn_no_full_package = True
        elif python_version >= LooseVersion("3.8"):
            # Reset the minimum Cython version if the Python version is 3.8 or higher
            # Due to internal changes in the CPython API, more recent versions of
            # Cython are necessary to build for Python 3.8. There is nothing Cantera
            # can do about this, the changes in CPython are handled by Cython. This
            # version bump is used to produce a more useful/actionable error message
            # for users than the compilation errors that result from using
            # Cython < 0.29.12.
            cython_min_version = LooseVersion("0.29.12")

        if numpy_version == LooseVersion('0.0.0'):
            print("NumPy not found.")
            warn_no_full_package = True
        elif numpy_version < numpy_min_test_version:
            print("WARNING: The installed version of Numpy is not tested and "
                  "support is not guaranteed. Found {0} but {1} or newer is preferred".format(
                  numpy_version, numpy_min_test_version))
        else:
            print('INFO: Using NumPy version {0}.'.format(numpy_version))

        if cython_version == LooseVersion('0.0.0'):
            print("Cython not found.")
            warn_no_full_package = True
        elif cython_version < cython_min_version:
            print("WARNING: Cython is an incompatible version: "
                  "Found {0} but {1} or newer is required.".format(
                    cython_version, cython_min_version))
            warn_no_full_package = True
        else:
            print('INFO: Using Cython version {0}.'.format(cython_version))

        if warn_no_full_package:
            if env['python_package'] == 'default':
                print('WARNING: Unable to build the full Python package because compatible '
                    'versions of Python, Numpy, and Cython could not be found.')
                print('INFO: Building the minimal Python package for Python {}'.format(python_version))
                env['python_package'] = 'minimal'
            else:
                print('ERROR: Unable to build the full Python package because compatible '
                    'versions of Python, Numpy, and Cython could not be found.')
                sys.exit(1)
        else:
            print('INFO: Building the full Python package for Python {0}'.format(python_version))
            env['python_package'] = 'full'

# Matlab Toolbox settings
if env['matlab_path'] != '' and env['matlab_toolbox'] == 'default':
    env['matlab_toolbox'] = 'y'

if env['matlab_toolbox'] == 'y':
    matPath = env['matlab_path']
    if matPath == '':
        print("ERROR: Unable to build the Matlab toolbox because 'matlab_path' has not been set.")
        sys.exit(1)

    if env['blas_lapack_libs']:
        print('ERROR: The Matlab toolbox is incompatible with external BLAS '
              'and LAPACK libraries. Unset blas_lapack_libs (e.g. "scons '
              'build blas_lapack_libs=") in order to build the Matlab '
              'toolbox, or set matlab_toolbox=n to use the specified BLAS/'
              'LAPACK libraries and skip building the Matlab toolbox.')
        sys.exit(1)

    if env['system_sundials'] == 'y':
        print('ERROR: The Matlab toolbox is incompatible with external '
              'SUNDIALS libraries. Set system_sundials to no (e.g., "scons build '
              'system_sundials=n") in order to build the Matlab '
              'toolbox, or set matlab_toolbox=n to use the specified '
              'SUNDIALS libraries and skip building the Matlab toolbox.')
        sys.exit(1)

    if not (os.path.isdir(matPath) and
            os.path.isdir(pjoin(matPath, 'extern'))):
        print("""ERROR: Path set for 'matlab_path' is not correct.""")
        print("""ERROR: Path was: '%s'""" % matPath)
        sys.exit(1)


# **********************************************
# *** Set additional configuration variables ***
# **********************************************

# On Debian-based systems, need to special-case installation to
# /usr/local because of dist-packages vs site-packages
env['debian'] = any(name.endswith('dist-packages') for name in sys.path)

# Directories where things will be after actually being installed. These
# variables are the ones that are used to populate header files, scripts, etc.
env['ct_installroot'] = env['prefix']
env['ct_libdir'] = pjoin(env['prefix'], env['libdirname'])
env['ct_bindir'] = pjoin(env['prefix'], 'bin')
env['ct_incdir'] = pjoin(env['prefix'], 'include', 'cantera')
env['ct_incroot'] = pjoin(env['prefix'], 'include')

if env['layout'] == 'compact':
    env['ct_datadir'] = pjoin(env['prefix'], 'data')
    env['ct_sampledir'] = pjoin(env['prefix'], 'samples')
    env['ct_mandir'] = pjoin(env['prefix'], 'man1')
    env['ct_matlab_dir'] = pjoin(env['prefix'], 'matlab', 'toolbox')
else:
    env['ct_datadir'] = pjoin(env['prefix'], 'share', 'cantera', 'data')
    env['ct_sampledir'] = pjoin(env['prefix'], 'share', 'cantera', 'samples')
    env['ct_mandir'] = pjoin(env['prefix'], 'share', 'man', 'man1')
    env['ct_matlab_dir'] = pjoin(env['prefix'], env['libdirname'],
                                 'cantera', 'matlab', 'toolbox')

# Always set the stage directory before building an MSI installer
if 'msi' in COMMAND_LINE_TARGETS:
    COMMAND_LINE_TARGETS.append('install')
    env['stage_dir'] = 'stage'
    env['prefix'] = '.'
    env['PYTHON_INSTALLER'] = 'binary'
elif env['layout'] == 'debian':
    COMMAND_LINE_TARGETS.append('install')
    env['stage_dir'] = 'stage/cantera'
    env['PYTHON_INSTALLER'] = 'debian'
    env['INSTALL_MANPAGES'] = False
else:
    env['PYTHON_INSTALLER'] = 'direct'


addInstallActions = ('install' in COMMAND_LINE_TARGETS or
                     'uninstall' in COMMAND_LINE_TARGETS)

# Directories where things will be staged for package creation. These
# variables should always be used by the Install(...) targets
if env['stage_dir']:
    instRoot = pjoin(os.getcwd(), env['stage_dir'],
                     stripDrive(env['prefix']).strip('/\\'))
    if env['python_prefix']:
        env['python_prefix'] = pjoin(os.getcwd(), env['stage_dir'],
                       stripDrive(env['python_prefix']).strip('/\\'))
    else:
        env['python_prefix'] = pjoin(os.getcwd(), env['stage_dir'])
else:
    instRoot = env['prefix']

# Prevent setting Cantera installation path to source directory
if os.path.abspath(instRoot) == Dir('.').abspath:
    print('ERROR: cannot install Cantera into source directory.')
    exit(1)

if env['layout'] == 'debian':
    base = pjoin(os.getcwd(), 'debian')

    env['inst_libdir'] = pjoin(base, 'cantera-dev', 'usr', env['libdirname'])
    env['inst_incdir'] = pjoin(base, 'cantera-dev', 'usr', 'include', 'cantera')
    env['inst_incroot'] = pjoin(base, 'cantera-dev', 'usr' 'include')

    env['inst_bindir'] = pjoin(base, 'cantera-common', 'usr', 'bin')
    env['inst_datadir'] = pjoin(base, 'cantera-common', 'usr', 'share', 'cantera', 'data')
    env['inst_docdir'] = pjoin(base, 'cantera-common', 'usr', 'share', 'cantera', 'doc')
    env['inst_sampledir'] = pjoin(base, 'cantera-common', 'usr', 'share', 'cantera', 'samples')
    env['inst_mandir'] = pjoin(base, 'cantera-common', 'usr', 'share', 'man', 'man1')

    env['inst_matlab_dir'] = pjoin(base, 'cantera-matlab', 'usr',
                                   env['libdirname'], 'cantera', 'matlab', 'toolbox')

    env['inst_python_bindir'] = pjoin(base, 'cantera-python', 'usr', 'bin')
    env['python_prefix'] = pjoin(base, 'cantera-python3')
else:
    env['inst_libdir'] = pjoin(instRoot, env['libdirname'])
    env['inst_bindir'] = pjoin(instRoot, 'bin')
    env['inst_python_bindir'] = pjoin(instRoot, 'bin')
    env['inst_incdir'] = pjoin(instRoot, 'include', 'cantera')
    env['inst_incroot'] = pjoin(instRoot, 'include')

    if env['layout'] == 'compact':
        env['inst_matlab_dir'] = pjoin(instRoot, 'matlab', 'toolbox')
        env['inst_datadir'] = pjoin(instRoot, 'data')
        env['inst_sampledir'] = pjoin(instRoot, 'samples')
        env['inst_docdir'] = pjoin(instRoot, 'doc')
        env['inst_mandir'] = pjoin(instRoot, 'man1')
    else: # env['layout'] == 'standard'
        env['inst_matlab_dir'] = pjoin(instRoot, env['libdirname'], 'cantera',
                                       'matlab', 'toolbox')
        env['inst_datadir'] = pjoin(instRoot, 'share', 'cantera', 'data')
        env['inst_sampledir'] = pjoin(instRoot, 'share', 'cantera', 'samples')
        env['inst_docdir'] = pjoin(instRoot, 'share', 'cantera', 'doc')
        env['inst_mandir'] = pjoin(instRoot, 'share', 'man', 'man1')

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

# Need to test all of these to see what platform.system() returns
configh['SOLARIS'] = 1 if env['OS'] == 'Solaris' else None
configh['DARWIN'] = 1 if env['OS'] == 'Darwin' else None

if env['OS'] == 'Solaris' or env['HAS_CLANG']:
    configh['NEEDS_GENERIC_TEMPL_STATIC_DECL'] = 1

configh['CT_SUNDIALS_VERSION'] = env['sundials_version'].replace('.','')

if env.get('has_sundials_lapack') and env['use_lapack']:
    configh['CT_SUNDIALS_USE_LAPACK'] = 1
else:
    configh['CT_SUNDIALS_USE_LAPACK'] = 0

cdefine('LAPACK_FTN_STRING_LEN_AT_END', 'lapack_ftn_string_len_at_end')
cdefine('LAPACK_FTN_TRAILING_UNDERSCORE', 'lapack_ftn_trailing_underscore')
cdefine('FTN_TRAILING_UNDERSCORE', 'lapack_ftn_trailing_underscore')
cdefine('LAPACK_NAMES_LOWERCASE', 'lapack_names', 'lower')
cdefine('CT_USE_LAPACK', 'use_lapack')
cdefine('CT_USE_SYSTEM_EIGEN', 'system_eigen')
cdefine('CT_USE_SYSTEM_FMT', 'system_fmt')
cdefine('CT_USE_SYSTEM_YAMLCPP', 'system_yamlcpp')
cdefine('CT_USE_DEMANGLE', 'has_demangle')

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

# preprocess input files (cti -> xml)
convertedInputFiles = set()
for cti in mglob(env, 'data/inputs', 'cti'):
    build(env.Command('build/data/%s' % cti.name, cti.path,
                      Copy('$TARGET', '$SOURCE')))
    outName = os.path.splitext(cti.name)[0] + '.xml'
    convertedInputFiles.add(outName)
    build(env.Command('build/data/%s' % outName, cti.path,
                      '$python_cmd_esc interfaces/cython/cantera/ctml_writer.py $SOURCE $TARGET'))

# Copy XML input files which are not present as cti:
for xml in mglob(env, 'data/inputs', 'xml'):
    dest = pjoin('build','data',xml.name)
    if xml.name not in convertedInputFiles:
        build(env.Command(dest, xml.path, Copy('$TARGET', '$SOURCE')))

for yaml in mglob(env, "data", "yaml"):
    dest = pjoin("build", "data", yaml.name)
    build(env.Command(dest, yaml.path, Copy("$TARGET", "$SOURCE")))

if addInstallActions:
    # Put headers in place
    headerBase = 'include/cantera'
    install(env.RecursiveInstall, '$inst_incdir', 'include/cantera')

    # Data files
    install('$inst_datadir', mglob(env, 'build/data', 'cti', 'xml', 'yaml'))


### List of libraries needed to link to Cantera ###
linkLibs = ['cantera']

### List of shared libraries needed to link applications to Cantera
linkSharedLibs = ['cantera_shared']

if env['system_sundials'] == 'y':
    env['sundials_libs'] = ['sundials_cvodes', 'sundials_ida', 'sundials_nvecserial']
    if env['use_lapack'] and sundials_ver >= LooseVersion('3.0'):
        if env.get('has_sundials_lapack'):
            env['sundials_libs'].extend(('sundials_sunlinsollapackdense',
                                         'sundials_sunlinsollapackband'))
        else:
            env['sundials_libs'].extend(('sundials_sunlinsoldense',
                                         'sundials_sunlinsolband'))
else:
    env['sundials_libs'] = []

linkLibs.extend(env['sundials_libs'])
linkSharedLibs.extend(env['sundials_libs'])

if env['system_yamlcpp']:
    linkLibs.append('yaml-cpp')
    linkSharedLibs.append('yaml-cpp')

#  Add LAPACK and BLAS to the link line
if env['blas_lapack_libs']:
    linkLibs.extend(env['blas_lapack_libs'])
    linkSharedLibs.extend(env['blas_lapack_libs'])

# Store the list of needed static link libraries in the environment
env['cantera_libs'] = linkLibs
env['cantera_shared_libs'] = linkSharedLibs
if not env['renamed_shared_libraries']:
    env['cantera_shared_libs'] = linkLibs

# Add targets from the SConscript files in the various subdirectories
Export('env', 'build', 'libraryTargets', 'install', 'buildSample')

# ext needs to come before src so that libraryTargets is fully populated
VariantDir('build/ext', 'ext', duplicate=0)
SConscript('build/ext/SConscript')

# Fortran needs to come before src so that libraryTargets is fully populated
if env['f90_interface'] == 'y':
    VariantDir('build/src/fortran/', 'src/fortran', duplicate=1)
    SConscript('build/src/fortran/SConscript')

VariantDir('build/src', 'src', duplicate=0)
SConscript('build/src/SConscript')

if env['python_package'] == 'full':
    SConscript('interfaces/cython/SConscript')
elif env['python_package'] == 'minimal':
    SConscript('interfaces/python_minimal/SConscript')

if env['CC'] != 'cl':
    VariantDir('build/platform', 'platform/posix', duplicate=0)
    SConscript('build/platform/SConscript')

if env['matlab_toolbox'] == 'y':
    SConscript('build/src/matlab/SConscript')

if env['doxygen_docs'] or env['sphinx_docs']:
    SConscript('doc/SConscript')

# Sample programs (also used from test_problems/SConscript)
VariantDir('build/samples', 'samples', duplicate=0)
sampledir_excludes = ['\\.o$', '^~$', '\\.in', 'SConscript']
SConscript('build/samples/cxx/SConscript')

# Install C++ samples
install(env.RecursiveInstall, '$inst_sampledir/cxx',
        'samples/cxx', exclude=sampledir_excludes)

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
    print("*******************************************************")
    print("Compilation completed successfully.\n")
    print("- To run the test suite, type 'scons test'.")
    if env['googletest'] == 'none':
        print("  WARNING: You set the 'googletest' to 'none' and all it's tests will be skipped.")
    if os.name == 'nt':
        print("- To install, type 'scons install'.")
        print("- To create a Windows MSI installer, type 'scons msi'.")
    else:
        print("- To install, type '[sudo] scons install'.")
    print("*******************************************************")

finish_build = env.Command('finish_build', [], postBuildMessage)
env.Depends(finish_build, buildTargets)
build_cantera = Alias('build', finish_build)

Default('build')

def postInstallMessage(target, source, env):
    # Only needed because Python 2 doesn't support textwrap.indent
    def indent(inp_str, indent):
        return '\n'.join([indent + spl for spl in inp_str.splitlines()])

    env_dict = env.Dictionary()
    install_message = textwrap.dedent("""
        Cantera has been successfully installed.

        File locations:

          applications                {ct_bindir!s}
          library files               {ct_libdir!s}
          C++ headers                 {ct_incroot!s}
          samples                     {ct_sampledir!s}
          data files                  {ct_datadir!s}""".format(**env_dict))

    if env['sphinx_docs'] or env['doxygen_docs']:
        install_message += indent(textwrap.dedent("""
            HTML documentation          {inst_docdir!s}""".format(**env_dict)), '  ')

    if env['python_package'] == 'full':
        env['python_example_loc'] = pjoin(env['python_module_loc'], 'cantera', 'examples')
        install_message += indent(textwrap.dedent("""
            Python package (cantera)    {python_module_loc!s}
            Python samples              {python_example_loc!s}""".format(**env_dict)), '  ')

    if env['matlab_toolbox'] == 'y':
        env['matlab_sample_loc'] = pjoin(env['ct_sampledir'], 'matlab')
        env['matlab_ctpath_loc'] = pjoin(env['ct_matlab_dir'], 'ctpath.m')
        install_message += textwrap.dedent("""
              Matlab toolbox              {ct_matlab_dir!s}
              Matlab samples              {matlab_sample_loc!s}

            An m-file to set the correct matlab path for Cantera is at:

              {matlab_ctpath_loc!s}
        """.format(**env_dict))

    if os.name != 'nt':
        install_message += textwrap.dedent("""
            Setup scripts to configure the environment for Cantera are at:

              setup script (bash)         {ct_bindir!s}/setup_cantera
              setup script (csh/tcsh)     {ct_bindir!s}/setup_cantera.csh

            It is recommended that you run the script for your shell by typing:

              source {ct_bindir!s}/setup_cantera

            before using Cantera, or else include its contents in your shell login script.
        """.format(**env_dict))

    print(install_message)

finish_install = env.Command('finish_install', [], postInstallMessage)
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

# Files installed by the Python installer
if os.path.exists('build/python-installed-files.txt'):
    with open('build/python-installed-files.txt', 'r') as f:
        file_list = f.readlines()

    install_base = os.path.dirname(file_list[0].strip())
    if os.path.exists(install_base):
        not_listed_files = [s for s in os.listdir(install_base) if not any(s in j for j in file_list)]
        for f in not_listed_files:
            f = pjoin(install_base, f)
            if not os.path.isdir(f) and os.path.exists(f):
                allfiles.append(File(f))
    for f in file_list:
        f = f.strip()
        if not os.path.isdir(f) and os.path.exists(f):
            allfiles.append(File(f))

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
    env['test_results'] = env.Command('test_results', [], testResults.printReport)

    if env['python_package'] == 'none':
        # copy scripts from the full Cython module
        test_py_int = env.Command('#build/python_local/cantera/__init__.py',
                                  '#interfaces/python_minimal/cantera/__init__.py',
                                  Copy('$TARGET', '$SOURCE'))
        for script in ['ctml_writer', 'ck2cti', 'ck2yaml', 'ctml2yaml']:
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
