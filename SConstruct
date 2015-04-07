"""
SCons build script for Cantera

Basic usage:
    'scons help' - print a description of user-specifiable options.

    'scons build' - Compile Cantera and the language interfaces using
                    default options.

    'scons clean' - Delete files created while building Cantera.

    '[sudo] scons install' - Install Cantera.

    '[sudo] scons uninstall' - Uninstall Cantera.

    'scons test' - Run full test suite.

    'scons test-clean' - Delete files created while running the tests.

    'scons test-help' - List available tests.

    'scons test-NAME' - Run the test named "NAME".

    'scons samples' - Compile the C++ and Fortran samples

    'scons msi' - Build a Windows installer (.msi) for Cantera.
"""

from buildutils import *

if not COMMAND_LINE_TARGETS:
    # Print usage help
    print __doc__
    sys.exit(0)

valid_commands = ('build','clean','install','uninstall',
                  'help','msi','samples')

for command in COMMAND_LINE_TARGETS:
    if command not in valid_commands and not command.startswith('test'):
        print 'ERROR: unrecognized command line target: %r' % command
        sys.exit(0)

extraEnvArgs = {}

if 'clean' in COMMAND_LINE_TARGETS:
    removeDirectory('build')
    removeDirectory('stage')
    removeDirectory('interfaces/python/build')
    removeDirectory('.sconf_temp')
    removeFile('.sconsign.dblite')
    removeFile('include/cantera/base/config.h')
    removeFile('interfaces/python/setup.py')
    removeFile('ext/f2c_libs/arith.h')
    for name in os.listdir('.'):
        if name.endswith('.msi'):
            removeFile(name)
    for name in os.listdir('interfaces/python/Cantera'):
        if name.startswith('_cantera') or name.startswith('cantera_shared'):
            removeFile('interfaces/python/Cantera/' + name)
    removeFile('interfaces/matlab/toolbox/cantera_shared.dll')
    for name in os.listdir('interfaces/matlab/toolbox'):
        if name.startswith('ctmethods.'):
            removeFile('interfaces/matlab/toolbox/' + name)
    print 'Done removing output files.'
    sys.exit(0)

# ******************************************************
# *** Set system-dependent defaults for some options ***
# ******************************************************

opts = Variables('cantera.conf')

windows_compiler_options = []
if os.name == 'nt':
    # On Windows, use the same version of Visual Studio that was used
    # to compile Python, and target the same architecture, unless
    # the user specified another option

    if 'MSC v.1400' in sys.version:
        msvc_version = '8.0' # Visual Studio 2005
    elif 'MSC v.1500' in sys.version:
        msvc_version = '9.0' # Visual Studio 2008
    elif 'MSC v.1600' in sys.version:
        msvc_version = '10.0' # Visual Studio 2010
    elif 'MSC v.1700' in sys.version:
        msvc_version = '11.0' # Visual Studio 2012
    else:
        msvc_version = None

    if '64 bit' in sys.version:
        target_arch = 'amd64'
    else:
        target_arch = 'x86'

    windows_compiler_options.extend([
        ('msvc_version',
         """Version of Visual Studio to use. The default is the same version
            that was used to compile the installed version of Python.""",
         msvc_version),
        ('target_arch',
         """Target architecture. The default is the same
            architecture as the installed version of Python""",
         target_arch),
        EnumVariable(
            'toolchain',
            """The preferred compiler toolchain.""",
            'msvc', ('msvc', 'mingw', 'intel'))])
    opts.AddVariables(*windows_compiler_options)

    pickCompilerEnv = Environment()
    opts.Update(pickCompilerEnv)

    if pickCompilerEnv['toolchain'] == 'msvc':
        toolchain = ['default']
        if msvc_version:
            extraEnvArgs['MSVC_VERSION'] = pickCompilerEnv['msvc_version']

    elif pickCompilerEnv['toolchain'] == 'mingw':
        toolchain = ['mingw', 'f90']
        extraEnvArgs['F77'] = None
        # Next line fixes http://scons.tigris.org/issues/show_bug.cgi?id=2683
        extraEnvArgs['WINDOWS_INSERT_DEF'] = 1

    elif pickCompilerEnv['toolchain'] == 'intel':
        toolchain = ['intelc'] # note: untested

    extraEnvArgs['TARGET_ARCH'] = pickCompilerEnv['target_arch']

else:
    toolchain = ['default']

env = Environment(tools=toolchain+['textfile', 'subst', 'recursiveInstall', 'wix'],
                  ENV={'PATH': os.environ['PATH']},
                  toolchain=toolchain,
                  **extraEnvArgs)

#
# To print the current environment
#
# print env.Dump()


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
    for key,val in env['ENV'].iteritems():
        env['ENV'][key] = str(val)

if 'FRAMEWORKS' not in env:
    env['FRAMEWORKS'] = []

add_RegressionTest(env)

class defaults: pass

if os.name == 'posix':
    defaults.prefix = '/usr/local'
    defaults.boostIncDir = ''
    defaults.boostLibDir = ''
    env['INSTALL_MANPAGES'] = True
elif os.name == 'nt':
    defaults.prefix = pjoin(os.environ['ProgramFiles'], 'Cantera')
    defaults.boostIncDir = ''
    defaults.boostLibDir = ''
    env['INSTALL_MANPAGES'] = False
else:
    print "Error: Unrecognized operating system '%s'" % os.name
    sys.exit(1)

compiler_options = [
    ('CXX',
     'The C++ compiler to use.',
     env['CXX']),
    ('CC',
     """The C compiler to use. This is only used to compile CVODE and
        the Python extension module.""",
     env['CC'])]
opts.AddVariables(*compiler_options)
opts.Update(env)

defaults.cxxFlags = ''
defaults.ccFlags = ''
defaults.noOptimizeCcFlags = ''
defaults.optimizeCcFlags = ''
defaults.debugCcFlags = ''
defaults.noDebugCcFlags = ''
defaults.debugLinkFlags = ''
defaults.noDebugLinkFlags = ''
defaults.warningFlags = ''

if env['CC'] == 'gcc' or env['CC'] == 'llvm-gcc':
    defaults.cxxFlags = '-ftemplate-depth-128'
    defaults.debugCcFlags = '-g'
    defaults.noOptimizeCcFlags = '-O0 -fno-inline'
    defaults.optimizeCcFlags = '-O3 -DNDEBUG -finline-functions -Wno-inline'
    defaults.warningFlags = '-Wall'

elif env['CC'] == 'cl': # Visual Studio
    defaults.cxxFlags = ['/EHsc']
    defaults.ccFlags = ['/MD', '/nologo', '/Zc:wchar_t', '/Zc:forScope',
                        '/D_SCL_SECURE_NO_WARNINGS', '/D_CRT_SECURE_NO_WARNINGS']
    if env['MSVC_VERSION'] == '11.0':
        # Fix compatibility issue between VS2012 and Google Test
        defaults.cxxFlags.append('/D_VARIADIC_MAX=10')
    defaults.debugCcFlags = '/Zi /Fd${TARGET}.pdb'
    defaults.noOptimizeCcFlags = '/Od /Ob0'
    defaults.optimizeCcFlags = '/O2 /DNDEBUG'
    defaults.debugLinkFlags = '/DEBUG'
    defaults.warningFlags = '/W3'

elif env['CC'] == 'icc':
    defaults.cxxFlags = '-ftemplate-depth-128'
    defaults.ccFlags = '-vec-report0 -diag-disable 1478'
    defaults.debugCcFlags = '-g'
    defaults.noOptimizeCcFlags = '-O0 -fno-inline'
    defaults.optimizeCcFlags = '-O3 -finline-functions -DNDEBUG'
    defaults.warningFlags = '-Wcheck'

elif env['CC'] == 'clang':
    defaults.cxxFlags = ''
    defaults.ccFlags = '-fcolor-diagnostics'
    defaults.debugCcFlags = '-g'
    defaults.noOptimizeCcFlags = '-O0'
    defaults.optimizeCcFlags = '-O3 -DNDEBUG'
    defaults.warningFlags = '-Wall'

else:
    print "WARNING: Unrecognized C compiler '%s'" % env['CC']

if env['OS'] in ('Windows', 'Darwin'):
    defaults.threadFlags = ''
else:
    defaults.threadFlags = '-pthread'

defaults.fsLayout = 'compact' if env['OS'] == 'Windows' else 'standard'
defaults.env_vars = 'LD_LIBRARY_PATH,PYTHONPATH'
defaults.python_prefix = '$prefix' if env['OS'] != 'Windows' else ''

# Transform lists into strings to keep cantera.conf clean
for key,value in defaults.__dict__.items():
    if isinstance(value, (list, tuple)):
        defaults.__dict__[key] = ' '.join(value)

# **************************************
# *** Read user-configurable options ***
# **************************************

# This check prevents requiring root permissions to build when the
# installation directory is not writable by the current user.
if 'install' in COMMAND_LINE_TARGETS:
    installPathTest = PathVariable.PathIsDirCreate
else:
    installPathTest = PathVariable.PathAccept

config_options = [
    PathVariable(
        'prefix',
        'Set this to the directory where Cantera should be installed.',
        defaults.prefix, installPathTest),
    EnumVariable(
        'python_package',
        """If you plan to work in Python, or you want to use the graphical
           MixMaster application, then you need either the 'new' or 'full'
           Cantera Python Package. If, on the other hand, you will only use
           Cantera from some other language (e.g. MATLAB or Fortran 90/95) and
           only need Python to process .cti files, then you only need a
           'minimal' subset of the package (actually, only one file). The
           default behavior is to build the Python package if the required
           prerequisites (numpy) are installed. NOTE: The legacy 'full' option
           is deprecated in favor of the 'new' Python package. The legacy
           Python package will be removed in Cantera 2.2 """,
        'default', ('new', 'full', 'minimal', 'none', 'default')),
    PathVariable(
        'python_cmd',
        """Cantera needs to know where to find the Python interpreter. If
           PYTHON_CMD is not set, then the configuration process will use the
           same Python interpreter being used by SCons.""",
        sys.executable),
    EnumVariable(
        'python_array',
        """The Cantera Python interface requires one of the Python array
           packages listed. Support for the legacy 'numeric' and 'numarray'
           packages is deprecated, and will be removed in a future version
           of Cantera.""",
        'numpy', ('numpy', 'numarray', 'numeric')),
    PathVariable(
        'python_array_home',
        """If numpy was installed using the --home option, set this to
           the home directory for numpy.""",
        '', PathVariable.PathAccept),
    PathVariable(
        'python_prefix',
        """Use this option if you want to install the Cantera Python package to
           an alternate location. On Unix-like systems, the default is the same
           as the $prefix option. If this option is set to the empty string (the
           default on Windows), then the Package will be installed to the system
           default 'site-packages' directory. To install to the current user's
           site-packages directory, use 'python_prefix=USER'.""",
        defaults.python_prefix, PathVariable.PathAccept),
    EnumVariable(
        'python3_package',
        """Controls whether or not the Python 3 module will be built. By
            default, the module will be built if the Python 3 interpreter can
            be found""",
        'default', ('y','n','default')),
    PathVariable(
        'python3_cmd',
        """ The name (full path if necessary) of the Python 3 interpreter.
        Required to build the Python 3 module.""",
        'python3', PathVariable.PathAccept),
    PathVariable(
        'python3_array_home',
        """"If numpy was installed to a custom location (e.g. using the --home
            option, set this to the directory for numpy""",
        '', PathVariable.PathAccept),
    PathVariable(
        'python3_prefix',
        """Use this option if you want to install the Cantera Python 3 package to
           an alternate location. On Unix-like systems, the default is the same
           as the $prefix option. If this option is set to the empty string (the
           default on Windows), then the Package will be installed to the system
           default 'site-packages' directory. To install to the current user's
           site-packages directory, use 'python3_prefix=USER'.""",
        defaults.python_prefix, PathVariable.PathAccept),
    PathVariable(
        'python_compiler',
        """ Compiler to use while building the Python extension module. By default,
            the compiler will be selected by distutils. Applies to both Python 2
            and Python 3.""",
        '', PathVariable.PathAccept),
    EnumVariable(
        'matlab_toolbox',
        """This variable controls whether the Matlab toolbox will be built. If
           set to 'y', you will also need to set the value of the 'matlab_path'
           variable. If set to 'default', the Matlab toolbox will be built if
           'matlab_path' is set.""",
        'default', ('y', 'n', 'default')),
    PathVariable(
        'matlab_path',
        """Path to the Matlab install directory. This should be the directory
           containing the 'extern', 'bin', etc. subdirectories. Typical values
           are: "C:/Program Files/MATLAB/R2011a" on Windows,
           "/Applications/MATLAB_R2011a.app" on OS X, or
           "/opt/MATLAB/R2011a" on Linux.""",
        '', PathVariable.PathAccept),
    EnumVariable(
        'f90_interface',
        """This variable controls whether the Fortran 90/95 interface will be
           built. If set to 'default', the builder will look for a compatible
           Fortran compiler in the $PATH, and compile the Fortran 90 interface
           if one is found.""",
        'default', ('y', 'n', 'default')),
    PathVariable(
        'F90',
        """The Fortran 90 compiler. If unspecified, the builder will look for a
           compatible compiler (gfortran, ifort, g95) in the $PATH.""",
        '', PathVariable.PathAccept),
    ('F90FLAGS',
     'Compilation options for the Fortran 90 compiler.',
     '-O3'),
    BoolVariable(
        'debug_verbose',
        """Enable extra printing to aid in debugging. This code is marked
            by the preprocessor macros DEBUG_MODE and DEBUG_MODE_ENABLED.""",
        False),
    BoolVariable(
        'coverage',
        """Enable collection of code coverage information with gcov.
           Available only when compiling with gcc.""",
        False),
    BoolVariable(
        'doxygen_docs',
        """Build HTML documentation for the C++ interface using Doxygen""",
        False),
    BoolVariable(
        'sphinx_docs',
        """Build HTML documentation for the Python module using Sphinx""",
        False),
    PathVariable(
        'sphinx_cmd',
        """Command to use for building the Sphinx documentation""",
        'sphinx-build', PathVariable.PathAccept),
    BoolVariable(
        'with_h298modify_capability',
        """Enable changing the 298K heats of formation directly via
           the C++ layer.""",
        False),
    BoolVariable(
        'with_html_log_files',
        """write HTML log files. Some multiphase equilibrium procedures can
           write copious diagnostic log messages. Deprecated -- to be removed
           in Cantera 2.2.""",
        False),
    EnumVariable(
        'use_sundials',
        """Cantera uses the CVODE or CVODES ODE integrator to time-integrate
           reactor network ODE's and for various other purposes. An older
           version of CVODE comes with Cantera, but it is possible to use the
           latest version as well, which now supports sensitivity analysis
           (CVODES). CVODES is a part of the 'sundials' package from Lawrence
           Livermore National Laboratory. Sundials is not distributed with
           Cantera, but it is free software that may be downloaded and
           installed separately. If you leave USE_SUNDIALS = 'default', then it
           will be used if you have it, and if not the older CVODE will be
           used. Or set USE_SUNDIALS to 'y' or 'n' to force using it or not.
           Note that sensitivity analysis with Cantera requires use of
           sundials. See: http://www.llnl.gov/CASC/sundials""",
        'default', ('default', 'y', 'n')),
    PathVariable(
        'sundials_include',
        """The directory where the Sundials header files are installed. This
           should be the directory that contains the "cvodes", "nvector", etc.
           subdirectories. Not needed if the headers are installed in a
           standard location, e.g. /usr/include.""",
        '', PathVariable.PathAccept),
    PathVariable(
        'sundials_libdir',
        """The directory where the sundials static libraries are installed.
           Not needed if the libraries are installed in a standard location,
           e.g. /usr/lib.""",
        '', PathVariable.PathAccept),
    PathVariable(
        'sundials_license',
        """Path to the sundials LICENSE file. Needed so that it can be included
           when bundling Sundials""",
        '', PathVariable.PathAccept),
    BoolVariable(
        'install_sundials',
        """Determines whether Sundials library and header files are installed
           alongside Cantera. Intended for use when installing on Windows.""",
        os.name == 'nt'),
    ('blas_lapack_libs',
     """Cantera comes with Fortran (or C) versions of those parts of BLAS and
        LAPACK it requires. But performance may be better if you use a version
        of these libraries optimized for your machine hardware. If you want to
        use your own libraries, set blas_lapack_libs to the the list of
        libraries that should be passed to the linker, separated by commas,
        e.g. "lapack,blas" or "lapack,f77blas,cblas,atlas". """,
     ''),
    PathVariable('blas_lapack_dir',
        """Directory containing the libraries specified by 'blas_lapack_libs'.""",
        '', PathVariable.PathAccept),
    EnumVariable(
        'lapack_names',
        """Set depending on whether the procedure names in the specified
           libraries are lowercase or uppercase. If you don't know, run 'nm' on
           the library file (e.g. 'nm libblas.a').""",
        'lower', ('lower','upper')),
    BoolVariable(
        'lapack_ftn_trailing_underscore', '', True),
    BoolVariable(
        'lapack_ftn_string_len_at_end', '', True),
    ('env_vars',
     """Environment variables to propagate through to SCons. Either the
        string "all" or a comma separated list of variable names, e.g.
        'LD_LIBRARY_PATH,HOME'""",
     defaults.env_vars),
    ('cxx_flags',
     'Compiler flags passed to the C++ compiler only.',
     defaults.cxxFlags),
    ('cc_flags',
     'Compiler flags passed to both the C and C++ compilers, regardless of optimization level',
     defaults.ccFlags),
    ('thread_flags',
     'Compiler and linker flags for POSIX multithreading support',
     defaults.threadFlags),
    BoolVariable(
        'optimize',
        """Enable extra compiler optimizations specified by the
           "optimize_flags" variable, instead of the flags specified by the
           "debug_flags" variable""",
        True),
    ('optimize_flags',
     'Additional compiler flags passed to the C/C++ compiler when optimize=yes.',
     defaults.optimizeCcFlags),
    ('no_optimize_flags',
     'Additional compiler flags passed to the C/C++ compiler when optimize=no.',
     defaults.noOptimizeCcFlags),
    BoolVariable(
        'debug',
        """Enable compiler debugging symbols.""",
        True),
    ('debug_flags',
     'Additional compiler flags passed to the C/C++ compiler when debug=yes.',
     defaults.debugCcFlags),
    ('no_debug_flags',
     'Additional compiler flags passed to the C/C++ compiler when debug=no.',
     defaults.noDebugCcFlags),
    ('debug_linker_flags',
     'Additional options passed to the linker when debug=yes',
     defaults.debugLinkFlags),
    ('no_debug_linker_flags',
     'Additional options passed to the linker when debug=yes',
     defaults.noDebugLinkFlags),
    ('warning_flags',
     """Additional compiler flags passed to the C/C++ compiler to enable
        extra warnings. Used only when compiling source code that part of
        Cantera (e.g. excluding code in the 'ext' directory).""",
     defaults.warningFlags),
    ('extra_inc_dirs',
     'Additional directories to search for header files (colon-separated list)',
     ''),
    ('extra_lib_dirs',
     'Additional directories to search for libraries (colon-separated list)',
     ''),
    BoolVariable(
        'build_thread_safe',
        """Cantera can be built so that it is thread safe. Doing so requires
           using procedures from the Boost library, so if you want thread
           safety then you need to get and install Boost (http://www.boost.org)
           if you don't have it.  This is turned off by default, in which case
           Boost is not required to build Cantera.""",
        False),
    PathVariable(
        'boost_inc_dir',
        'Location of the Boost header files',
        defaults.boostIncDir, PathVariable.PathAccept),
    PathVariable(
        'boost_lib_dir',
        'Directory containing the Boost.Thread library',
        defaults.boostLibDir, PathVariable.PathAccept),
    ('boost_thread_lib',
     'The name of the Boost.Thread library.',
     'boost_thread'),
    ('boost_windows_libs',
     """Comma-separated list containing the names of the Boost libraries
        required to link Cantera programs on Windows. These libraries will be
        copied to the Cantera installation directory.""",
     'thread,system,date_time,chrono'), # default is correct for Boost 1.54
    BoolVariable(
        'build_with_f2c',
        """For external procedures written in Fortran 77, both the
           original F77 source code and C source code generated by the
           'f2c' program are included.  Set this to "n" if you want to
           build Cantera using the F77 sources in the ext directory.""",
        True),
    ('F77',
     """Compiler used to build the external Fortran 77 procedures from
        the Fortran source code""",
     env.get('F77')),
    ('F77FLAGS',
     """Fortran 77 Compiler flags. Note that the Fortran compiler
      flags must be set to produce object code compatible with the
      C/C++ compiler you are using.""",
     '-O3'),
    PathVariable(
        'stage_dir',
        """ Directory relative to the Cantera source directory to be
            used as a staging area for building e.g. a Debian
            package. If specified, 'scons install' will install files
            to 'stage_dir/prefix/...' instead of installing into the
            local filesystem.""",
        '',
        PathVariable.PathAccept),
    BoolVariable(
        'legacy_headers',
        """Create symbolic links for headers that were installed to the
        'kernel' subdirectory in previous versions of Cantera.""",
        False),
    BoolVariable(
        'VERBOSE',
        """Create verbose output about what scons is doing.""",
        False),
    BoolVariable(
        'renamed_shared_libraries',
        """If this option is turned on, the shared libraries that are created
        will be renamed to have a "_shared" extension added to their base name.
        If not, the base names will be the same as the static libraries.
        In some cases this simplifies subsequent linking environments with
        static libaries and avoids a bug with using valgrind with
        the -static linking flag.""",
        True),
    BoolVariable(
        'install_gtest',
        """Determines whether to install the Google Test library and headers
           alongside Cantera's libraries. Not recommended.""",
        False),
    BoolVariable(
        'single_library',
        """If set to 'n', code from the 'ext' folder in the cantera library
        will be placed in separate libraries for ctlapack, ctf2c, etc., rather
        than being included in the cantera library. NOTE: This option is
        deprecated, and will be removed in Cantera 2.2.""",
        True),
    EnumVariable(
        'layout',
        """The layout of the directory structure. 'standard' installs files to
           several subdirectories under 'prefix', e.g. $prefix/bin,
           $prefix/include/cantera, $prefix/lib. This layout is best used in
           conjunction with 'prefix'='/usr/local'. 'compact' puts all installed
           files in the subdirectory define by 'prefix'. This layout is best for
           with a prefix like '/opt/cantera'. 'debian' installs to the stage
           directory in a layout used for generating Debian packages.""",
     defaults.fsLayout, ('standard','compact','debian')),
    PathVariable(
        'graphvizdir',
        """The directory location of the graphviz program, "dot". dot is used
           for creating the documentation, and for making reaction path
           diagrams. If "dot" is in your path, you can leave this unspecified.
           NOTE: Matlab comes with a stripped-down version of 'dot'. If 'dot'
           is on your path, make sure it is not the Matlab version!""",
        '', PathVariable.PathAccept),
    ('rpfont',
     """The font to use in reaction path diagrams. This must be a font
        name recognized by the 'dot' program. On linux systems, this
        should be lowercase 'helvetica'.""",
     'Helvetica'),
    ('cantera_version', '', '2.1.1')
]

opts.AddVariables(*config_options)
opts.Update(env)
opts.Save('cantera.conf', env)


if 'help' in COMMAND_LINE_TARGETS:
    ### Print help about configuration options and exit.
    print """
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
"""

    for opt in opts.options:
        print '\n'.join(formatOption(env, opt))
    sys.exit(0)

valid_arguments = (set(opt[0] for opt in windows_compiler_options) |
                   set(opt[0] for opt in compiler_options) |
                   set(opt[0] for opt in config_options))
for arg in ARGUMENTS:
    if arg not in valid_arguments:
        print 'Encountered unexpected command line argument: %r' % arg
        sys.exit(0)

# Print values of all build options:
print "Configuration variables read from 'cantera.conf' and command line:"
for line in open('cantera.conf'):
    print '   ', line.strip()
print

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
            env['ENV'][name] = os.environ[name]
        elif name not in defaults.env_vars:
            print 'WARNING: failed to propagate environment variable', name

env['extra_inc_dirs'] = [d for d in env['extra_inc_dirs'].split(':') if d]
env['extra_lib_dirs'] = [d for d in env['extra_lib_dirs'].split(':') if d]

env.Append(CPPPATH=env['extra_inc_dirs'],
           LIBPATH=env['extra_lib_dirs'])

# Try to find a Fortran compiler:
if env['f90_interface'] in ('y','default'):
    foundF90 = False
    if env['F90']:
        env['f90_interface'] = 'y'
        if which(env['F90']) is not None:
            foundF90 = True
        else:
            print "WARNING: Couldn't find specified Fortran compiler: '%s'" % env['F90']

    for compiler in ['gfortran', 'ifort', 'g95']:
        if foundF90:
            break
        if which(compiler) is not None:
            print "INFO: Using '%s' to build the Fortran 90 interface" % which(compiler)
            env['F90'] = compiler
            foundF90 = True

    if foundF90:
        env['f90_interface'] = 'y'
        env['FORTRAN'] = env['F90']
    elif env['f90_interface'] == 'y':
        print "ERROR: Couldn't find a suitable Fortran compiler to build the Fortran 90 interface"
        sys.exit(1)
    else:
        print "INFO: Skipping compilation of the Fortran 90 interface."

if 'gfortran' in env['F90']:
    env['FORTRANMODDIRPREFIX'] = '-J'
    env['FORTRANSYSLIBS'] = ['gfortran']
elif 'g95' in env['F90']:
    env['FORTRANMODDIRPREFIX'] = '-fmod='
    env['FORTRANSYSLIBS'] = ['f95']
elif 'ifort' in env['F90']:
    env['FORTRANMODDIRPREFIX'] = '-module '
    env['FORTRANSYSLIBS'] = []
else:
    env['FORTRANSYSLIBS'] = []

env['FORTRANMODDIR'] = '${TARGET.dir}'

if env['CC'] == 'cl':
    # embed manifest file
    env['LINKCOM'] = [env['LINKCOM'],
                      'if exist ${TARGET}.manifest mt.exe -nologo -manifest ${TARGET}.manifest -outputresource:$TARGET;1']
    env['SHLINKCOM'] = [env['SHLINKCOM'],
                        'if exist ${TARGET}.manifest mt.exe -nologo -manifest ${TARGET}.manifest -outputresource:$TARGET;2']

if env['boost_inc_dir']:
    env.Append(CPPPATH=env['boost_inc_dir'])

if env['boost_lib_dir']:
    env.Append(LIBPATH=[env['boost_lib_dir']])

if env['blas_lapack_dir']:
    env.Append(LIBPATH=[env['blas_lapack_dir']])

if (env['use_sundials'] == 'default' and
    (env['sundials_include'] or env['sundials_libdir'])):
    env['use_sundials'] = 'y'

if env['use_sundials'] in ('y','default'):
    if env['sundials_include']:
        env.Append(CPPPATH=[env['sundials_include']])
    if env['sundials_libdir']:
        env.Append(LIBPATH=[env['sundials_libdir']])


# BLAS / LAPACK configuration
if env['blas_lapack_libs'] != '':
    env['blas_lapack_libs'] = env['blas_lapack_libs'].split(',')
    env['BUILD_BLAS_LAPACK'] = False
elif env['OS'] == 'Darwin':
    env['blas_lapack_libs'] = []
    env['BUILD_BLAS_LAPACK'] = False
    env.Append(FRAMEWORKS=['Accelerate'])
else:
    # External BLAS/LAPACK were not given, so we need to compile them
    env['BUILD_BLAS_LAPACK'] = True
    env['blas_lapack_libs'] = [] # built into libcantera

# ************************************
# *** Compiler Configuration Tests ***
# ************************************

conf = Configure(env)

# Set up compiler options before running configuration tests
env['CXXFLAGS'] = listify(env['cxx_flags'])
env['CCFLAGS'] = listify(env['cc_flags']) + listify(env['thread_flags'])
env['LINKFLAGS'] += listify(env['thread_flags'])

env['warning_flags'] = listify(env['warning_flags'])

if env['optimize']:
    env['CCFLAGS'] += listify(env['optimize_flags'])
else:
    env['CCFLAGS'] += listify(env['no_optimize_flags'])

if env['debug']:
    env['CCFLAGS'] += listify(env['debug_flags'])
    env['LINKFLAGS'] += listify(env['debug_linker_flags'])
else:
    env['CCFLAGS'] += listify(env['no_debug_flags'])
    env['LINKFLAGS'] += listify(env['no_debug_linker_flags'])

if env['coverage']:
    if  env['CC'] == 'gcc':
        env.Append(CCFLAGS=['-fprofile-arcs', '-ftest-coverage'])
        env.Append(LINKFLAGS=['-fprofile-arcs', '-ftest-coverage'])

    else:
        print 'Error: coverage testing is only available with GCC'
        exit(0)

if env['toolchain'] == 'mingw':
    env.Append(LINKFLAGS=['-static-libgcc', '-static-libstdc++'])

def config_error(message):
    print 'ERROR:', message
    print "See 'config.log' for details."
    sys.exit(1)

# First, a sanity check:
if not conf.CheckCXXHeader('cmath', '<>'):
    config_error('The C++ compiler is not correctly configured.')

def get_expression_value(includes, expression):
    s = ['#include ' + i for i in includes]
    s.extend(('#include <iostream>',
              'int main(int argc, char** argv) {',
              '    std::cout << %s << std::endl;' % expression,
              '    return 0;',
              '}\n'))
    return '\n'.join(s)

env['HAS_TIMES_H'] = conf.CheckCHeader('sys/times.h', '""')
env['HAS_UNISTD_H'] = conf.CheckCHeader('unistd.h', '""')
env['HAS_MATH_H_ERF'] = conf.CheckDeclaration('erf', '#include <math.h>', 'C++')
env['HAS_BOOST_MATH'] = conf.CheckCXXHeader('boost/math/special_functions/erf.hpp', '<>')
boost_version_source = get_expression_value(['<boost/version.hpp>'], 'BOOST_LIB_VERSION')
retcode, boost_lib_version = conf.TryRun(boost_version_source, '.cpp')
env['BOOST_LIB_VERSION'] = boost_lib_version.strip()

import SCons.Conftest, SCons.SConf
ret = SCons.Conftest.CheckLib(SCons.SConf.CheckContext(conf),
                              ['sundials_cvodes'],
                              header='#include "cvodes/cvodes.h"',
                              language='C++',
                              call='CVodeCreate(CV_BDF, CV_NEWTON);',
                              autoadd=False,
                              extra_libs=env['blas_lapack_libs'])
env['HAS_SUNDIALS'] = not ret # CheckLib returns False to indicate success
env['NEED_LIBM'] = not conf.CheckLibWithHeader(None, 'math.h', 'C',
                                               'double x; log(x);', False)
env['LIBM'] = ['m'] if env['NEED_LIBM'] else []

if env['HAS_SUNDIALS'] and env['use_sundials'] != 'n':
    # Determine Sundials version
    sundials_version_source = get_expression_value(['"sundials/sundials_config.h"'],
                                                   'SUNDIALS_PACKAGE_VERSION')
    retcode, sundials_version = conf.TryRun(sundials_version_source, '.cpp')
    if retcode == 0:
        config_error("Failed to determine Sundials version.")

    # Ignore the minor version, e.g. 2.4.x -> 2.4
    env['sundials_version'] = '.'.join(sundials_version.strip().split('.')[:2])
    print """INFO: Using Sundials version %s""" % sundials_version.strip()

env = conf.Finish()


# Python 2 Package Settings
cython_min_version = LooseVersion('0.17')
env['install_python2_action'] = ''
if env['python_package'] in ('full','default','new'):
    # Check for Cython:
    try:
        import Cython
        cython_version = LooseVersion(Cython.__version__)
        print 'INFO: Found Cython version {0} for Python 2.x.'.format(cython_version)
    except ImportError:
        cython_version = LooseVersion('0.0.0')

    if cython_version >= cython_min_version:
        have_cython2 = True
    else:
        message = ("Cython not found or incompatible version: "
                   "Found {0} but {1} or newer is required".format(cython_version, cython_min_version))
        if env['python_package'] == 'new':
            print("ERROR: " + message)
            sys.exit(1)
        else:
            have_cython2 = False
            print ("WARNING: " + message)

    # Test to see if we can import the specified array module
    warnNoPython = False
    if env['python_array_home']:
        sys.path.append(env['python_array_home'])
    try:
        np = __import__(env['python_array'])
        try:
            env['python_array_include'] = np.get_include()
        except AttributeError:
            print """WARNING: Couldn't find include directory for Python array package"""
            env['python_array_include'] = ''

        if env['python_package'] == 'default':
            if have_cython2:
                env['python_package'] = 'new'
            else:
                env['python_package'] = 'full'
        package_desc = 'new' if env['python_package'] == 'new' else 'legacy'
        print """INFO: Building the %s Python package using %s.""" % (package_desc, env['python_array'])
    except ImportError:
        if env['python_package'] in ('full', 'new'):
            print ("""ERROR: Unable to find the array package """
                   """'%s' required by the Python package.""" % env['python_array'])
            sys.exit(1)
        else:
            print ("""WARNING: Not building the Python package """
                   """ because the array package '%s' could not be found.""" % env['python_array'])
            warnNoPython = True
            env['python_package'] = 'minimal'

    try:
        import site
        env['python_usersitepackages'] = site.getusersitepackages()
    except AttributeError: # getusersitepackages is only in Python 2.7+
        env['python_usersitepackages'] = '<user site-packages directory>'

    # Check for 3to2 if we're building the "new" Python module
    # See http://pypi.python.org/pypi/3to2
    if env['python_package'] == 'new':
        try:
            ret = getCommandOutput('3to2','-l')
        except OSError:
            ret = ''
        if 'print' in ret:
            env['python_convert_examples'] = True
        else:
            env['python_convert_examples'] = False
            print """WARNING: Couldn't find '3to2'. Python examples will not work correctly."""

    if env['python_package'] == 'full':
        print ("WARNING: The 'python_package=full' option is deprecated. "
               "This legacy Python package will be removed in Cantera 2.2. "
               "The new Python package may be build using 'python_package=new'.")

else:
    warnNoPython = False
    env['python_array_include'] = ''
    env['python_module_loc'] = ''

# Python 3 Package Settings
if env['python3_package'] in ('y', 'default'):
    # See if we can execute the Python 3 interpreter
    try:
        script = '\n'.join(("from distutils.sysconfig import *",
                            "import site",
                            "print(get_python_version())",
                            "try:",
                            "    print(site.getusersitepackages())",
                            "except AttributeError:",
                            "    print(site.USER_SITE)",
                            "try:",
                            "    import Cython",
                            "    print(Cython.__version__)",
                            "except ImportError:",
                            "    print('0.0.0')"))
        info = getCommandOutput(env['python3_cmd'], '-c', script)
        (env['python3_version'],
         env['python3_usersitepackages'],
         cython_version) = info.splitlines()[-3:]
    except OSError:
        info = False

    if not info:
        if env['python3_package'] == 'default':
            print ('INFO: Not building the Python 3 package because the Python '
                   '3 interpreter %r could not be found' % env['python3_cmd'])
            env['python3_package'] = 'n'
        else:
            print ('ERROR: Could not execute the Python 3 interpreter %r' %
                   env['python3_cmd'])
            sys.exit(1)
    elif cython_version < cython_min_version:
        message = ("Cython package for Python 3 not found or incompatible version: "
                   "Found {0} but {1} or newer is required".format(cython_version, cython_min_version))
        if env['python3_package'] == 'new':
            print("ERROR: " + message)
            sys.exit(1)
        else:
            print ("WARNING: " + message)
    else:
        print 'INFO: Using Cython version {0} for Python 3.x.'.format(cython_version)
        env['python3_package'] = 'y'


# Matlab Toolbox settings
if env['matlab_path'] != '' and env['matlab_toolbox'] == 'default':
    env['matlab_toolbox'] = 'y'

if env['matlab_toolbox'] == 'y':
    matPath = env['matlab_path']
    if matPath == '':
        print """ERROR: Unable to build the Matlab toolbox because 'matlab_path' has not been set."""
        sys.exit(1)

    if not (os.path.isdir(matPath) and
            os.path.isdir(pjoin(matPath, 'extern'))):
        print """ERROR: Path set for 'matlab_path' is not correct."""
        print """ERROR: Path was: '%s'""" % matPath
        sys.exit(1)


# Sundials Settings
if env['use_sundials'] == 'default':
    if env['HAS_SUNDIALS']:
        env['use_sundials'] = 'y'
    else:
        print "INFO: Sundials was not found. Building with minimal ODE solver capabilities."
        env['use_sundials'] = 'n'
elif env['use_sundials'] == 'y' and not env['HAS_SUNDIALS']:
    config_error("Unable to find Sundials headers and / or libraries.")
elif env['use_sundials'] == 'y' and env['sundials_version'] not in ('2.2','2.3','2.4','2.5'):
    print """ERROR: Sundials version %r is not supported.""" % env['sundials_version']
    sys.exit(1)

# Deprecation warnings for old Sundials versions
if env.get('sundials_version') in ('2.2', '2.3'):
    print 'WARNING: Support for Sundials %s is deprecated and will be removed.' % env['sundials_version']
    print 'WARNING: Upgrading to Sundials 2.5 is strongly recommended.'

# Deprecation warnings for numarray and numeric
if env.get('python_array') in ('numarray', 'numeric'):
    print 'WARNING: Support for "%s" is deprecated and will be removed.' % env['python_array']
    print 'WARNING: Upgrading to the "numpy" package is strongly recommended.'

# **********************************************
# *** Set additional configuration variables ***
# **********************************************

# Directories where things will be after actually being installed. These
# variables are the ones that are used to populate header files, scripts, etc.
env['ct_installroot'] = env['prefix']
env['ct_libdir'] = pjoin(env['prefix'], 'lib')
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
    env['ct_mandir'] = pjoin(env['prefix'], 'man', 'man1')
    env['ct_matlab_dir'] = pjoin(env['prefix'], 'lib', 'cantera', 'matlab', 'toolbox')

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
    env['prefix'] = '/usr/share/cantera'
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
    for k in ('python_prefix', 'python3_prefix'):
        if env[k]:
            env[k] = pjoin(os.getcwd(), env['stage_dir'],
                           stripDrive(env[k]).strip('/\\'))
        else:
            env[k] = pjoin(os.getcwd(), env['stage_dir'])
else:
    instRoot = env['prefix']

if env['layout'] == 'debian':
    base = pjoin(os.getcwd(), 'debian')

    env['inst_libdir'] = pjoin(base, 'cantera-dev', 'usr', 'lib')
    env['inst_incdir'] = pjoin(base, 'cantera-dev', 'usr', 'include', 'cantera')
    env['inst_incroot'] = pjoin(base, 'cantera-dev', 'usr' 'include')

    env['inst_bindir'] = pjoin(base, 'cantera-common', 'usr', 'bin')
    env['inst_datadir'] = pjoin(base, 'cantera-common', 'usr', 'share', 'cantera', 'data')
    env['inst_docdir'] = pjoin(base, 'cantera-common', 'usr', 'share', 'cantera', 'doc')
    env['inst_sampledir'] = pjoin(base, 'cantera-common', 'usr', 'share', 'cantera', 'samples')
    env['inst_mandir'] = pjoin(base, 'cantera-common', 'usr', 'share', 'man', 'man1')

    env['inst_matlab_dir'] = pjoin(base, 'cantera-matlab',
                                   'usr', 'lib', 'cantera', 'matlab', 'toolbox')

    env['inst_python_bindir'] = pjoin(base, 'cantera-python', 'usr', 'bin')
    env['python_prefix'] = pjoin(base, 'cantera-python', 'usr')
    env['ct_datadir'] = '/usr/share/cantera/data'
else:
    env['inst_libdir'] = pjoin(instRoot, 'lib')
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
        env['inst_matlab_dir'] = pjoin(instRoot, 'lib', 'cantera', 'matlab', 'toolbox')
        env['inst_datadir'] = pjoin(instRoot, 'share', 'cantera', 'data')
        env['inst_sampledir'] = pjoin(instRoot, 'share', 'cantera', 'samples')
        env['inst_docdir'] = pjoin(instRoot, 'share', 'cantera', 'doc')
        env['inst_mandir'] = pjoin(instRoot, 'man', 'man1')

# **************************************
# *** Set options needed in config.h ***
# **************************************

configh = {'CANTERA_VERSION': quoted(env['cantera_version'])}

# Conditional defines
def cdefine(definevar, configvar, comp=True, value=1):
    if env.get(configvar) == comp:
        configh[definevar] = value
    else:
        configh[definevar] = None

cdefine('DEBUG_MODE', 'debug_verbose')

# Need to test all of these to see what platform.system() returns
configh['SOLARIS'] = 1 if env['OS'] == 'Solaris' else None
configh['DARWIN'] = 1 if env['OS'] == 'Darwin' else None
cdefine('NEEDS_GENERIC_TEMPL_STATIC_DECL', 'OS', 'Solaris')

cdefine('HAS_NUMPY', 'python_array', 'numpy')
cdefine('HAS_NUMARRAY', 'python_array', 'numarray')
cdefine('HAS_NUMERIC', 'python_array', 'numeric')
if env['python_package'] == 'none' and env['python3_package'] == 'n':
    configh['HAS_NO_PYTHON'] = 1
else:
    configh['HAS_NO_PYTHON'] = None

cdefine('HAS_SUNDIALS', 'use_sundials', 'y')
if env['use_sundials'] == 'y':
    configh['SUNDIALS_VERSION'] = env['sundials_version'].replace('.','')
else:
    configh['SUNDIALS_VERSION'] = 0

cdefine('H298MODIFY_CAPABILITY', 'with_h298modify_capability')
cdefine('WITH_HTML_LOGS', 'with_html_log_files')
if env['with_html_log_files']:
    print 'WARNING: Support for HTML log files is deprecated and will be removed in Cantera 2.2.'

cdefine('LAPACK_FTN_STRING_LEN_AT_END', 'lapack_ftn_string_len_at_end')
cdefine('LAPACK_FTN_TRAILING_UNDERSCORE', 'lapack_ftn_trailing_underscore')
cdefine('FTN_TRAILING_UNDERSCORE', 'lapack_ftn_trailing_underscore')
cdefine('LAPACK_NAMES_LOWERCASE', 'lapack_names', 'lower')

configh['RXNPATH_FONT'] = quoted(env['rpfont'])
cdefine('THREAD_SAFE_CANTERA', 'build_thread_safe')

if not env['HAS_MATH_H_ERF']:
    if env['HAS_BOOST_MATH']:
        configh['USE_BOOST_MATH'] = 1
    else:
        config_error("Couldn't find 'erf' in either <math.h> or Boost.Math")
else:
    configh['USE_BOOST_MATH'] = None

config_h = env.Command('include/cantera/base/config.h',
                       'include/cantera/base/config.h.in',
                       ConfigBuilder(configh))
env.AlwaysBuild(config_h)
env['config_h_target'] = config_h

env['boost_libs'] = []
if env['build_thread_safe']:
    env['use_boost_libs'] = True
    if env['CC'] != 'cl':
        env['boost_libs'].append(env['boost_thread_lib'])
else:
    env['use_boost_libs'] = False

# *********************
# *** Build Cantera ***
# *********************

# Some options to speed up SCons
env.SetOption('max_drift', 2)
env.SetOption('implicit_cache', True)

buildTargets = []
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
                      '$python_cmd interfaces/python/ctml_writer.py $SOURCE $TARGET'))


# Copy input files which are not present as cti:
for xml in mglob(env, 'data/inputs', 'xml'):
    dest = pjoin('build','data',xml.name)
    if xml.name not in convertedInputFiles:
        build(env.Command(dest, xml.path, Copy('$TARGET', '$SOURCE')))

if addInstallActions:
    if env['VERBOSE']:
       print 'INFO 2: addInstallActions script started'

    # Put headers in place
    headerBase = 'include/cantera'
    install(env.RecursiveInstall, '$inst_incdir', 'include/cantera')

    # Make symlinks to replicate old header directory structure
    if env['legacy_headers']:
        print ('WARNING: Installation of legacy headers is deprecated. '
               'This installation option will be removed in Cantera 2.2.')
        install(env.Command, pjoin('$inst_incdir', 'kernel'), [], Mkdir("$TARGET"))
        install('$inst_incdir', 'platform/legacy/Cantera_legacy.h')

        if env['OS'] == 'Windows':
            cmd = Copy("$TARGET", "$SOURCE")
        else:
            def RelativeSymlink(target, source, env):
                if os.path.exists(target[0].path):
                    os.remove(target[0].path)
                srcpath = psplit(source[0].abspath)
                tgtpath = psplit(target[0].abspath)
                nCommon = max(i for i,(dir1,dir2) in enumerate(zip(srcpath, tgtpath))
                              if dir1 == dir2)
                relsrc = os.sep.join(['..'] + srcpath[nCommon-1:])
                os.symlink(relsrc, target[0].abspath)

            cmd = RelativeSymlink

        for name in os.listdir('include/cantera'):
            if not os.path.isdir(pjoin('include/cantera', name)):
                continue
            for filename in os.listdir(pjoin('include/cantera', name)):
                if not filename.endswith('.h'):
                    continue
                headerdir = pjoin(instRoot, 'include', 'cantera')
                install(env.Command, pjoin(headerdir, 'kernel', filename),
                        pjoin(headerdir, name, filename), cmd)

    # Data files
    install('$inst_datadir', mglob(env, pjoin('build','data'), 'cti', 'xml'))

    # Converter scripts
    pyExt = '.py' if env['OS'] == 'Windows' else ''
    install(env.InstallAs,
            pjoin('$inst_bindir','ck2cti%s' % pyExt),
            'interfaces/python/ck2cti.py')
    install(env.InstallAs,
            pjoin('$inst_bindir','ctml_writer%s' % pyExt),
            'interfaces/python/ctml_writer.py')

    # Copy external libaries for Windows installations
    if env['CC'] == 'cl' and env['use_boost_libs']:
        boost_suffix = '-vc%s-mt-%s.lib' % (env['MSVC_VERSION'].replace('.',''),
                                        env['BOOST_LIB_VERSION'])
        for lib in env['boost_windows_libs'].split(','):
            install('$inst_libdir', pjoin('$boost_lib_dir',
                                          'libboost_{0}{1}'.format(lib, boost_suffix)))

    # Copy sundials library and header files
    if env['install_sundials']:
        for subdir in ['cvode','cvodes','ida','idas','kinsol','nvector','sundials']:
            if os.path.exists(pjoin(env['sundials_include'], subdir)):
                install(env.RecursiveInstall, pjoin('$inst_incdir', '..', subdir),
                        pjoin(env['sundials_include'], subdir))
        if os.path.exists(env['sundials_license']):
            install(pjoin('$inst_incdir', '..', 'sundials'), env['sundials_license'])
        libprefix = '' if os.name == 'nt' else 'lib'
        install('$inst_libdir', mglob(env, env['sundials_libdir'],
                                      '^{0}sundials_*'.format(libprefix)))


### List of libraries needed to link to Cantera ###
linkLibs = ['cantera']
linkSharedLibs = ['cantera_shared']

if env['use_sundials'] == 'y':
    env['sundials_libs'] = ['sundials_cvodes', 'sundials_ida', 'sundials_nvecserial']
    linkLibs.extend(('sundials_cvodes', 'sundials_ida', 'sundials_nvecserial'))
    linkSharedLibs.extend(('sundials_cvodes', 'sundials_ida', 'sundials_nvecserial'))
else:
    env['sundials_libs'] = []
    if not env['single_library']:
        linkLibs.extend(['cvode'])
        linkSharedLibs.extend(['cvode_shared'])
    #print 'linkLibs = ', linkLibs

if not env['single_library']:
    linkLibs.append('ctmath')
    linkSharedLibs.append('ctmath_shared')

    # Add execstream to the link line
    linkLibs.append('execstream')
    linkSharedLibs.append('execstream_shared')

#  Add lapack and blas to the link line
if env['blas_lapack_libs']:
    linkLibs.extend(env['blas_lapack_libs'])
    linkSharedLibs.extend(env['blas_lapack_libs'])
elif not env['single_library']:
    linkLibs.extend(('ctlapack', 'ctblas'))
    linkSharedLibs.extend(('ctlapack_shared', 'ctblas_shared'))

if not env['build_with_f2c']:
    linkLibs.extend(env['FORTRANSYSLIBS'])
    linkSharedLibs.append(env['FORTRANSYSLIBS'])
elif not env['single_library']:
    # Add the f2c library when f2c is requested
    linkLibs.append('ctf2c')
    linkSharedLibs.append('ctf2c_shared')

linkLibs.extend(env['boost_libs'])
linkSharedLibs.extend(env['boost_libs'])

# Store the list of needed static link libraries in the environment
env['cantera_libs'] = linkLibs
env['cantera_shared_libs'] = linkSharedLibs
if env['renamed_shared_libraries'] == False :
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

if env['python_package'] in ('full','minimal'):
    VariantDir('build/src/python', 'src/python', duplicate=0)
    SConscript('build/src/python/SConscript')

if env['python3_package'] == 'y' or env['python_package'] == 'new':
    SConscript('interfaces/cython/SConscript')

SConscript('build/src/apps/SConscript')

if env['OS'] != 'Windows':
    VariantDir('build/platform', 'platform/posix', duplicate=0)
    SConscript('build/platform/SConscript')

if env['matlab_toolbox'] == 'y':
    SConscript('build/src/matlab/SConscript')

if env['doxygen_docs'] or env['sphinx_docs']:
    SConscript('doc/SConscript')

if 'samples' in COMMAND_LINE_TARGETS or addInstallActions:
    VariantDir('build/samples', 'samples', duplicate=0)
    sampledir_excludes = ['ct2ctml', '\\.o$', '^~$', 'xml$', '\\.in',
                          'SConscript', 'Makefile.am']
    SConscript('build/samples/cxx/SConscript')

    # Install C++ samples
    install(env.RecursiveInstall, pjoin('$inst_sampledir', 'cxx'),
            'samples/cxx', exclude=sampledir_excludes)

    if env['f90_interface'] == 'y':
        SConscript('build/samples/f77/SConscript')
        SConscript('build/samples/f90/SConscript')

        # install F90 / F77 samples
        install(env.RecursiveInstall, pjoin('$inst_sampledir', 'f77'),
                'samples/f77', sampledir_excludes)
        install(env.RecursiveInstall, pjoin('$inst_sampledir', 'f90'),
                'samples/f90', sampledir_excludes)

### Meta-targets ###
build_samples = Alias('samples', sampleTargets)

def postBuildMessage(target, source, env):
    print "*******************************************************"
    print "Compilation completed successfully.\n"
    print "- To run the test suite, type 'scons test'."
    if os.name == 'nt':
        print "- To install, type 'scons install'."
        print "- To create a Windows MSI installer, type 'scons msi'"
    else:
        print "- To install, type '[sudo] scons install'."
    print "*******************************************************"

finish_build = env.Command('finish_build', [], postBuildMessage)
env.Depends(finish_build, buildTargets)
build_cantera = Alias('build', finish_build)

Default('build')

def postInstallMessage(target, source, env):
    if env['python_package'] == 'none':
        env['python_module_loc'] = 'NONE'

    if env['python3_package'] == 'y':
        env['python3_example_loc'] = pjoin(env['python3_module_loc'], 'cantera', 'examples')

    env['python_example_loc'] = pjoin(env['python_module_loc'], 'cantera', 'examples')

    env['matlab_sample_loc'] = pjoin(env['ct_sampledir'], 'matlab')
    env['matlab_ctpath_loc'] = pjoin(env['ct_matlab_dir'], 'ctpath.m')
    print """
Cantera has been successfully installed.

File locations:

  applications                %(ct_bindir)s
  library files               %(ct_libdir)s
  C++ headers                 %(ct_incroot)s
  samples                     %(ct_sampledir)s
  data files                  %(ct_datadir)s""" % env,

    if env['python_package'] == 'full':
        print """
  Python 2 package (Cantera)  %(python_module_loc)s""" % env,
    elif warnNoPython:
        print """
    #################################################################
     WARNING: the Cantera Python package was not installed because a
     suitable array package (e.g. numpy) could not be found.
    #################################################################"""

    if env['python_package'] == 'new':
        print """
  Python 2 package (cantera)  %(python_module_loc)s
  Python 2 samples            %(python_example_loc)s""" % env,

    if env['python3_package'] == 'y':
        print """
  Python 3 package (cantera)  %(python3_module_loc)s
  Python 3 samples            %(python3_example_loc)s""" % env,

    if env['matlab_toolbox'] == 'y':
        print """
  Matlab toolbox              %(ct_matlab_dir)s
  Matlab samples              %(matlab_sample_loc)s

An m-file to set the correct matlab path for Cantera is at:

  %(matlab_ctpath_loc)s
    """ % env,

    if os.name != 'nt':
        print """
  setup script                %(ct_bindir)s/setup_cantera

The setup script configures the environment for Cantera. It is recommended that
you run this script by typing:

  source %(ct_bindir)s/setup_cantera

before using Cantera, or else include its contents in your shell login script.
    """ % env
    else:
        print ''

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

# Files installed by the Python installer(s)
pyFiles = ['build/python-installed-files.txt',
           'build/python2-installed-files.txt',
           'build/python3-installed-files.txt']

for filename in pyFiles:
    if os.path.exists(filename):
        allfiles.extend([File(f.strip()) for f in open(filename).readlines()])

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
                                  x64=env['TARGET_ARCH']=='amd64',
                                  includeMatlab=env['matlab_toolbox']=='y')
        wxs.make_wxs(str(target[0]))

    wxs_target = env.Command(pjoin('build', 'wix', 'cantera.wxs'),
                             [], build_wxs)
    env.AlwaysBuild(wxs_target)

    env.Append(WIXLIGHTFLAGS=['-ext', 'WixUIExtension'])
    msi_target = env.WiX('cantera.msi',
                         [pjoin('build', 'wix', 'cantera.wxs')])
    env.Depends(wxs_target, installTargets)
    env.Depends(msi_target, wxs_target)
    build_msi = Alias('msi', msi_target)

### Tests ###
if any(target.startswith('test') for target in COMMAND_LINE_TARGETS):
    env['testNames'] = []

    # Tests written using the gtest framework, the Python unittest module,
    # or the Matlab xunit package.
    VariantDir('build/test', 'test', duplicate=0)
    SConscript('build/test/SConscript')

    # Regression tests
    SConscript('test_problems/SConscript')

    if 'test-help' in COMMAND_LINE_TARGETS:
        print '\n*** Available tests ***\n'
        for name in env['testNames']:
            print 'test-%s' % name
        sys.exit(0)
