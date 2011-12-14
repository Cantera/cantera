from buildutils import *
import subst
import platform, sys, os

env = Environment(tools = ['default', 'textfile'])
env.AddMethod(RecursiveInstall)
subst.TOOL_SUBST(env)

# ******************************************************
# *** Set system-dependent defaults for some options ***
# ******************************************************

if os.name == 'posix':
    defaultPrefix = '/usr/local'
elif os.name == 'nt':
    defaultPrefix = os.environ['ProgramFiles']

# **************************************
# *** Read user-configurable options ***
# **************************************

opts = Variables('cantera.conf')
opts.AddVariables(
    PathVariable('prefix', 'Where to install Cantera',
                 defaultPrefix, PathVariable.PathIsDirCreate),
    EnumVariable('python_package', 'build python package?', 'default',
                 ('full', 'minimal', 'none','default')),
    PathVariable('python_cmd', 'Path to the python interpreter', sys.executable),
    EnumVariable('python_array', 'Which Python array package to use',
                 'numpy', ('numpy', 'numarray', 'numeric')),
    PathVariable('python_array_home',
                 'Location for array package (e.g. if installed with --home)',
                 '', PathVariable.PathAccept),
    PathVariable('cantera_python_home', 'where to install the python package',
                 '', PathVariable.PathAccept),
    EnumVariable('matlab_toolbox', '', 'default', ('y', 'n', 'default')),
    PathVariable('matlab_cmd', 'Path to the matlab executable',
                 'matlab', PathVariable.PathAccept),
    EnumVariable('f90_interface', 'Build Fortran90 interface?', 'default', ('y', 'n', 'default')),
    PathVariable('F90', 'Fortran compiler',
                 '', PathVariable.PathAccept),
    ('purify', '', ''),
    ('user_src_dir', '', 'Cantera/user'),
    BoolVariable('debug', '', False), # ok
    BoolVariable('with_lattice_solid', '', True), # ok
    BoolVariable('with_metal', '', True), # ok
    BoolVariable('with_stoich_substance', '', True), # ok
    BoolVariable('with_semiconductor', '', True), # ??
    BoolVariable('with_adsorbate', '', True),
    BoolVariable('with_spectra', '', True),
    BoolVariable('with_pure_fluids', '', True),
    BoolVariable('with_ideal_solutions', '', True), # ok
    BoolVariable('with_electrolytes', '', True), # ok
    BoolVariable('with_prime', '', False), # ok
    BoolVariable('with_h298modify_capability', '', False), # ok
    BoolVariable('enable_ck', '', True),
    BoolVariable('with_kinetics', '', True),
    BoolVariable('with_hetero_kinetics', '', True),
    BoolVariable('with_reaction_paths', '', True),
    BoolVariable('with_vcsnonideal', '', False),
    BoolVariable('enable_transport', '', True),
    BoolVariable('enable_equil', '', True),
    BoolVariable('enable_reactors', '', True),
    BoolVariable('enable_flow1d', '', True),
    BoolVariable('enable_solvers', '', True),
    BoolVariable('enable_rxnpath', '', True),
    BoolVariable('enable_tpx', '', True),
    BoolVariable('with_html_log_files', '', True),
    EnumVariable('use_sundials', '', 'default', ('default', 'y', 'n')),
    EnumVariable('sundials_version' ,'', '2.4', ('2.2','2.3','2.4')),
    PathVariable('sundials_include' ,'', ''),
    PathVariable('sundials_libdir', '', ''),
    ('blas_lapack_libs', '', ''), # 'lapack,blas' or 'lapack,f77blas,cblas,atlas' etc.
    ('blas_lapack_dir', '', ''), # '/usr/lib/lapack' etc
    EnumVariable('lapack_names', '', 'lower', ('lower','upper')),
    BoolVariable('lapack_ftn_trailing_underscore', '', True),
    BoolVariable('lapack_ftn_string_len_at_end', '', True),
    ('bitcompile', '', ''), # '32' or '64'
    ('CXX', '', env['CXX']),
    ('CC', '', env['CC']),
    ('CXXFLAGS', '', '-O3 -Wall'),
    BoolVariable('build_thread_safe', '', False),
    PathVariable('boost_inc_dir', '', '/usr/include/'),
    PathVariable('boost_lib_dir', '', '/usr/lib/'),
    ('boost_thread_lib', '', 'boost_thread'),
    BoolVariable('build_with_f2c', '', True),
    ('F77', '', env['F77']),
    ('F77FLAGS', '', '-O3'),
    ('F90FLAGS', '', '-O3'),
    ('install_bin', '', 'config/install-sh'),
    ('graphvisdir', '' ,''),
    ('ct_shared_lib', '', 'clib'),
    ('rpfont', '', 'Helvetica'),
    ('cantera_version', '', '1.8.x')
# These variables shouldn't be necessary any more...
#    ('exe_ext', '', ''),
#    ('lcxx_end_libs', '-lm'),
#    ('pic', '', '-fPIC'),
#    ('shared', '', '-dynamic'),
#    ('lfort_flags', '', '-L/usr/local/lib'),
#    ('AR', '', env['AR']),
#    ('ARFLAGS', '', env['ARFLAGS']), # ('archive', '', 'ar ruv'),
#    ('ranlib', '', 'ranlib'),
    )

opts.Update(env)
opts.Save('cantera.conf', env)

# ********************************************
# *** Configure system-specific properties ***
# ********************************************
env['OS'] = platform.system()

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
    elif env['f90_interface'] == 'y':
        print "ERROR: Couldn't find a suitable Fortran compiler to build the Fortran 90 interface"
        sys.exit(1)
    else:
        print "INFO: Skipping compilation of the Fortran 90 interface."

if env['F90'] == 'gfortran':
    env['FORTRANMODDIRPREFIX'] = '-J'
elif env['F90'] == 'g95':
    env['FORTRANMODDIRPREFIX'] = '-fmod='
elif env['F90'] == 'ifort':
    env['FORTRANMODDIRPREFIX'] = '-module '

env['FORTRANMODDIR'] = '${TARGET.dir}'

conf = Configure(env)

env['HAS_SSTREAM'] = conf.CheckCXXHeader('sstream', '<>')

env = conf.Finish()

if env['cantera_python_home'] == '' and env['prefix'] != defaultPrefix:
    env['cantera_python_home'] = env['prefix']

if env['python_package'] in ('full','default'):
    # Test to see if we can import the specified array module
    warnNoPython = False
    if env['python_array_home']:
        sys.path.append(env['python_array_home'])
    try:
        __import__(env['python_array'])
        print """INFO: Building the full Python package using %s.""" % env['python_array']
        env['python_package'] = 'full'
    except ImportError:
        if env['python_package'] == 'full':
            print ("""ERROR: Unable to find the array package """
                   """'%s' required by the full Python package.""" % env['python_array'])
            sys.exit(1)
        else:
            print ("""WARNING: Not building the full Python package """
                   """ because the array package '%s' could not be found.""" % env['python_array'])
            warnNoPython = True
            env['python_package'] = 'minimal'


if env['matlab_toolbox'] == 'y' and which(env['matlab_cmd']) is None:
    print """ERROR: Unable to find the Matlab executable '%s'""" % env['matlab_cmd']
    sys.exit(1)
elif env['matlab_toolbox'] == 'default':
    cmd = which(env['matlab_cmd'])
    if cmd is not None:
        env['matlab_toolbox'] = 'y'
        print """INFO: Building the Matlab toolbox using '%s'""" % cmd
    else:
        print """INFO: Skipping compilation of the Matlab toolbox. """


# **************************************
# *** Set options needed in config.h ***
# **************************************

configh = {'CANTERA_VERSION': quoted(env['cantera_version']),
           }

# Conditional defines
def cdefine(definevar, configvar, comp=True, value=1):
    if env.get(configvar) == comp:
        configh[definevar] = value
    else:
        configh[definevar] = None

cdefine('DEBUG_MODE', 'debug')
cdefine('PURIFY_MODE', 'purify')

# Need to test all of these to see what platform.system() returns
configh['SOLARIS'] = 1 if env['OS'] == 'Solaris' else None
configh['DARWIN'] = 1 if env['OS'] == 'Darwin' else None
configh['CYGWIN'] = 1 if env['OS'] == 'Cygwin' else None
configh['WINMSVC'] = 1 if env['OS'] == 'Windows' else None
cdefine('NEEDS_GENERIC_TEMPL_STATIC_DECL', 'OS', 'Solaris')

cdefine('HAS_NUMPY', 'python_array', 'numpy')
cdefine('HAS_NUMARRAY', 'python_array', 'numarray')
cdefine('HAS_NUMERIC', 'python_array', 'numeric')
cdefine('HAS_NO_PYTHON', 'python_package', 'none')
configh['PYTHON_EXE'] = quoted(env['python_cmd']) if env['python_package'] != 'none' else None

cdefine('HAS_SUNDIALS', 'use_sundials', 'y')
if env['use_sundials']:
    cdefine('SUNDIALS_VERSION_22', 'sundials_version', '2.2')
    cdefine('SUNDIALS_VERSION_23', 'sundials_version', '2.3')
    cdefine('SUNDIALS_VERSION_24', 'sundials_version', '2.4')

cdefine('WITH_ELECTROLYTES', 'with_electrolytes')
cdefine('WITH_IDEAL_SOLUTIONS', 'with_ideal_solutions')
cdefine('WITH_LATTICE_SOLID', 'with_lattice_solid')
cdefine('WITH_METAL', 'with_metal')
cdefine('WITH_STOICH_SUBSTANCE', 'with_stoich_substance')
cdefine('WITH_SEMICONDUCTOR', 'with_semiconductor')
cdefine('WITH_PRIME', 'with_prime')
cdefine('H298MODIFY_CAPABILITY', 'with_n298modify_capability')
cdefine('WITH_PURE_FLUIDS', 'with_pure_fluids')
cdefine('INCL_PURE_FLUIDS', 'with_pure_fluids') # TODO: fix redundancy
cdefine('WITH_HTML_LOGS', 'with_html_log_files')
cdefine('WITH_VCSNONIDEAL', 'with_vcsnonideal')

cdefine('LAPACK_FTN_STRING_LEN_AT_END', 'lapack_ftn_string_len_at_end')
cdefine('LAPACK_FTN_TRAILING_UNDERSCORE', 'lapack_ftn_trailing_underscore')
cdefine('FTN_TRAILING_UNDERSCORE', 'lapack_ftn_trailing_underscore')
cdefine('LAPACK_NAMES_LOWERCASE', 'lapack_names', 'lower')

configh['RXNPATH_FONT'] = quoted(env['rpfont'])
cdefine('THREAD_SAFE_CANTERA', 'build_thread_safe')
cdefine('HAS_SSTREAM', 'HAS_SSTREAM')
configh['CANTERA_DATA'] = quoted(os.path.join(env['prefix'], 'data'))

config_h = env.Command('config.h', 'config.h.in.scons', ConfigBuilder(configh))
env.AlwaysBuild(config_h)

# **********************************************
# *** Set additional configuration variables ***
# **********************************************
if env['blas_lapack_libs'] == '':
    # External BLAS/LAPACK were not given, so we need to compile them
    env['BUILD_BLAS_LAPACK'] = True
    env['blas_lapack_libs'] = ['ctlapack', 'ctblas']
else:
    ens['blas_lapack_libs'] = ','.split(env['blas_lapack_libs'])

if env['use_sundials'] == 'y' and env['sundials_include']:
    env.Append(CPPPATH=env['sundials_include'])
if env['use_sundials'] == 'y' and env['sundials_libdir']:
    env.Append(LIBPATH=env['sundials_libdir'])

env['ct_libdir'] = pjoin(env['prefix'], 'lib')
env['ct_bindir'] = pjoin(env['prefix'], 'bin')
env['ct_incdir'] = pjoin(env['prefix'], 'include', 'cantera')
env['ct_incroot'] = pjoin(env['prefix'], 'include')
env['ct_datadir'] = pjoin(env['prefix'], 'data')
env['ct_demodir'] = pjoin(env['prefix'], 'demos')
env['ct_templdir'] = pjoin(env['prefix'], 'templates')
env['ct_tutdir'] = pjoin(env['prefix'], 'tutorials')
env['ct_docdir'] = pjoin(env['prefix'], 'doc')
env['ct_dir'] = env['prefix']
env['ct_mandir'] = pjoin(env['prefix'], 'man1')
env['ct_matlab_dir'] = pjoin(env['prefix'], 'matlab', 'toolbox')

# *********************
# *** Build Cantera ***
# *********************

buildDir = 'build'
buildTargets = []
installTargets = []
demoTargets = []

env.SConsignFile()

env.Append(CPPPATH=[Dir(os.getcwd()),
                    Dir('build/include/cantera/kernel'),
                    Dir('build/include/cantera'),
                    Dir('build/include')],
           LIBPATH=[Dir('build/lib')],
           CCFLAGS=['-fPIC'],
           FORTRANFLAGS=['-fPIC'],
           F90FLAGS=['-fPIC'])

# Put headers in place
for header in mglob(env, 'Cantera/cxx/include', 'h'):
    header = env.Command('build/include/cantera/%s' % header.name, header,
                         Copy('$TARGET', '$SOURCE'))
    buildTargets.extend(header)
    inst = env.Install('$ct_incdir', header)
    installTargets.extend(inst)

for header in mglob(env, 'Cantera/clib/src', 'h'):
    hcopy = env.Command('build/include/cantera/clib/%s' % header.name, header,
                        Copy('$TARGET', '$SOURCE'))
    buildTargets.append(header)
    inst = env.Install(pjoin('$ct_incdir','clib'), header)
    installTargets.extend(inst)

### List of libraries needed to link to Cantera ###
linkLibs = ['clib','oneD','zeroD','equil','kinetics','transport',
            'thermo','ctnumerics','ctmath','tpx',
            'ctspectra','converters','ctbase']

if env['use_sundials']:
    linkLibs.extend(('sundials_cvodes','sundials_nvecserial'))

linkLibs.extend(env['blas_lapack_libs'])

if env['build_with_f2c']:
    linkLibs.append('ctf2c')
else:
    linkLibs.append('gfortran')

env['cantera_libs'] = linkLibs

configh = env.Command('build/include/cantera/config.h', 'config.h', Copy('$TARGET', '$SOURCE'))
inst = env.Install('$ct_incdir', configh)
installTargets.extend(inst)

# Add targets from the SConscript files in the various subdirectories
Export('env', 'buildDir', 'buildTargets', 'installTargets', 'demoTargets')

VariantDir('build/ext', 'ext', duplicate=0)
SConscript('build/ext/SConscript')

VariantDir('build/kernel', 'Cantera/src', duplicate=0)
SConscript('build/kernel/SConscript')

VariantDir('build/interfaces/clib', 'Cantera/clib', duplicate=0)
SConscript('build/interfaces/clib/SConscript')

VariantDir('build/interfaces/cxx', 'Cantera/cxx', duplicate=0)
SConscript('build/interfaces/cxx/SConscript')

if env['f90_interface'] == 'y':
    VariantDir('build/interfaces/fortran/', 'Cantera/fortran', duplicate=1)
    SConscript('build/interfaces/fortran/SConscript')

if env['python_package'] in ('full','minimal'):
    SConscript('Cantera/python/SConscript')

if env['matlab_toolbox'] == 'y':
    SConscript('Cantera/matlab/SConscript')

VariantDir('build/tools', 'tools', duplicate=0)
SConscript('build/tools/SConscript')

# Data files
inst = env.Install('$ct_datadir', mglob(env, pjoin('data','inputs'), 'cti', 'xml'))
installTargets.extend(inst)

# Install exp3to2.sh (used by some of the tests)
inst = env.Install('$ct_bindir', pjoin('bin', 'exp3to2.sh'))
installTargets.extend(inst)

### Meta-targets ###
build_demos = Alias('demos', demoTargets)

def postBuildMessage(target, source, env):
    print "**************************************************************"
    print "Compiliation complete. Type '[sudo] scons install' to install."
    print "**************************************************************"

finish_build = env.Command('finish_build', [], postBuildMessage)
env.Depends(finish_build, buildTargets)
build_cantera = Alias('build', finish_build)

Default('build')

def postInstallMessage(target, source, env):
    v = sys.version_info
    env['python_module_loc'] = pjoin(
        env['prefix'], 'lib', 'python%i.%i' % v[:2], 'site-packages')

    print """
Cantera has been successfully installed.

File locations:

    applications      %(ct_bindir)s
    library files     %(ct_libdir)s
    C++ headers       %(ct_incdir)s
    demos             %(ct_demodir)s
    data files        %(ct_datadir)s""" % env

    if env['python_package'] == 'full':
        print """
    Python package    %(python_module_loc)s""" % env
    elif warnNoPython:
        print """
    #################################################################
     WARNING: the Cantera Python package was not installed because a
     suitable array package (e.g. numpy) could not be found.
    #################################################################"""

    if env['matlab_toolbox'] == 'y':
        print """
    Matlab toolbox    %(ct_matlab_dir)s
    Matlab demos      %(ct_demodir)s/matlab
    Matlab tutorials  %(ct_tutdir)s/matlab

    An m-file to set the correct matlab path for Cantera is at:

        %(prefix)s/matlab/ctpath.m""" % env

    print """
    setup script      %(ct_bindir)s/setup_cantera

    The setup script configures the environment for Cantera. It is
    recommended that you run this script by typing:

        source %(ct_bindir)s/setup_cantera

    before using Cantera, or else include its contents in your shell
    login script.
    """ % env

finish_install = env.Command('finish_install', [], postInstallMessage)
env.Depends(finish_install, installTargets)
install_cantera = Alias('install', finish_install)
