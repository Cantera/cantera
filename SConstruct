from buildutils import *
import platform, sys, os

env = Environment()

# **************************************
# *** Read user-configurable options ***
# **************************************

opts = Variables('cantera.conf')
opts.AddVariables(
    PathVariable('prefix', 'Where to install Cantera',
                 '/usr/local', PathVariable.PathIsDirCreate),
    EnumVariable('python_package', 'build python package?', 'default',
                 ('full', 'minimal', 'none', 'default')),
    PathVariable('python_cmd', 'Path to the python interpreter', sys.executable),
    EnumVariable('matlab_toolbox', '', 'n', ('y', 'n', 'default')),
    BoolVariable('f90_interface', 'Build Fortran90 interface?', False),
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
    ('blas_lapack_libs', '', ''), # '-llapack -lblas' or '-llapack -lf77blas -lcblas -latlas' etc.
    ('blas_lapack_dir', '', ''), # '/usr/lib/lapack' etc
    EnumVariable('lapack_names', '', 'lower', ('lower','upper')),
    BoolVariable('lapack_ftn_trailing_underscore', '', True),
    BoolVariable('lapack_ftn_string_len_at_end', '', True),
    ('bitcompile', '', ''), # '32' or '64'
    ('CXX', '', env['CXX']),
    ('CC', '', env['CC']),
    ('CXXFLAGS', '', '-O3 -Wall'),
    BoolVariable('build_thread_safe', '', False),
    BoolVariable('build_with_f2c', '', True),
    ('F77', '', env['F77']),
    ('F77FLAGS', '', '-O3'),
    ('F90', '', env['F90']),
    ('F90FLAGS', '', '-O3'),
    ('install_bin', '', 'config/install-sh'),
    ('graphvisdir', '' ,''),
    ('ct_shared_lib', '', 'clib'),
    ('rpfont', '', 'Helvetica'),
    ('cantera_version', '', '1.8.x')
# These variables shouldn't be necessary any more...
#    ('cxx_ext', '', 'cpp'),
#    ('f77_ext', '', 'f'),
#    ('f90_ext', '', 'f90'),
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

# Additional options that apply only if building the full python package
if env['python_package'] in ('full', 'default'):
    opts.AddVariables(
        EnumVariable('python_array', 'Which Python array package to use',
                     'numpy', ('numpy', 'numarray', 'numeric')),
        BoolVariable('set_python_site_package_topdir', '', False),
        PathVariable('python_site_package_topdir', '', '/usr/local'),
        PathVariable('python_array_home',
                     'Location for array package (e.g. if installed with --home)',
                     None, PathVariable.PathAccept),
        PathVariable('cantera_python_home', 'where to install the python package',
                     None, PathVariable.PathAccept),
        )

# Options that apply only if building the Matlab interface
if env['matlab_toolbox'] != 'n':
    opts.AddVariables(
        PathVariable('matlab_cmd', 'Path to the matlab executable',
                     'default', PathVariable.PathAccept)
        )

# Options that apply only if building the Fortran interface
if env['f90_interface']:
    opts.AddVariables(
        PathVariable('f90', 'Fortran compiler', 'gfortran'),
        ('f90flags', '', '-O3')
        )

# Extra options for Sundials
if env['use_sundials'] != 'n':
    opts.AddVariables(
        EnumVariable('sundials_version' ,'', '2.4', ('2.2','2.3','2.4')),
        PathVariable('sundials_include' ,'', ''))

# Extra options for Boost.Thread
if env['build_thread_safe']:
    opts.AddVariables(
        PathVariable('boost_inc_dir', '', '/usr/include/'),
        PathVariable('boost_lib_dir', '', '/usr/lib/'),
        ('boost_thread_lib', '', 'boost_thread'))

opts.Update(env)
opts.Save('cantera.conf', env)

# ********************************************
# *** Configure system-specific properties ***
# ********************************************
env['OS'] = platform.system()

if env['F90'] == 'gfortran':
    env['FORTRANMODDIRPREFIX'] = '-J'
elif env['F90'] == 'g95':
    env['FORTRANMODDIRPREFIX'] = '-fmod='

env['FORTRANMODDIR'] = '${TARGET.dir}'

conf = Configure(env)

env['HAS_SSTREAM'] = conf.CheckCXXHeader('sstream', '<>')

env = conf.Finish()


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
    env['blas_lapack_libs'] = '-lctlapack -lctblas'

if env['use_sundials'] == 'y' and env['sundials_include']:
    env.Append(CPPPATH=env['sundials_include'])

# *********************
# *** Build Cantera ***
# *********************

# Put headers in place
for header in mglob(env, 'Cantera/cxx/include', 'h'):
    env.Command('build/include/cantera/%s' % header.name, header,
                Copy('$TARGET', '$SOURCE'))

for header in mglob(env, 'Cantera/clib/src', 'h'):
    env.Command('build/include/cantera/clib/%s' % header.name, header,
                Copy('$TARGET', '$SOURCE'))


build = 'build'
env.SConsignFile()
env.Append(CPPPATH=[Dir(os.getcwd()),
                    Dir('build/include/cantera/kernel'),
                    Dir('build/include/cantera')])

Export('env', 'build')

VariantDir('build/ext', 'ext', duplicate=0)
SConscript('build/ext/SConscript')

VariantDir('build/kernel', 'Cantera/src', duplicate=0)
SConscript('build/kernel/SConscript')

VariantDir('build/interfaces/', 'Cantera', duplicate=0)
SConscript('build/interfaces/SConscript')

VariantDir('build/interfaces/fortran/', 'Cantera/fortran', duplicate=1)
SConscript('build/interfaces/fortran/SConscript')
