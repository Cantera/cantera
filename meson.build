project('cantera', 'cpp', 'c',
  version: '4.0.0',
  license: 'BSD-3-Clause',
  default_options: [
    'cpp_std=c++20',
    'c_std=c99',
    'warning_level=2',
    'buildtype=release',
  ],
  meson_version: '>=1.0.0'
)

# Import modules
fs = import('fs')

# Get version components
version_parts = meson.project_version().split('.')
version_major = version_parts[0]
version_minor = version_parts[1]
cantera_short_version = f'@version_major@.@version_minor@'

# Compiler setup
cc = meson.get_compiler('c')
cpp = meson.get_compiler('cpp')

# Platform-specific settings
host_system = host_machine.system()
is_windows = host_system == 'windows'
is_darwin = host_system == 'darwin'

# Additional compiler flags
add_project_arguments('-DFMT_HEADER_ONLY', language: 'cpp')

if cpp.get_id() == 'msvc'
  add_project_arguments('/EHsc', '/utf-8', '/nologo', language: 'cpp')
  add_project_arguments('/D_SCL_SECURE_NO_WARNINGS', language: 'cpp')
  add_project_arguments('/D_CRT_SECURE_NO_WARNINGS', language: 'cpp')
  add_project_arguments('/MD', language: 'cpp')
else
  add_project_arguments('-fcolor-diagnostics', language: 'cpp')
endif

# Thread support
thread_dep = dependency('threads', required: true)

# Math library (needed on some platforms)
m_dep = cc.find_library('m', required: false)

# SUNDIALS dependency with fallback to subproject
sundials_dep = dependency('sundials_cvodes', version: '>=6.0', required: false)
sundials_idas_dep = dependency('sundials_idas', required: false)
sundials_nvecserial_dep = dependency('sundials_nvecserial', required: false)

use_system_sundials = false
if sundials_dep.found() and sundials_idas_dep.found() and sundials_nvecserial_dep.found()
  use_system_sundials = true
  message('Using system SUNDIALS')
else
  message('Using bundled SUNDIALS from ext/sundials')
  # TODO: Build SUNDIALS from submodule
  # For now, we'll add it as a subproject wrap
  sundials_dep = dependency('', required: false)
  sundials_idas_dep = dependency('', required: false)
  sundials_nvecserial_dep = dependency('', required: false)
endif

# Eigen dependency (header-only)
eigen_dep = dependency('eigen3', version: '>=3.4', required: false)
use_system_eigen = false
if eigen_dep.found()
  use_system_eigen = true
  message('Using system Eigen')
else
  message('Using bundled Eigen from ext/eigen')
  # Check if submodule is initialized
  eigen_include_dir = meson.current_source_dir() / 'ext' / 'eigen'
  if not fs.is_dir(eigen_include_dir)
    warning('Eigen submodule not initialized. Run: git submodule update --init ext/eigen')
  endif
  eigen_dep = declare_dependency(
    include_directories: include_directories('ext/eigen')
  )
endif

# fmt dependency with fallback
fmt_dep = dependency('fmt', version: '>=9.1.0', required: false)
use_system_fmt = false
if fmt_dep.found()
  use_system_fmt = true
  message('Using system fmt')
else
  message('Using bundled fmt from ext/fmt')
  # Check if submodule is initialized
  fmt_include_dir = meson.current_source_dir() / 'ext' / 'fmt' / 'include'
  if not fs.is_dir(fmt_include_dir)
    warning('fmt submodule not initialized. Run: git submodule update --init ext/fmt')
  endif
  # Use header-only mode for bundled fmt
  fmt_dep = declare_dependency(
    include_directories: include_directories('ext/fmt/include'),
    compile_args: ['-DFMT_HEADER_ONLY']
  )
endif

# yaml-cpp dependency with fallback
yamlcpp_dep = dependency('yaml-cpp', version: '>=0.6', required: false)
use_system_yamlcpp = false
if yamlcpp_dep.found()
  use_system_yamlcpp = true
  message('Using system yaml-cpp')
else
  message('Using bundled yaml-cpp from ext/yaml-cpp')
  # Check if submodule is initialized
  yamlcpp_include_dir = meson.current_source_dir() / 'ext' / 'yaml-cpp' / 'include'
  if not fs.is_dir(yamlcpp_include_dir)
    warning('yaml-cpp submodule not initialized. Run: git submodule update --init ext/yaml-cpp')
  endif
  # Build yaml-cpp from submodule
  yamlcpp_proj = subproject('yaml-cpp')
  yamlcpp_dep = yamlcpp_proj.get_variable('yaml-cpp_dep')
endif

# Boost dependency (header-only, required)
boost_dep = dependency('boost', version: '>=1.83', required: true,
  modules: []  # Header-only
)

# HDF5 dependency (optional)
hdf5_dep = dependency('hdf5', required: false)
use_hdf5 = hdf5_dep.found()
if use_hdf5
  message('HDF5 support enabled')
else
  message('HDF5 support disabled')
endif

# HighFive dependency (header-only, optional)
highfive_dep = dependency('HighFive', required: false)
use_system_highfive = false
if highfive_dep.found()
  use_system_highfive = true
  message('Using system HighFive')
elif use_hdf5
  message('Using bundled HighFive from ext/HighFive')
  highfive_dep = declare_dependency(
    include_directories: include_directories('ext/HighFive/include')
  )
else
  highfive_dep = dependency('', required: false)
endif

# BLAS/LAPACK dependency (optional)
lapack_dep = dependency('lapack', required: false)
blas_dep = dependency('blas', required: false)
use_lapack = false

if lapack_dep.found() and blas_dep.found()
  use_lapack = true
  message('Using system BLAS/LAPACK')
elif is_darwin
  # Try Accelerate framework on macOS
  accelerate_dep = dependency('Accelerate', required: false)
  if accelerate_dep.found()
    lapack_dep = accelerate_dep
    blas_dep = dependency('', required: false)
    use_lapack = true
    message('Using macOS Accelerate framework')
  endif
endif

if not use_lapack
  message('BLAS/LAPACK not found, will use Eigen for linear algebra')
  lapack_dep = dependency('', required: false)
  blas_dep = dependency('', required: false)
endif

# Git commit information
git = find_program('git', required: false)
git_commit = 'unknown'
if git.found()
  git_result = run_command(git, 'rev-parse', '--short', 'HEAD', 
                           check: false, capture: true)
  if git_result.returncode() == 0
    git_commit = git_result.stdout().strip()
  endif
endif

# Data directory
datadir = get_option('datadir')
if datadir == ''
  datadir = get_option('prefix') / 'share' / 'cantera' / 'data'
endif

# Generate config.h
conf_data = configuration_data()
conf_data.set_quoted('CANTERA_VERSION', meson.project_version())
conf_data.set_quoted('CANTERA_SHORT_VERSION', cantera_short_version)
conf_data.set10('CT_USE_SYSTEM_EIGEN', use_system_eigen)
conf_data.set10('CT_USE_SYSTEM_EIGEN_PREFIXED', false)  # Assume unprefixed
conf_data.set10('CT_USE_SYSTEM_FMT', use_system_fmt)
conf_data.set10('CT_USE_SYSTEM_YAMLCPP', use_system_yamlcpp)
conf_data.set10('CT_USE_LAPACK', use_lapack)
conf_data.set10('CT_USE_HDF5', use_hdf5)
conf_data.set10('CT_USE_SYSTEM_HIGHFIVE', use_system_highfive)
conf_data.set10('CT_SUNDIALS_USE_LAPACK', false)  # Conservative default

# Fortran settings
conf_data.set('FTN_TRAILING_UNDERSCORE', '#define FTN_TRAILING_UNDERSCORE')
conf_data.set('LAPACK_FTN_STRING_LEN_AT_END', '#define LAPACK_FTN_STRING_LEN_AT_END')
conf_data.set('LAPACK_FTN_TRAILING_UNDERSCORE', '#define LAPACK_FTN_TRAILING_UNDERSCORE')

config_h = configure_file(
  input: 'include/cantera/base/config.h.in',
  output: 'config.h',
  configuration: conf_data,
  install: true,
  install_dir: get_option('includedir') / 'cantera' / 'base'
)

# Include directories
cantera_include_dirs = include_directories(
  'include',
  'src',
  '.'  # For generated config.h
)

# Build the library
subdir('src')

# Summary
summary({
  'Version': meson.project_version(),
  'Git commit': git_commit,
  'System Eigen': use_system_eigen,
  'System fmt': use_system_fmt,
  'System yaml-cpp': use_system_yamlcpp,
  'System SUNDIALS': use_system_sundials,
  'BLAS/LAPACK': use_lapack,
  'HDF5 support': use_hdf5,
  'System HighFive': use_system_highfive,
}, section: 'Configuration')
