# Source files for the Cantera library
# Using run_command to automatically discover source files via find

# Helper to get cpp files from a directory
find_cpp = find_program('find')

# Gather all .cpp files from each module subdirectory
base_files = run_command(find_cpp, 'base', '-maxdepth', '1', '-name', '*.cpp', '-printf', '%f\n', check: true, capture: true).stdout().strip().split('\n')
base_sources = []
foreach f : base_files
  if f != ''
    base_sources += files('base' / f)
  endif
endforeach

thermo_files = run_command(find_cpp, 'thermo', '-maxdepth', '1', '-name', '*.cpp', '-printf', '%f\n', check: true, capture: true).stdout().strip().split('\n')
thermo_sources = []
foreach f : thermo_files
  if f != ''
    thermo_sources += files('thermo' / f)
  endif
endforeach

tpx_files = run_command(find_cpp, 'tpx', '-maxdepth', '1', '-name', '*.cpp', '-printf', '%f\n', check: true, capture: true).stdout().strip().split('\n')
tpx_sources = []
foreach f : tpx_files
  if f != ''
    tpx_sources += files('tpx' / f)
  endif
endforeach

equil_files_cpp = run_command(find_cpp, 'equil', '-maxdepth', '1', '-name', '*.cpp', '-printf', '%f\n', check: true, capture: true).stdout().strip().split('\n')
equil_files_c = run_command(find_cpp, 'equil', '-maxdepth', '1', '-name', '*.c', '-printf', '%f\n', check: false, capture: true).stdout().strip().split('\n')
equil_sources = []
foreach f : equil_files_cpp
  if f != ''
    equil_sources += files('equil' / f)
  endif
endforeach
foreach f : equil_files_c
  if f != ''
    equil_sources += files('equil' / f)
  endif
endforeach

numerics_files = run_command(find_cpp, 'numerics', '-maxdepth', '1', '-name', '*.cpp', '-printf', '%f\n', check: true, capture: true).stdout().strip().split('\n')
numerics_sources = []
foreach f : numerics_files
  if f != ''
    numerics_sources += files('numerics' / f)
  endif
endforeach

kinetics_files = run_command(find_cpp, 'kinetics', '-maxdepth', '1', '-name', '*.cpp', '-printf', '%f\n', check: true, capture: true).stdout().strip().split('\n')
kinetics_sources = []
foreach f : kinetics_files
  if f != ''
    kinetics_sources += files('kinetics' / f)
  endif
endforeach

transport_files = run_command(find_cpp, 'transport', '-maxdepth', '1', '-name', '*.cpp', '-printf', '%f\n', check: true, capture: true).stdout().strip().split('\n')
transport_sources = []
foreach f : transport_files
  if f != ''
    transport_sources += files('transport' / f)
  endif
endforeach

oned_files = run_command(find_cpp, 'oneD', '-maxdepth', '1', '-name', '*.cpp', '-printf', '%f\n', check: true, capture: true).stdout().strip().split('\n')
oned_sources = []
foreach f : oned_files
  if f != ''
    oned_sources += files('oneD' / f)
  endif
endforeach

zerod_files = run_command(find_cpp, 'zeroD', '-maxdepth', '1', '-name', '*.cpp', '-printf', '%f\n', check: true, capture: true).stdout().strip().split('\n')
zerod_sources = []
foreach f : zerod_files
  if f != ''
    zerod_sources += files('zeroD' / f)
  endif
endforeach

extensions_files = run_command(find_cpp, 'extensions', '-maxdepth', '1', '-name', '*.cpp', '-printf', '%f\n', check: true, capture: true).stdout().strip().split('\n')
extensions_sources = []
foreach f : extensions_files
  if f != ''
    extensions_sources += files('extensions' / f)
  endif
endforeach

# Combine all sources
cantera_sources = base_sources + thermo_sources + tpx_sources + equil_sources + numerics_sources + kinetics_sources + transport_sources + oned_sources + zerod_sources + extensions_sources

# Additional compile definitions for specific files
application_cpp_args = ['-DCANTERA_DATA="' + ct_datadir + '"']
global_cpp_args = ['-DGIT_COMMIT="' + git_commit + '"']

# Build application.cpp with special define
application_lib = static_library('cantera_application',
  'base/application.cpp',
  include_directories: [cantera_include_dirs, cantera_build_include],
  cpp_args: application_cpp_args,
  dependencies: [
    eigen_dep,
    fmt_dep,
    yamlcpp_dep,
    boost_dep,
    thread_dep,
    m_dep,
  ],
  pic: true,
)

# Build global.cpp with special define
global_lib = static_library('cantera_global',
  'base/global.cpp',
  include_directories: [cantera_include_dirs, cantera_build_include],
  cpp_args: global_cpp_args,
  dependencies: [
    eigen_dep,
    fmt_dep,
    yamlcpp_dep,
    boost_dep,
    thread_dep,
    m_dep,
  ],
  pic: true,
)

# Collect dependencies
cantera_deps = [
  eigen_dep,
  fmt_dep,
  yamlcpp_dep,
  boost_dep,
  thread_dep,
  m_dep,
]

if sundials_dep.found()
  cantera_deps += [sundials_dep, sundials_idas_dep, sundials_nvecserial_dep]
endif

if use_lapack
  cantera_deps += [lapack_dep]
  if blas_dep.found()
    cantera_deps += [blas_dep]
  endif
endif

if use_hdf5
  cantera_deps += [hdf5_dep, highfive_dep]
endif

# Build the shared library
cantera_lib = shared_library('cantera',
  cantera_sources,
  include_directories: [cantera_include_dirs, cantera_build_include],
  dependencies: cantera_deps,
  link_whole: [application_lib, global_lib],
  install: true,
  version: meson.project_version(),
  soversion: version_major,
)

# Install headers
install_subdir('include/cantera',
  install_dir: get_option('includedir'),
  strip_directory: false,
)

# Declare dependency for use by other subprojects
cantera_dep = declare_dependency(
  link_with: cantera_lib,
  include_directories: [cantera_include_dirs, cantera_build_include],
  dependencies: cantera_deps,
)

# pkg-config file
pkg = import('pkgconfig')
pkg.generate(cantera_lib,
  name: 'Cantera',
  description: 'Chemical kinetics, thermodynamics, and transport library',
  url: 'https://cantera.org',
  subdirs: 'cantera',
)
