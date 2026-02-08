# Meson Build System for Cantera

This directory contains the initial Meson build configuration for Cantera's C++ library.

## Status

The Meson build system is under development and currently builds most of the C++ library (~60% of targets). The build is blocked by the SUNDIALS dependency which needs to be added as a subproject or installed system-wide.

## Quick Start

### Prerequisites

- Meson >= 1.0.0
- Ninja build system
- C++20 compatible compiler (GCC, Clang, or MSVC)
- System dependencies:
  - Eigen3 >= 3.4 (or will use ext/eigen submodule)
  - Boost >= 1.83 (required, must be system-installed)
  - Optional: HDF5, BLAS/LAPACK

### Building

1. Initialize required submodules:
   ```bash
   git submodule update --init ext/fmt ext/yaml-cpp
   ```

2. Create a compatibility symlink for fmt (required for older fmt references):
   ```bash
   cd ext/fmt/include/fmt && ln -s ranges.h join.h
   ```

3. Configure the build:
   ```bash
   meson setup builddir
   ```

4. Compile:
   ```bash
   meson compile -C builddir
   ```

5. Install (optional):
   ```bash
   meson install -C builddir
   ```

## Dependencies

### Required
- **Boost** (>= 1.83): Header-only, must be system-installed
- **Eigen3** (>= 3.4): Uses system if available, otherwise ext/eigen submodule
- **fmt** (>= 9.1.0): Uses system if available, otherwise ext/fmt (header-only mode)
- **yaml-cpp** (>= 0.6): Uses system if available, otherwise auto-discovers and builds from ext/yaml-cpp

### Optional
- **SUNDIALS** (>= 6.0): Required for ODE/DAE integration (not yet implemented)
- **HDF5**: For HDF5 data file support
- **BLAS/LAPACK**: For optimized linear algebra (falls back to Eigen)
- **HighFive**: C++ wrapper for HDF5

## Configuration Options

- `cantera_datadir`: Directory for Cantera data files (default: `prefix/share/cantera/data`)

Example:
```bash
meson setup builddir -Dcantera_datadir=/opt/cantera/data
```

## Build Features

### Automatic Source Discovery

The Meson build automatically discovers source files using glob patterns, similar to how SCons uses `multi_glob()`:

```meson
# Automatically find all .cpp files in each module directory
base_files = run_command(find, 'base', '-name', '*.cpp', ...)
foreach f : base_files
  base_sources += files('base' / f)
endforeach
```

This eliminates the need to manually list individual source files.

### Dependency Management

Dependencies are tried in order:
1. **System packages** (via pkg-config/CMake) 
2. **Bundled sources** (from ext/ submodules with automatic file discovery)

For yaml-cpp, Meson automatically discovers all source files from `ext/yaml-cpp/src` instead of using the CMake build (which has issues with symlinked subprojects).

## Known Limitations

1. **SUNDIALS not yet supported**: The SUNDIALS dependency needs to be implemented as a Meson subproject or made available system-wide
2. **No Python bindings**: Only C++ library is currently supported
3. **No Fortran interface**: F90 interface not yet implemented
4. **Limited testing**: Build system needs more testing across platforms

## Comparison with SCons

The Meson build aims to eventually replace the SCons build system with these benefits:
- Faster build times with better parallelization
- Better IDE integration
- Standard pkg-config support
- Simpler dependency management
- Cross-compilation support
- **Automatic source discovery** like SCons `multi_glob()`

## Contributing

This is initial work to replace the SCons build system. Contributions are welcome to:
- Add SUNDIALS subproject support
- Improve dependency detection
- Add testing infrastructure
- Port more configuration options from SCons

## Files

- `meson.build`: Root build configuration
- `meson_options.txt`: Build options
- `src/meson.build`: C++ source compilation with automatic file discovery
- `include/cantera/base/config.h.meson.in`: Configuration header template
- `include/cantera/ext/`: Symlinks to submodule headers (fmt, yaml-cpp, Eigen)
