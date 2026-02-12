# Meson Build System for Cantera

This directory contains the Meson build configuration for Cantera's C++ library.

## Status

The Meson build system builds the complete C++ library using system-installed dependencies only.

## Quick Start

### Prerequisites

- Meson >= 1.0.0
- Ninja build system
- C++20 compatible compiler (GCC, Clang, or MSVC)
- System dependencies (all required):
  - Eigen3 >= 3.4
  - Boost >= 1.83
  - fmt >= 9.1.0
  - yaml-cpp >= 0.6
  - SUNDIALS >= 6.0 (cvodes, idas, nvecserial)
  - Optional: HDF5, BLAS/LAPACK, HighFive

### Building

1. Configure the build:
   ```bash
   meson setup builddir
   ```

2. Compile:
   ```bash
   meson compile -C builddir
   ```

3. Install (optional):
   ```bash
   meson install -C builddir
   ```

## Dependencies

All dependencies must be installed as system packages. The Meson build does not support fallback to bundled libraries in ext/ submodules.

### Required System Packages

- **Boost** (>= 1.83): Header-only library
- **Eigen3** (>= 3.4): Header-only library for linear algebra
- **fmt** (>= 9.1.0): Formatting library
- **yaml-cpp** (>= 0.6): YAML parser/emitter
- **SUNDIALS** (>= 6.0): ODE/DAE solvers (cvodes, idas, nvecserial components required)

### Optional System Packages

- **HDF5**: For HDF5 data file support
- **BLAS/LAPACK**: For optimized linear algebra (falls back to Eigen if not available)
- **HighFive**: C++ wrapper for HDF5 (required if HDF5 support is desired)

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

### Installation

On Ubuntu/Debian:
```bash
sudo apt install meson ninja-build libboost-dev libeigen3-dev \
  libfmt-dev libyaml-cpp-dev libsundials-dev
```

On macOS with Homebrew:
```bash
brew install meson ninja boost eigen fmt yaml-cpp sundials
```

On Fedora/RHEL:
```bash
sudo dnf install meson ninja-build boost-devel eigen3-devel \
  fmt-devel yaml-cpp-devel sundials-devel
```

## Known Limitations

1. **No Python bindings**: Only C++ library is currently supported
2. **No Fortran interface**: F90 interface not yet implemented
3. **Limited testing**: Build system needs more testing across platforms
4. **System packages required**: No support for bundled ext/ submodules

## Comparison with SCons

The Meson build aims to eventually replace the SCons build system with these benefits:
- Faster build times with better parallelization
- Better IDE integration
- Standard pkg-config support
- Simpler dependency management (system packages only)
- Cross-compilation support
- **Automatic source discovery** like SCons `multi_glob()`

## Contributing

Contributions are welcome to:
- Add Python bindings support
- Add testing infrastructure
- Port more configuration options from SCons
- Test on more platforms

## Files

- `meson.build`: Root build configuration
- `meson_options.txt`: Build options
- `src/meson.build`: C++ source compilation with automatic file discovery
- `include/cantera/base/config.h.meson.in`: Configuration header template
