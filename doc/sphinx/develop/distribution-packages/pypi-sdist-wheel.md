# Distributing sdist and Wheel Packages via PyPI

Cantera has been [distributed on PyPI](https://pypi.org/project/Cantera/) as a source
distribution (sdist) and a wheel since v2.6.0 in 2022. Wheels are Python's binary
distribution format, where source code is compiled to platform-specific libraries to
eliminate the need for user machines to have compilers and dependencies installed.
Source distributions are also uploaded so that platforms where a wheel is not available
can still attempt to build Cantera from the source code.

This page describes how to build sdists and wheels for Cantera, and how these are built
for distribution via PyPI. Some of the trade-offs in selecting the current build system
are documented in the relevant [enhancement
issue](https://github.com/Cantera/enhancements/issues/205). The source code to build and
publish the wheels is in a [specialized repository on
GitHub](https://github.com/Cantera/pypi-packages).

The overall process is conducted in two steps:

1. SCons is used to build an sdist from the git checkout with `scons sdist`
2. cibuildwheel/`build` is used to build wheels from the sdist, **not the git checkout**

There are a few reasons for this two-stage process. First, it ensures that our sdist can
be used to build a wheel on users' systems. Second, it tremendously simplifies the
cibuildwheel config, because installation of build dependencies is handled by
cibuildwheel at build time.

## TLDR

This set of steps should build an sdist and wheel on a developer's machine for local
testing.

1. Clone the repository (submodules are optional)
1. Install SCons and [build](https://pypi.org/project/build/)
1. Run `scons sdist`
   a. `scons clean sdist` can also be used to automatically remove the
   `build` folder prior to building the `sdist`.
1. Building the wheel requires that Boost and libhdf5 are available on the build system.
   If these are not installed in standard locations, you can set the `Boost_ROOT` and
   `HDF5_ROOT` environment variables to point to the correct location. On Linux, you
   will also need to install OpenBLAS (macOS uses Apple's Accelerate framework to
   provide LAPACK).
1. Change to the `build/python_sdist/dist` folder
1. Unpack the sdist: `tar xf cantera-3.1.0a2.tar.gz`. Note the version `3.1.0a2` may be
   different.
1. Change into the unpacked directory
1. Build the wheel: `python -m build --wheel .` Note the trailing `.` which indicates
   the current directory. You can also add the `-v` or `-vv` flags to increase build
   verbosity.

The wheel should be located in the `build/python_sdist/dist/cantera-3.1.0a2/dist`
folder. You can inspect the contents using any Zip-file reader. The CLI tool `unzip` can
list the contents with the `-l` flag. You can install the wheel by running `python -m
pip install <path to the wheel>`.

## Source Distributions

Sdists are built from the git checkout using SCons. Starting with v3.1, the only
dependencies required to build the sdist are SCons and
`build` (also known as `pyproject-build`). Dependencies such as NumPy, Boost, or any of
the dependencies available as submodules are not necessary to build the sdist.

To build the sdist, run `scons sdist` from the repository root.

```{note}
For versions 2.6 and 3.0, all the usual dependencies to build Cantera from source
are required. Starting with 3.1, the build process for the sdist was simplified to
pull in external dependencies at wheel-build time, rather than copying code from the
Cantera external source tree into the sdist.
```

The code for the sdist is located in the `interfaces/python_sdist` folder. Building the
sdist takes two steps:

1. SCons runs the `build_sdist.py` file to copy all the relevant source files into the
   `build/python_sdist` folder.
2. SCons runs `build` to turn the source files into a compliant sdist tarball using the
   `scikit-build-core` package as the build backend.

Both of these steps are done at "configure" time in SCons, when SCons is first parsing
the SConstruct and SConscript files. If the build is successful, SCons will exit before
performing any configuration or platform checks. This is why the build is done at
configure time, to avoid these time-consuming checks that aren't necessary since the
sdist just copies files into the `build` folder.

If successful, the `.tar.gz` file will be located in the `build/python_sdist/dist`
folder.

### Sdist configuration

The configuration for the sdist is done in the `pyproject.toml.in` file. This is a
template file that is filled in at build time by the `build_sdist.py` script. It needs
to be filled in by SCons because that's where we define the version of Cantera, in the
`SConstruct` file.

You can read more about the standard fields in the `pyproject.toml` file in the [PyPA
documentation](https://packaging.python.org/en/latest/specifications/pyproject-toml/).
You can read more about the `scikit-build` fields in the [`scikit-build-core`
documentation](https://scikit-build-core.readthedocs.io/en/latest/).

## Wheels

As mentioned above, wheels (`.whl` files) are Python's standard binary distribution
artifact. Similar to Conda packages, wheels include compiled libraries specific to a
particular platform, so that end users don't need to have compilers installed to use a
package.

Building Cantera's wheels is done using `build`, similar to the sdist. However, the
wheel is built from the unpacked sdist tarball, and you must pass the `--wheel` flag to
`build`. Then, `build` will install `scikit-build-core` to manage building the wheel
from Cantera's C++ and Cython source files.

`scikit-build-core` uses CMake to manage the compiling and linking process for the
wheel. Most of the configuration for the wheel build is therefore done in the
`CMakeLists.txt` files in the `interfaces/python_sdist` folder. These are copied into
the `build/python_sdist` folder and from there into the actual sdist tarball.

```{note}
The minimum required version of CMake is 3.27 to support all the options we use with
[`FetchContent`](https://cmake.org/cmake/help/latest/module/FetchContent.html). A
compatible version should be installed by `scikit-build-core` if it isn't available
on the system. Linux and macOS hosts also require the Ninja build system, version 1.11
or greater. This should also be installed by `scikit-build-core` at build time if isn't
available on the system.
```

There are essentially 4 stages to the wheel build:

1. External dependencies (`eigen`, `fmt`, `SUNDIALS`, `HighFive`, `yaml-cpp`) are either
   located on the system, or downloaded and compiled using the `FetchContent` mechanism
1. Cantera's C++ source code is built and linked to the appropriate dependencies
1. Cython runs and generates C++ source code for the Python interface
1. The generated C++ is compiled, linked to `libcantera`, and packaged into the `.whl`
   file

```{note}
Cantera's other dependencies, specifically Boost and libhdf5, are neither installed nor
built during the wheel build and must be pre-installed on the build system. I couldn't
find a way to have eigen and HighFive, respectively, find those dependencies if they
were also installed using `FetchContent`. If Boost and libhdf5 are installed in
non-standard directories, you can use the `Boost_ROOT` or `Boost_INCLUDE_DIRS` and
`HDF5_ROOT` variables to specify the correct locations. See the CMake documentation for
[`find_package()`](https://cmake.org/cmake/help/latest/command/find_package.html).
```

Steps 1 and 2 above are set up in the `interfaces/python_sdist/src/CMakeLists.txt` file.
Steps 3 and 4 are set up in the `interfaces/python_sdist/cantera/CMakeLists.txt` file.

## Publishing the sdist and Wheels

When a version of Cantera is ready to be published to PyPI, the workflow in the
[pypi-packages](https://github.com/Cantera/pypi-packages) repository must be executed.
The workflow essentially conducts the steps in the [](#tldr) section above.

The one big addition to the GitHub Actions workflow is that there are a couple of shell
scripts in the `pypi-packages` repository to build libhdf5, HighFive, and SUNDIALS for
macOS and Windows.

For Linux builds, we have a [custom-built manylinux
image](https://github.com/Cantera/cantera-base-manylinux) that already has the
dependencies installed. This saves time compiling our dependencies on Linux/ARM builds
because those are emulated and already very slow.

As much of the configuration for `cibuildwheel` as possible is done statically in
`pyproject.toml.in`. Only values that are dynamically determined when the workflow runs
should be set in the workflow.

This workflow should be automatically triggered when a new version of Cantera is tagged,
identified by a tag matching the pattern `v*`. This is managed by the GitHub action
workflow in `.github/workflows/packaging.yml` in the main Cantera repository. You can
[check the status](https://github.com/Cantera/pypi-packages/actions/workflows/python-package.yml)
of these package builds or trigger the workflow manually (in case it fails and requires
updates after the tag is pushed).

After the workflow successfully builds the packages, it will prompt the Cantera
maintainers to approve deploying the release to PyPI. Before doing so, download and test
some of the wheels that are uploaded as "artifacts" on GitHub.

## Tips

* You can manually run the packaging build in the `pypi-packages` repository by going
  to the _Actions_ tab in the repository and manually running the workflow. Choose
  your branch from the dropdown and paste in a commit from the `Cantera/cantera`
  repository. You can paste any valid commit associated with Cantera, including commits
  from branches associated with pull requests in `Cantera/cantera`.
* With CMake 3.27 or later, `*_ROOT` locations for dependencies can be spelled as either
  the canonical spelling defined by the library (for example, `Boost_ROOT`) or as all
  upper-case (for example, `BOOST_ROOT`). This behavior is configured by [CMake Policy
  0144](https://cmake.org/cmake/help/latest/policy/CMP0144.html).

[manylinux]:
