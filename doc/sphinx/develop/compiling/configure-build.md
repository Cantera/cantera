# Configure and Build Cantera

(sec-determine-config)=
## Determine configuration options

- Run `scons help --options` to see a list of all of the configuration options for
  Cantera, or see all of the options on the [Configuration Options](scons-config)
  page.

- Configuration options are specified as additional arguments to the `scons`
  command. For example:

  ```bash
  scons command option_name=value
  ```

  where `scons` is the program that manages the build steps, and `command`
  is most commonly one of

  - `build`
  - `test`
  - `clean`

  Other commands are explained in the [Build Commands](sec-build-commands) section.

- SCons saves configuration options specified on the command line in the file
  `cantera.conf` in the root directory of the source tree, so generally it is
  not necessary to respecify configuration options when rebuilding Cantera. To
  unset a previously set configuration option, either remove the corresponding
  line from `cantera.conf` or use the syntax:

  ```bash
  scons command option_name=
  ```

- Sometimes, changes in your environment can cause SCons's configuration tests
  (for example, checking for libraries or compiler capabilities) to unexpectedly fail.
  To force SCons to re-run these tests rather than trusting the cached results,
  run scons with the option `--config=force`.

- The following lists of options are not complete, they show only some commonly
  used options. The entire list of options can be found on the
  [Configuration options](scons-config) page.

### Common Options

- [`debug`](sconsopt-debug)
- [`optimize`](sconsopt-optimize)
- [`prefix`](sconsopt-prefix)

### Specifying Paths for Cantera's Dependencies

- [`blas_lapack_libs`](sconsopt-blas-lapack-libs)

  - On macOS, the Accelerate framework is automatically used to provide optimized
    versions of BLAS and LAPACK, so the `blas_lapack_libs` option should generally be
    left unspecified.

- [`blas_lapack_dir`](sconsopt-blas-lapack-dir)

- [`boost_inc_dir`](sconsopt-boost-inc-dir)

- [`sundials_include`](sconsopt-sundials-include)

- [`sundials_libdir`](sconsopt-sundials-libdir)

- [`hdf_include`](sconsopt-hdf-include)

- [`hdf_libdir`](sconsopt-hdf-libdir)

- [`extra_inc_dirs`](sconsopt-extra-inc-dirs)

- [`extra_lib_dirs`](sconsopt-extra-lib-dirs)

### Python Module Options

Compiling the Cantera Python module requires that NumPy and Cython are installed
for the target installation of Python. The following SCons options control how
the Python module is built:

- [`python_package`](sconsopt-python-package)

- [`python_cmd`](sconsopt-python-cmd)

  - By default, SCons will try to build the full Python interface for copy of
    Python that is running SCons. Use this option if you wish to build Cantera
    for a different Python installation.

- [`python_prefix`](sconsopt-python-prefix)

### Windows Only Options

:::{note}
The `cantera.conf` file uses the backslash character `\` as an escape
character. When modifying this file, backslashes in paths need to be escaped
like this: `boost_inc_dir = 'C:\\Program Files (x86)\\boost\\include'`
This does not apply to paths specified on the command line. Alternatively,
you can use forward slashes (`/`) in paths.
:::

- In Windows there aren't any proper default locations for many of the packages
  that Cantera depends on, so you will need to specify these paths explicitly.

- Remember to put double quotes around any paths with spaces in them, such as
  `"C:\Program Files"`.

- By default, SCons attempts to use the same architecture as the copy of Python
  that is running SCons, and the most recent installed version of the Visual
  Studio compiler. If you aren't building the Python module, you can override
  this with the configuration options `target_arch` and `msvc_toolset_version`.

- To compile with MinGW, specify the [`toolchain`](sconsopt-toolchain) option:

  ```
  toolchain=mingw
  ```

- [`msvc_toolset_version`](sconsopt-msvc-toolset-version)

- [`msvc_version`](sconsopt-msvc-version)

- [`target_arch`](sconsopt-target-arch)

- [`toolchain`](sconsopt-toolchain)

### Fortran Module Options

Building the Fortran module requires a compatible Fortran compiler. SCons will
attempt to find a compatible compiler by default in the `PATH` environment
variable. The following options control how the Fortran module is built:

- [`f90_interface`](sconsopt-f90-interface)
- [`FORTRAN`](sconsopt-fortran)

### Documentation Options

The following options control if the documentation is built:

- [`doxygen_docs`](sconsopt-doxygen-docs)
- [`sphinx_docs`](sconsopt-sphinx-docs)

### Less Common Options

- [`CC`](sconsopt-cc)
- [`CXX`](sconsopt-cxx)
- [`env_vars`](sconsopt-env-vars)
- [`layout`](sconsopt-layout)
- [`logging`](sconsopt-logging)
- [`gtest_flags`](sconsopt-gtest-flags)

(sec-build-commands)=
## Build Commands

The following *commands* are possible as arguments to SCons:

```bash
scons command
```

- `scons help`
  : Print a list of available SCons commands.
- `scons help --options`
  : Print a description of user-specifiable options.
- `scons build`
  : Compile Cantera and the language interfaces using
    default options.
- `scons clean`
  : Delete files created while building Cantera.
- `scons install`
  : Install Cantera.
- `scons uninstall`
  : Uninstall Cantera.
- `scons test`
  : Run all tests which did not previously pass or for which the
    results may have changed.
- `scons test-reset`
  : Reset the passing status of all tests.
- `scons test-clean`
  : Delete files created while running the tests.
- `scons test-help`
  : List available tests.
- `scons test-NAME`
  : Run the test named `NAME`.
- `scons <command> dump`
  : Dump the state of the SCons environment to the
    screen instead of doing `<command>`, for example,
    `scons build dump`. For debugging purposes.
- `scons samples`
  : Compile the C++ and Fortran samples.
- `scons msi`
  : Build a Windows installer (.msi) for Cantera.
- `scons sphinx`
  : Build the Sphinx documentation
- `scons doxygen`
  : Build the Doxygen documentation

## Compile Cantera & Test

- Run SCons with the list of desired configuration options:

  ```bash
  scons build ...
  ```

- If Cantera compiles successfully, you should see a message that looks like:

  ```none
  ********************** Compilation completed successfully **********************

  - To run the test suite, type 'scons test'.
  - To list available tests, type 'scons test-help'.
  - To install, type 'scons install'.

  ********************************************************************************
  ```

- If you do not see this message, check the output for errors to see what went
  wrong. You may also need to examine the contents of `config.log`.

- Cantera has a series of tests that can be run with the command:

```bash
scons test
```

- When the tests finish, you should see a summary indicating the number of
  tests that passed and failed.

- If you have tests that fail, try looking at the following to determine the
  source of the error:

  - Messages printed to the console while running `scons test`
  - Output files generated by the tests

### Building Documentation

To build the Cantera HTML documentation, run the commands:

```bash
scons doxygen
scons sphinx
```

or append the options `sphinx_docs=y` and `doxygen_docs=y` to the build
command:

```bash
scons build doxygen_docs=y sphinx_docs=y
```

## Installing Cantera

- To install Cantera into default directories, run the SCons installer as:

  ```bash
  scons install
  ```

  which may require super-user permissions if the installation directory is protected.
  An installation into an active [Conda environment](sec-conda) is recommended.

  Alternatively, you can specify location and layout of the installation via

  - [`prefix`](sconsopt-prefix)
  - [`layout`](sconsopt-layout)
  - [`python_prefix`](sconsopt-python-prefix)

- If Cantera installs successfully, you should see a message that looks similar to

  ```none
  **************** Cantera 3.x.y has been successfully installed *****************

  File locations:

    library files               C:/path/to/prefix/Library/lib
    C++ headers                 C:/path/to/prefix/Library/include
    samples                     C:/path/to/prefix/share/cantera/samples
    data files                  C:/path/to/prefix/share/cantera/data
    input file converters       C:/path/to/prefix/Scripts
    Python package              C:/path/to/prefix/Lib/site-packages
    Python examples             C:/path/to/prefix/share/cantera/samples/python

  ********************************************************************************
  ```

  where slight variations may depend on operating system and configuration.
