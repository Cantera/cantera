(sec-compilation-reqs)=
# Compilation Requirements

Click the buttons below to see the required software that you must install to
compile Cantera on your operating system.

`````{grid} auto

````{grid-item}
```{button-ref} sec-conda
:color: primary
:shadow:
```
````

````{grid-item}
```{button-ref} sec-ubuntu-debian-reqs
:color: primary
:shadow:
```
````

````{grid-item}
```{button-ref} sec-opensuse-reqs
:color: primary
:shadow:
  OpenSUSE & SLE
```
````

````{grid-item}
```{button-ref} sec-fedora-reqs
:color: primary
:shadow:
  Fedora & RHEL
```
````


````{grid-item}
```{button-ref} sec-macos
:color: primary
:shadow:
```
````

````{grid-item}
```{button-ref} sec-windows
:color: primary
:shadow:
```
````

`````

(sec-conda)=
## Conda & Anaconda

### General Notes

- These instructions will set you up to build Cantera with the dependencies installed in
  a Conda environment
- You will need to install compilers for your system by following the instructions in
  the sections below to install the compiler for your operating system.
- By default, Cantera is installed into the active conda environment, where the layout
  of the directory structure corresponds to the [configuration option](configure-build)
  `layout=conda`.

(sec-conda-reqs)=
### Conda Requirements

- Install [Anaconda](https://www.anaconda.com/download/),
  [Miniconda](https://conda.io/miniconda.html), or
  [Miniforge](https://github.com/conda-forge/miniforge).

- Launch the command line interface:

  - On macOS and Linux, the installer should add the appropriate activation mechanism
    for your normal terminal by default. You can test this by running

    ```bash
    conda --version
    ```

    in the terminal. If there is no output or an error appears, locate your Conda
    installation and run the following code in the terminal:

    ```bash
    /path/to/conda/install/folder/bin/conda init --all
    ```

    Then restart your terminal or shell.

  - On Windows, use the Anaconda PowerShell to run the build process (available from
    the Start Menu). When using MSVC compilers, you also need to set environment
    variables for x64-native tools (see [Developer command file locations](https://docs.microsoft.com/en-us/cpp/build/building-on-the-command-line?view=msvc-170#developer_command_file_locations))
    by running

    ```bash
    . "C:\path\to\MSVC\Auxiliary\Build\vcvars64.bat"
    ```

    (note that the period `'.'` is part of the command). The path can be found as
    follows: locate the **x64 Native Tools Command Prompt** in the Start Menu,
    right-click, select **More > Open File Location**, right-click on the shortcut,
    select **Properties** and copy the **Target** command.

- Create an environment `ct-build` with the dependencies to build Cantera. Create a
  file called `environment.yaml` with the following content

  ```yaml
  name: ct-build
  channels:
  - conda-forge
  dependencies:
  - python  # Cantera supports Python 3.10 and up
  - scons  # build system
  - boost-cpp  # C++ dependency
  - hdf5  # optional C++ dependency
  # - highfive  # C++ dependency; uncomment to override Cantera default
  # - sundials  # uncomment to override Cantera default
  # - fmt  # uncomment to override Cantera default
  # - eigen  # uncomment to override Cantera default
  # - yaml-cpp  # uncomment to override Cantera default
  # - libgomp  # optional (OpenMP implementation when using GCC)
  - cython  # needed to build Python package
  - numpy  # needed to build Python package
  - pip  # needed to build Python package
  - wheel  # needed to build Python package
  - setuptools  # needed to build Python package
  - pytest  # needed for the Python test suite
  # - pytest-cov  # optional (needed if running with test coverage enabled)
  - ruamel.yaml  # needed for converter scripts
  # - pandas  # optional (needed for pandas interface and some examples)
  # - scipy  # optional (needed for some examples)
  # - matplotlib  # optional (needed for plots and some examples)
  # - python-graphviz  # optional (needed for reaction path diagrams and some examples)
  - ipython  # optional (needed for nicer interactive command line)
  # - jupyter  # optional (needed for Jupyter Notebook)
  # - sphinx  # optional (needed for documentation)
  # - pydata-sphinx-theme  # optional (needed for documentation)
  # - myst-nb  # optional (needed for documentation)
  # - myst-parser  # optional (needed for documentation)
  # - sphinx-gallery  # optional (needed for documentation)
  # - sphinx-argparse  # optional (needed for documentation)
  # - sphinx-copybutton  # optional (needed for documentation)
  # - sphinx-design  # optional (needed for documentation)
  # - sphinxcontrib-bibtex  # optional (needed for documentation)
  # - sphinxcontrib-matlabdomain  # optional (needed for documentation)
  # - sphinxcontrib-doxylink  # optional (needed for documentation)
  # - doxygen  # optional (needed for documentation)
  # - graphviz  # optional (needed for documentation)
  # - texlive-core  # optional (needed for documentation)
  # - perl  # optional (needed for documentation)
  # - coolprop  # optional (needed for some examples)
  # - pint  # optional (needed for some examples)
  # - jinja2  # optional (needed for code generation; example: .NET interface)
  # - pip:  # optional (list of PyPI managed packages)
  #   - "git+https://github.com/Cantera/sphinx-tags.git@main"  # optional (needed for documentation)
  ```

  The environment is then created and activated using

  ```bash
  conda env create -f environment.yaml
  conda activate ct-build
  ```

  After creating the environment, it can be updated from within `ct-build` using

  ```bash
  conda env update -f environment.yaml --prune
  ```

- (Optional) If you want to override external libraries packaged with Cantera
  (`sundials`, `fmt`, `eigen`, `yaml-cpp`), simply uncomment corresponding lines in the
  file `environment.yaml` above. Note that specific versions can be forced by providing
  version numbers (example: replace `sundials` by `sundials=5.8` to install version
  `5.8`).

- (Optional) If you want to build the documentation, make sure to uncomment lines for
  `pip`, `sphinx`, `doxygen`, and any other lines marked as "needed for documentation"
  in the `environment.yaml` sample above.

- (Cantera \< 2.6 only) On previous Cantera versions, the build process required
  configuration options `boost_inc_dir` and `prefix` (see [configuration
  options](configure-build)); starting with Cantera 2.6, these settings are detected
  automatically.

:::{note}
As the compiled code is based on the conda environment `ct-build`, it is only
usable from within that environment. This means that in order to use the compiled
Cantera package, you have to activate your `ct-build` environment first.
:::

```{button-ref} source-code
:color: primary
:shadow:
:align: right
Next: Download the Source Code
```

(sec-linux)=
## Linux

### General Notes

- To download the source code, installing `git` is highly recommended in addition
  to the requirements listed below.
- The following instructions use the system-installed versions of Python, but alternate
  installations such as the Anaconda distribution of Python can be used as well.
- Cython is only required to be installed for the version of Python that also has SCons
  installed; following the instructions below will install Cython for the version of
  Python installed in the system directories. The minimum compatible Cython version is
  0.29.31. If your distribution does not contain a suitable version, you may be able to
  install a more recent version using Pip.
- Users of other distributions should install the equivalent packages, which may have
  slightly different names.
- In addition to the operating systems below, Cantera should work on any Unix-like
  system where the necessary prerequisites are available, but some additional
  configuration may be required.

(sec-ubuntu-debian-reqs)=

### Ubuntu & Debian

- Ubuntu 20.04 LTS (Focal Fossa) or newer

- Debian 11.0 (Bullseye) or newer

- The following packages must be installed to build any of the Cantera modules using
  your choice of package manager:

  ```
  g++ python3 scons libboost-dev libhdf5-dev
  ```

  - The HDF5 headers and libraries are not installed to directories on the compiler's
    default search path. When building Cantera, these paths need to be specified as
    options to `scons`, for example `extra_inc_dirs=/usr/include/hdf5/serial` and
    `extra_lib_dirs=/usr/lib/x86_64-linux-gnu/hdf5/serial`.

- If you want to use system system packages to provide the following dependencies,
  instead of the versions bundled with Cantera, you should also install:

  ```
  libsundials-dev libeigen3-dev libyaml-cpp-dev libfmt-dev
  ```

- In addition to the general packages, building the Python 3 module also requires:

  ```
  cython3 python3-setuptools python3-wheel python3-numpy python3-ruamel.yaml python3-pytest
  ```

  - Debian 12.0 (Bookworm) and Ubuntu 23.04 (Lunar Lobster) provide compatible Cython
    versions. For older releases, install Cython using Pip.

- In addition to the general packages, building the Fortran module also requires:

  ```
  gfortran
  ```

- In addition to the general packages, building the MATLAB toolbox also requires:

  - MATLAB version later than 2009a

    - Typically installed to:

      ```
      /opt/MATLAB/R20YYn
      ```

      where `YY` is a two digit year and `n` is either `a` or `b`

```{button-ref} source-code
:color: primary
:shadow:
:align: right
Next: Download the Source Code
```

(sec-fedora-reqs)=
### Fedora & RHEL

- The following packages must be installed to build any of the Cantera modules using
  your choice of package manager:

  ```
  gcc-c++ python3 scons boost-devel hdf5-devel
  ```

- If you want to use system system packages to provide the following dependencies,
  instead of the versions bundled with Cantera, you should also install:

  > sundials-devel eigen3-devel yaml-cpp-devel fmt-devel highfive-devel

- In addition to the general packages, building the Python 3 module also requires:

  ```
  python3-devel Cython python3-numpy python3-ruamel-yaml python3-pytest
  ```

- In addition to the general packages, building the Fortran module also requires:

  ```
  gcc-gfortran
  ```

- In addition to the general packages, building the MATLAB toolbox also requires:

  - MATLAB version later than 2009a

    - Typically installed to:

      ```
      /opt/MATLAB/R20YYn
      ```

      where `YY` is a two digit year and `n` is either `a` or `b`

```{button-ref} source-code
:color: primary
:shadow:
:align: right
Next: Download the Source Code
```

(sec-opensuse-reqs)=
### OpenSUSE & SUSE Linux Enterprise

- OpenSUSE Leap 15.5 or newer recommended

- The following packages must be installed to build any of the Cantera modules using
  your choice of package manager:

  ```
  gcc11-c++ python311 libboost_headers1_75_0-devel hdf5-devel python311-pip
  ```

  You can specify other version numbers for GCC, Python, and Boost, as long as they meet
  Cantera's minimum requirements.

- You will also need to install `scons` using the Pip version installed above.

- In addition to the general packages, building the Python module also requires:

  ```
  python3-devel
  ```

  as well as the following packages installed using Pip:

  > numpy wheel cython ruamel.yaml pytest

- In addition to the general packages, building the Fortran module also requires:

  ```
  gcc11-fortran
  ```

- In addition to the general packages, building the MATLAB toolbox also requires:

  - MATLAB version later than 2009a

    - Typically installed to:

      ```
      /opt/MATLAB/R20YYn
      ```

      where `YY` is a two digit year and `n` is either `a` or `b`

```{button-ref} source-code
:color: primary
:shadow:
:align: right
Next: Download the Source Code
```

(sec-windows)=
## Windows

### General Notes

- The build process will produce a Python module compatible with the version of Python
  used for the compilation. To generate different modules for other versions of Python,
  you will need to install those versions of Python and recompile.

- The following instructions use the versions of Python downloaded from
  <https://www.python.org/downloads/>, but alternate installations such as the Anaconda
  distribution of Python can be used as well.

- If you want to build the Matlab toolbox and you have a 64-bit copy of Windows, by
  default you will be using a 64-bit copy of Matlab, and therefore you need to compile
  Cantera in 64-bit mode. For simplicity, it is highly recommended that you use a 64-bit
  version of Python to handle this automatically.

- It is generally helpful to have SCons and Python in your `PATH` environment variable.
  This can be done by checking the appropriate box during the installation of Python or
  can be accomplished by adding the top-level Python directory and the `Scripts`
  subdirectory (for example, `C:\Python311;C:\Python311\Scripts`) to your `PATH`. The
  dialog to change the `PATH` is accessible from:

  ```
  Control Panel > System and Security > System > Advanced System Settings > Environment Variables
  ```

  Make sure that the installation of Python that has SCons comes first on your `PATH`.

- In order to use SCons to install Cantera to a system folder (for example,
  `C:\Program Files\Cantera`) you must run the `scons install` command in a
  command prompt that has been launched by selecting the *Run as Administrator*
  option.

(sec-windows-reqs)=
### Windows Requirements

- Windows 7 or later; either 32-bit or 64-bit

- To build any of the Cantera modules, you will need to install

  - Python

    - <https://www.python.org/downloads/>
    - Cantera supports Python 3.10 and higher
    - Be sure to choose the appropriate architecture for your system - either
      32-bit or 64-bit
    - When installing, make sure to choose the option to add to your `PATH`

  - SCons

    - <https://pypi.org/project/SCons/>
    - Be sure to choose the appropriate architecture for your system - either
      32-bit or 64-bit

  - One of the following supported compilers

    - Microsoft compilers

      - <https://visualstudio.microsoft.com/downloads/>
      - Known to work with Visual Studio 2017 (MSVC 14.1), Visual Studio 2019
        (MSVC 14.2), and Visual Studio 2022 (MSVC 14.3).

    - MinGW compilers

      - <http://mingw-w64.org/>
      - <http://tdm-gcc.tdragon.net/>
      - Known to work with Mingw-w64 12.2.

  - The Boost headers

    - <https://www.boost.org/doc/libs/1_82_0/more/getting_started/windows.html#get-boost>
    - It is not necessary to compile the Boost libraries since Cantera only uses
      the headers from Boost

- In addition to the general software, building the Python module also requires
  several Python packages: Cython, NumPy, setuptools, wheel, ruamel.yaml, and pytest.
  All of these can be installed using `pip`:

  ```bash
  py -m pip install setuptools wheel cython numpy ruamel.yaml pytest
  ```

```{button-ref} source-code
:color: primary
:shadow:
:align: right
Next: Download the Source Code
```

(sec-macos)=
## macOS

### General Notes

- Cantera 2.5.0 and higher do not support Python 2, which may be installed by default
  on your computer. You must install Python 3 from another source to be able to build
  Cantera. The instructions below use Homebrew.

(sec-mac-os-reqs)=
### macOS Requirements

- macOS 10.15 (Catalina) is required; Homebrew requires 11.0 or newer.

- To build any of the Cantera interfaces, you will need to install

  - Xcode

    - Download and install from the App Store

    - From a Terminal, run:

      ```bash
      sudo xcode-select --install
      ```

      and agree to the Xcode license agreement

  - Homebrew

    - <https://brew.sh>

    - From a Terminal, run:

      ```bash
      /bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/HEAD/install.sh)"
      ```

  - Once Homebrew is installed, the rest of the dependencies can be installed with:

    ```bash
    brew install python scons boost git hdf5 libomp
    ```

    Note that brew installs Python 3 by default, but does not over-write the existing
    system Python. When you want to use the brew-installed Python, check to make sure
    that `python3` and `pip3` refer to the Homebrew installation by running:

    ```bash
    which python3
    which pip3
    ```

    If these commands do not include the Homebrew path, you can run the correct ones as
    `$(brew --prefix)/bin/python3` and `$(brew --prefix)/bin/pip3`.

- In addition to the general software, building the Python module also requires:

  ```bash
  $(brew --prefix)/bin/pip3 install cython numpy wheel setuptools ruamel.yaml pytest
  ```

- In addition to the general software, building the Fortran module also requires:

  ```bash
  brew install gcc
  ```

- In addition to the general software, building the MATLAB toolbox also requires:

  - MATLAB version later than 2009a

    - Typically installed to:

      ```
      /Applications/MATLAB_R20YYn.app
      ```

      where `YY` is a two digit year and `n` is either `a` or `b`

- The Homebrew header and library directories will not be on the path for the system
  compiler (Xcode), so when compiling Cantera, you will need to provide the command line
  options `extra_inc_dirs=$(brew --prefix)/include` and
  `extra_lib_dirs=$(brew --prefix)/lib`.

```{button-ref} source-code
:color: primary
:shadow:
:align: right
Next: Download the Source Code
```
