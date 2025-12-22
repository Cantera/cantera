# Installing with Conda

(sec-install-conda)=

[Conda](https://docs.conda.io/projects/conda/en/stable/) is a package manager that can
be used to install Python packages and other software. Cantera packages are available
via the `conda-forge` channel, for use with Conda distributions that use that channel.
We highly recommend using the [Miniforge](https://conda-forge.org/download/)
distribution, which configures `conda-forge` as the default channel (in which case you
can omit the argument `--channel conda-forge` from any of the commands below).

For more details on how to use Conda, see the
[Conda documentation](https://docs.conda.io/projects/conda/en/latest/user-guide/index.html).

For instructions on upgrading an existing Conda-based installation of Cantera, see
[Upgrading from an earlier Cantera version](sec-conda-python-upgrade).

(sec-conda-python-interface)=

## Python interface

Cantera's Python interface can be installed from the popular `conda-forge` channel.
Packages are available for the following platforms:

- Windows (64-bit Intel)
- Linux (64-bit Intel, 64-bit ARM, and 64-bit PPCLE)
- macOS (64-bit ARM (Apple Silicon))

### Option 1: Create a new environment for Cantera

The following instructions will create a Conda environment where you can use Cantera
from Python. For this example, the environment is named `ct-env`. From the command
line (or the Conda Prompt on Windows), run:

```shell
conda create --name ct-env --channel conda-forge cantera ipython matplotlib jupyter
```

This will create an environment named `ct-env` with Cantera, IPython, Matplotlib, and
all their dependencies installed.

To use the scripts and modules installed in the `ct-env` environment, including Jupyter,
you must activate it it by running:

```shell
conda activate ct-env
```

### Option 2: Create a new environment using an environment file

This option is similar to **Option 1** but includes a few other packages that you may
find helpful as you're working with Cantera. Copy and paste the contents of the file
shown below into a file called `environment.yaml`. Then, save the the file somewhere and
remember that location.

```yaml
name: ct-env
channels:
- conda-forge
dependencies:
- python  # Cantera supports Python 3.10 and up
- cantera
- ipython  # optional (needed for nicer interactive command line)
- jupyter  # optional (needed for Jupyter Notebook)
- matplotlib  # optional (needed for plots)
- python-graphviz  # optional (needed for reaction path diagrams)
- pandas  # optional (needed for pandas interface)
```

From the command line (or the Conda Prompt on Windows), change directory into the
folder where you saved `environment.yaml`:

```shell
cd folder/where/you/saved
```

and then run:

```shell
conda env create -f environment.yaml
```

This will create an environment called `ct-env`. Once you've done that, you need to
activate the environment before using any scripts or modules that you just installed:

```shell
conda activate ct-env
```

### Option 3: Install the development version of Cantera

To install a recent development snapshot (that is, an alpha or beta version) of Cantera,
use the `conda-forge/label/cantera_dev` channel. Assuming you have an environment named
`ct-dev`, you can type:

```shell
conda activate ct-dev
conda install --channel conda-forge/label/cantera_dev cantera
```

If you later want to revert back to the stable version in that environment, first remove
and then reinstall Cantera:

```shell
conda activate ct-dev
conda remove cantera
conda install --channel conda-forge cantera
```

Alternatively, you can remove the `ct-dev` environment and follow Options 1 or 2 above
to create a new environment.

(sec-conda-python-upgrade)=

## Upgrading from an earlier Cantera version

If you already have Cantera installed in a Conda environment (named, for example,
`ct-dev`), you can upgrade it to the latest version available by running the commands:

```shell
conda activate ct-dev
conda update --channel conda-forge cantera
```

:::{attention}
This upgrade option will only work if the previous version of Cantera was also installed
from the `conda-forge` channel. If you used the packages from the `cantera` channel that
were also provided for Cantera 3.0 and earlier, this upgrade path will not work and you
should install Cantera in a new Conda environment.
:::

(sec-conda-development-interface)=

## Development (C++, Fortran 90, & MATLAB) Interface

The Cantera development interface provides header files and libraries needed to compile
your own C++, C, Fortran, and MATLAB applications that link to Cantera. It also
provides several sample programs and build scripts that you can adapt for
your own applications.

From the command line (or the Conda Prompt on Windows), create a new Conda
environment named `ct-dev` using:

```shell
conda create --name ct-dev --channel conda-forge libcantera-devel
```

C++ header and libraries are installed within the `ct-dev` environment folder, which
itself depends on the type of `conda` installation, and is abbreviated as
`path/to/conda/envs` below. Within the `ct-dev` folder, locations follow `conda`
recommendations for a given operating system.

### Linux and macOS Systems

Installation folders are:

```shell
library files         path/to/conda/envs/ct-dev/lib
pkg-config            path/to/conda/envs/ct-dev/lib/pkgconfig
C++ headers           path/to/conda/envs/ct-dev/include
Fortran module files  path/to/conda/envs/ct-dev/include/cantera
samples               path/to/conda/envs/ct-dev/share/cantera/samples
data files            path/to/conda/envs/ct-dev/share/cantera/data
```

In addition to `libcantera-devel`, installation of additional packages is recommended:

```shell
conda activate ct-dev
conda install --channel conda-forge cmake scons pkg-config
```

C++ programs can be compiled according to instructions outlined in the [C++
Guide](/userguide/compiling-cxx). Sample folders for C, C++ and Fortran include
pre-configured instruction files to facilitate compilation using the build tools `SCons`
and `CMake`, for example:

```shell
cd /path/to/conda/envs/ct-dev/share/cantera/samples/cxx/demo
scons  # uses SConstruct; or
cmake . && cmake --build .  # uses CMakeLists.txt
```

In addition, individual C++ Cantera sample programs can also be compiled using the
`pkg-config` build system:

```shell
g++ demo.cpp -o demo $(pkg-config --cflags --libs cantera)
```

In all cases, the build process yields the executable `demo`, which is run as:

```shell
./demo
```

### Windows Systems

Installation folders are:

```shell
library files               path\to\conda\envs\ct-dev\Library\lib
C++ headers                 path\to\conda\envs\ct-dev\Library\include
samples                     path\to\conda\envs\ct-dev\share\cantera\samples
data files                  path\to\conda\envs\ct-dev\share\cantera\data
```

C++ programs can be compiled according to instructions outlined in the
[C++ Guide](/userguide/compiling-cxx). Sample folders for C and C++ programs include
preconfigured instruction files to facilitate compilation using the build tools `SCons`
and `CMake`, for example:

```pwsh
cd path\to\conda\envs\ct-dev\share\cantera\samples\cxx\demo
scons  # uses SConstruct; or
cmake . && cmake --build . --config Release  # uses CMakeLists.txt
```

Fortran 90 support is not provided for Windows.

### Upgrading from an earlier Cantera version

If you already have the Cantera development interface installed in a Conda environment
(named, for example, `ct-dev`), you can upgrade it to the latest version available by
running the commands:

```shell
conda activate ct-dev
conda update --channel conda-forge libcantera-devel
```

## Matlab Interface

:::{attention}
The *legacy* Matlab Cantera interface was discontinued and removed in Cantera 3.1. Users
requiring support of legacy Matlab Cantera code should continue using Cantera 3.0
packages, or migrate their code base to the
[experimental Matlab toolbox](../matlab/index) that is currently under development.
:::
