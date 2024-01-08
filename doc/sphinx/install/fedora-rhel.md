(sec-install-fedora-rhel)=
# Fedora Packages

RPM packages are provided for supported versions of Fedora Linux. Stable builds are
available in the official repositories and development builds are in a Community
Projects (COPR) repository.

:::{attention}
The Matlab interface is not available from this archive. To install the Matlab interface
on Fedora, you must install it using [conda](sec-conda-matlab-interface) or
[compile the source code](sec-compiling).
:::

As of Cantera 3.0.0, packages are available for currently supported releases of Fedora
Linux and Fedora Rawhide as well as Enterprise Linux 8.

The available packages are:

- `python3-cantera` - The Cantera Python module for Python 3.
- `cantera-devel` - Shared object libraries and header files for compiling your own C++
  and Fortran 90 programs that use Cantera.
- `cantera-common` - Cantera data files and example programs.
- `cantera-static` - Static libraries for C++ and Fortran 90 development.

Cantera is available in the official repositories for Fedora; no configuration changes
are required.

## Basic Installation

On Enterprise Linux, if not already enabled, add the "Extra Packages for Enterprise
Linux" (EPEL) repository:

```bash
$ dnf install epel-release
```

To install all of the Cantera packages:

```bash
$ dnf install python3-cantera cantera-devel
```

or install whichever subset you need by adjusting the above command. The
`cantera-common` package is installed as a dependency if any other Cantera packages are
selected.

:::{tip}
If you plan on using Cantera from Python, you may also want to install IPython (an
advanced interactive Python interpreter) and Matplotlib (a plotting library). Matplotlib
is required to run some of the Python examples. These packages can be installed with:

```bash
$ dnf install python3-matplotlib python3-ipython
```
:::

## Installing pre-release Cantera versions

To access the development builds of Cantera (alpha or beta versions of the next version
of Cantera), enable the COPR:

```bash
$ dnf copr enable fuller/cantera-test
```

then install as described above.
