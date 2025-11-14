(sec-install-ubuntu)=
# Ubuntu Packages

As of Cantera 3.2.0, packages are available for Ubuntu 24.04 (Noble Numbat) and Ubuntu
25.10 (Questing Quokka). Generally, packages are available for the most-recent LTS
release and the latest non-LTS release. Packages for other Ubuntu releases supported by
[Launchpad](https://launchpad.net/ubuntu) are provided on a best-effort basis. To see
which Ubuntu releases and Cantera versions are currently supported, visit the [Cantera
PPA](https://launchpad.net/~cantera-team/+archive/ubuntu/cantera).

The available packages are:

- `python3-cantera` - The Cantera Python module for Python 3.
- `libcantera-dev` - Libraries and header files for compiling your own C, C++ and
  Fortran 90 programs that use Cantera.
- `cantera-common` - Cantera data files and example programs
  - Samples for all language interfaces are installed to `/usr/share/cantera/samples/`
  - Data files are installed in `/usr/share/cantera/data/`. This location is
    automatically added to the Cantera data path.
- `libcantera3.2` - The Cantera C++ library, for use by packaged C++ applications.
- `libcantera-fortran3.2` - The Cantera Fortran 90 library, for use by packaged
  Fortran 90 applications.
- `cantera` - A metapackage that will install everything except for the development
  files.

## Installing

To add the Cantera PPA:

```bash
sudo apt install software-properties-common
sudo apt-add-repository ppa:cantera-team/cantera
```

To install all of the Cantera packages:

```bash
sudo apt install python3-cantera libcantera-dev
```

or install whichever subset you need by adjusting the above command. The
`cantera-common` package is installed as a dependency if any other Cantera packages are
selected.

:::{tip}
If you plan on using Cantera from Python, you may also want to install IPython (an
advanced interactive Python interpreter) and Matplotlib (a plotting library). Matplotlib
is required to run some of the Python examples. These packages can be installed with:

```bash
sudo apt install python3-pip
pip3 install ipython matplotlib
```
:::

## Upgrading from an earlier Cantera version

If you already have Cantera installed from the `cantera-team` PPA, you can ensure that
you have the latest available version installed by running:

```bash
sudo apt update
sudo apt install python3-cantera
```

If you also have the `libcantera-dev` package installed, it should also be included on
the `apt install` command line.

## Installing pre-release Cantera versions

Sometimes, pre-release (alpha or beta) versions of Cantera which represent work toward
the next Cantera release will be available for users who want to use cutting-edge
features or test compatibility with the new version before it is released. To see the
latest Cantera versions available from this PPA, visit
<https://launchpad.net/~cantera-team/+archive/ubuntu/cantera-unstable>.

These packages can be installed by additionally enabling the
`cantera-team/cantera-unstable` PPA and then upgrading Cantera:

```bash
sudo apt-add-repository ppa:cantera-team/cantera-unstable
sudo apt install cantera-python3 libcantera-dev
```

You should also have the `cantera-team/cantera` PPA enabled, since the
`cantera-unstable` PPA *only* includes development versions.

If you later want to remove the development version and return to the latest stable
version, run the commands:

```bash
sudo apt-add-repository --remove ppa:cantera-team/cantera-unstable
sudo apt remove cantera cantera-common libcantera-dev python3-cantera
sudo apt install python3-cantera libcantera-dev
```
