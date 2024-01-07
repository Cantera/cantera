(sec-install-pip)=
# Installing with Pip

[Pip](https://pip.pypa.io/en/stable/) is a package installer for Python that can be used
to install the Cantera Python module from [PyPI](https://pypi.org/project/Cantera/).

:::{caution}
There are a few important limitations to note when Cantera is installed from PyPI:

- These packages are compiled without native HDF5 support. The new options for saving
  and restoring `SolutionArray` and flame objects to/from HDF5 files is not available.
- These packages are compiled using single-threaded implementations of LAPACK functions,
  and cannot make use of multiple cores to speed up reactor network or flame
  simulations.

If you want either of these features, you can install the [Conda](conda) packages
instead.
:::

## Prerequisites

The first step in installing the Cantera Python module using `pip` is to make sure you
have a compatible version of Python installed and are able to run `pip` from the command
line. Packages for Cantera 3.0.0 are available for Python versions 3.8, 3.9, 3.10, and
3.11.

If you don't already have Python installed, it can be downloaded from
[python.org](https://www.python.org/) or installed using your operating system's package
manager.

To check that you are able run `pip`, open a terminal / command prompt and run the
following command:

:::::{tab-set}
::::{tab-item} Linux / macOS
:sync: posix
```shell
python3 -m pip --version
```
::::

::::{tab-item} Windows
:sync: windows
```bat
py -m pip --version
```
::::
:::::

If the above command doesn't work, see the instructions at
[packaging.python.org](https://packaging.python.org/en/latest/tutorials/installing-packages/)
for how to get `pip` working with your Python installation.

## Virtual Environments

Virtual environments provide a way keeping separate sets of Python packages installed
for different projects, where different environments can have different versions of
packages that might otherwise conflict. To create and activate a virtual environment
named `ct-env` to be used with Cantera, run the commands:

:::::{tab-set}
::::{tab-item} Linux / macOS
:sync: posix
```shell
python3 -m venv ct-env
source ct-env/bin/activate
```
::::

::::{tab-item} Windows
:sync: windows
```bat
py -m venv ct-env
ct-env\Scripts\activate
```
::::
:::::

The second command should be run in the terminal each time you want to use the specified
environment.

## Installing Cantera

To install the Cantera Python module, first activate your virtual environment, if you're
using one. Then, run the command:

:::::{tab-set}
::::{tab-item} Linux / macOS
:sync: posix
```shell
python3 -m pip install cantera
```
::::

::::{tab-item} Windows
:sync: windows
```bat
py -m pip install cantera
```
::::
:::::

You can test that your installation is working by importing the Cantera module and
creating a `Solution` object:

:::::{tab-set}
::::{tab-item} Linux / macOS
:sync: posix
```shell
python3 -c 'import cantera; gas = cantera.Solution("gri30.yaml"); gas()'
```
::::

::::{tab-item} Windows
:sync: windows
```bat
py -c "import cantera; gas = cantera.Solution('gri30.yaml'); gas()"
```
::::
:::::

You should get the following output:

```
 gri30:

      temperature   300 K
         pressure   1.0133e+05 Pa
          density   0.081894 kg/m^3
 mean mol. weight   2.016 kg/kmol
  phase of matter   gas

                         1 kg             1 kmol
                    ---------------   ---------------
         enthalpy             26469             53361  J
  internal energy       -1.2108e+06        -2.441e+06  J
          entropy             64910        1.3086e+05  J/K
   Gibbs function       -1.9447e+07       -3.9204e+07  J
heat capacity c_p             14311             28851  J/K
heat capacity c_v             10187             20536  J/K

                     mass frac. Y      mole frac. X     chem. pot. / RT
                    ---------------   ---------------   ---------------
               H2                 1                 1           -15.717
    [  +52 minor]                 0                 0
```

### Installing a Pre-release

Sometimes, a pre-release (alpha or beta) version of Cantera may be available to install.
However, `pip` defaults to installing the latest stable version. To allow installation
of a pre-release, add the `--pre` flag:

:::::{tab-set}
::::{tab-item} Linux / macOS
:sync: posix
```shell
python3 -m pip install --pre cantera
```
::::

::::{tab-item} Windows
:sync: windows
```bat
py -m pip install --pre cantera
```
::::
:::::

You can check the version that was installed by running:

:::::{tab-set}
::::{tab-item} Linux / macOS
:sync: posix
```shell
python3 -c 'import cantera; print(cantera.__version__)'
```
::::

::::{tab-item} Windows
:sync: windows
```bat
py -c "import cantera; print(cantera.__version__)"
```
::::
:::::
