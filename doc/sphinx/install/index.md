# Install

The following instructions explain how to install the different Cantera interfaces on a
variety of platforms.

```{tip}
We highly recommend that all new users install the Python interface using
[Conda](conda).
```

## Installing the Cantera Python Interface

- If you don't already have Python installed (or already use Conda), the easiest way to
  install the Cantera Python interface on all operating systems is by using
  [Conda](sec-install-conda).
- If you already have a different Python installation, Cantera can be installed using
  [Pip](sec-install-pip).
- Ubuntu users can install the `cantera-python3` package from the
  [Cantera PPA](sec-install-ubuntu).
- Fedora / Enterprise Linux users can install `python3-cantera` using
  [dnf](sec-install-fedora-rhel).
- OpenSUSE users can install packages using [zypper](sec-install-opensuse).
- Gentoo users can install using [emerge](sec-install-gentoo).
- FreeBSD users can install using [pkg](sec-install-freebsd).
- If you want to use the current development version, or add features of your own, you
  should [compile Cantera from source](sec-compiling).

## Installing the Cantera Matlab Toolbox

- The Cantera Matlab toolbox can be installed on all operating systems using
  [Conda](sec-conda-matlab-interface).
- Windows users can use [MSI installer packages](sec-install-windows).
- macOS users can also install using a [.pkg installer](sec-install-macos).
- If you want to use the current development version, or add features of your own, you
  should [compile Cantera from source](sec-compiling).

## Installing the Cantera C++ Interface & Fortran 90 Module

- The Cantera development interface can be installed on all operating systems using
  [Conda](sec-conda-development-interface).
- Ubuntu users can install the `cantera-dev` package from the
  [Cantera PPA](sec-install-ubuntu).
- Fedora / Enterprise Linux users can install packages using
  [yum/dnf](sec-install-fedora-rhel).
- OpenSUSE users can install packages using [zypper](sec-install-opensuse).
- Gentoo users can install using [emerge](sec-install-gentoo).
- FreeBSD users can install using [pkg](sec-install-freebsd).
- Users of other Linux distributions should
  [compile Cantera from source](sec-compiling).

```{toctree}
:maxdepth: 1
:hidden:

conda
pip
ubuntu
fedora-rhel
opensuse
gentoo
freebsd
windows
macos
```

```{seealso}
Check the [FAQ](sec-faq-installation) for solutions to some common installation
problems.
```
