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
- Gentoo users can install using [emerge](sec-install-gentoo).
- FreeBSD users can install using [pkg](sec-install-freebsd).
- If you want to use the current development version, or add features of your own, you
  should [compile Cantera from source](sec-compiling).

## Installing the Cantera C++ Interface & Fortran 90 Module

- The Cantera development interface can be installed on all operating systems using
  [Conda](sec-conda-development-interface).
- Ubuntu users can install the `cantera-dev` package from the
  [Cantera PPA](sec-install-ubuntu).
- Gentoo users can install using [emerge](sec-install-gentoo).
- FreeBSD users can install using [pkg](sec-install-freebsd).
- Users of other Linux distributions should
  [compile Cantera from source](sec-compiling).

## Installing the Cantera Matlab Toolbox

:::{attention}
The *legacy* Matlab Cantera interface was discontinued and removed in Cantera 3.1. Users
requiring support of legacy Matlab Cantera code should continue using Cantera 3.0
packages, or migrate their code base to the
[experimental Matlab toolbox](../matlab/index) that is currently under development.
:::

```{toctree}
:maxdepth: 1
:hidden:

conda
pip
ubuntu
gentoo
freebsd
```

## Troubleshooting

```{seealso}
Check the [FAQ](sec-faq-installation) for solutions to some common installation
problems.
```

```{attention}
Packaged versions of Cantera for Fedora, Enterprise Linux (RHEL), and OpenSUSE are no
longer available due to the lack of a package maintainer for these distributions. If you
have experience with any of these packaging systems and are interested in maintaining
Cantera packages, the Cantera developers would appreciate your help.
See [Issue #1985](https://github.com/Cantera/cantera/issues/1985).

Installation using Conda and Pip on these distributions is still supported, along with
compilation from source.
```
