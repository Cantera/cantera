(sec-install-gentoo)=
# Gentoo Packages

The Gentoo [sci-libs/cantera](https://packages.gentoo.org/packages/sci-libs/cantera)
package is provided using a main portage tree. Additionally the
[app-doc/cantera-docs](https://packages.gentoo.org/packages/app-doc/cantera-docs)
package is provided for offline Documentation API reference for Cantera package
libraries.

:::{attention}
The Matlab interface is not available from this archive. To install the Matlab interface
on Gentoo, you must install it using [conda](sec-conda-matlab-interface) or
[compile the source code](sec-compiling).
:::

The following interfaces and tools are installed by default:

- C++ Libraries and header files for compiling your own programs that use Cantera.
- [YAML tools](../userguide/ck2yaml-tutorial).
- Python module for Python 3 (`python` USE flag with appropriate `PYTHON_SINGLE_TARGET`;
  optional).

The following additional interface is available:

- Fortran library and module files for compiling your own programs that use Cantera
  (`fortran` USE flag; optional)

More information about `USE` flags can be found in the
[Gentoo Handbook](https://wiki.gentoo.org/wiki/Handbook:Parts/Working/USE).
To learn about per-package control of `USE` flags, please refer to the
[/etc/portage/package.use](https://wiki.gentoo.org/wiki//etc/portage/package.use)
article.

To install the `sci-libs/cantera` and `app-doc/cantera-docs` packages:

```bash
emerge --ask cantera cantera-docs
```

Most likely, the latest versions of these packages and/or some of their dependencies
still have unstable status in the Gentoo portage tree. In this case, you have to
`unmask` (allow installation within stable system) them preliminarily using
[/etc/portage/package.accept_keywords](https://wiki.gentoo.org/wiki//etc/portage/package.accept_keywords).

If `/etc/portage/package.accept_keywords` is present in your system as file then (for
64-bit architecture) you could unmask `sci-libs/cantera` package by running the
following command (as root):

```bash
echo "sci-libs/cantera ~amd64" >> /etc/portage/package.accept_keywords
```

Otherwise, if `/etc/portage/package.accept_keywords` is present in your system as a
directory, run the following command (as root):

```bash
echo "sci-libs/cantera ~amd64" >> /etc/portage/package.accept_keywords/cantera
```

:::{tip}
If you plan on using Cantera from Python, you may also want to install IPython
([dev-python/ipython](https://packages.gentoo.org/packages/dev-python/ipython), an
advanced interactive Python interpreter) and Matplotlib
([dev-python/matplotlib](https://packages.gentoo.org/packages/dev-python/matplotlib), a
plotting library). Matplotlib is required to run some of the Python examples. These
packages can be installed with:

```bash
emerge --ask ipython matplotlib
```

:::
