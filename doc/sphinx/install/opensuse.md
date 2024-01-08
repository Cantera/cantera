(sec-install-opensuse)=
# openSUSE Packages

Cantera 3.0.0 packages are available for openSUSE Tumbleweed from a
[community repository](https://build.opensuse.org/package/show/home:fuller/cantera).

Installation is as follows:

```bash
$ zypper addrepo https://download.opensuse.org/repositories/home:fuller/openSUSE_Tumbleweed/home:fuller.repo
$ zypper refresh
$ zypper install cantera
```

:::{attention}
The Matlab interface is not available from this archive. To install the Matlab interface
on openSUSE, you must install it using [conda](sec-conda-matlab-interface) or
[compile the source code](sec-compiling).
:::
