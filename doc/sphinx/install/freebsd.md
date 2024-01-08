(sec-install-freebsd)=
# FreeBSD Packages

A community-maintained FreeBSD port ``science/cantera`` and package are available. This
port provides the Octave interface when the OCTAVE option is turned on.

:::{attention}
The Matlab interface is not available from this port; to install the Matlab interface on
FreeBSD, you must install it using [conda](sec-conda-matlab-interface) or
[compile the source code](sec-compiling).
:::

The package can be installed by running

```shell
pkg install cantera
```

Further information about the Cantera package can be found on
[FreeBSD.org](https://www.freebsd.org/cgi/ports.cgi?query=cantera&stype=all) and [FreshPorts.org](https://www.freshports.org/science/cantera/).
