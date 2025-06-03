# CLib Interface Documentation

```{caution}
The CLib API is an experimental part of Cantera and may be changed without notice.
```

Cantera has two CLib interfaces:

- The [](../develop/clib-extensions) replaces the _legacy_ CLib starting in Cantera 3.2.

- The _legacy_ CLib is hand-coded and is not feature-complete.

## Auto-Generated CLib

The auto-generated CLib is {ct}`fully documented <clibPage>`. The API is installed by
default when building Cantera:

```bash
scons build
```

While the auto-generated CLib is not feature-complete compared to the C++ API, it
re-implements all features of the legacy CLib and is easily extensible. See
[](../develop/clib-extensions) for further information.

## Legacy CLib

The (undocumented) legacy CLib is available when Cantera is built with the
`clib_legacy=y` option:

```bash
scons build clib_legacy=y
```

```{toctree}
:hidden:
:maxdepth: 1
```
