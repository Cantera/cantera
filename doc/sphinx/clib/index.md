# CLib Interface Documentation

```{caution}
The CLib API is an experimental part of Cantera and may be changed without notice.
```

Cantera has two CLib interfaces:

- The _traditional_ CLib is hand-coded and is not feature-complete.

- The [](../develop/clib-extensions) is currently in preview and will replace the
  _traditional_ CLib starting in Cantera 3.2.

## Traditional CLib

The (undocumented) traditional CLib is installed by default when building Cantera:

```bash
scons build
```

## Auto-generated CLib

The auto-generated CLib is {ct}`fully documented <clibPage>`. The API is available when
Cantera is built with the `clib_experimental=y` option:

```bash
scons build clib_experimental=y
```

See [](../develop/clib-extensions) for further information.

```{toctree}
:hidden:
:maxdepth: 1
```
