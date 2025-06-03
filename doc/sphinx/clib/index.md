# C Interface Documentation

```{caution}
The CLib API is an experimental part of Cantera and may be changed without notice.
```

Cantera has two C interfaces (CLib):

- The [](../develop/clib-extensions).

- The _legacy_ CLib is hand-coded and is not feature-complete.

Both C interfaces consists of a set of modules that are designed to encapsulate
Cantera functionality and to make it available for use in languages and applications
other than C++. Each modules provides a library of functions that are declared
`extern "C"`. All Cantera objects are stored and referenced by integers - no
pointers are passed to or from the calling application.

## Generated C Interface

:::{versionadded} 3.2
Replaces the _legacy_ C interface starting in Cantera 3.2.
:::

The generated CLib is {ct}`fully documented in Doxygen <clibPage>` and is installed by
default.

While the generated CLib is not feature-complete compared to the C++ API, it
re-implements all features of the legacy CLib and is easily extensible. See
[](../develop/clib-extensions) for further information.

## Legacy C Interface

:::{deprecated} 3.2
The legacy CLib is superseded by the generated CLib and will be removed in Cantera 3.3.
:::

The (undocumented) legacy CLib is available when Cantera is built with the
`clib_legacy=y` option:

```bash
scons build clib_legacy=y
```

```{toctree}
:hidden:
:maxdepth: 1
```
