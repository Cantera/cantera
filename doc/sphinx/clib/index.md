# C Interface Documentation

```{caution}
The CLib API is an experimental part of Cantera and may be changed without notice.
```

The [](../develop/clib-extensions) consists of a set of libraries that are designed to
encapsulate Cantera functionality and to make it available for use in languages and
applications other than C++. Each library provides functions that are declared
`extern "C"`. All Cantera objects are stored and referenced by integers - no pointers
are passed to or from the calling application.

## Generated C Interface

:::{versionadded} 3.2
Replaces the _legacy_ C interface starting in Cantera 3.2.
:::

The generated CLib is {ct}`fully documented in Doxygen <CAPImain>` and is installed by
default.

While the generated CLib is not feature-complete compared to the C++ API, it
re-implements all features of the legacy CLib and is easily extensible. See
[](../develop/clib-extensions) for further information.

```{toctree}
:hidden:
:maxdepth: 1
```
