@page clibPage Auto-Generated CLib Interface

# %Cantera – Auto-Generated CLib Interface

> **Note:** The auto-generated CLib API is a re-implementation of %Cantera's
> _traditional_ CLib interface. It is currently in preview and will replace the
> traditional CLib API in %Cantera 3.2.
>
> While many function names remain unchanged, the auto-generated CLib is not meant to
> be a drop-in replacement. Breaking changes involve improved consistency with the C++
> code base in terms of nomenclature and function/method signatures.

## CLib Modules

%Cantera classes and CLib modules have a one-to-one relationship. The @ref clibGroup
page provides a list of currently available modules.

## Code Generation

The auto-generated CLib interface is available if %Cantera is built with the
`clib_experimental=y` option:

```
scons build clib_experimental=y
```

A test suite ensures that code performs as expected:

```
scons test-clib
scons test-clib-demo
```
