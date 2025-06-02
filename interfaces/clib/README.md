@page clibPage Experimental CLib Interface

# %Cantera – Experimental CLib Interface

This directory and the associated `src/clib_experimental` folder contain an
experimental re-implementation of Cantera's traditional CLib interface.

## Code Generation

The auto-generated CLib interface is available if %Cantera is built with the the
`clib_experimental=y` option:

```
scons build clib_experimental=y
```

A test suite ensures that code performs as expected:

```
scons test-clib
scons test-clib-demo
```

## CLib Modules

%Cantera classes and CLib modules have a one-to-one relationship. The @ref clibGroup
page provides a list of currently available classes.
