# Cantera – Experimental CLib Interface

This directory and the associated `src/clib_experimental` folder contain an
experimental re-implementation of Cantera's traditional CLib interface.

## Code Generation

Run the following command from the Cantera root folder:

```
scons doxygen
python -m pip install -e interfaces/sourcegen
sourcegen --api=clib --output=.
scons build clib_experimental=y
```

A rudimentary test suite ensures that code performs as expected:

```
scons test-clib-experimental
scons test-clib-experimental-demo
```

## Status

The experimental CLib Interface is in preview and still missing many features
needed for parity with the traditional CLib interface.
