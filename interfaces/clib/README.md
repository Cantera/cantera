@page CAPImain Generated C Interface

# Generated C Interface (CLib)

@remark  The generated CLib API is a re-implementation of %Cantera's _legacy_ C
interface. It replaces the legacy CLib API starting in %Cantera 3.2.

@remark  While many function names remain unchanged, the generated CLib is not meant to
be a drop-in replacement. Breaking changes involve improved consistency with the C++
code base in terms of nomenclature and function/method signatures.

The generated C interface consists of a set of libraries that are designed to
encapsulate %Cantera functionality and to make it available for use in languages and
applications other than C++. Each library provides functions that are declared
`extern "C"`. All %Cantera objects are stored and referenced by integers - no pointers
are passed to or from the calling application.

## CLib Libraries

@warning  The generated CLib API is an experimental part of %Cantera and may be changed
without notice.

The @ref CAPIindex page provides a list of currently available libraries.

## Code Generation

The generated CLib interface is installed by default if %Cantera was built with default options. A test suite ensures that code performs as expected:

```
scons test-clib
scons test-clib-demo
```
