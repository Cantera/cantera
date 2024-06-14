# User Guide

## Introductory Tutorials

For those new to Cantera, we present here a set of short tutorials to familiarize you
with Cantera's basic functionality and capabilities, give some examples of how to work
Cantera within your preferred interface language, and demonstrate some basic
troubleshooting.

- [Getting Started with Python](python-tutorial)
- [](input-tutorial)
- [](reactor-tutorial)
- [](cxx-tutorial)
- [](compiling-cxx)

## Frequently asked questions

See the [FAQ](faq) for answers to some common issues that arise when using
Cantera. If your question isn't answered here, consider asking us on the
<a href="https://cantera.org/community.html#the-cantera-users-group">Cantera Users' Group</a>.

## Task Guides

The tutorials in this section are designed to help you accomplish a specific task
using Cantera, such as evaluating the ignition delay time for a fuel under different
conditions, or calculating the voltage of a Lithium-ion battery as it is discharged.

### Working with Input Data

- [](ck2yaml-tutorial)
- [](creating-mechanisms)
- [](thermobuild)
- [](input-errors)
- [](legacy2yaml-tutorial)

### Combustion Calculations
- [](flame-temperature)
- [](heating-value)

### Electrochemistry Calculations
- [](/examples/python/kinetics/lithium_ion_battery)

### Implementing Custom Models

- [](extensible-reactor)

## Advanced Resources

- For intermediate and advanced users, the [](/reference/index) section is an
  easily-searchable repository that describes the scientific models implemented by
  Cantera and documents the classes and functions used to access these models.
- The [](/examples/index) section provides demonstrations of how many Cantera features
  can be used to solve a range of different problems and often provide a good starting
  point for writing your own code.
- Finally, if you have trouble using Cantera and can't find an answer here in the
  documentation, please visit the
  <a href="https://cantera.org/community.html#the-cantera-users-group">Cantera Users' Group</a>.

```{toctree}
:hidden:

python-tutorial
input-tutorial
reactor-tutorial
cxx-tutorial
compiling-cxx

faq

ck2yaml-tutorial
creating-mechanisms
thermobuild
input-errors
legacy2yaml-tutorial

flame-temperature
heating-value

/examples/python/kinetics/lithium_ion_battery

extensible-reactor
```
