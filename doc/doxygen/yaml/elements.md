# Element Entries {#sec-yaml-elements}

[TOC]

## Overview

%Cantera provides built-in definitions for the chemical elements, including values for
their atomic weights taken from IUPAC / CIAAW. These elements can be used by specifying
the corresponding atomic symbols when specifying the composition of species. `element`
entries are needed only when defining custom elements that are not standard chemical
elements, or defining specific isotopes.

In order to give a name to a particular isotope or a virtual element
representing a surface site, a custom `element` entry can be used. The
default location for `element` entries is the `elements` section of the
input file. Elements defined in this section will automatically be
considered for addition to phases defined in the same file. Elements can
be defined in other sections of the input file if those sections are
named explicitly in the `elements` field of the phase definition.

**Example:**

``` yaml
elements:
- symbol: C13
  atomic-weight: 13.003354826
  atomic-number: 12
- symbol: O-18
  atomic-weight: 17.9991603
```

## Element API Reference

The fields of an `element` entry are:

<b>`symbol`</b>

The symbol used for the element, as used when specifying the
composition of species.

<b>`atomic-weight`</b>

The atomic weight of the element, in unified atomic mass units
(dalton).

<b>`atomic-number`</b>

The atomic number of the element. Optional.

<b>`entropy298`</b>

The standard molar entropy of the element at 298.15 K. Optional.
