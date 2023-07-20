# General File Structure {#sec-yaml-general}

The top level of a %Cantera YAML input file is a mapping that defines different input
file sections. Each section consists of a list of mappings that define objects of the
same type, such as reactions, species, phases, or elements. The `phases` section of
an input file contains all of the phase definitions. Multiple sections
containing reaction, species, or element definitions can be used. The
specific names `reactions`, `species`, and `elements` are used as defaults when looking
for @ref sec-yaml-reactions, @ref sec-yaml-species, and @ref sec-yaml-elements to add to
a phase. %Cantera uses SI input units by default, although input values can be
provided using a number of different units, as described in the @ref sec-yaml-units
section.

A simple input file has the following structure:

``` yaml
description: |-
  Optional description of the input file.

# units: optional unit settings used for input data

phases:
- name: spam
  thermo: ideal-gas
  # Additional fields come after this
- name: green-eggs
  thermo: model-name
  # Additional fields come after this

species:
- name: A
  # Additional fields come after this
- name: B
  # Additional fields come after this
- name: C
  # Additional fields come after this

reactions:
- equation: A + B <=> C + D
  # Additional fields come after this
- equation: A + C <=> 2 D
  # Additional fields come after this
```
