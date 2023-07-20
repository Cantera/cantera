@page sec-yaml-format-tutorial YAML Format Tutorial

[TOC]

Here we describe the syntax and structure of %Cantera YAML files, how
dimensional values in %Cantera YAML files are handled, and how to
understand some of the error messages that may be encountered when
reading these input files.

# Syntax

%Cantera YAML files use a subset of the [YAML 1.2](https://yaml.org/spec/1.2/spec.html)
specification. %Cantera YAML files consist of individual values, which may be strings,
numbers or booleans, that are then composed as elements of nested mappings and
sequences.

# Strings

Strings may be generally written without quotes, but may be enclosed in
single quotes or double quotes if needed in order to avoid certain
parsing ambiguities.

```yaml
A string
Another 'string'
"A string: that requires quotes"
```

# Numbers

Numbers can be written as integers, decimal values, or using E-notation

```yaml
3
3.14
6.022e23
```

# Booleans

Boolean values in YAML are written as the words `true` or `false`.

# Sequences

A sequence of multiple items is specified by separating the items by
commas and enclosing them in square brackets. The individual items can
have any type \--strings, integers, floating-point numbers, mappings, or
sequences.

```yaml
elements: [O, H, C, N, Ar]
temperature-ranges: [200.0, 1000.0, 3500.0]
```

The syntax above, using square brackets to define a list, is called
**flow style** in YAML. Sequences can also be written in **block style**,
using one line for each item in the sequence, with each line starting with a dash:

```yaml
elements:
- O
- H
- C
- N
- Ar
```

Sequences can also be nested. The following examples are all equivalent:

```yaml
data: [[1, 2], [3, 4]]

data:
-
  - 1
  - 2
-
  - 3
  - 4

data:
- - 1
  - 2
- - 3
  - 4
```

# Mappings

A mapping is a container consisting of key\--value pairs. The keys in a
mapping must be unique. Like sequences, there are two ways to write a
mapping. In the **flow style**, the mapping is enclosed in curly
brackets, colons (followed by spaces) are used to separate keys and
values, and key\--value pairs are separated by commas:

```yaml
composition: {H: 2, C: 1, O: 1}
```

In the **block style**, each key is written on a new line, followed by a
colon. The value can be placed either on the same line, or on the
following line, indented one level:

```yaml
composition:
  H: 2
  C:
    1
  O: 1
```

All keys in %Cantera YAML files are treated as strings. A %Cantera YAML
file is itself a mapping, usually in the **block style**. We refer to
the keys in this top-level mapping as the **sections** of the input
file.

# Sequences of Mappings

A common structure in %Cantera input files is a nested sequence of
mappings. This can be written in the **block style** as:

```yaml
- equation: O2 + CO <=> O + CO2
  rate-constant: {A: 2.5e+12, b: 0, Ea: 47800}
- equation: O2 + CH2O <=> HO2 + HCO
  rate-constant: {A: 1.0e+14, b: 0, Ea: 40000}
- equation: H + O2 + M <=> HO2 + M
  type: three-body
  rate-constant: {A: 2.8e+18, b: -0.86, Ea: 0}
  efficiencies: {AR: 0, C2H6: 1.5, CO: 0.75, CO2: 1.5, H2O: 0, N2: 0, O2: 0}
```

The keys in each mapping need not be the same. In this example, each of
the three mappings in the sequence has `equation` and `rate-constant`
keys, while only the third entry has `type` and `efficiencies` keys.

# Comments

The character `#` is the comment character. Everything to the right of
this character on a line is ignored:

```yaml
# set the default units
units:
  length: cm  # use centimeters for length
  quantity: mol  # use moles for quantity
```

# Dimensional Values

Many fields have numerical values that represent dimensional
quantities\-\--a pressure, or a density, for example. If these are
entered without specifying the units, the default units (set by the
`units` directive) will be used. However, it is also possible to specify
the units for each individual dimensional quantity, unless stated
otherwise. All that is required is to write the units after the value,
separated by a space:

```yaml
pressure: 1.0e5  # default is Pascals
pressure: 1.0 bar  # this is equivalent
density: 4.0 g/cm^3
density: 4000.0  # kg/m^3
```

Compound unit strings may be used, as long as a few rules are followed:

1.  Units in the denominator follow `/`.
2.  Units in the numerator follow `*`, except for the first one.
3.  Numerical exponents follow the unit string with a `^` character.

Examples of compound units:

```yaml
A: 1.0e20 cm^6/mol^2/s  # OK
h: 6.626e-34 J*s  # OK
density: 3.0 g*cm^-3  # OK
A: 1.0e20 cm6/mol/s  # error (missing '^')
A: 1.0e20 cm^6/mol^2-s  # error ('s' should be in denominator)
density: 3.0g/cm^3  # error (missing space between value and units)
```

See the [Units
API](%7B%7B%%20ct_docs%20sphinx/html/yaml/general.html#units%20%%7D%7D)
documentation for additional details, including the full set of
supported units.

# Default Units

Default units that apply to a whole input file or some portion thereof
can be set using `units` mapping. A `units` mapping placed at the top
level of an input file applies to the entire file. A `units` mapping
placed as a member of another mapping applies to that mapping and any
nested mappings or sequences, and overrides higher-level `units`
mappings:

```yaml
units: {length: cm, mass: kg}
section1:
  units: {length: m}
  density: 4000  # interpreted as 4000 kg/m^3
section2:
  density: 0.1  # interpreted as 0.1 kg/cm^3
section3:
- units: {mass: mg}  # must be the first item in the list
- name: species1
  density: 5e4  # interpreted as 5e4 mg/cm^3
```

Default units may be set for `mass`, `length`, `time`, `quantity`,
`pressure`, `energy`, and `activation-energy`.

# Error Handling

During processing of an input file, errors may be encountered. These
could be syntax errors, or could be ones that are flagged as errors by
%Cantera due to some apparent inconsistency in the data\-\--an unphysical
value, a species that contains an undeclared element, a reaction that
contains an undeclared species, missing species or element definitions,
multiple definitions of elements, species, or reactions, and so on.

# Syntax Errors

Syntax errors are caught by the YAML parser, and must be corrected
before proceeding further. If a syntax error is encountered, %Cantera
will raise an exception which includes the location of the error.
Additional information such as a traceback showing where in the code the
input file was being read may be printed as well.

For example, consider the following input file, which is intended to
create a gas with the species and reactions of GRI-Mech 3.0, but is
missing the colon which is needed after the `thermo` key:

```yaml
phases:
- name: gas
  thermo ideal-gas
  kinetics: gas
  elements: [H, O]
  species: [{gri30.yaml/species: all}]
  reactions: [gri30.yaml/reactions]
```

When this definition is imported into an application, an error message
like the following would be printed to the screen, and execution of the
program or script would terminate:

``` python
Traceback (most recent call last):
  File "<stdin>", line 1, in <module>
  File "/some/path/cantera/base.pyx", line 25, in cantera._cantera._SolutionBase.__cinit__
    self._init_yaml(infile, phaseid, phases, yaml)
  File "/some/path/cantera/base.pyx", line 49, in cantera._cantera._SolutionBase._init_yaml
    root = AnyMapFromYamlFile(stringify(infile))
cantera._cantera.CanteraError:
***********************************************************************
InputFileError thrown by AnyMap::fromYamlFile:
Error on line 4 of ./gas.yaml:
illegal map value
|  Line |
|     1 | phases:
|     2 | - name: gas
|     3 |   thermo ideal-gas
>     4 >   kinetics: gas
                    ^
|     5 |   elements: [H, O]
|     6 |   species: [{gri30.yaml/species: all}]
|     7 |   reactions: [gri30.yaml/reactions]
***********************************************************************
```

The top part of the error message shows the chain of functions that were
called before the error was encountered. For the most part, these are
internal %Cantera functions not of direct concern here. The relevant part
of this error message is the part between the lines of asterisks. This
message says that the YAML parser ran into a problem on line 4 of
`gas.yaml`. In many cases, including this one, the parser will fail
somewhere *after* the actual problem with the input file, since it must
continue parsing until it finds something that cannot possibly be valid
YAML syntax. In this case, the problem from the parser's perspective is
that the key which started on line 3 continues across a new line before
it finds a colon that can be considered as the separator. Since a key
can't be broken across lines like this, the parser indicates the error
at the point where it found the colon. By looking back from the
indicated point of the error, we can see that the problem is the missing
colon in the previous line.

# Cantera Errors

Now let's consider the other class of errors, ones that %Cantera itself
detects. Continuing the example above, suppose that the missing colon is
corrected, and the input file processed again. Again an error message
results, but this time it is from %Cantera:

``` python
Traceback (most recent call last):
  File "<stdin>", line 1, in <module>
  File "/some/path/cantera/base.pyx", line 25, in cantera._cantera._SolutionBase.__cinit__
    self._init_yaml(infile, phaseid, phases, yaml)
  File "/some/path/cantera/base.pyx", line 49, in cantera._cantera._SolutionBase._init_yaml
    root = AnyMapFromYamlFile(stringify(infile))
cantera._cantera.CanteraError:
***********************************************************************
CanteraError thrown by Phase::addSpecies:
Species 'C' contains an undefined element 'C'.
***********************************************************************
```

The problem is that the phase definition specifies that all species are
to be imported from the `gri30` mechanism, but only the elements H and O
are declared. The `gri30` mechanism contains species composed of the
elements H, O, C, N, and Ar. If the definition is modified to declare
these additional elements:

```yaml
phases:
- name: gas
  thermo: ideal-gas
  kinetics: gas
  elements: [H, O, C, N, Ar]
  species: [{gri30.yaml/species: all}]
  reactions: [gri30.yaml/reactions]
```

it may be imported successfully.
