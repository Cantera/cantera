# Introduction to YAML Input Files

All calculations in Cantera require an input file to describe the properties of the
relevant phase(s) of matter. These input files can be provided via one of several
methods:

- *Use an existing input file*. The input files included with Cantera should suffice for
  tutorials and introductory work with thermodynamic phases and reaction kinetics. A
  number of additional input files are available from various sources on the web.

- *Convert a mechanism from Chemkin (CK) format to YAML format*. Many reaction
  mechanisms are published in Chemkin (CK) format, and can be converted to Cantera's
  YAML format using the `ck2yaml` conversion tool. see [](ck2yaml-tutorial) for more
  information.

- *Create your own YAML file from scratch or by editing an existing file*. Advanced
  users may need to edit an existing input file in order to define additional species,
  reactions, or entirely new phases. See [](creating-mechanisms) for more information.

```{tip}
Whenever you edit a Cantera input file, it is *highly advised* that you begin by copying
the existing file and saving it under a new name, before editing the new file. Editing a
file under its original name can easily lead to errors, if one forgets that this file
does not represent the original mechanism.
```

## Understanding Cantera's Input File Syntax

With any input file, it can be helpful to understand the syntax requirements. Clearly,
anyone writing directly in the YAML formats must conform to these standards. However,
even when importing an externally-provided file or converting from CK format,
understanding the input file syntax can help diagnose and correct any errors (although
many/most of the CK conversion errors will be related to errors in the CK syntax
formatting).

Cantera YAML files use a subset of the [YAML 1.2](https://yaml.org/spec/1.2/spec.html)
specification. Cantera YAML files consist of individual values, which may be strings,
numbers or booleans, that are then composed as elements of nested mappings and
sequences.

### Strings

Strings may be generally written without quotes, but may be enclosed in single
quotes or double quotes if needed in order to avoid certain parsing ambiguities.

```yaml
A string
Another 'string'
"A string: that requires quotes"
```

### Numbers

Numbers can be written as integers, decimal values, or using E-notation

```yaml
3
3.14
6.022e23
```

### Booleans

Boolean values in YAML are written as the words `true` or `false`.

### Sequences

A sequence of multiple items is specified by separating the items by commas and
enclosing them in square brackets. The individual items can have any type -- strings,
integers, floating-point numbers, mappings, or sequences.

```yaml
elements: [O, H, C, N, Ar]
temperature-ranges: [200.0, 1000.0, 3500.0]
```

The syntax above, using square brackets to define a list, is called **flow style** in
YAML. Sequences can also be written in **block style**, using one line for each item in
the sequence, with each line starting with a dash:

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

### Mappings

A mapping is a container consisting of key--value pairs. The keys in a mapping must be
unique. Like sequences, there are two ways to write a mapping. In the **flow style**,
the mapping is enclosed in curly brackets, colons (followed by spaces) are used to
separate keys and values, and key--value pairs are separated by commas:

```yaml
composition: {H: 2, C: 1, O: 1}
```

In the **block style**, each key is written on a new line, followed by a colon. The
value can be placed either on the same line, or on the following line, indented one
level:

```yaml
composition:
  H: 2
  C:
    1
  O: 1
```

All keys in Cantera YAML files are treated as strings. A Cantera YAML file is itself a
mapping, usually in the **block style**. We refer to the keys in this top-level mapping
as the **sections** of the input file.

### Sequences of Mappings

A common structure in Cantera input files is a nested sequence of mappings. This can be
written in the **block style** as:

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

The keys in each mapping need not be the same. In this example, each of the three
mappings in the sequence has `equation` and `rate-constant` keys, while only the third
entry has `type` and `efficiencies` keys.

### Comments

The character `#` is the comment character. Everything to the right of this character on
a line is ignored:

```yaml
# set the default units
units:
  length: cm  # use centimeters for length
  quantity: mol  # use moles for quantity
```

### Dimensional Values

Many fields have numerical values that represent dimensional quantities---a pressure, or
a density, for example. If these are entered without specifying the units, the default
units (set by the `units` directive) will be used. However, it is also possible to
specify the units for each individual dimensional quantity, unless stated otherwise. All
that is required is to write the units after the value, separated by a space:

```yaml
pressure: 1.0e5  # default is Pascals
pressure: 1.0 bar  # this is equivalent
density: 4.0 g/cm^3
density: 4000.0  # kg/m³
```

Compound unit strings may be used, as long as a few rules are followed:

1. Units in the denominator follow `/`.
2. Units in the numerator follow `*`, except for the first one.
3. Numerical exponents follow the unit string with a `^` character.

Examples of compound units:

```yaml
A: 1.0e20 cm^6/mol^2/s  # OK
h: 6.626e-34 J*s  # OK
density: 3.0 g*cm^-3  # OK
A: 1.0e20 cm6/mol/s  # error (missing '^')
A: 1.0e20 cm^6/mol^2-s  # error ('s' should be in denominator)
density: 3.0g/cm^3  # error (missing space between value and units)
```

See the [Units API](sec-yaml-units) documentation for additional details, including the
full set of supported units.

### Default units

Default units that apply to a whole input file or some portion thereof can be set using
`units` mapping. A `units` mapping placed at the top level of an input file applies to
the entire file. A `units` mapping placed as a member of another mapping applies to that
mapping and any nested mappings or sequences, and overrides higher-level `units`
mappings:

```yaml
units: {length: cm, mass: kg}
section1:
  units: {length: m}
  density: 4000  # interpreted as 4000 kg/m³
section2:
  density: 0.1  # interpreted as 0.1 kg/cm^3
section3:
- units: {mass: mg}  # must be the first item in the list
- name: species1
  density: 5e4  # interpreted as 5e4 mg/cm^3
```

Default units may be set for `mass`, `length`, `time`, `quantity`, `pressure`, `energy`,
and `activation-energy`.

## Input Files Distributed with Cantera

Several reaction mechanism files are included in the Cantera distribution, including
ones that model natural gas combustion [`gri30.yaml`](../examples/input/gri30),
high-temperature air (`air.yaml`), a hydrogen/oxygen reaction mechanism
[`h2o2.yaml`](../examples/input/h2o2), some pure fluids in the liquid-vapor region
(`liquidvapor.yaml`), and a few surface reaction mechanisms (such as `ptcombust.yaml`,
[`diamond.yaml`](../examples/input/diamond), etc.), among others. A subset of these data
files are shown in the [Examples](../examples/input/index) section to demonstrate the
YAML format and some of the commonly-used models and options.

Additional data used in various examples are included in the `example_data` subdirectory
of the default Cantera data directory. These files can be imported using relative paths,
for example (in Python):

```py
gas = ct.Solution("example_data/ammonia-CO-H2-Alzueta-2023.yaml")
```

```{caution}
The input files included with Cantera are provided for convenience, and may not be
suited for research purposes.
```

On Windows, these files may be located in `C:\Program Files\Cantera\data` depending on
how you installed Cantera and the options you specified. On a Unix/Linux/macOS machine,
they are usually kept in the `data` subdirectory within the Cantera installation
directory. You can also browse the list of data files included in the base
[`data`](https://github.com/Cantera/cantera/tree/main/data) directory and the
[`example_data`](https://github.com/Cantera/cantera-example-data/tree/main/)
subdirectory in the respective source repositories.

% TODO: add link to the Matlab tutorial
Please see the tutorials for [Python](python-tutorial) and Matlab for
instructions on how to import from these pre-existing files.
