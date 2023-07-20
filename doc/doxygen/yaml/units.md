
@page sec-yaml-units Input Units

[TOC]

While %Cantera generally works internally in SI units, input values can
be provided using a number of different units.

# Overview

Compound units are written using the asterisk (`*`) to indicate
multiplication, the forward slash (`/`) to indicate division, and the
caret (`^`) to indicate exponentiation. Exponents can include negative
and decimal values. Standard one-letter metric prefixes can be applied
to any unit. Supported base units are:

-   Mass: `g`
-   Length: `m`, `micron`, `angstrom`, `Å`
-   Time: `s`, `min`, `hr`
-   Temperature: `K`, `C`
-   Current: `A`
-   Quantity: `mol` (gram mole), `gmol`, `mole`, `kmol`, `kgmol`,  `molec`

Supported compound units are:

-   Energy: `J`, `cal`, `erg`, `eV`
-   Activation Energy: `K`, any unit of energy per quantity (`J/kmol`,
    `cal/mol`, etc.), or any unit of energy (such as [eV])
-   Force: `N`, `dyn`
-   Pressure: `Pa`, `atm`, `bar`, `dyn/cm^2`
-   Volume: `m^3`, `liter`, `L`, `l`, `cc`
-   Other electrical units: `ohm`, `V`, `coulomb`

Units can be specified on individual input values by placing them after
the value, separated by a space:

```yaml
{A: 1.45e9 cm^3/kmol, b: 0.4, Ea: 21033 kJ/kmol}
```

or by using a `units` mapping:

```yaml
units: {mass: g, quantity: mol, pressure: atm, activation-energy: cal/mol}
```

A `units` mapping will set the default units for all values within the
same YAML list or mapping, including any nested lists and mappings.
Units not specified by a mapping use the values from higher level
mappings, or the %Cantera defaults if no `units` mapping specifies
applicable units. If a `units` mapping appears in a list, it must be the
first item in that list.

Default units may be set for `mass`, `length`, `time`, `temperature`,
`current`, `quantity`, `pressure`, `energy`, and `activation-energy`.
The units `pressure` and `energy` are used when these units appear
explicitly in the units that a value is being converted to within
%Cantera. For example, a conversion to `N/m^2` will use the default units
for mass, length, and time, while a conversion to `Pa` will use the
default units for pressure.

Activation energies given in temperature units will be implicitly
converted to energy per quantity by dividing by the gas constant.
Activation energies given in pure energy units such as eV will be
converted to energy per quantity by multiplying by the Avogadro
constant. Setting default units for `energy` and `quantity` will
determine the default units of `activation-energy`, which can be
overridden by explicitly giving the desired units of
`activation-energy`.

# Default Units {#sec-yaml-units-default}

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
