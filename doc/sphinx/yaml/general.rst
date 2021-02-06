.. highlight:: yaml

*****************
General Structure
*****************

Sections
--------

The top level of a Cantera `YAML <https://yaml.org/spec/1.2/spec.html#Introduction>`__
input file is a mapping that defines different input file sections. Each
section consists of a list of mappings that define objects of the same type,
e.g., reactions, species, phases, or elements. The ``phases`` section of an input
file contains all of the phase definitions. Multiple sections containing
reaction, species, or element definitions can be used. The specific names
``reactions``, ``species``, and ``elements`` are used as defaults when looking
for :ref:`sec-yaml-reactions`, :ref:`sec-yaml-species`, and
:ref:`sec-yaml-elements` to add to a phase. A simple input file has the
following structure::

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

Units
-----

While Cantera generally works internally in SI units, input values can be
provided using a number of different units.

Compound units are written using the asterisk (``*``) to indicate
multiplication, the forward slash (``/``) to indicate division, and the caret
(``^``) to indicate exponentiation. Exponents can include negative and decimal
values. Standard one-letter metric prefixes can be applied to any unit.
Supported base units are:

- Mass: ``g``
- Length: ``m``, ``micron``, ``angstrom``, ``Å``
- Time: ``s``, ``min``, ``hr``
- Temperature: ``K``, ``C``
- Current: ``A``
- Quantity: ``mol`` (gram mole), ``gmol``, ``mole``, ``kmol``, ``kgmol``, ``molec``

Supported compound units are:

- Energy: ``J``, ``cal``, ``erg``, ``eV``
- Activation Energy: ``K``, or any unit of energy per quantity (``J/kmol``,
  ``cal/mol``, etc.)
- Force: ``N``, ``dyn``
- Pressure: ``Pa``, ``atm``, ``bar``, ``dyn/cm^2``
- Volume: ``m^3``, ``liter``, ``L``, ``l``, ``cc``
- Other electrical units: ``ohm``, ``V``, ``coulomb``

Units can be specified on individual input values by placing them after the
value, separated by a space::

    {A: 1.45e9 cm^3/kmol, b: 0.4, Ea: 21033 kJ/kmol}

or by using a ``units`` mapping::

    units: {mass: g, quantity: mol, pressure: atm, activation-energy: cal/mol}

A ``units`` mapping will set the default units for all values within the same
YAML list or mapping, including any nested lists and mappings. Units not
specified by a mapping use the values from higher level mappings, or the Cantera
defaults if no ``units`` mapping specifies applicable units. If a ``units``
mapping appears in a list, it must be the first item in that list.

Default units may be set for ``mass``, ``length``, ``time``, ``temperature``,
``current``, ``quantity``, ``pressure``, ``energy``, and ``activation-energy``.
The units ``pressure`` and ``energy`` are used when these units appear
explicitly in the units that a value is being converted to within Cantera. For
example, a conversion to ``N/m^2`` will use the default units for mass, length,
and time, while a conversion to ``Pa`` will use the default units for pressure.

Conversions of activation energies implicitly include scaling by the gas
constant where necessary. Setting default units for ``energy`` and ``quantity``
will determine the default units of ``activation-energy``, which can be
overridden by explicitly giving the desired units of ``activation-energy``.
