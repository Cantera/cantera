# Reaction Entries {#sec-yaml-reactions}

[TOC]

## Common Fields

The fields common to all `reaction` entries are:

<b>`equation`</b>

The stoichiometric equation for the reaction. Each term (that is,
stoichiometric coefficient, species name, `+` or `<=>`) in the
equation must be separated by a space.

Reversible reactions may be written using `<=>` or `=` to separate
reactants and products. Irreversible reactions are written using `=>`.

<b>`type`</b>

A string specifying the type of reaction or rate coefficient
parameterization. The default is `elementary`. Reaction types are:

-   <tt>@ref sec-yaml-elementary</tt>
-   <tt>@ref sec-yaml-three-body</tt>
-   <tt>@ref sec-yaml-Blowers-Masel</tt>
-   <tt>@ref sec-yaml-two-temperature-plasma</tt>
-   <tt>@ref sec-yaml-falloff</tt>
-   <tt>@ref sec-yaml-chemically-activated</tt>
-   <tt>@ref sec-yaml-pressure-dependent-Arrhenius</tt>
-   <tt>@ref sec-yaml-Chebyshev</tt>

Reactions without a specified `type` on surfaces or edges are automatically treated as
<tt>@ref sec-yaml-interface-Arrhenius</tt> reactions, unless a `sticking-coefficient`
entry implies a <tt>@ref sec-yaml-sticking-Arrhenius</tt> reaction. Interface reactions
that involve charge transfer between phases are automatically treated as
<tt>@ref sec-yaml-electrochemical-reaction</tt> reactions.

Reactions on surfaces or edges specifying `type` as `Blowers-Masel`
are treated as <tt>@ref sec-yaml-interface-Blowers-Masel</tt> or
<tt>@ref sec-yaml-sticking-Blowers-Masel</tt>.

<b>`duplicate`</b>

Boolean indicating whether the reaction is a known duplicate of
another reaction. The default is `false`.

<b>`orders`</b>

An optional mapping of species to explicit reaction orders to use.
Reaction orders for reactant species not explicitly mentioned are
taken to be their respective stoichiometric coefficients. See
[Reaction orders](https://cantera.org/science/kinetics.html#reaction-orders)
for additional information.

<b>`negative-orders`</b>

Boolean indicating whether negative reaction orders are allowed. The
default is `false`.

<b>`nonreactant-orders`</b>

Boolean indicating whether orders for non-reactant species are
allowed. The default is `false`.

Depending on the reaction `type`, other fields may be necessary to
specify the rate of the reaction.

## Reaction Rate Expressions {#sec-yaml-reaction-rates}

### Arrhenius {#sec-yaml-Arrhenius-rate}

Arrhenius rate expressions are specified as a mapping with fields:

<b>`A`</b>

The pre-exponential factor \f$ A \f$

<b>`b`</b>

The temperature exponent \f$ b \f$

<b>`Ea`</b>

The activation energy \f$ E_a \f$

or a corresponding three-element list. The following are equivalent:

``` yaml
{A: -2.70000E+13 cm^3/mol/s, b: 0, Ea: 355 cal/mol}
[-2.70000E+13 cm^3/mol/s, 0, 355 cal/mol]
```

### Blowers-Masel {#sec-yaml-Blowers-Masel-rate}

Blowers-Masel rate expressions calculate the rate constant based on the
Blowers Masel approximation as [described
here](https://cantera.org/science/kinetics.html#sec-blowers-masel). The
rate parameters are specified as a mapping with fields:

<b>`A`</b>

The pre-exponential factor \f$ A\f$

<b>`b`</b>

The temperature exponent \f$ b\f$

<b>`Ea0`</b>

The intrinsic activation energy \f$ E_{a0}\f$

<b>`w`</b>

The average of the bond dissociation energy of the bond breaking and
that being formed in the reaction \f$ w\f$

or a corresponding four-element list. The following are equivalent:

``` yaml
{A: 3.87e+04 cm^3/mol/s, b: 2.7, Ea0: 6260.0 cal/mol, w: 1e9 cal/mol}
[3.87e+04 cm^3/mol/s, 2.7, 6260.0 cal/mol, 1e9 cal/mol]
```

### Two-Temperature Plasma {#sec-yaml-two-temperature-plasma-rate}

Two-temperature plasma reactions involve an electron as one of the
reactants, where the electron temperature may differ from the gas
temperature as [described
here](https://cantera.org/science/kinetics.html#two-temperature-plasma-reactions).
The rate parameters are specified as a mapping with fields:

<b>`A`</b>

The pre-exponential factor

<b>`b`</b>

The temperature exponent, which is applied to the electron temperature

<b>`Ea-gas`</b>

The activation energy term \f$ E_{a,g}\f$  that is related to the gas temperature

<b>`Ea-electron`</b>

The activation energy term \f$ E_{a,e}\f$  that is related to the electron temperature

or a corresponding four-element list. The following are equivalent:

``` yaml
{A: 17283, b: -3.1, Ea-gas: -5820 J/mol, Ea-electron: 1081 J/mol}
[17283, -3.1, -5820 J/mol, 1081 J/mol]
```

## Efficiencies {#sec-yaml-efficiencies}

Some reaction types include parameters for the \"efficiency\" of
different species as third-body colliders. For these reactions, the
following additional fields are supported:

<b>`efficiencies`</b>

A mapping of species names to efficiency values

<b>`default-efficiency`</b>

The efficiency for use for species not included in the
`efficiencies` mapping. Defaults to 1.0.
