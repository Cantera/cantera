# Reaction Rate Blocks {#sec-yaml-reaction-rate-blocks}

[TOC]

The following reaction rate expressions are used as building blocks for most
@ref sec-yaml-reaction-models.

## Arrhenius {#sec-yaml-Arrhenius-rate}

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

## Blowers-Masel {#sec-yaml-Blowers-Masel-rate}

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

## Two-Temperature Plasma {#sec-yaml-two-temperature-plasma-rate}

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
