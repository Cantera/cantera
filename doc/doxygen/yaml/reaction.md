# Reaction Entries {#sec-yaml-reactions}

[TOC]

## Overview

%Cantera supports a number of different types of reactions, including
several types of homogeneous reactions, surface reactions, and
electrochemical reactions. The reaction entries for all reaction types
some common features. These general fields of a reaction entry are
described first, followed by fields used for specific reaction types.

### Reaction Equation {#sec-yaml-reactions-equation}

The reaction equation, specified in the `equation` field of the reaction
entry, determines the reactant and product stoichiometry. All tokens
(species names, stoichiometric coefficients, `+`, and `<=>`) in the
reaction equation must be separated with spaces. Some examples of
correctly and incorrectly formatted reaction equations are shown below:

``` yaml
- equation: 2 CH2 <=> CH + CH3  # OK
- equation: 2 CH2<=>CH + CH3  # error - spaces required around '<=>''
- equation: 2CH2 <=> CH + CH3  # error - space required between '2' and 'CH2'
- equation: CH2 + CH2 <=> CH + CH3  # OK
- equation: 2 CH2 <=> CH+CH3  # error - spaces required around '+'
```

Whether the reaction is reversible or not is determined by the form of
the equality sign in the reaction equation. If either `<=>` or `=` is
found, then the reaction is regarded as reversible, and the reverse rate
will be computed based on the equilibrium constant. If, on the other
hand, `=>` is found, the reaction will be treated as irreversible.

### Reaction Type {#sec-yaml-reactions-type}

The type of the rate coefficient parameterization may be specified in
the `type` field of the `reaction` entry. Available reaction types are:

-   <tt>@ref sec-yaml-elementary</tt>:
    A reaction with a rate constant parameterized by a modified Arrhenius expression
-   <tt>@ref sec-yaml-three-body</tt>:
    A reaction involving a third-body collision
-   <tt>@ref sec-yaml-Blowers-Masel</tt>:
    A reaction rate constant parameterized as a modified Arrhenius
    reaction with one additional bond energy parameter to scale the
    activation energy according to the enthalpy of the reaction
-   <tt>@ref sec-yaml-two-temperature-plasma</tt>
-   <tt>@ref sec-yaml-falloff</tt>:
    A pressure-dependent reaction where the rate depends on the third-body
    concentration at low pressure but not at high pressure
-   <tt>@ref sec-yaml-chemically-activated</tt>:
    A pressure-dependent reaction where the rate depends on the third-body
    concentration at high pressure but not at low pressure
-   <tt>@ref sec-yaml-pressure-dependent-Arrhenius</tt>:
    A reaction rate parameterized by logarithmically interpolating between
    modified Arrhenius expressions at different pressures
-   <tt>@ref sec-yaml-Chebyshev</tt>:
    A reaction rate parameterized by a bivariate Chebyshev polynomial in
    pressure and temperature

Additional parameters defining the rate model for each of these
reaction types are described in the documentation linked above
(see also @subpage sec-yaml-reaction-models).

The default parameterization is `elementary`. Reactions without a specified
`type` on surfaces or edges are automatically treated as
<tt>@ref sec-yaml-interface-Arrhenius</tt> reactions, unless a `sticking-coefficient`
entry implies a <tt>@ref sec-yaml-sticking-Arrhenius</tt> reaction. Interface reactions
that involve charge transfer between phases are automatically treated as
<tt>@ref sec-yaml-electrochemical-reaction</tt> reactions.

Reactions on surfaces or edges specifying `type` as `Blowers-Masel`
are treated as <tt>@ref sec-yaml-interface-Blowers-Masel</tt> or
<tt>@ref sec-yaml-sticking-Blowers-Masel</tt>.

### Rate Parameterization

The following reaction rate building blocks are supported:

-   <tt>@ref sec-yaml-Arrhenius-rate</tt>
-   <tt>@ref sec-yaml-Blowers-Masel-rate</tt>
-   <tt>@ref sec-yaml-two-temperature-plasma-rate</tt>

Additional parameters defining these building blocks are described in the
documentation linked above (see also @subpage sec-yaml-reaction-rate-blocks).

Most reaction types in %Cantera are parameterized by one or more modified
<tt>@ref sec-yaml-Arrhenius-rate</tt> expressions, such as

@f[
    A T^b e^{-E_a / RT}
@f]

where @f$A@f$ is the pre-exponential factor, @f$T@f$ is the temperature, @f$b@f$ is
the temperature exponent, @f$E_a@f$ is the activation energy, and @f$R@f$ is the
gas constant. Rates in this form can be written as YAML mappings, where the key `E`
is used to specify @f$E_a@f$. For example:

``` yaml
{A: 1.0e13, b: 0, E: 7.3 kcal/mol}
```

The units of @f$A@f$ can be specified explicitly if desired. If not
specified, they will be determined based on the `quantity`, `length`,
and `time` units specified in the governing `units` fields. Since the
units of @f$A@f$ depend on the reaction order, the units of each reactant
concentration (dependent on phase type and dimensionality), and the
units of the rate of progress (different for homogeneous and
heterogeneous reactions), it is usually best not to specify units for
@f$A@f$, in which case they will be computed taking all of these factors
into account.

*Note:* if @f$b \ne 0@f$, then the term @f$T^b@f$ should have units of
@f$\mathrm{K}^b@f$, which would change the units of @f$A@f$. This is not done,
however, so the units associated with @f$A@f$ are really the units for
@f$k_f@f$. One way to formally express this is to replace @f$T^b@f$ by the
non-dimensional quantity @f$[T/(1\;\mathrm{K})]^b@f$.

### Third-Body Colliders {#sec-yaml-reactions-efficiencies}

Some reaction types include parameters for the \"efficiency\" of
different species as third-body colliders `M`.

**Example:**

``` yaml
- equation: O + H + M <=> OH + M
  type: three-body
  rate-constant: {A: 5.0e+5, b: -1.0, Ea: 0.0}
  default-efficiency: 0.
  efficiencies: {H2: 2.0, H2O: 6.0, AR: 0.7}
```

### Duplicate Reactions {#sec-yaml-reaction-duplicates}

When a reaction is imported into a phase, it is checked to see that it
is not a duplicate of another reaction already present in the phase, and
normally an error results if a duplicate is found. But in some cases, it
may be appropriate to include duplicate reactions, for example if a
reaction can proceed through two distinctly different pathways, each
with its own rate expression. Another case where duplicate reactions can
be used is if it is desired to implement a reaction rate coefficient of
the form:

@f[
    k_f(T) = \sum_{n=1}^{N} A_n T^{b_n} \exp(-E_n/RT)
@f]

While %Cantera does not provide such a form for reaction rates, it can be
implemented by defining @f$N@f$ duplicate reactions, and assigning one rate
coefficient in the sum to each reaction. By adding the field:

``` yaml
duplicate: true
```

to a reaction entry, then the reaction not only *may* have a duplicate,
it *must*. Any reaction that specifies that it is a duplicate, but
cannot be paired with another reaction in the phase that qualifies as
its duplicate generates an error.

### Negative Pre-exponential Factors

If some of the terms in the above sum have negative @f$A_n@f$, this scheme
fails, since %Cantera normally does not allow negative pre-exponential
factors. But if there are duplicate reactions such that the total rate
is positive, then the fact that negative @f$A@f$ parameters are acceptable
can be indicated by adding the field:

``` yaml
negative-A: true
```

### Reaction Orders {#sec-yaml-reaction-orders}

Explicit reaction orders different from the stoichiometric coefficients
are sometimes used for non-elementary reactions. For example, consider
the global reaction:

@f[
    \mathrm{C_8H_{18} + 12.5 O_2 \rightarrow 8 CO_2 + 9 H_2O}
@f]

the forward rate constant might be given as[^1]:

@f[
    k_f = 4.6 \times 10^{11} [\mathrm{C_8H_{18}}]^{0.25} [\mathrm{O_2}]^{1.5}
        \exp\left(\frac{30.0\,\mathrm{kcal/mol}}{RT}\right)
@f]

This reaction could be defined as:

``` yaml
- equation: C8H18 + 12.5 O2 => 8 CO2 + 9 H2O
  rate-constant: {A: 4.6e11, b: 0.0, Ea: 30.0 kcal/mol}
  orders: {C8H18: 0.25, O2: 1.5}
```

Special care is required in this case since the units of the
pre-exponential factor depend on the sum of the reaction orders, which
may not be an integer.

Note that you can change reaction orders only for irreversible
reactions.

**Negative Reaction Orders**

Normally, reaction orders are required to be positive. However, in some
cases negative reaction orders provide better fits for experimental
data. In these cases, the default behavior may be overridden by adding
the `negative-orders` field to the reaction entry. For example:

``` yaml
- equation: C8H18 + 12.5 O2 => 8 CO2 + 9 H2O
  rate-constant: {A: 4.6e11, b: 0.0, Ea: 30.0 kcal/mol}
  orders: {C8H18: -0.25, O2: 1.75}
  negative-orders: true
```

**Non-Reactant Orders**

Some global reactions could have reactions orders for non-reactant
species. In this case, the `nonreactant-orders` field must be added to
the reaction entry:

``` yaml
- equation: C8H18 + 12.5 O2 => 8 CO2 + 9 H2O
  rate-constant: {A: 4.6e11, b: 0.0, Ea: 30.0 kcal/mol}
  orders: {C8H18: -0.25, CO: 0.15}
  negative-orders: true
  nonreactant-orders: true
```

## Reaction API Reference

The fields common to most `reaction` entries are:

<b>`equation`</b>

The stoichiometric equation for the reaction. Each term (that is,
stoichiometric coefficient, species name, `+` or `<=>`) in the
equation must be separated by a space.

Reversible reactions may be written using `<=>` or `=` to separate
reactants and products. Irreversible reactions are written using `=>`.
See also @ref sec-yaml-reactions-equation.

<b>`type`</b>

A string specifying the type of reaction or rate coefficient
parameterization, see @ref sec-yaml-reactions-type.

<b>`duplicate`</b>

Boolean indicating whether the reaction is a known duplicate of
another reaction. The default is `false`.
See @ref sec-yaml-reaction-duplicates.

<b>`orders`</b>

An optional mapping of species to explicit reaction orders to use.
Reaction orders for reactant species not explicitly mentioned are
taken to be their respective stoichiometric coefficients. See
@ref sec-yaml-reaction-orders for additional information.

<b>`negative-orders`</b>

Boolean indicating whether negative reaction orders are allowed. The
default is `false`.

<b>`nonreactant-orders`</b>

Boolean indicating whether orders for non-reactant species are
allowed. The default is `false`.

<b>`efficiencies`</b>

Some reaction types include parameters for the \"efficiency\" of
different species as third-body colliders. The `efficiencies` field
holds a mapping of species names to efficiency values. See
@ref sec-yaml-reactions-efficiencies.

<b>`default-efficiency`</b>

The efficiency of third-body collider species not included in the
`efficiencies` mapping. Defaults to 1.0.

Depending on the reaction `type`, other fields may be necessary to
specify the rate of the reaction.
