# Reaction Rates

Here, we describe how Cantera calculates chemical reaction rates for various reaction
types.

(sec-elementary)=
## Elementary Reactions

The basic reaction type is a homogeneous reaction with a pressure-independent
rate coefficient and mass action kinetics. For example:

$$  \mathrm{A + B \rightleftharpoons C + D}  $$

The forward reaction rate is then calculated as:

$$  R_f = [\mathrm{A}] [\mathrm{B}] k_f  $$

where $k_f$ is the forward rate constant, calculated using one of the available rate
parameterizations such as [modified Arrhenius](sec-arrhenius-rate) form.

```{admonition} YAML Usage
:class: tip
An elementary reaction with an Arrhenius reaction rate can be defined in the YAML format
using the [`elementary`](sec-yaml-elementary) reaction `type`, or by omitting the
reaction `type` entry, as it represents the default. In case the `type` entry is omitted
and a species occurs on both sides of the reaction equation, Cantera infers that the
reaction type is [`three-body`](sec-yaml-three-body).
```

(sec-three-body-reaction)=
## Three-Body Reactions

A three-body reaction is a gas-phase reaction of the form:

$$  \mathrm{A + B + M \rightleftharpoons AB + M}  $$

Here $\mathrm{M}$ is an unspecified collision partner that carries away excess energy to
stabilize the $\mathrm{AB}$ molecule (forward direction) or supplies energy to break the
$\mathrm{AB}$ bond (reverse direction). In addition to the generic collision partner
$\mathrm{M}$, it is also possible to explicitly specify a colliding species. In this
case, the reaction type is automatically inferred by Cantera.

Different species may be more or less effective in acting as the collision partner. A
species that is much lighter than $\mathrm{A}$ and $\mathrm{B}$ may not be able to
transfer much of its kinetic energy, and so would be inefficient as a collision partner.
On the other hand, a species with a transition from its ground state that is nearly
resonant with one in the $\mathrm{AB^*}$ activated complex may be much more effective at
exchanging energy than would otherwise be expected.

These effects can be accounted for by defining a collision efficiency $\epsilon$ for
each species, defined such that the forward reaction rate is

$$  R_f = [\mathrm{A}][\mathrm{B}][\mathrm{M}]k_f(T)  $$

where

$$  [\mathrm{M}] = \sum_{k} \epsilon_k C_k  $$

where $C_k$ is the concentration of species $k$. Since any constant collision efficiency
can be absorbed into the rate coefficient $k_f(T)$, the default collision efficiency is
1.0.

:::{versionadded} 3.0
The rate coefficient $k_f(T)$ may be implemented using any rate parameterization
supported by Cantera, not just the modified Arrhenius form.
:::

Sometimes, accounting for a particular third body's collision efficiency may require an
alternate set of rate parameters entirely. In this case, two reactions are written:

$$
\mathrm{A + B + M \rightleftharpoons AB + M \quad (R1)}

\mathrm{A + B + C \rightleftharpoons AB + C \quad (R2)}
$$

where the third-body efficiency for C in the first reaction is set to zero.

```{admonition} YAML Usage
:class: tip
A three-body reaction may be defined in the YAML format using the
[`three-body`](sec-yaml-three-body) reaction `type`, or identified automatically if no
`type` is specified by the presence of the generic third body M or a specific
non-reactive species (for example, C in R2 above).
```

## Pressure-dependent Reactions

For pressure-dependent reactions where the behavior is more complex than described
by the three-body form, the pressure dependency is folded into the calculation of the
rate constant. Cantera supports several ways of representing pressure-dependent
reactions:

- [](sec-falloff-rate)
- [](sec-chemically-activated-rate)
- [](sec-plog-rate)
- [](sec-chebyshev-rate)

(sec-reaction-orders)=
## Reaction Orders

Explicit reaction orders different from the stoichiometric coefficients are sometimes
used for non-elementary reactions. For example, consider the global reaction:

$$
\mathrm{C_8H_{18} + 12.5 O_2 \rightarrow 8 CO_2 + 9 H_2O}
$$

the forward rate constant might be given as {cite:p}`westbrook1981`:

$$
k_f = 4.6 \times 10^{11} [\mathrm{C_8H_{18}}]^{0.25} [\mathrm{O_2}]^{1.5}
       \exp\left(\frac{30.0\,\mathrm{kcal/mol}}{RT}\right)
$$

Special care is required in this case since the units of the pre-exponential factor
depend on the sum of the reaction orders, which may not be an integer.

Note that you can change reaction orders only for irreversible reactions.

Normally, reaction orders are required to be positive. However, in some cases negative
reaction orders are found to be better fits for experimental data. In these cases, the
default behavior may be overridden in the input file.
