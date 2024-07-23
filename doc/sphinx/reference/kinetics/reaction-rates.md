# Reaction Rates

Here, we describe how Cantera calculates chemical reaction rates for various reaction
types.

(sec-elementary)=
## Elementary Reactions

The basic reaction type is a homogeneous reaction with a pressure-independent
rate coefficient and mass action kinetics. For example:

$$  a\t{A} + b\t{B} \rightleftharpoons c\t{C} + d\t{D}  $$

where A and B are reactant species, C and D are product species, and $a, b, c, $ and $d$
are stoichiometric coefficients.

The forward reaction rate is then calculated as:

$$  R_f = [\t{A}]^a [\t{B}]^b k_f  $$

where $k_f$ is the forward rate constant, calculated using one of the available rate
parameterizations such as the [modified Arrhenius](sec-arrhenius-rate) form.

```{admonition} YAML Usage
:class: tip
An elementary reaction with an Arrhenius reaction rate can be defined in the YAML format
using the [`elementary`](sec-yaml-elementary) reaction `type`, or by omitting the
reaction `type` entry, as it represents the default. An exception to this default is
when the same species occurs on both sides of the reaction equation, in which case the
reaction is treated as a
[three-body reaction for a specific collider](sec-three-body-specific-collider).
```

(sec-three-body-reaction)=
## Three-Body Reactions

A three-body reaction is a gas-phase reaction of the form:

$$  \t{A + B + M \rightleftharpoons AB + M}  $$

Here $\t{M}$ is an unspecified collision partner that carries away excess energy to
stabilize the $\t{AB}$ molecule (forward direction) or supplies energy to break the
$\t{AB}$ bond (reverse direction). In addition to the generic collision partner
$\t{M}$, it is also possible to explicitly specify a colliding species. In both
cases, the reaction type can be automatically inferred by Cantera and does not need to
be explicitly specified by the user.

Different species may be more or less effective in acting as the collision partner. A
species that is much lighter than $\t{A}$ and $\t{B}$ may not be able to
transfer much of its kinetic energy, and so would be inefficient as a collision partner.
On the other hand, a species with a transition from its ground state that is nearly
resonant with one in the $\t{AB^*}$ activated complex may be much more effective at
exchanging energy than would otherwise be expected.

These effects can be accounted for by defining a collision efficiency $\epsilon$ for
each species, defined such that the forward reaction rate is

$$  R_f = [\t{A}][\t{B}][\t{M}] k_f(T)  $$

where

$$  [\t{M}] = \sum_{k} \epsilon_k C_k  $$

where $C_k$ is the concentration of species $k$. Since any constant collision efficiency
can be absorbed into the rate coefficient $k_f(T)$, the default collision efficiency is
1.0.

:::{versionadded} 3.0
The rate coefficient $k_f(T)$ may be implemented using any rate parameterization
supported by Cantera, not just the modified Arrhenius form.
:::

(sec-three-body-specific-collider)=
### Collider-specific rate parameterizations

Sometimes, accounting for a particular third body's collision efficiency may require an
alternate set of rate parameters entirely. In this case, two reactions are written:

$$
\t{A + B + M \rightleftharpoons AB + M \quad (R1)}

\t{A + B + C \rightleftharpoons AB + C \quad (R2)}
$$

where the third-body efficiency for C in the first reaction should be explicitly set to
zero. For the second reaction, the efficiencies will automatically be set to one for C
and zero for all other colliders.

```{admonition} YAML Usage
:class: tip
A three-body reaction may be defined in the YAML format using the
[`three-body`](sec-yaml-three-body) reaction `type` or, if no `type` is specified,
identified automatically by the presence of the generic third body M or a specific
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


(sec-electrochemical-reactions)=
## Electrochemical Reactions

In an electrochemical reaction (one that moves electrical charge from one phase of matter to another), the electric potential difference $\Delta\phi$ at the phase boundary exerts an additional "force" on the reaction that must be accounted for in the rate expression.

The free energy of the reaction equals the electrochemical potential change, $\Delta\tilde{\mu}_{\rm rxn}$:

$$\Delta\tilde{\mu}_{\rm rxn} = \Delta\mu_{\rm rxn} + n_{\rm elec}F\Delta\phi$$

where $\mu_{\rm rxn}$ is the chemical potential, $n_{\rm elec}$ is the total electrical charge moved across the phase boundary (e.g. from "phase 1" to "phase 2", where $\Delta\phi = \phi_2 - \phi_1$), and $F$ is Faraday's constant (96,485 Coulombs per mole of charge).

Cantera's charge transfer treatment assumes a reversible reaction with a linear energy profile in the region of the transition state.  From above, for any $\Delta\phi$ the free energy of the reaction changes by $n_{\rm elec}F\Delta\phi$.  The transition state energy, meanwhile, changes by a fraction of this, $\beta n_{\rm elec}F\Delta\phi$, where the "symmetry parameter" $\beta$ is a number between 0 and 1.

This means that the activation energy for the reaction changes:
- The barrier height for the forward reaction increases by $\beta n_{\rm elec}F\Delta\phi$.
- The reverse reaction barrier height decreases by $\left(1-\beta\right) n_{\rm elec}F\Delta\phi$.

Note that $n_{\rm elec}$ and $\Delta \phi$ both have a sign, so the terms "increase" and "decrease" are relative; the forward barrier height might increase by a negative amount (i.e. decrease), for instance.

From transition state theory, the forward and reverse reaction rates are therefore calculated as:

$$ R_f = k_f\exp\left(-\frac{\beta n_{\rm elec}F\Delta\phi}{RT} \right)\Pi_k C_{{\rm ac},\,k}^{\nu_k^\prime}\$$

and

$$ R_r = k_r\exp\left(\frac{\left(1-\beta\right)n_{\rm elec}F\Delta\phi}{RT} \right)\Pi_k C_{{\rm ac},\,k}^{\nu_k^{\prime\prime}},$$

respectively, where $k_f$ and $k_r$ are the normal chemical rate coefficients in the absence of any electric potential difference (e.g., calculated using Arrhenius coefficients), $C_{{\rm ac},\,k}$ is the activity concentration of species $k$, $\nu_k^\prime$ and $\nu_k^{\prime\prime}$ are the forward and reverse stoichiometric coefficients, respectively, for species $k$ for this reaction, and $R$ and $T$ are the universal gas constant and temperature, respectively.

Note that Cantera's actual software implementation looks quite different from the description above, which is meant solely to give a clearer understanding of the science behind Cantera's calculations.
```{admonition} YAML Usage
:class: tip
- Electrochemical reactions only occur at phase boundaries and therefore use the standard [``interface``](sec-yaml-interface-Arrhenius) reaction rate implementation.
- Charge transfer is automatically detected, and $n_{\rm elec}$ automatically calculated.  If no value for `beta` is provided, an [``electrochemical``](sec-yaml-electrochemical-reaction) reaction assumes a default of ``beta = 0.5``.
```

(sec-butler-volmer)=
### The Butler-Volmer Form

Cantera's electrochemical reaction rate calculation is equivalent to the commonly-used Butler-Volmer rate form. In Butler-Volmer, the net rate of progress, $R_{\rm net} = R_f - R_r$, can be written as:

$$ R_{\rm net} = \frac{i_\circ}{n_{\rm elec}F}\left[\exp\left(-\frac{\beta n_{\rm elec}F\eta}{RT} \right) - \exp\left( \frac{\left(1-\beta\right)n_{\rm elec}F\eta}{RT}\right) \right]$$

where the kinetic rate constant $i_\circ$ is known as the "exchange current density" and $\eta$ the "overpotential" -- the difference between the actual electric potential difference and that which would set the reaction to equilibrium:

$$ \eta = \Delta\phi - \Delta\phi_{\rm equil}$$

To convert between the two forms, the exchange current density varies with the chemical state and can be calculated as:

$$i_\circ = n_{\rm elec}Fk_f^{\left(1-\beta\right)}k_r^\beta\Pi_k C_{{\rm ac},\,k}^{\left(1-\beta\right)\nu_k^\prime}\Pi_k C_{{\rm ac},\,k}^{\beta\nu_k^{\prime\prime}}.$$


````{admonition} YAML Usage
:class: tip
- One can explicitly provide an exchange current density, rather than the $k_f$ value, by setting the optional ``exchange-current-density-formulation`` field to ``true``.

```yaml
- equation: LiC6 <=> Li+(e) + C6
  rate-constant: [5.74, 0.0, 0.0]
  beta: 0.4
  exchange-current-density-formulation: true
```
Here, the rate constant Arrhenius parameters will be used to calculate the exchange current density.
````

(sec-reaction-orders)=
## Reaction Orders

Explicit reaction orders different from the stoichiometric coefficients are sometimes
used for non-elementary reactions. For example, consider the global reaction:

$$
\t{C_8H_{18} + 12.5 O_2 \rightarrow 8 CO_2 + 9 H_2O}
$$

the forward rate constant might be given as {cite:p}`westbrook1981`:

$$
k_f = 4.6 \times 10^{11} [\t{C_8H_{18}}]^{0.25} [\t{O_2}]^{1.5}
       \exp\left(\frac{30.0\,\t{kcal/mol}}{RT}\right)
$$

Special care is required in this case since the units of the pre-exponential factor
depend on the sum of the reaction orders, which may not be an integer.

Note that you can change reaction orders only for irreversible reactions.

Normally, reaction orders are required to be positive. However, in some cases negative
reaction orders are found to be better fits for experimental data. In these cases, the
default behavior may be overridden in the input file.

````{admonition} YAML Usage
:class: tip
To include explicit orders for the reaction above, it can be written in the YAML format as:

```yaml
- equation: C8H18 + 12.5 O2 => 8 CO2 + 9 H2O
  units: {length: cm, quantity: mol, activation-energy: kcal/mol}
  rate-constant: {A: 4.5e+11, b: 0.0, Ea: 30.0}
  orders: {C8H18: 0.25, O2: 1.5}
```
````
