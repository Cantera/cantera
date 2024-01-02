# Rate Constant Parameterizations

This page describes the different parameterizations available in Cantera for
calculating the forward rate constant $k_f$ for a reaction.

(sec-arrhenius-rate)=
## Arrhenius Rate Expressions

An Arrhenius rate is described by the
[modified Arrhenius function](https://en.wikipedia.org/wiki/Arrhenius_equation#Modified_Arrhenius_equation):

$$  k_f = A T^b e^{-E_a / RT}  $$

where $A$ is the pre-exponential factor, $T$ is the temperature, $b$ is the temperature
exponent, $E_a$ is the activation energy, and $R$ is the gas constant.

:::{admonition} YAML Usage
:class: tip
An Arrhenius rate can be specified for a reaction in the YAML format by providing an
[Arrhenius](sec-yaml-Arrhenius-rate) rate expression for the reaction's `rate-constant`
field.
:::

(sec-falloff-rate)=
## Falloff Reactions

A falloff reaction is one that has a rate that is first-order in the total concentration
of third-body colliders $\def\MM{[\mathrm{M}]} \MM$ at low pressure, like a
[three-body reaction](sec-three-body-reaction), but becomes zero-order in $\MM$ as $\MM$
increases. Dissociation/association reactions of polyatomic molecules often exhibit this
behavior.

The simplest expression for the rate coefficient for a falloff reaction is the Lindemann
form {cite:p}`lindemann1922`:

$$  k_f(T, \MM) = \frac{k_0 \MM}{1 + \frac{k_0 \MM}{k_\infty}}  $$

In the low-pressure limit, this approaches $k_0 \MM$, and in the high-pressure limit it
approaches $k_\infty$.

Defining the non-dimensional reduced pressure:

$$  P_r = \frac{k_0 \MM}{k_\infty}  $$

The rate constant may be written as

$$  k_f(T, P_r) = k_\infty \left(\frac{P_r}{1 + P_r}\right)  $$

More accurate models for unimolecular processes lead to other, more complex,
forms for the dependence on reduced pressure. These can be accounted for by
multiplying the Lindemann expression by a function $F(T, P_r)$:

$$  k_f(T, P_r) = k_\infty \left(\frac{P_r}{1 + P_r}\right) F(T, P_r)  $$

This expression is used to compute the rate coefficient for falloff reactions. The
function $F(T, P_r)$ is the falloff function.

```{admonition} YAML Usage
:class: tip
A falloff reaction may be defined in the YAML format using the
[`falloff`](sec-yaml-falloff) reaction `type`.
```

(sec-troe-falloff)=
### The Troe Falloff Function

A widely-used falloff function is the one proposed by {cite:t}`gilbert1983`:

\begin{gather*}
\log_{10} F(T, P_r) = \frac{\log_{10} F_{cent}(T)}{1 + f_1^2} \\

F_{cent}(T) = (1-A) \exp(-T/T_3) + A \exp (-T/T_1) + \exp(-T_2/T) \\

f_1 = (\log_{10} P_r + C) / (N - 0.14 (\log_{10} P_r + C)) \\

C = -0.4 - 0.67\; \log_{10} F_{cent} \\

N = 0.75 - 1.27\; \log_{10} F_{cent}
\end{gather*}

```{admonition} YAML Usage
:class: tip
A Troe falloff function may be specified in the YAML format using the
[`Troe`](sec-yaml-falloff) field in the reaction entry. The first three parameters, $(A,
T_3, T_1)$, are required. The fourth parameter, $T_2$, is optional; if omitted, the last
term of the falloff function is not used.
```

(sec-tsang-falloff)=
### Tsang's Approximation to $F_{cent}$

Wing Tsang presented approximations for the value of $F_{cent}$ for Troe falloff in
databases of reactions, for example, {cite:t}`tsang1991`. Tsang's approximations are
linear in temperature:

$$  F_{cent} = A + BT  $$

where $A$ and $B$ are constants. The remaining equations for $C$, $N$, $f_1$, and $F$
from the [Troe](sec-troe-falloff) falloff function are not affected.

```{admonition} YAML Usage
:class: tip
A Tsang falloff function may be specified in the YAML format using the
[`Tsang`](sec-yaml-falloff) field in the reaction entry.
```

```{versionadded} 2.6
```

(sec-sri-falloff)=
### The SRI Falloff Function

This falloff function is based on the one originally due to {cite:t}`stewart1989`, which
required three parameters $a$, $b$, and $c$. {cite:t}`kee1989` generalized this function
slightly by adding two more parameters $d$ and $e$. The original form corresponds to
$d = 1$ and $e = 0$. Cantera supports the extended 5-parameter form, given by:

$$  F(T, P_r) = d \bigl[a \exp(-b/T) + \exp(-T/c)\bigr]^{1/(1+\log_{10}^2 P_r )} T^e  $$

In keeping with the nomenclature of {cite:t}`kee1989`, we will refer to this as the
**SRI falloff function**.

```{admonition} YAML Usage
:class: tip
An SRI falloff function may be specified in the YAML format using the
[`SRI`](sec-yaml-falloff) field in the entry.
```

(sec-chemically-activated-rate)=
## Chemically-Activated Reactions

For these reactions, the rate falls off as the pressure increases, due to collisional
stabilization of a reaction intermediate. Example:

$$  \mathrm{Si + SiH_4 (+M) \leftrightarrow Si_2H_2 + H_2 (+M)}  $$

which competes with:

$$  \mathrm{Si + SiH_4 (+M) \leftrightarrow Si_2H_4 (+M)}  $$

Like falloff reactions, chemically-activated reactions are described by blending between
a low-pressure and a high-pressure rate expression. The difference is that the forward
rate constant is written as being proportional to the *low-pressure* rate constant:

$$  k_f(T, P_r) = k_0 \left(\frac{1}{1 + P_r}\right) F(T, P_r)  $$

and the optional blending function $F$ may described by any of the parameterizations
allowed for falloff reactions.

```{admonition} YAML Usage
:class: tip
Chemically-activated reactions can be defined in the YAML format using the
[`chemically-activated`](sec-yaml-chemically-activated) reaction `type`.
```

(sec-plog-rate)=
## Pressure-Dependent Arrhenius Rate Expressions (P-Log)

This parameterization represents pressure-dependent reaction rates by logarithmically
interpolating between Arrhenius rate expressions at various pressures. Given two rate
expressions at two specific pressures:

$$
P_1: k_1(T) = A_1 T^{b_1} e^{-E_1 / RT}

P_2: k_2(T) = A_2 T^{b_2} e^{-E_2 / RT}
$$

The rate at an intermediate pressure $P_1 < P < P_2$ is computed as

$$
\log k(T,P) = \log k_1(T) + \bigl(\log k_2(T) - \log k_1(T)\bigr)
    \frac{\log P - \log P_1}{\log P_2 - \log P_1}
$$

Multiple rate expressions may be given at the same pressure, in which case the rate used
in the interpolation formula is the sum of all the rates given at that pressure. For
pressures outside the given range, the rate expression at the nearest pressure is used.

```{caution}
Negative A-factors can be used for any of the rate expressions at a given pressure.
However, the sum of all of the rates at a given pressure **must** be positive, due to
the logarithmic interpolation of the rate for intermediate pressures. When a P-log type
reaction is initialized, Cantera does a validation check for a range of temperatures
that the sum of the reaction rates at each pressure is positive. Unfortunately, if these
checks fail, the only options are to remove the reaction or contact the author of the
reaction/mechanism in question, because the reaction is mathematically unsound.
```

```{admonition} YAML Usage
:class: tip
P-log reactions can be defined in the YAML format using the
[`pressure-dependent-Arrhenius`](sec-yaml-pressure-dependent-Arrhenius) reaction `type`.
```

(sec-chebyshev-rate)=
## Chebyshev Reaction Rate Expressions

Chebyshev rate expressions represent a phenomenological rate coefficient $k(T,P)$ in
terms of a bivariate Chebyshev polynomial. The rate constant can be written as:

$$
\log k(T,P) = \sum_{t=1}^{N_T} \sum_{p=1}^{N_P} \alpha_{tp}
                         \phi_t(\tilde{T}) \phi_p(\tilde{P})
$$

where $N_T$ is the order of the polynomial in the temperature dimension, $N_P$ is the
order of the polynomial in the pressure dimension, $\alpha_{tp}$ are the constants
defining the rate, $\phi_n(x)$ is the Chebyshev polynomial of the first kind of degree
$n$ evaluated at $x$, and

$$
\tilde{T} \equiv \frac{2T^{-1} - T_\mathrm{min}^{-1} - T_\mathrm{max}^{-1}}
                       {T_\mathrm{max}^{-1} - T_\mathrm{min}^{-1}}

\tilde{P} \equiv \frac{2 \log P - \log P_\mathrm{min} - \log P_\mathrm{max}}
                       {\log P_\mathrm{max} - \log P_\mathrm{min}}
$$

are reduced temperatures and reduced pressures which map the ranges $(T_\mathrm{min},
T_\mathrm{max})$ and $(P_\mathrm{min}, P_\mathrm{max})$ to $(-1, 1)$.

A Chebyshev rate expression is specified in terms of the coefficient matrix $\alpha$ and
the temperature and pressure ranges.

```{caution}
The Chebyshev polynomials are not defined outside the interval $(-1,1)$, and therefore
extrapolation of rates outside the range of temperatures and pressure for which they are
defined is strongly discouraged.
```

```{admonition} YAML Usage
:class: tip
Chebyshev reactions can be defined in the YAML format using the
[`Chebyshev`](sec-yaml-Chebyshev) reaction `type`.
```

(sec-blowers-masel)=
## Blowers-Masel Reactions

In some circumstances like thermodynamic sensitivity analysis, or modeling heterogeneous
reactions from one catalyst surface to another, the enthalpy change of a reaction
($\Delta H$) can be modified. Due to the change in $\Delta H$, the activation energy of
the reaction must be adjusted accordingly to provide accurate simulation results. To
adjust the activation energy due to changes in the reaction enthalpy, the Blowers-Masel
rate expression is available. This approximation was proposed by {cite:t}`blowers2000`
to automatically scale activation energy as the reaction enthalpy is changed. The
activation energy estimation can be written as:

$$
E_a = \begin{cases}
   0 & \text{if } \Delta H \leq -4 E_a^0 \\
   \Delta H & \text{if } \Delta H \geq 4 E_a^0 \\
   \frac{\left( w + \frac{\Delta H }{2} \right)  (V_P - 2 w + \Delta H) ^2}
            {V_P^2 - 4 w^2 + \Delta H^2} & \text{otherwise}
   \end{cases}
$$

where

$$  V_P = 2 w \frac{w + E_a^0}{w - E_a^0},  $$

$w$ is the average of the bond dissociation energy of the bond breaking and that being
formed, $E_a^0$ is the intrinsic activation energy, and $\Delta H$ is the enthalpy
change of the reaction. Note that the expression is insensitive to $w$ as long as $w \ge
2 E_a^0$, so we can use an arbitrarily high value of $w = 1000\text{ kJ/mol}$.

After $E_a$ is evaluated, the reaction rate can be calculated using the modified
Arrhenius expression

$$  k_f = A T^b e^{-E_a / RT}.  $$

```{versionadded} 2.6
```

```{admonition} YAML Usage
:class: tip
Blowers Masel reactions can be defined in the YAML format using the
[Blowers-Masel](sec-yaml-Blowers-Masel-rate) reaction `type`.
```

(sec-surface-rate)=
## Surface Reactions

Heterogeneous reactions on surfaces are represented by an extended Arrhenius- like rate
expression, which combines the modified Arrhenius rate expression with further
corrections dependent on the fractional surface coverages $\theta_{k}$ of one or more
surface species. The forward rate constant for a reaction of this type is:

$$
k_f = A T^b \exp \left( - \frac{E_a}{RT} \right)
   \prod_k 10^{a_k \theta_k}
   \theta_k^{m_k}
   \exp \left( \frac{- E_k \theta_k}{RT} \right)
$$

where $A$, $b$, and $E_a$ are the modified Arrhenius parameters and $a_k$, $m_k$, and
$E_k$ are the coverage dependencies from species $k$.

::::{admonition} YAML Usage
:class: tip
In the YAML format, surface reactions are identified by the presence of surface species
and support several [additional options](sec-yaml-interface-Arrhenius).

The surface reaction `type` defaults to `interface-Arrhenius`, where the rate
expression uses the [`Arrhenius`](sec-elementary) parameterization (see [YAML
documentation](sec-yaml-interface-Arrhenius)).

:::{versionadded} 2.6
As an alternative, Cantera also supports the `interface-Blowers-Masel` surface reaction
`type`, which uses the [`Blowers-Masel`](sec-Blowers-Masel) parameterization (see [YAML
documentation](sec-yaml-interface-Blowers-Masel)).
:::
::::

(sec-sticking-rate)=
## Sticking Reactions

Sticking reactions represent a special case of surface reactions, where collisions
between gas-phase molecules and surfaces result in the gas-phase molecule sticking to
the surface. This process can be described as a reaction which is parameterized by a
sticking coefficient:

$$  \gamma = a T^b e^{-c/RT}  $$

where $a$, $b$, and $c$ are constants specific to the reaction. The values of these
constants must be specified so that the sticking coefficient $\gamma$ is between 0 and 1
for all temperatures.

The sticking coefficient is related to the forward rate constant by the formula:

$$  k_f = \frac{\gamma}{\Gamma_\mathrm{tot}^m} \sqrt{\frac{RT}{2 \pi W}}  $$

where $\Gamma_\mathrm{tot}$ is the total molar site density, $m$ is the sum of all the
surface reactant stoichiometric coefficients, and $W$ is the molecular weight of the gas
phase species.

::::{admonition} YAML Usage
:class: tip
Sticking reactions can be defined in the YAML format by specifying the rate constant
in the reaction's
[sticking-coefficient](sec-yaml-sticking-Arrhenius) field.

The sticking reaction `type` defaults to `sticking-Arrhenius`, where the rate expression
uses the [`Arrhenius`](sec-elementary) parameterization (see [YAML
documentation](sec-yaml-sticking-Arrhenius)).

:::{versionadded} 2.6
As an alternative, Cantera also supports the `sticking-Blowers-Masel` surface reaction
`type`, which uses the [`Blowers-Masel`](sec-Blowers-Masel) parameterization (see
[YAML documentation](sec-yaml-sticking-Blowers-Masel)).
:::
::::

(sec-two-temperature-plasma-rate)=
## Two-Temperature-Plasma Reactions

The two-temperature-plasma reaction is commonly used for non-equilibrium plasmas. The
reaction rate of a two-temperature-plasma reaction depends on both gas and electron
temperature {cite:p}`kossyi1992`, and can be expressed as:

$$
k_f = A T_e^b \exp \left( - \frac{E_{a,g}}{RT} \right)
   \exp \left(\frac{E_{a,e}(T_e - T)}{R T T_e}\right),
$$

where $A$ is the pre-exponential factor, $T$ is the temperature, $T_e$ is the electron
temperature, $b$ is the electron temperature exponent, $E_{a,g}$ is the activation
energy for gas, $E_{a,e}$ is the activation energy for electron and $R$ is the gas
constant.

```{versionadded} 2.6
```

:::{admonition} YAML Usage
:class: tip

Two-temperature plasma reactions can be defined in the YAML format by specifying
[`two-temperature-plasma`](sec-yaml-two-temperature-plasma) as the reaction `type` and
providing the two activation energies as part of the `rate-constant`.
:::
