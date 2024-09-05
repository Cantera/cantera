# Rate Constant Parameterizations

This page describes the different parameterizations available in Cantera for
calculating the forward rate constant $k_f$ for a reaction.

(sec-arrhenius-rate)=
## Arrhenius Rate Expressions

An Arrhenius rate is described by the
[modified Arrhenius function](https://en.wikipedia.org/wiki/Arrhenius_equation#Modified_Arrhenius_equation):

$$  k_f = A T^b e^{-E_a / RT}  $$

where $A$ is the pre-exponential factor, $T$ is the temperature, $b$ is the temperature
exponent, $E_a$ is the activation energy, and $R$ is the universal gas constant.

:::{admonition} YAML Usage
:class: tip
An Arrhenius rate can be specified for a reaction in the YAML format by providing an
[Arrhenius](sec-yaml-Arrhenius-rate) rate expression for the reaction's `rate-constant`
field.
:::

(sec-falloff-rate)=
## Falloff Reactions

A falloff reaction is one that has a rate that is first-order in the total concentration
of third-body colliders $[\t{M}]$ at low pressure, like a
[three-body reaction](sec-three-body-reaction), but becomes zero-order in $[\t{M}]$ as
$[\t{M}]$ increases. Dissociation/association reactions of polyatomic molecules often
exhibit this behavior.

The simplest expression for the rate coefficient for a falloff reaction is the Lindemann
form {cite:p}`lindemann1922`:

$$  k_f(T, [\t{M}]) = \frac{k_0 [\t{M}]}{1 + \frac{k_0 [\t{M}]}{k_\infty}}  $$

In the low-pressure limit, this approaches $k_0 [\t{M}]$, and in the high-pressure limit
it approaches $k_\infty$.

Defining the non-dimensional reduced pressure:

$$  P_r = \frac{k_0 [\t{M}]}{k_\infty}  $$

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
\log_{10} F(T, P_r) = \frac{\log_{10} F_\t{cent}(T)}{1 + f_1^2} \\

F_\t{cent}(T) = (1-A) \exp(-T/T_3) + A \exp (-T/T_1) + \exp(-T_2/T) \\

f_1 = (\log_{10} P_r + C) / (N - 0.14 (\log_{10} P_r + C)) \\

C = -0.4 - 0.67\; \log_{10} F_\t{cent} \\

N = 0.75 - 1.27\; \log_{10} F_\t{cent}
\end{gather*}

```{admonition} YAML Usage
:class: tip
A Troe falloff function may be specified in the YAML format using the
[`Troe`](sec-yaml-falloff) field in the reaction entry. The first three parameters, $(A,
T_3, T_1)$, are required. The fourth parameter, $T_2$, is optional; if omitted, the last
term of the falloff function is not used.
```

(sec-tsang-falloff)=
### Tsang's Approximation to $F_\t{cent}$

Wing Tsang presented approximations for the value of $F_\t{cent}$ for Troe falloff in
databases of reactions, for example, {cite:t}`tsang1991`. Tsang's approximations are
linear in temperature:

$$  F_\t{cent} = A + BT  $$

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

$$  \t{Si + SiH_4 (+M) \leftrightarrow Si_2H_2 + H_2 (+M)}  $$

which competes with:

$$  \t{Si + SiH_4 (+M) \leftrightarrow Si_2H_4 (+M)}  $$

Like falloff reactions, chemically-activated reactions are described by blending between
a low-pressure and a high-pressure rate expression. The difference is that the forward
rate constant is written as proportional to the *low-pressure* rate constant:

$$  k_f(T, P_r) = k_0 \left(\frac{1}{1 + P_r}\right) F(T, P_r)  $$

and the optional blending function $F$ may be described by any of the parameterizations
allowed for falloff reactions.

```{admonition} YAML Usage
:class: tip
Chemically-activated reactions can be defined in the YAML format using the
[`chemically-activated`](sec-yaml-chemically-activated) reaction `type`.
```

(sec-plog-rate)=
## Pressure-Dependent Arrhenius Rate Expressions (P-Log)

This parameterization represents pressure-dependent reaction rates by logarithmically
interpolating between Arrhenius rate expressions at various pressures {cite:p}`gou2011`.
Given two rate expressions at two specific pressures:

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
\tilde{T} \equiv \frac{2T^{-1} - T_\t{min}^{-1} - T_\t{max}^{-1}}
                       {T_\t{max}^{-1} - T_\t{min}^{-1}}

\tilde{P} \equiv \frac{2 \log P - \log P_\t{min} - \log P_\t{max}}
                       {\log P_\t{max} - \log P_\t{min}}
$$

are reduced temperatures and reduced pressures which map the ranges $(T_\t{min},
T_\t{max})$ and $(P_\t{min}, P_\t{max})$ to $(-1, 1)$.

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

(sec-linear-burke)=
## Linear Burke Rate Expressions

Linear Burke rate expressions employ the reduced-pressure linear mixture rule (LMR-R). This mixture rule is used to evaluate the rate constants of complex-forming reactions, and is a mole-fraction-weighted sum of the bath gas rate constants (when pure) evaluated at the reduced pressure ($R$) and temperature ($T$) of the mixture.

$$
k_{\text{LMR-R}}(T,P,\textit{\textbf{X}}) = \sum_{i} k_{i}(T,R_{\text{LMR}})\tilde{X}_{i,\text{LMR}}
$$

where the reduced pressure, $R$, in its most general form

$$
R_{\text{LMR}}(T,P,\textit{\textbf{X}}) = \frac{\sum_{i} \Lambda_{0,i}(T)X_i[M]}{\Lambda_{\infty}(T)}
$$

and the fractional contribution of each component to the reduced pressure, $\tilde{X}_{i}$

$$
\tilde{X}_{i,\text{LMR}}(T,P,\textit{\textbf{X}})=\frac{\Lambda_{0,i}(T)X_i}{\sum_{j} \Lambda_{0,j}(T)X_j}
$$

can be cast in terms of the absolute value of the least negative chemically significant eigenvalue of the master equation for the $i^{th}$ collider (when pure) in the low-pressure limit, $\Lambda_{0,i}(T)[M]$, and high-pressure limit, $\Lambda_{\infty}(T)$, and $[M]$ is the total concentration.

Evaluating all rate constants at the reduced pressure ($R$)---instead of the pressure ($P$)---of the mixture takes advantage of the fact that rate constants (and their chemically significant eigenvectors) for different colliders are usually far more similar at the same $R$ than the same $P$. In practice, since rate constants are usually expressed with respect to pressure $P$ (which has units of Pa, Torr, bar, atm, etc.) rather than reduced pressure $R$ (which is dimensionless), one needs to find the effective pressure for the $i^{th}$ collider, $P_{i}^{\text{ eff}}$ (with units of $P$), such that the reduced pressure of pure collider $i$ is equal to the reduced pressure of the mixture, which can be shown to be

$$
P_{i,\text{LMR}}^{\text{ eff}}(T,P,\textit{\textbf{X}}) = \frac{\sum_{j} \Lambda_{0,j}(T)X_j}{\Lambda_{0,i}(T)}P
$$

such that an alternate version of the generalized LMR-R equation can be written as

\begin{equation}
k_{\text{LMR-R}}(T,P,\textit{\textbf{X}}) = \sum_{i} k_{i}(T,P_{i,\text{LMR}}^{\text{ eff}})\tilde{X}_{i,\text{LMR}}
\label{eq:LMRR_k_P_i_eff}
\end{equation}

The user can either specify $k_i(T,P)$ and $\Lambda_{0,i}(T)$ for each collider, or specify the third-body efficiency, $\epsilon_{0,i}(T)$, for each non-M collider and assign $\epsilon_{0,\text{M}}(T)=1$.  The pressure-dependent aspect of each i-th collider ($k_i(T,P)$) can be specified in any combination of Troe, PLOG, or Chebyshev formats for colliders for which data are available and $\Lambda_{0,i}(T)$ for each collider (or $\epsilon_{0,i}(T)$ for each non-M collider) can be specified in modified Arrhenius format for colliders for which data are available.

In fact, if $\Lambda_{0,i}(T)$ (or $\epsilon_{0,i}(T)$) for any colliders have available data or can be estimated using typical values (as is typically done in kinetic models for reactions in modified Lindemann expressions) but no data for $k_i(T,P)$ are available, this implementation also allows some colliders to be specified with unique $\Lambda_{0,i}(T)$ (or $\epsilon_{0,i}(T)$) without $k_i(T,P)$, by assuming the same reduced-pressure dependence as M (i.e. $k_{i}(T,R)=k_{M}(T,R)$). Cantera employs the following derived form of the generalized LMR-R equation, where the sum over $n$ is only for the colliders for which unique $k_n(T,P)$ are available. If unique $k_i(T,P)$ data are available for all colliders, then the second term in the above equation effectively disappears.

$$
k_{\text{LMR-R}}(T,P,\textit{\textbf{X}}) &= \sum_{n} k_{n}(T,P_{n,\text{LMR}}^{\text{ eff}})\tilde{X}_{n,\text{LMR}} + k_{M}(T,P_{M,\text{LMR}}^{\text{ eff}}) \left(1-\sum_{n}\tilde{X}_{n,\text{LMR}}\right)
$$

If the user has limited or incomplete access to parameter inputs (a likely scenario, given the scarcity of puplished third-body efficiencies and master equation eigenvalues), this computational implementation allowes them much greater flexibility and power to make educated assumptions. Further description of the LMR-R theory and computational method is available here {cite:singal}`singal2024` [citation not yet added].

```{admonition} YAML Usage
:class: tip
Linear Burke rate expressions can be defined in the YAML format using the
[`linear-burke`](sec-yaml-linear-burke) reaction `type`.
```

```{versionadded} 3.1
```

(sec-blowers-masel)=
## Blowers-Masel Reactions

In some circumstances like thermodynamic sensitivity analysis, or modeling heterogeneous
reactions from one catalyst surface to another, the enthalpy change of a reaction
($\Delta H$) can be modified. Due to the change in $\Delta H$, the activation energy of
the reaction must be adjusted accordingly to provide accurate simulation results. To
adjust the activation energy  due to changes in the reaction enthalpy, the Blowers-Masel
rate expression is available. This approximation was proposed by {cite:t}`blowers2000`
to automatically scale activation energy as the reaction enthalpy is changed. The
_intrinsic activation energy_ $E_a^0$ is defined as the activation energy when
$\Delta H = 0$. The activation energy can then be written as a function of $\Delta H$:

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

and $w$ is the average of the bond dissociation energy of the bond breaking and that
being formed. Note that the expression is insensitive to $w$ as long as $w \ge 2 E_a^0$,
so we can use an arbitrarily high value of $w = 1000\text{ kJ/mol}$.

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
corrections dependent on the fractional surface coverages $\theta_k$ of one or more
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

$$  k_f = \frac{\gamma}{\Gamma_\t{tot}^m} \sqrt{\frac{RT}{2 \pi W}}  $$

where $\Gamma_\t{tot}$ is the total molar site density, $m$ is the sum of all the
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
