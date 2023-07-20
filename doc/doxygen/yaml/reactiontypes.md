
@page sec-yaml-reaction-models Reaction Models

[TOC]

# elementary {#sec-yaml-elementary}

A homogeneous reaction with a pressure-independent rate coefficient and
mass action kinetics, as [described
here](https://cantera.org/science/kinetics.html#reactions-with-a-pressure-independent-rate).

Additional fields are:

<b>`rate-constant`</b>

An @ref sec-yaml-Arrhenius-rate list or mapping.

<b>`negative-A`</b>

A boolean indicating whether a negative value for the
pre-exponential factor is allowed. The default is `false`.

**Example:**

```yaml
equation: N + NO <=> N2 + O
rate-constant: {A: -2.70000E+13 cm^3/mol/s, b: 0, Ea: 355 cal/mol}
negative-A: true
```

# three-body {#sec-yaml-three-body}

A three body reaction as [described
here](https://cantera.org/science/kinetics.html#three-body-reactions).

The reaction equation should include the third body collision partner `M`.

Includes the fields of an <tt>@ref sec-yaml-elementary</tt> reaction, plus the fields
for specifying @ref sec-yaml-reactions-efficiencies.

**Example:**

```yaml
equation: 2 O + M = O2 + M
type: three-body
rate-constant: [1.20000E+17 cm^6/mol^2/s, -1, 0]
efficiencies: {AR: 0.83, H2O: 5}
```

# Blowers-Masel {#sec-yaml-Blowers-Masel}

Includes the fields of an <tt>@ref sec-yaml-elementary</tt> reaction, except that the
`rate-constant` field is a @ref sec-yaml-Blowers-Masel-rate list or mapping.

**Example:**

```yaml
equation: O + H2 <=> H + OH
type: Blowers-Masel
rate-constant: {A: 3.87e+04 cm^2/mol/s, b: 2.7, Ea0: 6260.0 cal/mol, w: 1e9 cal/mol}
```

# two-temperature-plasma {#sec-yaml-two-temperature-plasma}

Includes the fields of an <tt>@ref sec-yaml-elementary</tt> reaction, except that the
`rate-constant` field is a @ref sec-yaml-two-temperature-plasma-rate list or mapping.

**Example:**

```yaml
equation: O + H => O + H
type: two-temperature-plasma
rate-constant: {A: 17283, b: -3.1, Ea-gas: -5820 J/mol, Ea-electron: 1081 J/mol}
```

# falloff {#sec-yaml-falloff}

A falloff reaction as [described
here](https://cantera.org/science/kinetics.html#falloff-reactions).

The reaction equation should include the pressure-dependent third body
collision partner `(+M)` or `(+name)` where `name` is the name of a
species. The latter case is equivalent to setting the efficiency for
`name` to 1 and the efficiency for all other species to 0.

Includes field for specifying @ref sec-yaml-reactions-efficiencies as well as:

<b>`high-P-rate-constant`</b>

An @ref sec-yaml-Arrhenius-rate expression for the high-pressure limit

<b>`low-P-rate-constant`</b>

An @ref sec-yaml-Arrhenius-rate expression for the low-pressure limit

<b>`Troe`</b>

Parameters for the
[Troe](https://cantera.org/science/kinetics.html#the-troe-falloff-function)
falloff function. A mapping containing the keys `A`, `T3`, `T1` and
optionally `T2`. The default value for `T2` is 0.

<b>`SRI`</b>

Parameters for the
[SRI](https://cantera.org/science/kinetics.html#the-sri-falloff-function)
falloff function. A mapping containing the keys `A`, `B`, `C`, and
optionally `D` and `E`. The default values for `D` and `E` are 1.0
and 0.0, respectively.

<b>`Tsang`</b>

Parameters for the
[Tsang](https://cantera.org/science/kinetics.html#tsang-s-approximation-to-f-cent)
falloff function. A mapping containing the keys `A` and `B`. The
default value for `B` is 0.0.

**Example:**

```yaml
equation: H + CH2 (+ N2) <=> CH3 (+N2)
type: falloff
high-P-rate-constant: [6.00000E+14 cm^3/mol/s, 0, 0]
low-P-rate-constant: {A: 1.04000E+26 cm^6/mol^2/s, b: -2.76, Ea: 1600}
Troe: {A: 0.562, T3: 91, T1: 5836}
```

# chemically-activated {#sec-yaml-chemically-activated}

A chemically activated reaction as [described
here](https://cantera.org/science/kinetics.html#chemically-activated-reactions).

The parameters are the same as for <tt>@ref sec-yaml-falloff</tt> reactions.

**Example:**

```yaml
equation: CH3 + OH (+M) <=> CH2O + H2 (+M)
type: chemically-activated
high-P-rate-constant: [5.88E-14, 6.721, -3022.227]
low-P-rate-constant: [282320.078, 1.46878, -3270.56495]
```

# pressure-dependent-Arrhenius {#sec-yaml-pressure-dependent-Arrhenius}

A pressure-dependent reaction using multiple Arrhenius expressions as
[described
here](https://cantera.org/science/kinetics.html#pressure-dependent-arrhenius-rate-expressions-p-log).

The only additional field in this reaction type is:

<b>`rate-constants`</b>

A list of mappings, where each mapping is the mapping form of an
@ref sec-yaml-Arrhenius-rate expression with the addition of a pressure `P`.

**Example:**

```yaml
equation: H + CH4 <=> H2 + CH3
type: pressure-dependent-Arrhenius
rate-constants:
- {P: 0.039474 atm, A: 2.720000e+09 cm^3/mol/s, b: 1.2, Ea: 6834.0}
- {P: 1.0 atm, A: 1.260000e+20, b: -1.83, Ea: 15003.0}
- {P: 1.0 atm, A: 1.230000e+04, b: 2.68, Ea: 6335.0}
- {P: 1.01325 MPa, A: 1.680000e+16, b: -0.6, Ea: 14754.0}
```

# Chebyshev {#sec-yaml-Chebyshev}

A reaction parameterized as a bivariate Chebyshev polynomial as
[described
here](https://cantera.org/science/kinetics.html#chebyshev-reaction-rate-expressions).

Additional fields are:

<b>`temperature-range`</b>

A list of two values specifying the minimum and maximum temperatures
at which the rate constant is valid

<b>`pressure-range`</b>

A list of two values specifying the minimum and maximum pressures at
which the rate constant is valid

<b>`data`</b>

A list of lists containing the Chebyshev coefficients

**Example:**

```yaml
equation: CH4 <=> CH3 + H
type: Chebyshev
temperature-range: [290, 3000]
pressure-range: [0.0098692326671601278 atm, 98.692326671601279 atm]
data: [[-1.44280e+01,  2.59970e-01, -2.24320e-02, -2.78700e-03],
       [ 2.20630e+01,  4.88090e-01, -3.96430e-02, -5.48110e-03],
       [-2.32940e-01,  4.01900e-01, -2.60730e-02, -5.04860e-03],
       [-2.93660e-01,  2.85680e-01, -9.33730e-03, -4.01020e-03],
       [-2.26210e-01,  1.69190e-01,  4.85810e-03, -2.38030e-03],
       [-1.43220e-01,  7.71110e-02,  1.27080e-02, -6.41540e-04]]
```

# interface-Arrhenius {#sec-yaml-interface-Arrhenius}

A reaction occurring on a surface between two bulk phases, or along an
edge at the intersection of two surfaces, as [described
here](https://cantera.org/science/kinetics.html#sec-surface).

Includes the fields of an <tt>@ref sec-yaml-elementary</tt> reaction plus:

<b>`coverage-dependencies`</b>

A mapping of species names to coverage dependence parameters, where
these parameters are contained in either a mapping with the fields:

-   `a`: Coefficient for exponential dependence on the coverage
-   `m`: Power-law exponent of coverage dependence
-   `E`: Activation energy dependence on coverage, which uses the same
    sign convention as the leading-order activation energy term.
    This can be a scalar value for the linear dependency or a list
    of four values for the polynomial dependency given in the order
    of 1st, 2nd, 3rd, and 4th-order coefficients

or a list containing the three elements above, in the given order.

Note that parameters `a`, `m` and `E` correspond to parameters
@f$ \eta_{ki}@f$ , @f$ \mu_{ki}@f$  and @f$ \epsilon_{ki}@f$  in Eq 11.113 of
\[Kee, R. J., Coltrin, M. E., & Glarborg, P.(2003). Chemically
reacting flow: theory and practice. John Wiley & Sons\],
respectively.

**Examples:**

```yaml
- equation: 2 H(s) => H2 + 2 Pt(s)
  rate-constant: {A: 3.7e21 cm^2/mol/s, b: 0, Ea: 67400 J/mol}
  coverage-dependencies: {H(s): {a: 0, m: 0, E: -6000 J/mol}}

- equation: 2 O(S) => O2 + 2 Pt(S)
  rate-constant: {A: 3.7e+21, b: 0, Ea: 213200 J/mol}
  coverage-dependencies: {O(S): {a: 0.0, m: 0.0,
    E: [1.0e3 J/mol, 3.0e3 J/mol , -7.0e4 J/mol , 5.0e3 J/mol]}

- equation: CH4 + PT(S) + O(S) => CH3(S) + OH(S)
  rate-constant: {A: 5.0e+18, b: 0.7, Ea: 4.2e+04}
  coverage-dependencies:
    O(S): [0, 0, 8000]
    PT(S): [0, -1.0, 0]

- equation: 2 O(S) => O2 + 2 Pt(S)
  rate-constant: {A: 3.7e+21, b: 0, Ea: 213200 J/mol}
  coverage-dependencies:
    O(S): [0, 0, [1.0e6, 3.0e6, -7.0e7, 5.0e6]]
```

# interface-Blowers-Masel {#sec-yaml-interface-Blowers-Masel}

Includes the same fields as <tt>@ref sec-yaml-interface-Arrhenius</tt>, while using the
@ref sec-yaml-Blowers-Masel-rate parameterization for the rate constant.

**Example:**

```yaml
equation: 2 H(s) => H2 + 2 Pt(s)
type: Blowers-Masel
rate-constant: {A: 3.7e21 cm^2/mol/s, b: 0, Ea0: 67400 J/mol, w: 1000000 J/mol}
coverage-dependencies: {H(s): {a: 0, m: 0, E: -6000 J/mol}}
```

# sticking-Arrhenius {#sec-yaml-sticking-Arrhenius}

A sticking reaction occurring on a surface adjacent to a bulk phase, as
[described here](https://cantera.org/science/kinetics.html#sec-sticking).

Includes the fields of an <tt>@ref sec-yaml-interface-Arrhenius</tt> reaction plus:

<b>`sticking-coefficient`</b>

An @ref sec-yaml-Arrhenius-rate expression for the sticking coefficient

<b>`Motz-Wise`</b>

A boolean indicating whether to use the Motz-Wise correction factor
for sticking coefficients near unity. Defaults to `false`.

<b>`sticking-species`</b>

The name of the sticking species. Required if the reaction includes
multiple non-surface species.

**Example:**

```yaml
equation: OH + PT(S) => OH(S)
sticking-coefficient: {A: 1.0, b: 0, Ea: 0}
```

# sticking-Blowers-Masel {#sec-yaml-sticking-Blowers-Masel}

Includes the same fields as <tt>@ref sec-yaml-sticking-Arrhenius</tt>, while using the
@ref sec-yaml-Blowers-Masel-rate parameterization for the sticking coefficient.

**Example:**

```yaml
equation: OH + PT(S) => OH(S)
type: Blowers-Masel
sticking-coefficient: {A: 1.0, b: 0, Ea0: 0, w: 100000}
Motz-Wise: true
```

# electrochemical {#sec-yaml-electrochemical-reaction}

Interface reactions involving charge transfer between phases.

Includes the fields of an <tt>@ref sec-yaml-interface-Arrhenius</tt> reaction, plus:

<b>`beta`</b>

The symmetry factor for the reaction. Default is 0.5.

<b>`exchange-current-density-formulation`</b>

Set to `true` if the rate constant parameterizes the exchange
current density. Default is `false`.

**Example:**

```yaml
equation: LiC6 <=> Li+(e) + C6
rate-constant: [5.74, 0.0, 0.0]
beta: 0.4
```
