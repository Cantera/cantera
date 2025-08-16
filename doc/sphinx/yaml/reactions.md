(sec-yaml-reactions)=
# Reactions

The fields common to all `reaction` entries are:

`equation`
: The stoichiometric equation for the reaction. Each term (that is, stoichiometric
  coefficient, species name, `+` or `<=>`) in the equation must be separated by a space.

  Reversible reactions may be written using `<=>` or `=` to separate reactants and
  products. Irreversible reactions are written using `=>`.

`type`
: A string specifying the type of reaction or rate coefficient parameterization. The
  default is `elementary`. Reaction types are:
  - [`elementary`](sec-yaml-elementary)
  - [`three-body`](sec-yaml-three-body)
  - [`Blowers-Masel`](sec-yaml-Blowers-Masel)
  - [`two-temperature-plasma`](sec-yaml-two-temperature-plasma)
  - [`electron-collision-plasma`](sec-yaml-electron-collision-plasma)
  - [`electron-collisions`](sec-yaml-electron-collisions)
  - [`falloff`](sec-yaml-falloff)
  - [`chemically-activated`](sec-yaml-chemically-activated)
  - [`pressure-dependent-Arrhenius`](sec-yaml-pressure-dependent-Arrhenius)
  - [`Chebyshev`](sec-yaml-Chebyshev)
  - [`linear-Burke`](sec-yaml-linear-Burke)

  Reactions without a specified `type` on surfaces or edges are automatically treated as
  [`interface-Arrhenius`](sec-yaml-interface-Arrhenius) reactions, unless a
  `sticking-coefficient` implies a [`sticking-Arrhenius`](sec-yaml-sticking-Arrhenius)
  reaction. Interface reactions that involve charge transfer between phases are
  automatically treated as [`electrochemical`](sec-yaml-electrochemical-reaction)
  reactions.

  Reactions on surfaces or edges specifying `type` as `Blowers-Masel` are treated as
  [`interface-Blowers-Masel`](sec-yaml-interface-Blowers-Masel) or
  [`sticking-Blowers-Masel`](sec-yaml-sticking-Blowers-Masel).

`duplicate`
: Boolean indicating whether the reaction is a known duplicate of another reaction. The
  default is `false`.

`orders`
: An optional mapping of species to explicit reaction orders to use. Reaction orders for
  reactant species not explicitly mentioned are taken to be their respective
  stoichiometric coefficients. See [](sec-reaction-orders) for additional information.

`negative-orders`
: Boolean indicating whether negative reaction orders are allowed. The default is
  `false`.

`nonreactant-orders`
: Boolean indicating whether orders for non-reactant species are allowed. The default is
  `false`.

Depending on the reaction `type`, other fields may be necessary to specify the rate of
the reaction.

## Reaction rate expressions

(sec-yaml-arrhenius-rate)=
### Arrhenius

Arrhenius rate expressions are specified as a mapping with fields:

`A`
: The pre-exponential factor $A$

`b`
: The temperature exponent $b$

`Ea`
: The activation energy $E_a$

or a corresponding three-element list. The following are equivalent:

```yaml
{A: 2.70E+13 cm^3/mol/s, b: 0, Ea: 355 cal/mol}
[2.70E+13 cm^3/mol/s, 0, 355 cal/mol]
```

(sec-yaml-blowers-masel-rate)=
### Blowers-Masel

Blowers-Masel rate expressions calculate the rate constant based on the Blowers Masel
approximation as [described here](sec-blowers-masel). The rate parameters are specified
as a mapping with fields:

`A`
: The pre-exponential factor $A$

`b`
: The temperature exponent $b$

`Ea0`
: The intrinsic activation energy $E_{a0}$

`w`
: The average of the bond dissociation energy of the bond breaking and that being formed
  in the reaction $w$

or a corresponding four-element list. The following are equivalent:

```yaml
{A: 3.87e+04 cm^3/mol/s, b: 2.7, Ea0: 6260.0 cal/mol, w: 1e9 cal/mol}
[3.87e+04 cm^3/mol/s, 2.7, 6260.0 cal/mol, 1e9 cal/mol]
```

(sec-yaml-two-temperature-plasma-rate)=
### Two-Temperature Plasma

Two-temperature plasma reactions involve an electron as one of the reactants, where the
electron temperature may differ from the gas temperature as
[described here](sec-two-temperature-plasma-rate). The rate parameters are specified as
a mapping with fields:

`A`
: The pre-exponential factor

`b`
: The temperature exponent, which is applied to the electron temperature

`Ea-gas`
: The activation energy term $E_{a,g}$ that is related to the gas temperature

`Ea-electron`
: The activation energy term $E_{a,e}$ that is related to the electron temperature

or a corresponding four-element list. The following are equivalent:

```yaml
{A: 17283, b: -3.1, Ea-gas: -5820 J/mol, Ea-electron: 1081 J/mol}
[17283, -3.1, -5820 J/mol, 1081 J/mol]
```

(sec-yaml-efficiencies)=
## Efficiencies

Some reaction types include parameters for the "efficiency" of different species as
third-body colliders. For these reactions, the following additional fields are
supported:

`efficiencies`
: A mapping of species names to efficiency values

`default-efficiency`
: The efficiency for use for species not included in the `efficiencies` mapping.
  Defaults to 1.0.

Examples:

```yaml
- equation: O + H + M <=> OH + M
  rate-constant: {A: 5.0e+5, b: -1.0, Ea: 0.0}
  default-efficiency: 0.0
  efficiencies: {H2: 2.0, H2O: 6.0, AR: 0.7}
- equation: O + H + M <=> OH + M
  rate-constant: {A: 3.0e+5, b: 0.0, Ea: 123.0}
  # Efficiencies for the species in the previous reaction must be explicitly zero
  # to avoid double counting and being considered a duplicate reaction
  efficiencies: {H2: 0.0, H2O: 0.0, AR: 0.0, CO: 3.0}
```

(sec-yaml-rate-types)=
## Reaction types

(sec-yaml-elementary)=
### `elementary`

A homogeneous reaction with a pressure-independent rate coefficient and mass action
kinetics, as [described here](sec-arrhenius-rate).

Additional fields are:

`rate-constant`
: An [Arrhenius-type](sec-yaml-Arrhenius-rate) list or mapping.

`negative-A`
: A boolean indicating whether a negative value for the pre-exponential factor is
  allowed. The default is `false`.

Example:

```yaml
- equation: N + NO <=> N2 + O
  rate-constant: {A: -2.70000E+13 cm^3/mol/s, b: 0, Ea: 355 cal/mol}
  negative-A: true
```

(sec-yaml-three-body)=
### `three-body`

A three body reaction as [described here](sec-three-body-reaction).

The reaction equation must include a third body collision partner, which may be either a
specific species or the generic third body `M`.

Includes the fields of an `elementary` reaction, plus the fields for specifying
[efficiencies](sec-yaml-efficiencies).

Example:

```yaml
- equation: 2 O + M = O2 + M
  type: three-body
  rate-constant: [1.20000E+17 cm^6/mol^2/s, -1, 0]
  efficiencies: {AR: 0.83, H2O: 5}
```

The `type` field of the YAML entry may be omitted. Reactions containing the generic
third body M are automatically identified as three-body reactions. Reactions are also
identified as three-body reactions if all of the following conditions are met:

- There is exactly one species appearing as both a reactant and product
- All reactants and products have integral stoichiometric coefficients
- The sum of the stoichiometric coefficients for either the reactants or products is 3.

Examples:

```yaml
- equation: H + 2 O2 <=> HO2 + O2  # Reaction 34
  rate-constant: {A: 2.08e+19, b: -1.24, Ea: 0.0}
- equation: H + O2 + N2 <=> HO2 + N2  # Reaction 36
  rate-constant: {A: 2.6e+19, b: -1.24, Ea: 0.0}
```

````{caution}
If a corresponding reaction with the generic third body M also appears in the
mechanism, such as:

```yaml
- equation: H + O2 + M <=> HO2 + M  # Reaction 33
  rate-constant: {A: 2.8e+18, b: -0.86, Ea: 0.0}
  efficiencies: {O2: 0.0, H2O: 0.0, CO: 0.75, CO2: 1.5, C2H6: 1.5, N2: 0.0, AR: 0.0}
```

then the third body efficiency for any third bodies that are given in the explicit
form of Reaction 34 or Reaction 35 above must be set to zero, as shown here for
O2 and N2, or the reactions must be marked as duplicate.
````

:::{versionchanged} 3.0
Three body reactions are detected automatically and the the `type` field may be omitted.
Reactions with explicit third bodies are required to be marked as duplicates of
reactions with the generic third body if the corresponding efficiency is not zero.
:::

:::{versionadded} 3.1
Reactions with explicit third bodies and the corresponding reaction with "M" issue
warnings instead of raising errors by default. The
[](sec-yaml-phase-explicit-third-body-duplicates) field of the phase entry can be used
to control how these reactions are handled.
:::

(sec-yaml-blowers-masel)=
### `Blowers-Masel`

Includes the fields of an [`elementary`](sec-yaml-elementary) reaction, except that the
`rate-constant` field is a [Blowers-Masel-type](sec-yaml-Blowers-Masel-rate) list or
mapping.

Example:

```yaml
- equation: O + H2 <=> H + OH
  type: Blowers-Masel
  rate-constant: {A: 3.87e+04 cm^2/mol/s, b: 2.7, Ea0: 6260.0 cal/mol, w: 1e9 cal/mol}
```

(sec-yaml-two-temperature-plasma)=
### `two-temperature-plasma`

Includes the fields of an [`elementary`](sec-yaml-elementary) reaction, except that the
`rate-constant` field is a
[`two-temperature-plasma`](sec-yaml-two-temperature-plasma-rate) list or mapping.

Example:

```yaml
- equation: O + H => O + H
  type: two-temperature-plasma
  rate-constant: {A: 17283, b: -3.1, Ea-gas: -5820 J/mol, Ea-electron: 1081 J/mol}
```

(sec-yaml-electron-collision-plasma)=
### `electron-collision-plasma`

Electron collision plasma reactions involve an electron as one of the reactants, and are
parameterized by the collision cross section as a function of the electron energy. The
rate calculation is [described here](sec-electron-collision-plasma-rate). The rate
parameters are specified using the following additional fields in the reaction entry:

`energy-levels`
: A list of electron energy levels [V]

`cross-sections`
: A list of collision cross sections [m²] for the reaction at the specified energy
  levels.

Example:

```yaml
- equation: O2 + e => e + e + O2+
  type: electron-collision-plasma
  energy-levels: [13.0, 15.5, 18, 23]
  cross-sections: [1.17e-22, 7.3e-22, 1.64e-21, 3.66e-21]
```

:::{versionadded} 3.1
:::

(sec-yaml-electron-collisions)=
### `electron-collisions`

The `electron-collisions` field defines a list of cross-section datasets for
electron-impact processes that are used in plasma-phase simulations. These entries
are not formal reactions (they are not added to `Kinetics` objects), but serve
as data inputs for computing the electron energy distribution function.

Each entry includes:

`target`
: The name of the species that is the target of the collision

`energy-levels`
: A list of electron energy values [eV] at which the cross-section is provided

`cross-sections`
: Corresponding cross-section values [m²] for each energy level

`kind`
: A string indicating the process type. Options include:
  - `"effective"` – lumped or total effect of several channels
  - `"excitation"` – electronic excitation
  - `"ionization"` – electron-impact ionization
  - `"attachment"` – electron attachment processes

Example:

```yaml
electron-collisions:
- target: N2
  energy-levels: [0.0, 0.015, 0.03, 0.05, 0.1, 0.15, 0.2, 0.3, 0.4, 0.7, 1.2, 1.5, 1.9,
    2.2, 2.8, 3.3, 4.0, 5.0, 7.0, 10.0, 15.0, 20.0, 30.0, 75.0, 150.0]
  cross-sections: [1.1e-20, 2.55e-20, 3.4e-20, 4.33e-20, 5.95e-20, 7.1e-20, 7.9e-20,
    9e-20, 9.7e-20, 1e-19, 1.04e-19, 1.2e-19, 1.96e-19, 2.85e-19, 2.8e-19, 1.72e-19,
    1.26e-19, 1.09e-19, 1.01e-19, 1.04e-19, 1.1e-19, 1.02e-19, 9e-20, 6.6e-20, 4.9e-20]
  kind: effective
```

:::{versionadded} 3.2
:::

(sec-yaml-falloff)=
### `falloff`

A falloff reaction as [described here](sec-falloff-rate).

The reaction equation should include the pressure-dependent third body collision partner
`(+M)` or `(+name)` where `name` is the name of a species. The latter case is equivalent
to setting the efficiency for `name` to 1 and the efficiency for all other species to 0.

Includes field for specifying {ref}`efficiencies <sec-yaml-efficiencies>` as well as:

`high-P-rate-constant`
: An [](sec-yaml-Arrhenius-rate) expression for the high-pressure limit

`low-P-rate-constant`
: An [](sec-yaml-Arrhenius-rate) expression for the low-pressure limit

`Troe`
: Parameters for the [Troe](sec-troe-falloff) falloff function. A mapping containing the
  keys `A`, `T3`, `T1` and optionally `T2`. The default value for `T2` is 0.

`SRI`
: Parameters for the [SRI](sec-sri-falloff) falloff function. A mapping containing the
  keys `A`, `B`, `C`, and optionally `D` and `E`. The default values for `D` and `E` are
  1.0 and 0.0, respectively.

`Tsang`
: Parameters for the [Tsang](sec-tsang-falloff) falloff function. A mapping containing
  the keys `A` and `B`. The default value for `B` is 0.0.

Examples:

```yaml
- equation: H + CH2 (+ N2) <=> CH3 (+N2)
  type: falloff
  high-P-rate-constant: [6.00000E+14 cm^3/mol/s, 0, 0]
  low-P-rate-constant: {A: 1.04000E+26 cm^6/mol^2/s, b: -2.76, Ea: 1600}
  Troe: {A: 0.562, T3: 91, T1: 5836}
- equation: O + CO (+M) <=> CO2 (+M)
  type: falloff  # Lindemann form since no additional falloff function parameters given
  low-P-rate-constant: {A: 6.02e+14, b: 0.0, Ea: 3000.0}
  high-P-rate-constant: {A: 1.8e+10, b: 0.0, Ea: 2385.0}
  efficiencies: {H2: 2.0, O2: 6.0, H2O: 6.0, CH4: 2.0, CO: 1.5, CO2: 3.5,
    C2H6: 3.0, AR: 0.5}
```

(sec-yaml-chemically-activated)=
### `chemically-activated`

A chemically activated reaction as [described here](sec-chemically-activated-rate).

The parameters are the same as for {ref}`sec-yaml-falloff` reactions.

Example:

```yaml
- equation: CH3 + OH (+M) <=> CH2O + H2 (+M)
  type: chemically-activated
  high-P-rate-constant: [5.88E-14, 6.721, -3022.227]
  low-P-rate-constant: [282320.078, 1.46878, -3270.56495]
```

(sec-yaml-pressure-dependent-arrhenius)=
### `pressure-dependent-Arrhenius`

A pressure-dependent reaction using multiple Arrhenius expressions as
[described here](sec-plog-rate).

The only additional field in this reaction type is:

`rate-constants`
: A list of mappings, where each mapping is the mapping form of an
  [](sec-yaml-Arrhenius-rate) expression with the addition of a pressure `P`.

Example:

```yaml
- equation: H + CH4 <=> H2 + CH3
  type: pressure-dependent-Arrhenius
  rate-constants:
  - {P: 0.039474 atm, A: 2.720000e+09 cm^3/mol/s, b: 1.2, Ea: 6834.0}
  - {P: 1.0 atm, A: 1.260000e+20, b: -1.83, Ea: 15003.0}
  - {P: 1.0 atm, A: 1.230000e+04, b: 2.68, Ea: 6335.0}
  - {P: 1.01325 MPa, A: 1.680000e+16, b: -0.6, Ea: 14754.0}
```

(sec-yaml-chebyshev)=
### `Chebyshev`

A reaction parameterized as a bivariate Chebyshev polynomial as
[described here](sec-chebyshev-rate).

Additional fields are:

`temperature-range`
: A list of two values specifying the minimum and maximum temperatures at which the rate
  constant is valid

`pressure-range`
: A list of two values specifying the minimum and maximum pressures at which the rate
  constant is valid

`data`
: A list of lists containing the Chebyshev coefficients

Example:

```yaml
- equation: CH4 <=> CH3 + H
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

(sec-yaml-linear-burke)=
### `linear-Burke`

A complex-forming reaction (one that depends on both P and X) parameterized according to
the reduced-pressure linear mixture rule as [described here](sec-linear-Burke).

`efficiency` and `eig0` comprise the two acceptable ways to represent the contribution
of each bath gas component (collider) to the reduced pressure. All explicitly defined
colliders must include either `efficiency` or `eig0`, but the choice must remain
consistent throughout a single reaction (either all colliders are defined with
`efficiency`, or all are defined with `eig0`).

The pressure-dependent aspect of each collider rate constant can be parameterized in the
user's choice of [`Troe`](sec-yaml-falloff),
[`pressure-dependent-Arrhenius`](sec-yaml-pressure-dependent-Arrhenius), or
[`Chebyshev`](sec-yaml-Chebyshev) representations. The same parameters used for a
standalone Troe, PLOG, or Chebyshev reaction are then inserted directly below
`efficiency` or `eig0` for a given collider. At minimum, this treatment must be applied
to `M`. However, additional colliders may also be given their own Troe, PLOG, or
Chebyshev parameterization if so desired. Mixing and matching of types within the same
reaction is allowed (e.g., a PLOG table for `M`, Troe parameters for `H2`, and Chebyshev
data for `NH3`).

A mathematical description of this YAML implementation can be found in Eq. 8 of
{cite:t}`singal2024`.

Additional fields are:

`colliders`
: A list of dictionaries, where each entry contains parameters corresponding to
  individual colliders (species in the bath gas). Each entry within the `colliders` list
  may contain the following fields:

  `name`
  : The name of the collider species (e.g., `H2O`). The first collider defined must be
    `M`, which represents the generic reference collider (often `Ar` or `N2`) that
    represents all species lacking their own explicit parameterization.

  `eig0`
  : The absolute value of the least negative chemically significant eigenvalue of the
    master equation for the $i^{th}$ collider (when pure), evaluated at the low-pressure
    limit, $\Lambda_{0,i}(T)[M]$. The user must explicitly assign an `eig0` for `M`.
    This parameter is entered in modified Arrhenius format to enable consideration of
    temperature dependence.

  `efficiency`
  : The third-body efficiency of the collider relative to that of the reference collider
    `M`, defined as $\epsilon_{0,i}(T)=\Lambda_{0,i}(T)/\Lambda_{0,\text{M}}(T)$. The
    user does not need to assign an `efficiency` for `M`, as it is always equal to 1 by
    definition. However, they are free to do so, as long as it takes the form
    `efficiency: {A: 1, b: 0, Ea: 0}` (no variations are permitted). This parameter is
    entered in modified Arrhenius format to enable consideration of temperature
    dependence. If the user wishes to specify a temperature-independent value, then `A`
    can be set to this value and `b` and `Ea` can be set to 0 (such as
    `H2O: {A: 10, b: 0, Ea: 0}`).

  A [`Troe`](sec-yaml-falloff) implementation also requires: `high-P-rate-constant`,
  `low-P-rate-constant`, `Troe` (do not use the Troe `efficiencies` key).

  A [`pressure-dependent-Arrhenius`](sec-yaml-pressure-dependent-Arrhenius)
  implementation also requires: `rate-constants`.

  A [`Chebyshev`](sec-yaml-Chebyshev) implementation also requires: `temperature-range`,
  `pressure-range`, `data`.

Examples:

`linear-Burke` *rate with Troe format for the reference collider (N2)*:

```yaml
- equation: H + OH <=> H2O
  type: linear-Burke
  colliders:
  - name: M
    type: falloff
    low-P-rate-constant: {A: 4.530000e+21, b: -1.820309e+00, Ea: 4.987000e+02}
    high-P-rate-constant: {A: 2.510000e+13, b: 2.329303e-01, Ea: -1.142000e+02}
    Troe: {A: 9.995044e-01, T3: 1.0e-30, T1: 1.0e+30}
  - name: AR
    efficiency: {A: 2.20621e-02, b: 4.74036e-01, Ea: -1.13148e+02}
  - name: H2O
    efficiency: {A: 1.04529e-01, b: 5.50787e-01, Ea: -2.32675e+02}
```

`linear-Burke` *rate with PLOG format for the reference collider (Ar)*:

```yaml
- equation: H + O2 (+M) <=> HO2 (+M) # Adding '(+M)' is optional
  type: linear-Burke
  colliders:
  - name: M
    type: pressure-dependent-Arrhenius
    rate-constants:
    - {P: 1.316e-02 atm, A: 9.39968e+14, b: -2.14348e+00, Ea: 7.72730e+01}
    - {P: 1.316e-01 atm, A: 1.07254e+16, b: -2.15999e+00, Ea: 1.30239e+02}
    - {P: 3.947e-01 atm, A: 3.17830e+16, b: -2.15813e+00, Ea: 1.66994e+02}
    - {P: 1.000e+00 atm, A: 7.72584e+16, b: -2.15195e+00, Ea: 2.13473e+02}
    - {P: 3.000e+00 atm, A: 2.11688e+17, b: -2.14062e+00, Ea: 2.79031e+02}
    - {P: 1.000e+01 atm, A: 6.53093e+17, b: -2.13213e+00, Ea: 3.87493e+02}
    - {P: 3.000e+01 atm, A: 1.49784e+18, b: -2.10026e+00, Ea: 4.87579e+02}
    - {P: 1.000e+02 atm, A: 3.82218e+18, b: -2.07057e+00, Ea: 6.65984e+02}
  - name: HE
    efficiency: {A: 3.37601e-01, b: 1.82568e-01, Ea: 3.62408e+01}
  - name: N2
    efficiency: {A: 1.24932e+02, b: -5.93263e-01, Ea: 5.40921e+02}
  - name: H2
    efficiency: {A: 3.13717e+04, b: -1.25419e+00, Ea: 1.12924e+03}
  - name: CO2
    efficiency: {A: 1.62413e+08, b: -2.27622e+00, Ea: 1.97023e+03}
  - name: NH3
    efficiency: {A: 4.97750e+00, b: 1.64855e-01, Ea: -2.80351e+02}
  - name: H2O
    efficiency: {A: 3.69146e+01, b: -7.12902e-02, Ea: 3.19087e+01}
```

`linear-Burke` *rate with Chebyshev format for the reference collider (Ar)*:

```yaml
- equation: H2O2 <=> 2 OH
  type: linear-Burke
  colliders:
  - name: M
    type: Chebyshev
    temperature-range: [200.0, 2000.0]
    pressure-range: [1.000e-01 atm, 1.000e+02 atm]
    data:
    - [-1.58e+01, 8.71e-01, -9.44e-02, -2.81e-03, -4.48e-04, 1.58e-03, -2.51e-04]
    - [2.32e+01, 5.27e-01, 2.89e-02, -5.46e-03, 7.08e-04, -3.03e-03, 7.81e-04]
    - [-3.80e-01, 8.63e-02, 4.03e-02, -7.23e-03, 5.76e-04, 2.79e-03, -1.49e-03]
    - [-1.48e-01, -7.18e-03, 2.21e-02, 6.23e-03, -5.98e-03, -8.22e-06, 1.92e-03]
    - [-6.06e-02, -1.42e-02, 1.34e-03, 9.62e-03, 1.70e-03, -3.65e-03, -4.32e-04]
    - [-2.46e-02, -9.71e-03, -5.88e-03, 3.05e-03, 5.87e-03, 1.50e-03, -2.01e-03]
    - [-1.54e-02, -5.24e-03, -6.91e-03, -5.94e-03, -1.22e-03, 2.17e-03, 1.59e-03]
  - name: N2
    efficiency: {A: 1.14813e+00, b: 4.60090e-02, Ea: -2.92413e+00}
  - name: CO2
    efficiency: {A: 8.98839e+01, b: -4.27974e-01, Ea: 2.41392e+02}
  - name: H2O2
    efficiency: {A: 6.45295e-01, b: 4.26266e-01, Ea: 4.28932e+01}
  - name: H2O
    efficiency: {A: 1.36377e+00, b: 3.06592e-01, Ea: 2.10079e+02}
```

(sec-yaml-interface-arrhenius)=
### `interface-Arrhenius`

A reaction occurring on a surface between two bulk phases, or along an edge at the
intersection of two surfaces, as [described here](sec-surface-rate).

Includes the fields of an [](sec-yaml-elementary) reaction plus:

`coverage-dependencies`
: A mapping of species names to coverage dependence parameters, where these
  parameters are contained in either a mapping with the fields:

  `a`
  : Coefficient for exponential dependence on the coverage

  `m`
  : Power-law exponent of coverage dependence

  `E`
  : Activation energy dependence on coverage, which uses the same sign convention as the
    leading-order activation energy term. This can be a scalar value for the linear
    dependency or a list of four values for the polynomial dependency given in the order
    of 1st, 2nd, 3rd, and 4th-order coefficients

  or a list containing the three elements above, in the given order.

  Note that parameters `a`, `m` and `E` correspond to parameters $\eta_{ki}$, $\mu_{ki}$
  and $\epsilon_{ki}$ in Eq 11.113 of {cite:t}`kee2003`, respectively.

Examples:

```yaml
- equation: 2 H(s) => H2 + 2 Pt(s)
  rate-constant: {A: 3.7e21 cm^2/mol/s, b: 0, Ea: 67400 J/mol}
  coverage-dependencies: {H(s): {a: 0, m: 0, E: -6000 J/mol}}

- equation: 2 O(S) => O2 + 2 Pt(S)
  rate-constant: {A: 3.7e+21, b: 0, Ea: 213200 J/mol}
  coverage-dependencies: {O(S): {a: 0.0, m: 0.0,
    E: [1.0e3 J/mol, 3.0e3 J/mol , -7.0e4 J/mol , 5.0e3 J/mol]}}

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

(sec-yaml-interface-blowers-masel)=
### `interface-Blowers-Masel`

Includes the same fields as [`interface-Arrhenius`](sec-yaml-interface-Arrhenius), while
using the [Blowers-Masel](sec-yaml-Blowers-Masel-rate) parameterization for the rate
constant.

Example:

```yaml
- equation: 2 H(s) => H2 + 2 Pt(s)
  type: Blowers-Masel
  rate-constant: {A: 3.7e21 cm^2/mol/s, b: 0, Ea0: 67400 J/mol, w: 1000000 J/mol}
  coverage-dependencies: {H(s): {a: 0, m: 0, E: -6000 J/mol}}
```

(sec-yaml-sticking-arrhenius)=
### `sticking-Arrhenius`

A sticking reaction occurring on a surface adjacent to a bulk phase, as
[described here](sec-sticking-rate).

Includes the fields of an [](sec-yaml-interface-Arrhenius) reaction plus:

`sticking-coefficient`
: An [Arrhenius-type](sec-yaml-Arrhenius-rate) expression for the sticking coefficient

`Motz-Wise`
: A boolean indicating whether to use the Motz-Wise correction factor for sticking
  coefficients near unity. Defaults to `false`.

`sticking-species`
: The name of the sticking species. Required if the reaction includes multiple
  non-surface species.

Example:

```yaml
- equation: OH + PT(S) => OH(S)
  sticking-coefficient: {A: 1.0, b: 0, Ea: 0}
```

(sec-yaml-sticking-blowers-masel)=
### `sticking-Blowers-Masel`

Includes the same fields as [`sticking-Arrhenius`](sec-yaml-sticking-Arrhenius), while
using the [Blowers-Masel](sec-yaml-Blowers-Masel-rate) parameterization
for the sticking coefficient.

Example:

```yaml
- equation: OH + PT(S) => OH(S)
  type: Blowers-Masel
  sticking-coefficient: {A: 1.0, b: 0, Ea0: 0, w: 100000}
  Motz-Wise: true
```

(sec-yaml-electrochemical-reaction)=
### `electrochemical`

Interface reactions involving [charge transfer](sec-electrochemical-reactions) between
phases.

Includes the fields of an [](sec-yaml-interface-Arrhenius) reaction, plus:

`beta`
: The symmetry factor for the reaction. Default is 0.5.

`exchange-current-density-formulation`
: Set to `true` if the rate constant parameterizes the exchange current density. Default
  is `false`.

Example:

```yaml
- equation: LiC6 <=> Li+(e) + C6
  rate-constant: [5.74, 0.0, 0.0]
  beta: 0.4
```
