@page sec-yaml-species-coverage-models Species Coverage Models

[TOC]

Species coverage models are specified using the `coverage-dependencies` field of
@ref sec-yaml-species.

# linear {#sec-yaml-species-coverage-linear}

<b>`enthalpy`</b> and <b>`entropy`</b>

Slope of the coverage-dependent enthalpy or entropy.

**Example:**

```yaml
coverage-dependencies:
  O_Pt:
    model: linear
    units: {energy: eV, quantity: molec}
    enthalpy: 0.48
    entropy: -0.031
  # + other entries (optional)
```

# polynomial {#sec-yaml-species-coverage-polynomial}

<b>`enthalpy-coefficients`</b> and <b>`entropy-coefficients`</b>

Array of polynomial coefficients in order of 1st, 2nd, 3rd, and
4th-order used in coverage-dependent enthalpy or entropy calculation.

**Example:**

```yaml
coverage-dependencies:
  OC_Pt:
    model: polynomial
    units: {energy: J, quantity: mol}
    enthalpy-coefficients: [0.0, -3.86e4, 0.0, 4.2e5]
    entropy-coefficients: [0.8e3, 0.0, -1.26e4, 0.0]
  # + other entries (optional)
```

# piecewise-linear {#sec-yaml-species-coverage-piecewise-linear}

<b>`enthalpy-low`</b> and <b>`entropy-low`</b>

Slope of the coverage-dependent enthalpy or entropy for the lower coverage region.

<b>`enthalpy-high`</b> and <b>`entropy-high`</b>

Slope of the coverage-dependent enthalpy or entropy for the higher coverage region.

<b>`enthalpy-change`</b> and <b>`entropy-change`</b>

Coverage that separates the lower and higher coverage regions of the
coverage-dependent enthalpy or entropy.

<b>`heat-capacity-a`</b> and <b>`heat-capacity-b`</b>

Coefficient @f$c^{(a)}@f$ or @f$c^{(b)}@f$ used in the coverage-dependent
`heat capacity` model.

**Example:**

```yaml
coverage-dependencies:
  CO2_Pt:
    model: piecewise-linear
    units: {energy: kJ, quantity: mol}
    enthalpy-low: 0.5e2
    enthalpy-high: 1.0e2
    enthalpy-change: 0.4
    entropy-low: 0.1e2
    entropy-high: -0.2e2
    entropy-change: 0.4
    heat-capacity-a: 0.02e-1
    heat-capacity-b: -0.156e-1
  # + other entries (optional)
```

# interpolative {#sec-yaml-species-coverage-interpolative}

<b>`enthalpy-coverages`</b> and <b>`entropy-coverages`</b>

Array of discrete coverage values used in coverage-dependent enthalpy or entropy.

<b>`enthalpies`</b> and <b>`entropies`</b>

Array of discrete enthalpy or entropy values corresponding to the
coverages in `enthalpy-coverages` and `entropy-coverages`, respectively.

**Example:**

```yaml
coverage-dependencies:
  C_Pt:
    model: interpolative
    units: {energy: kcal, quantity: mol}
    enthalpy-coverages: [0.0, 0.2, 0.4, 0.7, 0.9, 1.0]
    entropy-coverages: [0.0, 0.5, 1.0]
    enthalpies: [0.0, 0.5, 1.0, 2.7, 3.5, 4.0]
    entropies: [0.0, -0.7, -2.0]
  # + other entries (optional)
```
