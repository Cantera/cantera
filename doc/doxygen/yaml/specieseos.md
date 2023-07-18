# Species Equation of State Models {#sec-yaml-species-eos-models}

[TOC]

## constant-volume {#sec-yaml-eos-constant-volume}

A constant volume model, see @ref Cantera.PDSS_ConstVol.

Any one of the following may be specified:

`molar-volume`

The molar volume of the species.

`molar-density`

The molar density of the species.

`density`

The mass density of the species.

**Example:**

``` yaml
equation-of-state:
  model: constant-volume
  molar-volume: 1.3 cm^3/mol
```

## coverage-dependent-surface {#sec-yaml-coverage-dependent-surface-species}

A model where species thermodynamic properties are calculated as a
function coverage, see @ref Cantera.CoverageDependentSurfPhase.

Additional fields:

`coverage-dependencies`

Mapping where keys are the name of species whose coverage affects
thermodynamic properties of the node-owner species. The map values
are the dependency entries including `model`, model-specific
parameters, `heat-capacity-a`, and `heat-capacity-b` that correspond
to an individual dependency between the node-owner species and keyed species.

`model`

Dependency model for coverage-dependent enthalpy or entropy. It
should be one of the four: `linear`, `polynomial`,
`piecewise-linear` or `interpolative`. The `model` and
model-specific parameters are grouped as follow.

-   `linear`: `enthalpy`, `entropy`
-   `polynomial`: `enthalpy-coefficients`, `entropy-coefficients`
-   `piecewise-linear`: `enthalpy-low`, `enthalpy-high`,
    `enthalpy-change`, `entropy-low`, `entropy-high`, `entropy-change`
-   `interpolative`: `enthalpy-coverages`, `enthalpies`,
    `entropy-coverages`, `entropies`

`enthalpy` or `entropy`

Slope of the coverage-dependent enthalpy or entropy used in the
`linear` model.

`enthalpy-coefficients` or `entropy-coefficients`

Array of polynomial coefficients in order of 1st, 2nd, 3rd, and
4th-order used in coverage-dependent enthalpy or entropy calculation
with the `polynomial` model.

`enthalpy-low` or `entropy-low`

Slope of the coverage-dependent enthalpy or entropy for the lower
coverage region used in the `piecewise-linear` model.

`enthalpy-high` or `entropy-high`

Slope of the coverage-dependent enthalpy or entropy for the higher
coverage region used in the `piecewise-linear` model.

`enthalpy-change` or `entropy-change`

Coverage that separates the lower and higher coverage regions of the
coverage-dependent enthalpy or entropy used in the
`piecewise-linear` model.

`enthalpy-coverages` or `entropy-coverages`

Array of discrete coverage values used in coverage-dependent
enthalpy or entropy used in the `interpolative` model.

`enthalpies` or `entropies`

Array of discrete enthalpy or entropy values corresponding to the
coverages in `enthalpy-coverages` or `entropy-coverages`,
respectively, used in the `interpolative` model.

`heat-capacity-a` or `heat-capacity-b`

Coefficient \f$c^{(a)}\f$ or \f$c^{(b)}\f$ used in the coverage-dependent
`heat capacity` model.

**Example:**

``` yaml
coverage-dependencies:
  OC_Pt: {model: linear,
          units: {energy: eV, quantity: molec},
          enthalpy: 0.48, entropy: -0.031}
  C_Pt: {model: polynomial,
         units: {energy: J, quantity: mol},
         enthalpy-coefficients: [0.0, -3.86e4, 0.0, 4.2e5],
         entropy-coefficients: [0.8e3, 0.0, -1.26e4, 0.0]}
  CO2_Pt: {model: piecewise-linear,
           units: {energy: kJ, quantity: mol},
           enthalpy-low: 0.5e2, enthalpy-high: 1.0e2,
           enthalpy-change: 0.4,
           entropy-low: 0.1e2, entropy-high: -0.2e2,
           entropy-change: 0.4,
           heat-capacity-a: 0.02e-1, heat-capacity-b: -0.156e-1}
  O_Pt: {model: interpolative,
         units: {energy: kcal, quantity: mol},
         enthalpy-coverages: [0.0, 0.2, 0.4, 0.7, 0.9, 1.0],
         enthalpies: [0.0, 0.5, 1.0, 2.7, 3.5, 4.0],
         entropy-coverages: [0.0, 0.5, 1.0],
         entropies: [0.0, -0.7, -2.0]}
```

## density-temperature-polynomial {#sec-yaml-eos-density-temperature-polynomial}

A model in which the density varies with temperature, see @ref Cantera.PDSS_SSVol.

Additional fields:

`data`

Vector of 4 coefficients for a cubic polynomial in temperature

**Example:**

``` yaml
equation-of-state:
  model: density-temperature-polynomial
  units: {mass: g, length: cm}
  data: [0.536504, -1.04279e-4, 3.84825e-9, -5.2853e-12]
```

## HKFT {#sec-yaml-eos-hkft}

The Helgeson-Kirkham-Flowers-Tanger model, see @ref Cantera.PDSS_HKFT.

Additional fields:

`h0`

Enthalpy of formation at the reference temperature and pressure

`s0`

Entropy of formation at the reference temperature and pressure

`a`

4-element vector containing the coefficients \f$a_1, \ldots , a_4\f$

`c`

2-element vector containing the coefficients \f$c_1\f$ and \f$c_2\f$

`omega`

The \f$\omega\f$ parameter at the reference temperature and pressure

**Example:**

``` yaml
equation-of-state:
  model: HKFT
  h0: -57433. cal/gmol
  s0: 13.96 cal/gmol/K
  a: [0.1839 cal/gmol/bar, -228.5 cal/gmol,
     3.256 cal*K/gmol/bar, -27260. cal*K/gmol]
  c: [18.18 cal/gmol/K, -29810. cal*K/gmol]
  omega: 33060 cal/gmol
```

## ideal-gas {#sec-yaml-eos-ideal-gas}

A species using the ideal gas equation of state, see @ref Cantera.PDSS_IdealGas.

@deprecated This species thermo model is deprecated and will be removed after %Cantera 3.0.

## ions-from-neutral-molecule {#sec-yaml-eos-ions-from-neutral}

A species equation of state model used with the
`ions-from-neutral-molecule` phase model, see @ref Cantera.PDSS_IonsFromNeutral.

@deprecated This species thermo model is deprecated and will be removed after %Cantera 3.0.

Additional fields:

`special-species`

Boolean indicating whether the species is the \"special species\" in
the phase. Default is `false`.

`multipliers`

A dictionary mapping species to neutral species multiplier values.

**Example:**

``` yaml
equation-of-state:
  model: ions-from-neutral-molecule
  multipliers: {KCl(l): 1.2}
```

## liquid-water-IAPWS95 {#sec-yaml-eos-liquid-water-iapws95}

A detailed equation of state for liquid water, see @ref Cantera.PDSS_Water.  as [described
here](https://cantera.org/documentation/dev/doxygen/html/de/d64/classCantera_1_1PDSS__Water.html#details).

## molar-volume-temperature-polynomial {#sec-yaml-eos-molar-volume-temperature-polynomial}

A model in which the molar volume varies with temperature, see @ref Cantera.PDSS_SSVol.

Additional fields:

`data`

Vector of 4 coefficients for a cubic polynomial in temperature

## Peng-Robinson {#sec-yaml-eos-peng-robinson}

A model where species follow the Peng-Robinson equation of state, see
@ref Cantera.PengRobinson.

Additional fields:

`a`

Pure-species `a` coefficient \[Pa\*m\^6/kmol\^2\]

`b`

Pure-species `b` coefficient \[m\^3/kmol\]

`acentric-factor`

Pitzer's acentric factor \[-\]

`binary-a`

Optional mapping where the keys are species names and the values are
the `a` coefficients for binary interactions between the two species.

**Example:**

``` yaml
equation-of-state:
  model: Peng-Robinson
  units: {length: cm, quantity: mol}
  a: 5.998873E+11
  b: 18.9714
  acentric-factor: 0.344
  binary-a:
    H2: 4 bar*cm^6/mol^2
    CO2: 7.897e7 bar*cm^6/mol^2
```

## Redlich-Kwong {#sec-yaml-eos-redlich-kwong}

A model where species follow the Redlich-Kwong equation of state, see
@ref Cantera.RedlichKwongMFTP.

Additional fields:

`a`

Pure-species `a` coefficient. Scalar or list of two values for a
temperature-dependent expression.

`b`

Pure-species `b` coefficient.

`binary-a`

Mapping where the keys are species and the values are the `a`
coefficients for binary interactions between the two species.
