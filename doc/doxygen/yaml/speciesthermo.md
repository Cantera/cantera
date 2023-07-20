# Species Thermodynamic Models {#sec-yaml-species-thermo-models}

[TOC]

Species thermodynamic models are specified using the `thermo` field of
@ref sec-yaml-species.

## constant-cp {#sec-yaml-constant-cp}

The constant heat capacity model, see @ref Cantera.ConstCpPoly or [described
here](https://cantera.org/science/species-thermo.html#constant-heat-capacity).

Additional fields of a `constant-cp` thermo entry are:

<b>`T0`</b>

The reference temperature. Defaults to 298.15 K.

<b>`h0`</b>

The molar enthalpy at the reference temperature. Defaults to 0.0.

<b>`s0`</b>

The molar entropy at the reference temperature. Defaults to 0.0.

<b>`cp0`</b>

The heat capacity at constant pressure. Defaults to 0.0.

<b>`T-min`</b>

The minimum temperature at which this thermo data should be used.
Defaults to 0.0.

<b>`T-max`</b>

The maximum temperature at which this thermo data should be used.
Defaults to infinity.

**Example:**

``` yaml
thermo:
  model: constant-cp
  T0: 1000 K
  h0: 9.22 kcal/mol
  s0: -3.02 cal/mol/K
  cp0: 5.95 cal/mol/K
```

## NASA7 {#sec-yaml-nasa7}

The NASA 7-coefficient polynomial form given for one (@ref Cantera.NasaPoly1) or two
temperature regions (@ref Cantera.NasaPoly2); also [described
here](https://cantera.org/science/species-thermo.html#the-nasa-7-coefficient-polynomial-parameterization).

Additional fields of a `NASA7` thermo entry are:

<b>`temperature-ranges`</b>

A list giving the temperature intervals on which the polynomials are
valid. For one temperature region, this list contains the minimum
and maximum temperatures for the polynomial. For two temperature
regions, this list contains the minimum, intermediate, and maximum
temperatures.

<b>`data`</b>

A list with one item per temperature region, where that item is a 7
item list of polynomial coefficients. The temperature regions are
arranged in ascending order. Note that this is different from the
standard CHEMKIN formulation that uses two temperature regions
listed in descending order.

**Example:**

``` yaml
thermo:
  model: NASA7
  temperature-ranges: [300.0, 1000.0, 5000.0]
  data:
  - [3.298677, 0.0014082404, -3.963222e-06, 5.641515e-09,
    -2.444854e-12, -1020.8999, 3.950372]
  - [2.92664, 0.0014879768, -5.68476e-07, 1.0097038e-10,
    -6.753351e-15, -922.7977, 5.980528]
```

## NASA9 {#sec-yaml-nasa9}

The NASA 9-coefficient polynomial form given for a single (@ref Cantera.Nasa9Poly1) or
multiple temperature regions (@ref Cantera.Nasa9PolyMultiTempRegion); also [described
here](https://cantera.org/science/species-thermo.html#the-nasa-9-coefficient-polynomial-parameterization).

Additional fields of a `NASA9` thermo entry are:

<b>`temperature-ranges`</b>

A list giving the temperature intervals on which the polynomials are
valid. This list contains the minimum temperature, the intermediate
temperatures between each set pair of regions, and the maximum
temperature.

<b>`data`</b>

A list with one item per temperature region, where that item is a 9
item list of polynomial coefficients. The temperature regions are
arranged in ascending order.

**Example:**

``` yaml
thermo:
  model: NASA9
  temperature-ranges: [200.00, 1000.00, 6000.0, 20000]
  reference-pressure: 1 bar
  data:
  - [2.210371497E+04, -3.818461820E+02, 6.082738360E+00, -8.530914410E-03,
     1.384646189E-05, -9.625793620E-09, 2.519705809E-12, 7.108460860E+02,
     -1.076003744E+01]
  - [5.877124060E+05, -2.239249073E+03, 6.066949220E+00, -6.139685500E-04,
     1.491806679E-07,  -1.923105485E-11, 1.061954386E-15, 1.283210415E+04,
     -1.586640027E+01]
  - [8.310139160E+08, -6.420733540E+05, 2.020264635E+02, -3.065092046E-02,
     2.486903333E-06, -9.705954110E-11, 1.437538881E-15, 4.938707040E+06,
     -1.672099740E+03]
```

## piecewise-Gibbs {#sec-yaml-piecewise-gibbs}

A model based on piecewise interpolation of the Gibbs free energy, see
@ref Cantera.Mu0Poly or [described
here](https://cantera.org/documentation/dev/doxygen/html/d4/d9e/classCantera_1_1Mu0Poly.html#details).

Additional fields of a `piecewise-Gibbs` entry are:

<b>`h0`</b>

The molar enthalpy at the reference temperature of 298.15 K.
Defaults to 0.0.

<b>`dimensionless`</b>

A boolean flag indicating whether the values of the Gibbs free
energy are given in a dimensionless form, that is, divided by @f$RT@f$.
Defaults to `false`.

<b>`data`</b>

A mapping of temperatures to values of the Gibbs free energy. The
Gibbs free energy can be either in molar units (if `dimensionless`
is `false`) or nondimensionalized by the corresponding temperature
(if `dimensionless` is `true`). A value must be provided at
@f$T^\circ = 298.15@f$ K.

<b>`T-min`</b>

The minimum temperature at which this thermo data should be used.
Defaults to 0.0.

<b>`T-max`</b>

The maximum temperature at which this thermo data should be used.
Defaults to infinity.

**Example:**

``` yaml
thermo:
  model: piecewise-Gibbs
  h0: -230.015 kJ/mol
  dimensionless: true
  data: {298.15: -91.50963, 333.15: -85.0}
```

## Shomate {#sec-yaml-shomate}

The Shomate polynomial form given for one (@ref Cantera.ShomatePoly) or two temperature
regions (@ref Cantera.ShomatePoly2); also [described
here](https://cantera.org/science/species-thermo.html#the-shomate-parameterization).

Additional fields of a `Shomate` thermo entry are:

<b>`temperature-ranges`</b>

A list giving the temperature intervals on which the polynomials are
valid. For one temperature region, this list contains the minimum
and maximum temperatures for the polynomial. For two temperature
regions, this list contains the minimum, intermediate, and maximum
temperatures.

<b>`data`</b>

A list with one item per temperature region, where that item is a 7
item list of polynomial coefficients. The temperature regions are
arranged in ascending order.

**Example:**

``` yaml
thermo:
  model: Shomate
  temperature-ranges: [298, 1300, 6000]
  reference-pressure: 1 bar
  data:
  - [25.56759, 6.096130, 4.054656, -2.671301, 0.131021,
    -118.0089, 227.3665]
  - [35.15070, 1.300095, -0.205921, 0.013550, -3.282780,
    -127.8375, 231.7120]
```
