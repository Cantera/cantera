@page sec-yaml-species-transport-models Species Transport Models

Species transport models are specified using the `transport` field of
@ref sec-yaml-species.

Species transport properties are a rare exception to %Cantera's use of
SI units, and use the units in which these properties are customarily
reported. No conversions are supported.

# gas {#sec-yaml-species-transport-gas}

The fields of a `gas` transport entry are:

<b>`geometry`</b>

A string specifying the geometry of the molecule. One of `atom`,
`linear`, or `nonlinear`.

<b>`diameter`</b>

The Lennard-Jones collision diameter \[Å\]

<b>`well-depth`</b>

The Lennard-Jones well depth \[K\]

<b>`dipole`</b>

The permanent dipole moment \[Debye\]. Default 0.0.

<b>`polarizability`</b>

The dipole polarizability \[Å\^3\]. Default 0.0.

<b>`rotational-relaxation`</b>

The rotational relaxation collision number at 298 K \[-\]. Default 0.0.

<b>`acentric-factor`</b>

Pitzer's acentric factor \[-\]. Default 0.0. This value may also be
specified as part of the `critical-parameters` field,
in which case the value provided there supersedes this one.

<b>`dispersion-coefficient`</b>

The dispersion coefficient, normalized by @f$e^2@f$ \[Å\^5\]. Default 0.0.

<b>`quadrupole-polarizability`</b>

The quadrupole polarizability \[Å\^5\]. Default 0.0.

**Example:**

```yaml
transport:
  model: gas
  geometry: linear
  well-depth: 107.4
  diameter: 3.458
  polarizability: 1.6
  rotational-relaxation: 3.8
```
