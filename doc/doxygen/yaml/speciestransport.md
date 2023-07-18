# Species Transport Models {#sec-yaml-species-transport-models}

## gas {#sec-yaml-species-transport-gas}

Species transport properties are a rare exception to %Cantera's use of
SI units, and use the units in which these properties are customarily
reported. No conversions are supported.

The additional fields of a `gas` transport entry are:

`geometry`

A string specifying the geometry of the molecule. One of `atom`,
`linear`, or `nonlinear`.

`diameter`

The Lennard-Jones collision diameter \[Å\]

`well-depth`

The Lennard-Jones well depth \[K\]

`dipole`

The permanent dipole moment \[Debye\]. Default 0.0.

`polarizability`

The dipole polarizability \[Å\^3\]. Default 0.0.

`rotational-relaxation`

The rotational relaxation collision number at 298 K \[-\]. Default 0.0.

`acentric-factor`

Pitzer's acentric factor \[-\]. Default 0.0. This value may also be
specified as part of the @ref sec-yaml-species-crit-props field,
in which case the value provided there supersedes this one.

`dispersion-coefficient`

The dispersion coefficient, normalized by $e^2$ \[Å\^5\]. Default 0.0.

`quadrupole-polarizability`

The quadrupole polarizability \[Å\^5\]. Default 0.0.

**Example:**

``` yaml
transport:
  model: gas
  geometry: linear
  well-depth: 107.4
  diameter: 3.458
  polarizability: 1.6
  rotational-relaxation: 3.8
```
