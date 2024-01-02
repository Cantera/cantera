# Species Thermodynamic Models

The phase models discussed in the [](phase-thermo) section implement specific models for
the thermodynamic properties appropriate for the type of phase or interface they
represent. Although each one may use different expressions to compute the properties,
they all require thermodynamic property information for the individual species.

Generally, the phase models require a parameterization of the standard state heat
capacity, enthalpy, and entropy for each species at a fixed pressure $p^\circ$ as a
function of $T$. In addition, phase models may require information describing how each
species affects the equation of state, either in terms of the species standard molar
volume, or through parameters that are specific to the equation of state.

(sec-standard-state-species-thermo)=
## Models for Species Standard State Enthalpy, Entropy, and Heat Capacity

Many of Cantera's phase thermodynamic models are formulated to make use of the
[standard state](https://goldbook.iupac.org/terms/view/S05925) thermodynamic properties
for individual species, defined at a standard pressure $p^\circ$ and for the composition
specified by the phase model. For example, this could include a pure gas in the case of
the ideal gas model, or an ion at infinite dilution in water in the case of aqueous
solutions.

```{caution}
In some parts of the Cantera documentation, properties calculated at the standard
pressure $p^\circ$ are referred to as *reference-state* thermodynamic properties, while
properties calculated using the composition defining the standard state but at any
pressure are referred to as *standard state* properties. This nomenclature is fairly
unique to Cantera, based on the desire to distinguish these different steps in the
calculation of the full thermodynamic properties, and is not often seen in other
descriptions of solution thermodynamics.
```

The necessary properties are:

1. $\hat{c}^\circ_p(T)$: the standard-state molar heat capacity at constant pressure as
   a function of temperature;
2. $\hat{h}^\circ(T)$, the standard-state molar enthalpy as a function of temperature;
3. $\hat{s}^\circ(T)$: the standard-state molar entropy as a function of temperature.

While each of these functions is implemented explicitly in Cantera for computational
efficiency, $\hat{h}^\circ(T)$ and $\hat{s}^\circ(T)$ can be expressed in terms of
$\hat{c}^\circ_p(T)$ using the relations

$$  \hat{h}^\circ(T) = \hat{h}^\circ(T_\mathrm{ref}) +
                       \int_{T_\mathrm{ref}}^T  \hat{c}^\circ_p(T) \; dT $$

and

$$  \hat{s}^\circ(T) = \hat{s}^\circ(T_\mathrm{ref}) +
                       \int_{T_\mathrm{ref}}^T  \frac{\hat{c}^\circ_p(T)}{T} \; dT $$

respectively. This means that a parameterization of $\hat{c}_p^\circ(T)$ plus the
constants $\hat{h}^\circ(T_\mathrm{ref})$ and $\hat{s}^\circ(T_\mathrm{ref})$ at a
reference temperature $T_\mathrm{ref}$ is sufficient to define the standard state
properties for a species.

The models described in this section can be used to provide standard state thermodynamic
data for each species in a phase. They are implemented by classes deriving from
{ct}`SpeciesThermoInterpType`.

```{note}
There is no requirement that all species in a phase use the same parameterization; each
species can use the one most appropriate to represent how the heat capacity depends on
temperature.
```

(sec-thermo-nasa7)=
### The NASA 7-Coefficient Polynomial Parameterization

The NASA 7-coefficient polynomial parameterization is used to compute the species
standard-state thermodynamic properties $\hat{c}^\circ_p(T)$, $\hat{h}^\circ(T)$, and
$\hat{s}^\circ(T)$.

The NASA parameterization represents $\hat{c}^\circ_p(T)$ with a fourth-order
polynomial:

$$
\frac{\hat{c}_p^\circ(T)}{R} &= a_0 + a_1 T + a_2 T^2 + a_3 T^3 + a_4 T^4

\frac{\hat{h}^\circ (T)}{R T} &= a_0 + \frac{a_1}{2} T + \frac{a_2}{3} T^2 +
    \frac{a_3}{4} T^3 + \frac{a_4}{5} T^4 + \frac{a_5}{T}

\frac{\hat{s}^\circ(T)}{R} &= a_0 \ln T + a_1 T + \frac{a_2}{2} T^2 +
    \frac{a_3}{3} T^3 + \frac{a_4}{4} T^4 + a_6
$$

This model is implemented by the C++ classes {ct}`NasaPoly1` and {ct}`NasaPoly2`.

:::{note}
This is sometimes referred to as the *NASA7 polynomial* within Cantera. It was used in
the original NASA equilibrium program and in Chemkin, which uses 7 coefficients in each
of two temperature regions. It is not compatible with the form used in more recent
versions of the NASA equilibrium program, which uses 9 coefficients for each temperature
region.
:::

:::{admonition} YAML Usage
:class: tip
A NASA7 parameterization can be defined in the YAML format by specifying
[`NASA7`](sec-yaml-nasa7) as the `model` in the species `thermo` field.
:::

(sec-thermo-nasa9)=
### The NASA 9-Coefficient Polynomial Parameterization

The NASA 9-coefficient polynomial parameterization {cite:p}`mcbride2002`, often called
*NASA9* within Cantera for short, is an extension of the NASA 7-coefficient polynomial
parameterization that includes two additional terms in each temperature region and
supports an arbitrary number of temperature regions.

The NASA9 parameterization represents the species thermodynamic properties with the
following equations:

$$
\frac{\hat{c}_p^\circ(T)}{R} &= a_0 T^{-2} + a_1 T^{-1} + a_2 + a_3 T
    + a_4 T^2 + a_5 T^3 + a_6 T^4

\frac{\hat{h}^\circ(T)}{R T} &= - a_0 T^{-2} + a_1 \frac{\ln T}{T} + a_2
    + \frac{a_3}{2} T + \frac{a_4}{3} T^2  + \frac{a_5}{4} T^3 +
    \frac{a_6}{5} T^4 + \frac{a_7}{T}

\frac{\hat{s}^\circ(T)}{R} &= - \frac{a_0}{2} T^{-2} - a_1 T^{-1} + a_2 \ln T
   + a_3 T + \frac{a_4}{2} T^2 + \frac{a_5}{3} T^3  + \frac{a_6}{4} T^4 + a_8
$$

A common source for species data in the NASA9 format is the
[NASA ThermoBuild](/userguide/thermobuild) tool. This model is implemented by the C++
classes {ct}`Nasa9Poly1` and {ct}`Nasa9PolyMultiTempRegion`.

:::{admonition} YAML Usage
:class: tip
A NASA-9 parameterization can be defined in the YAML format by specifying
[`NASA9`](sec-yaml-nasa9) as the `model` in the species `thermo` field.
:::

(sec-thermo-shomate)=
### The Shomate Parameterization

The Shomate parameterization is:

$$
\hat{c}_p^\circ(T) &= A + Bt + Ct^2 + Dt^3 + \frac{E}{t^2}

\hat{h}^\circ(T) &= At + \frac{Bt^2}{2} + \frac{Ct^3}{3} + \frac{Dt^4}{4} -
    \frac{E}{t} + F

\hat{s}^\circ(T) &= A \ln t + B t + \frac{Ct^2}{2} + \frac{Dt^3}{3} - \frac{E}{2t^2} + G
$$

where $t = T / 1000\textrm{ K}$. It requires 7 coefficients $A$, $B$, $C$, $D$, $E$,
$F$, and $G$. This parameterization is used to represent standard-state properties in
the [NIST Chemistry WebBook](http://webbook.nist.gov/chemistry). The values of the
coefficients $A$ through $G$ should be entered precisely as shown there, with no units
attached. Unit conversions to SI is handled internally. This model is implemented by the
C++ classes {ct}`ShomatePoly` and {ct}`ShomatePoly2`.

:::{admonition} YAML Usage
:class: tip
A Shomate parameterization can be defined in the YAML format by specifying
[`Shomate`](sec-yaml-shomate) as the `model` in the species `thermo` field.
:::

(sec-thermo-const-cp)=
### Constant Heat Capacity

In some cases, species properties may only be required at a single temperature or over a
narrow temperature range. In such cases, the heat capacity can be approximated as
constant, and simple expressions can be used for the thermodynamic properties:

$$
\hat{c}_p^\circ(T) &= \hat{c}_p^\circ(T_\mathrm{ref})

\hat{h}^\circ(T) &= \hat{h}^\circ\left(T_\mathrm{ref}\right) + \hat{c}_p^\circ \left(T-T_\mathrm{ref}\right)

\hat{s}^\circ(T) &= \hat{s}^\circ(T_\mathrm{ref}) + \hat{c}_p^\circ \ln{\left(\frac{T}{T_\mathrm{ref}}\right)}
$$

The parameterization uses four constants: $T_\mathrm{ref}$,
$\hat{c}_p^\circ(T_\mathrm{ref})$, $\hat{h}^\circ(T_\mathrm{ref})$, and
$\hat{s}^\circ(T)$. The default value of $T_\mathrm{ref}$ is 298.15 K; the default value
for the other parameters is 0.0. This model is implemented by the C++ class
{ct}`ConstCpPoly`.

:::{admonition} YAML Usage
:class: tip
A constant heat capacity parameterization can be defined in the YAML format by
specifying [`constant-cp`](sec-yaml-constcp) as the `model` in the species `thermo`
field.
:::

(sec-thermo-piecewise-Gibbs)=
### Piecewise Gibbs Free Energy

This parameterization uses a piecewise interpolation of the Gibbs free energy between specified values, with a constant heat capacity between each pair of
interpolation points. This model is implemented by the C++ class {ct}`Mu0Poly`.

:::{admonition} YAML Usage
:class: tip
A piecewise Gibbs parameterization can be defined in the YAML format by specifying
[`piecewise-Gibbs`](sec-yaml-constcp) as the `model` in the species `thermo` field.
:::

## Models for Species Contributions to the Equation of State

Besides enthalpy, entropy, and heat capacity data, some phase models require additional
parameters that describe how each species affects the mechanical equation of state.
These models are used in combination with one of the above parameterizations for the
standard state enthalpy, entropy, and heat capacity. This applies to most models that
represent condensed phases.

### Constant Volume Model

This species equation of state model can be used for phase models that require standard
state density information for each species, and simply specifies a constant specific
volume, density or molar density for the species. Implemented by the C++ class
{ct}`PDSS_ConstVol`.

:::{admonition} YAML Usage
:class: tip
A constant volume parameterization can be defined in the YAML format by specifying
[`constant-volume`](sec-yaml-eos-constant-volume) as the `model` in the species
`equation-of-state` field.
:::

### Temperature-dependent Density

This species equation of state model can be used for phase models that require the
standard state molar volume $V_{m,k}^\circ$ for each species $k$. The standard state
density $\rho^\circ_k$ is described using a third-order polynomial:

$$  \rho^\circ_k(T) = \frac{W_k}{V_{m,k}^\circ(T)} = a_0 + a_1 T + a_2 T^2 + a_3 T^3  $$

where $W_k$ is the molecular weight of the species. This model is implemented by class
{ct}`PDSS_SSVol`.

:::{admonition} YAML Usage
:class: tip
A temperature-dependent density parameterization can be defined in the YAML format by
specifying [`density-temperature-polynomial`](sec-yaml-eos-density-temperature-polynomial)
as the `model` in the species `equation-of-state` field.
:::

### Temperature-dependent Molar Volume

This species equation of state model can be used for phase models that require the
standard state molar volume $V_{m,k}^\circ$ for each species $k$. The standard state
molar volume is described using a third-order polynomial:

$$  V_{m,k}^\circ(T) = a_0 + a_1 T + a_2 T^2 + a_3 T^3  $$

This model is implemented by class {ct}`PDSS_SSVol`.

:::{admonition} YAML Usage
:class: tip
A temperature-dependent molar volume parameterization can be defined in the YAML format
by specifying [`molar-volume-temperature-polynomial`](sec-yaml-eos-molar-volume-temperature-polynomial)
as the `model` in the species `equation-of-state` field.
:::

### Peng-Robinson Equation of State

This species equation of state model is parameterized using the coefficients $a$, $b$,
and $\omega$ for a pure species that follows the Peng-Robinson equation of state:

$$  P = \frac{RT}{V_m - b} - \frac{a\alpha}{V_m^2 + 2bV_m - b^2}  $$

where $V_m$ is the molar volume,

$$  \alpha = \left[ 1 + \kappa \left(1 - \sqrt{T_r}\right) \right]^2  $$

$$ \kappa =
\begin{cases}
0.37464 + 1.54226\omega - 0.26992\omega^2,  & \omega \le 0.491 \\
0.379642 + 1.487503\omega - 0.164423\omega^2 + 0.016667\omega^3 , & \omega > 0.491
\end{cases}
$$

$T_r = T / T_c$ is the reduced temperature, and $T_c$ is the critical temperature.

These pure-species properties are combined in the multi-species Peng-Robinson phase
model, implemented by class {ct}`PengRobinson`.

:::{admonition} YAML Usage
:class: tip
A Peng-Robinson parameterization for a species can be defined in the YAML format by
specifying [`Peng-Robinson`](sec-yaml-eos-peng-robinson) as the `model` in the species
`equation-of-state` field.
:::

### Redlich-Kwong Equation of State

This species equation of state model is parameterized using the coefficients $a$ and $b$
for a pure species that follows the Redlich-Kwong equation of state:

$$  P = \frac{RT}{V_m-b} - \frac{a}{\sqrt{T} V_m (V_m + b) }  $$

where $V_m$ is the molar volume.

These pure-species properties are combined in the multi-species Redlich-Kwong phase
model, implemented by class {ct}`RedlichKwongMFTP`.

:::{admonition} YAML Usage
:class: tip
A Redlich-Kwong parameterization for a species can be defined in the YAML format by
specifying [`Redlich-Kwong`](sec-yaml-eos-redlich-kwong) as the `model` in the species
`equation-of-state` field.
:::

## Other Species Equation of State Models

The following models provide complete standard state parameterizations for a species,
including molar volume, enthalpy, entropy, and heat capacity.

### Helgeson-Kirkham-Flowers-Tanger Model

The Helgeson-Kirkham-Flowers-Tanger model for aqueous species. Implemented by the C++
class {ct}`PDSS_HKFT`.

:::{admonition} YAML Usage
:class: tip
An HKFT parameterization can be defined in the YAML format by specifying
[`HKFT`](sec-yaml-eos-hkft) as the `model` in the species `equation-of-state` field.
:::

### Liquid Water

A model for liquid water that uses the IAPWS95 formulation {cite:p}`wagner2002`.
Implemented by the C++ class {ct}`PDSS_Water`.

:::{admonition} YAML Usage
:class: tip
A liquid water parameterization can be defined in the YAML format by specifying
[`liquid-water-IAPWS95`](sec-yaml-eos-liquid-water-iapws95) as the `model` in the
species `equation-of-state` field.
:::
