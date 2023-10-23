# Species Thermodynamic Models

The phase models discussed in the [](phase-thermo) section implement specific models for
the thermodynamic properties appropriate for the type of phase or interface they
represent. Although each one may use different expressions to compute the properties,
they all require thermodynamic property information for the individual species. For the
phase types implemented at present, the properties needed are:

1. the molar heat capacity at constant pressure $\hat{c}^\circ_p(T)$ for a range of
   temperatures and a reference pressure $p^\circ$;
2. the molar enthalpy $\hat{h}(T^\circ, p^\circ)$ at $p^\circ$ and a reference
   temperature $T^\circ$;
3. the absolute molar entropy $\hat{s}(T^\circ, p^\circ)$ at $(T^\circ, p^\circ)$.

The superscript $^\circ$ here represents the *reference state*--a specified state (that
is, a set of conditions $T^\circ$ and $p^\circ$ and fixed chemical composition) at which
thermodynamic properties are known.

The models described in this section can be used to provide thermodynamic data for each
species in a phase. Each model implements a different *parameterization* (functional
form) for the heat capacity. Note that there is no requirement that all species in a
phase use the same parameterization; each species can use the one most appropriate to
represent how the heat capacity depends on temperature.

Currently, several types are implemented that provide species properties appropriate for
models of ideal gas mixtures, ideal solutions, and pure compounds.

(sec-thermo-nasa7)=
## The NASA 7-Coefficient Polynomial Parameterization

The NASA 7-coefficient polynomial parameterization is used to compute the species
reference-state thermodynamic properties $\hat{c}^\circ_p(T)$, $\hat{h}^\circ(T)$, and
$\hat{s}^\circ(T)$.

The NASA parameterization represents $\hat{c}^\circ_p(T)$ with a fourth-order
polynomial:

$$
\frac{\hat{c}_p^\circ(T)}{\overline{R}} &= a_0 + a_1 T + a_2 T^2 + a_3 T^3 + a_4 T^4

\frac{\hat{h}^\circ (T)}{\overline{R} T} &= a_0 + \frac{a_1}{2} T + \frac{a_2}{3} T^2 +
    \frac{a_3}{4} T^3 + \frac{a_4}{5} T^4 + \frac{a_5}{T}

\frac{\hat{s}^\circ(T)}{\overline{R}} &= a_0 \ln T + a_1 T + \frac{a_2}{2} T^2 +
    \frac{a_3}{3} T^3 + \frac{a_4}{4} T^4 + a_6
$$

:::{note}
This is the "old" NASA polynomial form, used in the original NASA equilibrium program
and in Chemkin, which uses 7 coefficients in each of two temperature regions. It is not compatible with the form used in more recent versions of the NASA equilibrium program,
which uses 9 coefficients for each temperature region.
:::

:::{admonition} YAML Usage
:class: tip
A NASA-7 parameterization can be defined in the YAML format by specifying
[`NASA7`](sec-yaml-nasa7) as the `model` in the species `thermo` field.
:::

(sec-thermo-nasa9)=
## The NASA 9-Coefficient Polynomial Parameterization

The NASA 9-coefficient polynomial parameterization {cite:p}`mcbride2002` ("NASA9" for
short) is an extension of the NASA 7-coefficient polynomial parameterization which
includes two additional terms in each temperature region, as well as supporting an
arbitrary number of temperature regions.

The NASA9 parameterization represents the species thermodynamic properties with the
following equations:

$$
\frac{\hat{c}_p^\circ(T)}{\overline{R}} &= a_0 T^{-2} + a_1 T^{-1} + a_2 + a_3 T
    + a_4 T^2 + a_5 T^3 + a_6 T^4

\frac{\hat{h}^\circ(T)}{\overline{R} T} &= - a_0 T^{-2} + a_1 \frac{\ln T}{T} + a_2
    + \frac{a_3}{2} T + \frac{a_4}{3} T^2  + \frac{a_5}{4} T^3 +
    \frac{a_6}{5} T^4 + \frac{a_7}{T}

\frac{\hat{s}^\circ(T)}{\overline{R}} &= - \frac{a_0}{2} T^{-2} - a_1 T^{-1} + a_2 \ln T
   + a_3 T + \frac{a_4}{2} T^2 + \frac{a_5}{3} T^3  + \frac{a_6}{4} T^4 + a_8
$$

A common source for species data in the NASA9 format is the
[NASA ThermoBuild](/userguide/thermobuild) tool.

:::{admonition} YAML Usage
:class: tip
A NASA-9 parameterization can be defined in the YAML format by specifying
[`NASA9`](sec-yaml-nasa9) as the `model` in the species `thermo` field.
:::

(sec-thermo-shomate)=
## The Shomate Parameterization

The Shomate parameterization is:

$$
\hat{c}_p^\circ(T) &= A + Bt + Ct^2 + Dt^3 + \frac{E}{t^2}

\hat{h}^\circ(T) &= At + \frac{Bt^2}{2} + \frac{Ct^3}{3} + \frac{Dt^4}{4} -
    \frac{E}{t} + F

\hat{s}^\circ(T) &= A \ln t + B t + \frac{Ct^2}{2} + \frac{Dt^3}{3} - \frac{E}{2t^2} + G
$$

where $t = T / 1000\textrm{ K}$. It requires 7 coefficients $A$, $B$, $C$, $D$, $E$,
$F$, and $G$. This parameterization is used to represent reference-state properties in
the [NIST Chemistry WebBook](http://webbook.nist.gov/chemistry). The values of the
coefficients $A$ through $G$ should be entered precisely as shown there, with no units
attached. Unit conversions to SI is handled internally.

:::{admonition} YAML Usage
:class: tip
A Shomate parameterization can be defined in the YAML format by specifying
[`Shomate`](sec-yaml-shomate) as the `model` in the species `thermo` field.
:::

(sec-thermo-const-cp)=
## Constant Heat Capacity

In some cases, species properties may only be required at a single temperature or over a
narrow temperature range. In such cases, the heat capacity can be approximated as
constant, and simple expressions can be used for the thermodynamic properties:

$$
\hat{c}_p^\circ(T) &= \hat{c}_p^\circ(T^\circ)

\hat{h}^\circ(T) &= \hat{h}^\circ\left(T^\circ\right) + \hat{c}_p^\circ \left(T-T^\circ\right)

\hat{s}^\circ(T) &= \hat{s}^\circ(T^\circ) + \hat{c}_p^\circ \ln{\left(\frac{T}{T^\circ}\right)}
$$

The parameterization uses four constants: $T^\circ, \hat{c}_p^\circ(T^\circ),
\hat{h}^\circ(T^\circ)$, and $\hat{s}^\circ(T)$. The default value of $T^\circ$ is
298.15 K; the default value for the other parameters is 0.0.

:::{admonition} YAML Usage
:class: tip
A constant heat capacity parameterization can be defined in the YAML format by
specifying [`constant-cp`](sec-yaml-constcp) as the `model` in the species `thermo`
field.
:::
