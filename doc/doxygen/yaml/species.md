# Species Entries {#sec-yaml-species}

[TOC]

## Overview

A species entry in %Cantera is used to specify the name, composition,
thermodynamic, and transport properties of an individual species.

The default location for species entries is in the `species` section of
the input file. Species defined in this section will automatically be
considered for addition to phases defined in the same file. Species can
be defined in other sections of the input file or in other input files,
and these species definitions can be used in phase definitions by
explicitly referencing the section name.

### Species Name

The name of a species is given in the `name` field of a `species` entry.
Names may include almost all printable characters, with the exception of
spaces. The use of some characters such as `[`, `]`, and `,` may require
that species names be enclosed in quotes when written in YAML. Some
valid species names given in a YAML list include:

``` yaml
[CH4, methane, argon_2+, "C[CH2]", CH2(singlet), "H2O,l"]
```

### Elemental Composition

The elemental composition of a species is specified as a mapping in the
`composition` entry.

For gaseous species, the elemental composition is well-defined, since
the species represent distinct molecules. For species in solid or liquid
solutions, or on surfaces, there may be several possible ways of
defining the species. For example, an aqueous species might be defined
with or without including the water molecules in the solvation cage
surrounding it.

For surface species, it is possible for the `composition` mapping to be
empty, in which case the species is composed of nothing, and represents
an empty surface site. This can also be done to represent vacancies in
solids. A charged vacancy can be defined to be composed solely of
electrons.

The number of atoms of an element must be non-negative, except for the
special \"element\" `E` that represents an electron.

**Examples:**

``` yaml
composition: {C: 1, O: 2}  # carbon dioxide
composition: {Ar: 1, E: -2}  # Ar++
composition: {Y: 1, Ba: 2, Cu: 3, O: 6.5}  # stoichiometric YBCO
composition: {}  # A surface species representing an empty site
```

### Thermodynamic Properties {#sec-yaml-species-thermo-properties}

In addition to the thermodynamic model used at the phase level for
computing properties, parameterizations are usually required for the
enthalpy, entropy, and specific heat capacities of individual species
under standard conditions. These parameterizations are provided in the
`thermo` field of each `species` entry.

The parameterization used to provide this information is specified by
the `model` field of the `thermo` field. The models available are:

-   <tt>@ref sec-yaml-constant-cp</tt>: Constant heat capacity
-   <tt>@ref sec-yaml-nasa7</tt>: 7-coefficient NASA polynomials in one or two
    temperature regions
-   <tt>@ref sec-yaml-nasa9</tt>: 9-coefficient NASA polynomials in one or more
    temperature regions
-   <tt>@ref sec-yaml-piecewise-gibbs</tt>: Interpolation between tabulated Gibbs free
    energies using a constant heat capacity in each temperature interval
-   <tt>@ref sec-yaml-shomate</tt>: Shomate polynomials in one or two temperature
    regions

For descriptions of fields used by each model as well as examples,
see @subpage sec-yaml-species-thermo-models.

### Equation of State {#sec-yaml-species-equation-of-state}

For some phase thermodynamic models, additional equation of state
parameterizations are needed for each species. This information is
provided in the `equation-of-state` field of each `species` entry, with
the type of parameterization used specified by the `model` field of the
`equation-of-state` field. The models available are:

-   <tt>@ref sec-yaml-eos-constant-volume</tt>: A fixed value of mass density, molar
    density, or molar volume
-   <tt>@ref sec-yaml-eos-density-temperature-polynomial</tt>: Mass density
    parameterized using a cubic polynomial in temperature
-   <tt>@ref sec-yaml-eos-hkft</tt>: The Helgeson-Kirkham-Flowers-Tanger model for
    aqueous species
-   <tt>@ref sec-yaml-eos-ideal-gas</tt>: A species following the ideal gas law
-   <tt>@ref sec-yaml-eos-ions-from-neutral</tt>: Used with the
    ions-from-neutral-molecule phase model
-   <tt>@ref sec-yaml-eos-liquid-water-iapws95</tt>: The IAPWS95 equation of state for
    water, applied only in the liquid region
-   <tt>@ref sec-yaml-eos-molar-volume-temperature-polynomial</tt>: Molar volume
    parameterized using a cubic polynomial in temperature
-   <tt>@ref sec-yaml-eos-peng-robinson</tt>: A species which follows the Peng-Robinson
    equation of state
-   <tt>@ref sec-yaml-eos-redlich-kwong</tt>: A species which follows the Redlich-Kwong
    equation of state

For descriptions of fields used by each model as well as examples,
see @subpage sec-yaml-species-eos-models.

### Transport Coefficients {#sec-yaml-species-transport-coefficients}

Transport-related parameters for each species are needed in order to
calculate transport properties of a phase. These parameters are provided
in the `transport` field of each `species` entry, with the type of the
parameterization used specified by the `model` field of the `transport`
field. The parameters used depend on the transport model specified at the phase level.

Currently, there is a single model type that is specifically handled:

-   <tt>@ref sec-yaml-species-transport-gas</tt>: Gas phase transport

For descriptions of fields used by the model as well as an example,
see @subpage sec-yaml-species-transport-models.

### Coverage Dependencies {#sec-yaml-species-coverage}

The coverage-dependent surface species formulation calculates coverage-dependent
correction factors to the ideal surface phase properties. Used in conjunction with the
[coverage-dependent surface phase model](@ref sec-yaml-coverage-dependent-surface).
Supported algebraic models are:

-   <tt>@ref sec-yaml-species-coverage-linear</tt>: Linear dependency model
-   <tt>@ref sec-yaml-species-coverage-polynomial</tt>: Polynomial dependency model
-   <tt>@ref sec-yaml-species-coverage-piecewise-linear</tt>:
    Piecewise-linear dependency model
-   <tt>@ref sec-yaml-species-coverage-interpolative</tt>:
    Interpolative dependency model

For descriptions of fields used by each model as well as examples,
see @subpage sec-yaml-species-coverage-models.

## Species API Reference

The fields of a `species` entry are:

<b>`name`</b>

String identifier used for the species. Required.

<b>`composition`</b>

Mapping that specifies the elemental composition of the species, for
example, `{C: 1, H: 4}`. Required.

<b>`thermo`</b>

Mapping containing the reference state thermodynamic model
specification and parameters.

-   <b>`model`</b>: String specifying the model to be used. Required.
    See @ref sec-yaml-species-thermo-properties for supported models.

-   <b>`reference-pressure`</b>: The reference pressure at which the given
    thermodynamic properties apply. Defaults to 1 atm.

-   Additional fields are specific to the species thermodynamics model.
    (see @ref sec-yaml-species-thermo-models).

<b>`equation-of-state`</b>

A mapping or list of mappings. Each mapping contains an equation of
state model specification for the species, any parameters for that
model, and any parameters for interactions with other species. See
@ref sec-yaml-species-eos-models. If this field
is absent and a model is required, the `ideal-gas` model is assumed.

-   <b>`model`</b>: String specifying the model to be used. Required.
    See @ref sec-yaml-species-equation-of-state for supported models.

-   Additional fields are specific to the species equation of state model
    (see @ref sec-yaml-species-eos-models).

<b>`critical-parameters`</b>

Mapping containing parameters related to the critical state of a
species. Used in models that incorporate "real gas" effects, such
as @ref sec-yaml-eos-redlich-kwong.

-   <b>`critical-temperature`</b>: The critical temperature of the species \[K\]

-   <b>`critical-pressure`</b>: The critical pressure of the species \[Pa\]

-   <b>`acentric-factor`</b>: Pitzer's acentric factor @f$\omega@f$ \[-\]

<b>`transport`</b>

Mapping containing the species transport model specification and
parameters.

-   <b>`model`</b>: String specifying the model type. See
    @ref sec-yaml-species-transport-coefficients and
    @ref sec-yaml-species-transport-models.

<b>`sites`</b>

The number of sites occupied by a surface or edge species. Default is 1.

<b>`Debye-Huckel`</b>

Additional model parameters used in the Debye-Hückel model. See
@ref sec-yaml-Debye-Huckel for more information.

<b>`coverage-dependencies`</b>

Mappings where keys are the name of species whose coverage affects thermodynamic
properties of the node-owner species, see @ref Cantera.CoverageDependentSurfPhase.

-   <b>`model`</b>: String specifying the model to be used. Required.
    See @ref sec-yaml-species-coverage for supported models.

-   Additional map entries that correspond to an individual dependency between the
    node-owner species and keyed species (see @ref sec-yaml-species-coverage-models).
