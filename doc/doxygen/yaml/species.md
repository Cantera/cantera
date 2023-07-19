# Species Entries {#sec-yaml-species}

[TOC]

## Fields

The fields of a `species` entry are:

<b>`name`</b>

String identifier used for the species. Required.

<b>`composition`</b>

Mapping that specifies the elemental composition of the species, for
example, `{C: 1, H: 4}`. Required.

<b>`thermo`</b>

Mapping containing the reference state thermodynamic model
specification and parameters. See <tt>@ref sec-yaml-species-thermo</tt>.

<b>`equation-of-state`</b>

A mapping or list of mappings. Each mapping contains an equation of
state model specification for the species, any parameters for that
model, and any parameters for interactions with other species. See
@ref sec-yaml-species-eos-models. If this field
is absent and a model is required, the `ideal-gas` model is assumed.
See <tt>@ref sec-yaml-species-eos</tt>.

<b>`critical-parameters`</b>

Mapping containing parameters related to the critical state of a
species. Used in models that incorporate "real gas" effects, such
as @ref sec-yaml-eos-redlich-kwong. See <tt>@ref sec-yaml-species-crit-props</tt>.

<b>`transport`</b>

Mapping containing the species transport model specification and
parameters. See <tt>@ref sec-yaml-species-transport</tt>.

<b>`sites`</b>

The number of sites occupied by a surface or edge species. Default is 1.

<b>`Debye-Huckel`</b>

Additional model parameters used in the Debye-Hückel model. See
@ref sec-yaml-Debye-Huckel for more information.

<b>`coverage-dependencies`</b>

Used in conjunction with the [coverage-dependent surface phase
model](@ref sec-yaml-coverage-dependent-surface). See
<tt>@ref sec-yaml-species-coverage</tt>.

### thermo {#sec-yaml-species-thermo}

Fields of a species `thermo` entry used by all models are:

<b>`model`</b>

String specifying the model to be used. Required. Supported model strings are:

-   <tt>@ref sec-yaml-constant-cp</tt>
-   <tt>@ref sec-yaml-nasa7</tt>
-   <tt>@ref sec-yaml-nasa9</tt>
-   <tt>@ref sec-yaml-piecewise-gibbs</tt>
-   <tt>@ref sec-yaml-shomate</tt>

<b>`reference-pressure`</b>

The reference pressure at which the given thermodynamic properties apply.
Defaults to 1 atm.

Additional fields are specific to the species thermodynamics model.

### equation-of-state {#sec-yaml-species-eos}

A species `equation-of-state` entry is identified by the field

<b>`model`</b>

String specifying the model to be used. Required. Supported model strings are:

-   <tt>@ref sec-yaml-eos-constant-volume</tt>
-   <tt>@ref sec-yaml-eos-density-temperature-polynomial</tt>
-   <tt>@ref sec-yaml-eos-hkft</tt>
-   <tt>@ref sec-yaml-eos-ideal-gas</tt>
-   <tt>@ref sec-yaml-eos-ions-from-neutral</tt>
-   <tt>@ref sec-yaml-eos-liquid-water-iapws95</tt>
-   <tt>@ref sec-yaml-eos-molar-volume-temperature-polynomial</tt>
-   <tt>@ref sec-yaml-eos-peng-robinson</tt>
-   <tt>@ref sec-yaml-eos-redlich-kwong</tt>

Additional fields are specific to the species equation of state model
(see @ref sec-yaml-species-eos-models).

### critical-parameters {#sec-yaml-species-crit-props}

<b>`critical-temperature`</b>

The critical temperature of the species \[K\]

<b>`critical-pressure`</b>

The critical pressure of the species \[Pa\]

<b>`acentric-factor`</b>

Pitzer's acentric factor \f$\omega\f$ \[-\]

### transport {#sec-yaml-species-transport}

<b>`model`</b>

String specifying the model type. The only model that is
specifically handled is <tt>@ref sec-yaml-species-transport-gas</tt>.

### coverage-dependencies {#sec-yaml-species-coverage}

Species `coverage-dependencies` entries are mappings where keys are the name of species
whose coverage affects thermodynamic properties of the node-owner species, see
@ref Cantera.CoverageDependentSurfPhase.

Map entries of `coverage-dependencies` are identified by the field

<b>`model`</b>

String specifying the model to be used. Required. Supported model strings are:

-   <tt>@ref sec-yaml-species-coverage-linear</tt>
-   <tt>@ref sec-yaml-species-coverage-polynomial</tt>
-   <tt>@ref sec-yaml-species-coverage-piecewise-linear</tt>
-   <tt>@ref sec-yaml-species-coverage-interpolative</tt>

Additional map entries that correspond to an individual dependency between the
node-owner species and keyed species (see @ref sec-yaml-species-coverage-models).
