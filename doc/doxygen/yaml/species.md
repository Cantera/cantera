# Species Entries {#sec-yaml-species}

[TOC]

## Fields

The fields of a `species` entry are:

`name`

String identifier used for the species. Required.

`composition`

Mapping that specifies the elemental composition of the species, for
example, `{C: 1, H: 4}`. Required.

`thermo`

Mapping containing the reference state thermodynamic model
specification and parameters. See <tt>@ref sec-yaml-species-thermo</tt>.

`equation-of-state`

A mapping or list of mappings. Each mapping contains an equation of
state model specification for the species, any parameters for that
model, and any parameters for interactions with other species. See
@ref sec-yaml-species-eos-models. If this field
is absent and a model is required, the `ideal-gas` model is assumed.
See <tt>@ref sec-yaml-species-eos</tt>.

`critical-parameters`

Mapping containing parameters related to the critical state of a
species. Used in models that incorporate "real gas" effects, such
as @ref sec-yaml-eos-redlich-kwong. See <tt>@ref sec-yaml-species-crit-props</tt>.

`transport`

Mapping containing the species transport model specification and
parameters. See @ref sec-yaml-species-transport.

`sites`

The number of sites occupied by a surface or edge species. Default is 1.

`Debye-Huckel`

Additional model parameters used in the Debye-Hückel model. See
@ref sec-yaml-Debye-Huckel for more information.

### thermo {#sec-yaml-species-thermo}

Fields of a species `thermo` entry of @ref sec-yaml-species used by all models are:

`model`

String specifying the model to be used. Required. Supported model strings are:

-   <tt>@ref sec-yaml-constant-cp</tt>
-   <tt>@ref sec-yaml-nasa7</tt>
-   <tt>@ref sec-yaml-nasa9</tt>
-   <tt>@ref sec-yaml-piecewise-gibbs</tt>
-   <tt>@ref sec-yaml-shomate</tt>

`reference-pressure`

The reference pressure at which the given thermodynamic properties apply.
Defaults to 1 atm.

Additional fields are specific to the species thermodynamics model.

### equation-of-state {#sec-yaml-species-eos}

An `equation-of-state` entry of a @ref sec-yaml-species is identified as:

`model`

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

`critical-temperature`

The critical temperature of the species \[K\]

`critical-pressure`

The critical pressure of the species \[Pa\]

`acentric-factor`

Pitzer's acentric factor $omega$ \[-\]

### transport {#sec-yaml-species-transport}

`model`

String specifying the model type. The only model that is
specifically handled is <tt>@ref sec-yaml-species-transport-gas</tt>.
