# Transport Properties

Here, we describe how Cantera uses species and phase information to calculate transport
properties and rates. Similar to Cantera's approach to
[thermodynamic properties](/reference/thermo/index), transport property calculations in
Cantera depend on information at both the species and phase levels. The user must
specify transport models for both levels, and these selections must be compatible with
one another.

- The user must specify a transport model for each species and provide inputs that
  inform how species properties are calculated. For example, the user provides inputs
  that allow Cantera to calculate species collision integrals based on species-specific
  Lennard-Jones parameters.
- The user also selects a phase transport model. This model describes how the species
  interact with one another to determine properties such as viscosity, thermal
  conductivity, and diffusion coefficients for a given thermodynamic state.

## Species Transport Parameters

Transport property models in general require parameters that express the effect of each
species on the transport properties of the phase. Currently, most transport models
available in Cantera are applicable to gaseous phases.

```{admonition} YAML Usage
:class: tip
Gas transport properties can be defined in the YAML format using the
[`transport`](sec-yaml-species-transport) field of a `species` entry.
```

(sec-phase-transport-models)=
## Phase Transport Models

Multicomponent
: A multicomponent transport model for ideal gases, based on the model described by
  Dixon-Lewis {cite:t}`dixon-lewis1968`; See also Kee et al. {cite:t}`kee2017`. The
  multicomponent transport model can be specified in the YAML format by setting the
  [`transport`](sec-yaml-phase-transport) field of the phase entry to `multicomponent`.
  Implemented by class {ct}`MultiTransport`.

Mixture-averaged
: A mixture-averaged transport model for ideal gases, as described in Kee et al.
  {cite:t}`kee2017`. The mixture-averaged transport model can be specified in the YAML
  format by setting the [`transport`](sec-yaml-phase-transport) field of the phase entry
  to `mixture-averaged`. Implemented by class {ct}`MixTransport`.

High-pressure Gas
: A model for high-pressure gas transport properties based on a method of corresponding
  states {cite:p}`takahashi1975,poling2001`. The high-pressure gas transport model can
  be specified in the YAML format by setting the [`transport`](sec-yaml-phase-transport)
  field of the phase entry to `high-pressure`. Implemented by
  class {ct}`HighPressureGasTransport`.

Ionized Gas
: A model implementing the Stockmayer-(n,6,4) model for transport of ions in a gas. The
  ionized gas transport model can be specified in the YAML format by setting the
  [`transport`](sec-yaml-phase-transport) field of the phase entry to `ionized-gas`.
  Implemented by class {ct}`IonGasTransport`.

Unity Lewis Number
: A transport model for ideal gases where diffusion coefficients for all species are set
  so that the [Lewis number](https://en.wikipedia.org/wiki/Lewis_number) is 1. The unity
  Lewis number transport model can be specified in the YAML format by setting the
  [`transport`](sec-yaml-phase-transport) field of the phase entry to
  `unity-Lewis-number`. Implemented by class {ct}`UnityLewisTransport`.

Water
: A transport model for pure water applicable in both liquid and vapor phases. The water
  transport model can be specified in the YAML format by setting the
  [`transport`](sec-yaml-phase-transport) field of the phase entry to `water`.
  Implemented by class {ct}`WaterTransport`.
