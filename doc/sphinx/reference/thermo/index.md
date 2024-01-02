# Thermodynamic Properties

In this section, we describe how Cantera uses species and phase information to calculate
thermodynamic properties.

Thermodynamic properties typically depend on information at both the species and phase
levels. Generally, the species thermodynamic model and accompanying coefficient data
specifies how the reference enthalpy and entropy values for each species are calculated
as a function of temperature. The phase model then describes how the species interact
with one another to determine phase properties and species specific properties for a
given thermodynamic state. This includes both the mechanical equation of state
($p$-$\hat{v}$-$T$ relationship) as well as how species-specific properties, such as
internal energy, entropy, and others, depend on the state variables.

The user must specify the thermodynamic models and provide input data to be used for
both levels, and these selections must be compatible with one another. For instance: one
cannot pair certain non-ideal species thermodynamic models with an ideal phase model.

The following sections describe the species and phase thermodynamic models available
in Cantera.

[](species-thermo)
: The models and equations that Cantera uses to calculate species thermodynamic
  properties, such as the NASA 7-parameter polynomial form.

[](phase-thermo)
: The theory behind some of Cantera's phase models, such as the ideal gas law.

```{toctree}
:hidden:
:maxdepth: 1
:caption: Thermodynamic Models

species-thermo
phase-thermo
```
