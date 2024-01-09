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
internal energy, entropy, and others, depend on the state variables. The following
sections describe the species and phase thermodynamic models available in Cantera.

[](species-thermo)
: The models and equations that Cantera uses to calculate species thermodynamic
  properties, such as the NASA 7-parameter polynomial form.

[](phase-thermo)
: The theory behind some of Cantera's phase models, such as the ideal gas law.

The user must specify the thermodynamic models and provide input data to be used for
both levels, and these selections must be compatible with one another. For instance: one
cannot pair certain non-ideal species thermodynamic models with an ideal phase model.

```{toctree}
:hidden:
:maxdepth: 1
:caption: Thermodynamic Models

species-thermo
phase-thermo
```

### The Intensive Thermodynamic State

Cantera's phase thermodynamic model, implemenented by the C++ {ct}`ThermoPhase` class
and classes derived from it, works only with the intensive thermodynamic state. That is,
all extensive properties (enthalpy, entropy, internal energy, volume, etc.) are computed
for a unit quantity (on a mass or mole basis). For example, there is a method
{ct}`ThermoPhase::enthalpy_mole` that returns the molar enthalpy (J/kmol), and a method
{ct}`ThermoPhase::enthalpy_mass` that returns the specific enthalpy (J/kg), but no
method `ThermoPhase::enthalpy()` that would return the total enthalpy (J). This is
because class {ct}`ThermoPhase` does not store the total amount (mass or mole) of the
phase.

The intensive state of a single-component phase in equilibrium is fully specified by the
values of any $r+1$ independent thermodynamic properties, where $r$ is the number of
reversible work modes. If the only reversible work mode is compression (a "simple
compressible substance"), then two properties suffice to specify the intensive state. By
default, class {ct}`ThermoPhase` stores internally the values of the *temperature*, the
*mass density*, and the *mass fractions* of all species. These values are sufficient to
fix the intensive thermodynamic state of the phase and to compute any other intensive
properties. This choice is arbitrary, and for most purposes you can't tell which
properties are stored and which are computed. For some phase models, other choices of
the *intrinsic state variables* are chosen, such as incompressible phases where the
pressure replaces the mass density as an independent variable.
