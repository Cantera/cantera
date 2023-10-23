# Phase Thermodynamic Models

On this page, we list the phase thermodynamic models implemented in Cantera, with
links to the documentation for their YAML input parameters and the documentation for
the C++ classes which implement these models. This API documentation may also provide
references or a mathematical description of the model.

## Ideal Gas Mixture

A mixture which follows the ideal gas law. Defined in the YAML format by specifying
[`ideal-gas`](sec-yaml-ideal-gas) in the `thermo` field. Implemented by class
{ct}`IdealGasPhase`.

## Stoichiometric Solid

A *stoichiometric solid* is one that is modeled as having a precise, fixed composition,
given by the composition of the one species present. Defined in the YAML format by
specifying [`fixed-stoichiometry`](sec-yaml-fixed-stoichiometry) in the `thermo` field.
Implemented by class {ct}`StoichSubstance`.

## Ideal Surface

An interface between two bulk phases where the species behave as an ideal solution and
the composition is described by the coverage of each species on the surface. Defined in
the YAML format by specifying [`ideal-surface`](sec-yaml-ideal-surface) in the `thermo`
field. Implemented by class {ct}`SurfPhase`.
