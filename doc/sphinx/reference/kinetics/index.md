# Chemical Reactions

Calculation of reaction rates in Cantera is done in two steps. First, a [*forward rate
constant*](https://en.wikipedia.org/wiki/Reaction_rate_constant) is calculated, which is
typically a function of temperature, but can also depend on other properties of the
mixture state such as pressure and composition. This rate constant is then used along
with the reactant and product concentrations to determine the forward and reverse rates
of reactions, as well as other quantities such as the net rate of production for each
species.

[](reaction-rates)
: This page describes how forward and reverse reaction rates and other quantities are
  calculated by Cantera.

[](rate-constants)
: This page describes the different parameterizations available in Cantera for
  calculating forward rate constants.


## Plasma Reactions

For electron-impact processes in plasma simulations, forward rate constants
calculations are done in a more complex way. The following
sections provide information to the interested user:

- [Electron-collision rates of progress](sec-electron-collision-rates-of-progress):
  How electron-impact reactions contribute to forward and reverse
  reaction rates, and how collision-only processes differ from
  kinetic reactions.
- [Electron-collision rate coefficients](sec-electron-collision-plasma-rate):
  Calculation of rate coefficients from the EEDF and tabulated
  cross sections, including superelastic reactions.
- [Electron-collision YAML format](sec-yaml-electron-collisions):
  Definition of collision-only processes and kinetic reactions
  using the root-level `electron-collisions` section.
- [Plasma phase configuration](sec-yaml-plasma):
  Configuration of the EEDF model and cross-section extrapolation.


```{toctree}
:hidden:

reaction-rates
rate-constants
```
