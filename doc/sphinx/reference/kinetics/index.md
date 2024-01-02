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

```{toctree}
:hidden:

reaction-rates
rate-constants
```
