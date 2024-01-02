# One-dimensional Flames

Cantera includes a set of models for representing steady-state, quasi-one-dimensional
reacting flows.

These models can be used to simulate a number of common flames, such as:

- freely-propagating premixed laminar flames
- burner-stabilized premixed flames
- counterflow diffusion flames
- counterflow (strained) premixed flames

Additional capabilities include simulation of surface reactions, which can be used to
represent processes such as combustion on a catalytic surface or chemical vapor
deposition processes.

All of these configurations are simulated using a common set of governing equations
within a 1D flow domain, with the differences between the models being represented by
differences in the boundary conditions applied.

[](governing-equations)
: This page describes the governing equations and the various boundary conditions that
  can be applied.

[](discretization)
: This page describes the finite difference schemes used to discretize the 1D governing
  equations and the criteria used for refining the grid.

[](nonlinear-solver)
: This page describes the hybrid time-stepping--steady-state damped Newton solver that
  is used to solve the discretized governing equations.

```{toctree}
:hidden:

governing-equations
discretization
nonlinear-solver
```
