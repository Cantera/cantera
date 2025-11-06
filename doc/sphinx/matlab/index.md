
# MATLAB Toolbox Documentation

```{caution}
The MATLAB toolbox is an experimental part of Cantera and may be changed without notice.
```

The _experimental_ MATLAB toolbox for Cantera is currently in preview and replaces the
_legacy_ MATLAB toolbox that was discontinued after Cantera 3.0. The
replacement introduces modern MATLAB object-oriented programming structure and syntax,
while replacing a custom-compiled MEX interface with MATLAB's
[`clibgen`](https://www.mathworks.com/help/matlab/ref/clibgen-package.html) to
automatically wrap functions exposed by [Cantera CLib](../clib/index).

```{note}
For current information on the state of the _experimental_ MATLAB toolbox, refer to
GitHub issues and enhancement requests, specifically:

- [**Open Pull Requests**](https://github.com/Cantera/cantera/pulls?q=is%3Apr+state%3Aopen+label%3Amatlab)
- [**Open Issues**](https://github.com/Cantera/cantera/issues?q=is%3Aissue%20state%3Aopen%20label%3Amatlab)
- [**Open Enhancements**](https://github.com/Cantera/enhancements/issues?q=is%3Aissue%20state%3Aopen%20label%3Amatlab)

Installation instructions are provided in
[`interfaces/matlab_experimental/README.md`](https://github.com/Cantera/cantera/blob/main/interfaces/matlab_experimental/README.md)
within the Cantera source code.
```

## [Objects Representing Phases](./phases)

The most frequently used class in Cantera is the {mat:class}`Solution`. It can represent
a mixture of gases, a liquid solution, or a solid solution and provides access to the
solution's thermodynamic, kinetic, and transport properties. The {mat:class}`Interface`
class represents surfaces formed by two adjacent phases or edges where three phases
meet. Several [special constructors](sec-matlab-purefluid) are provided to instantiate
{mat:class}`Solution` objects that implement pure fluid equations of state for certain
substances.

## [Thermodynamic Properties](./thermo)

Class {mat:class}`ThermoPhase` is one of the base classes for {mat:class}`Solution`
objects. It represents the intensive thermodynamic state using one of the [phase
thermodynamic models](/reference/thermo/phase-thermo) implemented by Cantera and
provides methods for computing thermodynamic properties.

The {mat:class}`Mixture` class provides an interface for computing equilibrium
properties of mixtures composed of multiple phases.

## [Chemical Kinetics](./kinetics)

Class {mat:class}`Kinetics` is a base class of {mat:class}`Solution` that provides access
to [reaction rates](/reference/kinetics/reaction-rates) of progress, species production
rates, and other quantities pertaining to a reaction mechanism.

## [Transport Properties](./transport)

Class {mat:class}`Transport` is a base class of {mat:class}`Solution` that provides
access to transport properties such as viscosity and species diffusivities, using one of
the available [transport models](sec-phase-transport-models).

## [Zero-Dimensional Reactor Networks](./zerodim)

A reactor network consists of one or more interconnected reactors. Several reactor types
are implemented by [classes derived from `Reactor`](sec-matlab-reactors), each with its
own set of [governing equations](sec-homogenous-reactor-types).

Reactors can be connected to each other and upstream or downstream
{mat:class}`Reservoir`s using {mat:class}`Valve`s, {mat:class}`MassFlowController`s, and
{mat:class}`Wall`s, which introduce [additional terms to the governing
equations](sec-reactor-interactions). Heterogeneous reactions are handled by
{mat:class}`ReactorSurface`.

Time integration of reactor networks is handled by the {mat:class}`ReactorNet` class.

## [One-dimensional Reacting Flows](./onedim)

A [1D simulation](/reference/onedim/index) in Cantera is set up as a [flow
domain](sec-matlab-flow-domains) surrounded by two [boundary
domains](sec-matlab-boundary-domains) that define inlet/outlet boundary conditions.

These domain objects are combined in a {mat:class}`Sim1D` object that manages the
numerical solution. A specialization {mat:class}`CounterFlowDiffusionFlame` handles
setup and initial condition definitions for diffusion flame simulations.

## [Physical Constants](./constants.rst)

Cantera provides definitions for a number of frequently used physical constants. The
values are consistent with the 2018 CODATA recommendations.

## [Utilities](./utilities.rst)

The Matlab interface also includes a number of other methods for controlling
interactions with the Cantera C++ library, getting information about the C++ library,
and working with other global library state.


```{toctree}
:hidden:
:maxdepth: 1
phases
thermo
kinetics
transport
zerodim
onedim
constants
utilities
```
