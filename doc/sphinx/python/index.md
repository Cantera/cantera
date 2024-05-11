```{py:currentmodule} cantera
```

(sec-python-documentation)=

# Python Module Documentation

## [Objects Representing Phases](./importing)

The most frequently used class in Cantera is the {py:class}`Solution`. It can represent
a mixture of gases, a liquid solution, or a solid solution and provides access to the
solution's thermodynamic, kinetic, and transport properties. The {py:class}`Interface`
class represents surfaces formed by two adjacent phases or edges where three phases
meet. Several [special constructors](sec-python-purefluid) are provided to instantiate
{py:class}`Solution` objects that implement pure fluid equations of state for certain
substances.

The {py:class}`Quantity` class represents a specific quantity (mass) of a
{py:class}`Solution`. It provides methods for accessing the extensive properties of the
quantity and computing the state resulting from mixing two substances.

The {py:class}`SolutionArray` class provides a convenient interface for representing a
multidimensional array of states and accessing properties of those states as arrays.

## [Thermodynamic Properties](./thermo)

Class {py:class}`ThermoPhase` is one of the base classes for {py:class}`Solution`
objects. It represents the intensive thermodynamic state using one of the [phase
thermodynamic models](/reference/thermo/phase-thermo) implemented by Cantera and
provides methods for computing thermodynamic properties.

A phase is composed of a number of {py:class}`Species` objects. Each species has a
{py:class}`SpeciesThermo` object that implements a particular [species thermodynamic
model](/reference/thermo/species-thermo). {py:class}`Element` objects can be used to
access information about the elements comprising each species or to define custom
isotopes.

The {py:class}`Mixture` class provides an interface for computing equilibrium properties
of mixtures composed of multiple phases.

## [Chemical Kinetics](./kinetics)

Class {py:class}`Kinetics` is a base class of {py:class}`Solution` that provides access
to [reaction rates](/reference/kinetics/reaction-rates) of progress, species production
rates, and other quantities pertaining to a reaction mechanism.
{py:class}`InterfaceKinetics` provides this information for {py:class}`Interface`
objects.

A reaction mechanism consists of a set of {py:class}`Reaction` objects. The
{py:class}`Reaction` object specifies the reactants and products of the reaction and
has a {py:class}`ReactionRate` object that handles the evaluation of the rate constant
as a function of the mixture composition, using one of the [rate
parameterizations](/reference/kinetics/rate-constants) implemented by Cantera or a
user-defined rate implemented using the {py:class}`ExtensibleRate` class.

The {py:class}`ReactionPathDiagram` class can be used to analyze reaction pathways.

## [Transport Properties](./transport)

Class {py:class}`Transport` is a base class of {py:class}`Solution` that provides access
to transport properties such as viscosity and species diffusivities, using one of the
available [transport models](sec-phase-transport-models). The
{py:class}`DustyGasTransport` provides a special transport model applicable to porous
media. Species properties needed to compute transport properties are handled by the
{py:class}`GasTransportData` class.

## [Zero-Dimensional Reactor Networks](./zerodim)

A reactor network consists of one or more interconnected reactors. Several reactor types
are implemented by [classes derived from `Reactor`](sec-python-reactors), each with its
own set of [governing equations](sec-homogenous-reactor-types).

Reactors can be connected to each other and upstream or downstream
{py:class}`Reservoir`s using {py:class}`Valve`s, {py:class}`MassFlowController`s,
{py:class}`PressureController`s, and {py:class}`Wall`s, which introduce [additional
terms to the governing equations](sec-reactor-interactions). Heterogeneous reactions are
handled by {py:class}`ReactorSurface`. User-defined reactor governing equations can be
implemented using the {py:class}`ExtensibleReactor` class or any of the classes derived
from it.

Time integration of reactor networks is handled by the {py:class}`ReactorNet` class. For
[certain types of networks](sec-reactor-preconditioning), integration can be accelerated
by using the {py:class}`AdaptivePreconditioner` class.

## [One-dimensional Reacting Flows](./onedim)

The main way of setting up a [1D simulation](/reference/onedim/index) in Cantera is
through one of the specializations of the {py:class}`FlameBase` class:
{py:class}`FreeFlame`, {py:class}`BurnerFlame`, {py:class}`CounterflowDiffusionFlame`,
{py:class}`CounterflowPremixedFlame`, {py:class}`CounterflowTwinPremixedFlame`, or
{py:class}`ImpingingJet`. These classes all consist of a [flow
domain](sec-python-flow-domains) with two [boundaries](sec-python-boundary-domains)
defining inlet/outlet boundary conditions.

## [Mechanism Conversion](./scripts)

Cantera includes modules implementing conversion of mechanisms between the YAML and
Chemkin formats, conversion from the LXCat format, and conversion from the legacy CTI
and CTML (XML) formats that were used prior to Cantera 3.0. The preferred interfaces for
these conversions are the executable scripts [`ck2yaml`](/yaml/ck2yaml),
[`yaml2ck`](/yaml/yaml2ck), [`lxcat2yaml`](/yaml/lxcat2yaml),
[`cti2yaml`](/yaml/cti2yaml), and [`ctml2yaml`](/yaml/ctml2yaml).

## [Python Interface With Units](./units.rst)

This subpackage provides an alternative interface to thermodynamic properties where
physical units are associated with all values, using the
[pint](https://pint.readthedocs.io/en/stable/) library. This capability is implemented by the {py:class}`with_units.Solution` and {py:class}`with_units.PureFluid` classes.

## [Physical Constants](./constants.rst)

Cantera provides definitions for a number of frequently used physical constants. The
values are consistent with the 2018 CODATA recommendations.

## [Utilities](./utilities.rst)

Classes {py:class}`AnyMap` and {py:class}`YamlWriter` provide functionality for
interacting with input data defined using Cantera's [YAML input format](/yaml/index).
Classes {py:class}`UnitSystem` and {py:class}`Units` are used for expressing dimensional
quantities in input files.

A number of [global functions](sec-python-global-funcs) are provided for managing
locations where Cantera looks for data files, for setting how Cantera handles warnings
and errors, and how some transitional/legacy features are handled.

The {py:func}`extension` function is used as a decorator for implementations of
{py:class}`ExtensibleRate` to allow user-defined rate types to be specified from YAML
input files.

Classes {py:class}`CanteraError` and {py:class}`ThermoModelMethodError` are raised for
types of errors that are not mapped onto a built-in Python exception type.

```{toctree}
:hidden:

importing
thermo
kinetics
transport
zerodim
onedim
scripts
units
constants
utilities
```
