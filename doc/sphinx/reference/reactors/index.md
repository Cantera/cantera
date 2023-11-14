```{py:currentmodule} cantera
```

# Reactors and Reactor Networks

```{caution}
This page is a work in progress. For more complete documentation of the current Cantera
release (Cantera 3.0), please see <a href="/science/reactors/reactors.html">this page</a>.
```

In Cantera, a *reactor network* represents a set of one or more homogeneous reactors and
reacting surfaces that may be connected to each other and to the environment through
devices representing mass flow, heat transfer, and moving walls. The system is generally
unsteady -- that is, all states are functions of time. In particular, transient state
changes due to chemical reactions are possible.

## Reactor Types and Governing Equations

A Cantera *reactor* represents the simplest form of a chemically reacting system. It
corresponds to an extensive thermodynamic control volume $V$, in which all state
variables are homogeneously distributed. By default, these reactors are closed (no
inlets or outlets) and have adiabatic, chemically-inert walls. These properties may all
be changed by adding components such as walls, surfaces, mass flow controllers, and
valves, as described [below](sec-reactor-interactions).

The specific governing equations defining Cantera's reactor models are derived and
described below. These models represent different combinations of whether pressure or
volume are held constant, whether they support any equation of state or are specialized
for ideal gas mixtures, and whether mass fractions or moles of each species are used as
the state variables representing the composition.

[Control Volume Reactor](controlreactor)
: A reactor where the volume is prescribed by the motion of the reactor's walls.

[Constant Pressure Reactor](constant-pressure-reactor)
: A reactor where the pressure is held constant.

[Ideal Gas Control Volume Reactor](ideal-gas-reactor)
: A reactor where the volume is prescribed by the motion of the reactor's walls,
  specialized for ideal gas mixtures.

[Ideal Gas Constant Pressure Reactor](ideal-gas-constant-pressure-reactor)
: A reactor where the pressure is held constant, specialized for ideal gas mixtures.

[Control Volume Mole Reactor](mole-reactor)
: A reactor where the volume is prescribed by the motion of the reactor's walls, with
  the composition stored in moles.

[Constant Pressure Mole Reactor](constant-pressure-mole-reactor)
: A reactor where the pressure is held constant and the composition is stored in moles.

[Ideal Gas Control Volume Mole Reactor](ideal-gas-mole-reactor)
: A reactor where the volume is prescribed by the motion of the reactor's walls,
  specialized for ideal gas mixtures and with the composition stored in moles.

[Ideal Gas Constant Pressure Mole Reactor](ideal-gas-constant-pressure-mole-reactor)
: A reactor where the pressure is held constant, specialized for ideal gas mixtures and
  with the composition stored in moles.

(sec-reactor-interactions)=
## Reactor Interactions

Reactors can interact with each other and the surrounding environment in multiple ways.
Mass can flow from one reactor into another can be incorporated, heat can be
transferred, and the walls between reactors can move. In addition, reactions can occur
on surfaces within a reactor. The models used to establish these interconnections are
described in the following sections:

All of these interactions can vary as a function of time or system state. For example,
heat transfer can be described as a function of the temperature difference between the
reactor and the environment, or wall movement can be modeled as depending on the
pressure difference. Interactions of the reactor with the environment are defined using
the following models:

[Reservoirs](sec-reservoir)
: Reservoirs are used to represent constant conditions defining the inlets, outlets, and
  surroundings of a reactor network.

[Flow Devices](sec-flow-device)
: Flow devices are used to define mass transfer between two reactors, or between
  reactors and the surrounding environment as defined by a reservoir.

[Walls](sec-wall)
: Walls between reactors are used to allow heat transfer between reactors. By moving the
  walls of the reactor, its volume can be changed and expansion or compression work can
  be done by or on the reactor.

[Reacting Surfaces](sec-reactor-surface)
: Reactions may occur on the surfaces of a reactor. These reactions may include both
  catalytic reactions and reactions resulting in net mass transfer between the surface
  and the fluid.

```{seealso}
Cantera comes with a broad variety of well-commented example scripts for reactor
networks. Please see the [Cantera Examples](/examples/python/reactors/index) for further
information.
```

## Plug Flow Reactor

A *plug flow reactor* (PFR) represents a steady-state flow in a channel. The fluid is
considered to be homogeneous perpendicular to the flow direction, while the state of the
gas is allowed to change in the axial direction. However, all diffusion processes are
neglected.

These assumptions result in a system of equations that is similar to those used to model
homogeneous reactors, with the time variable replaced by the axial coordinate. Because
of this mathematical similarity, PFRs are also solved by Cantera's reactor network
model. However, they can only be simulated alone, and not part of a network containing
time-dependent reactors.

[Plug Flow Reactor](pfr)
: A reactor modeling one-dimensional steady-state flow in a channel that may contain
  catalytically active surfaces where heterogeneous reactions occur.

```{toctree}
:hidden:
:caption: Reactor models
:maxdepth: 1

controlreactor
constant-pressure-reactor
ideal-gas-reactor
ideal-gas-constant-pressure-reactor
mole-reactor
constant-pressure-mole-reactor
ideal-gas-mole-reactor
ideal-gas-constant-pressure-mole-reactor

interactions
pfr
```
