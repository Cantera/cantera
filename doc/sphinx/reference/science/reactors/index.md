# Reactors and Reactor Networks

```{caution}
This page is a work in progress. For more complete documentation of the current Cantera
release (Cantera 3.0), please see <a href="/science/reactors/reactors.html">this page</a>.
```

## Reactors

A Cantera *reactor* represents the simplest form of a chemically reacting system. It
corresponds to an extensive thermodynamic control volume $V$, in which all state
variables are homogeneously distributed. The system is generally unsteady -- that is,
all states are functions of time. In particular, transient state changes due to chemical
reactions are possible. However, thermodynamic (but not chemical) equilibrium is assumed
to be present throughout the reactor at all instants of time.

Reactors can interact with the surrounding environment in multiple ways:

- Expansion/compression work: By moving the walls of the reactor, its volume can be
  changed and expansion or compression work can be done by or on the reactor.

- Heat transfer: An arbitrary heat transfer rate can be defined to cross the boundaries
  of the reactor.

- Mass transfer: The reactor can have multiple inlets and outlets. For the inlets,
  arbitrary states can be defined. Fluid with the current state of the reactor exits the
  reactor at the outlets.

- Surface interaction: One or multiple walls can influence the chemical reactions in the
  reactor. This is not just restricted to catalytic reactions, but mass transfer between
  the surface and the fluid can also be modeled.

All of these interactions do not have to be constant, but can vary as a function of time
or state. For example, heat transfer can be described as a function of the temperature
difference between the reactor and the environment, or the wall movement can be modeled
depending on the pressure difference. Interactions of the reactor with the environment
are defined on one or more walls, inlets, and outlets.

In addition to single reactors, Cantera is also able to interconnect reactors into a
*reactor network*. Each reactor in a network may be connected so that the contents of one
reactor flow into another. Reactors may also be in contact with one another or the
environment via walls that conduct heat or move to do work.

## Reactor Types and Governing Equations

All reactor types are modelled using combinations of Cantera's governing equations of
state. The specific governing equations defining Cantera's supported reactor models are
derived and described below.


````{grid} 2
:gutter: 3

```{grid-item-card} Control Volume Reactor
:link: controlreactor
:link-type: doc

A reactor where the volume is prescribed by the motion of the reactor's walls. The state variables are the volume, mass, total internal energy, and species mass fractions.
```

````

```{toctree}
:hidden:

controlreactor
```