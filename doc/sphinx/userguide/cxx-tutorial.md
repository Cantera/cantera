# Using Cantera from C++

One of the main reasons to use the Cantera C++ interface directly is to be able to
couple it with a higher-level simulation code, where Cantera is used to calculate
thermodynamic properties, reaction rates, and transport properties for the host code.
This tutorial demonstrates how to initialize Cantera objects from
[YAML input files](/userguide/input-tutorial) and use those objects to calculate
properties for mixtures with different thermodynamic states.

## A very simple C++ program

A short C++ program that uses Cantera is shown below. This program reads in a
specification of a gas mixture from an input file, and then builds a new object
representing the mixture. It then sets the thermodynamic state and composition of the
gas mixture, and prints out a summary of its properties.

```{literalinclude} demo1a.cpp
:language: c++
```

:::{tip}
Before you can run this program, it first needs to be compiled. On a Linux system using
the GCC compiler, a typical command line for compiling this program might look like
this:

```bash
g++ combustor.cpp -o combustor -O3 $(pkg-config --cflags --libs cantera)
```

This example relies on the `pkg-config` tool to determine the appropriate compiler
flags, such as those specifying the Cantera header and library files. For more advanced
and flexible methods of compiling programs that use the Cantera C++ library, see
[](compiling-cxx).
:::

Running this program produces the output below:

```{literalinclude} ../../../../test_problems/cxx_userguide/demo1a_blessed.txt
:language: none
```

As C++ programs go, this one is quite short. It is the Cantera equivalent of the "Hello,
World" program most programming textbooks begin with. But it illustrates some important
points in writing Cantera C++ programs.

### Making Cantera classes and functions available

Cantera provides several header files designed for use in C++ application programs.
These are designed to include those portions of Cantera needed for particular types of
calculations. The headers and their functions are:

`core.h`
: Base classes and functions for creating {ct}`Solution` and {ct}`Interface` objects
  from input files, as well as associated classes accessed through these objects such as
  {ct}`ThermoPhase`, {ct}`Kinetics`, and {ct}`Transport`. *(New in Cantera 3.0)*.

`zerodim.h`
: Zero-dimensional reactor networks.

`onedim.h`
: One-dimensional reacting flows.

`reactionpaths.h`
: Reaction path diagrams.

### Creating a `Solution` object from an input file

Here, {ct}`newSolution` imports all information held by a YAML input file into a Cantera
{ct}`Solution` object, which is accessed by the pointer `sol`. The thermodynamic
information is accessible via `sol->thermo()`, which itself returns a pointer to a
{ct}`ThermoPhase` object.

Class {ct}`ThermoPhase` is the base class for Cantera classes that represent phases of
matter. It defines the public interface for all classes that represent phases. For
example, it specifies that they all have a method `temperature()` that returns the current
temperature, a method `setTemperature(double T)` that sets the temperature, a method
`getChemPotentials(double* mu)` that writes the species chemical potentials into array
`mu`, and so on.

{ct}`ThermoPhase` objects can represent the intensive state of any single-phase solution
of multiple species. The phase may be a bulk, three-dimensional phase (a gas, a liquid,
or a solid), or it may be a two-dimensional surface phase, or even a one-dimensional
"edge" phase. Many different [phase thermodynamic models](/reference/thermo/phase-thermo)
are implemented as classes derived from {ct}`ThermoPhase`.

### The `report` function

The {ct}`ThermoPhase::report` method generates a nicely-formatted report of the
properties of a phase, including its composition in both mole (X) and mass (Y) units.
For each species present, the non-dimensional chemical potential is also printed. This
is handy particularly when doing equilibrium calculations. This function is very useful
to see at a glance the state of some phase.

### Handling errors

The entire body of the program is put inside a function that is invoked within a `try`
block in the main program. In this way, exceptions thrown in the function or in any
procedure it calls may be caught. In this program, a `catch` block is defined for
exceptions of type {ct}`CanteraError`. Cantera throws exceptions of this type, so it is
always a good idea to catch them.

## Accessing thermodynamic properties

In the program below, a gas mixture object is created, and a few thermodynamic
properties are computed and printed out:

```{literalinclude} thermodemo.cpp
:language: c++
```

Methods that compute scalar properties take no input parameters. The properties are
computed for the state that has been previously set and stored internally within the
object. Methods that return *molar* properties have names that end in `_mole`, while
methods that return properties *per unit mass* have names that end in `_mass`.

Methods that compute an array of values have names beginning with `get` and write the
computed values into an array provided as an input argument. For example, the method
`getChemPotentials(double* mu)` writes the species chemical potentials into the output
array `mu`.

All thermodynamic property methods are declared in class {ct}`ThermoPhase`, which is the
base class from which all classes that represent any type of phase of matter derive. See
the documentation for class {ct}`ThermoPhase` for the full list of available
thermodynamic properties.

## Computing chemical equilibrium

In the program below, the {ct}`ThermoPhase::equilibrate` method is called to set a gas
to a state of chemical equilibrium, holding the temperature and pressure fixed.

```{literalinclude} demoequil.cpp
:language: c++
```

The program output is:

```{literalinclude} ../../../../test_problems/cxx_userguide/demoequil_blessed.txt
:language: none
```

How can we tell that this is really a state of chemical equilibrium? Well, by applying
the equation of reaction equilibrium to formation reactions from the elements, it is
straightforward to show that

$$  \mu_k = \sum_m \lambda_m a_{km}  $$

where $\mu_k$ is the chemical potential of species $k$, $a_{km}$ is the number of atoms
of element $m$ in species $k$, and $\lambda_m$ is the chemical potential of the
elemental species per atom (the so-called "element potential"). In other words, the
chemical potential of each species in an equilibrium state is a linear sum of
contributions from each atom. We can see that this is true in the output above---the
chemical potential of H2 is exactly twice that of H, the chemical potential for OH is
the sum of the values for H and O, the value for H2O2 is twice as large as the value for
OH, and so on.

## Computing reaction rates and transport properties

The following program demonstrates the typical method for accessing the following object
types from a {ct}`Solution` object:

- {ct}`ThermoPhase`: Represents the thermodynamic properties of mixtures containing one
  or more species. Accessed using {ct}`Solution::thermo` method.
- {ct}`Kinetics`: Represents a kinetic mechanism involving one or more phases. Accessed
  using the {ct}`Solution::kinetics`.
- {ct}`Transport`: Computes transport properties for a phase. Accessed using the
  {ct}`Solution::transport` method.

```{literalinclude} kinetics_transport.cpp
:language: c++
```

This program produces the output below:

```{literalinclude} ../../../../test_problems/cxx_userguide/kinetics_transport_blessed.txt
:language: none
```

## Examples of additional C++ functionality

- [`combustor.cpp`](/examples/cxx/combustor): An example that shows how to set up a
  reactor network model, using classes {ct}`Reactor`, {ct}`ReactorNet`, {ct}`Reservoir`,
  {ct}`MassFlowController`, and {ct}`Valve`.
- [`flamespeed.cpp`](/examples/cxx/flamespeed): An example simulating a freely
  propagating laminar flame using the {ct}`Sim1D` solver.
