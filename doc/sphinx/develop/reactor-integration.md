```{py:currentmodule} cantera
```

# How Time Integration of Reactor Networks Works

This section provides a developer-oriented description of how time integration of
reactor networks works in Cantera.

A {ct}`ReactorNet` object doesn't perform time integration on its own. Instead, the role
of the {ct}`ReactorNet` object is to assemble a system of governing equations from
multiple reactors. It then uses the CVODES or IDAS integrators from SUNDIALS to
integrate the governing equations in time or space, respectively.

## Creating a Reactor and Reactor Network

First, let's take a look at what happens when creating a reactor network by setting
up an isolated reactor in Python.

```python
import cantera as ct
gas = ct.Solution("gri30.yaml")
gas.TPX = 1000.0, ct.one_atm, "H2:2,O2:1,N2:4"
reac = ct.IdealGasReactor(gas, clone=False)
sim = ct.ReactorNet([reac])
```

The constructor for the {py:class}`ReactorNet` class builds the network from the
provided list of {py:class}`Reactor` and {py:class}`ReactorSurface` objects and creates
a new {ct}`Integrator` used to integrate the governing equations. If there are
connections such as walls and mass flow controllers between different reactors, the
adjacent reactors are automatically added to the reactor network.

The {ct}`Integrator` class is Cantera's interface for ODE/DAE system integrators.
{ct}`Integrator` is a [polymorphic base class](http://www.cplusplus.com/doc/tutorial/polymorphism/);
it defines a set of *virtual methods* that derived classes (the actual system
integrators) provide implementations for.

The {ct}`newIntegrator` factory function creates and returns an {ct}`Integrator` object
of the specified type. Calling `newIntegrator("CVODE")` creates a new
{ct}`CVodesIntegrator` object for integrating an ODE system, while calling
`newIntegrator("IDA")` creates a new {ct}`IdasIntegrator` object for integrating a DAE
system. The appropriate integrator type is determined by the {ct}`ReactorNet` class
based on the types of the installed reactors. A {ct}`FlowReactor` defines a DAE system
and uses the IDAS integrator, while the other reactor types define ODE systems and use
the CVODES integrator. In this guide, the implementation of {ct}`CVodesIntegrator` is
described; {ct}`IdasIntegrator` works in a similar way, although the IDAS function names
are different.

## Communicating with SUNDIALS using Wrapper Functions

Because SUNDIALS is written in C, the {ct}`CVodesIntegrator` C++ wrapper is used to
access the `CVODES` solver and make the appropriate calls such as those to the
integrator function `CVode`. For details on CVODES API, see the
[CVODES User Guide](https://sundials.readthedocs.io/en/latest/cvodes/Introduction_link.html).

## Calling the `ReactorNet.advance` Method

After configuring a reactor network and its components in Cantera, a call to the
{py:meth}`ReactorNet.advance` method can be used to obtain the state of the network at a
specified time. The initial condition information is passed off to the {ct}`Integrator`
when calling `advance`. Transient physical and chemical interactions are simulated by
integrating the network's system of ODE governing equations through time.

```python
sim.advance(1)  # advance the simulation to the specified absolute time, t = 1 sec
gas()  # view the updated state of the mixture, reflecting properties at t = 1 sec
```

Calling the {ct}`ReactorNet::advance` method invokes the method
{ct}`CVodesIntegrator::integrate` to integrate the system to the specified time. It is
most efficient to let CVODES determine the actual integration step sizes on its own.
Therefore, we take individual time steps using CVODES using the `CVode` function until
we have reached or passed the specified time. We then interpolate the system state back
to the specified time using the `CVodeGetDky` function. With some additional handling of
special cases and errors, the implementation of {ct}`CVodesIntegrator::integrate` is:

```{literalinclude} ../../../../src/numerics/CVodesIntegrator.cpp
:language: c++
:start-at: "void CVodesIntegrator::integrate("
:end-before: " CVodesIntegrator::"
```

The arguments taken by the SUNDIALS
[`CVode()`](https://sundials.readthedocs.io/en/latest/cvode/Usage/index.html#cvode-solver-function)
function are:

- {ct}`CVodesIntegrator::m_cvode_mem`, a pointer to the block of memory that was
  allocated and configured during initialization.
- `tout` is the desired integrator output time. CVODES will not necessarily reach this
  time when operating in "one step" mode, but it is used in the selection of the initial
  step size.
- After execution, {ct}`CVodesIntegrator::m_y` will contain the computed system state
  at the time reached by the integrator, and will later be used to update the
  {ct}`ReactorNet` to its time-integrated state.
- After execution, {ct}`CVodesIntegrator::m_tInteg` will contain the time reached by the
  integrator.
- The `CV_ONE_STEP` option tells the solver to take a single internal step.

The return value of the `CVode()` function is assigned to the `flag` variable. `CVode()`
returns the constant `CV_SUCCESS` to indicate success an error code if integration was
unsuccessful.

## Specifying ODEs using Class `FuncEval`

How does `CVODES` know what ODE system it should be solving?

In the example above, the ODE system was actually already specified using `CVodeInit()`,
one of the functions automatically invoked by the {ct}`ReactorNet::initialize` method.
CVODES requires a C function with a specific signature that defines the ODE system by
computing the right-hand side of the ODE system $dy/dt$ for a given value of the
independent variable, $t$, and the state vector $y$. For more information about ODE
right-hand side function requirements, see the
[CVODES User Guide](https://sundials.readthedocs.io/en/latest/cvodes/Usage/SIM.html#user-supplied-functions).

The {ct}`CVodesIntegrator` wrapper class provides a useful C++ interface for configuring
this C function by pairing with {ct}`FuncEval`, an abstract base class for ODE and DAE
right-hand-side function evaluators. Classes derived from {ct}`FuncEval` implement the
evaluation of the provided ODE system.

An ODE right-hand-side evaluator is always needed in the ODE solution process to provide
the definition of the system, and for that reason a {ct}`FuncEval` object is a required
parameter when initializing any type of {ct}`Integrator`.

{ct}`ReactorNet` handles this requirement by inheriting the {ct}`FuncEval` class,
meaning that it provides the implementation for the ODE function and actually specifies
itself using the [`this`](https://en.cppreference.com/w/cpp/language/this) pointer
when calling {ct}`Integrator::initialize` in {ct}`ReactorNet::initialize`.

To be a valid {ct}`FuncEval` object, a {ct}`ReactorNet` needs to provide implementations
for several of {ct}`FuncEval`'s virtual methods, particularly the actual ODE
right-hand-side computation method, {ct}`FuncEval::eval`:

```C++
virtual void eval(double t, double* y, double* ydot, double* p);
```

where the arguments are:

- `t`, the current time in seconds.
- `y`, a pointer to the start of the current state vector.
- `ydot`, a pointer to the start of the array where the computed time derivatives should
  be stored.
- `p`, a pointer to the start of an array containing the (potentially perturbed)
  sensitivity parameters, if any have been defined.

To take a single timestep, CVODES will call `eval()` one or more times and use the
computed values of `ydot` to determine the value of `y` at the new time.

Within {ct}`ReactorNet::eval`, the governing equations for each each reactor in the
network are evaluated and assembled to form the full `ydot` vector for the system. To
handle some complexities of the `Reactor` model, namely the fact that multiple
`Reactor`s can share the same `ThermoPhase` object while having different states, the
governing equation evaluation takes place in two steps. First, the state of each
`Reactor` is set according to the values in the global state vector `y` using the
{ct}`ReactorNet::updateState` method, which calls the {ct}`Reactor::updateState` method
for each reactor:

```{literalinclude} ../../../../src/zeroD/ReactorNet.cpp
:start-at: "void ReactorNet::updateState"
:end-before: " ReactorNet::"
:language: c++
```

To simplify implementation of some modifications of the governing equations, for
example, using the {py:class}`ExtensibleReactor` class, the governing equations for each
reactor are written in the form:

$$
\t{LHS}_i \frac{dy_i}{dt} = \t{RHS}_i
$$

where the {ct}`Reactor::eval` method or the `eval()` method of any class derived from
{ct}`Reactor` calculates the values for the LHS (left-hand side) and RHS (right-hand
side) vectors, whose values default to 1 and 0, respectively, by implementing a method
with the signature:

```c++
void eval(double time, double* LHS, double* RHS);
```

These values are then assembled into the global `ydot` vector by `ReactorNet::eval()`:

```{literalinclude} ../../../../src/zeroD/ReactorNet.cpp
:start-at: "void ReactorNet::eval("
:end-before: " ReactorNet::"
:language: c++
```
