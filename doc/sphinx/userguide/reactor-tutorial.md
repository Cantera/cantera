---
file_format: mystnb
kernelspec:
  name: python3
---

```{py:currentmodule} cantera
```

# Using Reactor Networks for Time Integration

This guide explains the basics of setting up and integrating the governing equations of
a transient reactor network problem.

## Setting up a Reactor Network

First, let's take a look at a basic example to see how we might utilize Cantera's time
integration functionality. We'll simulate an isolated reactor in Python that is filled
with a homogeneous gas mixture. The gas state used in this example is arbitrary, but
interesting because it's explosive.

```{code-cell} python
# import the Cantera Python module
import cantera as ct

# Create a gas mixture using the GRI-Mech 3.0 mechanism.
# This mechanism is included with Cantera for use in examples, but is not
# recommended for research use.
gas = ct.Solution("gri30.yaml")

# set gas to an interesting state
gas.TPX = 1000.0, ct.one_atm, "H2:2,O2:1,N2:4"

# create a reactor containing the gas
reac = ct.IdealGasReactor(gas)

# Add the reactor to a new ReactorNet simulator
sim = ct.ReactorNet([reac])

# View the initial state of the mixture
reac.phase()
```

Now, let's advance the simulation in time to an absolute time of 1 second, and examine
how the state of the gas has changed. There isn't anything special about choosing 1
second as the end time for the integration; any end time can be chosen. We are choosing
1 second here because after that amount of time, the mixture will have almost certainly
ignited for a wide range of initial conditions.

```{code-cell} python
# advance the simulation to the specified absolute time, t = 1 sec
sim.advance(1)

# view the updated state of the mixture, reflecting properties at t = 1 sec
reac.phase()
```

## Methods for Time Integration

The state of a {py:class}`ReactorNet` can be advanced in time by one of the following
three methods:

- `step()`: The `step()` method computes the state of the system after one time step.
  The size of the step is determined by SUNDIALS when the method is called. SUNDIALS
  determines the step size by estimating the local error, which must satisfy tolerance
  conditions. The step is redone with reduced step size whenever that error test fails.
  SUNDIALS also periodically checks if the maximum step size is being used. The time
  step must not be larger than a predefined maximum time step $\Delta t_\t{max}$. The
  new time $t_\t{new}$ at the end of the single step is returned by this function. This
  method produces the highest time resolution in the output data of the methods
  implemented in Cantera.

- `advance(t_new)`: This method computes the state of the system at the user-provided
  time $t_\t{new}$. $t_\t{new}$ is the absolute time from the initial time of the
  system. Although the user specifies the time when integration should stop, SUNDIALS
  chooses the time step size as the network is integrated. Many internal SUNDIALS time
  steps are usually required to reach $t_\t{new}$. As such, `advance(t_new)` preserves
  the accuracy of using `step()` but allows consistent spacing in the output data.

- `advance_to_steady_state(max_steps, residual_threshold, atol, write_residuals)`
  *Python interface only*: If the steady state solution of a reactor network is of
  interest, this method can be used. Internally, the steady state is approached by time
  stepping. The network is considered to be at steady state if the feature-scaled
  residual of the state vector is below a given threshold value (which by default is 10
  times the time step `rtol`).

### Limiting Timesteps

The `advance(t_new)` approach is typically used when consistent, regular, time steps are
required in the output data. This usually simplifies comparisons among many simulation
results at a single time point. However, some detail, for example, a fast ignition
process, might not be resolved in the output data if the user-provided time step is not
small enough.

To avoid losing this detail, the {ct}`Reactor::setAdvanceLimit` method (C++) or the
{py:func}`Reactor.set_advance_limit` method (Python) can be used to set the maximum
amount that a specified solution component can change between output times. For an
example of this feature's use, see the example
[reactor1.py](/examples/python/reactors/reactor1). However, as a tradeoff, the time step
sizes in the output data are no longer guaranteed to be uniform.

### Setting Tolerances

Even though Cantera comes pre-defined with reasonable default values for tolerances and
the maximum internal time step, the solver may not be correctly configured for the
specific system. In this case SUNDIALS may stop with an error. To solve this problem,
three parameters can be tuned:

1. the absolute tolerances
2. the relative tolerances
3. the maximum time step size

Reducing the third value is particularly useful when dealing with abrupt changes in the
boundary conditions. One example of this is opening or closing valves; see also the
[IC engine example](/examples/python/reactors/ic_engine).

## Adding Inlets and Outlets

The following examples demonstrate the use of [flow devices](sec-flow-device) such
as valves and mass flow controllers to model the inlets and outlets of reactors:

- [](/examples/python/reactors/mix1)
- [](/examples/python/reactors/fuel_injection)
- [](/examples/python/reactors/ic_engine)
- [](/examples/cxx/combustor) (C++)

## Modeling Volume Changes and Heat Transfer

The following examples demonstrate the use of [walls](sec-wall) to model volume changes
and heat transfer between reactors:

- [](/examples/python/reactors/piston)
- [](/examples/python/reactors/reactor2)
- [](/examples/python/reactors/ic_engine)

Additional examples can be found in the
[Python Reactor Network Examples](/examples/python/reactors/index) section.

## Accelerating Integration with Preconditioned Sparse Iterative Methods

The [](/reference/reactors/ideal-gas-mole-reactor) and the
[](/reference/reactors/ideal-gas-constant-pressure-mole-reactor) are able to be
solved using preconditioned sparse iterative methods, instead of the normal direct
linear methods. Using these solvers can greatly speed up integration for systems with
hundreds or thousands of species {cite:p}`walker2023`.

To enable the use of this solver, a preconditioner object needs to be created and
specified for use by the reactor network. The previous example can be modified as
follows:

```{code-cell} python
import cantera as ct
gas = ct.Solution("gri30.yaml")
gas.TPX = 1000.0, ct.one_atm, "H2:2,O2:1,N2:4"
reac = ct.IdealGasMoleReactor(gas)  # Use reactor type that supports preconditioning
sim = ct.ReactorNet([reac])
precon = ct.AdaptivePreconditioner()  # Create the preconditioner
sim.preconditioner = precon  # Add it to the network
```

The results of the time integration are the same as before:
```{code-cell} python
sim.advance(1)
reac.phase()
```

The approximate Jacobian matrices used to construct the preconditioners do not currently
account for terms that link multiple reactors in the same network. For networks of this
type, the iterative solver may exhibit convergence errors for systems where the default
direct solver does not. If the solver converges, the solution will satisfy the error
tolerances, even though the Jacobian matrix is not exact.

## Common Reactor Types and their Implementation in Cantera

### Batch Reactor at Constant Volume or at Constant Pressure

If you are interested in how a homogeneous chemical mixture changes state in time when
isolated from the surroundings, a simple batch reactor can be used. Two versions are
commonly considered: a rigid vessel with fixed volume but variable pressure, or a system
idealized at constant pressure but varying volume.

The initial state of the solution can be specified by composition and a set of intensive
thermodynamic parameters, like temperature and pressure, as a standard Cantera
{py:class}`Solution` object. From this base, a {py:class}`Reactor` or a
{py:class}`ConstPressureReactor` can be created, depending on if a constant volume or
constant pressure batch reactor should be considered, respectively. The behavior of the
solution in time can be simulated as a {py:class}`ReactorNet` containing only the single
reactor.

An example for such a batch reactor is given in
[`reactor1.py`](/examples/python/reactors/reactor1).

### Continuously Stirred Tank Reactor

A Continuously Stirred Tank Reactor (CSTR), also often referred to as Well-Stirred
Reactor (WSR), Perfectly Stirred Reactor (PSR), or Longwell Reactor, is essentially a
single Cantera reactor with an inlet, an outlet, and constant volume.

Steady state solutions to CSTRs are often of interest. In this case, the mass flow rate
$\dot{m}$ is constant and equal at inlet and outlet. The mass contained in the control
volume, $m$, divided by $\dot{m}$ defines the mean residence time of the fluid in the
control volume.

While Cantera always solves a transient problem, if you are interested in steady-state
conditions, you can run your simulation for a long time until the states are converged.
This integration process can be managed automatically using the
{py:meth}`ReactorNet.advance_to_steady_state` method. See
[`combustor.py`](/examples/python/reactors/combustor) for an example using this method.

:::{tip}
A problem can be the ignition of a CSTR: If the reactants are not reactive enough, the
simulation can result in the trivial solution that inflow and outflow states are
identical. To solve this problem, the reactor can be initialized with a high temperature
and/or radical concentration. A good approach is to use the equilibrium composition of
the reactants (which can be computed using Cantera's {py:meth}`ThermoPhase.equilibrate`
function) as an initial guess. This method is demonstrated in the
[`combustor.py`](/examples/python/reactors/combustor) example.
:::
