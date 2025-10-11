---
file_format: mystnb
kernelspec:
  name: python3
---

```{py:currentmodule} cantera
```

# Extending Reactor Models with Custom Equations

In some cases, Cantera's built-in reactor types are insufficient to model a problem.
In this situation, the {py:class}`ExtensibleReactor` family of classes can be used to
modify and governing equations, starting from one of the built-in reactor types.

These reactor types allow the methods that implement the governing equations and related
reactor configuration to be augmented or replaced with user-defined Python functions.
New state variables can be added to the reactor, and existing ones can be redefined.
User-defined `ExtensibleReactor` implementations can be used alongside the built-in
reactor types to build reactor networks that can be integrated using the
{py:class}`ReactorNet` integrator provided by Cantera.

This guide introduces the use of the `ExtensibleReactor` model to modify the governing
equations for a reactor.

## Modifying Existing Governing Equations

The variables in the governing equations for a reactor that are differentiated with
respect to time are known as the *state variables*. The state variables depend on the
type of reactor base class chosen. For example, choosing an
[ideal gas constant pressure reactor](/reference/reactors/ideal-gas-constant-pressure-reactor)
allows the user to modify the governing equations for the following state variables:

- $m$, the mass of the reactor's contents (in kg)
- $T$, the temperature (in K)
- $Y_k$, the mass fractions for each species (dimensionless)

As shown in the derivations of the governing equations, the equations implemented by
the {ct}`IdealGasConstPressureReactor` class are:

$$
\frac{dm}{dt} = \sum_\t{in} \dot{m}_\t{in} - \sum_\t{out} \dot{m}_\t{out}
    + \dot{m}_\t{wall}
$$

$$
m c_p \frac{dT}{dt} = \dot{Q} - \sum_k h_k \dot{m}_{k,\t{gen}}
     + \sum_\t{in} \dot{m}_\t{in} \left(h_\t{in} - \sum_k h_k Y_{k,\t{in}} \right)
$$

$$
\frac{d(mY_k)}{dt} = \sum_\t{in} \dot{m}_\t{in} Y_{k,\t{in}}
                     - \sum_\t{out} \dot{m}_\t{out} Y_k + \dot{m}_{k,gen}
$$

Each of these equations is written with an expression on the left-hand side (LHS)
multiplied by the time derivative of a state variable which equals an expression on the
right-hand side (RHS). The interface of the reactor network model is defined so that the
reactor governing equations are evaluated by calculating these LHS and RHS expressions
for each governing equation. The extensible reactor model then allows calculation of
these LHS and RHS terms to be replaced or augmented by user-provided methods.

For example, to add a term for a large mass, say a rock, inside the reactor that affects
the thermal mass, the energy equation would become:

$$
\left(m c_p + m_\t{rock} c_{p,\t{rock}}\right) \frac{dT}{dt} = \dot{Q}
    - \sum_k h_k \dot{m}_{k,\t{gen}}
    + \sum_\t{in} \dot{m}_\t{in} \left(h_\t{in} - \sum_k h_k Y_{k,\t{in}} \right)
$$

Here, the LHS coefficient has changed from $m c_p$ to
$m c_p + m_{\t{rock}} c_{p,\t{rock}}$. Since the rock does not change the
composition of the species in the reactor and does not change the mass flow rate of any
inlets or outlets, the other governing equations defining the ideal gas constant
pressure reactor can be left unmodified. To implement this change, we define a new class
derived from {py:class}`ExtensibleIdealGasConstPressureReactor`. In this new class, we
then define an implementation for the `after_eval` method which has the signature

```py
def after_eval(self, t, LHS, RHS) -> None
```

where `t` is the current simulation time and `LHS` and `RHS` are arrays with one element
for each state variable. In this case, the length of each array is two plus the number of species.

Optional user-defined methods for the `ExtensibleReactor` class start with the prefix
`before_`, `after_`, or `replace_`, indicating whether the this method will called
before, after, or instead of the corresponding method from the base class. In this case,
our method will be called after the `eval` method of {ct}`IdealGasConstPressureReactor`,
with the values of `LHS` and `RHS` already set based on the original governing
equations. Instead of `after_eval`, we could implement either `before_eval` or
`replace_eval` methods. However, in the case of `before_eval`, any element of `LHS` that
we set would be overwritten by the initial implementation, and if we use `replace_eval`,
we would have needed to calculate values for all elements of `RHS` and `LHS` instead of
just modifying a single value.

To obtain the above modification to the governing equations, the implementation of
`after_eval` inside a new reactor class is then:

```{code-cell} python
import cantera as ct

class RockReactor(ct.ExtensibleIdealGasConstPressureReactor):
    def __init__(self, *args, mass_rock, cp_rock, **kwargs):
        super().__init__(*args, **kwargs)
        self.mass_rock = mass_rock
        self.cp_rock = cp_rock

    def after_eval(self, t, LHS, RHS):
        # The index 1 corresponds to the governing equation for the second
        # state variable (temperature)
        LHS[1] += self.mass_rock * self.cp_rock
```

We can then use this new reactor class in a reactor network with some different values
of `mass_rock` to see the impact of the modified governing equations. To do this, we
can define a reactor network with an imposed heat transfer rate.


```{code-cell} python
cp_rock = 790  # J/kg/K

for mass_rock in [0.0, 1.0, 3.0]:
    # Define objects, properties, and initial conditions.
    gas = ct.Solution('h2o2.yaml')
    gas.TPX = 300, ct.one_atm, {'O2': 1, 'N2': 3.76}
    T0 = gas.T

    r1 = RockReactor(gas, mass_rock=mass_rock, cp_rock=cp_rock)
    r1.volume = 2.0  # volume (not including 'rock')

    env = ct.Reservoir(gas)
    wall = ct.Wall(env, r1, Q=1e4, A=1.0)  # Heat transfer of 10 kW
    net = ct.ReactorNet([r1])

    # Integrate for 10 seconds
    net.advance(10)

    print(f'mass_rock = {mass_rock} kg; ΔT = {r1.phase.T - T0:.2f} K')
```

As expected, we see that the temperature rise decreases as the thermal mass of the
reactor's contents increases.

Other modifications to the reactor's governing equations and behavior may be made by
implementing `before_*`, `after_`, and `replace_*` versions of any of the reactor
methods listed in the documentation for the {py:class}`ExtensibleReactor` class.

## Additional Examples

[Using `ExtensibleReactor` to implement wall inertia](/examples/python/reactors/custom2)
: This example shows how to introduce additional governing equations to a reactor. Here,
  in addition to the `eval` method, overrides are introduced that modify the behavior of
  the `initialize`, `get_state`, `update_state`, `component_index`, and `component_name`
  methods. The result is a reactor network with a wall between two reactors that has
  inertia and takes time to respond to changes in the reactor's pressure.

[Reactor cascade model for reactive flows in inert porous media](/examples/python/reactors/porous_media_burner)
: This example showcases the use of {py:class}`ExtensibleReactor` to add a temperature
  equation for a solid phase and custom heat transfer/radiation models. Each reactor in
  the network represents a different spatial locations along the length of a porous
  media burner.
