```{py:currentmodule} cantera
```

# Constant Pressure Mole Reactor

A constant pressure mole reactor is implemented by the C++ class
{ct}`ConstPressureMoleReactor` and is available in Python as the
{py:class}`ConstPressureMoleReactor` class. It is defined by the state variables:

- $H$, the total enthalpy of the reactor's contents (in J)
- $n_k$, the number of moles for each species (in kmol)

Equations 1 and 2 below are the governing equations for a constant pressure mole
reactor.

## Species Equations

The moles of each species in the reactor changes as a result of flow through the
reactor's [inlets and outlets](sec-flow-device) and production of gas phase species
through homogeneous reactions and reactions on the reactor
[surfaces](sec-reactor-surface). The rate at which species $k$ is generated through
homogeneous phase reactions is $V \dot{\omega}_k$, and the total rate at which moles of
species $k$ changes is:

$$
\frac{dn_k}{dt} = V \dot{\omega}_k + \sum_{in} \dot{n}_{k, in}
                  - \sum_{out} \dot{n}_{k, out} + \dot{n}_{k, wall}
$$ (const-pressure-mole-reactor-species)

Where the subscripts *in* and *out* refer to the sum of the corresponding property over
all inlets and outlets respectively. A dot above a variable signifies a time derivative.

## Energy Equation

Writing the first law for an open system gives:

$$
\frac{dU}{dt} = - p \frac{dV}{dt} + \dot{Q} +
                \sum_{in} \dot{n}_{in} \hat{h}_{in} - \hat{h} \sum_{out} \dot{n}_{out}
$$

where positive $\dot{Q}$ represents heat addition to the system and $h$ is the specific
enthalpy of the reactor's contents.

Differentiating the definition of the total enthalpy, $H = U + pV$, with respect to time
gives:

$$  \frac{dH}{dt} = \frac{dU}{dt} + p \frac{dV}{dt} + V \frac{dp}{dt}  $$

Noting that $dp/dt = 0$ and substituting into the energy equation yields:

$$
\frac{dH}{dt} = \dot{Q} + \sum_{in} \dot{n}_{in} \hat{h}_{in}
                - \hat{h} \sum_{out} \dot{n}_{out}
$$ (const-pressure-mole-reactor-energy)
