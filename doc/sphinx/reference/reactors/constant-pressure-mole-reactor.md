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
reactor's [inlets and outlets](sec-flow-device), and production of homogeneous gas phase
species and reactions on the reactor [surfaces](sec-reactor-surface). The rate at which
species $k$ is generated through homogeneous phase reactions is $V \dot{\omega}_k$, and
the total rate at which moles of species $k$ changes is:

$$
\frac{dn_k}{dt} = V \dot{\omega}_k + \sum_{in} \dot{n}_{k, in}
                  - \sum_{out} \dot{n}_{k, out} + \dot{n}_{k, wall}
$$ (const-pressure-mole-reactor-species)

Where the subscripts *in* and *out* refer to the sum of the corresponding property over
all inlets and outlets respectively. A dot above a variable signifies a time derivative.

## Energy Equation

Taking the definition of the total enthalpy and differentiating with respect to time
gives:

$$
H &= U + pV

\frac{dH}{dt} &= \frac{d U}{d t} + p \frac{dV}{dt} + V \frac{dp}{dt}
$$

Noting that $dp/dt = 0$ and substituting into the control volume mole reactor energy
equation {eq}`molereactor-energy` yields:

$$
\frac{dH}{dt} = \dot{Q} + \sum_{in} \dot{n}_{in} \hat{h}_{in}
                - \hat{h} \sum_{out} \dot{n}_{out}
$$ (const-pressure-mole-reactor-energy)

where positive $\dot{Q}$ represents heat addition to the system.
