```{py:currentmodule} cantera
```

# Constant Pressure Reactor

For this reactor model, the pressure is held constant and the energy equation is defined
in terms of the total enthalpy. This model is implemented by the C++ class
{ct}`ConstPressureReactor` and available in Python as the
{py:class}`ConstPressureReactor` class. A constant pressure reactor is defined by the
state variables:

- $m$, the mass of the reactor's contents (in kg)
- $H$, the total enthalpy of the reactor's contents (in J)
- $Y_k$, the mass fractions for each species (dimensionless)

Equations 1-3 below are the governing equations for a constant pressure reactor.

## Mass Conservation

The total mass of the reactor's contents changes as a result of flow through the
reactor's [inlets and outlets](sec-flow-device), and production of homogeneous phase
species on [surfaces](sec-reactor-surface):

$$
\frac{dm}{dt} = \sum_{in} \dot{m}_{in} - \sum_{out} \dot{m}_{out} + \dot{m}_{wall}
$$ (constpressurereactor-mass)

Where the subscripts *in* and *out* refer to the sum of the superscripted property over
all inlets and outlets respectively. A dot above a variable signifies a time derivative.

## Species Equations

The rate at which species $k$ is generated through homogeneous phase reactions is $V
\dot{\omega}_k W_k$, and the total rate at which species $k$ is generated is:

$$  \dot{m}_{k,gen} = V \dot{\omega}_k W_k + \dot{m}_{k,wall}  $$

The rate of change in the mass of each species is:

$$
\frac{d(mY_k)}{dt} = \sum_{in} \dot{m}_{in} Y_{k,in} - \sum_{out} \dot{m}_{out} Y_k +
                     \dot{m}_{k,gen}
$$

Expanding the derivative on the left hand side and substituting the equation for
$dm/dt$, the equation for each homogeneous phase species is:

$$
m \frac{dY_k}{dt} = \sum_{in} \dot{m}_{in} (Y_{k,in} - Y_k) +
                    \dot{m}_{k,gen} - Y_k \dot{m}_{wall}
$$ (constpressurereactor-species)

## Energy Equation

Writing the first law for an open system gives:

$$
\frac{dU}{dt} = - p \frac{dV}{dt} + \dot{Q} +
                \sum_{in} \dot{m}_{in} h_{in} - h \sum_{out} \dot{m}_{out}
$$

where positive $\dot{Q}$ represents heat addition to the system and $h$ is the specific
enthalpy of the reactor's contents.

Differentiating the definition of the total enthalpy, $H = U + pV$, with respect to time
gives:

$$  \frac{dH}{dt} = \frac{dU}{dt} + p \frac{dV}{dt} + V \frac{dp}{dt}  $$

Noting that $dp/dt = 0$ and substituting into the energy equation yields:

$$
\frac{dH}{dt} = \dot{Q} + \sum_{in} \dot{m}_{in} h_{in} - h \sum_{out} \dot{m}_{out}
$$ (constpressurereactor-energy)
