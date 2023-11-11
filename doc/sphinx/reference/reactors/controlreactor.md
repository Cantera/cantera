```{py:currentmodule} cantera
```

# Control Volume Reactor

This model represents a homogeneous zero-dimensional reactor, as implemented by the C++
class {ct}`Reactor` and available in Python as the {py:class}`Reactor` class. A control
volume reactor is defined by the state variables:

- $m$, the mass of the reactor's contents (in kg)
- $V$, the reactor volume (in m$^3$)
- $U$, the total internal energy of the reactors contents (in J)
- $Y_k$, the mass fractions for each species (dimensionless)

Equations 1-4 below are the governing equations for a control volume reactor.

## Mass Conservation

The total mass of the reactor's contents changes as a result of flow through the
reactor's [inlets and outlets](sec-flow-device), and production of homogeneous phase
species on [surfaces](sec-reactor-surface):

$$
\frac{dm}{dt} = \sum_{in} \dot{m}_{in} - \sum_{out} \dot{m}_{out} + \dot{m}_{wall}
$$ (mass)

Where the subscripts *in* and *out* refer to the sum of the corresponding property over
all inlets and outlets respectively. A dot above a variable signifies a time derivative.

## Volume Equation

The reactor volume changes as a function of time due to the motion of one or more
[walls](sec-wall):

$$
\frac{dV}{dt} = \sum_w f_w A_w v_w(t)
$$ (volume)

where $f_w = \pm 1$ indicates the facing of the wall (whether moving the wall increases
or decreases the volume of the reactor), $A_w$ is the surface area of the wall, and
$v_w(t)$ is the velocity of the wall as a function of time.

## Species Equations

The rate at which species $k$ is generated through homogeneous phase reactions is
$V \dot{\omega}_k W_k$, and the total rate at which species $k$ is generated is:

$$
\dot{m}_{k,gen} = V \dot{\omega}_k W_k + \dot{m}_{k,wall}
$$

The rate of change in the mass of each species is:

$$
\frac{d(mY_k)}{dt} = \sum_{in} \dot{m}_{in} Y_{k,in} - \sum_{out} \dot{m}_{out} Y_k +
                     \dot{m}_{k,gen}
$$

Expanding the derivative on the left hand side and substituting the equation
for $dm/dt$, the equation for each homogeneous phase species is:

$$
m \frac{dY_k}{dt} = \sum_{in} \dot{m}_{in} (Y_{k,in} - Y_k) +
                    \dot{m}_{k,gen} - Y_k \dot{m}_{wall}
$$ (species)

## Energy Equation

The equation for the total internal energy is found by writing the first law for an open
system:

$$
\frac{dU}{dt} = - p \frac{dV}{dt} + \dot{Q} +
                \sum_{in} \dot{m}_{in} h_{in} - h \sum_{out} \dot{m}_{out}
$$ (cv-energy)

Where $\dot{Q}$ is the net rate of heat addition to the system.
