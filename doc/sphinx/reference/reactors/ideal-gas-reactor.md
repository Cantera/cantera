```{py:currentmodule} cantera
```

# Ideal Gas Control Volume Reactor

An ideal gas control volume reactor, as implemented by the C++ class
{ct}`IdealGasReactor` and available in Python as the {py:class}`IdealGasReactor` class.
It is defined by the state variables:

- $m$, the mass of the reactor's contents (in kg)
- $V$, the reactor volume (in m{sup}`3`)
- $T$, the temperature (in K)
- $Y_k$, the mass fractions for each species (dimensionless)

Equations 1-4 below are the governing equations for an ideal gas reactor.

## Mass Conservation

The total mass of the reactor's contents changes as a result of flow through the
reactor's [inlets and outlets](sec-flow-device), and production of gas phase species on
[surfaces](sec-reactor-surface):

$$
\frac{dm}{dt} = \sum_{in} \dot{m}_{in} - \sum_{out} \dot{m}_{out} + \dot{m}_{wall}
$$ (igr-mass)

where the subscripts *in* and *out* refer to the sum of the corresponding property over
all inlets and outlets respectively. A dot above a variable signifies a time derivative.

## Volume Equation

The reactor volume can change as a function of time due to the motion of one or more
[walls](sec-wall):

$$
\frac{dV}{dt} = \sum_w f_w A_w v_w(t)
$$ (igr-volume)

where $f_w = \pm 1$ indicates the facing of the wall (whether moving the wall increases
or decreases the volume of the reactor), $A_w$ is the surface area of the wall, and
$v_w(t)$ is the velocity of the wall as a function of time.

## Species Equations

The rate at which species $k$ is generated through homogeneous phase reactions is
$V \dot{\omega}_k W_k$, and the total rate at which species $k$ is generated is:

$$  \dot{m}_{k,gen} = V \dot{\omega}_k W_k + \dot{m}_{k,wall}  $$

The rate of change in the mass of each species is:

$$
\frac{d(mY_k)}{dt} = \sum_{in} \dot{m}_{in} Y_{k,in} - \sum_{out} \dot{m}_{out} Y_k +
                     \dot{m}_{k,gen}
$$

Expanding the derivative on the left hand side and substituting the equation
for $dm/dt$, the equation for each homogeneous phase species is:

$$
m \frac{dY_k}{dt} = \sum_{in} \dot{m}_{in} (Y_{k,in} - Y_k)+ \dot{m}_{k,gen}
                    - Y_k \dot{m}_{wall}
$$ (igr-species)

## Energy Equation

In the case of the ideal gas control volume reactor model, the reactor temperature $T$
is used instead of the total internal energy $U$ as a state variable. For an ideal gas,
we can rewrite the total internal energy in terms of the mass fractions and temperature:

$$  U = m \sum_k Y_k u_k(T)  $$

and differentiate it with respect to time to obtain:

$$
\frac{dU}{dt} = u \frac{dm}{dt} + m c_v \frac{dT}{dt} + m \sum_k u_k \frac{dY_k}{dt}
$$

Substituting this into the energy equation for the control volume reactor
{eq}`cv-energy` yields an equation for the temperature:

$$
m c_v \frac{dT}{dt} =& - p \frac{dV}{dt} + \dot{Q} + \sum_{in} \dot{m}_{in} \left( h_{in} - \sum_k u_k Y_{k,in} \right) \\
    &- \frac{p V}{m} \sum_{out} \dot{m}_{out} - \sum_k \dot{m}_{k,gen} u_k
$$ (igr-energy)

While this form of the energy equation is somewhat more complicated, it significantly
reduces the cost of evaluating the system Jacobian, since the derivatives of the species
equations are taken at constant temperature instead of constant internal energy.
