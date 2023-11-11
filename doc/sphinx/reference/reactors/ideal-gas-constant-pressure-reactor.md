```{py:currentmodule} cantera
```

# Ideal Gas Constant Pressure Reactor

An ideal gas constant pressure reactor, as implemented by the C++ class
{ct}`IdealGasConstPressureReactor` and available in Python as the
{py:class}`IdealGasConstPressureReactor` class. It is defined by the state variables:

- $m$, the mass of the reactor's contents (in kg)
- $T$, the temperature (in K)
- $Y_k$, the mass fractions for each species (dimensionless)

Equations 1-3 below are the governing equations for an ideal gas constant pressure
reactor.

## Mass Conservation

The total mass of the reactor's contents changes as a result of flow through the
reactor's [inlets and outlets](sec-flow-device), and production of homogeneous phase
species on [surfaces](sec-reactor-surface):

$$
\frac{dm}{dt} = \sum_{in} \dot{m}_{in} - \sum_{out} \dot{m}_{out} + \dot{m}_{wall}
$$ (igcpr-mass)

Where the subscripts *in* and *out* refer to the sum of the corresponding property over
all inlets and outlets respectively. A dot above a variable signifies a time derivative.

## Species Equations

The rate at which species $k$ is generated through homogeneous phase reactions is
$V \dot{\omega}_k W_k$, and the total rate at which species $k$ is generated is:

$$  \dot{m}_{k,gen} = V \dot{\omega}_k W_k + \dot{m}_{k,wall}  $$

The rate of change in the mass of each species is:

$$
\frac{d(mY_k)}{dt} = \sum_{in} \dot{m}_{in} Y_{k,in} - \sum_{out} \dot{m}_{out} Y_k
                     + \dot{m}_{k,gen}
$$

Expanding the derivative on the left hand side and substituting the equation
for $dm/dt$, the equation for each homogeneous phase species is:

$$
m \frac{dY_k}{dt} = \sum_{in} \dot{m}_{in} (Y_{k,in} - Y_k) + \dot{m}_{k,gen}
                    - Y_k \dot{m}_{wall}
$$ (igcpr-species)

## Energy Equation

As for the [ideal gas reactor](ideal-gas-reactor), we replace the total enthalpy as a
state variable with the temperature by writing the total enthalpy in terms of the mass
fractions and temperature and differentiating with respect to time:

$$
H &= m \sum_k Y_k h_k(T)

\frac{dH}{dt} &= h \frac{dm}{dt} + m c_p \frac{dT}{dt}
                + m \sum_k h_k \frac{dY_k}{dt}
$$

Substituting the corresponding derivatives into the constant pressure reactor energy
equation {eq}`constpressurereactor-energy` yields an equation for the temperature:

$$
m c_p \frac{dT}{dt} = \dot{Q} - \sum_k h_k \dot{m}_{k,gen}
     + \sum_{in} \dot{m}_{in} \left(h_{in} - \sum_k h_k Y_{k,in} \right)
$$ (igcpr-energy)
