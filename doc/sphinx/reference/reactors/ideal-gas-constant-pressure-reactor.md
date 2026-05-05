```{py:currentmodule} cantera
```

# Ideal Gas Constant Pressure Reactor

A constant pressure reactor using temperature as a state variable, as implemented by
the C++ class {ct}`IdealGasConstPressureReactor` and available in Python as the
{py:class}`IdealGasConstPressureReactor` class. It is defined by the state variables:

- $m$, the mass of the reactor's contents (in kg)
- $T$, the temperature (in K)
- $Y_k$, the mass fractions for each species (dimensionless)

Equations 1-3 below are the governing equations for this reactor model. While the
class name is historical, this formulation is applicable to non-ideal equations of
state as well.

## Mass Conservation

The total mass of the reactor's contents changes as a result of flow through the
reactor's [inlets and outlets](sec-flow-device), and production of homogeneous phase
species on [surfaces](sec-reactor-surface):

$$
\frac{dm}{dt} = \sum_\t{in} \dot{m}_\t{in} - \sum_\t{out} \dot{m}_\t{out}
              + \dot{m}_\t{wall}
$$ (igcpr-mass)

Where the subscripts *in* and *out* refer to the sum of the corresponding property over
all inlets and outlets respectively. A dot above a variable signifies a time derivative.

## Species Equations

The rate at which species $k$ is generated through homogeneous phase reactions is
$V \dot{\omega}_k W_k$, and the total rate at which species $k$ is generated is:

$$  \dot{m}_{k,\t{gen}} = V \dot{\omega}_k W_k + \dot{m}_{k,\t{wall}}  $$

The rate of change in the mass of each species is:

$$
\frac{d(mY_k)}{dt} = \sum_\t{in} \dot{m}_\t{in} Y_{k,\t{in}}
                   - \sum_\t{out} \dot{m}_\t{out} Y_k + \dot{m}_{k,gen}
$$

Expanding the derivative on the left hand side and substituting the equation
for $dm/dt$, the equation for each homogeneous phase species is:

$$
m \frac{dY_k}{dt} = \sum_\t{in} \dot{m}_\t{in} (Y_{k,\t{in}} - Y_k)
                  + \dot{m}_{k,\t{gen}} - Y_k \dot{m}_\t{wall}
$$ (igcpr-species)

## Energy Equation

In this reactor model, the reactor temperature $T$ is used as a state variable instead
of the total enthalpy $H$. For the mass-based form, write:

$$
H = H(T, P, m_1, \ldots, m_K), \qquad m_k = mY_k
$$

At constant pressure, applying the chain rule gives:

$$
\frac{dH}{dt} = m c_p \frac{dT}{dt}
    + \sum_k \frac{\bar{h}_k}{W_k} \frac{dm_k}{dt}
$$

where $\bar{h}_k$ are the partial molar enthalpies and $W_k$ are the molecular weights.
Substituting the species and mass equations into the constant
pressure reactor energy equation {eq}`constpressurereactor-energy` yields:

$$
m c_p \frac{dT}{dt} = \dot{Q}
     - \sum_k \frac{\bar{h}_k}{W_k} \dot{m}_{k,\t{gen}}
     + \sum_\t{in} \left(\dot{m}_\t{in} h_\t{in}
        - \sum_k \frac{\bar{h}_k}{W_k} \dot{m}_{k,\t{in}} \right)
$$ (igcpr-energy)
