```{py:currentmodule} cantera
```

# Ideal Gas Constant Pressure Mole Reactor

A constant pressure mole reactor using temperature as a state variable, as implemented
by the C++ class {ct}`IdealGasConstPressureMoleReactor` and available in Python as the
{py:class}`IdealGasConstPressureMoleReactor` class. It is defined by the state
variables:

- $T$, the temperature (in K)
- $n_k$, the number of moles for each species (in kmol)

Equations 1 and 2 below are the governing equations for this reactor model. While the
class name is historical, this formulation is applicable to non-ideal equations of
state as well.

## Species Equations

The moles of each species in the reactor changes as a result of flow through the
reactor's [inlets and outlets](sec-flow-device), and production of homogeneous gas phase
species and reactions on the reactor [surfaces](sec-reactor-surface). The rate at which
species $k$ is generated through homogeneous phase reactions is $V \dot{\omega}_k$, and
the total rate at which moles of species $k$ changes is:

$$
\frac{dn_k}{dt} = V \dot{\omega}_k + \sum_\t{in} \dot{n}_{k, \t{in}}
                  - \sum_\t{out} \dot{n}_{k, \t{out}} + \dot{n}_{k, \t{wall}}
$$ (ig-const-pressure-mole-reactor-species)

Where the subscripts *in* and *out* refer to the sum of the corresponding property over
all inlets and outlets respectively. A dot above a variable signifies a time derivative.

## Energy Equation

Writing the first law for an open system gives:

$$
\frac{dU}{dt} = - p \frac{dV}{dt} + \dot{Q} + \sum_\t{in} \hat{h}_\t{in} \dot{n}_\t{in}
              - \hat{h} \sum_\t{out} \dot{n}_\t{out}
$$

where positive $\dot{Q}$ represents heat addition to the system and $\hat{h}$ is the
molar enthalpy of the reactor's contents.

Differentiating the definition of the total enthalpy, $H = U + pV$, with respect to time
gives:

$$  \frac{dH}{dt} = \frac{dU}{dt} + p \frac{dV}{dt} + V \frac{dp}{dt}  $$

Noting that $dp/dt = 0$ and substituting into the energy equation yields:

$$
\frac{dH}{dt} = \dot{Q} + \sum_\t{in} \hat{h}_\t{in} \dot{n}_\t{in}
                - \hat{h} \sum_\t{out} \dot{n}_\t{out}
$$

In this reactor model, the reactor temperature $T$ is used as a state variable instead
of the total enthalpy $H$. For a general equation of state, write:

$$
H = H(T, P, n_1, \ldots, n_K)
$$

At constant pressure, applying the chain rule gives:

$$
\frac{dH}{dt} = N \hat{c}_p \frac{dT}{dt}
    + \sum_k \bar{h}_k \frac{dn_k}{dt}
$$

where $N$ is the total number of moles, $\hat{c}_p$ is the mixture molar heat
capacity, and $\bar{h}_k$ are the partial molar enthalpies. Making this substitution
yields:

$$
N \hat{c}_p \frac{dT}{dt} = & \dot{Q}
        + \sum_\t{in} \dot{n}_\t{in} \hat{h}_\t{in}
        - \hat{h} \sum_\t{out} \dot{n}_\t{out} \\
        & - \sum_k \bar{h}_k \left(V \dot{\omega}_k + \dot{n}_{k,\t{wall}}
                                   + \sum_\t{in} \dot{n}_{k,\t{in}}
                                   - \sum_\t{out} \dot{n}_{k,\t{out}} \right)
$$

Rearranging and simplifying gives the final energy equation:

$$
N \hat{c}_p \frac{dT}{dt} = \dot{Q}
        - \sum_k \bar{h}_k \left(V \dot{\omega}_k + \dot{n}_{k,\t{wall}} \right)
        + \sum_\t{in} \left(\dot{n}_\t{in} \hat{h}_\t{in}
            - \sum_k \bar{h}_k \dot{n}_{k,\t{in}} \right)
$$  (ig-const-pressure-mole-reactor-energy)
