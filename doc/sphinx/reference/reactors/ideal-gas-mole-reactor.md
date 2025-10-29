```{py:currentmodule} cantera
```

# Ideal Gas Control Volume Mole Reactor

An ideal gas control volume mole reactor, as implemented by the C++ class
{ct}`IdealGasMoleReactor` and available in Python as the {py:class}`IdealGasMoleReactor`
class. It is defined by the state variables:

- $T$, the temperature (in K)
- $V$, the reactor volume (in m{sup}`3`)
- $n_k$, the number of moles for each species (in kmol)

Equations 1-3 are the governing equations for an ideal gas control volume mole reactor.

## Volume Equation

The reactor volume can change as a function of time due to the motion of one or more
[walls](sec-wall):

$$
\frac{dV}{dt} = \sum_w f_w A_w v_w(t)
$$ (ig-mole-reactor-volume)

Where $f_w = \pm 1$ indicates the facing of the wall (whether moving the wall increases
or decreases the volume of the reactor), $A_w$ is the surface area of the wall, and
$v_w(t)$ is the velocity of the wall as a function of time.

## Species Equations

The moles of each species in the reactor changes as a result of flow through the
reactor's [inlets and outlets](sec-flow-device), and production of homogeneous gas phase
species and reactions on the reactor [surfaces](sec-reactor-surface). The rate at which
species $k$ is generated through homogeneous phase reactions is $V \dot{\omega}_k$, and
the total rate at which moles of species $k$ changes is:

$$
\frac{dn_k}{dt} = V \dot{\omega}_k + \sum_\t{in} \dot{n}_{k, \t{in}}
                  - \sum_\t{out} \dot{n}_{k, \t{out}} + \dot{n}_{k, \t{wall}}
$$ (ig-mole-reactor-species)

## Energy Equation

In the case of the ideal gas control volume mole reactor model, the reactor temperature
$T$ is used instead of the total internal energy $U$ as a state variable. For an ideal
gas, we can rewrite the total internal energy in terms of the species moles and
temperature:

$$  U = \sum_k \hat{u}_k(T) n_k  $$

and differentiate it with respect to time to obtain:

$$  \frac{dU}{dt} = \frac{dT}{dt}\sum_k n_k \hat{c}_{v,k}
                    + \sum_k \hat{u}_k \frac{dn_k}{dt}
$$

Substituting the total energy equation for the general control volume mole reactor
{eq}`molereactor-energy` for $dU/dt$ yields an equation for the temperature:

$$
\sum_k n_k \hat{c}_{v,k} \frac{dT}{dt} = -p \frac{dV}{dt} + \dot{Q}
        + \sum_\t{in} \dot{n}_\t{in} \hat{h}_\t{in}
        - \hat{h} \sum_\t{out} \dot{n}_\t{out}
        - \sum_k \hat{u}_k \frac{dn_k}{dt}
$$

Substituting the species equation {eq}`ig-mole-reactor-species` for $dn_k/dt$ and making
some conversions between mass and moles gives the final form of the energy equation:

$$
m c_v \frac{dT}{dt} = & -p \frac{dV}{dt} + \dot{Q}
        - \sum_k \hat{u}_k \left(\dot{\omega}_k V + \dot{n}_{k, \t{wall}} \right) \\
        & - \frac{pV}{m} \sum_\t{out} \dot{m}_\t{out}
        + \sum_\t{in} \left(h_\t{in} \dot{m}_\t{in}
                            - \sum_k \hat{u}_k \dot{n}_{k,\t{in}} \right)
$$ (ig-mole-reactor-energy)
