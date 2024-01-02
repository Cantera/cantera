```{py:currentmodule} cantera
```

# Control Volume Mole Reactor

A control volume mole reactor, as implemented by the C++ class {ct}`MoleReactor` and
available in Python as the {py:class}`MoleReactor` class. It is defined by the state
variables:

- $U$, the total internal energy of the reactor's contents (in J)
- $V$, the reactor volume (in m{sup}`3`)
- $n_k$, the number of moles for each species (in kmol)

Equations 1-3 are the governing equations for a control volume mole reactor.

## Volume Equation

The reactor volume changes as a function of time due to the motion of one or
more [walls](sec-wall):

$$
\frac{dV}{dt} = \sum_w f_w A_w v_w(t)
$$ (molereactor-volume)

where $f_w = \pm 1$ indicates the facing of the wall (whether moving the wall increases
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
$$ (molereactor-species)

where the subscripts *in* and *out* refer to the sum of the corresponding property over
all inlets and outlets respectively. A dot above a variable signifies a time derivative.

## Energy Equation

The equation for the total internal energy is found by writing the first law for an open
system:

$$
\frac{dU}{dt} = - p \frac{dV}{dt} + \dot{Q} + \sum_\t{in} \dot{n}_\t{in} \hat{h}_\t{in}
                - \hat{h} \sum_\t{out} \dot{n}_\t{out}
$$ (molereactor-energy)

where $\dot{Q}$ is the net rate of heat addition to the system and $\hat{h}$ is the
molar enthalpy.
