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

Equations 1-4 below are the governing equations for this reactor model. While the class
name is historical, this formulation is now applicable to non-ideal equations of state
as well.

## Mass Conservation

The total mass of the reactor's contents changes as a result of flow through the
reactor's [inlets and outlets](sec-flow-device), and production of gas phase species on
[surfaces](sec-reactor-surface):

$$
\frac{dm}{dt} = \sum_\t{in} \dot{m}_\t{in} - \sum_\t{out} \dot{m}_\t{out}
              + \dot{m}_\t{wall}
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

$$  \dot{m}_{k,\t{gen}} = V \dot{\omega}_k W_k + \dot{m}_{k,\t{wall}}  $$

The rate of change in the mass of each species is:

$$
\frac{d(mY_k)}{dt} = \sum_\t{in} \dot{m}_\t{in} Y_{k,\t{in}} - \sum_\t{out} \dot{m}_\t{out} Y_k +
                     \dot{m}_{k,\t{gen}}
$$

Expanding the derivative on the left hand side and substituting the equation
for $dm/dt$, the equation for each homogeneous phase species is:

$$
m \frac{dY_k}{dt} = \sum_\t{in} \dot{m}_\t{in} (Y_{k,\t{in}} - Y_k)+ \dot{m}_{k,\t{gen}}
                    - Y_k \dot{m}_\t{wall}
$$ (igr-species)

## Energy Equation

In this reactor model, the reactor temperature $T$ is used as a state variable instead
of the total internal energy $U$. For the mass-based form, write:

$$
U = U(T, V, m_1, \ldots, m_K), \qquad m_k = mY_k
$$

and apply the chain rule:

$$
\frac{dU}{dt} =
    \frac{\partial U}{\partial T}\Bigg|_{V,m_k}\frac{dT}{dt}
    + \frac{\partial U}{\partial V}\Bigg|_{T,m_k}\frac{dV}{dt}
    + \sum_k \frac{\partial U}{\partial m_k}\Bigg|_{T,V,m_{j\ne k}}
    \frac{dm_k}{dt}
$$

Here we can make the substitutions $(\partial U / \partial T)_{V,m_k} = m c_v$;
$\pi_T \equiv (\partial U / \partial V)_{T,m_k}$ is the internal pressure (see
; and
$(\partial U / \partial m_k)_{T,V,m_{j\ne k}} = \tilde{u}_k /W_k $ where $\tilde{u}_k$
are the partial molar internal energies at constant temperature and volume (see
). This gives:

$$
\frac{dU}{dt} = m c_v \frac{dT}{dt} + \pi_T \frac{dV}{dt}
    + \sum_k \frac{\tilde{u}_k}{W_k}\frac{dm_k}{dt}
$$

Substituting into the control-volume energy equation {eq}`cv-energy` gives:

$$
m c_v \frac{dT}{dt} =& -(p + \pi_T) \frac{dV}{dt} + \dot{Q}
    - \sum_k \dot{m}_{k,\t{gen}} \frac{\tilde{u}_k}{W_k}
    - \frac{p V}{m} \sum_\t{out} \dot{m}_\t{out} \\
    &+ \sum_\t{in} \dot{m}_\t{in}
       \left( h_\t{in} - \sum_k Y_{k,\t{in}} \frac{\tilde{u}_k}{W_k} \right)
$$ (igr-energy)

While this form of the energy equation is somewhat more complicated, it significantly
reduces the cost of evaluating the system Jacobian, since the derivatives of the species
equations are taken at constant temperature instead of constant internal energy. It
also eliminates the need to iteratively solve the equation of state to find $T$
given $u$ and $v$.

In the case of an ideal gas, $\pi_T = 0$ and $\tilde{u}_k$ are equal to the usual
partial molar internal energies $\hat{u}_k$.

```{seealso}
- {ct}`ThermoPhase::internalPressure`
- {ct}`ThermoPhase::getPartialMolarIntEnergies_TV`
```