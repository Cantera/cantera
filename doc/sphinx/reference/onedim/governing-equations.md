# Governing Equations for One-dimensional Flow

Cantera models flames that are stabilized in an axisymmetric stagnation flow, and
computes the solution along the stagnation streamline ($r=0$), using a similarity
solution to reduce the three-dimensional governing equations to a single dimension.

## Axisymmetric Stagnation Flow

The governing equations for a steady axisymmetric stagnation flow follow those derived
in Section 7.2 of {cite:t}`kee2017` and are implemented by class {ct}`StFlow`.

*Continuity*:

$$  \pxpy{u}{z} + 2 \rho V = 0  $$

*Radial momentum*:

$$
\rho u \pxpy{V}{z} + \rho V^2 = - \Lambda + \pxpy{}{z}\left(\mu \pxpy{V}{z}\right)
$$

*Energy*:

$$
\rho c_p u \pxpy{T}{z} = \pxpy{}{z}\left(\lambda \pxpy{T}{z}\right)
    - \sum_k j_k \pxpy{h_k}{z} - \sum_k h_k W_k \dot{\omega}_k
$$

*Species*:

$$
\rho u \pxpy{Y_k}{z} = - \pxpy{j_k}{z} + W_k \dot{\omega}_k
$$

where the following variables are used:

- $z$ is the axial coordinate
- $r$ is the radial coordinate
- $\rho$ is the density
- $u$ is the axial velocity
- $v$ is the radial velocity
- $V = v/r$ is the scaled radial velocity
- $\Lambda$ is the pressure eigenvalue (independent of $z$)
- $\mu$ is the dynamic viscosity
- $c_p$ is the heat capacity at constant pressure
- $T$ is the temperature
- $\lambda$ is the thermal conductivity
- $Y_k$ is the mass fraction of species $k$
- $j_k$ is the diffusive mass flux of species $k$
- $c_{p,k}$ is the specific heat capacity of species $k$
- $h_k$ is the enthalpy of species $k$
- $W_k$ is the molecular weight of species $k$
- $\dot{\omega}_k$ is the molar production rate of species $k$.

The tangential velocity $w$ has been assumed to be zero. The model is applicable to both
ideal and non-ideal fluids, which follow ideal-gas or real-gas (Redlich-Kwong and
Peng-Robinson) equations of state.

```{versionadded} 3.0
Support for real gases in the flame models was introduced in Cantera 3.0.
```

To help in the solution of the discretized problem, it is useful to write a
differential equation for the scalar $\Lambda$:

$$  \frac{d\Lambda}{dz} = 0  $$

When discretized, the Jacobian terms introduced by this equation match the block
diagonal structure produced by the other governing equations, rather than creating a
column of entries that would cause fill-in when factorizing as part of the Newton
solver.

## Diffusive Fluxes

The species diffusive mass fluxes $j_k$ are computed according to either a
mixture-averaged or multicomponent formulation. If the mixture-averaged formulation is
used, the calculation performed is:

$$
j_k^* = - \rho \frac{W_k}{\overline{W}} D_{km}^\prime \pxpy{X_k}{z}

j_k = j_k^* - Y_k \sum_i j_i^*
$$

where $\overline{W}$ is the mean molecular weight of the mixture, $D_{km}^\prime$ is the
mixture-averaged diffusion coefficient for species $k$, and $X_k$ is the mole fraction
for species $k$. The diffusion coefficients used here are those computed by the method
{ct}`GasTransport::getMixDiffCoeffs`. The correction applied by the second equation
ensures that the sum of the mass fluxes is zero, a condition which is not inherently
guaranteed by the mixture-averaged formulation.

When using the multicomponent formulation, the mass fluxes are computed according to:

$$
j_k = \frac{\rho W_k}{\overline{W}^2} \sum_i W_i D_{ki} \pxpy{X_i}{z}
      - \frac{D_k^T}{T} \pxpy{T}{z}
$$

where $D_{ki}$ is the multicomponent diffusion coefficient and $D_k^T$ is the Soret
diffusion coefficient. Inclusion of the Soret calculation must be explicitly enabled
when setting up the simulation, on top of specifying a multicomponent transport model,
for example by using the {ct}`StFlow::enableSoret` method (C++) or setting the
{py:attr}`~cantera.FlameBase.soret_enabled` property (Python).

## Boundary Conditions

### Inlet boundary

For a boundary located at a point $z_0$ where there is an inflow, values are supplied
for the temperature $T_0$, the species mass fractions $Y_{k,0}$ the scaled radial
velocity $V_0$, and the mass flow rate $\dot{m}_0$. In the case of the
freely-propagating flame, the mass flow rate is not an input but is determined
indirectly by holding the temperature fixed at an intermediate location within the
domain; see [](discretization) for details.

The following equations are solved at the point $z = z_\t{in}$:

$$
T(z_\t{in}) &= T_0

V(z_\t{in}) &= V_0

\dot{m}_0 Y_{k,\t{in}} - j_k(z_\t{in}) - \rho(z_\t{in}) u(z_\t{in}) Y_k(z_\t{in}) &= 0
$$

If the mass flow rate is specified, we also solve:

$$
\rho(z_\t{in}) u(z_\t{in}) = \dot{m}_0
$$

Otherwise, we solve:

$$  \Lambda(z_\t{in}) = 0  $$

These equations are implemented by class {ct}`Inlet1D`.

### Outlet boundary

For a boundary located at a point $z_\t{out}$ where there is an outflow, we solve:

$$
\Lambda(z_\t{out}) = 0

\left.\pxpy{T}{z}\right|_{z_\t{out}} = 0

\left.\pxpy{Y_k}{z}\right|_{z_\t{out}} = 0

V(z_\t{out}) = 0
$$

These equations are implemented by class {ct}`Outlet1D`.

### Symmetry boundary

For a symmetry boundary located at a point $z_\t{symm}$, we solve:

$$
\rho(z_\t{symm}) u(z_\t{symm}) = 0

\left.\pxpy{V}{z}\right|_{z_\t{symm}} = 0

\left.\pxpy{T}{z}\right|_{z_\t{symm}} = 0

j_k(z_\t{symm}) = 0
$$

These equations are implemented by class {ct}`Symm1D`.

### Reacting surface

For a surface boundary located at a point $z_\t{surf}$ on which reactions may
occur, the temperature $T_\t{surf}$ is specified. We solve:

$$
\rho(z_\t{surf}) u(z_\t{surf}) &= 0

V(z_\t{surf}) &= 0

T(z_\t{surf}) &= T_\t{surf}

j_k(z_\t{surf}) + \dot{s}_k W_k &= 0
$$

where $\dot{s}_k$ is the molar production rate of the gas-phase species $k$ on the
surface. In addition, the surface coverages $\theta_i$ for each surface species $i$ are
computed such that $\dot{s}_i = 0$.

These equations are implemented by class {ct}`ReactingSurf1D`.

## The Drift-Diffusion Model

To account for the transport of charged species in a flame, class {ct}`IonFlow` adds the
drift term to the diffusive fluxes of the mixture-average formulation according to
{cite:t}`pedersen1993`,

$$
j_k^* = \rho \frac{W_k}{\overline{W}} D_{km}^\prime \pxpy{X_k}{z} + s_k \mu_k E Y_k,
$$

where $s_k$ is the sign of charge (1,-1, and 0 respectively for positive, negative, and
neutral charge), $\mu_k$ is the mobility, and $E$ is the electric field. The diffusion
coefficients and mobilities of charged species can be more accurately calculated by
{ct}`IonGasTransport::getMixDiffCoeffs` and {ct}`IonGasTransport::getMobilities`. The
following correction is applied instead to preserve the correct fluxes of charged
species:

$$
j_k = j_k^* - \frac {1 - |s_k|} {1 - \sum_i |s_i| Y_i} Y_k \sum_i j_i^*.
$$

In addition, Gauss's law is solved simultaneously with the species and energy equations,

$$
\pxpy{E}{z} &= \frac{e}{\epsilon_0}\sum_k Z_k n_k ,

n_k &= N_a \rho Y_k / W_k,

E|_{z=0} &= 0,
$$

where $Z_k$ is the charge number, $n_k$ is the number density, and $N_a$ is the Avogadro
number.

## Counterflow Two-Point Flame Control

A two-point temperature control feature is available for counterflow diffusion flames. This
feature allows users to set a control points on both sides of a flame and incrementally lower
the flame temperature. This allows for the simulation of the stable burning branch as well
as the unstable burning branch of the standard flamelet "S-curve". The implementation is based
on the method discussed in {cite:t}`nishioka1996` and {cite:t}`huo2014`. The diagram below shows
the general concept of the two-point flame control method, with control points located on either
side of the peak flame temperature. An initial flame solution is used as a starting point, and
the temperatures at the control points are lowered to produce a new flame solution that satisfies
the governing equations and passes through the new temperatures at the control points.

```{image} /_static/images/two_point_control_diagram.svg
:width: 50%
:alt: Two-Point Flame Control Diagram
:align: center
```

For the two-point control method, one governing equation was modified ($\Lambda$), and a new
governing equation for the axial oxidizer velocity was added ($U_o$). The fuel and oxidizer velocity boundary
conditions are modified when the two-point control is active. These equations allow
for the temperature reduction to be performed in a numerically consistent manner (preventing
any issues of over-defining the system of governing equations). The two equations that are activated
when two-point control is turned on are:

$$
\frac{d\Lambda}{dz} = 0
$$

and

$$
\frac{dU_o}{dz} = 0
$$

These equations are zero everywhere, except at their respective control points. At the left control point
the residual for the $\Lambda$ equation is:

$$
residual = T(z=Z_L) - T_{L, control}
$$

At the left control point the residual for the $U_o$ equation is:

$$
residual = T(z=Z_R) - T_{R, control}
$$

Where T(z=Z_L) is the temperature of the flowfield at the left control point, T(z=Z_R) is the temperature of
the flowfield at the right control point, T_{L, control} is the left control point desired temperature, and
T_{R, control} is the right control point desired temperature.


The values of $\Lambda$ and $U_o$ are influenced by the left and right control points, respectively.
A residual error is induced because of the difference between the flow's temperature at that point and
the desired control point temperature.
In order to drive this error to zero, the solver adjusts the flow rates at the boundaries, which changes
the temperature distribution, which in turn affects the values of $\Lambda$ and $U_o$.

At the right boundary, the boundary condition for the continuity equation is imposed by using
the solution from the oxidizer velocity equation. At the left boundary, the boundary condition for
the continuity equation is imposed by using the value of the axial velocity at the left boundary.

$$
mdot(z=0) = rho(z=0)*U(z=0)
$$

$$
mdot(z=L) = rho(z=L)*U_0(z=L)
$$


