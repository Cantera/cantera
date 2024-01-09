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
