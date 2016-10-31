.. default-role:: math

.. py:currentmodule:: cantera

**********************
One-Dimensional Flames
**********************

Cantera includes a set of models for representing steady-state, quasi-one-
dimensional reacting flows, which can be used to simulate a number of common
flames, such as:

- freely-propagating premixed laminar flames
- burner-stabilized premixed flames
- counterflow diffusion flames
- counterflow (strained) premixed flames

Additional capabilities include simulation of surface reactions, which can be
used to represent processes such as combustion on a catalytic surface or
chemical vapor deposition processes.

All of these configurations are simulated using a common set of governing
equations within a 1D "flow" domain, with the differences between the models
being represented by differences in the boundary conditions applied. Here, we
describe the governing equations and the various boundary conditions which can
be applied.

Stagnation Flow Governing Equations
===================================

Cantera models flames which are stabilized in an axisymmetric stagnation flow,
and computes the solution along the stagnation streamline (`r=0`), using a
similarity solution to reduce the three-dimensional governing equations to a
single dimension.

The governing equations for a steady axisymmetric stagnation flow follow those
derived in Section 6.2 of [KCG2003]_:

*Continuity*:

.. math::

    \frac{\partial\rho u}{\partial z} + 2 \rho V = 0

*Radial momentum*:

.. math::

    \rho u \frac{\partial V}{\partial z} + \rho V^2 =
        - \Lambda
        + \frac{\partial}{\partial z}\left(\mu \frac{\partial V}{\partial z}\right)


*Energy*:

.. math::

    \rho c_p u \frac{\partial T}{\partial z} =
        \frac{\partial}{\partial z}\left(\lambda \frac{\partial T}{\partial z}\right)
        - \sum_k j_k c_{p,k} \frac{\partial T}{\partial z}
        - \sum_k h_k W_k \dot{\omega}_k

*Species*:

.. math::

    \rho u \frac{\partial Y_k}{\partial z} = - \frac{\partial j_k}{\partial z}
        + W_k \dot{\omega}_k

where `\rho` is the density, `u` is the axial velocity, `v` is the radial
velocity, `V = v/r` is the scaled radial velocity, `\Lambda` is the pressure
eigenvalue (independent of `z`), `\mu` is the dynamic viscosity, `c_p` is the
heat capacity at constant pressure, `T` is the temperature, `\lambda` is the
thermal conductivity, `Y_k` is the mass fraction of species `k`, `j_k` is the
diffusive mass flux of species `k`, `c_{p,k}` is the specific heat capacity of
species `k`, `h_k` is the enthalpy of species `k`, `W_k` is the molecular weight
of species `k`, and `\dot{\omega}_k` is the molar production rate of species
`k`.

The tangential velocity `w` has been assumed to be zero, and the fluid has been
assumed to behave as an ideal gas.

To help in the solution of the discretized problem, it is convenient to write a
differential equation for the scalar `\Lambda`:

.. math::

    \frac{d\Lambda}{dz} = 0

Diffusive Fluxes
----------------

The species diffusive mass fluxes `j_k` are computed according to either a
mixture-averaged or multicomponent formulation. If the mixture-averaged
formulation is used, the calculation performed is:

.. math::

    j_k^* = \rho \frac{W_k}{\overline{W}} D_{k,m} \frac{\partial X_k}{\partial z}

    j_k = j_k^* - Y_k \sum_i j_i^*

where `\overline{W}` is the mean molecular weight of the mixture, `D_{k,m}` is the
mixture-averaged diffusion coefficient for species `k`, and `X_k` is the mole
fraction for species `k`. The diffusion coefficients used here are those
computed by the method :ct:`GasTransport::getMixDiffCoeffs`. The correction
applied by the second equation ensures that the sum of the mass fluxes is zero,
a condition which is not inherently guaranteed by the mixture-averaged
formulation.

When using the multicomponent formulation, the mass fluxes are computed
according to:

.. math::

    j_k = \frac{\rho W_k}{\overline{W}^2} \sum_i W_i D_{ki} \frac{\partial X_i}{\partial z}
          - \frac{D_k^T}{T} \frac{\partial T}{\partial z}

where `D_{ki}` is the multicomponent diffusion coefficient and `D_k^T` is the
Soret diffusion coefficient (used only if calculation of this term is
specifically enabled).

Boundary Conditions
===================

Inlet boundary
--------------

For a boundary located at a point `z_0` where there is an inflow, values are
supplied for the temperature `T_0`, the species mass fractions `Y_{k,0}` the
scaled radial velocity `V_0`, and the mass flow rate `\dot{m}_0` (except in the
case of the freely-propagating flame).

The following equations are solved at the point `z = z_0`:

.. math::

    T(z_0) = T_0

    V(z_0) = V_0

    \dot{m}_0 Y_{k,0} - j_k(z_0) - \rho(z_0) u(z_0) Y_k(z_0) = 0

If the mass flow rate is specified, we also solve:

.. math::

    \rho(z_0) u(z_0) = \dot{m}_0

Otherwise, we solve:

.. math::

    \Lambda(z_0) = 0

Outlet boundary
---------------

For a boundary located at a point `z_0` where there is an outflow, we solve:

.. math::

    \Lambda(z_0) = 0

    \left.\frac{\partial T}{\partial z}\right|_{z_0} = 0

    \left.\frac{\partial Y_k}{\partial z}\right|_{z_0} = 0

    V(z_0) = 0


Symmetry boundary
-----------------

For a symmetry boundary located at a point `z_0`, we solve:

.. math::

    \rho(z_0) u(z_0) = 0

    \left.\frac{\partial V}{\partial z}\right|_{z_0} = 0

    \left.\frac{\partial T}{\partial z}\right|_{z_0} = 0

    j_k(z_0) = 0

Reacting surface
----------------

For a surface boundary located at a point `z_0` on which reactions may occur,
the temperature `T_0` is specified. We solve:

.. math::

    \rho(z_0) u(z_0) = 0

    V(z_0) = 0

    T(z_0) = T_0

    j_k(z_0) + \dot{s}_k W_k = 0

where `\dot{s}_k` is the molar production rate of the gas-phase species `k` on
the surface. In addition, the surface coverages `\theta_i` for each surface
species `i` are computed such that `\dot{s}_i = 0`.


References
==========

.. [KCG2003] Kee, Coltrin, Glarborg: *Chemically Reacting Flow*.
             Wiley-Interscience, 2003

