# Discretization of 1D Equations

# Governing Equations

This section outlines the discretizations used for the governing equations, including the continuity, radial momentum, energy, and species equations, see(governingEquationsSection) for more details on the governing equations. The discretization methods applied to these equations include upwinding for terms involving \( u \), central differences for second derivative terms, and one-sided discretizations for boundary conditions.

The domain can be thought of in terms of the following diagram.

```{image} /_static/images/1d_domain_diagram_1.svg
:width: 50%
:alt: Diagram of the 1D domain
:align: center
```

The 1D solver sweeps over all points in a domain and computes a residual for each grid point and each equation. This residual vector has entries
in a form similar to what is shown below.

$$
F(x) = F(u_0, V_0, T_0, \lambda_0, u_1, V_1, T_1, \lambda_1, ..., u_{N-1}, V_{N-1}, T_{N-1}, \lambda_{N-1}, u_{N}, V_{N}, T_{N}, \lambda_{N})
$$

For each governing equation, there is a default boundary condition applied at the left and right boundaries. During the sweep, the last portions that
are evaluated are the actual domain boundary points, which allows for the default residual to be augmented or overwritten by a boundary condition that is
different from the default boundary condition.

## Continuity Equation

The continuity equation is given by:

$$
 \frac{\partial (\rho u)}{ \partial z} + 2\rho V = 0
$$

**Discretization:**
The continuity equation does not need any upwinding in the standard sense, but because the equation is first order and the
boundary condition is applied at the right boundary in the 1D domain, the first derivative approximation uses a one-sided
form. For consistency, the radial momentum term also uses the same points and uses an average value between two nodes.


For the interior points in the domain, the discretized equation in residual form (all terms moved to one side) is:
$$
F = \frac{\rho_j u_j - \rho_{j+1} u_{j+1}}{z_{j+1} - z_j} + 2 \left( \frac{\rho_j V_j + \rho_{j+1} V_{j+1}}{2} \right)
$$

**Boundary Conditions:**
At the right boundary, the default boundary condition is a zero velocity condition, which is representative of a stagnation
surface. The radial momentum at the stagnation surface is also zero $ V=0 $.

At the right boundary(`j=N`)
$$
F = \rho_j u_j - \dot m
$$

Expressing the boundary residual in this form will cause the Newton root finding algorithm to find the value of $ \rho_j u_j $ that
minimizes the difference between it and the imposed(say from a user-defined boundary mass flow rate specification) boundary condition.

There is no imposed boundary condition at the left boundary, and so the residual equation at the left boundary is the same as the interior points.

At the left boundary(`j=0`)
$$
F = \frac{\rho_j u_j - \rho_{j+1} u_{j+1}}{z_{j+1} - z_j} + 2 \left( \frac{\rho_j V_j + \rho_{j+1} V_{j+1}}{2} \right)
$$


## Radial Momentum Equation

The radial momentum equation is:

$$
\rho u \frac{\partial V}{\partial z} + \rho V^2 = -\Lambda + \frac{\partial}{\partial z} \left( \mu \frac{\partial V}{\partial z} \right)
$$

**Discretization:**
- The term $ \rho u \frac{\partial V}{\partial z} $ uses upwinding.
- The second derivative term $ \frac{\partial}{\partial z} \left( \mu \frac{\partial V}{\partial z} \right) $ uses central differences.

In general, the upwinding formula for a variable, A, is:

$$
\frac{\partial A}{\partial z} \bigg|_{j} \approx \frac{A(x, j_{\text{loc}}) -
  A(x, j_{\text{loc} - 1})}{z_{jloc} - z_{jloc-1}}
$$

Where the value of `loc` is determined by the sign of the axial velocity. If the axial velocity is positive, the value of `loc` is j. If the axial velocity
is negative, the value of `loc` is `j+1`. A positive velocity means that the flow is moving left-to-right.


For the second derivative term (shear term in the momentum equation), a three-point central difference formula is used. The term we are discretizing is:

$$
\frac{d}{dz}\left(\mu \frac{dV}{dz}\right)
$$

Let $ A = \mu \frac{dV}{dz} $ for simplicity. We'll call this the inner term or inner discretization. In this situation, the inner term is evaluated
using a central difference formula, but evaluated at `j+1/2` and `j-1/2` (halfway between the grid points around point j).

The value of $ A $ at point $ j-1/2 $ is estimating using a central difference formula:

$$
A_{j-1/2} = \mu_{j-1/2} \frac{V_j - V_{j-1}}{z_j - z_{j-1}}
$$

The value of $ A $ at point $ j+1/2 $ is:

$$
A_{j+1/2} = \mu_{j+1/2} \frac{V_{j+1} - V_j}{z_{j+1} - z_j}
$$

$ \mu_{j+1/2} $ is the viscosity, estimated at the midpoint between grid points.

The outer discretization uses a central difference between the `j+1/2` and `j-1/2`
locations.

$$
\frac{dA}{dz} \approx \frac{A_{j+1/2} - A_{j-1/2}}{z_{j+1/2} - z_{j-1/2}}
$$

Where the values of $ z $ are:
$ z_{j+1/2} = z_{j} + 0.5*(z_{j+1} - z_j) = 0.5*(z_{j} + z_{j+1}) $ and
$ z_{j-1/2} = z_j - 0.5*(z_j - z_{j-1}) = 0.5*(z_{j} + z_{j-1}) $. The
difference between these two values is
$ z_{j+1/2} - z_{j-1/2} = \frac{z_{j+1} - z_{j-1}}{2} $.

Substituting these values into the central difference formula gives:

$$
\frac{d}{dz}\left(\mu \frac{dV}{dz}\right) \approx \frac{\mu_{j+1/2} \frac{V_{j+1} - V_j}{z_{j+1} - z_j} -
  \mu_{j-1/2} \frac{V_j - V_{j-1}}{z_j - z_{j-1}}}{\frac{z_{j+1} - z_{j-1}}{2}}
$$


## Energy Equation

The steady-state energy equation is described by:

$$
\rho c_p u \frac{\partial T}{\partial z} = \frac{\partial}{\partial z} \left( \lambda \frac{\partial T}{\partial z} \right) - \sum_k j_k \frac{\partial h_k}{\partial z} - \sum_k h_k W_k \dot{\omega}_k
$$

**Discretization:**
- The term $ \rho c_p u \frac{\partial T}{\partial z} is discretized using upwinding.
- The second derivative term $ \frac{\partial}{\partial z} \left( \lambda \frac{\partial T}{\partial z} \right) $ uses central differences.

These terms are discretized in the same way as was described above for the momentum equation for the upwinding term($\rho c_p u \frac{\partial T}{\partial z}$) and the second derivative term ($\( \frac{\partial}{\partial z} \left( \lambda \frac{\partial T}{\partial z} \right) \)$). An additional term that needs
to be discretized in this equation is the quantity that involes the species diffusive mass fluxes, $ j_k $ and the gradient of enthalpy.


## Species Equation

The species equation is given by:

$$
\rho u \frac{\partial Y_k}{\partial z} = -\frac{\partial j_k}{\partial z} + W_k \dot{\omega}_k
$$

**Discretization:**
- The term $ \rho u \frac{\partial Y_k}{\partial z} $ uses upwinding.
- The second derivative term $ -\frac{\partial j_k}{\partial z} $ uses a central difference formula in the same way as has been described earlier.


For the interior points in the domain, moving all terms to the right-hand-side, the discretized equation is:
$$
F = -\rho_j u_j \left( \frac{Y(x, j_{\text{loc}}) - Y(x, j_{\text{loc} - 1})}{m_{\text{dz}}[j_{\text{loc} - 1}]} \right) -
  \frac{j_{k, j+1/2} - j_{k, j-1/2}}{z_{j+1/2} - z_{j-1/2}} + \dot \omega_{k, j} W_k
$$

As a reminder, the `loc` is a stand-in variable that depends on the actual upwinding direction. If the axial velocity is positive, the value of `loc` is j. If the axial velocity is negative, the value of `loc` is `j+1`. This discretization can be seen in {cite:t}`kee2003` , equation 16.106.


These discretizations are taken from the work of Kee[2003].

