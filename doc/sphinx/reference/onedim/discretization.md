# Discretization of 1D Equations

# Governing Equations

The governing equations outlined in the [governing equations](./governing-equations) must be discretized in order to obtain solutions for 1D flame configurations. Various discretization methods are used for different terms in the equations, and their details are outlined here.

For this discussion, consider a discretized 1D domain. Canter uses a non-uniform node spacing for the 1D domain.

```{image} /_static/images/1d_domain_diagram_1.svg
:width: 75%
:alt: Diagram of the 1D domain
:align: center
```

The 1D solver sweeps over all points in a domain and computes a residual for each grid point and each equation. The solver attempts to find the combination of the variables that minimizes the residual vector. This residual vector has entries in a form similar to what is shown below.

$$
F(x) = F(u_0, V_0, T_0, \lambda_0, u_1, V_1, T_1, \lambda_1, ..., u_{N-1}, V_{N-1}, T_{N-1}, \lambda_{N-1}, u_{N}, V_{N}, T_{N}, \lambda_{N})
$$

For each governing equation, there is a default boundary condition applied at the left and right boundaries. During the sweep, the final evaluations are the  done for the domain boundary points, which allows for the default residual to be augmented or overwritten by a boundary condition that is
different from the default boundary condition.


## Continuity Equation

The continuity equation is:

$$
 \frac{\partial (\rho u)}{ \partial z} + 2\rho V = 0
$$


### Discretization

The continuity equation does not need any upwinding in the standard sense, but because the equation is first order and the
boundary condition is applied at the right boundary in the 1D domain, the first derivative approximation uses a one-sided
form which uses a stencil that pulls from points to right of a grid point. For consistency, the radial momentum term, \($2\rho V$\), also uses the same points and averages the value between two nodes.


The discretized equation in residual form (all terms moved to one side) at the interior points in the domain is given below.

$$
F = \frac{\rho_j u_j - \rho_{j+1} u_{j+1}}{z_{j+1} - z_j} + 2 \left( \frac{\rho_j V_j + \rho_{j+1} V_{j+1}}{2} \right)
$$


### Boundary Conditions

At the right boundary, the default boundary condition is a zero velocity condition, which is representative of a stagnation
surface.

At the right boundary (`j=N`):

$$
F = \rho_j u_j - \dot m
$$

Expressing the boundary residual in this form will drive the Newton root finding algorithm to find the value of $ \rho_j u_j $ that
minimizes the difference between it and the imposed (from a user-defined boundary mass flow rate specification, for example) boundary condition.

There is no imposed boundary condition at the left boundary because only one boundary condition can be enforced for a first-order differential equation. As such, the residual equation at the left boundary is the same as the interior points.

At the left boundary (`j=0`):

$$
F = \frac{\rho_j u_j - \rho_{j+1} u_{j+1}}{z_{j+1} - z_j} + 2 \left( \frac{\rho_j V_j + \rho_{j+1} V_{j+1}}{2} \right)
$$


## Radial Momentum Equation

The radial momentum equation is:

$$
\rho u \frac{\partial V}{\partial z} + \rho V^2 = -\Lambda + \frac{\partial}{\partial z} \left( \mu \frac{\partial V}{\partial z} \right)
$$

### Discretization

- The term $ \rho u \frac{\partial V}{\partial z} $ uses upwinding.
- The second derivative term $ \frac{\partial}{\partial z} \left( \mu \frac{\partial V}{\partial z} \right) $ uses a central difference approximation.


The discretized equation in residual form (all terms moved to one side) at the interior points in the domain is given below.

$$
-\rho_j u_j \left( \frac{V_{j_{\text{loc}}} - V_{j_{\text{loc} - 1}}}{z_{jloc} - z_{jloc-1}} \right) - \rho_j V_j^2 + \Lambda_j +
  \frac{\mu_{j+1/2} \frac{V_{j+1} - V_j}{z_{j+1} - z_j} - \mu_{j-1/2} \frac{V_j - V_{j-1}}{z_j - z_{j-1}}}{\frac{z_{j+1} - z_{j-1}}{2}}
$$

#### Upwinding

The upwinding formula for the radial velocity derivative term \( $ \rho u \frac{\partial V}{\partial z} $ \) is:

$$
\left( \rho u \frac{\partial A}{\partial z} \right) \bigg|_{j} \approx \rho_j u_j \frac{A_{j_{\text{loc}}} -
  A_{j_{\text{loc} - 1}}}{z_{jloc} - z_{jloc-1}}
$$

Where the value of `loc` is determined by the sign of the axial velocity \(u\). If the axial velocity is positive, the value of `loc` is j. If the axial velocity
is negative, the value of `loc` is `j+1`. A positive velocity means that the flow is moving left-to-right.


#### Second Derivative Term

For the second derivative term (the shear term in the momentum equation), a three-point central difference formula is used. The term being discretized is:

$$
\frac{d}{dz}\left(\mu \frac{dV}{dz}\right)
$$

For simplicity, let $ A = \mu \frac{dV}{dz} $ for simplicity. This will be called the inner term. In this situation, the inner term is evaluated using a central difference formula, but instead of using the `j+1` and `j-1` points, the derivative is estimated using `j+1/2` and `j-1/2` (halfway between the grid points around point j).

The values of $A$ at point `j-1/2` and `j+1/2` are estimating using a central difference formula:

$$
A \big|_{j-1/2} = \mu_{j-1/2} \frac{V_j - V_{j-1}}{z_j - z_{j-1}}
$$


$$
A \big|_{j+1/2} = \mu_{j+1/2} \frac{V_{j+1} - V_j}{z_{j+1} - z_j}
$$

$ \mu_{j+1/2} $ is the viscosity, estimated at the midpoint between grid points.

The outer discretization uses a central difference between the `j+1/2` and `j-1/2`
locations.

$$
\frac{dA}{dz} \bigg|_{j} \approx \frac{A \big|_{j+1/2} - A \big|_{j-1/2}}{z_{j+1/2} - z_{j-1/2}}
$$

Where the values of $ z $ are:

$$
z_{j+1/2} = z_{j} + \frac{1}{2}(z_{j+1} - z_j) = \frac{1}{2}(z_{j} + z_{j+1})
$$

and,

$$
z_{j-1/2} = z_j - \frac{1}{2}(z_j - z_{j-1}) = \frac{1}{2}(z_{j} + z_{j-1})
$$

The difference between these two values is:

$$
z_{j+1/2} - z_{j-1/2} = \frac{z_{j+1} - z_{j-1}}{2}
$$

Substituting these values into the central difference formula gives:

$$
\frac{d}{dz}\left(\mu \frac{dV}{dz}\right) \approx \frac{\mu_{j+1/2} \frac{V_{j+1} - V_j}{z_{j+1} - z_j} -
  \mu_{j-1/2} \frac{V_j - V_{j-1}}{z_j - z_{j-1}}}{\frac{z_{j+1} - z_{j-1}}{2}}
$$

This formula is used to discretize the second derivative term in the radial momentum equation. It takes this form due to the
construction of a conservative scheme, meaning that point `j` sees the same value of the flux \( $\mu \frac{dV}{dz}$ \) on its
right side as the point `j+1` sees on its left side. This is a conservative scheme because the fluxes are balanced at each point.

### Boundary Conditions

At the right boundary, the default boundary condition is a zero radial velocity.

At the right boundary (`j=N`):

$$
F = V_j
$$

Expressing the boundary residual in this form will drive the Newton root finding algorithm to find the value of $ V_j $ that
minimizes the residual, which in this case, the value of $ V_j $ that minimizes the relation above is $V=0$ .

The same boundary condition is used at the left boundary.

At the left boundary (`j=0`):

$$
F = V_j
$$



## Energy Equation

The steady-state energy equation is described by:

$$
\rho c_p u \frac{\partial T}{\partial z} = \frac{\partial}{\partial z} \left( \lambda \frac{\partial T}{\partial z} \right) - \sum_k j_k \frac{\partial h_k}{\partial z} - \sum_k h_k W_k \dot{\omega}_k
$$


### Discretization

- The term $ \rho c_p u \frac{\partial T}{\partial z} $ uses upwinding.
- The second derivative term $ \frac{\partial}{\partial z} \left( \lambda \frac{\partial T}{\partial z} \right) $ uses a central difference approximation.

These terms are discretized in the same way as was described above for the momentum equation for the upwinded term \( $\rho c_p u \frac{\partial T}{\partial z}$ \) and the second derivative term \( $ \frac{\partial}{\partial z} \left( \lambda \frac{\partial T}{\partial z} \right) $ \).
An additional term that needs to be discretized in this equation is the quantity that involves the species diffusive mass fluxes, $ j_k $ and the gradient of enthalpy, \( $ \sum_k j_k \frac{\partial h_k} $ \).


The discretized equation in residual form (all terms moved to one side) at the interior points in the domain is given below.

$$
-\rho_j c_p u_j \left( \frac{T_{j_{\text{loc}}} - T_{j_{\text{loc} - 1}}}{z_{jloc} - z_{jloc-1}} \right) +
  \frac{\lambda_{j+1/2} \frac{T_{j+1} - T_j}{z_{j+1} - z_j} - \lambda_{j-1/2} \frac{T_j - T_{j-1}}{z_j - z_{j-1}}}{\frac{z_{j+1} - z_{j-1}}{2}} -
  \sum_k j_{k, j} \left( \frac{h_{k, j_{\text{loc}}} - h_{k, j_{\text{loc} - 1}}}{z_{jloc} - z_{jloc-1}} \right)
$$

The upwinding on the enthalpy gradient term uses upwinding, but it does not need to use it because it is not multiplied by the axial velocity \( u \). This discretization may change in the future.

### Boundary Conditions

At the right boundary, the default boundary condition is a zero temperature.

At the right boundary (`j=N`):

$$
F = T_j
$$

Expressing the boundary residual in this form will drive the Newton root finding algorithm to find the value of $ T_j $ that
minimizes the residual, which in this case, the value of $ T_j $ that minimizes the relation above is $T=0$ .

The same boundary condition is used at the left boundary.

At the left boundary (`j=0`):

$$
F = T_j
$$


## Species Equation

The species equation is given by:

$$
\rho u \frac{\partial Y_k}{\partial z} = -\frac{\partial j_k}{\partial z} + W_k \dot{\omega}_k
$$

**Discretization:**
- The term $ \rho u \frac{\partial Y_k}{\partial z} $ uses upwinding.
- The diffusive mass flux term $ \frac{\partial j_k}{\partial z} $ uses a conservative central difference formula.


For the interior points in the domain, moving all terms to the right-hand-side, the discretized equation is:

$$
F = -\rho_j u_j \left( \frac{Y_{k, j_{\text{loc}}} - Y_{j_{\text{loc} - 1}}}{z_{jloc} - z_{jloc-1}} \right) -
  \frac{j_{k, j+1/2} - j_{k, j-1/2}}{z_{j+1/2} - z_{j-1/2}} + \dot \omega_{k, j} W_k
$$

This discretization can be seen in {cite:t}`kee2003` , equation 16.106.

## Boundary Conditions

At the right boundary, the default boundary condition is a zero species mass-flux. This is different from a zero species
mass fraction condition, because species may diffuse towards the boundaries during the solution process.

At the right boundary (`j=N`):

$$
F = \rho_j u_j Y_{k, j} - j_{k, j-1/2}
$$

The same boundary condition is used at the left boundary.

At the left boundary (`j=0`):

$$
F = \rho_j u_j Y_{k, j} - j_{k, j+1/2}
$$