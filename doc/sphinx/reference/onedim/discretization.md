# Discretization of 1D Equations

# Governing Equations

This section outlines the discretizations used for the governing equations, including the continuity, radial momentum, energy, and species equations, see(governingEquationsSection) for more details on the governing equations. The discretization methods applied to these equations include upwinding for terms involving \( u \), central differences for second derivative terms, and one-sided discretizations for boundary conditions.

## Continuity Equation

The continuity equation is given by:

$$
\[ \frac{\partial (\rho u)}{ \partial z} + 2\rho V = 0 \]
$$

**Discretization:**
The continuity equation does not need any upwinding in the standard sense, but because the equation is first order and the
boundary condition is applied at the right boundary in the 1D domain, the first derivative approximation uses a one-sided
form.


## Radial Momentum Equation

The radial momentum equation is:

\[ \rho u \frac{\partial V}{\partial z} + \rho V^2 = -\Lambda + \frac{\partial}{\partial z} \left( \mu \frac{\partial V}{\partial z} \right) \]

**Discretization:**
- The term \( \rho u \frac{\partial V}{\partial z} \) uses upwinding for \( u \).
- The second derivative term \( \frac{\partial}{\partial z} \left( \mu \frac{\partial V}{\partial z} \right) \) uses central differences.

## Energy Equation

The energy equation is described by:

\[ \rho c_p u \frac{\partial T}{\partial z} = \frac{\partial}{\partial z} \left( \lambda \frac{\partial T}{\partial z} \right) - \sum_k j_k \frac{\partial h_k}{\partial z} - \sum_k h_k W_k \dot{\omega}_k \]

**Discretization:**
- The term \( \rho c_p u \frac{\partial T}{\partial z} \) is discretized using upwinding for \( u \).
- The second derivative term \( \frac{\partial}{\partial z} \left( \lambda \frac{\partial T}{\partial z} \right) \) uses central differences.

## Species Equation

The species equation is given by:

\[ \rho u \frac{\partial Y_k}{\partial z} = -\frac{\partial j_k}{\partial z} + W_k \dot{\omega}_k \]

**Discretization:**
- The term \( \rho u \frac{\partial Y_k}{\partial z} \) uses upwinding for \( u \).
- The second derivative term \( -\frac{\partial j_k}{\partial z} \) uses central differences.

---

These discretizations are taken from the work of Kee[2017].

