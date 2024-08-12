# Nonlinear Solver for One-dimensional Flows

## Overview

Cantera uses a hybrid time stepping / steady-state algorithm to solve the discretized
1-dimensional flame equations. For both the time stepping and steady-state problems,
a damped Newton's method solver is used. The general principles of the solver used in
Cantera are described in {cite:t}`kee2003`(Chapter 15).

## Problem Definition

The solution to the 1-dimensional set of governing equations is expressed in the form
of a root finding equation for use with a Newton solver. The equation to be solved
takes the form of $F(x) = 0$. The function $F(x)$ is a vector-value function of the
solution vector, $x$, which is a vector of all solution components at all grid points.

$$
x =
\begin{pmatrix}
u_0 \\
V_0 \\
T_0 \\
\Lambda_0 \\
Y_{1,0} \\
\vdots \\
u_1 \\
V_1 \\
\vdots \\
Y_{m,N}
\end{pmatrix}
$$

The vector-value function $F(x)$ is a vector of the residuals for each of the governing
equations at all of the grid points.

$$
F(x) =
\begin{pmatrix}
F_{u,0}(x) \\
F_{V,0}(x) \\
F_{T,0}(x) \\
F_{\Lambda,0}(x) \\
F_{Y_1,0}(x) \\
\vdots \\
F_{u,1}(x) \\
F_{V,1}(x) \\
\vdots \\
F_{Y_m,N}(x)
\end{pmatrix}
$$

Residuals in this context, at interior grid
points, are the difference between the left-hand side and right-hand side of the
governing equations. If the perfect solution was obtained, then the difference
between the left-hand side and right-hand side of the governing equations would be
zero, and this is what the solver is trying to achieve by examining the residuals
during each attempt at solving the system of equations.

One of the key components of the solver is the Jacobian matrix, which is the matrix of
partial derivatives of the residuals with respect to the solution vector. The Jacobian
matrix is used to determine the direction of the correction vector that will drive the
solution towards zero error.

$$
J(x) =
\frac{\partial F(x)}{\partial x}
$$

For the vector-value residual vector $F(x)$, the Jacobian matrix looks like,

$$
\begin{pmatrix}
\frac{\partial F_{U,0}(x)}{\partial U_0} &
\frac{\partial F_{U,0}(x)}{\partial V_0} &
\frac{\partial F_{U,0}(x)}{\partial T_0} &
\frac{\partial F_{U,0}(x)}{\partial Y_{m,0}} &
\cdots &
\frac{\partial F_{U,0}(x)}{\partial U_N} &
\frac{\partial F_{U,0}(x)}{\partial V_N} &
\frac{\partial F_{U,0}(x)}{\partial T_N} &
\frac{\partial F_{U,0}(x)}{\partial Y_{m,N}}

\\

\frac{\partial F_{V,0}(x)}{\partial U_0} &
\frac{\partial F_{V,0}(x)}{\partial V_0} &
\frac{\partial F_{V,0}(x)}{\partial T_0} &
\frac{\partial F_{V,0}(x)}{\partial Y_{m,0}} &
\cdots &
\frac{\partial F_{V,0}(x)}{\partial U_N} &
\frac{\partial F_{V,0}(x)}{\partial V_N} &
\frac{\partial F_{U,0}(x)}{\partial T_N} &
\frac{\partial F_{V,0}(x)}{\partial Y_{m,N}}

\\

\vdots & \vdots & \vdots & \vdots & & \vdots & \vdots & \vdots & \vdots

\\

\frac{\partial F_{U,N}(x)}{\partial U_0} &
\frac{\partial F_{U,N}(x)}{\partial V_0} &
\frac{\partial F_{U,N}(x)}{\partial T_0} &
\frac{\partial F_{U,N}(x)}{\partial Y_{m,0}} &
\cdots &
\frac{\partial F_{U,N}(x)}{\partial U_N} &
\frac{\partial F_{U,N}(x)}{\partial V_N} &
\frac{\partial F_{U,N}(x)}{\partial T_N} &
\frac{\partial F_{U,N}(x)}{\partial Y_{m,N}}

\\

\frac{\partial F_{V,N}(x)}{\partial U_0} &
\frac{\partial F_{V,N}(x)}{\partial V_0} &
\frac{\partial F_{V,N}(x)}{\partial T_0} &
\frac{\partial F_{V,N}(x)}{\partial Y_{m,0}} &
\cdots &
\frac{\partial F_{V,N}(x)}{\partial U_N} &
\frac{\partial F_{V,N}(x)}{\partial V_N} &
\frac{\partial F_{V,N}(x)}{\partial T_N} &
\frac{\partial F_{V,N}(x)}{\partial Y_{m,N}}

\end{pmatrix}
$$

This is approximated numerically in the 1D solver instead of having analytical
relations derived for each governing equation.

### Damped Newton Method

The damped Newton method starts with an initial guess for the solution, $x^{(0)}$, and
performs a series of iterations until the solution converges with the help of a
damping parameter.

For each iteration, k, the solution is updated using the following relation:

$$
J(x^{(k)}) \left( x^{(k+1)} - x^{(k)} \right) = -\lambda^{(k)} F(x^{(k)})
  \quad (k = 0, 1, 2, 3, \ldots)
$$

Here, $ J(x^{(k)}) $ is the Jacobian matrix of $ F(x^{(k)}) $.

Another way to looking at the equation is:

$$
 \left( x^{(k+1)} - x^{(k)} \right) = -\lambda^{(k)} J(x^{(k)})^{-1} F(x^{(k)})
   = \Delta x^{(k)}
$$

Where $ \Delta x^{(k)} $ is a vector that represents a correction to the current
solution that will take the solution from $ x^{(k)} $ to $ x^{(k+1)} $. During each
iteration this correction vector changes and during the iteration process it
should reduce in size until it becomes zero. This correction vector points in the
direction that will drive the solution towards a zero error.

The damping parameter, $\lambda$ is a value that is between 0 and 1. This damping
parameter is selected to satisfy two conditions:
  1. Each component of $x$ must stay within a trust region, which is the bounds that
     are assigned to each solution component. These are bounds such as limitations on
     the magnitude or sign of the velocity, mass fractions, etc.
  2. The norms of succeeding undamped steps decrease in magnitude.

The following image visually illustrates the damped Newton method. In it, the undamped
Newton step \( $ \Delta x^{(k)} $ \) is shown. The second vector is the undamped step
that would occur if the full initial step were to be taken. The shorter vector
representing the damped correction \( $ \lambda^{(k)} \Delta x^{(k)} $ \) is a trial
step. The method takes this trial step to get to a new solution at \( $ x^{(k+1)} $ \).
A new step vector is then computed using the trial solution, which gives
$ \Delta x^{(k+1)} $. The length of this vector is compared to the length of the
original undamped step \( $ \Delta x^{(k)} $ \). If the length of the new step is less
than the length of the original step, then the trial step is accepted, and the damping
value is accepted. The step is then taken, and the process is repeated until the
solution converges.

```{figure} /_static/images/damped_newton_diagram.svg
:width: 75%
:alt: Representation of the damped Newton method. Adapted from {cite:t}`kee2003`.
:align: center

Representation of the damped Newton method. Adapted from {cite:t}`kee2003`.
```

For a more mathematical representation of the damped Newton method, we consider:

$$
x^{(m+1)} = x^{(m)} + \lambda^{(m)} \Delta x^{(m)}
$$

A value of $\lambda$ needs to be picked that satisfies:

$$
|\Delta x^{(m+1)}| < |\Delta x^{(m)}|
$$

Where:

$$
\Delta x^{(m)} = J(x^{(m)})^{-1} F(x^{(m)})
$$

and,

$$
\Delta x^{(m+1)} = J(x^{(m)})^{-1} F(x^{(m+1)})
$$

During the search for the correct value of $\lambda$, the value of $\lambda$ starts at 1, it is adjusted down to a value that keeps the solution within
the trust region. The process then begins for finding $\lambda$, failures result in the damping factor being reduced by a constant factor. The current
factor in  Cantera is the $\sqrt{2}$.

During the damped Newton method, the Jacobian is kept at the $x^{(m)}$ value. This sometimes can cause issues with convergence if the Jacobian becomes out of date (due to sensitivity to the solution). In Cantera, the Jacobian is updated if too many attempted steps fail to take any damped step. The balance of how long to wait before updating the Jacobian is a trade-off between the cost of updating the Jacobian and the cost of failing to converge.


### Convergence Criteria

As was discussed earlier, the Newton method is an iterative method, and it's important to assess when the method has reached a point where the iterations
can be stopped. This point is called convergence. Cantera's implementation uses a **weighted norm** of the step vector to determine convergence, rather than a simple absolute norm. A damped Newton step is considered to be converged when the **weighted norm** of the correction vector is less than 1. During the solution, the process of
finding and taking a damped Newton step is repeated until the **weighted norm** of the correction vector is less than 1, if it is not, then the process continues.

In a multivariate system, different variables may have vastly different magnitudes and units. A simple absolute norm could either be dominated by large components or fail to account for smaller components effectively. By normalizing the step vector components using \( w_n \), the weighted norm ensures that the convergence criterion is meaningful across all solution components, regardless of their individual scales.

This approach provides a more robust and scale-invariant method for assessing convergence, making it especially useful in systems with diverse variables.

#### Definition of the Weighted Norm

The weighted norm of the step vector \( \mathbf{s} \) is calculated as:

$$
\text{Weighted Norm} = \sum_{n,j} \left(\frac{s_{n,j}}{w_n}\right)^2
$$

where:
- \( s_{n,j} \) is the Newton step vector component for the \( n \)-th solution variable at the \( j \)-th grid point.
- \( w_n \) is the error weight for the \( n \)-th solution component, given by:

$$
w_n = \epsilon_{r,n} \cdot \frac{\sum_j |x_{n,j}|}{J} + \epsilon_{a,n}
$$

Here:
- \( \epsilon_{r,n} \) is the relative error tolerance for the \( n \)-th solution component.
- \( \frac{\sum_j |x_{n,j}|}{J} \) is the average magnitude of the \( n \)-th solution component over all grid points, and $J$ is the total number of grid points.
- \( \epsilon_{a,n} \) is the absolute error tolerance for the \( n \)-th solution component.

#### Interpretation of the Weighted Norm

The weighted norm is a relative measure that helps bring all components of the step vector into a comparable range, taking into account the scales of the different solution components. Here's how to interpret it:

- **Relative Error Term** \( (\epsilon_{r,n}) \): Scales the step size relative to the average magnitude of the corresponding solution component. This means that larger components can tolerate larger steps.
- **Absolute Error Term** \( (\epsilon_{a,n}) \): Ensures that even very small solution components are considered in the convergence check by providing a minimum threshold.

#### Convergence Criterion

The Newton iteration is considered converged when the weighted norm is less than 1:

\[
\sum_{n,j} \left(\frac{s_{n,j}}{w_n}\right)^2 < 1
\]

This criterion indicates that each component of the step vector \( s_{n,j} \) is sufficiently small relative to the expected precision (as defined by the weights \( w_n \)).

## Transient Solution

There will be times when the solution of the steady-state problem can not be found using the damped Newton method. In this case, a transient solution is solved and a specified number of time steps is taken before the steady-state damped Newton method is attempted again. The transient equation for stepping forward in time (using backward euler method) is given by:

$$
\left( J(x_{(n)}) - \frac{I}{\Delta t} \right) \Delta x_{(n+1)} = -F(x_{(n)})
$$

Here the `n+1` is the solution at the next time step, `n` is the solution at the current time step. This problem is analogous to the method that was just described above with the exception that the Jacobian matrix is modified to include the time step size. In Cantera, this is exactly what is done. The Jacobian matrix is modified to include the time step size, and the damped Newton method is then used to solve the transient problem because it has the same solution procedure as the steady-state problem, with the exception of that modified Jacobian matrix.

### Differential-Algebraic Equation (DAE) Form

A generic DAE system can be expressed as:

$$
F\left(t, x, \frac{dx}{dt}\right) = 0
$$

Where:
- $x$ is the solution vector.
- $\frac{dx}{dt}$ represents the time derivatives of the differential components of $x$.
- $F(t, x, \frac{dx}{dt})$ is the residual vector, which combines both differential equations and algebraic constraints.

### Backward Euler Time Discretization

Using the backward Euler method to discretize the time derivative:

$$
\frac{dx}{dt} \approx \frac{x_{n+1} - x_n}{\Delta t}
$$

The residual equation for the time step $n+1$ becomes:

$$
F\left(t_{n+1}, x_{n+1}, \frac{x_{n+1} - x_n}{\Delta t}\right) = 0
$$

### Newton's Method for Solving the Nonlinear System

To solve the nonlinear system using Newton's method:

1. **Newton Iteration**:

   $$
   J\left(x_{n+1}^{(k)}\right) \Delta x_{n+1}^{(k)} = -F\left(t_{n+1}, x_{n+1}^{(k)}, \frac{x_{n+1}^{(k)} - x_n}{\Delta t}\right)
   $$

   Where:
   - $x_{n+1}^{(k)}$ is the solution at the $k$-th Newton iteration for the time step $n+1$.
   - $J$ is the Jacobian matrix, which is the derivative of the residual function with respect to $x$.

2. **Jacobian Structure**:
   The Jacobian matrix is composed of:

   $$
   J = \frac{1}{\Delta t} I + \frac{\partial F}{\partial x_{n+1}^{(k)}}
   $$

   - The term $\frac{1}{\Delta t} I$ arises from the derivative of the backward Euler time discretization.
   - The second term is the Jacobian of the steady-state residual with respect to $x_{n+1}^{(k)}$.

3. **Update Solution**:

   $$
   x_{n+1}^{(k+1)} = x_{n+1}^{(k)} + \Delta x_{n+1}^{(k)}
   $$

4. **Convergence**:
   - Iterate until the Newton correction $\Delta x_{n+1}^{(k)}$ is sufficiently small.
   - Once converged, $x_{n+1}^{(k+1)}$ is accepted as $x_{n+1}$, and the process advances to the next time step.


