# Nonlinear Solver for One-dimensional Flows

## Overview

Cantera uses a hybrid time stepping / steady-state algorithm to solve the discretized
1-dimensional flame equations. For both the time stepping and steady-state problems,
a damped Newton's method solver is used.

## Problem Definition

The solution to the 1-dimensional set of governing equations is expressed in the form
of a root finding equation for use with a Newton solver. The equation to be solved
takes the form of $F(x) = 0$. The function $F(x)$ is a vector-value function of the
solution vector, $x$, which is a vector of all solution components at all grid points.

$$
\mathbf{x} =
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
\mathbf{F}(\mathbf{x}) =
\begin{pmatrix}
F_{u,0}(\mathbf{x}) \\
F_{V,0}(\mathbf{x}) \\
F_{T,0}(\mathbf{x}) \\
F_{\Lambda,0}(\mathbf{x}) \\
F_{Y_1,0}(\mathbf{x}) \\
\vdots \\
F_{u,1}(\mathbf{x}) \\
F_{V,1}(\mathbf{x}) \\
\vdots \\
F_{Y_m,N}(\mathbf{x})
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
\mathbf{J}(\mathbf{x}^{(k)}) =
\frac{\partial \mathbf{F}(\mathbf{x}^{(k)})}{\partial \mathbf{x}} =
\begin{pmatrix}
\frac{\partial F_{U,0}(\mathbf{x}^{(k)})}{\partial U_0} &
\frac{\partial F_{U,0}(\mathbf{x}^{(k)})}{\partial V_0} &
\frac{\partial F_{U,0}(\mathbf{x}^{(k)})}{\partial T_0} &
\frac{\partial F_{U,0}(\mathbf{x}^{(k)})}{\partial Y_{m,0}} &
\cdots &
\frac{\partial F_{U,0}(\mathbf{x}^{(k)})}{\partial U_N} &
\frac{\partial F_{U,0}(\mathbf{x}^{(k)})}{\partial V_N} &
\frac{\partial F_{U,0}(\mathbf{x}^{(k)})}{\partial T_N} &
\frac{\partial F_{U,0}(\mathbf{x}^{(k)})}{\partial Y_{m,N}}

\\

\frac{\partial F_{V,0}(\mathbf{x}^{(k)})}{\partial U_0} &
\frac{\partial F_{V,0}(\mathbf{x}^{(k)})}{\partial V_0} &
\frac{\partial F_{V,0}(\mathbf{x}^{(k)})}{\partial T_0} &
\frac{\partial F_{V,0}(\mathbf{x}^{(k)})}{\partial Y_{m,0}} &
\cdots &
\frac{\partial F_{V,0}(\mathbf{x}^{(k)})}{\partial U_N} &
\frac{\partial F_{V,0}(\mathbf{x}^{(k)})}{\partial V_N} &
\frac{\partial F_{U,0}(\mathbf{x}^{(k)})}{\partial T_N} &
\frac{\partial F_{V,0}(\mathbf{x}^{(k)})}{\partial Y_{m,N}}

\\

\vdots & \vdots & \vdots & \vdots & & \vdots & \vdots & \vdots & \vdots

\\

\frac{\partial F_{U,N}(\mathbf{x}^{(k)})}{\partial U_0} &
\frac{\partial F_{U,N}(\mathbf{x}^{(k)})}{\partial V_0} &
\frac{\partial F_{U,N}(\mathbf{x}^{(k)})}{\partial T_0} &
\frac{\partial F_{U,N}(\mathbf{x}^{(k)})}{\partial Y_{m,0}} &
\cdots &
\frac{\partial F_{U,N}(\mathbf{x}^{(k)})}{\partial U_N} &
\frac{\partial F_{U,N}(\mathbf{x}^{(k)})}{\partial V_N} &
\frac{\partial F_{U,N}(\mathbf{x}^{(k)})}{\partial T_N} &
\frac{\partial F_{U,N}(\mathbf{x}^{(k)})}{\partial Y_{m,N}}

\\

\frac{\partial F_{V,N}(\mathbf{x}^{(k)})}{\partial U_0} &
\frac{\partial F_{V,N}(\mathbf{x}^{(k)})}{\partial V_0} &
\frac{\partial F_{V,N}(\mathbf{x}^{(k)})}{\partial T_0} &
\frac{\partial F_{V,N}(\mathbf{x}^{(k)})}{\partial Y_{m,0}} &
\cdots &
\frac{\partial F_{V,N}(\mathbf{x}^{(k)})}{\partial U_N} &
\frac{\partial F_{V,N}(\mathbf{x}^{(k)})}{\partial V_N} &
\frac{\partial F_{V,N}(\mathbf{x}^{(k)})}{\partial T_N} &
\frac{\partial F_{V,N}(\mathbf{x}^{(k)})}{\partial Y_{m,N}}

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

Here, $ \mathbf{J}(x^{(k)}) $ is the Jacobian matrix of $ F(x^{(k)}) $.

Another way to looking at the equation is:

$$
 \left( x^{(k+1)} - x^{(k)} \right) = -\lambda^{(k)} J(x^{(k)})^{-1} F(x^{(k)})
   = \Delta x^{(k)}
$$

Where $ \Delta x^{(k)} $ is a vector that represents a correction to the current
solution that will take the solution from $ x^{(k)} $ to $ x^{(k+1)} $. During each
iteration this correction vector changes and during during the iteration process it
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

A damped newton step is considered to be converged when the norm of the correction vector is less than 1. During the steady-state solution, the process of
finding and taking a damped newton step is repeated until the norm of the correction vector is less than 1, if it is not, then the process continues.


## Transient Solution

There will be times when the solution of the steady-state problem can not be found using the damped newton method. In this case, a transient solution is solved and a specified number of time steps is taken before the steady-state damped newton method is attempted again. The transient equation for stepping forward in time (using backward euler method) is given by:

$$
\left( J(x^{(n)}) - \frac{I}{\Delta t} \right) \Delta x^{(n+1)} = -F(x^{(n)})
$$

Here the `n+1` is the solution at the next time step, `n` is the solution at the current time step. This problem is analogous to the method that was just described above with the exception that the Jacobian matrix is modified to include the time step size. In Cantera, this is exactly what is done. The Jacobian matrix is modified to include the time step size, and the damped newton method is then used to solve the transient problem because it has the same solution proceedure as the steady-state problem, with the exception of that modified Jacobian matrix.