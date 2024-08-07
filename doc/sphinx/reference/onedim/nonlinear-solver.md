# Nonlinear Solver for One-dimensional Flows

# Newton Solver for \( F(x) = 0 \)

## Overview

The Newton solver is a numerical method used to find the roots of a nonlinear equation $ F(x) = 0 $. Cantera uses this method for solving the 1-dimensional
flame equations.

## Method

The function F(x) is a vector function where x is a large vector of all solution components at all grid points.

$$
\mathbf{x} = \begin{pmatrix}
u_0 \\
V_0 \\
T_0 \\
\lambda_0 \\
u_1 \\
V_1 \\
T_1 \\
\vdots
\end{pmatrix}
$$

### Initial Guess
Start with an initial guess $ x^{(0)} $ for the solution.

### Iteration
For each iteration, k, the solution is updated using the following relation:
$$
J(x^{(k)}) \left( x^{(k+1)} - x^{(k)} \right) = -F(x^{(k)}) \quad (k = 0, 1, 2, 3, \ldots)
$$

Here, $ \mathbf{J}(x^{(k)}) $ is the Jacobian matrix of $ F(x^{(k)}) $ .

Another way to looking at the equation is in the following:

$$
 \left( x^{(k+1)} - x^{(k)} \right) = -J(x^{(k)})^{-1} F(x^{(k)})  = \Delta x^{(k)}
$$

Where $ \Delta x^{(k)} $ is a vector that represents the correction that the equations want to
perform to take the solution from $ x^{(k)} $ to $ x^{(k+1)} $. During each iteration, this
correction vector changes and during convergence, should reduce in size until it becomes zero.
This correction vector points in the direction that will drive the solution towards a zero
error.

### Jacobian Matrix
The Jacobian matrix $ \mathbf{J} $ is the matrix of partial derivatives:

$$
\mathbf{J}(x^{(k)}) = \frac{\partial F(x^{(k)})}{\partial x}
$$

This is approximated numerically in the 1D solver instead of having analytical relations derived for
each governing equation.

### Damped Newton Method

This method is essentially the same as the previous method, but with a new damping parameter introduced.

$$
 \Delta x^{(k)} = -\lambda^{(m)} J(x^{(k)})^{-1} F(x^{(k)})
$$

The damping parameter, $\lambda$ is a value that is between 0 and 1. This damping parameter is selected to satisfy two conditions:
  1. Each component of x must stay within a trust region, which is the bounds that are assigned to each solution component. Things such as
      limitations on the magnitude or sign of the velocity, mass fractions, etc.
  2. The norms of succeeding undamped steps decrease in magnitude.


```{image} /_static/images/damped_newton_diagram.svg
:width: 75%
:alt: Representation of the damped Newton method. Adapted from {cite:t}`kee2003`.
:align: center

Representation of the damped Newton method. Adapted from {cite:t}`kee2003`.
```

The algorithm is that given:

$$
x^{(m+1)} = x^{(m)} + \lambda^{(m)} \Delta x^{(m)}
$$

We need to pick a value of $\lambda$ that satisfies:

$$
|\Delta x^{(m+1)}| < |\Delta x^{(m)}|
$$

Where $ \Delta x^{(m)} = J(x^{(m)})^{-1} F(x^{(m)}) $ and $ \Delta x^{(m+1)} = J(x^{(m)})^{-1} F(x^{(m+1)}) $.

During the search for the correct value of $\lambda$, the value of $\lambda$ starts at 1, it is adjusted down to a value that keeps the solution within
the trust region. The process then begins for finding $\lambda$, failures result in the damping factor being reduced by a constant factor. The current
factor in  Cantera is the $\sqrt{2}$.


## Convergence Criteria

A damped newton step is considered to be converged when the norm of the correction vector is less than 1. During the steady-state solution, the process of
finding a damped newton step and taking it is repeated until the norm of the correction vector is less than 1, if it is not, then the process continues.


## Transient Solution

Sometimes the solution of the steady-state problem can not be found from the damped newton method. In this case, a transient solution is solve for a specified
number of time steps before the damped newton method is attempted again. The transient equation for stepping forward in time (using backward euler method) is given by:

$$
\left( J(x^{(k)}) - \frac{I}{\Deltat} \right) \Delta x^{(n+1)} = -F(x^{(n)})
$$

Here the `n+1` is the solution at the next time step, `n` is the solution at the current time step. One can see that this problem is analogous to the method that was just described above with the exception that the Jacobian matrix is modified to include the time step size. In Cantera, this is exactly what is done. The Jacobian matrix is modified to include the time step size, and the damped newton method is then used to solve the transient problem because it is the same as the steady-state problem, with the exception of that modified Jacobian matrix.