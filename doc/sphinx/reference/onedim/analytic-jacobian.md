# Analytic Jacobian for 1D Flames

## Overview

The [nonlinear solver](nonlinear-solver) for 1D flames needs the Jacobian $\mathbf{J}$
of the residual function to compute each Newton step. The Jacobian is reused across
Newton steps and only re-evaluated when the steps stop making progress (see the
[nonlinear solver](nonlinear-solver) page). Even so, each re-evaluation is one of the
more expensive operations in the solution process, so reducing its cost is worthwhile.

By default, the Jacobian is built by finite differences: each of the $N_v$ solution
components is perturbed in turn and the resulting change in the residual gives one
column of $\mathbf{J}$. Because the discretization is block tridiagonal --- the residual
at grid point $j$ depends only on the solution at points $j-1$, $j$, and $j+1$ ---
perturbing a variable at point $p$ changes the residual only at points $p-1$, $p$, and
$p+1$. Each column therefore costs a *local* residual evaluation over three grid points,
not a sweep of the whole domain.

The analytic Jacobian mode (`domain.jacobian_mode = "analytic"`) replaces the
finite-difference perturbation for the species mass-fraction columns at interior grid
points with formulas derived from the kinetics composition derivatives. For a domain
with $K$ species and $N$ grid points, the solution vector has length

$$
N_v = (K + c) N,
$$

where $c$ is the number of non-species components per grid point (the axial velocity
$u$, scaled radial velocity $V$, temperature $T$, pressure eigenvalue $\Lambda$, and,
for two-point-controlled flames, the axial mass flux). Of the $N_v$ columns, the
analytic mode handles the $K\,(N-2)$ species columns at the $N-2$ interior points
directly. The remaining

$$
N_v - K (N-2) = c N + 2K
$$

columns --- every non-species column, plus all columns at the two boundary points ---
are still evaluated by finite differences. The relative savings therefore grow with the
species fraction $K/(K+c)$ of the system, i.e. with the size of the mechanism.

The derivation below covers the interior column points $p = 1, \ldots, N-2$ (the two
boundary points always use finite differences). All transport coefficients are treated
as frozen at the base state, matching the frozen-transport approximation of the
finite-difference Jacobian, so the two modes make equivalent approximations.

## Notation

The notation follows the [discretization](discretization) page: grid points are denoted
by a subscript, species mass fractions by $Y_{k,j}$ (species $k$ at grid point $j$) and
mole fractions by $X_{k,j}$, with $W_k$ the molecular weight of species $k$ and
$\overline{W}_j$ the mean molecular weight at point $j$. A Jacobian entry involves two
grid points: the row (residual) point $j$ and the column (perturbed-variable) point $p$.
We write the column we are computing as $\partial / \partial Y_{m,p}$ --- the derivative
with respect to the mass fraction of species $m$ at point $p$. As in the
finite-difference solver, this is an *unnormalized* perturbation: a single $Y_{m,p}$ is
varied while the other mass fractions at point $p$ are held fixed (the mixture is not
renormalized to sum to one). Summation indices over species are written as $n$, and
$\delta_{km}$ is the Kronecker delta ($1$ if $k = m$, $0$ otherwise).

## Species residual

Using the mixture-averaged diffusive flux and an upwind convection scheme, the
steady-state species residual that Cantera assembles at interior point $j$ is

$$
F_{Y_k,j} = \frac{1}{\rho_j}\left[
  W_k \dot\omega_{k,j}
  - \rho_j u_j \left.\frac{\partial Y_k}{\partial z}\right|_j
  - \zeta_j \bigl(j_{k,j+1/2} - j_{k,j-1/2}\bigr)
\right],
$$

where $\zeta_j \equiv 2 / (z_{j+1} - z_{j-1})$ and $j_{k,j+1/2}$ is the diffusive mass
flux of species $k$ at the midpoint between points $j$ and $j+1$, $\dot\omega_{k,j}$ is
the net molar production rate, and the convection derivative $\partial Y_k/\partial
z|_j$ is upwinded. (The transient term that the time-stepping solver adds to the
diagonal is handled separately and is not part of the analytic species column.) Compared
with the [discretization](discretization) page, this is the same residual divided
through by $\rho_j$, which is the form actually assembled in {ct}`Flow1D::evalSpecies`
and the reason the $1/\rho_j$ factors appear in the Jacobian below.

With frozen transport, the midpoint flux is

$$
j_{k,j+1/2} =
  \frac{\tilde{D}_{k,j+1/2}}{\Delta z_j}\bigl(X_{k,j} - X_{k,j+1}\bigr)
  + Y_{k,j} S_{j+1/2},
$$

where

$$
S_{j+1/2} = -\sum_n \frac{\tilde{D}_{n,\,j+1/2}}{\Delta z_j}
                    \bigl(X_{n,j} - X_{n,j+1}\bigr),
$$

and $\Delta z_j = z_{j+1} - z_j$. The notation $\tilde{D}_{k,\,j+1/2}$ is a reminder
that this is not a bare diffusion coefficient but the composite *flux prefactor* that
multiplies the mole-fraction gradient. For the default molar-gradient, mixture-averaged
case it is

$$
\tilde{D}_{k,\,j+1/2} = \frac{\rho\, W_k}{\overline{W}}\, D_{k,\mathrm{mix}},
$$

evaluated at the midpoint state, where $D_{k,\mathrm{mix}}$ is the mixture-averaged
diffusion coefficient of species $k$; it therefore carries units of a mass-flux density
coefficient rather than $\mathrm{m^2/s}$. This is the quantity stored in
{ct}`Flow1D::m_diff`, and (being a transport property) it is held frozen here. The term
$Y_{k,j}\,S_{j+1/2}$ is the correction flux that enforces $\sum_k j_{k,\,j+1/2} = 0$.

The mass-gradient form replaces each $X_{n,j}$ by $Y_{n,j}$, sets
$\tilde{D}_{k,\,j+1/2} = \rho\, D_{k,\mathrm{mix}}^{\,\mathrm{mass}}$, and drops the
$\overline{W}/W$ factors that appear in the derivatives below.

## Jacobian column for $Y_{m,p}$

The residual $F_{Y_k,j}$ depends on $Y_{m,p}$ through the two midpoint fluxes adjacent
to point $p$ --- $j_{k,\,p-1/2}$ (which has point $p$ as its right endpoint) and
$j_{k,\,p+1/2}$ (which has point $p$ as its left endpoint) --- and, when $j = p$,
additionally through the reaction term, the upwinded convection term, and the $1/\rho_j$
prefactor. The column therefore has nonzero entries only in rows $j = p-1,\, p,\, p+1$:

$$
\frac{\partial F_{Y_k,j}}{\partial Y_{m,p}} =
\begin{cases}
\displaystyle
-\frac{\zeta_j}{\rho_j} \frac{\partial j_{k,\,p-1/2}}{\partial Y_{m,p}}
  & j = p-1 \\[12pt]
\displaystyle
-\frac{\zeta_j}{\rho_j}
  \left(\frac{\partial j_{k,\,p+1/2}}{\partial Y_{m,p}}
        - \frac{\partial j_{k,\,p-1/2}}{\partial Y_{m,p}}\right)
  + \frac{W_k}{\rho_j}\frac{\partial \dot\omega_k}{\partial Y_m}
  + R_{km}^{(p)} + C_{km}^{(p)}
  & j = p \\[12pt]
\displaystyle
+\frac{\zeta_j}{\rho_j} \frac{\partial j_{k,\,p+1/2}}{\partial Y_{m,p}}
  & j = p+1
\end{cases}
$$

Here $R_{km}^{(p)}$ is the *density chain rule* term: because
$\rho_j = P\overline{W}_j/(RT_j)$ depends on composition through $\overline{W}_j$,
the $1/\rho_j$ prefactor contributes

$$
R_{km}^{(p)} = \bigl(W_k\dot\omega_{k,p} - d_{k,p}\bigr)\,
  \frac{\overline{W}_p}{W_m\,\rho_p},
\qquad
d_{k,p} = \frac{2\bigl(j_{k,\,p+1/2} - j_{k,\,p-1/2}\bigr)}{z_{p+1} - z_{p-1}},
$$

obtained from $\partial(1/\rho_p)/\partial Y_{m,p} = \overline{W}_p / (W_m\,\rho_p)$
acting on the reaction and diffusion parts of the residual numerator. (The convection
part has no density chain rule contribution because the $\rho_j$ in $\rho_j u_j$ cancels
the $1/\rho_j$ prefactor.) The convection term $C_{km}^{(p)}$ is the derivative of the
upwinded $-u_j\,\partial Y_k/\partial z|_j$: it is species-diagonal ($k=m$) and
contributes to whichever rows have point $p$ in their upwind difference stencil.

Boundary rows $j = 0$ and $j = N-1$ are skipped: those residuals are fixed-value
constraints with no dependence on interior species values.

## Diffusive-flux derivatives

The two flux derivatives above are evaluated at the endpoints of the relevant interval.
With frozen transport, differentiating $j_{k,\,j+1/2}$ with respect to the mass fraction
at its _left_ endpoint (point $j$) gives

$$
\frac{\partial j_{k,\,j+1/2}}{\partial Y_{m,j}} =
  \frac{\tilde{D}_{k,\,j+1/2}}{\Delta z_j}\, f_j\,\bigl(\delta_{km} - X_{k,j}\bigr)
  + Y_{k,j}\,\frac{\partial S_{j+1/2}}{\partial Y_{m,j}}
  + \delta_{km}\, S_{j+1/2},
$$

and with respect to the _right_ endpoint (point $j+1$),

$$
\frac{\partial j_{k,\,j+1/2}}{\partial Y_{m,j+1}} =
  -\frac{\tilde{D}_{k,\,j+1/2}}{\Delta z_j}\, f_{j+1}\,\bigl(\delta_{km} - X_{k,j+1}\bigr)
  + Y_{k,j}\,\frac{\partial S_{j+1/2}}{\partial Y_{m,j+1}},
$$

where the $f_j = \overline{W}_j / W_{m,j}$ is the $\partial X/\partial Y$ scale factor
and the correction-flux derivatives are

$$
\frac{\partial S_{j+1/2}}{\partial Y_{m,j}} & =
  -f_j\!\left(\frac{\tilde{D}_{m,\,j+1/2}}{\Delta z_j} - a_j\right),

\frac{\partial S_{j+1/2}}{\partial Y_{m,j+1}} & =
  +f_{j+1}\!\left(\frac{\tilde{D}_{m,\,j+1/2}}{\Delta z_j} - a_{j+1}\right),
$$

with $a_j = \sum_n (\tilde{D}_{n,\,j+1/2}/\Delta z_j)\, X_{n,j}$. The right-endpoint
form has no $\delta_{km}\,S$ term because the $Y_{k,j}\,S_{j+1/2}$ part of the flux is
anchored at the left endpoint. The implementation computes only the endpoint actually
needed for each row.

## Reaction derivatives: $\partial\dot\omega_k / \partial Y_m$

At constant $T$ and $P$, the chain rule from molar concentrations $C_n = \rho Y_n / W_n$
to the unnormalized mass fractions gives

$$
\frac{\partial \dot\omega_k}{\partial Y_m} =
  \frac{\rho}{W_m}\frac{\partial \dot\omega_k}{\partial C_m}
  - \frac{\overline{W}}{W_m} \sum_n \frac{\partial \dot\omega_k}{\partial C_n}\, C_n.
$$

The first term is the direct response to perturbing $C_m$; the second arises because
$\rho$ and $\overline{W}$ themselves depend on the unnormalized $Y_m$. The sparse matrix
$\partial \dot\omega / \partial C$ is supplied by
{ct}`Kinetics::netProductionRates_ddCi()`, whose nonzero pattern is fixed for a given
mechanism and is reused across all grid points.

## Energy and flow rows

The energy-row block $\partial F_{T,j} / \partial Y_{m,p}$
(for $j \in \{p-1,\,p,\,p+1\}$) reuses the same diffusive-flux derivatives through the
enthalpy-flux term $\sum_k h_k\, j_k/W_k$, and adds a density chain rule correction
together with a reaction-heat term and (when enabled) a radiation term at $j = p$. The
continuity and radial-momentum rows depend on the species mass fractions only through
the density, so they receive only a density chain rule correction
($\partial \rho_p / \partial Y_{m,p} = -\rho_p\,\overline{W}_p / W_m$) from the column
at point $p$.
