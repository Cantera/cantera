# Analytic Jacobian for 1D Flames

## Overview

The [nonlinear solver](nonlinear-solver) for 1D flames evaluates the Jacobian
$\mathbf{J}$ of the residual function at each Newton iteration. By default this is done
by finite differences (perturbing each solution variable in turn), which requires $N_v$
residual evaluations for a solution vector of length $N_v$. The analytic Jacobian mode
(`domain.jacobian_mode = "analytic"`) instead computes the species-column blocks of
$\mathbf{J}$ directly from the kinetics concentration derivatives, reducing the number of
residual evaluations per Jacobian to $N_v - K \cdot (N-2)$, where $K$ is the number of
species and $N$ is the number of grid points.

The derivation below covers the interior grid points $j = 1, \ldots, N-2$ (the two
boundary points always use finite differences). All transport coefficients are treated as
frozen at the base state, matching the frozen-transport approximation of the FD Jacobian.

## Species residual

The discretized species mass-fraction residual at interior point $j$ (using an upwind
convection scheme and a central-difference diffusion flux) is:

$$
F_k^{(j)} = \rho_j \frac{\partial Y_k}{\partial t}
  - \frac{2(\bar{F}_k^{(j)} - \bar{F}_k^{(j-1)})}{z_{j+1} - z_{j-1}}
  - W_k \dot\omega_k^{(j)}
  - [\text{convection}]_k^{(j)}
$$

where $\bar{F}_k^{(j)} = \tfrac{1}{2}(F_k^{(j)} + F_k^{(j+1)})$ is the mean
diffusive flux over the interval $(z_j, z_{j+1})$ and $F_k^{(j)}$ is the mixture-
averaged diffusive mass flux at the right boundary of interval $j$:

$$
F_k^{(j)} = \frac{1}{\Delta z_j}\Bigl[
  \bar{c}_j D_{kj} \bigl(X_k^{(j)} - X_k^{(j+1)}\bigr)
  + Y_k^{(j)} S_j
\Bigr]
$$

with $S_j = -\sum_n \tfrac{\bar{c}_j D_{nj}}{\Delta z_j}(X_n^{(j)} - X_n^{(j+1)})$.

## Jacobian column for $Y_m$ at point $p$

The column $\partial F / \partial Y_m(p)$ has nonzero blocks only in rows $j = p-1, p,
p+1$ (because only the two adjacent diffusive fluxes $F^{(p-1)}$ and $F^{(p)}$ depend
on $Y(p)$, plus the reaction term at $j = p$). The three blocks are:

$$
\frac{\partial F_k^{(j)}}{\partial Y_m(p)} =
\begin{cases}
-\frac{2}{\rho_{j} (z_{j+1}-z_{j-1})} \frac{\partial F_k^{(p-1)}}{\partial Y_m(p)}
  & j = p-1 \\[4pt]
-\frac{2}{\rho_{j} (z_{j+1}-z_{j-1})}
  \Bigl(\frac{\partial F_k^{(p)}}{\partial Y_m(p)}
        - \frac{\partial F_k^{(p-1)}}{\partial Y_m(p)}\Bigr)
  + \frac{W_k}{\rho_j}\frac{\partial \dot\omega_k}{\partial Y_m}
  + [\text{density chain}]_{km}^{(p)}
  & j = p \\[4pt]
+\frac{2}{\rho_{j} (z_{j+1}-z_{j-1})} \frac{\partial F_k^{(p)}}{\partial Y_m(p)}
  & j = p+1
\end{cases}
$$

Boundary rows $j = 0$ and $j = N-1$ are skipped because those residuals are boundary
conditions that do not depend on interior species values.

## Diffusive-flux Jacobian

With frozen transport coefficients, the molar-gradient flux Jacobian at interval $q$
with respect to $Y_m$ at endpoint $q$ (left, $\text{AtQ}=\text{true}$) is:

$$
\frac{\partial F_k^{(q)}}{\partial Y_m(q)} =
  \frac{\bar{c} D_{kq}}{\Delta z_q} f_q \bigl(\delta_{km} - X_k^{(q)}\bigr)
  + Y_k^{(q)} \frac{\partial S_q}{\partial Y_m(q)}
  + \delta_{km} S_q
$$

where $f_q = \bar{w}_q / W_m$ is the $\partial X / \partial Y$ scale factor,
$\partial S_q / \partial Y_m(q) = -f_q(D_{mq}/\Delta z_q - a_q)$, and
$a_q = \sum_n (D_{nq}/\Delta z_q) X_n^{(q)}$.

The right-endpoint derivative (at $q+1$) has an analogous formula with $f_{q+1}$,
$a_{q+1}$, and a sign flip. The implementation uses a template parameter to compute
only the needed endpoint, avoiding redundant work.

## Reaction derivatives: $\partial\dot\omega_k / \partial Y_m$

At constant $T$ and $P$, the chain rule from mole fractions to unnormalized mass
fractions gives:

$$
\frac{\partial \dot\omega_k}{\partial Y_m} =
  \frac{\rho}{W_m}\frac{\partial \dot\omega_k}{\partial C_m}
  - \frac{\bar{w}}{W_m} \sum_n \frac{\partial \dot\omega_k}{\partial C_n} C_n
$$

where $C_n = \rho Y_n / W_n$ is the molar concentration. The sparse matrix
$\partial \dot\omega / \partial C$ is provided by
`Kinetics::netProductionRates_ddCi()`, whose nonzero pattern is fixed for a given
mechanism and is reused across grid points.

## Energy and flow residual columns

The energy-row block $\partial F_T^{(j)} / \partial Y_m(p)$ for $j \in \{p-1, p, p+1\}$
is computed from the same diffusive-flux Jacobian blocks (via the enthalpy-flux term
$\sum_k h_k F_k / W_k$) plus a density-chain correction and a reaction-heat term at
$j = p$. The continuity and momentum rows receive only density-chain corrections
($\partial \rho / \partial Y_m = -\rho \bar{w} / W_m$) from the same column $p$.
