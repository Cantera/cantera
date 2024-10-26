//! @file StFlow.cpp

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#include "cantera/oneD/StFlow.h"
#include "cantera/base/global.h"

using namespace std;

namespace Cantera
{

StFlow::StFlow(ThermoPhase* ph, size_t nsp, size_t points) :
    Flow1D(ph, nsp, points)
{
    warn_deprecated("StFlow::StFlow",
        "To be removed after Cantera 3.1. Class replaced by Flow1D.");
}

StFlow::StFlow(shared_ptr<ThermoPhase> th, size_t nsp, size_t points)
    : Flow1D(th, nsp, points)
{
    warn_deprecated("StFlow::StFlow",
        "To be removed after Cantera 3.1. Class replaced by Flow1D.");
}

StFlow::StFlow(shared_ptr<Solution> sol, const string& id, size_t points)
    : Flow1D(sol, id, points)
{
    warn_deprecated("StFlow::StFlow",
        "To be removed after Cantera 3.1. Class replaced by Flow1D.");
}

void StFlow::eval(size_t jg, double* xg, double* rg, integer* diagg, double rdt)
{
    // if evaluating a Jacobian, and the global point is outside the domain of
    // influence for this domain, then skip evaluating the residual
    if (jg != npos && (jg + 1 < firstPoint() || jg > lastPoint() + 1)) {
        return;
    }

    // start of local part of global arrays
    double* x = xg + loc();
    double* rsd = rg + loc();
    integer* diag = diagg + loc();

    size_t jmin, jmax;
    if (jg == npos) { // evaluate all points
        jmin = 0;
        jmax = m_points - 1;
    } else { // evaluate points for Jacobian
        size_t jpt = (jg == 0) ? 0 : jg - firstPoint();
        jmin = std::max<size_t>(jpt, 1) - 1;
        jmax = std::min(jpt+1,m_points-1);
    }

    updateProperties(jg, x, jmin, jmax);
    evalResidual(x, rsd, diag, rdt, jmin, jmax);
    evalUo(x, rsd, diag, rdt, jmin, jmax);
}

void StFlow::evalResidual(double* x, double* rsd, int* diag,
                          double rdt, size_t jmin, size_t jmax)
{
    //----------------------------------------------------
    // evaluate the residual equations at all required
    // grid points
    //----------------------------------------------------

    // calculation of qdotRadiation (see docstring of enableRadiation)
    if (m_do_radiation) {
        // variable definitions for the Planck absorption coefficient and the
        // radiation calculation:
        double k_P_ref = 1.0*OneAtm;

        // polynomial coefficients:
        const double c_H2O[6] = {-0.23093, -1.12390, 9.41530, -2.99880,
                                     0.51382, -1.86840e-5};
        const double c_CO2[6] = {18.741, -121.310, 273.500, -194.050,
                                     56.310, -5.8169};

        // calculation of the two boundary values
        double boundary_Rad_left = m_epsilon_left * StefanBoltz * pow(T(x, 0), 4);
        double boundary_Rad_right = m_epsilon_right * StefanBoltz * pow(T(x, m_points - 1), 4);

        // loop over all grid points
        for (size_t j = jmin; j < jmax; j++) {
            // helping variable for the calculation
            double radiative_heat_loss = 0;

            // calculation of the mean Planck absorption coefficient
            double k_P = 0;
            // absorption coefficient for H2O
            if (m_kRadiating[1] != npos) {
                double k_P_H2O = 0;
                for (size_t n = 0; n <= 5; n++) {
                    k_P_H2O += c_H2O[n] * pow(1000 / T(x, j), (double) n);
                }
                k_P_H2O /= k_P_ref;
                k_P += m_press * X(x, m_kRadiating[1], j) * k_P_H2O;
            }
            // absorption coefficient for CO2
            if (m_kRadiating[0] != npos) {
                double k_P_CO2 = 0;
                for (size_t n = 0; n <= 5; n++) {
                    k_P_CO2 += c_CO2[n] * pow(1000 / T(x, j), (double) n);
                }
                k_P_CO2 /= k_P_ref;
                k_P += m_press * X(x, m_kRadiating[0], j) * k_P_CO2;
            }

            // calculation of the radiative heat loss term
            radiative_heat_loss = 2 * k_P *(2 * StefanBoltz * pow(T(x, j), 4)
            - boundary_Rad_left - boundary_Rad_right);

            // set the radiative heat loss vector
            m_qdotRadiation[j] = radiative_heat_loss;
        }
    }

    for (size_t j = jmin; j <= jmax; j++) {
        //----------------------------------------------
        //         left boundary
        //----------------------------------------------

        if (j == 0) {
            // these may be modified by a boundary object

            // Continuity. This propagates information right-to-left, since
            // rho_u at point 0 is dependent on rho_u at point 1, but not on
            // mdot from the inlet.
            rsd[index(c_offset_U,0)] =
                -(rho_u(x,1) - rho_u(x,0))/m_dz[0]
                -(density(1)*V(x,1) + density(0)*V(x,0));

            // the inlet (or other) object connected to this one will modify
            // these equations by subtracting its values for V, T, and mdot. As
            // a result, these residual equations will force the solution
            // variables to the values for the boundary object
            rsd[index(c_offset_V,0)] = V(x,0);
            rsd[index(c_offset_T,0)] = T(x,0);
            if (m_usesLambda) {
                rsd[index(c_offset_L, 0)] = -rho_u(x, 0);
            } else {
                rsd[index(c_offset_L, 0)] = lambda(x, 0);
                diag[index(c_offset_L, 0)] = 0;
            }

            // The default boundary condition for species is zero flux. However,
            // the boundary object may modify this.
            double sum = 0.0;
            for (size_t k = 0; k < m_nsp; k++) {
                sum += Y(x,k,0);
                rsd[index(c_offset_Y + k, 0)] =
                    -(m_flux(k,0) + rho_u(x,0)* Y(x,k,0));
            }
            rsd[index(c_offset_Y + leftExcessSpecies(), 0)] = 1.0 - sum;

            // set residual of poisson's equ to zero
            rsd[index(c_offset_E, 0)] = x[index(c_offset_E, j)];
        } else if (j == m_points - 1) {
            evalRightBoundary(x, rsd, diag, rdt);
        } else { // interior points
            evalContinuity(j, x, rsd, diag, rdt);
            // set residual of poisson's equ to zero
            rsd[index(c_offset_E, j)] = x[index(c_offset_E, j)];

            //------------------------------------------------
            //    Radial momentum equation
            //
            //    \rho dV/dt + \rho u dV/dz + \rho V^2
            //       = d(\mu dV/dz)/dz - lambda
            //-------------------------------------------------
            if (m_usesLambda) {
                rsd[index(c_offset_V,j)] =
                    (shear(x, j) - lambda(x, j) - rho_u(x, j) * dVdz(x, j)
                    - m_rho[j] * V(x, j) * V(x, j)) / m_rho[j]
                    - rdt * (V(x, j) - V_prev(j));
                diag[index(c_offset_V, j)] = 1;
            } else {
                rsd[index(c_offset_V, j)] = V(x, j);
                diag[index(c_offset_V, j)] = 0;
            }

            //-------------------------------------------------
            //    Species equations
            //
            //   \rho dY_k/dt + \rho u dY_k/dz + dJ_k/dz
            //   = M_k\omega_k
            //-------------------------------------------------
            getWdot(x,j);
            for (size_t k = 0; k < m_nsp; k++) {
                double convec = rho_u(x,j)*dYdz(x,k,j);
                double diffus = 2.0*(m_flux(k,j) - m_flux(k,j-1))
                                / (z(j+1) - z(j-1));
                rsd[index(c_offset_Y + k, j)]
                = (m_wt[k]*(wdot(k,j))
                   - convec - diffus)/m_rho[j]
                  - rdt*(Y(x,k,j) - Y_prev(k,j));
                diag[index(c_offset_Y + k, j)] = 1;
            }

            //-----------------------------------------------
            //    energy equation
            //
            //    \rho c_p dT/dt + \rho c_p u dT/dz
            //    = d(k dT/dz)/dz
            //      - sum_k(\omega_k h_k_ref)
            //      - sum_k(J_k c_p_k / M_k) dT/dz
            //-----------------------------------------------
            if (m_do_energy[j]) {

                setGas(x,j);
                double dtdzj = dTdz(x,j);
                double sum = 0.0;

                grad_hk(x, j);
                for (size_t k = 0; k < m_nsp; k++) {
                    double flxk = 0.5*(m_flux(k,j-1) + m_flux(k,j));
                    sum += wdot(k,j)*m_hk(k,j);
                    sum += flxk * m_dhk_dz(k,j) / m_wt[k];
                }

                rsd[index(c_offset_T, j)] = - m_cp[j]*rho_u(x,j)*dtdzj
                                            - conduction(x,j) - sum;
                rsd[index(c_offset_T, j)] /= (m_rho[j]*m_cp[j]);
                rsd[index(c_offset_T, j)] -= rdt*(T(x,j) - T_prev(j));
                rsd[index(c_offset_T, j)] -= (m_qdotRadiation[j] / (m_rho[j] * m_cp[j]));
                diag[index(c_offset_T, j)] = 1;
            } else {
                // residual equations if the energy equation is disabled
                rsd[index(c_offset_T, j)] = T(x,j) - T_fixed(j);
                diag[index(c_offset_T, j)] = 0;
            }

            if (m_usesLambda) {
                rsd[index(c_offset_L, j)] = lambda(x, j) - lambda(x, j - 1);
            } else {
                rsd[index(c_offset_L, j)] = lambda(x, j);
            }
            diag[index(c_offset_L, j)] = 0;
        }
    }
}

void StFlow::evalRightBoundary(double* x, double* rsd, int* diag, double rdt)
{
    size_t j = m_points - 1;

    // the boundary object connected to the right of this one may modify or
    // replace these equations. The default boundary conditions are zero u, V,
    // and T, and zero diffusive flux for all species.

    rsd[index(c_offset_V,j)] = V(x,j);
    diag[index(c_offset_V,j)] = 0;
    double sum = 0.0;
    // set residual of poisson's equ to zero
    rsd[index(c_offset_E, j)] = x[index(c_offset_E, j)];
    for (size_t k = 0; k < m_nsp; k++) {
        sum += Y(x,k,j);
        rsd[index(k+c_offset_Y,j)] = m_flux(k,j-1) + rho_u(x,j)*Y(x,k,j);
    }
    rsd[index(c_offset_Y + rightExcessSpecies(), j)] = 1.0 - sum;
    diag[index(c_offset_Y + rightExcessSpecies(), j)] = 0;
    if (m_usesLambda) {
        rsd[index(c_offset_U, j)] = rho_u(x, j);
    } else {
        rsd[index(c_offset_U, j)] = rho_u(x, j) - rho_u(x, j-1);
    }

    rsd[index(c_offset_L, j)] = lambda(x, j) - lambda(x, j-1);
    diag[index(c_offset_L, j)] = 0;
    rsd[index(c_offset_T, j)] = T(x, j);
}

void StFlow::evalContinuity(size_t j, double* x, double* rsd, int* diag, double rdt)
{
    //algebraic constraint
    diag[index(c_offset_U, j)] = 0;
    //----------------------------------------------
    //    Continuity equation
    //
    //    d(\rho u)/dz + 2\rho V = 0
    //----------------------------------------------
    if (m_usesLambda) {
        // Note that this propagates the mass flow rate information to the left
        // (j+1 -> j) from the value specified at the right boundary. The
        // lambda information propagates in the opposite direction.
        rsd[index(c_offset_U,j)] =
            -(rho_u(x,j+1) - rho_u(x,j))/m_dz[j]
            -(density(j+1)*V(x,j+1) + density(j)*V(x,j));
    } else if (m_isFree) {
        // terms involving V are zero as V=0 by definition
        if (z(j) > m_zfixed) {
            rsd[index(c_offset_U,j)] =
                - (rho_u(x,j) - rho_u(x,j-1))/m_dz[j-1];
        } else if (z(j) == m_zfixed) {
            if (m_do_energy[j]) {
                rsd[index(c_offset_U,j)] = (T(x,j) - m_tfixed);
            } else {
                rsd[index(c_offset_U,j)] = (rho_u(x,j)
                                            - m_rho[0]*0.3); // why 0.3?
            }
        } else if (z(j) < m_zfixed) {
            rsd[index(c_offset_U,j)] =
                - (rho_u(x,j+1) - rho_u(x,j))/m_dz[j];
        }
    } else {
        // unstrained with fixed mass flow rate
        rsd[index(c_offset_U, j)] = rho_u(x, j) - rho_u(x, j - 1);
    }
}

} // namespace
