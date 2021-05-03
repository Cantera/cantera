//! @file Newton.cpp: A damped & bounded quasi-Newton solver

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#include "cantera/numerics/Newton.h"
#include "cantera/base/utilities.h"

#include <ctime>

using namespace std;

namespace Cantera
{

// constants
const double damp_factor = sqrt(2.0);
const double damp_min = 0.1;

Newton::Newton(FuncEval& func) {
    m_residfunc = &func;
    m_nv = m_residfunc->neq();
    m_x.resize(m_nv);
    m_x1.resize(m_nv);
    m_stp.resize(m_nv);
    m_stp1.resize(m_nv);
    m_upper_bounds.resize(m_nv, 0.0);
    m_lower_bounds.resize(m_nv, 0.0);
    m_rtol_ss.resize(m_nv, 1.0e-4);
    m_atol_ss.resize(m_nv, 1.0e-9);
    // m_rtol_ts.resize(m_nv, 1.0e-4);
    // m_atol_ts.resize(m_nv, 1.0e-11);
    m_xlast.resize(m_nv);
    m_xsave.resize(m_nv);

    m_converge_tol = 1.0e-14;
    m_rdt = 0;

    m_jacobian = DenseMatrix(m_nv, m_nv);
    m_jacAge = npos;
    m_jacMaxAge = 1;
    m_jacRtol = 1.0e-15;
    m_jacAtol = sqrt(std::numeric_limits<double>::epsilon());
}

void Newton::evalJacobian(doublereal* x, doublereal* xdot) {

    // calculate unperturbed residual
    m_residfunc->eval(0, x, xdot, 0);

    for (size_t n = 0; n < m_nv; n++) {
        // calculate the nth Jacobian column
        double xsave = x[n];

        // calculate the perturbation amount, preserving the sign of x[n]
        double dx;
        if (xsave >= 0) {
            dx = xsave*m_jacRtol + m_jacAtol;
        } else {
            dx = xsave*m_jacRtol - m_jacAtol;
        }

        // perturb the solution vector
        x[n] = xsave + dx;
        dx = x[n] - xsave;

        // calculate perturbed residual
        vector_fp xdotPerturbed(m_nv); //make this member for speed?
        m_residfunc->eval(0, x, xdotPerturbed.data(), 0);

        // compute nth column of Jacobian
        for (size_t m = 0; m < m_nv; m++) {
            m_jacobian.value(m,n) = (xdotPerturbed[m] - xdot[m])/dx;
        }
        // restore solution vector
        x[n] = xsave;

        // for constant-valued components: 1 in diagonal position, all 0's in row and column
        for (size_t i : m_constantComponents) {
            for (size_t j = 0; j < m_nv; j++) {
                m_jacobian.value(i,j) = 0;
                m_jacobian.value(j,i) = 0;
            }
            m_jacobian.value(i,i) = 1;
        }
    }

    // writelog("\nnew jac:\n");
    // for (int i = 0; i < m_nv; i++) {
    //     writelog("ROW {} | ", i);
    //     for (int j = 0; j < m_nv; j++) {
    //         writelog("{:.5} ", m_jacobian.value(i,j));
    //     }
    //     writelog("\n\n");
    // }

    m_jacAge = 0;

    m_jacFactored = m_jacobian;
    Cantera::factor(m_jacFactored);
}

// RMSD (weighted)
doublereal Newton::weightedNorm(const doublereal* x, const doublereal* step) const
{
    double square = 0.0;
    for (size_t i = 0; i < m_nv; i++) {
        square += pow(step[i]/(x[i] + m_atol_ss[i]), 2);
    }
    return sqrt(square/m_nv);
}

void Newton::step(doublereal* x, doublereal* step, int loglevel)
{
    m_residfunc->eval(0, x, step, 0);
    for (size_t n = 0; n < m_nv; n++) {
        step[n] = -step[n] + m_rdt*(x[n]-m_xlast[n]);
    }

    DenseMatrix solvejac = m_jacFactored;
    try {
        //Note: this function takes an unfactored jacobian, then finds its LU factorization before
        // solving. Optimization is possible by saving the factored jacobian, since it can be reused.
        // Also, the DenseMatrix provided here will be overwritten with the LU factored version, so
        // a copy is passed instead in order to preserve the original for reuse.
        // Cantera::solve(solvejac, step);
        Cantera::solveFactored(solvejac, step);

    } catch (CanteraError&) {
        // int iok = m_jac->info() - 1;
        int iok = -1; //TODO: enable error info
        if (iok >= 0) {
            throw CanteraError("Newton::step",
                "Jacobian is singular for component {} (Matrix row {})",
                iok, iok); //TODO: add component name
        } else {
            throw;
        }
    }
}

int Newton::hybridSolve() {
    int MAX = 100;
    int newtonsolves = 0;
    int timesteps = 0;

    m_residfunc->getState(m_x.data());
    copy(m_x.begin(), m_x.end(), m_xsave.begin());

    for(int i = 0; i < MAX; i++) {
        copy(m_xsave.begin(), m_xsave.end(), m_x.begin());
        for(int j = 0; j < MAX; j++) {
            timesteps++;
            timestep();
        }
        copy(m_x.begin(), m_x.end(), m_xsave.begin());
        newtonsolves++;
        m_rdt = 0;
        if(solve() == 1) {
            writelog("\nConverged in {} newton solves, {} timesteps.", newtonsolves, timesteps);
            return 1;
        }
    }
    writelog("Solver failure...");
    return 0;
}

int Newton::timestep() {
    m_rdt = 1.0/(1.0e-5);

    // calculate time-integration Jacobian
    // m_residfunc->getState(m_x.data());
    copy(m_x.begin(), m_x.end(), m_xlast.begin());
    evalJacobian(&m_x[0], &m_stp[0]);
    for (size_t i = 0; i < m_nv; i++) {
        m_jacobian.value(i,i) -= m_rdt;
    }
    return solve();
}

int Newton::solve(int loglevel)
{
    m_residfunc->getState(m_x.data());

    while (true) {
        // Check whether the Jacobian should be re-evaluated.
        if (m_jacAge > m_jacMaxAge) {
            evalJacobian(&m_x[0], &m_stp[0]);
        }

        // compute the undamped Newton step
        step(&m_x[0], &m_stp[0], loglevel-1);
        m_jacAge++;

        // compute the weighted norm of the undamped step size step0
        double step_rms = weightedNorm(&m_x[0], &m_stp[0]);

        // compute the multiplier to keep all components in bounds.
        double bound_factor = 1.0;
        for (size_t i = 0; i < m_nv; i++) {
            double upper_bound = m_upper_bounds[i];
            double lower_bound = m_lower_bounds[i];
            double val = m_x[i];
            if (val > upper_bound + 1.0e-12 || val < lower_bound - 1.0e-12) {
                throw CanteraError("Newton::dampStep", "solution out of bounds");
            }
            double newval = val + m_stp[i];

            if (newval > upper_bound) {
                bound_factor = max(0.0, min(bound_factor, (upper_bound - val)/(newval - val)));
            } else if (newval < lower_bound) {
                bound_factor = min(bound_factor, (val - lower_bound)/(val - newval));
            }
        }
        // if bound factor is very small, then x0 is already close to the boundary and
        // step0 points out of the allowed domain. In this case, the Newton
        // algorithm fails, so return an error condition.
        if (bound_factor < 1.e-10) {
            throw CanteraError("Newton::dampStep", "solution at limits");
        }

        int m = -1;
        // damped step - attempt to find a damping coefficient such that the next
        // undamped step would have a RMSD smaller than that of step0
        for (double damp = bound_factor; damp > damp_min; damp /= damp_factor) {
            // step the solution by the damped step size
            for (size_t j = 0; j < m_nv; j++) {
                m_x1[j] = damp*m_stp[j] + m_x[j];
            }
            // compute the next undamped step that would result if x1 is accepted
            step(&m_x1[0], &m_stp1[0], loglevel-1);

            // compute the weighted norm of step1
            double nextstep_rms = weightedNorm(&m_x1[0], &m_stp1[0]);

            // converged solution criteria
            if (nextstep_rms < m_converge_tol) {
                return 1;
            }
            // Also accept it if this step would result in a
            // converged solution - successful step, but not converged yet.
            // Take the damped step and try again.
            if (nextstep_rms < step_rms) {
                m = 0;
                copy(m_x1.begin(), m_x1.end(), m_x.begin());
                break;
            }
        }

        if (m < 0) {
            // If dampStep fails, first try a new Jacobian if an old one was
            // being used. If it was a new Jacobian, then return -1 to signify
            // failure.
            if (m_jacAge > 1) {
                m_jacAge = npos;
            } else {
                return m;
            }
        }
    }
}

} // end namespace Cantera
