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
const doublereal DampFactor = sqrt(2.0);
const size_t NDAMP = 7;

Newton::Newton(FuncEval& func) {
    m_residfunc = &func;
    m_nv = m_residfunc->neq();
    m_constant.resize(m_nv, false);
    m_x.resize(m_nv);
    m_x1.resize(m_nv);
    m_stp.resize(m_nv);
    m_stp1.resize(m_nv);
    m_max.resize(m_nv, 0.0);
    m_min.resize(m_nv, 0.0);
    m_rtol_ss.resize(m_nv, 1.0e-4);
    m_atol_ss.resize(m_nv, 1.0e-9);
    // m_rtol_ts.resize(m_nv, 1.0e-4);
    // m_atol_ts.resize(m_nv, 1.0e-11);

    m_jacobian = DenseMatrix(m_nv, m_nv);
    m_jacAge = 10000;
    m_jacMaxAge = 5;
    m_jacRtol = 1.0e-5;
    m_jacAtol = sqrt(std::numeric_limits<double>::epsilon());
}

void Newton::evalJacobian(doublereal* x, doublereal* xdot) {

    // // calculate unperturbed residual
    m_residfunc->eval(0, x, xdot, 0);

    for (size_t n = 0; n < m_nv; n++) {
        // calculate the nth Jacobian column, unless component n is constant
        if (!m_constant[n]) {
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

        // for constant components: Jacobian column of 0's, with 1 on diagonal
        // note: Jacobian row must also be 0's w/ 1 on diagonal
        } else {
            for (size_t m = 0; m < m_nv; m++) {
                m_jacobian.value(m,n) = 0;
            }
            m_jacobian.value(n,n) = 1;
        }
    }

    // writelog("\nnew jac:\n");
    // for (int i = 0; i < m_nv; i++) {
    //     for (int j = 0; j < m_nv; j++) {
    //         writelog("{:14.5} ", m_jacobian.value(i,j));
    //     }
    //     writelog("\n");
    // }

    m_jacAge = 0;
}

doublereal Newton::weightedNorm(const doublereal* x, const doublereal* step) const
{
    double sum = 0.0;
    for (size_t n = 0; n < m_nv; n++) {
        double weight = m_rtol_ss[n]*fabs(x[n]) + m_atol_ss[n]; //TODO: add transient tolerances if rdt!=0
        double f = step[n]/weight;
        sum += f*f;
    }
    return sqrt(sum/m_nv);
}

void Newton::step(doublereal* x, doublereal* step, int loglevel)
{
    m_residfunc->eval(0, x, step, 0);
    for (size_t n = 0; n < m_nv; n++) {
        step[n] = -step[n];
    }

    DenseMatrix solvejac = m_jacobian;
    try {
        //Note: this function takes an unfactored jacobian, then finds its LU factorization before
        // solving. Optimization is possible by saving the factored jacobian, since it can be reused.
        // Also, the DenseMatrix provided here will be overwritten with the LU factored version, so
        // a copy is passed instead in order to preserve the original for reuse.
        Cantera::solve(solvejac, step);

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

doublereal Newton::boundStep(const doublereal* x, const doublereal* step, int loglevel)
{
    doublereal boundFactor = 1.0;
    bool wroteTitle = false;
    for (size_t n = 0; n < m_nv; n++) {
        double upperBound = m_max[n];
        double lowerBound = m_min[n];
        double val = x[n];
        if (loglevel > 0 && (val > upperBound + 1.0e-12 || val < lowerBound - 1.0e-12)) {
            writelog("\nERROR: solution out of bounds.\n");
            writelog("Component {}: {:10.3e} with bounds ({:10.3e}, {:10.3e})\n",
                        n, val, lowerBound, upperBound);
            // writelog("domain {:d}: {:>20s}({:d}) = {:10.3e} ({:10.3e}, {:10.3e})\n",
            //             r.domainIndex(), r.componentName(m), j, val, below, above);
        }
        double newval = val + step[n];

        if (newval > upperBound) {
            boundFactor = std::max(0.0, std::min(boundFactor, (upperBound - val)/(newval - val)));
        } else if (newval < lowerBound) {
            boundFactor = std::min(boundFactor, (val - lowerBound)/(val - newval));
        }
        if (loglevel > 1 && (newval > upperBound || newval < lowerBound)) {
            if (!wroteTitle) {
                writelog("\nNewton step takes solution out of bounds.\n\n");
                // writelog("  {:>12s}  {:>12s}  {:>4s}  {:>10s}  {:>10s}  {:>10s}  {:>10s}\n",
                //          "domain","component","pt","value","step","min","max");
                wroteTitle = true;
            }
            writelog("Component {}: {:10.3e} with bounds ({:10.3e}, {:10.3e}), step = {:10.3e}\n",
                        n, val, lowerBound, upperBound, step[n]);
            // writelog("          {:4d}  {:>12s}  {:4d}  {:10.3e}  {:10.3e}  {:10.3e}  {:10.3e}\n",
            //          r.domainIndex(), r.componentName(m), j,
            //          val, step[index(m,j)], below, above);
        }
    }
    return boundFactor;
}

int Newton::dampStep(const doublereal* x0, const doublereal* step0,
                          doublereal* x1, doublereal* step1, doublereal& s1,
                          int loglevel, bool writetitle)
{
    // write header
    if (loglevel > 0 && writetitle) {
        writelog("\n\nDamped Newton iteration:\n");
        writeline('-', 65, false);

        writelog("\n{}  {:>9s}   {:>9s}     {:>9s}   {:>9s}   {:>9s}  {:>5s} {:>5s}\n",
                "m","F_damp","F_bound","log10(ss)",
                "log10(s0)","log10(s1)","N_jac","Age");
        writeline('-', 65);
    }

    // compute the weighted norm of the undamped step size step0
    doublereal s0 = weightedNorm(x0, step0);

    // compute the multiplier to keep all components in bounds
    doublereal boundFactor = boundStep(x0, step0, loglevel-1);

    // if bound factor is very small, then x0 is already close to the boundary and
    // step0 points out of the allowed domain. In this case, the Newton
    // algorithm fails, so return an error condition.
    if (boundFactor < 1.e-10) {
        debuglog("\nAt limits.\n", loglevel);
        return -3;
    }

    // ---------- Attempt damped step ----------

    // damping coefficient starts at 1.0
    doublereal damp = 1.0;
    size_t m;
    for (m = 0; m < NDAMP; m++) {
        double ff = boundFactor*damp;

        // step the solution by the damped step size
        for (size_t j = 0; j < m_nv; j++) {
            x1[j] = ff*step0[j] + x0[j];
        }

        // compute the next undamped step that would result if x1 is accepted
        step(x1, step1, loglevel-1);

        // compute the weighted norm of step1
        s1 = weightedNorm(x1, step1);

        // write log information
        if (loglevel > 0) {
            doublereal ss = weightedNorm(x1,step1);
            writelog("\n{:d}  {:9.5f}   {:9.5f}   {:9.5f}   {:9.5f}   {:9.5f}  {:d}/{:d}",
                     m, damp, boundFactor, log10(ss+SmallNumber),
                     log10(s0+SmallNumber), log10(s1+SmallNumber),
                     m_jacAge, m_jacMaxAge);
        }

        // if the norm of s1 is less than the norm of s0, then accept this
        // damping coefficient. Also accept it if this step would result in a
        // converged solution. Otherwise, decrease the damping coefficient and
        // try again.
        if (s1 < 1.0 || s1 < s0) {
            break;
        }
        damp /= DampFactor;
    }

    // If a damping coefficient was found, return 1 if the solution after
    // stepping by the damped step would represent a converged solution, and
    // return 0 otherwise. If no damping coefficient could be found, return -2.
    if (m < NDAMP) {
        if (s1 > 1.0) {
            return 0;
        } else {
            return 1;
        }
    } else {
        return -2;
    }
}

int Newton::solve(int loglevel)
{
    int m = 0;
    bool forceNewJac = false;
    doublereal s1=1.e30;

    m_residfunc->getState(m_x.data());

    bool frst = true;
    int nJacReeval = 0;

    while (true) {
        // Check whether the Jacobian should be re-evaluated.
        if (m_jacAge > m_jacMaxAge) {
            if (loglevel > 0) {
                writelog("\nMaximum Jacobian age reached ({})\n", m_jacMaxAge);
            }
            forceNewJac = true;
        }

        if (forceNewJac) {
            evalJacobian(&m_x[0], &m_stp[0]);
            forceNewJac = false;
        }

        // compute the undamped Newton step
        step(&m_x[0], &m_stp[0], loglevel-1);

        // increment the Jacobian age
        m_jacAge++;

        // damp the Newton step
        m = dampStep(&m_x[0], &m_stp[0], &m_x1[0], &m_stp1[0], s1, loglevel-1, frst);
        if (loglevel == 1 && m >= 0) {
            if (frst) {
                writelog("\n\n    {:>10s}    {:>10s}   {:>5s}",
                         "log10(ss)","log10(s1)","N_jac");
                writelog("\n    ------------------------------------");
            }
            doublereal ss = weightedNorm(&m_x[0], &m_stp[0]);
            writelog("\n    {:10.4f}    {:10.4f}",
                     log10(ss),log10(s1));
        }
        frst = false;

        // Successful step, but not converged yet. Take the damped step, and try
        // again.
        if (m == 0) {
            copy(m_x1.begin(), m_x1.end(), m_x.begin());
        } else if (m == 1) {
            // convergence
            // if (rdt == 0) {
            //     jac.setAge(0); // for efficient sensitivity analysis
            // }
            break;
        } else if (m < 0) {
            // If dampStep fails, first try a new Jacobian if an old one was
            // being used. If it was a new Jacobian, then return -1 to signify
            // failure.
            if (m_jacAge > 1) {
                forceNewJac = true;
                if (nJacReeval > 3) {
                    break;
                }
                nJacReeval++;
                debuglog("\nRe-evaluating Jacobian, since no damping "
                         "coefficient\ncould be found with this Jacobian.\n",
                         loglevel);
            } else {
                break;
            }
        }
    }

    if (m < 0) {
        //TODO: add get solution method? vs copy into provided vector
        //copy(m_x.begin(), m_x.end(), x1);
    }
    // if (m > 0 && m_jac->nEvals() == j0) {
    //     m = 100;
    // }
    return m;
}

} // end namespace Cantera
