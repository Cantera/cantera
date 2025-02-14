//! @file MultiNewton.cpp Damped Newton solver for 1D multi-domain problems

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#include "cantera/oneD/MultiNewton.h"
#include "cantera/base/utilities.h"

#include <ctime>

using namespace std;

namespace Cantera
{

MultiNewton::MultiNewton(int sz)
    : m_n(sz)
{
}

void MultiNewton::resize(size_t sz)
{
    m_n = sz;
    m_x.resize(m_n);
    m_stp.resize(m_n);
    m_stp1.resize(m_n);
}

void MultiNewton::step(double* x, double* step, SteadyStateSystem& r, int loglevel)
{
    r.eval(npos, x, step);
    for (size_t n = 0; n < r.size(); n++) {
        step[n] = -step[n];
    }

    auto jac = r.linearSolver();
    try {
        jac->solve(r.size(), step, step);
    } catch (CanteraError&) {
        if (jac->info() > 0) {
            // Positive value for "info" indicates the row where factorization failed
            size_t row = static_cast<size_t>(jac->info() - 1);
            throw CanteraError("MultiNewton::step",
                "Jacobian is singular for matrix row {}:\n{}",
                row, r.componentName(row));
        }
        throw;
    }
}

double MultiNewton::boundStep(const double* x0, const double* step0, const SteadyStateSystem& r,
                              int loglevel)
{
    const static string separator = fmt::format("\n     {:=>71}", ""); // equals sign separator
    double fbound = 1.0;
    bool wroteTitle = false;

    for (size_t i = 0; i < size(); i++) {
        double above = r.upperBound(i);
        double below = r.lowerBound(i);

        double val = x0[i];
        if (loglevel > 0 && (val > above + 1.0e-12 || val < below - 1.0e-12)) {
            writelog("\nERROR: solution component {} out of bounds.\n", i);
            writelog("{}: value = {:10.3e} (lower = {:10.3e}, upper = {:10.3e})\n",
                     r.componentName(i), val, below, above);
        }

        double newval = val + step0[i];

        if (newval > above) {
            fbound = std::max(0.0, std::min(fbound,
                                            (above - val)/(newval - val)));
        } else if (newval < below) {
            fbound = std::min(fbound, (val - below)/(val - newval));
        }

        if (loglevel > 1 && (newval > above || newval < below)) {
            if (!wroteTitle){
                string header = fmt::format("     {:=>10}", "") +
                                " Undamped Newton step takes solution out of bounds " +
                                fmt::format("{:=>10}", "");
                writelog("\n{}", header);
                writelog(separator);
                const auto& [custom1, custom2] = r.componentTableHeader();
                // Split header across 2 lines to shorten the line length
                writelog("\n     {:<24s}    {:<10s}  {:<10s}  {:<9s}  {:<9s}",
                         custom1, "", "Value", "Min", "Max");
                writelog("\n     {:<24s}    {:<10s}  {:<10s}  {:<9s}  {:<9s}",
                         custom2, "Value", "Change", "Bound", "Bound");
                writelog(separator);
                wroteTitle = true;
            }
            string comp_info = r.componentTableLabel(i);
            writelog("\n     {:<24s}  {:>10.3e}  {:>10.3e}  {:>9.2e}  {:>9.2e}",
                     comp_info, val, step0[i], below, above);
        }
    }
    if (loglevel > 1 && wroteTitle) { // If a title was written, close up the table
        writelog(separator);
    }
return fbound;
}

int MultiNewton::dampStep(const double* x0, const double* step0,
                          double* x1, double* step1, double& s1,
                          SteadyStateSystem& r, int loglevel, bool writetitle)
{
    // write header
    if (loglevel > 0 && writetitle) {
        writelog("\n\n  {:-^70}", " Damped Newton iteration ");
        writelog("\n  {:<4s}  {:<10s}   {:<10s}  {:<7s}  {:<7s}  {:<7s}  {:<5s}  {:<3s}\n",
                 "Iter", "F_damp", "F_bound", "log(ss)",
                 "log(s0)", "log(s1)", "N_jac", "Age");
        writelog("  {:->70}", "");
    }

    // compute the weighted norm of the undamped step size step0
    double s0 = r.norm2(step0);

    // compute the multiplier to keep all components in bounds
    double fbound = boundStep(x0, step0, r, loglevel-1);

    // if fbound is very small, then x0 is already close to the boundary and
    // step0 points out of the allowed domain. In this case, the Newton
    // algorithm fails, so return an error condition.
    if (fbound < 1.e-10) {
        debuglog("\n  No damped step can be taken without violating solution component bounds.", loglevel);
        return -3;
    }

    // ---------- Attempt damped step ----------

    // damping coefficient starts at 1.0, but must be scaled by the
    // fbound factor to ensure that the solution remains within bounds.
    double alpha = fbound*1.0;
    size_t m;
    auto jac = r.linearSolver();
    for (m = 0; m < m_maxDampIter; m++) {
        // step the solution by the damped step size
        // x_{k+1} = x_k + alpha_k*J(x_k)^-1 F(x_k)
        for (size_t j = 0; j < m_n; j++) {
            x1[j] = x0[j] + alpha*step0[j];
        }

        // compute the next undamped step that would result if x1 is accepted
        // J(x_k)^-1 F(x_k+1)
        step(x1, step1, r, loglevel-1);

        // compute the weighted norm of step1
        s1 = r.norm2(step1);

        if (loglevel > 0) {
            double ss = r.ssnorm(x1,step1);
            writelog("\n  {:<4d}  {:<9.3e}   {:<9.3e}   {:>6.3f}   {:>6.3f}   {:>6.3f}    {:<5d}  {:d}/{:d}",
                     m, alpha, fbound, log10(ss+SmallNumber),
                     log10(s0+SmallNumber), log10(s1+SmallNumber),
                     jac->nEvals(), jac->age(), m_maxAge);
        }

        // If the norm of s1 is less than the norm of s0, then accept this
        // damping coefficient. Also accept it if this step would result in a
        // converged solution. Otherwise, decrease the damping coefficient and
        // try again.
        if (s1 < 1.0 || s1 < s0) {
            break;
        }
        alpha /= m_dampFactor;
    }

    // If a damping coefficient was found, return 1 if the solution after
    // stepping by the damped step would represent a converged solution, and
    // return 0 otherwise. If no damping coefficient could be found, return -2.
    if (m < m_maxDampIter) {
        if (s1 > 1.0) {
                debuglog("\n  Damping coefficient found (solution has not converged yet)", loglevel);
            return 0;
        } else {
                debuglog("\n  Damping coefficient found (solution has converged)", loglevel);
            return 1;
        }
    } else {
            debuglog("\n  No damping coefficient found (max damping iterations reached)", loglevel);
        return -2;
    }
}

int MultiNewton::solve(double* x0, double* x1, SteadyStateSystem& r, int loglevel)
{
    clock_t t0 = clock();
    int status = 0;
    bool forceNewJac = false;
    bool write_header = true;
    double s1=1.e30;

    copy(x0, x0 + m_n, &m_x[0]);

    double rdt = r.rdt();
    int nJacReeval = 0;
    auto jac = r.linearSolver();
    while (true) {
        // Check whether the Jacobian should be re-evaluated.
        if (jac->age() > m_maxAge) {
            if (loglevel > 1) {
                writelog("\n  Maximum Jacobian age reached ({}), updating it.", m_maxAge);
            }
            forceNewJac = true;
        }

        if (forceNewJac) {
            r.evalJacobian(&m_x[0]);
            jac->updateTransient(rdt, r.transientMask().data());
            forceNewJac = false;
        }

        // compute the undamped Newton step
        step(&m_x[0], &m_stp[0], r, loglevel-1);

        // increment the Jacobian age
        jac->incrementAge();

        // damp the Newton step
        status = dampStep(&m_x[0], &m_stp[0], x1, &m_stp1[0], s1, r, loglevel-1, write_header);
        write_header = false;

        // Successful step, but not converged yet. Take the damped step, and try
        // again.
        if (status == 0) {
            copy(x1, x1 + m_n, m_x.begin());
        } else if (status == 1) { // convergence
            if (rdt == 0) {
                jac->setAge(0); // for efficient sensitivity analysis
            }
            break;
        } else if (status < 0) {
            // If dampStep fails, first try a new Jacobian if an old one was
            // being used. If it was a new Jacobian, then return -1 to signify
            // failure.
            if (jac->age() > 1) {
                forceNewJac = true;
                if (nJacReeval > 3) {
                    break;
                }
                nJacReeval++;
                if (loglevel > 1) {
                    writelog("\n  Re-evaluating Jacobian (damping coefficient not found"
                            " with this Jacobian)");
                }
            } else {
                break;
            }
        }
    }
    // Close off the damped iteration table that is written by the dampedStep() method
    if (loglevel > 1) {
        writelog("\n  {:->70}", "");
    }

    if (status < 0) { // Reset x1 to x0 if the solution failed
        copy(m_x.begin(), m_x.end(), x1);
    }
    m_elapsed += (clock() - t0)/(1.0*CLOCKS_PER_SEC);
    return status;
}

} // end namespace Cantera
