//! @file MultiNewton.cpp Damped Newton solver for 1D multi-domain problems

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#include "cantera/oneD/MultiNewton.h"
#include "cantera/base/utilities.h"

#include <ctime>

using namespace std;

namespace Cantera
{

// unnamed-namespace for local helpers
namespace
{

class Indx
{
public:
    Indx(size_t nv, size_t np) : m_nv(nv), m_np(np) {}
    size_t m_nv, m_np;
    size_t operator()(size_t m, size_t j) {
        return j*m_nv + m;
    }
};

/**
 * Return a damping coefficient that keeps the solution after taking one
 * Newton step between specified lower and upper bounds. This function only
 * considers one domain.
 */
double bound_step(const double* x, const double* step, Domain1D& r, int loglevel)
{
    size_t np = r.nPoints();
    size_t nv = r.nComponents();
    Indx index(nv, np);
    double fbound = 1.0;
    bool wroteTitle = false;
    for (size_t m = 0; m < nv; m++) {
        double above = r.upperBound(m);
        double below = r.lowerBound(m);

        for (size_t j = 0; j < np; j++) {
            double val = x[index(m,j)];
            if (loglevel > 0 && (val > above + 1.0e-12 || val < below - 1.0e-12)) {
                writelog("\nERROR: solution out of bounds.\n");
                writelog("domain {:d}: {:>20s}({:d}) = {:10.3e} ({:10.3e}, {:10.3e})\n",
                         r.domainIndex(), r.componentName(m), j, val, below, above);
            }

            double newval = val + step[index(m,j)];

            if (newval > above) {
                fbound = std::max(0.0, std::min(fbound,
                                                (above - val)/(newval - val)));
            } else if (newval < below) {
                fbound = std::min(fbound, (val - below)/(val - newval));
            }

            if (loglevel > 1 && (newval > above || newval < below)) {
                if (!wroteTitle) {
                    writelog("\n   Undamped Newton step takes solution out of bounds:\n");
                    writelog("   {:<10s}  {:<12s}  {:<4s}  {:<10s}   {:<10s}   {:<10s}   {:<10s}\n",
                             "domain","component","pt","value","step","min","max");
                    wroteTitle = true;
                }
                writelog("   {:<10d}  {:<12s}  {:<4d}  {:<10.3e}   {:<10.3e}   {:<10.3e}   {:<10.3e}\n",
                         r.domainIndex(), r.componentName(m), j,
                         val, step[index(m,j)], below, above);
            }
        }
    }
    return fbound;
}

/**
 * This function computes the square of a weighted norm of a step vector for one
 * domain.
 *
 * @param x     Solution vector for this domain.
 * @param step  Newton step vector for this domain.
 * @param r     Object representing the domain. Used to get tolerances,
 *              number of components, and number of points.
 *
 * The return value is
 * @f[
 *    \sum_{n,j} \left(\frac{s_{n,j}}{w_n}\right)^2
 * @f]
 * where the error weight for solution component @f$ n @f$ is given by
 * @f[
 *     w_n = \epsilon_{r,n} \frac{\sum_j |x_{n,j}|}{J} + \epsilon_{a,n}.
 * @f]
 * Here @f$ \epsilon_{r,n} @f$ is the relative error tolerance for component n,
 * and multiplies the average magnitude of solution component n in the domain.
 * The second term, @f$ \epsilon_{a,n} @f$, is the absolute error tolerance for
 * component n.
 */
double norm_square(const double* x, const double* step, Domain1D& r)
{
    double sum = 0.0;
    double f2max = 0.0;
    size_t nv = r.nComponents();
    size_t np = r.nPoints();

    for (size_t n = 0; n < nv; n++) {
        double esum = 0.0;
        for (size_t j = 0; j < np; j++) {
            esum += fabs(x[nv*j + n]);
        }
        double ewt = r.rtol(n)*esum/np + r.atol(n);
        for (size_t j = 0; j < np; j++) {
            double f = step[nv*j + n]/ewt;
            sum += f*f;
            f2max = std::max(f*f, f2max);
        }
    }
    return sum;
}

} // end unnamed-namespace


// constants
const double DampFactor = sqrt(2.0);
const size_t NDAMP = 7;

// ---------------- MultiNewton methods ----------------

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

double MultiNewton::norm2(const double* x, const double* step, OneDim& r) const
{
    double sum = 0.0;
    size_t nd = r.nDomains();
    for (size_t n = 0; n < nd; n++) {
        double f = norm_square(x + r.start(n), step + r.start(n), r.domain(n));
        sum += f;
    }
    sum /= r.size();
    return sqrt(sum);
}

void MultiNewton::step(double* x, double* step, OneDim& r, MultiJac& jac, int loglevel)
{
    r.eval(npos, x, step);
    for (size_t n = 0; n < r.size(); n++) {
        step[n] = -step[n];
    }

    try {
        jac.solve(step, step);
    } catch (CanteraError&) {
        if (jac.info() > 0) {
            // Positive value for "info" indicates the row where factorization failed
            size_t row = static_cast<size_t>(jac.info() - 1);
            // Find the domain, grid point, and solution component corresponding
            // to this row
            for (size_t n = 0; n < r.nDomains(); n++) {
                Domain1D& dom = r.domain(n);
                size_t nComp = dom.nComponents();
                if (row >= dom.loc() && row < dom.loc() + nComp * dom.nPoints()) {
                    size_t offset = row - dom.loc();
                    size_t pt = offset / nComp;
                    size_t comp = offset - pt * nComp;
                    throw CanteraError("MultiNewton::step",
                        "Jacobian is singular for domain {}, component {} at point {}\n"
                        "(Matrix row {})",
                        dom.id(), dom.componentName(comp), pt, row);
                }
            }
        }
        throw;
    }
}

double MultiNewton::boundStep(const double* x0, const double* step0, const OneDim& r,
                              int loglevel)
{
    double fbound = 1.0;
    for (size_t i = 0; i < r.nDomains(); i++) {
        fbound = std::min(fbound,
                          bound_step(x0 + r.start(i), step0 + r.start(i),
                                     r.domain(i), loglevel));
    }
    return fbound;
}

int MultiNewton::dampStep(const double* x0, const double* step0,
                          double* x1, double* step1, double& s1,
                          OneDim& r, MultiJac& jac, int loglevel, bool writetitle)
{
    // write header
    if (loglevel > 0) {
        writelog("\n\n{}\n", "Damped Newton iteration:");
        writeline('-', 83, false);

        writelog("\n{:<6s}  {:<9s}    {:<9s}      {:<9s}   {:<9s}   {:<9s}  {:<5s}  {:<5s}\n",
                 "Iter", "F_damp", "F_bound", "log10(ss)",
                "log10(s0)", "log10(s1)", "N_jac", "Age");
        writeline('-', 83, false);
    }

    // compute the weighted norm of the undamped step size step0
    double s0 = norm2(x0, step0, r);

    // compute the multiplier to keep all components in bounds
    double fbound = boundStep(x0, step0, r, loglevel-1);

    // if fbound is very small, then x0 is already close to the boundary and
    // step0 points out of the allowed domain. In this case, the Newton
    // algorithm fails, so return an error condition.
    if (fbound < 1.e-10) {
        debuglog("\nNo damping step can be taken without leaving bounds on solution.\n", loglevel);
        return -3;
    }

    // ---------- Attempt damped step ----------

    // damping coefficient starts at 1.0, but must be scaled by the
    // fbound factor to ensure that the solution remains within bounds.
    double alpha = fbound*1.0;
    size_t m;
    for (m = 0; m < NDAMP; m++) {
        // step the solution by the damped step size
        // x_{k+1} = x_k + alpha_k*J(x_k)^-1 F(x_k)
        for (size_t j = 0; j < m_n; j++) {
            x1[j] = x0[j] + alpha*step0[j];
        }

        // compute the next undamped step that would result if x1 is accepted
        // J(x_k)^-1 F(x_k+1)
        step(x1, step1, r, jac, loglevel-1);

        // compute the weighted norm of step1
        s1 = norm2(x1, step1, r);

        if (loglevel > 0) {
            double ss = r.ssnorm(x1,step1);
            writelog("\n{:<6d}  {:<9.5g}   {:<9.5g}     {:<9.5f}   {:<9.5f}   {:<9.5f}  {:<5d}  {:d}/{:d}",
                     m, alpha, fbound, log10(ss+SmallNumber),
                     log10(s0+SmallNumber), log10(s1+SmallNumber),
                     jac.nEvals(), jac.age(), m_maxAge);
        }

        // If the norm of s1 is less than the norm of s0, then accept this
        // damping coefficient. Also accept it if this step would result in a
        // converged solution. Otherwise, decrease the damping coefficient and
        // try again.
        if (s1 < 1.0 || s1 < s0) {
            break;
        }
        alpha /= DampFactor;
    }

    // If a damping coefficient was found, return 1 if the solution after
    // stepping by the damped step would represent a converged solution, and
    // return 0 otherwise. If no damping coefficient could be found, return -2.
    if (m < NDAMP) {
        if (s1 > 1.0) {
            if (loglevel > 0) {
                writelog("\nDamping coefficient found, but solution does not converge.\n");
            }
            return 0;
        } else {
            if (loglevel > 0) {
                writelog("\nDamping coefficient found and solution converges.\n");
            }
            return 1;
        }
    } else {
        if (loglevel > 0) {
            writelog("\nNo damping coefficient found.\n");
        }
        return -2;
    }
}

int MultiNewton::solve(double* x0, double* x1, OneDim& r, MultiJac& jac, int loglevel)
{
    clock_t t0 = clock();
    int status = 0;
    bool forceNewJac = false;
    double s1=1.e30;

    copy(x0, x0 + m_n, &m_x[0]);

    double rdt = r.rdt();
    int nJacReeval = 0;
    while (true) {
        // Check whether the Jacobian should be re-evaluated.
        if (jac.age() > m_maxAge) {
            if (loglevel > 0) {
                writelog("\nMaximum Jacobian age reached ({}), updating it.", m_maxAge);
            }
            forceNewJac = true;
        }

        if (forceNewJac) {
            r.eval(npos, &m_x[0], &m_stp[0], 0.0, 0);
            jac.eval(&m_x[0], &m_stp[0], 0.0);
            jac.updateTransient(rdt, r.transientMask().data());
            forceNewJac = false;
        }

        // compute the undamped Newton step
        step(&m_x[0], &m_stp[0], r, jac, loglevel-1);

        // increment the Jacobian age
        jac.incrementAge();

        // damp the Newton step
        status = dampStep(&m_x[0], &m_stp[0], x1, &m_stp1[0], s1, r, jac, loglevel-1, true);

        // Successful step, but not converged yet. Take the damped step, and try
        // again.
        if (status == 0) {
            copy(x1, x1 + m_n, m_x.begin());
        } else if (status == 1) {
            // convergence
            if (rdt == 0) {
                jac.setAge(0); // for efficient sensitivity analysis
            }
            break;
        } else if (status < 0) {
            // If dampStep fails, first try a new Jacobian if an old one was
            // being used. If it was a new Jacobian, then return -1 to signify
            // failure.
            if (jac.age() > 1) {
                forceNewJac = true;
                if (nJacReeval > 3) {
                    break;
                }
                nJacReeval++;
                if (loglevel > 0) {
                    writelog("\nRe-evaluating Jacobian(damping coefficient not found"
                            " with this Jacobian)");
                }
            } else {
                break;
            }
        }
    }

    if (status < 0) { // Why is this done?
        copy(m_x.begin(), m_x.end(), x1);
    }
    m_elapsed += (clock() - t0)/(1.0*CLOCKS_PER_SEC);
    return status;
}

} // end namespace Cantera
