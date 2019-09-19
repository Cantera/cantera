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
doublereal bound_step(const doublereal* x, const doublereal* step,
                      Domain1D& r, int loglevel)
{
    size_t np = r.nPoints();
    size_t nv = r.nComponents();
    Indx index(nv, np);
    doublereal fbound = 1.0;
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
                    writelog("\nNewton step takes solution out of bounds.\n\n");
                    writelog("  {:>12s}  {:>12s}  {:>4s}  {:>10s}  {:>10s}  {:>10s}  {:>10s}\n",
                             "domain","component","pt","value","step","min","max");
                    wroteTitle = true;
                }
                writelog("          {:4d}  {:>12s}  {:4d}  {:10.3e}  {:10.3e}  {:10.3e}  {:10.3e}\n",
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
 * \f[
 *    \sum_{n,j} \left(\frac{s_{n,j}}{w_n}\right)^2
 * \f]
 * where the error weight for solution component \f$n\f$ is given by
 * \f[
 *     w_n = \epsilon_{r,n} \frac{\sum_j |x_{n,j}|}{J} + \epsilon_{a,n}.
 * \f]
 * Here \f$\epsilon_{r,n} \f$ is the relative error tolerance for component n,
 * and multiplies the average magnitude of solution component n in the domain.
 * The second term, \f$\epsilon_{a,n}\f$, is the absolute error tolerance for
 * component n.
 */
doublereal norm_square(const doublereal* x,
                       const doublereal* step, Domain1D& r)
{
    double sum = 0.0;
    doublereal f2max = 0.0;
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
const doublereal DampFactor = sqrt(2.0);
const size_t NDAMP = 7;

// ---------------- MultiNewton methods ----------------

MultiNewton::MultiNewton(int sz)
    : m_maxAge(5)
{
    m_n = sz;
    m_elapsed = 0.0;
}

void MultiNewton::resize(size_t sz)
{
    m_n = sz;
    m_x.resize(m_n);
    m_stp.resize(m_n);
    m_stp1.resize(m_n);
}

doublereal MultiNewton::norm2(const doublereal* x,
                              const doublereal* step, OneDim& r) const
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

void MultiNewton::step(doublereal* x, doublereal* step,
                       OneDim& r, MultiJac& jac, int loglevel)
{
    r.eval(npos, x, step);
    for (size_t n = 0; n < r.size(); n++) {
        step[n] = -step[n];
    }

    try {
        jac.solve(step, step);
    } catch (CanteraError&) {
        int iok = jac.info() - 1;
        if (iok >= 0) {
            size_t nd = r.nDomains();
            size_t n;
            for (n = nd-1; n != npos; n--) {
                if (iok >= static_cast<int>(r.start(n))) {
                    break;
                }
            }
            Domain1D& dom = r.domain(n);
            size_t offset = iok - r.start(n);
            size_t pt = offset/dom.nComponents();
            size_t comp = offset - pt*dom.nComponents();
            throw CanteraError("MultiNewton::step",
                "Jacobian is singular for domain {}, component {} at point {}\n"
                "(Matrix row {})",
                dom.id(), dom.componentName(comp), pt, iok);
        } else {
            throw;
        }
    }
}

doublereal MultiNewton::boundStep(const doublereal* x0,
                                  const doublereal* step0, const OneDim& r, int loglevel)
{
    doublereal fbound = 1.0;
    for (size_t i = 0; i < r.nDomains(); i++) {
        fbound = std::min(fbound,
                          bound_step(x0 + r.start(i), step0 + r.start(i),
                                     r.domain(i), loglevel));
    }
    return fbound;
}

int MultiNewton::dampStep(const doublereal* x0, const doublereal* step0,
                          doublereal* x1, doublereal* step1, doublereal& s1,
                          OneDim& r, MultiJac& jac, int loglevel, bool writetitle)
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
    doublereal s0 = norm2(x0, step0, r);

    // compute the multiplier to keep all components in bounds
    doublereal fbound = boundStep(x0, step0, r, loglevel-1);

    // if fbound is very small, then x0 is already close to the boundary and
    // step0 points out of the allowed domain. In this case, the Newton
    // algorithm fails, so return an error condition.
    if (fbound < 1.e-10) {
        debuglog("\nAt limits.\n", loglevel);
        return -3;
    }

    // ---------- Attempt damped step ----------

    // damping coefficient starts at 1.0
    doublereal damp = 1.0;
    size_t m;
    for (m = 0; m < NDAMP; m++) {
        double ff = fbound*damp;

        // step the solution by the damped step size
        for (size_t j = 0; j < m_n; j++) {
            x1[j] = ff*step0[j] + x0[j];
        }

        // compute the next undamped step that would result if x1 is accepted
        step(x1, step1, r, jac, loglevel-1);

        // compute the weighted norm of step1
        s1 = norm2(x1, step1, r);

        // write log information
        if (loglevel > 0) {
            doublereal ss = r.ssnorm(x1,step1);
            writelog("\n{:d}  {:9.5f}   {:9.5f}   {:9.5f}   {:9.5f}   {:9.5f} {:4d}  {:d}/{:d}",
                     m, damp, fbound, log10(ss+SmallNumber),
                     log10(s0+SmallNumber), log10(s1+SmallNumber),
                     jac.nEvals(), jac.age(), m_maxAge);
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

int MultiNewton::solve(doublereal* x0, doublereal* x1,
                       OneDim& r, MultiJac& jac, int loglevel)
{
    clock_t t0 = clock();
    int m = 0;
    bool forceNewJac = false;
    doublereal s1=1.e30;

    copy(x0, x0 + m_n, &m_x[0]);

    bool frst = true;
    doublereal rdt = r.rdt();
    int j0 = jac.nEvals();
    int nJacReeval = 0;

    while (true) {
        // Check whether the Jacobian should be re-evaluated.
        if (jac.age() > m_maxAge) {
            if (loglevel > 0) {
                writelog("\nMaximum Jacobian age reached ({})\n", m_maxAge);
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
        m = dampStep(&m_x[0], &m_stp[0], x1, &m_stp1[0], s1, r, jac, loglevel-1, frst);
        if (loglevel == 1 && m >= 0) {
            if (frst) {
                writelog("\n\n    {:>10s}    {:>10s}   {:>5s}",
                         "log10(ss)","log10(s1)","N_jac");
                writelog("\n    ------------------------------------");
            }
            doublereal ss = r.ssnorm(&m_x[0], &m_stp[0]);
            writelog("\n    {:10.4f}    {:10.4f}       {:d}",
                     log10(ss),log10(s1),jac.nEvals());
        }
        frst = false;

        // Successful step, but not converged yet. Take the damped step, and try
        // again.
        if (m == 0) {
            copy(x1, x1 + m_n, m_x.begin());
        } else if (m == 1) {
            // convergence
            if (rdt == 0) {
                jac.setAge(0); // for efficient sensitivity analysis
            }
            break;
        } else if (m < 0) {
            // If dampStep fails, first try a new Jacobian if an old one was
            // being used. If it was a new Jacobian, then return -1 to signify
            // failure.
            if (jac.age() > 1) {
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
        copy(m_x.begin(), m_x.end(), x1);
    }
    if (m > 0 && jac.nEvals() == j0) {
        m = 100;
    }
    m_elapsed += (clock() - t0)/(1.0*CLOCKS_PER_SEC);
    return m;
}

} // end namespace Cantera
