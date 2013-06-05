/**
 *  @file MultiNewton.cpp Damped Newton solver for 1D multi-domain problems
 */

/*
 *  Copyright 2001 California Institute of Technology
 */

#include "cantera/oneD/MultiNewton.h"

#include "cantera/base/ctexceptions.h"
#include "cantera/base/vec_functions.h"
#include "cantera/base/stringUtils.h"

#include <cstdio>

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
    char buf[100];
    size_t np = r.nPoints();
    size_t nv = r.nComponents();
    Indx index(nv, np);
    doublereal above, below, val, newval;
    size_t m, j;
    doublereal fbound = 1.0;
    bool wroteTitle = false;
    for (m = 0; m < nv; m++) {
        above = r.upperBound(m);
        below = r.lowerBound(m);

        for (j = 0; j < np; j++) {
            val = x[index(m,j)];
            if (loglevel > 0) {
                if (val > above + 1.0e-12 || val < below - 1.0e-12) {
                    sprintf(buf, "domain %s: %20s(%s) = %10.3e (%10.3e, %10.3e)\n",
                            int2str(r.domainIndex()).c_str(),
                            r.componentName(m).c_str(), int2str(j).c_str(),
                            val, below, above);
                    writelog(string("\nERROR: solution out of bounds.\n")+buf);
                }
            }

            newval = val + step[index(m,j)];

            if (newval > above) {
                fbound = std::max(0.0, std::min(fbound,
                                                (above - val)/(newval - val)));
            } else if (newval < below) {
                fbound = std::min(fbound, (val - below)/(val - newval));
            }

            if (loglevel > 1 && (newval > above || newval < below)) {
                if (!wroteTitle) {
                    writelog("\nNewton step takes solution out of bounds.\n\n");
                    sprintf(buf,"  %12s  %12s  %4s  %10s  %10s  %10s  %10s\n",
                            "domain","component","pt","value","step","min","max");
                    wroteTitle = true;
                    writelog(buf);
                }
                sprintf(buf, "          %4s  %12s  %4s  %10.3e  %10.3e  %10.3e  %10.3e\n",
                        int2str(r.domainIndex()).c_str(),
                        r.componentName(m).c_str(), int2str(j).c_str(),
                        val, step[index(m,j)], below, above);
                writelog(buf);
            }
        }
    }
    return fbound;
}

/**
 * This function computes the square of a weighted norm of a step
 * vector for one domain.
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
 * Here \f$\epsilon_{r,n} \f$ is the relative error tolerance for
 * component n, and multiplies the average magnitude of
 * solution component n in the domain. The second term,
 * \f$\epsilon_{a,n}\f$, is the absolute error tolerance for component
 * n.
 */
doublereal norm_square(const doublereal* x,
                       const doublereal* step, Domain1D& r)
{
    doublereal f, ewt, esum, sum = 0.0;
    size_t n, j;
    doublereal f2max = 0.0;
    size_t nv = r.nComponents();
    size_t np = r.nPoints();

    for (n = 0; n < nv; n++) {
        esum = 0.0;
        for (j = 0; j < np; j++) {
            esum += fabs(x[nv*j + n]);
        }
        ewt = r.rtol(n)*esum/np + r.atol(n);
        for (j = 0; j < np; j++) {
            f = step[nv*j + n]/ewt;
            sum += f*f;
            if (f*f > f2max) {
                f2max = f*f;
            }
        }
    }
    return sum;
}

} // end unnamed-namespace

//-----------------------------------------------------------
//                  constants
//-----------------------------------------------------------

const string dashedline =
    "-----------------------------------------------------------------";

const doublereal DampFactor = sqrt(2.0);
const size_t NDAMP = 7;

//-----------------------------------------------------------
//                 MultiNewton methods
//-----------------------------------------------------------

MultiNewton::MultiNewton(int sz)
    : m_maxAge(5)
{
    m_n  = sz;
    m_elapsed = 0.0;
}

MultiNewton::~MultiNewton()
{
    for (size_t i = 0; i < m_workarrays.size(); i++) {
        delete[] m_workarrays[i];
    }
}

void MultiNewton::resize(size_t sz)
{
    m_n = sz;
    for (size_t i = 0; i < m_workarrays.size(); i++) {
        delete[] m_workarrays[i];
    }
    m_workarrays.clear();
}

doublereal MultiNewton::norm2(const doublereal* x,
                              const doublereal* step, OneDim& r) const
{
    doublereal f, sum = 0.0;//, fmx = 0.0;
    size_t nd = r.nDomains();
    for (size_t n = 0; n < nd; n++) {
        f = norm_square(x + r.start(n), step + r.start(n),
                        r.domain(n));
        sum += f;
    }
    sum /= r.size();
    return sqrt(sum);
}

void MultiNewton::step(doublereal* x, doublereal* step,
                       OneDim& r, MultiJac& jac, int loglevel)
{
    size_t iok;
    size_t sz = r.size();
    r.eval(npos, x, step);
#undef DEBUG_STEP
#ifdef DEBUG_STEP
    vector_fp ssave(sz, 0.0);
    for (size_t n = 0; n < sz; n++) {
        step[n] = -step[n];
        ssave[n] = step[n];
    }
#else
    for (size_t n = 0; n < sz; n++) {
        step[n] = -step[n];
    }
#endif

    iok = jac.solve(step, step);

    // if iok is non-zero, then solve failed
    if (iok != 0) {
        iok--;
        size_t nd = r.nDomains();
        size_t n;
        for (n = nd-1; n != npos; n--)
            if (iok >= r.start(n)) {
                break;
            }
        Domain1D& dom = r.domain(n);
        size_t offset = iok - r.start(n);
        size_t pt = offset/dom.nComponents();
        size_t comp = offset - pt*dom.nComponents();
        throw CanteraError("MultiNewton::step",
                           "Jacobian is singular for domain "+
                           dom.id() + ", component "
                           +dom.componentName(comp)+" at point "
                           +int2str(pt)+"\n(Matrix row "
                           +int2str(iok)+") \nsee file bandmatrix.csv\n");
    } else if (int(iok) < 0)
        throw CanteraError("MultiNewton::step",
                           "iok = "+int2str(iok));

#ifdef DEBUG_STEP
    bool ok = false;
    Domain1D* d;
    if (!ok) {
        for (size_t n = 0; n < sz; n++) {
            d = r.pointDomain(n);
            int nvd = d->nComponents();
            int pt = (n - d->loc())/nvd;
            cout << "step: " << pt << "  " <<
                 r.pointDomain(n)->componentName(n - d->loc() - nvd*pt)
                 << "    " << x[n] << "     " << ssave[n] << "   " << step[n] << endl;
        }
    }
#endif
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
        writelog(dashedline);

        sprintf(m_buf,"\n%s  %9s   %9s     %9s   %9s   %9s  %5s %5s\n",
                "m","F_damp","F_bound","log10(ss)",
                "log10(s0)","log10(s1)","N_jac","Age");
        writelog(m_buf);
        writelog(dashedline+"\n");
    }

    // compute the weighted norm of the undamped step size step0
    doublereal s0 = norm2(x0, step0, r);

    // compute the multiplier to keep all components in bounds
    doublereal fbound = boundStep(x0, step0, r, loglevel-1);

    // if fbound is very small, then x0 is already close to the
    // boundary and step0 points out of the allowed domain. In
    // this case, the Newton algorithm fails, so return an error
    // condition.
    if (fbound < 1.e-10) {
        writelog("\nAt limits.\n", loglevel);
        return -3;
    }


    //--------------------------------------------
    //           Attempt damped step
    //--------------------------------------------

    // damping coefficient starts at 1.0
    doublereal damp = 1.0;

    doublereal ff;

    size_t m;
    for (m = 0; m < NDAMP; m++) {

        ff = fbound*damp;

        // step the solution by the damped step size
        for (size_t j = 0; j < m_n; j++) {
            x1[j] = ff*step0[j] + x0[j];
        }

        // compute the next undamped step that would result if x1
        // is accepted
        step(x1, step1, r, jac, loglevel-1);

        // compute the weighted norm of step1
        s1 = norm2(x1, step1, r);

        // write log information
        if (loglevel > 0) {
            doublereal ss = r.ssnorm(x1,step1);
            sprintf(m_buf,"\n%s  %9.5f   %9.5f   %9.5f   %9.5f   %9.5f %4d  %d/%d",
                    int2str(m).c_str(), damp, fbound, log10(ss+SmallNumber),
                    log10(s0+SmallNumber),
                    log10(s1+SmallNumber),
                    jac.nEvals(), jac.age(), m_maxAge);
            writelog(m_buf);
        }

        // if the norm of s1 is less than the norm of s0, then
        // accept this damping coefficient. Also accept it if this
        // step would result in a converged solution. Otherwise,
        // decrease the damping coefficient and try again.

        if (s1 < 1.0 || s1 < s0) {
            break;
        }
        damp /= DampFactor;
    }

    // If a damping coefficient was found, return 1 if the
    // solution after stepping by the damped step would represent
    // a converged solution, and return 0 otherwise. If no damping
    // coefficient could be found, return -2.
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

    doublereal* x    = getWorkArray();
    doublereal* stp  = getWorkArray();
    doublereal* stp1 = getWorkArray();

    copy(x0, x0 + m_n, x);

    bool frst = true;
    doublereal rdt = r.rdt();
    int j0 = jac.nEvals();
    int nJacReeval = 0;

    while (1 > 0) {

        // Check whether the Jacobian should be re-evaluated.
        if (jac.age() > m_maxAge) {
            writelog("\nMaximum Jacobian age reached ("+int2str(m_maxAge)+")\n", loglevel);
            forceNewJac = true;
        }

        if (forceNewJac) {
            r.eval(npos, x, stp, 0.0, 0);
            jac.eval(x, stp, 0.0);
            jac.updateTransient(rdt, DATA_PTR(r.transientMask()));
            forceNewJac = false;
        }

        // compute the undamped Newton step
        step(x, stp, r, jac, loglevel-1);

        // increment the Jacobian age
        jac.incrementAge();

        // damp the Newton step
        m = dampStep(x, stp, x1, stp1, s1, r, jac, loglevel-1, frst);
        if (loglevel == 1 && m >= 0) {
            if (frst) {
                sprintf(m_buf,"\n\n    %10s    %10s   %5s ",
                        "log10(ss)","log10(s1)","N_jac");
                writelog(m_buf);
                sprintf(m_buf,"\n    ------------------------------------");
                writelog(m_buf);
            }
            doublereal ss = r.ssnorm(x, stp);
            sprintf(m_buf,"\n    %10.4f    %10.4f       %d ",
                    log10(ss),log10(s1),jac.nEvals());
            writelog(m_buf);
        }
        frst = false;

        // Successful step, but not converged yet. Take the damped
        // step, and try again.
        if (m == 0) {
            copy(x1, x1 + m_n, x);
        }

        // convergence
        else if (m == 1) {
            break;
        }

        // If dampStep fails, first try a new Jacobian if an old
        // one was being used. If it was a new Jacobian, then
        // return -1 to signify failure.
        else if (m < 0) {
            if (jac.age() > 1) {
                forceNewJac = true;
                if (nJacReeval > 3) {
                    break;
                }
                nJacReeval++;
                writelog("\nRe-evaluating Jacobian, since no damping "
                         "coefficient\ncould be found with this Jacobian.\n",
                         loglevel);
            } else {
                break;
            }
        }
    }

    if (m < 0) {
        copy(x, x + m_n, x1);
    }
    if (m > 0 && jac.nEvals() == j0) {
        m = 100;
    }
    releaseWorkArray(x);
    releaseWorkArray(stp);
    releaseWorkArray(stp1);
    m_elapsed += (clock() - t0)/(1.0*CLOCKS_PER_SEC);
    return m;
}

doublereal* MultiNewton::getWorkArray()
{
    doublereal* w = 0;

    if (!m_workarrays.empty()) {
        w = m_workarrays.back();
        m_workarrays.pop_back();
    } else {
        w = new doublereal[m_n];
    }
    return w;
}

void MultiNewton::releaseWorkArray(doublereal* work)
{
    m_workarrays.push_back(work);
}

} // end namespace Cantera
