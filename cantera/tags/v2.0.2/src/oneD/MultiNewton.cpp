/**
 *  @file MultiNewton.cpp
 *
 *  Damped Newton solver for 1D multi-domain problems
 */

/*
 *  Copyright 2001 California Institute of Technology
 */

#include <vector>
using namespace std;

#include "cantera/oneD/MultiNewton.h"

#include "cantera/base/ctexceptions.h"
#include "cantera/base/vec_functions.h"
#include "cantera/base/stringUtils.h"

#include <cstdio>
#include <cmath>
#include <ctime>

using namespace std;

namespace Cantera
{

//----------------------------------------------------------
//                  function declarations
//----------------------------------------------------------

// declarations for functions in newton_utils.h
doublereal bound_step(const doublereal* x,
                      const doublereal* step, Domain1D& r, int loglevel=0);
doublereal norm_square(const doublereal* x,
                       const doublereal* step, Domain1D& r);



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

/**
 * Prepare for a new solution vector length.
 */
void MultiNewton::resize(size_t sz)
{
    m_n = sz;
    for (size_t i = 0; i < m_workarrays.size(); i++) {
        delete[] m_workarrays[i];
    }
    m_workarrays.clear();
}


/**
 * Compute the weighted 2-norm of 'step'.
 */
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


/**
 * Compute the undamped Newton step.  The residual function is
 * evaluated at x, but the Jacobian is not recomputed.
 */
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


/**
 * Return the factor by which the undamped Newton step 'step0'
 * must be multiplied in order to keep all solution components in
 * all domains between their specified lower and upper bounds.
 */
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


/**
 * On entry, step0 must contain an undamped Newton step for the
 * solution x0. This method attempts to find a damping coefficient
 * such that the next undamped step would have a norm smaller than
 * that of step0. If successful, the new solution after taking the
 * damped step is returned in x1, and the undamped step at x1 is
 * returned in step1.
 */
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
        if (loglevel > 0) {
            writelog("\nAt limits.\n");
        }
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


/**
 * Find the solution to F(X) = 0 by damped Newton iteration.  On
 * entry, x0 contains an initial estimate of the solution.  On
 * successful return, x1 contains the converged solution.
 */
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
            if (loglevel > 0) {
                writelog("\nMaximum Jacobian age reached ("+int2str(m_maxAge)+")\n");
            }
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
            goto done;
        }

        // If dampStep fails, first try a new Jacobian if an old
        // one was being used. If it was a new Jacobian, then
        // return -1 to signify failure.
        else if (m < 0) {
            if (jac.age() > 1) {
                forceNewJac = true;
                if (nJacReeval > 3) {
                    goto done;
                }
                nJacReeval++;
                if (loglevel > 0)
                    writelog("\nRe-evaluating Jacobian, since no damping "
                             "coefficient\ncould be found with this Jacobian.\n");
            } else {
                goto done;
            }
        }
    }

done:
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


/**
 * Get a pointer to an array of length m_n for temporary work
 * space.
 */
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

/**
 * Release a work array by pushing its pointer onto the stack of
 * available arrays.
 */
void MultiNewton::releaseWorkArray(doublereal* work)
{
    m_workarrays.push_back(work);
}
}


// $Log: Newton.cpp,v
