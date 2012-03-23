/**
 *  @file newton_utils.cpp
 */

#include "cantera/base/ct_defs.h"
#include "cantera/oneD/Domain1D.h"

#include <cstdio>

using namespace std;

namespace Cantera
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
 *
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
}
