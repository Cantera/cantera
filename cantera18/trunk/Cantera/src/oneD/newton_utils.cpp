/**
 *  @file newton_utils.cpp
 */

#ifdef WIN32
#pragma warning(disable:4786)
#pragma warning(disable:4503)
#endif

#include "ct_defs.h"
#include "Domain1D.h"

using namespace std;

namespace Cantera {

    class Indx {
    public:
        Indx(int nv, int np) : m_nv(nv), m_np(np) {}
        int m_nv, m_np;
        int operator()(int m, int j) { return j*m_nv + m; }
    };

        
    /**
     * Return a damping coefficient that keeps the solution after taking one 
     * Newton step between specified lower and upper bounds. This function only
     * considers one domain.
     */ 
    doublereal bound_step(const doublereal* x, const doublereal* step, 
        Domain1D& r, int loglevel) {

        char buf[100];
        int np = r.nPoints();
        int nv = r.nComponents();
        Indx index(nv, np);
        doublereal above, below, val, newval;
        int m, j;
        doublereal fbound = 1.0;
        bool wroteTitle = false;
        for (m = 0; m < nv; m++) {
            above = r.upperBound(m);
            below = r.lowerBound(m);

            for (j = 0; j < np; j++) {
                val = x[index(m,j)];
                if (loglevel > 0) {
                    if (val > above + 1.0e-12 || val < below - 1.0e-12) {
                        sprintf(buf, "domain %d: %20s(%d) = %10.3e (%10.3e, %10.3e)\n",
                            r.domainIndex(), r.componentName(m).c_str(), j, val, below, above);
                        writelog(string("\nERROR: solution out of bounds.\n")+buf);
                    }
                }

                newval = val + step[index(m,j)];

                if (newval > above) {
                    fbound = fmaxx( 0.0, fminn( fbound, 
                                       (above - val)/(newval - val)));
                }
                else if (newval < below) {
                    fbound = fminn(fbound, (val - below)/(val - newval));                
                }

                if (loglevel > 1 && (newval > above || newval < below)) {
                    if (!wroteTitle) { 
                        writelog("\nNewton step takes solution out of bounds.\n\n");
                        sprintf(buf,"  %12s  %12s  %4s  %10s  %10s  %10s  %10s\n",
                            "domain","component","pt","value","step","min","max");
                        wroteTitle = true;
                        writelog(buf);
                    }
                    sprintf(buf, "          %4i  %12s  %4i  %10.3e  %10.3e  %10.3e  %10.3e\n",
                        r.domainIndex(), r.componentName(m).c_str(), j, val, 
                        step[index(m,j)], below, above);
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
        const doublereal* step, Domain1D& r) {
        doublereal f, ewt, esum, sum = 0.0;
        int n, j;
        doublereal f2max = 0.0;
        int nmax = 0;
        int jmax = 0;
        int nv = r.nComponents();
        int np = r.nPoints();

        for (n = 0; n < nv; n++) {
            esum = 0.0;
            for (j = 0; j < np; j++) esum += fabs(x[nv*j + n]);
            ewt = r.rtol(n)*esum/np + r.atol(n);
            for (j = 0; j < np; j++) {
                f = step[nv*j + n]/ewt;
                sum += f*f;
                if (f*f > f2max) {
                    jmax = j;
                    nmax = n;
                    f2max = f*f;
                }
            }
        }
#undef DEBUG_NORM
#ifdef DEBUG_NORM
        cout << "max step in domain " << r.id() << ": " << f2max << endl << 
            " for component " << r.componentName(nmax) << "  at point " << jmax << endl;
#endif
        return sum;
    }
}
