/**
 *  @file newton_utils.cpp
 */

#ifdef WIN32
#pragma warning(disable:4786)
#pragma warning(disable:4503)
#endif

#include "../ct_defs.h"
#include "Resid1D.h"

namespace Cantera {

    class Indx {
    public:
        Indx(int nv, int np) : m_nv(nv), m_np(np) {}
        int m_nv, m_np;
        int operator()(int m, int j) { return j*m_nv + m; }
    };

        
    doublereal bound_step(const doublereal* x, const doublereal* step, 
        Resid1D& r, int loglevel=0) {

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
                    if (val > above + Tiny || val < below - Tiny) {
                        sprintf(buf, "domain %d: %20s(%d) = %10.3e (%10.3e, %10.3e)\n",
                            r.domainIndex(), r.componentName(m).c_str(), j, val, below, above);
                        writelog(string("ERROR: solution out of bounds.\n")+buf);
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


    doublereal norm_square(const doublereal* x, 
        const doublereal* step, Resid1D& r) {
        doublereal f, ewt, esum, sum = 0.0;
        int n, j;

        int nv = r.nComponents();
        int np = r.nPoints();

        for (n = 0; n < nv; n++) {
            esum = 0.0;
            for (j = 0; j < np; j++) esum += fabs(x[nv*j + n]);
            ewt = r.rtol(n)*esum/np + r.atol(n);
            for (j = 0; j < np; j++) {
                f = step[nv*j + n]/ewt;
                sum += f*f;
//                 if (fabs(f) > fmx) {
//                     fmx = fabs(f);
//                     jmx = j;
//                     nmx = n;
//                 }
            }
        }
        return sum;
    }
}
