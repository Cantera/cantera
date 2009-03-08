
// miscellaneous functions

#ifdef WIN32
#pragma warning(disable:4786)
#pragma warning(disable:4503)
#endif

#include <vector>
#include <algorithm>

using namespace std;

#include "ct_defs.h"
#include "ctexceptions.h"
#include "stringUtils.h"

extern "C" {

    int dpolft_(integer* n, doublereal* x, doublereal* y, doublereal* w, 
        integer* maxdeg, integer* ndeg, doublereal* eps, doublereal* r,
        integer* ierr, doublereal* a);

    int dpcoef_(integer* l, doublereal* c, doublereal* tc, doublereal* a);
}

namespace Cantera {

    /**
     * Linearly interpolate a function defined on a discrete grid.
     * vector xpts contains a monotonic sequence of grid points, and 
     * vector fpts contains function values defined at these points.
     * The value returned is the linear interpolate at point x.
     * If x is outside the range of xpts, the value of fpts at the 
     * nearest end is returned.
     */

    doublereal linearInterp(doublereal x, const vector_fp& xpts, 
        const vector_fp& fpts) {
        if (x <= xpts[0]) 
            return fpts[0];
        if (x >= xpts.back()) 
            return fpts.back();
        vector_fp::const_iterator loc = 
            lower_bound(xpts.begin(), xpts.end(), x);
        int iloc = int(loc - xpts.begin()) - 1;
        doublereal ff = fpts[iloc] + 
            (x - xpts[iloc])*(fpts[iloc + 1] 
                - fpts[iloc])/(xpts[iloc + 1] - xpts[iloc]);
        return ff;
    }



    doublereal polyfit(int n, doublereal* x, doublereal* y, doublereal* w, 
        int maxdeg, int& ndeg, doublereal eps, doublereal* r) {
        integer nn = n;
        integer mdeg = maxdeg;
        integer ndg = ndeg;
        doublereal epss = eps;
        integer ierr;
        int worksize = 3*n + 3*maxdeg + 3;
        vector_fp awork(worksize,0.0);
        vector_fp coeffs(n+1, 0.0);
        doublereal zer = 0.0;

        dpolft_(&nn, x, y, w, &mdeg, &ndg, &epss, &coeffs[0],
            &ierr, &awork[0]);
        if (ierr != 1) throw CanteraError("polyfit",
            "DPOLFT returned error code IERR = " + int2str(ierr) + 
            "while attempting to fit " + int2str(n) + " data points "
            + "to a polynomial of degree " + int2str(maxdeg));
        ndeg = ndg;
        dpcoef_(&ndg, &zer, r, &awork[0]);
        return epss;
    }   


}
