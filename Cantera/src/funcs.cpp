
// miscellaneous functions

#ifdef WIN32
#pragma warning(disable:4786)
#pragma warning(disable:4503)
#endif

#include <vector>
using namespace std;

#include "ct_defs.h"
#include "ctexceptions.h"

extern "C" {

    int dpolft_(integer* n, doublereal* x, doublereal* y, doublereal* w, 
        integer* maxdeg, integer* ndeg, doublereal* eps, doublereal* r,
        integer* ierr, doublereal* a);

    int dpcoef_(integer* l, doublereal* c, doublereal* tc, doublereal* a);
}

namespace Cantera {

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

        dpolft_(&nn, x, y, w, &mdeg, &ndg, &epss, coeffs.begin(), 
            &ierr, awork.begin());
        if (ierr != 1) throw CanteraError("polyfit",
            "DPOLFT returned error code IERR = " + int2str(ierr) + 
            "while attempting to fit " + int2str(n) + " data points "
            + "to a polynomial of degree " + int2str(maxdeg));
        ndeg = ndg;
        dpcoef_(&ndg, &zer, r, awork.begin());
        return epss;
    }   


}
