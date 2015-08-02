//! @file polyfit.cpp

#include "cantera/numerics/polyfit.h"
#include "cantera/base/ctexceptions.h"
#include "cantera/base/stringUtils.h"

#ifndef FTN_TRAILING_UNDERSCORE
#define _DPOLFT_ dpolft
#define _DPCOEF_ dpcoef
#else
#define _DPOLFT_ dpolft_
#define _DPCOEF_ dpcoef_
#endif

extern "C" {
    int _DPOLFT_(integer* n, doublereal* x, doublereal* y, doublereal* w,
                 integer* maxdeg, integer* ndeg, doublereal* eps, doublereal* r,
                 integer* ierr, doublereal* a);

    int _DPCOEF_(integer* l, doublereal* c, doublereal* tc, doublereal* a);
}

namespace Cantera
{

doublereal polyfit(int n, doublereal* x, doublereal* y, doublereal* w,
                   int maxdeg, int& ndeg, doublereal eps, doublereal* r)
{
    integer nn = n;
    integer mdeg = maxdeg;
    integer ndg = ndeg;
    doublereal epss = eps;
    integer ierr;
    int worksize = 3*n + 3*maxdeg + 3;
    vector_fp awork(worksize,0.0);
    vector_fp coeffs(n+1, 0.0);
    doublereal zer = 0.0;

    _DPOLFT_(&nn, x, y, w, &mdeg, &ndg, &epss, &coeffs[0],
             &ierr, &awork[0]);
    if (ierr != 1) {
        throw CanteraError("polyfit",
                           "DPOLFT returned error code IERR = " + int2str(ierr) +
                           "while attempting to fit " + int2str(n) + " data points "
                           + "to a polynomial of degree " + int2str(maxdeg));
    }
    ndeg = ndg;
    _DPCOEF_(&ndg, &zer, r, &awork[0]);
    return epss;
}

}
