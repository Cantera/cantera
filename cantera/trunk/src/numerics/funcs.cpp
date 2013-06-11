/**
 *  @file funcs.cpp file containing miscellaneous
 *                numerical functions.
 */
/*
 *  Copyright 2001-2003 California Institute of Technology
 *  See file License.txt for licensing information
 */

#include "cantera/numerics/funcs.h"

#include "cantera/numerics/ctlapack.h"
#include "cantera/base/ctexceptions.h"
#include "cantera/base/stringUtils.h"
#include "cantera/numerics/polyfit.h"

#include <algorithm>

using namespace std;

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

doublereal linearInterp(doublereal x, const vector_fp& xpts,
                        const vector_fp& fpts)
{
    if (x <= xpts[0]) {
        return fpts[0];
    }
    if (x >= xpts.back()) {
        return fpts.back();
    }
    vector_fp::const_iterator loc =
        lower_bound(xpts.begin(), xpts.end(), x);
    int iloc = int(loc - xpts.begin()) - 1;
    doublereal ff = fpts[iloc] +
                    (x - xpts[iloc])*(fpts[iloc + 1]
                                      - fpts[iloc])/(xpts[iloc + 1] - xpts[iloc]);
    return ff;
}

//! Fits a polynomial function to a set of data points
/*!
 *     Given a collection of points X(I) and a set of values Y(I) which
 *     correspond to some function or measurement at each of the X(I),
 *     subroutine  DPOLFT  computes the weighted least-squares polynomial
 *     fits of all degrees up to some degree either specified by the user
 *     or determined by the routine.  The fits thus obtained are in
 *     orthogonal polynomial form.  Subroutine  DP1VLU  may then be
 *     called to evaluate the fitted polynomials and any of their
 *     derivatives at any point.  The subroutine  DPCOEF  may be used to
 *     express the polynomial fits as powers of (X-C) for any specified
 *     point C.
 *
 *   @param n   The number of data points.
 *
 *   @param x   A set of grid points on which the data is specified.
 *              The array of values of the independent variable.  These
 *              values may appear in any order and need not all be
 *              distinct. There are n of them.
 *
 *   @param y  array of corresponding function values. There are n of them
 *
 *   @param w   array of positive values to be used as weights.  If
 *              W[0] is negative,  DPOLFT  will set all the weights
 *              to 1.0, which means unweighted least squares error
 *              will be minimized.  To minimize relative error, the
 *              user should set the weights to:  W(I) = 1.0/Y(I)**2,
 *              I = 1,...,N .
 *
 *   @param maxdeg  maximum degree to be allowed for polynomial fit.
 *                  MAXDEG  may be any non-negative integer less than  N.
 *                  Note -- MAXDEG  cannot be equal to  N-1  when a
 *                  statistical test is to be used for degree selection,
 *                  i.e., when input value of  EPS  is negative.
 *
 *   @param ndeg    output degree of the fit computed.
 *
 *   @param eps     Specifies the criterion to be used in determining
 *                  the degree of fit to be computed.
 *                  (1)  If  EPS  is input negative,  DPOLFT chooses the
 *                       degree based on a statistical F test of
 *                       significance.  One of three possible
 *                       significance levels will be used:  .01, .05 or
 *                      .10.  If  EPS=-1.0 , the routine will
 *                       automatically select one of these levels based
 *                       on the number of data points and the maximum
 *                       degree to be considered.  If  EPS  is input as
 *                       -.01, -.05, or -.10, a significance level of
 *                       .01, .05, or .10, respectively, will be used.
 *                  (2)  If  EPS  is set to 0.,  DPOLFT  computes the
 *                       polynomials of degrees 0 through  MAXDEG .
 *                  (3)  If  EPS  is input positive,  EPS  is the RMS
 *                       error tolerance which must be satisfied by the
 *                       fitted polynomial.  DPOLFT  will increase the
 *                       degree of fit until this criterion is met or
 *                       until the maximum degree is reached.
 *
 *  @param r    Output vector containing the first LL+1 Taylor coefficients
 *                  where LL=ABS(ndeg).
 *                    P(X) = r[0] + r[1]*(X-C) + ... + r[ndeg] * (X-C)**ndeg
 *                   ( here C = 0.0)
 *
 * @return  returns the RMS error of the polynomial of degree ndeg .
 */
doublereal polyfit(int n, doublereal* x, doublereal* y, doublereal* w, int maxdeg, int& ndeg, doublereal eps, doublereal* r)
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
    if (ierr != 1) throw CanteraError("polyfit",
                                          "DPOLFT returned error code IERR = " + int2str(ierr) +
                                          "while attempting to fit " + int2str(n) + " data points "
                                          + "to a polynomial of degree " + int2str(maxdeg));
    ndeg = ndg;
    _DPCOEF_(&ndg, &zer, r, &awork[0]);
    return epss;
}

}
