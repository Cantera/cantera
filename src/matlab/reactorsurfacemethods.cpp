/**
 *  @file reactorsurfacemethods.cpp
 */

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#include "cantera/clib/ctreactor.h"
#include "cantera/clib/ct.h"
#include "ctmatutils.h"

void reactorsurfacemethods(int nlhs, mxArray* plhs[],
                           int nrhs, const mxArray* prhs[])
{
    int iok = 0, n;
    int job = getInt(prhs[1]);
    int i = getInt(prhs[2]);
    double r = Undef;
    double v = Undef;
    if (nrhs > 3) {
        v = getDouble(prhs[3]);
    }

    // constructor
    if (job == 0) {
        n = reactorsurface_new(i);
        plhs[0] = mxCreateNumericMatrix(1,1,mxDOUBLE_CLASS,mxREAL);
        double* h = mxGetPr(plhs[0]);
        *h = double(n);
        if (n < 0) {
            reportError();
        }
        return;
    }

    // options that do not return a value
    if (job < 20) {
        switch (job) {
        case 1:
            iok = reactorsurface_del(i);
            break;
        case 4:
            n = getInt(prhs[3]);
            iok = reactorsurface_install(i, n);
            break;
        case 5:
            iok = reactorsurface_setArea(i, v);
            break;
        case 12:
            n = getInt(prhs[3]);
            iok = reactorsurface_setkinetics(i, n);
            break;
        default:
            mexErrMsgTxt("unknown job parameter");
        }
        plhs[0] = mxCreateNumericMatrix(1,1,mxDOUBLE_CLASS,mxREAL);
        double* h = mxGetPr(plhs[0]);
        *h = double(iok);
        if (iok < 0) {
            reportError();
        }
        return;
    } else if (job < 40) {
        // options that return a value of type 'double'
        switch (job) {
        case 23:
            r = reactorsurface_area(i);
            break;
        default:
            mexErrMsgTxt("unknown job parameter");
        }
        plhs[0] = mxCreateNumericMatrix(1,1,mxDOUBLE_CLASS,mxREAL);
        double* h = mxGetPr(plhs[0]);
        *h = r;
        if (r == DERR) {
            reportError();
        }
        return;
    }
}
