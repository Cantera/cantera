/**
 * @file reactornetmethods.cpp
 */

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#include "cantera/clib/ctreactor.h"
#include "cantera/clib/ct.h"
#include "ctmatutils.h"

void reactornetmethods(int nlhs, mxArray* plhs[],
                       int nrhs, const mxArray* prhs[])
{
    int iok = 0, n;
    int job = getInt(prhs[1]);
    int i = getInt(prhs[2]);
    double r = Undef;
    double v = Undef;
    double v2 = -1.0;
    if (nrhs > 3) {
        v = getDouble(prhs[3]);
    }
    if (nrhs > 4) {
        v2 = getDouble(prhs[4]);
    }

    // constructor
    if (job == 0) {
        n = reactornet_new();
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
            iok = reactornet_del(i);
            break;
        case 4:
            iok = reactornet_addreactor(i, int(v));
            break;
        case 5:
            iok = reactornet_setInitialTime(i, v);
            break;
        case 6:
            iok = reactornet_setMaxTimeStep(i, v);
            break;
        case 7:
            iok = reactornet_setTolerances(i, v, v2);
            break;
        case 8:
            iok = reactornet_advance(i, v);
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
        case 21:
            r = reactornet_step(i);
            break;
        case 22:
            r = reactornet_time(i);
            break;
        case 23:
            r = reactornet_rtol(i);
            break;
        case 24:
            r = reactornet_atol(i);
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
