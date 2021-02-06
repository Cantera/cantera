/**
 *  @file wallmethods.cpp
 */

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#include "cantera/clib/ctreactor.h"
#include "cantera/clib/ct.h"
#include "ctmatutils.h"

void wallmethods(int nlhs, mxArray* plhs[],
                 int nrhs, const mxArray* prhs[])
{
    int m, iok = 0, n;
    int job = getInt(prhs[1]);
    double r = Undef;
    double v = Undef;
    if (nrhs > 3) {
        v = getDouble(prhs[3]);
    }

    // constructor
    if (job == 0) {
        char* type = getString(prhs[2]);
        n = wall_new2(type);
        plhs[0] = mxCreateNumericMatrix(1,1,mxDOUBLE_CLASS,mxREAL);
        double* h = mxGetPr(plhs[0]);
        *h = double(n);
        if (n < 0) {
            reportError();
        }
        return;
    }

    // options that do not return a value

    int i = getInt(prhs[2]);
    if (job < 20) {
        switch (job) {
        case 1:
            iok = wall_del(i);
            break;
        case 4:
            m = getInt(prhs[4]);
            iok = wall_install(i, int(v), m);
            break;
        case 5:
            iok = wall_setArea(i, v);
            break;
        case 6:
            iok = wall_setThermalResistance(i, v);
            break;
        case 7:
            iok = wall_setHeatTransferCoeff(i, v);
            break;
        case 8:
            iok = wall_setHeatFlux(i, int(v));
            break;
        case 9:
            iok = wall_setExpansionRateCoeff(i, v);
            break;
        case 10:
            iok = wall_setVelocity(i, int(v));
            break;
        case 11:
            iok = wall_ready(i);
            break;
        case 13:
            iok = wall_setEmissivity(i, v);
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
            r = wall_vdot(i, v);
            break;
        case 22:
            r = wall_Q(i, v);
            break;
        case 23:
            r = wall_area(i);
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
