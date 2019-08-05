// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#include "ctmatutils.h"
#include "cantera/clib/ctreactor.h"
#include "cantera/clib/ct.h"

void flowdevicemethods(int nlhs, mxArray* plhs[],
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

        n = flowdev_new2(type);
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
            iok = flowdev_del(i);
            break;
        case 2:
            m = getInt(prhs[4]);
            iok = flowdev_install(i, int(v), m);
            break;
        case 3:
            // @deprecated To be removed after Cantera 2.5.
            iok = flowdev_setMassFlowRate(i, v);
            break;
        case 4:
            iok = flowdev_setValveCoeff(i, v);
            break;
        case 5:
            // @deprecated To be removed after Cantera 2.5.
            iok = flowdev_setFunction(i, int(v));
            break;
        case 7:
            iok = flowdev_setMaster(i, int(v));
            break;
        case 8:
            iok = flowdev_setPressureFunction(i, int(v));
            break;
        case 9:
            iok = flowdev_setTimeFunction(i, int(v));
            break;
        case 10:
            iok = flowdev_setMassFlowCoeff(i, v);
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
            r = flowdev_massFlowRate(i, v);
            break;
        default:
            mexErrMsgTxt("unknown job parameter");
        }
        plhs[0] = mxCreateNumericMatrix(1,1,mxDOUBLE_CLASS,mxREAL);
        double* h = mxGetPr(plhs[0]);
        *h = r;
        if (r == Undef) {
            reportError();
        }
        return;
    }
}
