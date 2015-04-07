/**
 *  @file wallmethods.cpp
 */

#include "clib/ctreactor.h"
#include "clib/ct.h"
#include "ctmatutils.h"

void wallmethods(int nlhs, mxArray* plhs[],
                 int nrhs, const mxArray* prhs[])
{
    int m, iok = 0, n;
    int job = getInt(prhs[1]);
    int i = getInt(prhs[2]);

    double r = Undef;
    double v = Undef;
    if (nrhs > 3) {
        v = getDouble(prhs[3]);
    }

    // constructor
    if (job == 0) {
        n = wall_new(i);
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
            iok = wall_del(i);
            break;
        case 2:
            iok = wall_copy(i);
            break;
        case 3:
            iok = wall_assign(i,int(v));
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
        case 12:
            n = getInt(prhs[3]);
            m = getInt(prhs[4]);
            iok = wall_setkinetics(i, n, m);
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
    }


    // options that return a value of type 'double'

    else if (job < 40) {
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

