// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#include <fstream>

#include "cantera/clib/ct.h"
#include "ctmatutils.h"

void reportError();

void transportmethods(int nlhs, mxArray* plhs[],
                      int nrhs, const mxArray* prhs[])
{
    double vv = 0.0;
    int n = getInt(prhs[1]);
    int job = getInt(prhs[2]);
    double* h;
    int iok = 0;
    int nsp;

    if (job == -1) {
        char* model = getString(prhs[3]);
        int loglevel = getInt(prhs[4]);
        int m = -2;
        m = (int) trans_new(model, n, loglevel);
        if (m < 0) {
            reportError();
        }

        // Create matrix for the return argument.
        plhs[0] = mxCreateDoubleMatrix(1,1, mxREAL);
        double* x = mxGetPr(plhs[0]);
        *x = m;
        return;
    } else if (job == -2) {
        int loglevel = getInt(prhs[3]);
        int m = -2;
        m = (int) trans_newDefault(n, loglevel);
        if (m < 0) {
            reportError();
        }

        // Create matrix for the return argument.
        plhs[0] = mxCreateDoubleMatrix(1,1, mxREAL);
        double* x = mxGetPr(plhs[0]);
        *x = m;
        return;
    }

    if (job < 10) {
        switch (job) {
        case 0:
            trans_del(n);
            vv = 0.0;
            break;
        case 1:
            vv = trans_viscosity(n);
            break;
        case 2:
            vv = trans_thermalConductivity(n);
            break;
        case 3:
            vv = trans_electricalConductivity(n);
            break;
        default:
            mexErrMsgTxt("unknown Transport method");
        }
        if (vv < 0.0) {
            reportError();
        }
        plhs[0] = mxCreateNumericMatrix(1,1,mxDOUBLE_CLASS,mxREAL);
        h = mxGetPr(plhs[0]);
        *h = vv;
        return;
    } else if (job < 20) {
        nsp = getInt(prhs[3]);
        plhs[0] = mxCreateNumericMatrix(nsp,1,mxDOUBLE_CLASS,mxREAL);
        h = mxGetPr(plhs[0]);

        switch (job) {
        case 11:
            iok = trans_getMixDiffCoeffs(n, nsp, h);
            break;
        case 12:
            iok = trans_getThermalDiffCoeffs(n, nsp, h);
            break;
        default:
            mexErrMsgTxt("unknown Transport method");
        }
    } else if (job < 30) {
        nsp = getInt(prhs[3]);
        plhs[0] = mxCreateNumericMatrix(nsp,nsp,mxDOUBLE_CLASS,mxREAL);
        h = mxGetPr(plhs[0]);
        switch (job) {
        case 21:
            iok = trans_getBinDiffCoeffs(n, nsp, h);
            break;
        case 22:
            iok = trans_getMultiDiffCoeffs(n, nsp, h);
            break;
        default:
            mexErrMsgTxt("unknown Transport method");
        }
    } else if (job < 40) {
        // set parameters
        double* params;
        int typ, k;
        switch (job) {
        case 31:
            typ = getInt(prhs[3]);
            k = getInt(prhs[4]);
            params = mxGetPr(prhs[5]);
            iok = trans_setParameters(n, typ, k, params);
            break;
        default:
            mexErrMsgTxt("unknown Transport method");
        }
        plhs[0] = mxCreateNumericMatrix(1,1,mxDOUBLE_CLASS,mxREAL);
        h = mxGetPr(plhs[0]);
        *h = double(iok);
    } else {
        mexErrMsgTxt("unknown Transport method");
    }
    if (iok < 0) {
        reportError();
    }
}
