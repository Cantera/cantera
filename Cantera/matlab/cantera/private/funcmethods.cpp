
#include "mex.h"
#include "ctmatutils.h"
#include "../../../clib/src/ctfunc.h"


void funcmethods( int nlhs, mxArray *plhs[],
    int nrhs, const mxArray *prhs[] )
{
    int job = getInt(prhs[1]);
    int nn;

    // constructor
    if (job == 0) {
        int type = getInt(prhs[2]);
        int n = getInt(prhs[3]);
        double* ptr = mxGetPr(prhs[4]);
        int msize = mxGetM(prhs[4]);
        int nsize = mxGetN(prhs[4]);
        int lenp = msize*nsize;
        nn = func_new(type, n, lenp, ptr);
        plhs[0] = mxCreateNumericMatrix(1,1,mxDOUBLE_CLASS,mxREAL);
        double *h = mxGetPr(plhs[0]);
        *h = double(nn);
        if (nn < 0) reportError();
        return;
    }

    else {
        int nn = 0;
        double t;
        double v;
        int i = getInt(prhs[2]);
        if (job == 1) {
            nn = func_del(int i);
            if (nn < 0) reportError();
            v = double(nn);
        }
        else if (job == 2) {
            t = getDouble(prhs[3]);
            v = func_value(i, t);
            if (v == Undef) reportError();
        }
        else {
            mexErrMsgTxt("unknown job parameter");
        }
        plhs[0] = mxCreateNumericMatrix(1,1,mxDOUBLE_CLASS,mxREAL);
        double *h = mxGetPr(plhs[0]);
        *h = v;
        return;
    }
}

