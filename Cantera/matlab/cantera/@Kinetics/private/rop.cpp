
#include "mex.h"
#include "../../../../clib/src/ct.h"
#include "../../private/ctmatutils.h"

extern "C" {

    void mexFunction( int nlhs, mxArray *plhs[],
        int nrhs, const mxArray *prhs[] )
    {
        int kin  = getInt(prhs[0]);
        int job = getInt(prhs[1]);
        int nr = kin_nReactions(kin);
        plhs[0] = mxCreateNumericMatrix(nr,1,mxDOUBLE_CLASS,mxREAL);
        double *h = mxGetPr(plhs[0]);
        int ok = -10;
        switch (job) {
        case 0:
            ok = kin_getFwdRatesOfProgress(kin,nr,h); break;
        case 1:
            ok = kin_getRevRatesOfProgress(kin,nr,h); break;
        case 2:
            ok = kin_getNetRatesOfProgress(kin,nr,h); break;
        case 3:
            ok = kin_getEquilibriumConstants(kin,nr,h); break;
        default:
            ;
        }
        if (ok < 0)
            mexErrMsgTxt("error computing rates of progress");
    }
}
