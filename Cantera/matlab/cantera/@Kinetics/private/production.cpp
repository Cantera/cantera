
#include "mex.h"
#include "../../../../clib/src/ct.h"
#include "../../private/ctmatutils.h"

extern "C" {

    void mexFunction( int nlhs, mxArray *plhs[],
        int nrhs, const mxArray *prhs[] )
    {
        int kin  = getInt(prhs[0]);
        int nsp  = getInt(prhs[1]);
        int job = getInt(prhs[2]);
        plhs[0] = mxCreateNumericMatrix(nsp,1,mxDOUBLE_CLASS,mxREAL);
        double *h = mxGetPr(plhs[0]);
        int ok = -10;
        switch (job) {
        case 0:
            ok = kin_getCreationRates(kin,nsp,h); break;
        case 1:
            ok = kin_getDestructionRates(kin,nsp,h); break;
        case 2:
            ok = kin_getNetProductionRates(kin,nsp,h); break;
        case 3:
            ok = kin_getSourceTerms(kin, nsp, h); break;
        default:
            ;
        }
        if (ok < 0)
            mexErrMsgTxt("error computing production rates");
    }
}
