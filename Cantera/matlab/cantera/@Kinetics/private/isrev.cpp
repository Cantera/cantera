
#include "mex.h"
#include "../../../../clib/src/ct.h"
#include "../../private/ctmatutils.h"

extern "C" {

    void mexFunction( int nlhs, mxArray *plhs[],
        int nrhs, const mxArray *prhs[] )
    {
        double vv;
        int kin  = getInt(prhs[0]);
        int irxn = getInt(prhs[1]);
        vv = kin_isReversible(kin,irxn-1);
        plhs[0] = mxCreateNumericMatrix(1,1,mxDOUBLE_CLASS,mxREAL);
        double *h = mxGetPr(plhs[0]);
        *h = vv;
    }
}
