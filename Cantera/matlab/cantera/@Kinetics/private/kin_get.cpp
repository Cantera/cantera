
#include "mex.h"
#include "../../../../clib/src/ct.h"
#include "../../private/ctmatutils.h"

extern "C" {

    void mexFunction( int nlhs, mxArray *plhs[],
        int nrhs, const mxArray *prhs[] )
    {
        double vv;
        int kin  = getInt(prhs[0]);
        int job = getInt(prhs[1]);
        int irxn = getInt(prhs[2]);

        bool ok = true;
        switch (job) {
            
        case 1:
            vv = kin_nReactions(kin); break;
        case 2:
            vv = kin_multiplier(kin, irxn-1); break;
        case 3:
            vv = kin_nSpecies(kin); break;
        default:
            ok = false;
        }
        if (ok) {
            plhs[0] = mxCreateNumericMatrix(1,1,mxDOUBLE_CLASS,mxREAL);
            double *h = mxGetPr(plhs[0]);
            *h = vv;
            return;
        }
    }
}
