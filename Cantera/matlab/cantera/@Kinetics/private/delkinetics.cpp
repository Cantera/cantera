
#include "mex.h"
#include "../../../../clib/src/ct.h"
#include "../../private/ctmatutils.h"

extern "C" {

    /*
     * Create a Cantera 'Kinetics' object
     */
    void mexFunction( int nlhs, mxArray *plhs[],
        int nrhs, const mxArray *prhs[] )
    {
        int k = getInt(prhs[0]);
        int ok = delKinetics(k);

        // Create matrix for the return argument.
        plhs[0] = mxCreateDoubleMatrix(1,1, mxREAL);
        double* x = mxGetPr(plhs[0]);
        *x = ok;
    }
}
