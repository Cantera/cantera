
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
        // Check for proper number of arguments
        if(nrhs != 4) {
            mexErrMsgTxt("Four inputs required.");
        } 
        else if(nlhs > 1) {
            mexErrMsgTxt("Too many output arguments");
        }
        
        int root = getInt(prhs[0]);
        int iph = getInt(prhs[1]);
        int in1 = getInt(prhs[2]);
        int in2 = getInt(prhs[3]);
        int n = newKineticsFromXML(root, iph, in1, in2);

        // Create matrix for the return argument.
        plhs[0] = mxCreateDoubleMatrix(1,1, mxREAL);
        double* x = mxGetPr(plhs[0]);
        *x = n;
    }
}
