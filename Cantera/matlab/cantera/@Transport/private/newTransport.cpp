
#include "mex.h"
#include "../../../../clib/src/ct.h"
#include "../../private/ctmatutils.h"

extern "C" {

    static void reportError() {
        int buflen = 300;
        char* output_buf = (char*)mxCalloc(buflen, sizeof(char));
        getCanteraError(buflen, output_buf);
        mexErrMsgTxt(output_buf);
    }

    /*
     * Create a Cantera 'Transport' object
     */
    void mexFunction( int nlhs, mxArray *plhs[],
        int nrhs, const mxArray *prhs[] )
    {
        try {

            // Check for proper number of arguments
            if(nrhs != 3) {
                mexErrMsgTxt("Three inputs required.");
            } 
            else if(nlhs > 1) {
                mexErrMsgTxt("Too many output arguments");
            }
        
            //int model = getInt(prhs[0]);
            int ph = getInt(prhs[1]);
            int loglevel = getInt(prhs[2]);
            
            int buflen = (mxGetM(prhs[0]) * mxGetN(prhs[0])) + 1;
            char* model = (char*)mxCalloc(buflen, sizeof(char));
            int status = mxGetString(prhs[0],model,buflen);
            if (status != 0) 
                mexErrMsgTxt("error reading model.");
            
            int n = newTransport(model, ph, loglevel);
            if (n < 0) mexErrMsgTxt("Unknown transport model");
            
            // Create matrix for the return argument.
            plhs[0] = mxCreateDoubleMatrix(1,1, mxREAL);
            double* x = mxGetPr(plhs[0]);
            *x = n;
        }
        catch (...) {
            reportError();
            //mexErrMsgTxt("exception raised in newTransport.");
        }
    }
}
