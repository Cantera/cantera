
#include "mex.h"
#include "../../private/ctmatutils.h"
#include "../../../../clib/src/ct.h"

extern "C" {

    void reportError() {
        int buflen = 300;
        char* output_buf = (char*)mxCalloc(buflen, sizeof(char));
        getCanteraError(buflen, output_buf);
        mexErrMsgTxt(output_buf);
    }

    void mexFunction( int nlhs, mxArray *plhs[],
        int nrhs, const mxArray *prhs[] ) {
        // try {
            double vv;
            int kin = getInt(prhs[0]);
            int job = getInt(prhs[1]);
            int irxn = getInt(prhs[2]);
            double* ptr = mxGetPr(prhs[3]);
            int m = mxGetM(prhs[3]);
            int n = mxGetN(prhs[3]);
            
            // set scalar attributes
            int iok = -1;
            if (job < 10) {
                if (m != 1 || n != 1)
                    mexErrMsgTxt("value must be scalar.");
                
                switch (job) {
                case 1:
                    iok = kin_setMultiplier(kin,irxn-1,*ptr); break; 
                default:
                    iok = -1;
                }
            }
            if (iok < 0) mexErrMsgTxt("error in kin_set.");
            //       }
    //catch (...) {
    //      reportError();
            //            mexErrMsgTxt("exception in kin_set.");
    //      return;
    //  }
    }
}
