
#include "mex.h"
#include "../../../../clib/src/ctxml.h"
#include "../../private/ctmatutils.h"


extern "C" {

    /*
     * Create a Cantera 'XML_Node' object
     */
    void mexFunction( int nlhs, mxArray *plhs[],
        int nrhs, const mxArray *prhs[] )
    {
        // Check for proper number of arguments
        if(nrhs != 1) {
            mexErrMsgTxt("One input required.");
        } 
        else if(nlhs > 1) {
            mexErrMsgTxt("Too many output arguments");
        }        
        char* input_buf = getString(prhs[0]);
        int id = 0;
        //try {
            id = xml_new(input_buf);
            //}
            //catch (...) {
            //mexErrMsgTxt("exception occurred");
            //}
        plhs[0] = mxCreateNumericMatrix(1,1,mxDOUBLE_CLASS,mxREAL);
        double *h = mxGetPr(plhs[0]);
        *h = double(id);
    }
}
