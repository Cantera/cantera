
#include "mex.h"
#include "../../../../clib/src/ct.h"
#include "../../private/ctmatutils.h"
#include <fstream.h>

extern "C" {

    void mexFunction( int nlhs, mxArray *plhs[],
        int nrhs, const mxArray *prhs[] )
    {
        double vv;
        int n = getInt(prhs[0]);
        int job = getInt(prhs[1]);
        double* h;
        int iok = 0;
        int nsp; 

        if (job < 10) {

            bool ok = true;
            switch (job) {
            case 0:
                delTransport(n);
                vv = 0.0; 
                break;
            case 1:
                vv = trans_viscosity(n);
                break;
            case 2:
                vv = trans_thermalConductivity(n); break;
            default:
                mexErrMsgTxt("unknown Transport method");                
            }
            plhs[0] = mxCreateNumericMatrix(1,1,mxDOUBLE_CLASS,mxREAL);
            h = mxGetPr(plhs[0]);
            *h = vv;
            return;
        }
        else if (job < 20) {
            nsp = getInt(prhs[2]);
            plhs[0] = mxCreateNumericMatrix(nsp,1,mxDOUBLE_CLASS,mxREAL);
            h = mxGetPr(plhs[0]);

            switch (job) {
            case 11:
                iok = trans_getMixDiffCoeffs(n, nsp, h); break;
            case 12:
                iok = trans_getThermalDiffCoeffs(n, nsp, h); break;
            default:
                mexErrMsgTxt("unknown Transport method");
            }
        }
        else if (job < 30) {
            nsp = getInt(prhs[2]);
            plhs[0] = mxCreateNumericMatrix(nsp,nsp,mxDOUBLE_CLASS,mxREAL);
            h = mxGetPr(plhs[0]);
            switch (job) {
            case 21:
                iok = trans_getBinDiffCoeffs(n, nsp, h); break;
            case 22:
                iok = trans_getMultiDiffCoeffs(n, nsp, h); break;
            default:
                mexErrMsgTxt("unknown Transport method");
            }
        }
        else {
            mexErrMsgTxt("unknown Transport method");
        }
        if (iok < 0) {
            if (iok == -1) {
                int buflen = 80;
                char* buf = (char*)mxCalloc(buflen, sizeof(char));
                getCanteraError(buflen, buf);
                mexErrMsgTxt(buf);
            }
            else {
                mexErrMsgTxt("exception thrown.");
            }
        }
    }
}
