
#include "mex.h"
#include "../../../../clib/src/ct.h"
#include "../../private/ctmatutils.h"

extern "C" {

    void mexFunction( int nlhs, mxArray *plhs[],
        int nrhs, const mxArray *prhs[] )
    {
        int kin  = getInt(prhs[0]);
        int irxn = getInt(prhs[1]);
        char* buf;
        int buflen = 80;
        buf = (char*)mxCalloc(buflen, sizeof(char));
        int iok = kin_getReactionString(kin, irxn-1, buflen, buf);
        if (iok >= 0) {
            plhs[0] = mxCreateString(buf);
            return;
        }
        else {
            mexErrMsgTxt("error getting reaction string");
        }
    }
}
