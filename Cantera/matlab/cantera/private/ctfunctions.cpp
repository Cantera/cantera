
#include "mex.h"
#include "../../../clib/src/ct.h"
#include "ctmatutils.h"
#include <string.h>

void reportError() {
    int buflen = 300;
    char* output_buf = (char*)mxCalloc(buflen, sizeof(char));
    getCanteraError(buflen, output_buf);
    mexErrMsgTxt(output_buf);
}

void ctfunctions( int nlhs, mxArray *plhs[],
    int nrhs, const mxArray *prhs[] )
{
    int job = getInt(prhs[1]);
        
    int j, m, iok, id;
    char *file, *key, *val;

    char *infile, *dbfile, *trfile, *outfile, *idtag;
    int buflen;
    char* output_buf;

    switch (job) {

        // convert CK file to CTI
    case 1:
        if (nrhs < 6) {
            mexErrMsgTxt("Wrong number of inputs.");
            return;
        } 
        infile = getString(prhs[2]);
        dbfile = getString(prhs[3]);
        trfile = getString(prhs[4]);
        idtag = getString(prhs[5]);

        iok = ck_to_cti(infile, dbfile, trfile, idtag);
        break;  

        // get Cantera error
    case 2:
        buflen = 300;
        output_buf = (char*)mxCalloc(buflen, sizeof(char));
        iok = getCanteraError(buflen, output_buf);
        plhs[0] = mxCreateString(output_buf);
        return;

        // add directory
    case 3:
        infile = getString(prhs[2]);
        iok = addCanteraDirectory(strlen(infile), infile);
        break;

    default:
        mexErrMsgTxt("ctfunctions: unknown job");
    }
    plhs[0] = mxCreateNumericMatrix(1,1,mxDOUBLE_CLASS,mxREAL);
    double *h = mxGetPr(plhs[0]);
    *h = double(iok);
    if (iok < 0) reportError();
}
