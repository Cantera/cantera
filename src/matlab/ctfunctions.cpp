/**
 *   @file ctfunctions.cpp
 */

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#include <string.h>
#include <iostream>

#include "cantera/clib/ct.h"
#include "cantera/clib/ctreactor.h"
#include "cantera/clib/ctonedim.h"
#include "cantera/clib/ctmultiphase.h"
#include "cantera/clib/ctxml.h"
#include "cantera/clib/ctfunc.h"
#include "cantera/clib/ctrpath.h"
#include "ctmatutils.h"

using namespace std;

void reportError()
{
    int buflen = 0;
    char* output_buf = 0;
    buflen = ct_getCanteraError(buflen, output_buf) + 1;
    output_buf = (char*)mxCalloc(buflen, sizeof(char));
    ct_getCanteraError(buflen, output_buf);
    mexErrMsgTxt(output_buf);
}

void ctfunctions(int nlhs, mxArray* plhs[],
                 int nrhs, const mxArray* prhs[])
{
    int job = getInt(prhs[1]);
    int iok = 0, dbg, validate;
    char* infile, *dbfile, *trfile, *idtag, *sep;
    int buflen = 0;
    char* output_buf = 0;

    switch (job) {
        // convert CK file to CTI
    case 1:
        if (nrhs < 8) {
            mexErrMsgTxt("Wrong number of inputs.");
            return;
        }
        infile = getString(prhs[2]);
        dbfile = getString(prhs[3]);
        trfile = getString(prhs[4]);
        idtag = getString(prhs[5]);
        dbg = getInt(prhs[6]);
        validate = getInt(prhs[7]);
        iok = ct_ck2cti(infile, dbfile, trfile, idtag, dbg, validate);
        break;

        // get Cantera error
    case 2:
        buflen = ct_getCanteraError(buflen, output_buf) + 1;
        output_buf = (char*)mxCalloc(buflen, sizeof(char));
        iok = ct_getCanteraError(buflen, output_buf);
        plhs[0] = mxCreateString(output_buf);
        iok = 0;
        return;

        // add directory
    case 3:
        infile = getString(prhs[2]);
        iok = ct_addCanteraDirectory(strlen(infile), infile);
        break;

        // clear storage
    case 4:
        iok = ct_clearOneDim();
        iok = ct_clearMix();
        iok = ct_clearXML();
        iok = ct_clearFunc();
        iok = ct_clearStorage();
        iok = ct_clearReactors();
        iok = ct_clearReactionPath();
        break;

        // get string of data directories
    case 5:
        sep = getString(prhs[2]);
        buflen = ct_getDataDirectories(0, 0, sep);
        output_buf = (char*)mxCalloc(buflen, sizeof(char));
        iok = ct_getDataDirectories(buflen, output_buf, sep);
        plhs[0] = mxCreateString(output_buf);
        iok = 0;
        return;

        // get cantera version string
    case 6:
        buflen = ct_getCanteraVersion(0, 0);
        output_buf = (char*)mxCalloc(buflen, sizeof(char));
        iok = ct_getCanteraVersion(buflen, output_buf);
        plhs[0] = mxCreateString(output_buf);
        return;

        // get cantera git commit
    case 7:
        buflen = ct_getGitCommit(0, 0);
        output_buf = (char*)mxCalloc(buflen, sizeof(char));
        iok = ct_getGitCommit(buflen, output_buf);
        plhs[0] = mxCreateString(output_buf);
        return;

    default:
        mexErrMsgTxt("ctfunctions: unknown job");
    }
    plhs[0] = mxCreateNumericMatrix(1,1,mxDOUBLE_CLASS,mxREAL);
    double* h = mxGetPr(plhs[0]);
    *h = double(iok);
    if (iok < 0) {
        reportError();
    }
}
