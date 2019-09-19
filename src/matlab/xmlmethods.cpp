/**
 *  @file xmlmethods.cpp
 */

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#include "cantera/clib/ctxml.h"
#include "cantera/clib/ct.h"
#include "ctmatutils.h"

void reportError();

static bool nargs_ok(int job, int n)
{
    switch (n) {
    case 1:
    case 2:
    case 10:
    case 21:
    case 22:
        return (n == 2);
    case 3:
    case 4:
    case 6:
    case 7:
    case 8:
    case 9:
    case 13:
    case 14:
    case 20:
        return (n == 3);
    case 5:
    case 11:
    case 12:
        return (n == 4);
    default:
        return true;
    }
}

void xmlmethods(int nlhs, mxArray* plhs[],
                int nrhs, const mxArray* prhs[])
{
    int iok = 0;
    char* file, *key, *val, *nm;
    int job = getInt(prhs[1]);
    int i = getInt(prhs[2]);

    // Check for proper number of arguments
    if (!nargs_ok(job,nrhs-1)) {
        mexErrMsgTxt("Wrong number of inputs.");
        return;
    } else if (nlhs > 1) {
        mexErrMsgTxt("Too many output arguments");
    }

    // options that do not return a value
    if (job < 20) {
        switch (job) {
        case 0:
            nm = getString(prhs[3]);
            iok = xml_new(nm);
            break;
        case 1:
            iok = xml_del(i);
            break;
        case 4:
            file = getString(prhs[3]);
            iok = xml_build(i, file);
            break;
        case 6:
            key = getString(prhs[3]);
            iok = xml_child(i, key);
            break;
        case 8:
            key = getString(prhs[3]);
            iok = xml_findID(i, key);
            break;
        case 9:
            key = getString(prhs[3]);
            iok = xml_findByName(i, key);
            break;
        case 10:
            iok = xml_nChildren(i);
            break;
        case 11:
            key = getString(prhs[3]);
            val = getString(prhs[4]);
            iok = xml_addChild(i, key, val);
            break;
        case 13:
            file = getString(prhs[3]);
            iok = xml_write(i, file);
            break;
        case 15:
            file = getString(prhs[3]);
            iok = xml_get_XML_File(file, 0);
            break;
        default:
            mexErrMsgTxt("unknown job parameter");
        }
        plhs[0] = mxCreateNumericMatrix(1,1,mxDOUBLE_CLASS,mxREAL);
        double* h = mxGetPr(plhs[0]);
        *h = double(iok);
        if (iok < 0) {
            reportError();
        }
        return;
    }

    // options that return strings
    char* v;
    int buflen;
    iok = -1;
    switch (job) {
    case 20:
        // return an attribute
        key = getString(prhs[3]);
        buflen = xml_attrib(i, key, 0, 0);
        if (buflen > 0) {
            v = (char*) mxCalloc(buflen, sizeof(char));
            iok = xml_attrib(i, key, buflen, v);
        }
        break;
    case 21:
        // return the value of the node
        buflen = xml_value(i, 0, 0);
        if (buflen > 0) {
            v = (char*) mxCalloc(buflen, sizeof(char));
            iok = xml_value(i, buflen, v);
        }
        break;
    default:
        mexErrMsgTxt("unknown job parameter");
    }
    if (iok < 0) {
        plhs[0] = mxCreateNumericMatrix(1,1,mxDOUBLE_CLASS,mxREAL);
        double* h = mxGetPr(plhs[0]);
        *h = double(iok);
        if (iok < 0) {
            reportError();
        }
    } else {
        plhs[0] = mxCreateString(v);
    }
}
