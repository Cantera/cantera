/**
 *  @file xmlmethods.cpp
 */

#include "clib/ctxml.h"
#include "clib/ct.h"
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
    int j, m, iok = 0;
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
        case 2:
            iok = xml_copy(i);
            break;
        case 3:
            j = getInt(prhs[3]);
            iok = xml_assign(i,j);
            break;
        case 4:
            file = getString(prhs[3]);
            iok = xml_build(i, file);
            break;
        case 5:
            key = getString(prhs[3]);
            val = getString(prhs[4]);
            iok = xml_addAttrib(i, key, val);
            break;
        case 6:
            key = getString(prhs[3]);
            iok = xml_child(i, key);
            break;
        case 7:
            m = getInt(prhs[3]);
            iok = xml_child_bynumber(i, m);
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
        case 12:
            key = getString(prhs[3]);
            j = getInt(prhs[4]);
            iok = xml_addChildNode(i, j);
            break;
        case 13:
            file = getString(prhs[3]);
            iok = xml_write(i, file);
            break;
        case 14:
            j = getInt(prhs[3]);
            iok = xml_removeChild(i, j);
            break;
        case 15:
            file = getString(prhs[3]);
            iok = xml_get_XML_File(file);
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
    char* v = (char*)mxCalloc(80, sizeof(char));
    switch (job) {
    case 20:
        // return an attribute
        key = getString(prhs[3]);
        iok = xml_attrib(i, key, v);
        break;
    case 21:
        // return the value of the node
        iok = xml_value(i, v);
        break;
    case 22:
        iok = xml_tag(i, v);
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
