// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#ifndef CT_MATUTILS_H
#define CT_MATUTILS_H

const double Undef = -999.123;

#include <string>

// Workaround for VS2010 and Matlab 2010a.
// mex.h must be included after <string> or another include from
// the standard library which includes <yvals.h>.
#if (_MSC_VER >= 1600)
#define __STDC_UTF_16__
#endif
#include "mex.h"

void reportError();

void checkNArgs(const int n, const int nrhs);

template<class A>
inline int getInt(A* mxhndl)
{
    return int(mxGetScalar(mxhndl));
}

template<class A>
inline double getDouble(A* mxhndl)
{
    return double(mxGetScalar(mxhndl));
}

inline char* getString(const mxArray* p)
{
    char* input_buf = 0;
    int status;
    size_t m = mxGetM(p);
    size_t n = mxGetN(p);
    mwSize buflen = (mwSize)(m*n + 1);
    std::string msg;

    if (m == 1) {
        input_buf = (char*)mxCalloc(buflen, sizeof(char));
        status = mxGetString(p, input_buf, buflen);
        if (status != 0) {
            msg = std::string(input_buf)
                  + "\nNot enough space. String is truncated.";
            mexWarnMsgTxt(msg.c_str());
        }
    } else {
        mexErrMsgTxt("string must be a row vector");
    }
    return input_buf;
}

#endif
