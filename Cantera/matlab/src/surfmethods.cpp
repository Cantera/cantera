#include <iostream>
#include <string>

// Workaround for VS2010 and Matlab 2010a
#if (_MSC_VER >= 1600) && defined(_CHAR16T)
#define CHAR16_T char16_t
#endif
#include "mex.h"

#include "ctmatutils.h"
#include <cantera/clib/ctsurf.h>
#include <cantera/clib/ct.h>

using namespace std;


namespace Cantera {
    void writelog(const std::string& s);
}
using namespace Cantera;

void surfmethods( int nlhs, mxArray *plhs[],
    int nrhs, const mxArray *prhs[] ) {
    double vv;
    int job = getInt(prhs[2]); 
    int iok;
    double* ptr;
    char* str;
    size_t nsp, n, m;
    int surf = getInt(prhs[1]);

    // set parameters
    if (job < 100) {

        switch (job) {

        case 1:
            checkNArgs(4, nrhs);
            vv = getDouble(prhs[3]);
            iok = surf_setsitedensity(surf, vv);
            break;

        case 3:
            checkNArgs(4, nrhs);
            ptr = mxGetPr(prhs[3]);
            m = mxGetM(prhs[3]);
            n = mxGetN(prhs[3]);
            nsp = phase_nSpecies(surf);
            if ((m == nsp && n == 1) || (m == 1 && n == nsp)) {
                iok = surf_setcoverages(surf, ptr);
            }
            else {
                mexErrMsgTxt("wrong array size for coverages");
            }
            break;

        case 5:
            checkNArgs(4, nrhs);
            str = getString(prhs[3]);
            iok = surf_setcoveragesbyname(surf, str);
            break;

        default:
            mexErrMsgTxt("unknown job");
        }
        if (iok < 0) reportError();
        plhs[0] = mxCreateNumericMatrix(1,1,mxDOUBLE_CLASS,mxREAL);
        double *h = mxGetPr(plhs[0]);
        *h = double(iok);
        return;
    }

    // return array parameters
    else if (job < 200) {
        nsp = phase_nSpecies(surf);
        double* x = new double[nsp];

        switch (job) {
        case 101:
            checkNArgs(3,nrhs);
            iok = surf_getcoverages(surf,x);
            break;
        case 103:
            iok = surf_getconcentrations(surf,x);
            break;
        default:
            ;
        }
        plhs[0] = mxCreateNumericMatrix((mwSize) nsp,1,
            mxDOUBLE_CLASS,mxREAL);
        double *h = mxGetPr(plhs[0]);
        if (iok >= 0) {
            for (size_t i = 0; i < nsp; i++) h[i] = x[i];
            delete x;
            return;
        }
        else {
            for (size_t i = 0; i < nsp; i++) h[i] = -999.99;
            delete x;
            reportError();
            return;
        }
    }
    else {
        mexErrMsgTxt("unknown job");
    }
}
