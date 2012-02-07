/**
 * @file reactornetmethods.cpp
 */
/*
 *   $Id: reactornetmethods.cpp,v 1.4 2009/07/11 16:43:13 hkmoffa Exp $
 */

#include "mex.h"
#include "../../../clib/src/ctreactor.h"
#include "../../../clib/src/ct.h"
#include "ctmatutils.h"

//const double Undef = -999.123;

    void reactornetmethods( int nlhs, mxArray *plhs[],
        int nrhs, const mxArray *prhs[] )
    {
        int iok, n;
        int job = getInt(prhs[1]);
        int i = getInt(prhs[2]);

        double r = Undef;
        double v = Undef;
        double v2 = -1.0;
        if (nrhs > 3) v = getDouble(prhs[3]);
        if (nrhs > 4) v2 = getDouble(prhs[4]);

        // constructor
        if (job == 0) {
            n = reactornet_new();
            plhs[0] = mxCreateNumericMatrix(1,1,mxDOUBLE_CLASS,mxREAL);
            double *h = mxGetPr(plhs[0]);
            *h = double(n);
            if (n < 0) reportError();
            return;
        }

        // options that do not return a value

        if (job < 20) {
            switch (job) {

            case 1:
                iok = reactornet_del(i);
                break;
            case 2:
                iok = reactornet_copy(i);
                break;
            case 3:
                iok = reactornet_assign(i,int(v));
                break;
            case 4:
                iok = reactornet_addreactor(i, int(v));
                 break;
            case 5:
                iok = reactornet_setInitialTime(i, v);
                break;
            case 6:
                iok = reactornet_setMaxTimeStep(i, v);
                break;
            case 7:
                iok = reactornet_setTolerances(i, v, v2);
                break;
            case 8:
                iok = reactornet_advance(i, v);
                break;
            default:
                mexErrMsgTxt("unknown job parameter");
            }
            plhs[0] = mxCreateNumericMatrix(1,1,mxDOUBLE_CLASS,mxREAL);
            double *h = mxGetPr(plhs[0]);
            *h = double(iok);
            if (iok < 0) reportError();
            return;
        }

        // options that return a value of type 'double'

        else if (job < 40) {
            switch (job) {
            case 21:
                r = reactornet_step(i, v);
                break;
            case 22:
                r = reactornet_time(i);
                break;
            case 23:
                r = reactornet_rtol(i);
                break;
            case 24:
                r = reactornet_atol(i);
                break;
            default:
                mexErrMsgTxt("unknown job parameter");
            }
            plhs[0] = mxCreateNumericMatrix(1,1,mxDOUBLE_CLASS,mxREAL);
            double *h = mxGetPr(plhs[0]);
            *h = r;
            if (r == Undef) reportError();
            return;
        }
    }
