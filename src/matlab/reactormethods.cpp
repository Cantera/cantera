/**
 * @file reactormethods.cpp
 */

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#include "cantera/clib/ctreactor.h"
#include "cantera/clib/ct.h"
#include "ctmatutils.h"

void reactormethods(int nlhs, mxArray* plhs[],
                    int nrhs, const mxArray* prhs[])
{
    int iok = 0, n;
    int job = getInt(prhs[1]);
    double r = Undef;
    double v = Undef;
    if (nrhs > 3) {
        v = getDouble(prhs[3]);
    }

    // constructor
    if (job == 0) {
        char* type = getString(prhs[2]);

        n = reactor_new2(type);
        plhs[0] = mxCreateNumericMatrix(1,1,mxDOUBLE_CLASS,mxREAL);
        double* h = mxGetPr(plhs[0]);
        *h = double(n);
        if (n < 0) {
            reportError();
        }
        return;
    }

    // options that do not return a value

    int i = getInt(prhs[2]);
    if (job < 20) {
        switch (job) {
        case 1:
            iok = reactor_del(i);
            break;
        case 4:
            iok = reactor_setInitialVolume(i, v);
            break;
        case 6:
            iok = reactor_setThermoMgr(i, int(v));
            break;
        case 7:
            iok = reactor_setKineticsMgr(i, int(v));
            break;
        case 8:
            iok = reactor_setChemistry(i, bool(v));
            break;
        case 9:
            iok = reactor_setEnergy(i, int(v));
            break;
        case 10:
            iok = flowReactor_setMassFlowRate(i, v);
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
    } else if (job < 40) {
        // options that return a value of type 'double'
        switch (job) {
        case 23:
            r = reactor_mass(i);
            break;
        case 24:
            r = reactor_volume(i);
            break;
        case 25:
            r = reactor_density(i);
            break;
        case 26:
            r = reactor_temperature(i);
            break;
        case 27:
            r = reactor_enthalpy_mass(i);
            break;
        case 28:
            r = reactor_intEnergy_mass(i);
            break;
        case 29:
            r = reactor_pressure(i);
            break;
        case 30:
            r = reactor_massFraction(i, int(v));
            break;
        default:
            mexErrMsgTxt("unknown job parameter");
        }
        plhs[0] = mxCreateNumericMatrix(1,1,mxDOUBLE_CLASS,mxREAL);
        double* h = mxGetPr(plhs[0]);
        *h = r;
        if (r == DERR) {
            reportError();
        }
        return;
    }
}
