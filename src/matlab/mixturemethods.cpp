/**
 * @file mixturemethods.cpp
 */

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#include <iostream>
#include <vector>

#include "cantera/clib/ctmultiphase.h"
#include "cantera/clib/ct.h"
#include "ctmatutils.h"
using namespace std;

void mixturemethods(int nlhs, mxArray* plhs[],
                    int nrhs, const mxArray* prhs[])
{
    int m, iok = 0, n;
    int job = getInt(prhs[1]);
    int i = getInt(prhs[2]);

    double r = Undef;
    double v = Undef;
    if (nrhs > 3 && job != 8 && job != 9 && job != 22 && job != 23) {
        v = getDouble(prhs[3]);
    }

    // constructor
    if (job == 0) {
        n = mix_new();
        plhs[0] = mxCreateNumericMatrix(1,1,mxDOUBLE_CLASS,mxREAL);
        double* h = mxGetPr(plhs[0]);
        *h = double(n);
        if (n < 0) {
            reportError();
        }
        return;
    }

    // options that do not return a value
    double moles, err;
    char* nmstr, *XY, *nm;
    int maxiter, maxsteps, loglevel;
    if (job < 15) {
        switch (job) {
        case 1:
            iok = mix_del(i);
            break;
        case 4:
            checkNArgs(5, nrhs);
            moles = getDouble(prhs[4]);
            iok = mix_addPhase(i, int(v), moles);
            break;
        case 5:
            iok = mix_setTemperature(i, v);
            break;
        case 6:
            iok = mix_setPressure(i, v);
            break;
        case 7:
            checkNArgs(5, nrhs);
            moles = getDouble(prhs[4]);
            iok = mix_setPhaseMoles(i, int(v)-1, moles);
            break;
        case 8:
            checkNArgs(4, nrhs);
            nmstr = getString(prhs[3]);
            iok = mix_setMolesByName(i, nmstr);
            break;
        case 9:
            iok = mix_updatePhases(i);
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
        case 19:
            r = (double) mix_nPhases(i);
            break;
        case 21:
            r = (double) mix_nElements(i);
            break;
        case 22:
            checkNArgs(4, nrhs);
            nm = getString(prhs[3]);
            r = (double) mix_elementIndex(i, nm)+1;
            break;
        case 23:
            checkNArgs(5, nrhs);
            m = getInt(prhs[3]);
            n = getInt(prhs[4]);
            r = (double) mix_speciesIndex(i, m-1, n-1)+1;
            break;
        case 24:
            r = (double) mix_nSpecies(i);
            break;
        case 25:
            r = mix_temperature(i);
            break;
        case 26:
            r = mix_pressure(i);
            break;
        case 27:
            m = getInt(prhs[4]);
            r = mix_nAtoms(i,int(v), m);
            break;
        case 28:
            r = mix_phaseMoles(i, int(v)-1);
            break;
        case 29:
            r = mix_speciesMoles(i, int(v)-1);
            break;
        case 30:
            r = mix_elementMoles(i, int(v)-1);
            break;
        case 31:
            checkNArgs(8, nrhs);
            XY = getString(prhs[3]);
            err = getDouble(prhs[4]);
            maxsteps = getInt(prhs[5]);
            maxiter = getInt(prhs[6]);
            loglevel = getInt(prhs[7]);
            r = mix_equilibrate(i, XY, err, maxsteps,
                                maxiter, loglevel);
            break;
        default:
            mexErrMsgTxt("unknown job parameter");
        }
        plhs[0] = mxCreateNumericMatrix(1,1,mxDOUBLE_CLASS,mxREAL);
        double* h = mxGetPr(plhs[0]);
        *h = r;
        if (r == DERR || r == Undef || r == -1) {
            reportError();
        }
        return;
    } else if (job < 60) {
        // species properties
        int iok = 0;
        mwSize nsp = (mwSize) mix_nSpecies(i);
        std::vector<double> x(nsp);
        switch (job) {
        case 41:
            iok = mix_getChemPotentials(i,nsp, &x[0]);
            break;
        default:
            ;
        }
        plhs[0] = mxCreateNumericMatrix(nsp,1, mxDOUBLE_CLASS,mxREAL);
        double* h = mxGetPr(plhs[0]);
        if (iok >= 0) {
            for (mwSize i = 0; i < nsp; i++) {
                h[i] = x[i];
            }
            return;
        } else {
            for (mwSize i = 0; i < nsp; i++) {
                h[i] = -999.99;
            }
            mexErrMsgTxt("unknown attribute");
            return;
        }
    }
}
