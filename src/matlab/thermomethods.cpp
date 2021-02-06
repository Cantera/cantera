/**
 *  @file thermomethods.cpp
 */

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#include "cantera/clib/ct.h"
#include "ctmatutils.h"
#include <vector>

static void thermoset(int nlhs, mxArray* plhs[],
                      int nrhs, const mxArray* prhs[])
{
    int ierr = 0;
    int th = getInt(prhs[1]);
    int job = -getInt(prhs[2]);
    double* ptr = 0;
    if (mxIsDouble(prhs[3]) == 1) {
        ptr = mxGetPr(prhs[3]);
    }
    size_t m = mxGetM(prhs[3]);
    size_t n = mxGetN(prhs[3]);

    // scalar attributes
    if (job < 20) {
        if (m != 1 || n != 1) {
            mexErrMsgTxt("value must be scalar.");
        }
        switch (job) {
        case 10:
            ierr = thermo_del(th);
            break;
        case 1:
            ierr = thermo_setPressure(th,*ptr);
            break;
        case 2:
            ierr = thermo_setElectricPotential(th, *ptr);
            break;
        default:
            mexErrMsgTxt("unknown attribute.");
        }
    } else if (job < 40) {
        // property pairs
        if ((m == 2 && n == 1) || (m == 1 && n == 2)) {
            switch (job) {
            case 20:
                ierr = thermo_set_HP(th,ptr);
                break;
            case 21:
                ierr = thermo_set_UV(th,ptr);
                break;
            case 22:
                ierr = thermo_set_SV(th,ptr);
                break;
            case 23:
                ierr = thermo_set_SP(th,ptr);
                break;
            case 24:
                ierr = thermo_setState_Psat(th,ptr[0],ptr[1]);
                break;
            case 25:
                ierr = thermo_setState_Tsat(th,ptr[0],ptr[1]);
                break;
            case 26:
                ierr = thermo_set_RP(th,ptr);
                break;
            case 27:
                ierr = thermo_set_ST(th,ptr);
                break;
            case 28:
                ierr = thermo_set_TV(th,ptr);
                break;
            case 29:
                ierr = thermo_set_PV(th,ptr);
                break;
            case 30:
                ierr = thermo_set_UP(th,ptr);
                break;
            case 31:
                ierr = thermo_set_VH(th,ptr);
                break;
            case 32:
                ierr = thermo_set_TH(th,ptr);
                break;
            case 33:
                ierr = thermo_set_SH(th,ptr);
                break;
            default:
                mexErrMsgTxt("unknown pair attribute.");
            }
        } else {
            mexErrMsgTxt("wrong size");
        }
    } else if (job == 50) {
        // equilibrate
        char* xy = getString(prhs[3]);
        int solver = getInt(prhs[4]);
        double rtol = getDouble(prhs[5]);
        int maxsteps = getInt(prhs[6]);
        int maxiter = getInt(prhs[7]);
        int loglevel = getInt(prhs[8]);
        ierr = thermo_equilibrate(th, xy, solver, rtol, maxsteps, maxiter, loglevel);
    }
    if (ierr < 0) {
        reportError();
    }

    plhs[0] = mxCreateNumericMatrix(1,1,mxDOUBLE_CLASS,mxREAL);
    double* h = mxGetPr(plhs[0]);
    *h = 1.0*ierr;
    return;
}

static void thermoget(int nlhs, mxArray* plhs[],
                      int nrhs, const mxArray* prhs[])
{
    double vv, psat, tsat, TK;
    int n = getInt(prhs[1]);
    int job = getInt(prhs[2]);

    if (job == 0) {
        checkNArgs(5, nrhs);
        std::string fileName = getString(prhs[3]);
        std::string phaseName = getString(prhs[4]);
        if (phaseName == "-") {
            phaseName = "";
        }
        vv = (double) thermo_newFromFile(fileName.c_str(), phaseName.c_str());
        if (vv == DERR) {
            reportError();
        }
        plhs[0] = mxCreateNumericMatrix(1, 1, mxDOUBLE_CLASS, mxREAL);
        double* h = mxGetPr(plhs[0]);
        *h = vv;
        return;
    } else if (job < 30) {
        bool ok = true;
        switch (job) {
        case 2:
            vv = thermo_enthalpy_mole(n);
            break;
        case 3:
            vv = thermo_intEnergy_mole(n);
            break;
        case 4:
            vv = thermo_entropy_mole(n);
            break;
        case 5:
            vv = thermo_gibbs_mole(n);
            break;
        case 6:
            vv = thermo_cp_mole(n);
            break;
        case 7:
            vv = thermo_cv_mole(n);
            break;
        case 8:
            vv = thermo_pressure(n);
            break;
        case 9:
            vv = thermo_enthalpy_mass(n);
            break;
        case 10:
            vv = thermo_intEnergy_mass(n);
            break;
        case 11:
            vv = thermo_entropy_mass(n);
            break;
        case 12:
            vv = thermo_gibbs_mass(n);
            break;
        case 13:
            vv = thermo_cp_mass(n);
            break;
        case 14:
            vv = thermo_cv_mass(n);
            break;
        case 15:
            vv = thermo_refPressure(n);
            break;
        case 16:
            vv = thermo_minTemp(n, -1);
            break;
        case 17:
            vv = thermo_maxTemp(n, -1);
            break;
        case 19:
            vv = thermo_critTemperature(n);
            break;
        case 20:
            vv = thermo_critPressure(n);
            break;
        case 21:
            vv = thermo_critDensity(n);
            break;
        case 22:
            vv = thermo_vaporFraction(n);
            break;
        case 23:
            psat = getDouble(prhs[3]);
            vv = thermo_satTemperature(n, psat);
            break;
        case 24:
            tsat = getDouble(prhs[3]);
            vv = thermo_satPressure(n, tsat);
            break;
        case 25:
            vv = thermo_electricPotential(n);
                break;
        case 26:
            vv = thermo_isothermalCompressibility(n);
                break;
        case 27:
            vv = thermo_thermalExpansionCoeff(n);
                break;
        default:
            ok = false;
        }
        if (ok) {
            if (vv == DERR) {
                reportError();
            }
            plhs[0] = mxCreateNumericMatrix(1, 1, mxDOUBLE_CLASS, mxREAL);
            double* h = mxGetPr(plhs[0]);
            *h = vv;
            return;
        }
    } else if (job < 50) {
        int iok = 0;
        size_t nsp = thermo_nSpecies(n);
        std::vector<double> x(nsp);
        switch (job) {
        case 32:
            iok = thermo_getEnthalpies_RT(n, nsp, &x[0]);
            break;
        case 34:
            iok = thermo_chemPotentials(n, nsp, &x[0]);
            break;
        case 36:
            iok = thermo_getEntropies_R(n, nsp, &x[0]);
            break;
        case 38:
            iok = thermo_getCp_R(n, nsp, &x[0]);
            break;
        default:
            ;
        }
        plhs[0] = mxCreateNumericMatrix((mwSize) nsp, 1, mxDOUBLE_CLASS, mxREAL);
        double* h = mxGetPr(plhs[0]);
        if (iok >= 0) {
            for (size_t i = 0; i < nsp; i++) {
                h[i] = x[i];
            }
            return;
        } else {
            for (size_t i = 0; i < nsp; i++) {
                h[i] = -999.99;
            }
            mexErrMsgTxt("unknown attribute");
            return;
        }
    } else {
        mexErrMsgTxt("unknown attribute");
    }
}


void thermomethods(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[])
{
    int job = getInt(prhs[2]);
    if (job < 0) {
        thermoset(nlhs, plhs, nrhs, prhs);
    } else {
        thermoget(nlhs, plhs, nrhs, prhs);
    }
}
