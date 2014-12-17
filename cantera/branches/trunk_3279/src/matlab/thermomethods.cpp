/**
 *  @file thermomethods.cpp
 */

#include "clib/ct.h"
#include "ctmatutils.h"

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
            ierr = delThermo(th);
            break;
        case 1:
            ierr = th_setPressure(th,*ptr);
            break;
        case 2:
            ierr = th_setElectricPotential(th, *ptr);
            break;
        default:
            mexErrMsgTxt("unknown attribute.");
        }
    }

    // property pairs
    else if (job < 40) {
        if ((m == 2 && n == 1) || (m == 1 && n == 2)) {
            switch (job) {
            case 20:
                ierr = th_set_HP(th,ptr);
                break;
            case 21:
                ierr = th_set_UV(th,ptr);
                break;
            case 22:
                ierr = th_set_SV(th,ptr);
                break;
            case 23:
                ierr = th_set_SP(th,ptr);
                break;
            case 24:
                ierr = th_setState_Psat(th,ptr[0],ptr[1]);
                break;
            case 25:
                ierr = th_setState_Tsat(th,ptr[0],ptr[1]);
                break;
            default:
                mexErrMsgTxt("unknown pair attribute.");
            }
        } else {
            mexErrMsgTxt("wrong size");
        }
    }

    // equilibrate
    else if (job == 50) {
        char* xy = getString(prhs[3]);
        int solver = getInt(prhs[4]);
        double rtol = getDouble(prhs[5]);
        int maxsteps = getInt(prhs[6]);
        int maxiter = getInt(prhs[7]);
        int loglevel = getInt(prhs[8]);
        ierr = th_equil(th, xy, solver, rtol, maxsteps, maxiter, loglevel);
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

    if (job < 30) {
        bool ok = true;
        switch (job) {
        case 0:
            vv = (double) newThermoFromXML(n);
            break;
        case 2:
            vv = th_enthalpy_mole(n);
            break;
        case 3:
            vv = th_intEnergy_mole(n);
            break;
        case 4:
            vv = th_entropy_mole(n);
            break;
        case 5:
            vv = th_gibbs_mole(n);
            break;
        case 6:
            vv = th_cp_mole(n);
            break;
        case 7:
            vv = th_cv_mole(n);
            break;
        case 8:
            vv = th_pressure(n);
            break;
        case 9:
            vv = th_enthalpy_mass(n);
            break;
        case 10:
            vv = th_intEnergy_mass(n);
            break;
        case 11:
            vv = th_entropy_mass(n);
            break;
        case 12:
            vv = th_gibbs_mass(n);
            break;
        case 13:
            vv = th_cp_mass(n);
            break;
        case 14:
            vv = th_cv_mass(n);
            break;
        case 15:
            vv = th_refPressure(n);
            break;
        case 16:
            vv = th_minTemp(n);
            break;
        case 17:
            vv = th_maxTemp(n);
            break;
        case 18:
            vv = double(th_eosType(n));
            break;
        case 19:
            vv  = th_critTemperature(n);
            break;
        case 20:
            vv = th_critPressure(n);
            break;
        case 21:
            vv = th_critDensity(n);
            break;
        case 22:
            vv = th_vaporFraction(n);
            break;
        case 23:
            psat = getDouble(prhs[3]);
            vv = th_satTemperature(n, psat);
            break;
        case 24:
            tsat = getDouble(prhs[3]);
            vv = th_satPressure(n, tsat);
            break;
        case 25:
            vv = th_electricPotential(n);
                break;
        case 26:
            vv = th_isothermalCompressibility(n);
                break;
        case 27:
            vv = th_thermalExpansionCoeff(n);
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
        size_t nsp = th_nSpecies(n);
        std::vector<double> x(nsp);
        switch (job) {
        case 32:
            iok = th_getEnthalpies_RT(n, nsp, &x[0]);
            break;
        case 34:
            iok = th_chemPotentials(n, nsp, &x[0]);
            break;
        case 36:
            iok = th_getEntropies_R(n, nsp, &x[0]);
            break;
        case 38:
            iok = th_getCp_R(n, nsp, &x[0]);
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
    }

    else {
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
