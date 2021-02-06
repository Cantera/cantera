/**
 *  @file phasemethods.cpp
 */

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#include "ctmatutils.h"
#include "cantera/clib/ct.h"

#include <vector>

void phasemethods(int nlhs, mxArray* plhs[],
                  int nrhs, const mxArray* prhs[])
{
    double vv = 0.0;
    int iok=0, k;
    int ph = getInt(prhs[1]);
    int job = getInt(prhs[2]);

    char* input_buf;
    double* ptr = 0;
    size_t nsp, n, m;
    int mjob, show_thermo;

    // methods to set attributes
    if (job < 0) {
        mjob = -job;
        if (mxIsChar(prhs[3]) != 1) {
            ptr = mxGetPr(prhs[3]);
        }
        m = mxGetM(prhs[3]);
        n = mxGetN(prhs[3]);

        nsp = thermo_nSpecies(ph);

        // set scalar attributes
        if (mjob < 10) {
            if (m != 1 || n != 1) {
                mexErrMsgTxt("value must be scalar.");
            }

            switch (mjob) {
            case 1:
                iok = thermo_setTemperature(ph,*ptr);
                break;
            case 2:
                iok = thermo_setDensity(ph,*ptr);
                break;
            default:
                mexErrMsgTxt("Unknown job number");
            }
        } else if (mjob < 30) {
            // set array attributes
            if ((m == nsp && n == 1) || (m == 1 && n == nsp)) {
                int norm = 1;
                switch (mjob) {
                case 20:
                    iok = thermo_setMoleFractions(ph, nsp, ptr, norm);
                    break;
                case 21:
                    iok = thermo_setMassFractions(ph, nsp, ptr, norm);
                    break;
                case 22:
                    norm = 0;
                    iok = thermo_setMoleFractions(ph, nsp, ptr, norm);
                    break;
                case 23:
                    norm = 0;
                    iok = thermo_setMassFractions(ph, nsp, ptr, norm);
                    break;
                default:
                    mexErrMsgTxt("Unknown job number");
                }
            } else {
                mexErrMsgTxt("wrong array size");
            }
        } else {
            // set attributes from a string
            int status;
            mwSize buflen;
            char* input_buf;
            if (mxIsChar(prhs[3]) == 1) {
                if (mxGetM(prhs[3]) != 1) {
                    mexErrMsgTxt("Input must be a row vector.");
                }
                buflen = (mwSize)(mxGetM(prhs[3]) * mxGetN(prhs[3])) + 1;
                input_buf = (char*)mxCalloc(buflen, sizeof(char));
                status = mxGetString(prhs[3], input_buf, buflen);
                if (status != 0) {
                    mexWarnMsgTxt("Not enough space. " "String is truncated.");
                }

                switch (mjob) {
                case 30:
                    iok = thermo_setMoleFractionsByName(ph, input_buf);
                    break;
                case 31:
                    iok = thermo_setMassFractionsByName(ph, input_buf);
                    break;
                case 32:
                    iok = thermo_setName(ph, input_buf);
                    break;
                default:
                    mexErrMsgTxt("Unknown job number");
                }
            } else {
                mexErrMsgTxt("expected a string.");
            }
        }
    } else if (job < 20) {
        double threshold;
        switch (job) {
        case 0:
            vv = (double) thermo_newFromXML(ph);
            break;
        case 1:
            // floating-point attributes
            vv = thermo_temperature(ph);
            if (vv == DERR) {
                reportError();
            }
            break;
        case 2:
            vv = thermo_density(ph);
            if (vv == DERR) {
                reportError();
            }
            break;
        case 3:
            vv = thermo_molarDensity(ph);
            if (vv == DERR) {
                reportError();
            }
            break;
        case 4:
            vv = thermo_meanMolecularWeight(ph);
            if (vv == DERR) {
                reportError();
            }
            break;
        case 10:
            vv = static_cast<int>(thermo_nElements(ph));
            if (vv == -1) {
                reportError();
            }
            break;
        case 11:
            vv = static_cast<int>(thermo_nSpecies(ph));
            if (vv == -1) {
                reportError();
            }
            break;
        case 12:
            input_buf = getString(prhs[3]);
            vv = static_cast<int>(thermo_speciesIndex(ph, input_buf)) + 1;
            break;
        case 13:
            input_buf = getString(prhs[3]);
            vv = static_cast<int>(thermo_elementIndex(ph, input_buf)) + 1;
            break;
        case 14:
            k = getInt(prhs[3]);
            m = getInt(prhs[4]);
            vv = thermo_nAtoms(ph,k-1,m-1);
            if (vv == ERR) {
                reportError();
            }
            break;
        case 15:
            show_thermo = getInt(prhs[3]);
            threshold = getDouble(prhs[4]);
            vv = thermo_print(ph,show_thermo,threshold);
            if (vv == -1 || vv == ERR) {
                reportError();
            }
            break;
        default:
            mexErrMsgTxt("Unknown job number");
        }
        plhs[0] = mxCreateNumericMatrix(1,1,mxDOUBLE_CLASS,mxREAL);
        double* h = mxGetPr(plhs[0]);
        *h = vv;
        return;
    } else if (job < 30) {
        iok = 0;
        size_t nsp = thermo_nSpecies(ph);
        std::vector<double> x(nsp);
        switch (job) {
        case 20:
            iok = thermo_getMoleFractions(ph,nsp, &x[0]);
            break;
        case 21:
            iok = thermo_getMassFractions(ph,nsp, &x[0]);
            break;
        case 22:
            iok = thermo_getMolecularWeights(ph,nsp, &x[0]);
            break;
        case 23:
            iok = thermo_getCharges(ph, nsp, &x[0]);
            break;
        default:
            mexErrMsgTxt("Unknown job number");
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
    } else if (job < 40) {
        iok = 0;
        size_t nel = thermo_nElements(ph);
        std::vector<double> x(nel);
        switch (job) {
        case 30:
            iok = thermo_getAtomicWeights(ph,nel, &x[0]);
            break;
        default:
            ;
        }
        plhs[0] = mxCreateNumericMatrix((mwSize) nel, 1, mxDOUBLE_CLASS, mxREAL);
        double* h = mxGetPr(plhs[0]);
        if (iok >= 0) {
            for (size_t i = 0; i < nel; i++) {
                h[i] = x[i];
            }
            return;
        } else {
            for (size_t i = 0; i < nel; i++) {
                h[i] = -999.99;
            }
            mexErrMsgTxt("unknown attribute");
            return;
        }
    } else if (job < 50) {
        iok = -1;
        int ksp, mel;
        int buflen;
        char* output_buf;
        switch (job) {
        case 40:
            ksp = getInt(prhs[3]);
            buflen = thermo_getSpeciesName(ph, ksp-1, 0, 0);
            if (buflen > 0) {
                output_buf = (char*)mxCalloc(buflen, sizeof(char));
                iok = thermo_getSpeciesName(ph, ksp-1, buflen, output_buf);
            }
            break;
        case 41:
            mel = getInt(prhs[3]);
            buflen = thermo_getElementName(ph, mel-1, 0, 0);
            if (buflen > 0) {
                output_buf = (char*)mxCalloc(buflen, sizeof(char));
                iok = thermo_getElementName(ph, mel-1, buflen, output_buf);
            }
            break;
        case 42:
            buflen = thermo_getName(ph, 0, 0);
            if (buflen > 0) {
                output_buf = (char*)mxCalloc(buflen, sizeof(char));
                iok = thermo_getName(ph, buflen, output_buf);
            }
            break;
        case 43:
            buflen = thermo_getEosType(ph, 0, 0);
            if (buflen > 0) {
                output_buf = (char*)mxCalloc(buflen, sizeof(char));
                iok = thermo_getEosType(ph, buflen, output_buf);
            }
            break;
        default:
            iok = -1;
        }
        if (iok == 0) {
            plhs[0] = mxCreateString(output_buf);
            return;
        } else {
            reportError();
            return;
        }
    } else {
        mexErrMsgTxt("unimplemented method.");
        return;
    }
    if (iok < 0) {
        reportError();
    }
}
