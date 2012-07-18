/**
 *  @file phasemethods.cpp
 */

#include "ctmatutils.h"
#include "clib/ct.h"

void phasemethods(int nlhs, mxArray* plhs[],
                  int nrhs, const mxArray* prhs[])
{
    double vv = 0.0;
    int iok=0, k;
    int ph  = getInt(prhs[1]);
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

        nsp = phase_nSpecies(ph);

        // set scalar attributes
        if (mjob < 10) {
            if (m != 1 || n != 1) {
                mexErrMsgTxt("value must be scalar.");
            }

            switch (mjob) {
            case 1:
                iok = phase_setTemperature(ph,*ptr);
                break;
            case 2:
                iok = phase_setDensity(ph,*ptr);
                break;
            default:
                mexErrMsgTxt("Unknown job number");
            }
        }

        // set array attributes
        else if (mjob < 30) {
            if ((m == nsp && n == 1) || (m == 1 && n == nsp)) {
                int norm = 1;
                switch (mjob) {
                case 20:
                    iok = phase_setMoleFractions(ph, nsp, ptr, norm);
                    break;
                case 21:
                    iok = phase_setMassFractions(ph, nsp, ptr, norm);
                    break;
                case 22:
                    norm = 0;
                    iok = phase_setMoleFractions(ph, nsp, ptr, norm);
                    break;
                case 23:
                    norm = 0;
                    iok = phase_setMassFractions(ph, nsp, ptr, norm);
                    break;
                default:
                    mexErrMsgTxt("Unknown job number");
                }
            } else {
                mexErrMsgTxt("wrong array size");
            }
        }

        // set attributes from a string
        else {
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
                    iok = phase_setMoleFractionsByName(ph, input_buf);
                    break;
                case 31:
                    iok = phase_setMassFractionsByName(ph, input_buf);
                    break;
                case 32:
                    iok = phase_setName(ph, input_buf);
                    break;
                default:
                    mexErrMsgTxt("Unknown job number");
                }
            } else {
                mexErrMsgTxt("expected a string.");
            }
        }
    }

    else if (job < 20) {

        switch (job) {
        case 0:
            vv = (double) newThermoFromXML(ph);
            break;
            // floating-point attributes
        case 1:
            vv = phase_temperature(ph);
            break;
        case 2:
            vv = phase_density(ph);
            break;
        case 3:
            vv = phase_molarDensity(ph);
            break;
        case 4:
            vv = phase_meanMolecularWeight(ph);
            break;
        case 8:
            vv = 1.0/phase_density(ph);
            break;
        case 10:
            vv = static_cast<int>(phase_nElements(ph));
            break;
        case 11:
            vv = static_cast<int>(phase_nSpecies(ph));
            break;
        case 12:
            input_buf = getString(prhs[3]);
            vv = static_cast<int>(phase_speciesIndex(ph, input_buf)) + 1;
            break;
        case 13:
            input_buf = getString(prhs[3]);
            vv = static_cast<int>(phase_elementIndex(ph, input_buf)) + 1;
            break;
        case 14:
            k = getInt(prhs[3]);
            m = getInt(prhs[4]);
            vv = phase_nAtoms(ph,k-1,m-1);
            break;
        case 15:
            show_thermo = getInt(prhs[3]);
            vv = write_phase(ph,show_thermo);
            break;
        default:
            mexErrMsgTxt("Unknown job number");
        }
        if (vv == DERR || vv == -1 || vv == ERR) {
            reportError();
        }
        plhs[0] = mxCreateNumericMatrix(1,1,mxDOUBLE_CLASS,mxREAL);
        double* h = mxGetPr(plhs[0]);
        *h = vv;
        return;
    }

    else if (job < 30) {
        iok = 0;
        size_t nsp = phase_nSpecies(ph);
        std::vector<double> x(nsp);
        switch (job) {
        case 20:
            iok = phase_getMoleFractions(ph,nsp, &x[0]);
            break;
        case 21:
            iok = phase_getMassFractions(ph,nsp, &x[0]);
            break;
        case 22:
            iok = phase_getMolecularWeights(ph,nsp, &x[0]);
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
    }

    else if (job < 40) {
        iok = 0;
        size_t nel = phase_nElements(ph);
        std::vector<double> x(nel);
        switch (job) {
        case 30:
            iok = phase_getAtomicWeights(ph,nel, &x[0]);
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
    }

    else if (job < 50) {
        iok = -1;
        int ksp, mel;
        int buflen;
        char* output_buf;
        switch (job) {
        case 40:
            ksp = getInt(prhs[3]);
            buflen = 40;
            output_buf = (char*)mxCalloc(buflen, sizeof(char));
            iok = phase_getSpeciesName(ph, ksp-1, buflen, output_buf);
            break;
        case 41:
            mel = getInt(prhs[3]);
            buflen = 40;
            output_buf = (char*)mxCalloc(buflen, sizeof(char));
            iok = phase_getElementName(ph, mel-1, buflen, output_buf);
            break;
        case 42:
            buflen = 40;
            output_buf = (char*)mxCalloc(buflen, sizeof(char));
            iok = phase_getName(ph, buflen, output_buf);
            break;
        default:
            iok = -1;
        }
        if (iok >= 0) {
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
