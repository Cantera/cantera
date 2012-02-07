/**
 *  @file phasemethods.cpp
 */
/*
 *  $Id: phasemethods.cpp,v 1.5 2009/07/11 16:43:12 hkmoffa Exp $
 */

#include "mex.h"
#include "ctmatutils.h"
#include "../../../clib/src/ct.h"

    void phasemethods( int nlhs, mxArray *plhs[],
        int nrhs, const mxArray *prhs[] )
    {
        double vv;
        int iok=0, k, m;
        int ph  = getInt(prhs[1]);
        int job = getInt(prhs[2]);
            
        bool ok = true;
        char* input_buf;
        double* ptr = 0;
        int n, nsp, mjob, show_thermo;

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
            bool ok = true;
            if (mjob < 10) {
                if (m != 1 || n != 1)
                    mexErrMsgTxt("value must be scalar.");

                switch (mjob) {
                case 1:
                    iok = phase_setTemperature(ph,*ptr); break; 
                case 2:
                    iok = phase_setDensity(ph,*ptr); break;
                default:
                    ok = false;
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
                        ok = false;
                    }
                }
                else { 
                    mexErrMsgTxt("wrong array size");
                }
            }

            // set attributes from a string
            else {
                int buflen, status;
                char* input_buf;
                if (mxIsChar(prhs[3]) == 1) {
                    if(mxGetM(prhs[3]) != 1)
                        mexErrMsgTxt("Input must be a row vector.");
                    buflen = (mxGetM(prhs[3]) * mxGetN(prhs[3])) + 1;
                    input_buf = (char*)mxCalloc(buflen, sizeof(char));
                    status = mxGetString(prhs[3], input_buf, buflen);
                    if (status != 0) 
                        mexWarnMsgTxt("Not enough space. "
                            "String is truncated.");

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
                        mexErrMsgTxt("what?");
                    }   
                }
                else {
                    mexErrMsgTxt("expected a string.");
                }
            }
        }

        else if (job < 20) {

            switch (job) {
            case 0:
                vv = newThermoFromXML(ph); break;
                // floating-point attributes
            case 1:
                vv = phase_temperature(ph); break; 
            case 2:
                vv = phase_density(ph); break;
            case 3:
                vv = phase_molarDensity(ph); break;
            case 4:
                vv = phase_meanMolecularWeight(ph); break;
            case 8:
                vv = 1.0/phase_density(ph); break;            
            case 10:
                vv = phase_nElements(ph); break;
            case 11:
                vv = phase_nSpecies(ph); break;
            case 12:
                input_buf = getString(prhs[3]);
                vv = phase_speciesIndex(ph, input_buf) + 1;
                break;
            case 13:
                input_buf = getString(prhs[3]);
                vv = phase_elementIndex(ph, input_buf) + 1;
                break;
            case 14:
                k = getInt(prhs[3]);
                m = getInt(prhs[4]);
                vv = phase_nAtoms(ph,k-1,m-1); break;
            case 15:
                show_thermo = getInt(prhs[3]);
                vv = write_phase(ph,show_thermo); 
                break;
            default:
                ok = false;
            }
            if (ok) {
                if (vv == DERR) reportError();
                plhs[0] = mxCreateNumericMatrix(1,1,mxDOUBLE_CLASS,mxREAL);
                double *h = mxGetPr(plhs[0]);
                *h = vv;
                return;
            }
        }
        //ok = true;

        else if (job < 30) {

            iok = 0;
            int nsp = phase_nSpecies(ph);
            double* x = new double[nsp];
            switch (job) {
            case 20:
                iok = phase_getMoleFractions(ph,nsp,x);
                break;
            case 21:
                iok = phase_getMassFractions(ph,nsp,x);
                break;
            case 22:
                iok = phase_getMolecularWeights(ph,nsp,x);
                break;
            default:
                ;
            }
            plhs[0] = mxCreateNumericMatrix(nsp,1,
                mxDOUBLE_CLASS,mxREAL);
            double *h = mxGetPr(plhs[0]);
            if (iok >= 0) {
                for (int i = 0; i < nsp; i++) h[i] = x[i];
                delete x;
                return;
            }
            else {
                for (int i = 0; i < nsp; i++) h[i] = -999.99;
                delete x;
                mexErrMsgTxt("unknown attribute");
                return;
            }
        }

        else if (job < 40) {

            iok = 0;
            int nel = phase_nElements(ph);
            double* x = new double[nel];
            switch (job) {
            case 30:
                iok = phase_getAtomicWeights(ph,nel,x);
                break;
            default:
                ;
            }
            plhs[0] = mxCreateNumericMatrix(nel,1,
                mxDOUBLE_CLASS,mxREAL);
            double *h = mxGetPr(plhs[0]);
            if (iok >= 0) {
                for (int i = 0; i < nel; i++) h[i] = x[i];
                delete x;
                return;
            }
            else {
                for (int i = 0; i < nel; i++) h[i] = -999.99;
                delete x;
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
            }
            else {
                mexErrMsgTxt("error or unknown method.");
                reportError();
                return;
            }
        }
        else {
            mexErrMsgTxt("unimplemented method.");
            return;
        }
        if (iok < 0) reportError();
    }
