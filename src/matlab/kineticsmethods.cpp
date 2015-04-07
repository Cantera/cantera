#include "ctmatutils.h"
#include "clib/ct.h"

void checkNArgs(const int n, const int nrhs)
{
    if (n != nrhs) {
        mexErrMsgTxt("Wrong number of arguments.");
    }
}

void kineticsmethods(int nlhs, mxArray* plhs[],
                     int nrhs, const mxArray* prhs[])
{
    double vv = 0.0;
    int job = getInt(prhs[2]);
    int kin, irxn;

    // construct a new instance
    if (job == 0) {
        checkNArgs(8, nrhs);
        int root = getInt(prhs[1]);
        int iph = getInt(prhs[3]);
        int in1 = getInt(prhs[4]);
        int in2 = getInt(prhs[5]);
        int in3 = getInt(prhs[6]);
        int in4 = getInt(prhs[7]);
        vv = static_cast<int>(newKineticsFromXML(root, iph, in1, in2, in3, in4));
        plhs[0] = mxCreateNumericMatrix(1,1,mxDOUBLE_CLASS,mxREAL);
        double* h = mxGetPr(plhs[0]);
        *h = vv;
        return;
    }

    // methods
    else if (job > 0) {
        int isp = 1;
        if (job < 5 || job > 6) {
            checkNArgs(4,nrhs);
        } else {
            checkNArgs(5,nrhs);
            isp = getInt(prhs[4]);
        }
        kin = getInt(prhs[1]);
        irxn = getInt(prhs[3]);

        // get scalar attributes
        if (job < 10) {

            switch (job) {

            case 1:
                vv = (double) kin_nReactions(kin);
                break;
            case 2:
                vv = kin_multiplier(kin, irxn-1);
                break;
            case 3:
                vv = (double) kin_nSpecies(kin);
                break;
            case 4:
                vv = kin_isReversible(kin,irxn-1);
                break;
            case 5:
                vv = kin_reactantStoichCoeff(kin, isp - 1, irxn-1);
                break;
            case 6:
                vv = kin_productStoichCoeff(kin, isp - 1, irxn-1);
                break;
            default:
                mexErrMsgTxt("unknown job");
            }
            plhs[0] = mxCreateNumericMatrix(1,1,mxDOUBLE_CLASS,mxREAL);
            double* h = mxGetPr(plhs[0]);
            *h = vv;
            return;
        } else if (job < 20) {

            // get reaction array attributes
            mwSize nr = (mwSize) kin_nReactions(kin);
            plhs[0] = mxCreateNumericMatrix(nr,1,mxDOUBLE_CLASS,mxREAL);
            double* h = mxGetPr(plhs[0]);
            int ok = -10;
            switch (job) {
            case 11:
                ok = kin_getFwdRatesOfProgress(kin,nr,h);
                break;
            case 12:
                ok = kin_getRevRatesOfProgress(kin,nr,h);
                break;
            case 13:
                ok = kin_getNetRatesOfProgress(kin,nr,h);
                break;
            case 14:
                ok = kin_getEquilibriumConstants(kin,nr,h);
                break;
            default:
                ;
            }
            if (ok < 0) {
                mexErrMsgTxt("error computing rates of progress");
            }
        } else if (job < 30) {
            mwSize nsp = (mwSize) kin_nSpecies(kin);
            plhs[0] = mxCreateNumericMatrix(nsp,1,mxDOUBLE_CLASS,mxREAL);
            double* h = mxGetPr(plhs[0]);
            int ok = -10;
            switch (job) {
            case 21:
                ok = kin_getCreationRates(kin,nsp,h);
                break;
            case 22:
                ok = kin_getDestructionRates(kin,nsp,h);
                break;
            case 23:
                ok = kin_getNetProductionRates(kin,nsp,h);
                break;
            case 24:
                ok = kin_getSourceTerms(kin, nsp, h);
                break;
            default:
                ;
            }
            if (ok < 0) {
                mexErrMsgTxt("error computing production rates");
            }
        } else if (job < 40) {
            char* buf;
            int iok = -1, buflen = 80;
            switch (job) {
            case 31:
                buf = (char*)mxCalloc(buflen, sizeof(char));
                iok = kin_getReactionString(kin, irxn-1, buflen, buf);
                break;
            default:
                ;
            }
            if (iok >= 0) {
                plhs[0] = mxCreateString(buf);
                return;
            } else {
                reportError();
            }
        }
    }

    else {

        // set attributes
        int iok = -1;
        job = -job;
        kin = getInt(prhs[1]);
        irxn = getInt(prhs[3]);

        if (job < 10) {
            checkNArgs(5,nrhs);
            double v = getDouble(prhs[4]);
            switch (job) {
            case 1:
                iok = kin_setMultiplier(kin,irxn-1,v);
                break;
            case 3:
                iok = delKinetics(kin);
                break;
            case 5:
                iok = kin_advanceCoverages(kin,v);
                break;
            default:
                mexErrMsgTxt("unknown job");
            }
        }

        if (iok < 0) {
            reportError();
        }
    }
}
