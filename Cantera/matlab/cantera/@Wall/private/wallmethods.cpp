
#include "mex.h"
#include "../../../../clib/src/ctreactor.h"
#include "../../../../clib/src/ct.h"
#include "../../private/ctmatutils.h"

const double Undef = -999.123;

extern "C" {


    void reportError() {
        int buflen = 300;
        char* output_buf = (char*)mxCalloc(buflen, sizeof(char));
        getCanteraError(buflen, output_buf);
        mexErrMsgTxt(output_buf);
    }


    void mexFunction( int nlhs, mxArray *plhs[],
        int nrhs, const mxArray *prhs[] )
    {
        int j, m, iok, n;
        char *file, *key, *val;

        int job = getInt(prhs[0]);
        int i = getInt(prhs[1]);

        double r = Undef;
        double v = Undef;
        if (nrhs > 2) v = getDouble(prhs[2]);

        // constructor
        if (job == 0) {
            n = wall_new(i);
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
                iok = wall_del(i);
                break;
            case 2:
                iok = wall_copy(i);
                break;
            case 3:
                iok = wall_assign(i,int(v));
                break;
            case 4:
                m = getInt(prhs[3]);
                iok = wall_install(i, int(v), m);
                break;
            case 5:
                iok = wall_setArea(i, v);
                break;
            case 6:
                iok = wall_setThermalResistance(i, v);
                break;
            case 7:
                iok = wall_setHeatTransferCoeff(i, v);
                break;
            case 8:
                iok = wall_setHeatFlux(i, int(v));
                break;
            case 9:
                iok = wall_setExpansionRateCoeff(i, v);
                break;
            case 10:
                iok = wall_setExpansionRate(i, int(v));
                break;
            case 11:
                iok = wall_ready(i);
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
                r = wall_vdot(i, v);
                break;
            case 22:
                r = wall_Q(i, v);
                break;
            case 23:
                r = wall_area(i);
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
}
