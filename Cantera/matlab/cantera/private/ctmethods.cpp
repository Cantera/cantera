
#include "mex.h"
#include "../../../clib/src/ct.h"
#include "ctmatutils.h"

const int NO_CLASS = 0;
const int XML_CLASS = 10;
const int THERMO_CLASS = 20;
const int PHASE_CLASS = 30;
const int KINETICS_CLASS = 40;
const int TRANSPORT_CLASS = 50;
const int REACTOR_CLASS = 60;
const int WALL_CLASS = 70;
const int FLOWDEVICE_CLASS = 80;

void ctfunctions( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] );
void xmlmethods( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] );
void thermomethods( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] );
void phasemethods( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] );
void kineticsmethods( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] );
void transportmethods( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] );
void reactormethods( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] );
void wallmethods( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] );
void flowdevicemethods( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] );

extern "C" {

    void mexFunction( int nlhs, mxArray *plhs[],
        int nrhs, const mxArray *prhs[] )
    {
        int iclass = getInt(prhs[0]);

        switch (iclass) {
        case NO_CLASS:
            ctfunctions(nlhs, plhs, nrhs, prhs); break;
        case XML_CLASS:
            xmlmethods(nlhs, plhs, nrhs, prhs); break;
        case THERMO_CLASS:
            thermomethods(nlhs, plhs, nrhs, prhs); break;
        case PHASE_CLASS:
            phasemethods(nlhs, plhs, nrhs, prhs); break;
        case KINETICS_CLASS:
            kineticsmethods(nlhs, plhs, nrhs, prhs); break;
        case TRANSPORT_CLASS:
            transportmethods(nlhs, plhs, nrhs, prhs); break;
        case REACTOR_CLASS:
            reactormethods(nlhs, plhs, nrhs, prhs); break;
        case WALL_CLASS:
            wallmethods(nlhs, plhs, nrhs, prhs); break;
        case FLOWDEVICE_CLASS:
            flowdevicemethods(nlhs, plhs, nrhs, prhs); break;
        default:
            mexErrMsgTxt("unknown class");
        }
    }
}
